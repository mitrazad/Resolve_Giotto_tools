#' Convert Resolve Seurat object to Giotto object
#'
#' @param resolve.obj Input Resolve Seurat object to convert to Giotto object
#' @param spatial_assay Specify name of the spatial assay slot in Seurat. Default is \code{"Spatial"}.
#' @param dim_reduction Specify which dimensional reduction computations to fetch from
#'  input Seurat object. Default is \code{"c('pca', 'umap')"}.
#' @param subcellular_assay Specify name of the subcellular assay in input object. Default is \code{"Vizgen"}.
#' @param instructions Specify instructions during the \code {createGiottoInstructions} step before. Default is \code{NULL}.
#' @return A Giotto object converted from Seurat object with all computations stored in it.
#' @references Adapted from "https://github.com/drieslab/Giotto_site_suite"
#' @export

ResolveSeuratToGiotto <- function(resolve.obj,
                                  spatial_assay = "Spatial",
                                  dim_reduction = c("pca", "umap"),
                                  subcellular_assay = "Vizgen",
                                  instructions = NULL) {
  PackageCheck("Seurat")
  
  if(is.null(GetAssayData(object = resolve.obj, slot = "counts", assay = spatial_assay))) {
    wrap_msg("No raw expression values are provided in spatial_assay")
    return(resolve.obj)
    
  } else {
    
    exp = GetAssayData(object = resolve.obj, slot = "counts", assay = spatial_assay) 
    if(!is.null(resolve.obj@assays$SCT)){
      normexp = GetAssayData(object = resolve.obj, slot = "counts", assay = "SCT")
    }
    
    if(!is.null(slot(resolve.obj, "assays")[[spatial_assay]]@data)){
      normexp = GetAssayData(object = resolve.obj, slot = "data", assay = spatial_assay)
    }
    
    # Cell Metadata
    cell_metadata = resolve.obj@meta.data    
    
    # Dimension Reduction
    if(sum(sapply(dim_reduction,function(x) length(resolve.obj@reductions[[x]]))) == 0) {
      dimReduc_list = NULL
    } else {
      dimReduc_list = lapply(dim_reduction,function(x){
        dim_coord = as.matrix(Seurat::Embeddings(object = resolve.obj,reduction = x))
        dim_load = as.matrix(Seurat::Loadings(object = resolve.obj, reduction = x))
        dim_eig = Seurat::Stdev(resolve.obj, reduction = x)
        if (length(dim_eig) > 0){
          dim_eig = sapply(dim_eig, function(x) x ^ 2)
        }
        colnames(dim_coord) = paste0("Dim.",1:ncol(dim_coord))
        if (length(dim_load) > 0){
          colnames(dim_load) = paste0("Dim.",1:ncol(dim_load))
        }
        dimReduc = Giotto:::create_dim_obj(name = x,
                                           reduction = "cells",
                                           reduction_method = x,
                                           spat_unit = "cell",
                                           feat_type = "rna",
                                           provenance = NULL,
                                           coordinates = dim_coord,
                                           misc = list(eigenvalues = dim_eig, loadings = dim_load))
        
        return(dimReduc)
      })
      names(dimReduc_list) <- dim_reduction
    }
    
    # Polygons: requires images object
    if (!length(resolve.obj@images[[subcellular_assay]]@boundaries$segmentation@polygons) == 0) {
      spat_coord <- Seurat::GetTissueCoordinates(resolve.obj)
      spat_coord <- spat_coord[,c(3,1,2)]
      colnames(spat_coord) = c("poly_ID", "x", "y")
      gpolygons <- createGiottoPolygonsFromDfr(spat_coord, name = "cell")
    }
    
    # Subcellular/ Centroids
    if(length(resolve.obj@assays[[spatial_assay]]) == 0) {
      spat_loc = NULL
    } else {
      name = names(resolve.obj@images)
      spat_coord = cbind(data.frame(cell_id = resolve.obj@images[[subcellular_assay]]@boundaries$centroids@cells),
                         resolve.obj@images[[subcellular_assay]]@boundaries$centroids@coords)
      colnames(spat_coord) = c("cell_id", "sdimx", "sdimy")
      exp = exp[, c(intersect(spat_coord$cell_id, colnames(exp)))]
      spat_loc = spat_coord
    }
    
    
    if (!length(resolve.obj@images) == 0) {
      if ("molecules" %in% methods::slotNames(resolve.obj@images[[name]]) == TRUE) {
        if (!length(resolve.obj@images[[name]][["molecules"]]) == 0) {
          featnames = rownames(resolve.obj@assays[[spatial_assay]]@meta.features)
          mol_spatlocs = data.table::data.table()
          for (x in featnames) {
            df = (Seurat::FetchData(resolve.obj[[name]][["molecules"]],vars = x))
            mol_spatlocs = rbind(mol_spatlocs, df)
          }
          
          gpoints = createGiottoPoints(mol_spatlocs,feat_type = "rna")
        }
      }
    }
  }
  
  
  gobject = createGiottoObject(expression = exp,
                               spatial_locs = spat_loc,
                               dimension_reduction = dimReduc_list,
                               instructions = NULL)
  
  if (exists("normexp") == TRUE) {
    exprObj = createExprObj(expression_data = normexp,
                            name = "normalized",
                            spat_unit = "cell",
                            feat_type = "rna",
                            provenance = "cell")
    gobject = set_expression_values(gobject = gobject, values = exprObj, set_defaults = FALSE)
  }
  gobject = addCellMetadata(gobject = gobject, new_metadata = cell_metadata)
  
  if (exists("gpoints") == TRUE) {
    gobject = addGiottoPoints(gobject = gobject, gpoints = list(gpoints))
  }
  
  if(exists("gpolygons") == TRUE) {
    gobject = addGiottoPolygons(gobject = gobject, gpolygons = list(gpolygons))
  }
  
  return (gobject)
}