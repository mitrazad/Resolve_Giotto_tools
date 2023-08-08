# Part 03 of Experiment_01_RStG.R
rm(list = ls())
# Ensure Giotto Suite is installed.
if(!"Giotto" %in% installed.packages()) {
  devtools::install_github("drieslab/Giotto@suite")
}

# Ensure GiottoData, a small, helper module for tutorials, is installed.
if(!"GiottoData" %in% installed.packages()) {
  devtools::install_github("drieslab/GiottoData")
}
library(Giotto)
# Ensure the Python environment for Giotto has been installed.
genv_exists = checkGiottoEnvironment()
if(!genv_exists){
  # The following command need only be run once to install the Giotto environment.
  # installGiottoEnvironment()
}

# Enter commands in R (or R studio, if installed)
# Install the remotes package
install.packages('remotes')
remotes::install_github(repo = 'satijalab/seurat', ref = 'develop')
library(Seurat)
library(future)
plan("multisession", workers = 10)
library(data.table)
library(RImageJROI)
library(magick)

# Set working directory
results_folder = '~/Documents/Resolve/Resolve Analysis/'

# Upload resolve object
resolve.obj <- readRDS("~/Documents/Resolve/Resolve Analysis/cellpose.dragonfly.W6A1.rds")

# Michele: Quick test single cell analysis, not meaningful just to have something to plot
# Remove unwanted cells = skip cells with no transcripts
resolve.obj <- subset(resolve.obj, subset = nCount_Resolve > 10) 

# Apply a centered log ratio transformation
# resolve.obj <- NormalizeData(object = resolve.obj, normalization.method = "CLR", margin = 2)

# Standard pre-processing step to PCA and UMAP
# resolve.obj <- ScaleData(resolve.obj)

# Interpretation of scRNA-seq data requires effective pre-processing and normalization.
# SCTransform = modeling framework for normalization and variance stabilization of molecular count data from scRNA-seq experiment 
resolve.obj <- SCTransform(resolve.obj, assay = "Resolve", verbose = FALSE)

resolve.obj <- FindVariableFeatures(resolve.obj, assay = "Resolve")

resolve.obj <- RunPCA(resolve.obj, npcs = 20, verbose = FALSE)

resolve.obj <- RunUMAP(resolve.obj, dims = 1:20, verbose = FALSE)

resolve.obj <- FindNeighbors(resolve.obj, k.param = 20, dims = 1:20, verbose = TRUE)

resolve.obj <- FindClusters(resolve.obj, verbose = FALSE, resolution = 0.6, algorithm = 4)

# From seurat to Giotto (only useful if we want to use Seurat before,
# if not we can just make a giotto object from scratch)
# See: https://github.com/RubD/Giotto/blob/suite/R/interoperability.R
# seuratToGiotto

# The SeuratToGiotto function supports only one Assay and image in the Seurat object
# resolve.obj@assays <- resolve.obj@assays["Resolve"]
# resolve.obj@images <- resolve.obj@images["cellpose.dragonfly.W6A1"]

# Adjust the Seurat object so it works with SeuratToGiotto
ResolveSeuratToGiotto <- function(resolve.obj,
                                  spatial_assay = "Resolve",
                                  dim_reduction = c('pca','umap'),
                                  subcellular_assay = "cellpose.dragonfly.W6A1") {
  Giotto:::package_check('Seurat')
  
  if(is.null(Seurat::GetAssayData(object = resolve.obj, slot = "counts", assay = spatial_assay))) {
    wrap_msg('No raw expression values are provided in spatial_assay')
    return(resolve.obj)
    
  } else {
    
    exp = Seurat::GetAssayData(object = resolve.obj, slot = "counts", assay = spatial_assay) 
    if(!is.null(resolve.obj@assays$SCT)){
      normexp = Seurat::GetAssayData(object = resolve.obj, slot = "counts", assay = 'SCT')
    }
    
    if(!is.null(slot(resolve.obj, 'assays')[[spatial_assay]]@data)){
      normexp = Seurat::GetAssayData(object = resolve.obj, slot = "data", assay = spatial_assay)
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
        colnames(dim_coord) = paste0('Dim.',1:ncol(dim_coord))
        if (length(dim_load) > 0){
          colnames(dim_load) = paste0('Dim.',1:ncol(dim_load))
        }
        dimReduc = Giotto:::create_dim_obj(name = x,
                                           reduction = 'cells',
                                           reduction_method = x,
                                           spat_unit = 'cell',
                                           feat_type = 'rna',
                                           provenance = NULL,
                                           coordinates = dim_coord,
                                           misc = list(eigenvalues = dim_eig, loadings = dim_load))
        
        return(dimReduc)
      })
      names(dimReduc_list) <- dim_reduction
    }
    
    #Spatial locations
    if (length(resolve.obj@assays[[spatial_assay]]) == 0) {
      spat_loc = NULL
    } else {
      
      # requires image objects!
      spat_coord = Seurat::GetTissueCoordinates(resolve.obj)
      spat_coord <- spat_coord[,c(3,1,2)]
      colnames(spat_coord) = c("cell_ID", "sdimx", "sdimy")
      spat_loc = spat_coord
    }
    
    #Subcellular
    name = names(resolve.obj@images)
    if (length(resolve.obj@assays[[subcellular_assay]]) == 1) {
      spat_coord = Seurat::GetTissueCoordinates(resolve.obj)
      colnames(spat_coord) = c("sdimx", "sdimy", "cell_ID")
      spat_coord <- spat_coord[,c(3,1,2)]
      exp = exp[, c(intersect(spat_coord$cell_ID, colnames(exp)))]
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
  
  gobject = createGiottoObject(exp,
                               spatial_locs = spat_loc,
                               dimension_reduction = dimReduc_list)
  if (exists('normexp') == TRUE) {
    exprObj = Giotto:::create_expr_obj(name = 'normalized',
                                       exprMat = normexp,
                                       spat_unit = 'cell',
                                       feat_type = 'rna',
                                       provenance = 'cell')
    gobject = set_expression_values(gobject = gobject, values = exprObj, set_defaults = FALSE)
    # gobject@expression$cell$rna$normalized = normexp
  }
  gobject = addCellMetadata(gobject = gobject, new_metadata = cell_metadata)
  
  if (exists('gpoints') == TRUE) {
    gobject = addGiottoPoints(gobject = gobject,gpoints = list(gpoints))
  }
  return (gobject)
}

resolve.obj <- ResolveSeuratToGiotto(resolve.obj)

# UMAP
plotUMAP(resolve.obj)

#PCA
plotPCA(resolve.obj)

# Duplicate original gobject
resolve.obj_sNN = resolve.obj
resolve.obj_sNN <-  createNearestNetwork(resolve.obj, dimensions_to_use = 1:20, k = 20)
## Leiden clustering not available
# visium_brain <- doLeidenCluster(resolve.obj_sNN, resolution = 0.6, n_iterations = 1000)




# Part 02 of Experiment_01_RStG.R
rm(list = ls())
library(Seurat)
library(future)
plan("multisession", workers = 10)
library(data.table)
library(RImageJROI)
library(Giotto)
library(GiottoData)
library(magick)


resolve.obj <- readRDS("~/Documents/Resolve/Resolve Analysis/cellpose.dragonfly.W6A1.rds")

################
# From seurat to Giotto (only useful if we want to use Seurat before,
# if not we can just make a giotto object from scratch)
# See: https://github.com/RubD/Giotto/blob/suite/R/interoperability.R
# seuratToGiotto

# The SeuratToGiotto function supports only one Assay and image in the Seurat object
resolve.obj@assays <- resolve.obj@assays["Resolve"]
resolve.obj@images <- resolve.obj@images["cellpose.dragonfly.W6A1"]

ResolveSeuratToGiotto <- function(resolve.obj) {
  
  if (length(resolve.obj@assays[[spatial_assay]]) == 0) {
    spat_loc = NULL
  } else {
    
    # requires image objects!
    spat_coord = Seurat::GetTissueCoordinates(resolve.obj)
    spat_coord = cbind(rownames(spat_coord), data.frame(spat_coord, 
                                                        row.names = NULL))
    colnames(spat_coord) = c("cell_ID", "sdimy", "sdimx")
    spat_coord <- spat_coord[, -1]
    spat_coord <- spat_coord[,c(3,2,1)]
    colnames(spat_coord) = c("cell_ID", "sdimx", "sdimy")
    spat_loc = spat_coord
  }
  name = names(resolve.obj@images)
  if (length(resolve.obj@assays[[subcellular_assay]]) == 1) {
    spat_coord = Seurat::GetTissueCoordinates(resolve.obj)
    colnames(spat_coord) = c("sdimx", "sdimy", "cell_ID")
    spat_coord <- spat_coord[,c(3,1,2)]
    exp = exp[, c(intersect(spat_coord$cell_ID, colnames(exp)))]
    spat_loc = spat_coord
  }
  if (!length(resolve.obj@images) == 0) {
    if ("molecules" %in% methods::slotNames(resolve.obj@images[[name]]) == 
        TRUE) {
      if (!length(resolve.obj@images[[name]][["molecules"]]) == 
          0) {
        assay = names(resolve.obj@assays)
        featnames = rownames(resolve.obj@assays[[assay]]@meta.features)
        mol_spatlocs = data.table::data.table()
        for (x in featnames) {
          df = (Seurat::FetchData(resolve.obj[[name]][["molecules"]], 
                                  vars = x))
          mol_spatlocs = rbind(mol_spatlocs, df)
        }
        gpoints = createGiottoPoints(mol_spatlocs, 
                                     feat_type = "rna")
      }
    }
  }
}

# Part 01 of Experiment_01_RStG.R
giotto_res <- seuratToGiotto(resolve.obj, spatial_assay = "Resolve",
                             subcellular_assay =  "cellpose.dragonfly.W6A1")

# Adjust the Seurat object so it works with SeuratToGiotto
resolve.obj@images$cellpose.dragonfly.W6A1@boundaries$segmentation

ResolveSeuratToGiotto <- function(resolve.obj) {
  
  if (length(resolve.obj@assays[[spatial_assay]]) == 0) {
    spat_loc = NULL
  } else {
    
    # requires image objects!
    spat_coord = Seurat::GetTissueCoordinates(resolve.obj)
    spat_coord = cbind(rownames(spat_coord), data.frame(spat_coord, 
                                                        row.names = NULL))
    colnames(spat_coord) = c("cell_ID", "sdimy", "sdimx")
    spat_coord <- spat_coord[, -1]
    spat_coord <- spat_coord[,c(3,2,1)]
    colnames(spat_coord) = c("cell_ID", "sdimx", "sdimy")
    spat_loc = spat_coord
  }
  name = names(resolve.obj@images)
  if (length(resolve.obj@assays[[subcellular_assay]]) == 1) {
    spat_coord = Seurat::GetTissueCoordinates(resolve.obj)
    colnames(spat_coord) = c("sdimx", "sdimy", "cell_ID")
    spat_coord <- spat_coord[,c(3,1,2)]
    exp = exp[, c(intersect(spat_coord$cell_ID, colnames(exp)))]
    spat_loc = spat_coord
  }
  if (!length(resolve.obj@images) == 0) {
    if ("molecules" %in% methods::slotNames(resolve.obj@images[[name]]) == 
        TRUE) {
      if (!length(resolve.obj@images[[name]][["molecules"]]) == 
          0) {
        assay = names(resolve.obj@assays)
        featnames = rownames(resolve.obj@assays[[assay]]@meta.features)
        mol_spatlocs = data.table::data.table()
        for (x in featnames) {
          df = (Seurat::FetchData(resolve.obj[[name]][["molecules"]], 
                                  vars = x))
          mol_spatlocs = rbind(mol_spatlocs, df)
        }
        gpoints = createGiottoPoints(mol_spatlocs, 
                                     feat_type = "rna")
      }
    }
  }
}

gobject = createGiottoObject(exp, spatial_locs = spat_loc, 
                             dimension_reduction = dimReduc_list)
if (exists("normexp") == TRUE) {
  exprObj = create_expr_obj(name = "normalized", exprMat = normexp, 
                            spat_unit = "cell", feat_type = "rna", provenance = "cell")
  gobject = set_expression_values(gobject = gobject, values = exprObj, 
                                  set_defaults = FALSE)
}
gobject = addCellMetadata(gobject = gobject, new_metadata = cell_metadata)

if (exists("gpoints") == TRUE) {
  gobject = addGiottoPoints(gobject = gobject, gpoints = list(gpoints))
}
return(gobject)

# Giotto object with the included row name changes of lane 54 in the debugger
# gobject <- readRDS("gobject.rds")

# Centroids are sufficient, giotto does not support ROIs
giotto_res <- seuratToGiotto(resolve.obj, spatial_assay = "Resolve",
                             subcellular_assay =  "cellpose.dragonfly.W6A1")

# Extra notes

{
  package_check('Seurat')
  
  if(is.null(Seurat::GetAssayData(object = resolve.obj, slot = "counts", assay = spatial_assay))) {
    wrap_msg('No raw expression values are provided in spatial_assay')
    return(resolve.obj)
  }
  else {
    exp = Seurat::GetAssayData(object = resolve.obj, slot = "counts", assay = spatial_assay)
    if(!is.null(resolve.obj@assays$Resolve)){
      normexp = Seurat::GetAssayData(object = resolve.obj, slot = "counts", assay = 'Resolve')
    }
    if(!is.null(slot(resolve.obj, 'assays')[[spatial_assay]]@data)) {
      normexp = Seurat::GetAssayData(object = resolve.obj, slot = "data", assay = spatial_assay)
    }
    
    # Spatial Locations
    spatialLocs <- SpatialExperiment::spatialCoords(spe)
    if(ncol(spatialLocs) > 0){
      if(verbose) message("Copying spatial locations")
      spatialLocsDT <- data.table(sdimx = spatialLocs[, 1], sdimy = spatialLocs[, 2], cell_ID = rownames(spatialLocs))
      giottoObj <- set_spatial_locations(gobject = giottoObj, spatlocs = cbind(spatialLocsDT, cell_ID = colnames(spe)))
    }
    
    # Spatial Images
    spatialImages <- SpatialExperiment::imgData(spe)
    if(nrow(spatialImages) > 0){
      for(i in seq(nrow(spatialImages))){
        if(verbose) message("Copying spatial images")
        spImg <- SpatialExperiment::getImg(spe,
                                           spatialImages[i, "sample_id"],
                                           spatialImages[i, "image_id"])
        mObject <- magick::image_read(grDevices::as.raster(spImg))
        giottoImage <- createGiottoImage(gobject = giottoObj,
                                         mg_object = mObject,
                                         scale_factor = spatialImages[i, "scaleFactor"])
        giottoObj <- addGiottoImage(gobject = giottoObj,
                                    images = list(giottoImage))
        
        giotto_img <- createGiottoImage(giotto_res, spatial_locs = NULL,
                                        mg_object = DAPI_img, name = "image", xmax_adj = 0, xmin_adj = 0,
                                        ymax_adj = 0, ymin_adj = 0, scale_factor = DAPI_scale_factor,
                                        do_manual_adj = T)