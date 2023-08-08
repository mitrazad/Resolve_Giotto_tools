rm(list = ls())
# Ensure Giotto Suite is installed.
#if(!"Giotto" %in% installed.packages()) {
#  devtools::install_github("drieslab/Giotto@suite")
#}

# Ensure GiottoData, a small, helper module for tutorials, is installed.
#if(!"GiottoData" %in% installed.packages()) {
#  devtools::install_github("drieslab/GiottoData")
#}
library(Giotto)
library(GiottoData)
# Ensure the Python environment for Giotto has been installed.
genv_exists = checkGiottoEnvironment()
#if(!genv_exists){
  # The following command need only be run once to install the Giotto environment.
  # installGiottoEnvironment()
#}

library(remotes)
# remotes::install_github(repo = 'satijalab/seurat', ref = 'develop', force = TRUE)
library(Seurat)
library(future)
plan("multisession", workers = 10)
library(data.table)
library(RImageJROI)
library(magick)

# Set working directory
data_path = '~/Giotto_loader_Resolve/'

# Upload resolve object
resolve.obj <- readRDS("~/Giotto_loader_Resolve/Cellpose.dragonfly.W6A1/cellpose.dragonfly.W6A1.rds")

# From seurat to Giotto (only useful if we want to use Seurat before,
# if not we can just make a giotto object from scratch)
# See: https://github.com/RubD/Giotto/blob/suite/R/interoperability.R
# seuratToGiotto

#resolve.obj <- SCTransform(resolve.obj, assay = "Resolve", verbose = FALSE)

#resolve.obj <- FindVariableFeatures(resolve.obj, assay = "Resolve")

#resolve.obj <- RunPCA(resolve.obj, npcs = 20, verbose = FALSE)

#npcs <- 20
#n_neighbors <- c(10, 30, 50, 100)
#min_dist <- c(0.2, 0.3, 0.4, 0.6)

#for (n in n_neighbors) {
  #for (md in min_dist) {
    #resolve.obj <- RunUMAP(resolve.obj,dims = 1:npcs, verbose = TRUE, seed.use = 42, n.neighbors = n, min.dist = md)
    # DimPlot(resolve.obj, reduction = "umap", pt.size = 0.2, group.by = "orig.BioClassification") + ggtitle(paste("UMAP with", n, "neighbors and min.dist =", md))
    
    # ggsave(filename = paste0("UMAP_", n, "_neighbors_", md, "_min_dist.png"))
#  }
#}
#resolve.obj <- RunUMAP(resolve.obj, dims = 1:20, verbose = FALSE)

#resolve.obj <- FindNeighbors(resolve.obj, dims = 1:20, verbose = TRUE)

#resolve.obj <- FindClusters(resolve.obj, verbose = FALSE, resolution = 0.6)

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

# Duplicate resolve.obj
resolve.obj_red <- resolve.obj

## Normalize Giotto object
resolve.obj_red <- normalizeGiotto(gobject = resolve.obj_red, scalefactor = 6000, verbose = T)

## add gene & cell statistics
resolve.obj_red <- addStatistics(gobject = resolve.obj_red)

## Adjust Giotto matrix
resolve.obj_red <- adjustGiottoMatrix(gobject = resolve.obj_red,
                                expression_values = c('normalized'),
                                covariate_columns = c('nr_feats', 'total_expr'))
## New method?
resolve.obj_red <- normalizeGiotto(gobject = resolve.obj_red, norm_methods = 'pearson_resid', update_slot = 'pearson')

showGiottoExpression(resolve.obj_red)

# PCA
## highly variable features / genes (HVF)
resolve.obj_red <- calculateHVF(resolve.obj_red, HVFname = 'hvg_orig')

# new method based on variance of pearson residuals for each gene
subc_test <- calculateHVF(gobject = resolve.obj_red,
                          method = 'var_p_resid', expression_values = 'pearson',
                          show_plot = T)

## run PCA on expression values (default)
gene_metadata = fDataDT(resolve.obj_red)
featgenes = gene_metadata[hvf == 'yes']$feat_ID

## run PCA on expression values (default)
resolve.obj_red <- runPCA(resolve.obj_red, feats_to_use = featgenes)
screePlot(resolve.obj_red, ncp = 20)
dimPlot2D(resolve.obj_red, dim_reduction_to_use = "pca")

# UMAP
## run UMAP on PCA space
resolve.obj_red <- runUMAP(resolve.obj_red, dimensions_to_use = 1:11)
plotUMAP(resolve.obj_red)

resolve.obj_red <-  createNearestNetwork(resolve.obj_red, dimensions_to_use = 1:11, k = 20)
## Leiden clustering not available
# visium_brain <- doLeidenCluster(resolve.obj_red, resolution = 0.6, n_iterations = 1000)

# If set to NULL (default), the Python executable within the previously installed Giotto environment will be used.
my_python_path = NULL

# Create Giotto Instructions
instrs = createGiottoInstructions(save_dir = data_path,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  python_path = my_python_path)

original_DAPI_image = paste0("Documents/Resolve/Giotto_loader_Resolve/Cellpose.dragonfly.W6A1/MAX_Resolve_MBENrun1_pseudoHE_C2-1_2022-07-19_13.48.59_Fused_405nm_corr_cut_registered.png")

segmentation_mask = paste0('Documents/Resolve/Giotto_loader_Resolve/Cellpose.dragonfly.W6A1/MAX_Resolve_MBENrun1_pseudoHE_C2-1_2022-07-19_13.48.59_Fused_405nm_corr_cut_registered_cp_masks.pngg')

DAPI_image = createGiottoImage(gobject = resolve.obj_red,
                               name = 'DAPI',
                               do_manual_adj = T,
                               xmax_adj = 0,ymax_adj = 0,
                               xmin_adj = 0,ymin_adj = 0,
                               image_transformations = c("flip_x_axis", "flip_y_axis"),
                               mg_object = original_DAPI_image)

resolve.obj_red = addGiottoImage(resolve.obj_red,
                         images = list(DAPI_image))


spatPlot2D(gobject = resolve.obj_red, image_name = 'DAPI', point_size = 0.1)

# nr_feats = numerical features of cells
spatDimPlot(resolve.obj_red, cell_color = 'nr_feats', color_as_factor = F,
            dim_point_size = 1, spat_point_size = 2.5)
