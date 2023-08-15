# Adapted from "https://giottosuite.readthedocs.io/en/latest/subsections/datasets/resolve_bc_210928.html"

rm(list = ls())

## Ensure Giotto Suite is installed.
# if(!"Giotto" %in% installed.packages()) {
#   devtools::install_github("drieslab/Giotto@suite")
#   }

library(Giotto)

## Ensure GiottoData, a small, helper module for tutorials, is installed.
# if(!"GiottoData" %in% installed.packages()) {
#   devtools::install_github("drieslab/GiottoData")
#   }

library(GiottoData)

## Ensure the Python environment for Giotto has been installed.
genv_exists = checkGiottoEnvironment()

## The following packages should be installed
# if(!genv_exists){
#   installGiottoEnvironment(packages_to_install = c("pandas", "networkx", "python-igraph","leidenalg", "python-louvain", "python.app", "scikit-learn", "umap-learn"),
#                            force_miniconda = FALSE,
#                            force_environment = TRUE,
#                            verbose = TRUE)
#   }

## If packages are already present, the following command need only be run once to install the Giotto environment.
# if(!genv_exists){
#   installGiottoEnvironment()
#   }

library(Seurat)
library(magick)
library(data.table)
library(RImageJROI)
library(future)
library(dplyr)
library(ggdendro)

#################### INPUT FILES ####################
## 1. Set working directory
results_folder = "~/path/to/result"

## 2. Load Resolve object
resolve.obj <- readRDS("~/resolve_file.rds")

## 3. Specify a path to a Python executable within a conda or miniconda environment
## If set to NULL (default), the Python executable within the previously installed Giotto environment will be used
my_python_path = NULL # alternatively, "/local/python/path/python" if desired.

## Optional: Provide path to Resolve Bioscience folder
## This might be needed if one prefers "paste0" over "image_read" when loading the original and segmentation images in 5. and 6.
# data_path = "~/path/to/Resolve/Bioscience/folder"

## 4. Create Giotto instructions
instrs = createGiottoInstructions(save_dir = results_folder,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  python_path = my_python_path)

## 5. Original image as png
# original_DAPI_image = paste0(data_path, "/", "DAPI_image.jpg")
original_DAPI_image = image_read("/path/to/DAPI_image.jpg")

## 6. Input cell segmentation as mask file
## Can also be provided as a 3-column polygon file
## Used as image background AND to store segmentations as polygons
## Can be obtained through Fiji / QuPath / Ilastik / Cellpose / ...
# segmentation_mask = paste0(data_path, "/", "segm_mask_image.jpg")
segmentation_mask = image_read("/path/to/segm_mask_image.jpg")

#################### PART 1: RESOLVE FILE: SEURAT TO GIOTTO OBJECT ####################

source("~/Giotto_loader_Resolve/ResolveSeurattoGiotto/RStG_function.R")

resolve.gobj <- ResolveSeuratToGiotto(resolve.obj,
                                      spatial_assay = "Spatial",
                                      dim_reduction = c("pca", "umap"),
                                      subcellular_assay = "Vizgen",
                                      instructions = instrs)

#################### PART 2: ADD IMAGE INFORMATION ####################

## Create Giotto images
DAPI_image = createGiottoImage(gobject = resolve.gobj,
                               name = "DAPI",
                               do_manual_adj = TRUE,
                               xmax_adj = 0,ymax_adj = 0,
                               xmin_adj = 0,ymin_adj = 0,
                               mg_object = original_DAPI)

segm_image = createGiottoImage(gobject = resolve.gobj,
                               name = "segmentation",
                               do_manual_adj = TRUE,
                               xmax_adj = 0,ymax_adj = 0,
                               xmin_adj = 0,ymin_adj = 0,
                               mg_object = segmentation_mask)

## Add images to Giotto object
resolve.gobj = addGiottoImage(resolve.gobj, 
                              images = list(DAPI_image, segm_image))

## This provides an overview of available images
showGiottoImageNames(resolve.gobj)

#################### PART 3: VISUALIZE ORIGINAL IMAGES ####################

## Visualize overlay of calculated cell centroids with original image and segmentation mask file

DAPI_spatplot <- spatPlot2D(gobject = resolve.gobj,
                            image_name = "DAPI",
                            point_size = 1.0,
                            show_image = TRUE)
# print(DAPI_spatplot)

segm_spatplot <- spatPlot2D(gobject = resolve.gobj,
                            image_name = "segmentation",
                            point_size = 1.0,
                            show_image = TRUE)
# print(segm_spatplot)

#################### PART 4: CALCULATE CELL SHAPE OVERLAP ####################

filt_resolve.gobj = calculateOverlapParallel(resolve.gobj,
                                             spatial_info = "cell",
                                             feat_info = "rna")

#################### PART 5: FILTER DATA ####################

## Features can be filtered individually
## Cells will be filtered across features

## First, filter on RNA
filt_resolve.gobj <- filterGiotto(gobject = filt_resolve.gobj,
                                  expression_threshold = 1,
                                  feat_det_in_min_cells = 20, 
                                  min_det_feats_per_cell = 5)

filt_spatplot <- spatPlot2D(gobject = filt_resolve.gobj, 
                            image_name = "segmentation",
                            point_size = 1.0,
                            show_image = TRUE)

# print(filt_spatplot)

## RNA data (default)
## Other feature modalities can be processed and filtered in an analogous manner

## Normalize Giotto object
filt_resolve.gobj <- normalizeGiotto(gobject = filt_resolve.gobj, 
                                     scalefactor = 6000, 
                                     verbose = TRUE)

## Add gene & cell statistics
filt_resolve.gobj <- addStatistics(gobject = filt_resolve.gobj)

## Adjust Giotto matrix
filt_resolve.gobj <- adjustGiottoMatrix(gobject = filt_resolve.gobj,
                                        expression_values = c("normalized"),
                                        covariate_columns = c("nr_feats", "total_expr"))

## Normalize according to Pearson residuals
filt_resolve.gobj <- normalizeGiotto(gobject = filt_resolve.gobj, 
                                     norm_methods = "pearson_resid", 
                                     update_slot = "pearson")

showGiottoExpression(filt_resolve.gobj)

#################### PART 6: DIMENSION REDUCTION ####################

## Find highly variable features / genes (HVF)
## Typical way of calculating HVF
filt_resolve.gobj <- calculateHVF(filt_resolve.gobj, 
                                  HVFname = "hvg_orig")

## New method based on variance of pearson residuals for each gene
filt_resolve.gobj <- calculateHVF(gobject = filt_resolve.gobj,
                                  method = "var_p_resid", 
                                  expression_values = "pearson",
                                  show_plot = TRUE)

## Run PCA
filt_resolve.gobj <- runPCA(gobject = filt_resolve.gobj,
                            expression_values = "pearson",
                            scale_unit = FALSE,
                            center = FALSE)

## Scree plot to determine the principal components to keep in PCA
screePlot(gobject = filt_resolve.gobj, 
          ncp = 20)

## Plot PCA
plotPCA(gobject = filt_resolve.gobj,
        dim1_to_use = 1,
        dim2_to_use = 2)

## Run UMAP
filt_resolve.gobj <- runUMAP(gobject = filt_resolve.gobj, 
                             dimensions_to_use = 1:5, 
                             n_threads = 2)
## Plot UMAP
plotUMAP(gobject = filt_resolve.gobj)

#################### PART 7: CLUSTER ####################

filt_resolve.gobj <- createNearestNetwork(gobject = filt_resolve.gobj, 
                                          dimensions_to_use = 1:5, 
                                          k = 5)

filt_resolve.gobj <- doLeidenCluster(gobject = filt_resolve.gobj, 
                                     resolution = 0.05, 
                                     n_iterations = 1000, 
                                     name = "leiden_clus")

## Visualize UMAP cluster results
plotUMAP(gobject = filt_resolve.gobj, 
         cell_color = "leiden_clus",
         show_NN_network = TRUE, 
         point_size = 2.5)

## Visualize UMAP and spatial results
spatDimPlot2D(gobject = filt_resolve.gobj,
              show_image = TRUE, 
              image_name = "segmentation",
              cell_color = "leiden_clus",
              spat_point_size = 2.5)

## Plot a cluster heatmap
showClusterHeatmap(gobject = filt_resolve.gobj, 
                   cluster_column = "leiden_clus",
                   save_param = list(save_format = "pdf",
                                     base_height = 6, 
                                     base_width = 8, 
                                     units = "cm"))

## See cluster relationships in a dendogram
showClusterDendrogram(gobject = filt_resolve.gobj, 
                      h = 0.5, 
                      rotate = TRUE, 
                      cluster_column = "leiden_clus")

#################### PART 8: CREATE A SPATIAL NETWORK ####################

filt_resolve.gobj = createSpatialNetwork(gobject = filt_resolve.gobj,
                                         spat_loc_name = "raw",
                                         minimum_k = 3,
                                         maximum_distance_delaunay = 100)

spatPlot2D(gobject = filt_resolve.gobj,
           image_name = "segmentation", 
           show_image = TRUE,
           point_size = 1.5, 
           show_network = TRUE)

#################### PART 9: VISUALIZE SUBCELLULAR DATA ####################

## Visualize clustered cells
spatInSituPlotPoints(gobject = filt_resolve.gobj,
                     show_polygon = TRUE,
                     polygon_feat_type = "cell",
                     polygon_color = "white",
                     polygon_line_size = 0.1,
                     polygon_fill = "leiden_clus",
                     polygon_fill_as_factor = TRUE)

## Individual plotting of transcripts and polygon information

## All cells
spatInSituPlotPoints(gobject = filt_resolve.gobj,
                     feats = list("rna" = c("LAMA2", "MKI67", "NRXN3")),
                     point_size = 0.2,
                     show_polygon = TRUE,
                     polygon_feat_type = "cell",
                     polygon_color = "white",
                     polygon_line_size = 0.1)

## Filtered cells
spatInSituPlotPoints(gobject = filt_resolve.gobj,
                     feats = list("rna" = c("LAMA2", "MKI67", "NRXN3")),
                     point_size = 0.2,
                     show_polygon = TRUE,
                     polygon_feat_type = "cell",
                     polygon_color = "white",
                     polygon_line_size = 0.1)

## Faster plotting method if you have many points- millions of points
#spatInSituPlotPoints(gobject = filt_resolve.gobj,
#                     plot_method = "scattermore",
#                     feats = list("rna" = c("LAMA2", "MKI67", "NRXN3")),
#                     point_size = 0.2,
#                     show_polygon = TRUE,
#                     polygon_feat_type = "cell",
#                     polygon_color = "white",
#                     polygon_line_size = 0.1)

#################### PART 8.1: SUBSET BY LOCATION ####################
## Can be used to focus on specific spatial structures
## For example to zoom in on niche environments

subloc = subsetGiottoLocs(gobject = filt_resolve.gobj,
                          x_min = 2000, x_max = 4000,
                          y_min = 4000, y_max = 6000)

## Show subset of genes
spatInSituPlotPoints(gobject = subloc,
                     feats = list("rna" = c("LAMA2", "MKI67", "NRXN3")),
                     point_size = 0.6,
                     show_polygon = TRUE,
                     polygon_feat_type = "cell",
                     polygon_color = "white",
                     polygon_line_size = 0.1)

## Show subset of genes and color cells according to clusters
spatInSituPlotPoints(gobject = subloc,
                     feats = list("rna" = c("LAMA2", "NRXN3", "MKI67")),
                     point_size = 0.6,
                     show_polygon = TRUE,
                     polygon_feat_type = "cell",
                     polygon_color = "white",
                     polygon_line_size = 0.1,
                     polygon_fill = "leiden_clus",
                     polygon_fill_as_factor = TRUE)

## Show subset of genes and color cells according to total expression
## Use a faster and more efficient point plotting method = scattermore
spatInSituPlotPoints(gobject = subloc,
                     plot_method = "scattermore",
                     feats = list("rna" = c("LAMA2", "NRXN3", "MKI67")),
                     point_size = 0.6,
                     show_polygon = TRUE,
                     polygon_feat_type = "cell",
                     polygon_color = "white",
                     polygon_line_size = 0.1,
                     polygon_fill = "total_expr",
                     polygon_fill_as_factor = FALSE)

## Show cells and color them according to total expression
spatInSituPlotPoints(gobject = subloc,
                     show_polygon = TRUE,
                     polygon_feat_type = "cell",
                     polygon_color = "white",
                     polygon_line_size = 0.1,
                     polygon_fill = "total_expr",
                     polygon_fill_as_factor = FALSE)

## Show cells and color them according to total cluster information
spatInSituPlotPoints(gobject = subloc,
                     show_polygon = TRUE,
                     polygon_feat_type = "cell",
                     polygon_color = "white",
                     polygon_line_size = 0.1,
                     polygon_fill = "leiden_clus",
                     polygon_fill_as_factor = TRUE)

#################### PART 9: FIND INTERACTION CHANGED FEATURES ####################

## Find interaction changed Features
## In this case, features are genes whose expression difference is associated with a neighboring cell type

test = findInteractionChangedFeats(gobject = filt_resolve.gobj,
                                   cluster_column = "leiden_clus")

test$ICFscores[type_int == "hetero"]

spatInSituPlotPoints(gobject = filt_resolve.gobj,
                     feats = list("rna" = c("LAMA2", "MKI67", "NRXN3")),
                     point_size = 0.6,
                     show_polygon = TRUE,
                     polygon_feat_type = "cell",
                     polygon_color = "black",
                     polygon_line_size = 0.1,
                     polygon_fill = "leiden_clus",
                     polygon_fill_as_factor = TRUE)

