rm(list = ls())

#################### START GIOTTO ####################

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
#  installGiottoEnvironment(
#    packages_to_install = c("pandas", "networkx", "python-igraph",
#                            "leidenalg", "python-louvain", "python.app",
#                            "scikit-learn"),
#    force_miniconda = FALSE,
#    force_environment = TRUE,
#    verbose = TRUE)
#}


library(magick)
library(data.table)
library(RImageJROI)
library(future)
library(tidyr)
library(dplyr)
library(ggdendro)

# 1. set working directory
my_working_dir = '~/Giotto_loader_Resolve/'

# Optional: Specify a path to a Python executable within a conda or miniconda
# environment. If set to NULL (default), the Python executable within the previously
# installed Giotto environment will be used.
my_python_path = './.local/share/r-miniconda/envs/giotto_env/bin/python3'

#################### INPUT FILES ####################

## provide path to resolve bioscience folder
data_path = '~/Giotto_loader_Resolve/Cellpose.dragonfly.W6A1/'

# Create Giotto Instructions
instrs = createGiottoInstructions(save_dir = data_path,
                                  save_plot = TRUE,
                                  show_plot = FALSE,
                                  python_path = my_python_path)

# 1. original image as png
original_DAPI_image = '~/Giotto_loader_Resolve/Cellpose.dragonfly.W6A1/MAX_Resolve_MBENrun1_pseudoHE_C2-1_2022-07-19_13.48.59_Fused_405nm_corr_cut_registered.tif'
original_DAPI_image_magick = image_read('~/Giotto_loader_Resolve/Cellpose.dragonfly.W6A1/MAX_Resolve_MBENrun1_pseudoHE_C2-1_2022-07-19_13.48.59_Fused_405nm_corr_cut_registered.tif') %>% image_flop()

# 2. input cell segmentation as mask file
# can also be provided as a 3-column polygon file
# to be used as image background AND to store segmentations as polygons
# can be obtained through Fiji / QuPath / Ilastik / Cellpose / ...
segmentation_mask = '~/Giotto_loader_Resolve/Cellpose.dragonfly.W6A1/MAX_Resolve_MBENrun1_pseudoHE_C2-1_2022-07-19_13.48.59_Fused_405nm_corr_cut_registered_cp_masks.png'
segmentation_mask_magick = image_read('~/Giotto_loader_Resolve/Cellpose.dragonfly.W6A1/MAX_Resolve_MBENrun1_pseudoHE_C2-1_2022-07-19_13.48.59_Fused_405nm_corr_cut_registered_cp_masks.png') %>% image_flop()

# 3. input features coordinates
tx_coord = fread(paste0('~/Giotto_loader_Resolve/Cellpose.dragonfly.W6A1/Panorama_MBEN-slide1_W6A1_results_withFP.txt'), fill = TRUE)
colnames(tx_coord) = c('x', 'y', 'number_of_transcripts', 'gene_id', 'transcription_activity', 'blank')
tx_coord = tx_coord[,.(x, y, gene_id)]

#################### PART 1: CREATE SUBCELLULAR GIOTTO OBJECTS ####################

testobj = createGiottoObjectSubcellular(gpoints = list('rna' = tx_coord),
                                        gpolygons = list('cell' = segmentation_mask),
                                        #polygon_mask_list_params = list(flip_horizontal = FALSE),
                                        instructions = instrs,
                                        verbose = FALSE,
                                        cores = 32)

#################### PART 2: CREATE SPATIAL LOCATIONS ####################

# centroids are now used to provide the spatial locations (centroid of each cell)
# needed for certain downstream spatial analyses
testobj = addSpatialCentroidLocations(testobj,
                                      poly_info = 'cell')

#################### PART 3: ADD IMAGE INFORMATION ####################

# create Giotto images
DAPI_image = createGiottoImage(gobject = testobj,
                               name = 'DAPI',
                               do_manual_adj = TRUE,
                               xmax_adj = 0,ymax_adj = 0,
                               xmin_adj = 0,ymin_adj = 0,
                               mg_object = original_DAPI_image_magick)
print(DAPI_image)

segm_image = createGiottoImage(gobject = testobj,
                               name = 'segmentation',
                               do_manual_adj = TRUE,
                               xmax_adj = 0,ymax_adj = 0,
                               xmin_adj = 0,ymin_adj = 0,
                               mg_object = segmentation_mask_magick)
print(segm_image)

# add images to Giotto object
testobj = addGiottoImage(testobj, images = list(segm_image, DAPI_image))

# provides an overview of available images
showGiottoImageNames(testobj)

#################### PART 4: VISUALIZE ORIGINAL IMAGES ####################

# visualize overlay of calculated cell centroid with original image and segmentation mask file
# by setting show_plot to FALSE and save_plot to TRUE you can save quite some time when creating plots
# with big images it sometimes takes quite long for R/Rstudio to render them
DAPI_spatplot <- spatPlot2D(gobject = testobj,
                            image_name = 'DAPI',
                            point_size = 1.0,
                            show_plot = TRUE)
print(DAPI_spatplot)

segm_spatplot <- spatPlot2D(gobject = testobj,
                            image_name = 'segmentation',
                            point_size = 1.0,
                            show_plot = TRUE)
print(segm_spatplot)

#################### PART 5: CALCULATE CELL SHAPE OVERLAP ####################

tictoc::tic()
testobj = calculateOverlapParallel(testobj,
                           spatial_info = 'cell',
                           feat_info = 'rna')
tictoc::toc()

#convert overlap to matrix
testobj = overlapToMatrix(testobj,
                          poly_info = 'cell',
                          feat_info = 'rna',
                          name = 'raw')

#################### PART 6: FILTER DATA ####################

# features can be filtered individually
# cells will be filtered across features

# first filter on rna: 
# expression_threshold: threshold to consider a gene expressed
# feat_det_in_min_cells: minimum # of cells that need to express a feature
# min_det_feats_per_celL: minimum # of features that need to be detected in a cell

subc_test <- filterGiotto(gobject = testobj,
                          expression_threshold = 1,
                          feat_det_in_min_cells = 20, 
                          min_det_feats_per_cell = 5)

subc_spatplot <- spatPlot2D(gobject = subc_test, 
                            image_name = 'segmentation', 
                            show_image = TRUE, 
                            point_size = 1.0)

print(subc_spatplot)

# rna data, default.
# other feature modalities can be processed and filtered in an anologous manner
subc_test <- normalizeGiotto(gobject = subc_test, 
                             scalefactor = 6000, 
                             verbose = T)

subc_test <- addStatistics(gobject = subc_test)

subc_test <- adjustGiottoMatrix(gobject = subc_test,
                                expression_values = c('normalized'),
                                covariate_columns = c('nr_feats', 'total_expr'))

subc_test <- normalizeGiotto(gobject = subc_test, 
                             norm_methods = 'pearson_resid', 
                             update_slot = 'pearson')

showGiottoExpression(subc_test)

#################### PART 8: DIMENSION REDUCTION ####################

# Find highly valuable Features

# typical way of calculating HVF
subc_test <- calculateHVF(gobject = subc_test, HVFname= 'hvg_orig')

# new method based on variance of pearson residuals for each gene
subc_test <- calculateHVF(gobject = subc_test,
                          method = 'var_p_resid', 
                          expression_values = 'pearson',
                          show_plot = TRUE)

#run PCA
subc_test <- runPCA(gobject = subc_test,
                    expression_values = 'pearson',
                    scale_unit = FALSE,
                    center = FALSE)

screePlot(subc_test, ncp = 20)

plotPCA(subc_test,
        dim1_to_use = 1,
        dim2_to_use = 2)

# run UMAP
subc_test <- runUMAP(subc_test, 
                     dimensions_to_use = 1:5, 
                     n_threads = 2)

plotUMAP(gobject = subc_test)

#################### PART 9: CLUSTER ####################

subc_test <- createNearestNetwork(gobject = subc_test, 
                                  dimensions_to_use = 1:5, 
                                  k = 5)

subc_test <- doLeidenCluster(gobject = subc_test, 
                             resolution = 0.05, 
                             n_iterations = 1000, 
                             name = 'leiden_clus')

# Create color palettes, or proceed with Giotto defaults
#library(RSkittleBrewer)
#colorcode = RSkittleBrewer(flavor = c("smarties"))
#featcolor = RSkittleBrewer(flavor = c("M&M"))

# visualize UMAP cluster results
plotUMAP(gobject = subc_test, 
         cell_color = 'leiden_clus',
         show_NN_network = TRUE, 
         point_size = 2.5)

# visualize UMAP and spatial results
spatDimPlot2D(gobject = subc_test,
              show_image = TRUE, 
              image_name = 'segmentation',
              cell_color = 'leiden_clus',
              spat_point_size = 2.5)

# Plot a cluster heatmap
showClusterHeatmap(gobject = subc_test, cluster_column = 'leiden_clus',
                   save_param = list(save_format = 'pdf',
                                     base_height = 6, 
                                     base_width = 8, 
                                     units = 'cm'))

# See cluster relationships in a dendogram
showClusterDendrogram(subc_test, 
                      h = 0.5, 
                      rotate = TRUE, 
                      cluster_column = 'leiden_clus')

#################### PART 10: CREATE A SPATIAL NETWORK ####################

subc_test = createSpatialNetwork(gobject = subc_test,
                                 spat_loc_name = 'raw',
                                 minimum_k = 3,
                                 maximum_distance_delaunay = 100)

spatPlot2D(gobject = subc_test,
           image_name = 'segmentation', 
           show_image = TRUE,
           point_size = 1.5, 
           show_network = TRUE)

#################### PART 11: VISUALIZE SUBCELLULAR DATA ####################

# Visualize clustered cells
spatInSituPlotPoints(subc_test,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'white',
                     polygon_line_size = 0.1,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = TRUE)

# individual plotting of transcripts and polygon information

# all cells
spatInSituPlotPoints(testobj,
                     #feats = list('rna' = c("MMP2", "VEGFA", "IGF1R", 'CDH2', 'MKI67')),
                     feats = list('rna' = c("LAMA2", "MKI67", "NRXN3")),
                     point_size = 0.2,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'white',
                     polygon_line_size = 0.1)

# filtered cells
tictoc::tic()
spatInSituPlotPoints(subc_test,
                     #feats = list('rna' = c("MMP2", "VEGFA", "IGF1R", 'CDH2', 'MKI67')),
                     feats = list('rna' = c("LAMA2", "MKI67", "NRXN3")),
                     point_size = 0.2,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'white',
                     polygon_line_size = 0.1)
tictoc::toc()
# 17.019 sec elapsed

# faster plotting method if you have many points- millions of points
#tictoc::tic()
#spatInSituPlotPoints(subc_test,
#                     plot_method = 'scattermore',
#                     #feats = list('rna' = c("MMP2", "VEGFA", "IGF1R", 'CDH2', 'MKI67')),
#                     feats = list('rna' = c("LAMA2", "MKI67", "NRXN3")),
#                     point_size = 0.2,
#                     show_polygon = TRUE,
#                     polygon_feat_type = 'cell',
#                     polygon_color = 'white',
#                     polygon_line_size = 0.1)
#tictoc::toc()
# 17.311 sec elapsed

#################### PART 11.1: SUBSET BY LOCATION ####################
# can be used to focus on specific spatial structures
# to zoom in on niche environments

subloc = subsetGiottoLocs(subc_test,
                          x_min = 2000, 
                          x_max = 4000,
                          y_min = 4000, 
                          y_max = 6000)

# show subset of genes
spatInSituPlotPoints(subloc,
                     #feats = list('rna' = c("MMP2", "VEGFA", "IGF1R", 'CDH2', 'MKI67')),
                     feats = list('rna' = c("LAMA2", "MKI67", "NRXN3")),
                     point_size = 0.6,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'white',
                     polygon_line_size = 0.1)

# show subset of genes and color cells according to clusters
spatInSituPlotPoints(subloc,
                     #feats = list('rna' = c("MMP2", "VEGFA", "IGF1R", 'CDH2', 'MKI67')),
                     feats = list('rna' = c("LAMA2", "NRXN3", "MKI67")),
                     point_size = 0.6,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'white',
                     polygon_line_size = 0.1,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = TRUE)

# show subset of genes and color cells according to total expression
# use a faster and more efficient point plotting method = scattermore
spatInSituPlotPoints(subloc,
                     plot_method = 'scattermore',
                     #feats = list('rna' = c("MMP2", "VEGFA", "IGF1R", 'CDH2', 'MKI67')),
                     feats = list('rna' = c("LAMA2", "NRXN3", "MKI67")),
                     point_size = 0.6,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'white',
                     polygon_line_size = 0.1,
                     polygon_fill = 'total_expr',
                     polygon_fill_as_factor = FALSE)

# show cells and color them according to total expression
spatInSituPlotPoints(subloc,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'white',
                     polygon_line_size = 0.1,
                     polygon_fill = 'total_expr',
                     polygon_fill_as_factor = FALSE)

# show cells and color them according to total cluster information
spatInSituPlotPoints(subloc,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'white',
                     polygon_line_size = 0.1,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = TRUE)

#################### PART 12: FIND INTERACTION CHANGED FEATURES ####################
# find interaction changed Features
# In this case, features are genes whose expression difference is associated with a neighboring cell type
future::plan('multisession', workers = 4) # sometimes unstable, restart R session

test = findInteractionChangedFeats(gobject = subc_test,
                                   cluster_column = 'leiden_clus')

test$ICFscores[type_int == 'hetero']

spatInSituPlotPoints(subc_test,
                     feats = list('rna' = c("LAMA2", "MKI67", "NRXN3")),
                     point_size = 0.6,
                     show_polygon = TRUE,
                     polygon_feat_type = 'cell',
                     polygon_color = 'black',
                     polygon_line_size = 0.1,
                     polygon_fill = 'leiden_clus',
                     polygon_fill_as_factor = TRUE)

