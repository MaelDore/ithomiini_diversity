###


### Clean environment
rm(list = ls())

### Load libraries

library(raster)
library(tidyverse)

### Load Patricia shp files

# Chacoan_dominion_2 <- rgdal::readOGR(dsn = "./input_data/Map_stuff/Patricia_data", layer = "Chacoan_dominion_2") # Pour charger depuis le dossier de travail
# South_eastern_Amazonian_dominion_2 <- rgdal::readOGR(dsn = "./input_data/Map_stuff/Patricia_data", layer = "South_eastern_Amazonian_dominion_2") # Pour charger depuis le dossier de travail
# 
# plot(Chacoan_dominion_2)
# plot(South_eastern_Amazonian_dominion_2)
# 
# 
# Lowenberg_Neto_2014 <- rgdal::readOGR(dsn = "./input_data/Map_stuff/Patricia_data", layer = "Lowenberg_Neto_2014") # Pour charger depuis le dossier de travail
# 
# save(Lowenberg_Neto_2014, file = "./input_data/Map_stuff/Patricia_data/Lowenberg_Neto_2014.RData")
# 
# plot(test, col = factor(test@data$Dominions))
# plot(test, col = factor(test@data$Province_1))
# plot(test, col = factor(test@data$Subregio_1))

### Load recursively all shp files in a list

# Get shp names
path_to_shp_files <- str_replace(string = list.files(path = "./input_data/Map_stuff/Patricia_data/", pattern = ".shp$", full.names = F), pattern = ".shp", replacement = "")
# Filter out the one you do not want
Patricia_region_names <- path_to_shp_files <- path_to_shp_files[!str_detect(string = path_to_shp_files, pattern = "Morrone|Lowenberg_Neto")]
# Import all in a list
Patricia_list_shp_regions <- lapply(X = path_to_shp_files, FUN = rgdal::readOGR, dsn = "./input_data/Map_stuff/Patricia_data")

# Merge all
Patricia_all_shp_polygons <- do.call(raster::bind, Patricia_list_shp_regions)

str(Patricia_list_shp_regions[[2]])

color_palette_pastel <- RColorBrewer::brewer.pal(n = length(Patricia_list_shp_regions), name = "Set3")
plot(Patricia_all_shp_polygons, col = color_palette_pastel)
# See the issue with multiple polygon regions which take multiple values


### Load my own data for a test

load(file = paste0("./outputs/Indices_maps/sp_binary_stack.RData"))

Ithomiini_phylo <- readRDS(file = "./input_data/Phylogenies/Final_phylogeny.rds")

### Try on my own data

source(file = "./functions/compute_betadiversity.R")

# list_shp_regions <- list(Chacoan_dominion_2, South_eastern_Amazonian_dominion_2)
# region_names <- c("Chacoan_dominion_2", "South_eastern_Amazonian_dominion_2")

within_regional_betadiversity <- compute_within_regional_betadiversity(proba_stack = sp_binary_stack, 
                                                                       list_shp_regions = Patricia_list_shp_regions,
                                                                       region_names = Patricia_region_names,
                                                                       subsample_size = 10000,
                                                                       diversity_type = c("taxo", "phylo"),
                                                                       # diversity_type = c("taxo"),
                                                                       phylo = Ithomiini_phylo,
                                                                       index_family = "sorensen",
                                                                       beta_part = c("turnover"),
                                                                       aggreg_type = c("mean_pairwise"),
                                                                       regions_map = T,
                                                                       output_type = "region_df",
                                                                       quiet = F)

View(within_regional_betadiversity[[1]])
plot(within_regional_betadiversity[[2]], zlim = c(0, 0.6))



### Create Patricia raster's stack

path_to_tif_files <- list.files(path = "./input_data/Map_stuff/Patricia_data/", pattern = ".tif$", full.names = T)
Patricia_sp_binary_stack <- stack(path_to_tif_files)

# Assign CRS EPSG:4326
crs(Patricia_sp_binary_stack) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"  

names(Patricia_sp_binary_stack)

plot(Patricia_sp_binary_stack)

# Check if binary
table(getValues(Patricia_sp_binary_stack))

### Load Patrica's phylogeny

Actinote_phylo <- ape::read.tree(file = "./input_data/Map_stuff/Patricia_data/actinote_certa.tree")
plot(Actinote_phylo)

### Try with Patricia's data

source(file = "./functions/compute_betadiversity.R")

# list_shp_regions <- list(Chacoan_dominion_2, South_eastern_Amazonian_dominion_2)
# region_names <- c("Chacoan_dominion_2", "South_eastern_Amazonian_dominion_2")



### Within regional betadiversity

within_regional_betadiversity <- compute_within_regional_betadiversity(proba_stack = Patricia_sp_binary_stack,
                                                                       list_shp_regions = Patricia_list_shp_regions,
                                                                       region_names = Patricia_region_names,
                                                                       subsample_size = 10000,
                                                                       diversity_type = c("taxo", "phylo"),
                                                                       # diversity_type = c("taxo"),
                                                                       phylo = Actinote_phylo,
                                                                       index_family = "sorensen",
                                                                       beta_part = c("turnover"),
                                                                       aggreg_type = c("mean_pairwise"),
                                                                       regions_map = T,
                                                                       output_type = "region_df",
                                                                       quiet = F)

View(within_regional_betadiversity[[1]])
plot(within_regional_betadiversity[[2]])


### Compute pairwise betadiversity between regions ####

pairwise_regional_betadiversity <- compute_pairwise_regional_betadiversity(proba_stack = Patricia_sp_binary_stack,
                                                                           list_shp_regions, region_names,
                                                                           diversity_type = c("taxo", "phylo"),
                                                                           # diversity_type = c("taxo"),
                                                                           phylo = Actinote_phylo, index_family = "sorensen",
                                                                           beta_part = c("turnover"))

View(pairwise_regional_betadiversity[[1]])

### Map both within and between betadiversity
map_betadiversity_within_and_between_regions (region_df_within = within_regional_betadiversity[[1]],   # To provide directly the df with indices values within each region. 2 columns: the region names, the index values
                                              region_matrix_between = pairwise_regional_betadiversity[[1]],  # To provide directly the matrix of betadiversity values between regions
                                              list_shp_regions, # To provide the list of SpatialPolygons of regions to extract centroid coordinates
                                              region_labels = region_names, # To provide labels to use instead of full region names
                                              nod_cex = 80,
                                              edge_cex = 2.5,
                                              plot_type = "ggplot")

map_betadiversity_within_and_between_regions (region_df_within = within_regional_betadiversity[[1]],   # To provide directly the df with indices values within each region. 2 columns: the region names, the index values
                                              region_matrix_between = pairwise_regional_betadiversity[[1]],  # To provide directly the matrix of betadiversity values between regions
                                              list_shp_regions, # To provide the list of SpatialPolygons of regions to extract centroid coordinates
                                              region_labels = region_names, # To provide labels to use instead of full region names
                                              nod_cex = 120,
                                              edge_cex = 25,
                                              plot_type = "igraph")

### Map regional residuals of PBD ~ TBD

source(file = "./functions/compute_residuals.R")

PBD_TBD_regions_residuals <- map_residuals_within_regions (y_index = within_regional_betadiversity[[2]][[2]],
                                                           x_index = within_regional_betadiversity[[2]][[1]], 
                                                           list_shp_regions, # Provide region Spatial Polygon to define borders
                                                           method = "GAM", # Either LM for linear relationship, LM_quadratic to add quadratic effects, or GLS to account for spatial autocorrelation, or GAM to allow non-linear relationship.
                                                           corSpatial = "corExp", # To choose the type of spatial autocorrelation structure as in nlme::corSpatial() for GLS
                                                           include_model = F, # To include the fitted model in the output
                                                           scatterplot = T,  # To display the scatter plot with regression line and R²
                                                           plot_variogram_correlogram = F) # To display the variogram and correlogram in case of spatial GLS


plot(PBD_TBD_regions_residuals, col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))


### Test for significant PBD relative to TBD (species composition is fixed, only phylogeny is randomized)

source(file = "./functions/compute_SES_tests.R")

SES_PBD_vs_TBD_within_regions <- compute_SES_PBD_vs_TBD_within_regions(proba_stack = Patricia_sp_binary_stack, 
                                                                       list_shp_regions, region_names,
                                                                       subsample_size = 1000, # Provide number of pixels to subsample globally and regularly in order to save computing time
                                                                       # subsample_size = NA, # Provide number of pixels to subsample globally and regularly in order to save computing time
                                                                       phylo = Actinote_phylo, # Provide the phylogeny
                                                                       index_family = "sorensen", # Chose the family of Beta-diversity indices
                                                                       beta_part = "turnover",    # Chose the partition
                                                                       aggreg_type = "mean_pairwise",     # Chose the aggregation method
                                                                       seed = 1, randomizations = 99, # To select the number of randomization
                                                                       alpha = 0.05,  # Threshold for significance for a two-sided test
                                                                       export_null_df = F, # To include the df with indices values from all randomization in the output
                                                                       export_region_stats_df = T, # To include the df with statistics from the tests for each region
                                                                       plot_distri = T, # To plot the null distribution for an "average" community as an example
                                                                       example_region_name = "South_Brazilian_dominion_3") # To provide the name of the region to use as example. Otherwise it is randomly chosen.


