###### Function to compute betadiversity within regions ##########

library(ape)
library(picante)
library(geiger)
library(treespace)
library(maptools)
library(raster)
library(readr)
library(tmap)
library(sp)
library(rgdal)
library(tidyverse)
library(data.table)
library(sf)
library(lwgeom)

### Create function to aggregate probabilities to higher hierachical level (aggregate pixel, or go up on a phylogenetic tree)
aggreg_prob = function(x, na.rm) { 
  y <- 1-prod(1-x) 
  return(y) # Output
}

### Load phylogeny
actinote_tree <- "actinote_certa.tree"
phylogeny_object<- ape::read.tree(actinote_tree)

### Load the complete stack of species SDM outputs
tif <- list.files (pattern ="Actinote")
tif
sp_proba_stack <- raster::stack(tif)
is.na(sp_proba_stack) <- 0

#### Function that match species list in stack and phylogeny and clean them ####

match_stack_and_phylo <- function (proba_stack, phylo)
{
  # List species in the stack
  species_list <- names(proba_stack)
  
  # Remove species not in the phylogeny from the stack
  proba_stack_for_phylo <- subset(proba_stack, which(species_list %in% phylo$tip.label))
  
  # Print species remove from stack because they are not in the phylogeny
  not_in_phylo <- species_list[!(species_list %in% phylo$tip.label)]
  if (length(not_in_phylo) > 0)
  {
    cat(paste0("\n", length(not_in_phylo), " species removed from stack because they are absent from the phylogeny:\n\n"))
    print(not_in_phylo)
    cat("\n")
  }
  
  # Prune the phylogeny to keep only species in the stack
  pruned_tree <- ape::keep.tip(phylo, names(proba_stack_for_phylo))
  
  # Print species remove from phylogeny because they are not in the stack
  not_in_stack <- phylo$tip.label[!(phylo$tip.label %in% names(proba_stack_for_phylo))]
  if (length(not_in_stack) > 0)
  {
    cat(paste0("\n", length(not_in_stack), " species removed from the phylogeny because they are absent from the Raster Stack:\n\n"))
    print(not_in_stack)
    cat("\n")
  }
  
  # Return cleaned stack and phylogeny in a list
  cleaned_stack_and_phylo <- list(proba_stack_for_phylo, pruned_tree)
  return(cleaned_stack_and_phylo)
}

test <- match_stack_and_phylo(sp_proba_stack, phylogeny_object)

#### download and prepare regions shapefile  #####

Boreal_Brazilian_dominion <- readShapeSpatial("Boreal_Brazilian_dominion_2.shp")
Boreal_Brazilian_dominion$fid = NULL
Boreal_Brazilian_dominion$cat = NULL
Mesoamerican_dominion <- readShapeSpatial("Mesoamerican_dominion_2.shp")
Mesoamerican_dominion$fid = NULL
Mesoamerican_dominion$cat = NULL
Parana_dominion <- readShapeSpatial("Parana_dominion_2.shp")
South_Brazilian_dominion <- readShapeSpatial("South_Brazilian_dominion_3.shp")
South_Brazilian_dominion$fid = NULL
South_Brazilian_dominion$cat = NULL
Chacoan_dominion <- readShapeSpatial("Chacoan_dominion_2.shp")
Chacoan_dominion$fid = NULL
Chacoan_dominion$cat = NULL
Pacific_dominion <- readShapeSpatial("Pacific_dominion_2.shp")
Pacific_dominion$fid = NULL
Pacific_dominion$cat = NULL
South_eastern_Amazonian_dominion <- readShapeSpatial("South_eastern_Amazonian_dominion_2.shp")
South_andes <- readShapeSpatial("South_andes_4.shp")
South_andes$fid = NULL
South_andes$cat = NULL
Central_America <- readShapeSpatial("Central_America_3.shp")
Central_America$fid = NULL
Central_America$cat = NULL
North_Andes <- readShapeSpatial("North_Andes_2.shp")
plot(Central_America)

list_sp_regions<- list(Boreal_Brazilian_dominion, Mesoamerican_dominion, Parana_dominion,
                       South_Brazilian_dominion, Chacoan_dominion, Pacific_dominion,
                       South_eastern_Amazonian_dominion, South_andes, Central_America, North_Andes)

region_names <- c("Boreal_Brazilian_dominion", "Mesoamerican_dominion", "Parana_dominion",
                  "South_Brazilian_dominion", "Chacoan_dominion", "Pacific_dominion",
                  "South_eastern_Amazonian_dominion", "South_andes", "Central_America", "North_Andes" )

### 12.1/ Betadiversity within regions ####

### 12.1.1/ Function to compute any type of betadiversity from a sites x species matrix ####

compute_betadiversity <- function (regional_data, diversity_type, # "taxo" or "phylo"
                                   phylo, index_family, # "jaccard" or "sorensen"
                                   beta_part, # "total", "nestedness" or "turnover"
                                   aggreg_type, # "multi" or "mean_pairwise"
                                   quiet = F) # To display progess
{
  if (!quiet) { cat(paste0(Sys.time(), " - Betadiversity Computation - Start for ", diversity_type," betadiversity ; Family = ",index_family," ; Partition = ",beta_part," ; Aggregation = ",aggreg_type,"\n")) }
  
  ### Prepare data 
  
  # Initiate check for the presence of data
  no_regional_data <- F
  
  # Clean sites with NA
  regional_data_clean <- na.omit(regional_data) 
  
  # Check there are at least two communities with data in the region
  if (nrow(regional_data_clean) < 2)
  {
    no_regional_data <- T
    beta_value <- NA
  }
  
  ### Compute indices if needed 
  if (!no_regional_data)
  {
    # For taxonomic beta-diversity based on species identity
    if (diversity_type == "taxo")
    {
      if (aggreg_type == "multi")
      {
        # Compute multi-sites taxonomic beta-diversity for all sites
        beta_results <- betapart::beta.multi(x = regional_data_clean, index.family = index_family)
      } else { 
        # Compute pairwise taxonomic beta-diversity for all pairs of sites
        beta_results <- betapart::beta.pair(x = regional_data_clean, index.family = index_family)
      }
    }
    
    # For phylogenetic beta-diversity based on Faith's Phylogenetic Diversity
    if (diversity_type == "phylo")
    {
      if (aggreg_type == "multi")
      {
        # Compute multi-sites phylogenetic beta-diversity for all sites
        beta_results <- betapart::phylo.beta.multi(x = regional_data_clean, tree = phylo, index.family = index_family)
      } else { 
        # Compute pairwise phylogenetic beta-diversity for all pairs of sites
        beta_results <- betapart::phylo.beta.pair(x = regional_data_clean, tree = phylo, index.family = index_family)
      }
    }
    
    ### Extract results
    
    # Depending on the partition of beta-diversity
    if (beta_part == "turnover") {beta_dist <- beta_results[[1]]}
    if (beta_part == "nestedness") {beta_dist <- beta_results[[2]]}
    if (beta_part == "total") {beta_dist <- beta_results[[3]]}
    
    # For multi-site index, pixel value is the multi-site beta-diversity
    if (aggreg_type == "multi")
    {
      beta_value <- round(as.numeric(beta_dist), 5)
    }
    
    # For pairwise index, pixel value is the mean pairwise beta-diversity across all pairs of regional sites
    if (aggreg_type == "mean_pairwise")
    {
      beta_value <- round(mean(beta_dist, na.rm = T), 5)
    }
  }  
  
  if (!quiet) { cat(paste0(Sys.time(), " - Betadiversity Computation - Done for ", diversity_type," betadiversity ; Family = ",index_family," ; Partition = ",beta_part," ; Aggregation = ",aggreg_type,"\n")) }
  
  return(beta_value)
}

#source(file = "./functions/compute_betadiversity.R")

### 12.1.2/ Function to compute betadiversity within regions ####

# Compute betadiversity across all pixels/sites of each region, and eventually map them

compute_within_regional_betadiversity <- function (proba_stack, list_sp_regions, region_names,
                                                   subsample_size = NA, # Provide number of pixels to subsample globally and regularly in order to save computing time, especially if phylobetadiversity is computed
                                                   diversity_type = "taxo",
                                                   phylo = NULL, index_family = "sorensen",
                                                   beta_part = c("turnover", "total"), 
                                                   aggreg_type = "mean_pairwise",
                                                   regions_map = T, # To create and display a map of regions with "within" regional betadiversity values
                                                   output_type = "region_df", # To choose the type of df as output. "region_df" fo a region x indices df. "ggplot_df" for a row per value. "Or"raster_only" to get just the raster Stack with all indices.
                                                   quiet = F) # To display progress at each iteration of the purrr:pmap()
{
  if (!quiet) { cat(paste0(Sys.time(), " - Within regional Betadiversity Computation - Start\n")) }
  
  ### Prepare input data
  final_output <- list()
  
  # Subsample sites in a standardized fashion such as density of sampling remains the same for all regions
  if (!is.na(subsample_size))
  {
    proba_stack <- raster::sampleRegular(x = proba_stack, size = subsample_size, asRaster = T)
  }
  
  # Binarize output following probability ranks if needed
  if (any(!is.element(unique(as.vector(getValues(proba_stack))), c(NA, 0, 1))))
  {
    binary_stack <- binarize_output_with_ranks(proba_stack)
  } else {
    binary_stack <- proba_stack
  }
  
  # If using the phylogenetic beta-diversity, need to match Raster Stack and phylogeny species lists
  pruned_tree <- phylo
  if ("phylo" %in% diversity_type)
  {
    # Match raster Stack and Phylogeny species lists
    clean_stack_and_phylo <- match_stack_and_phylo(binary_stack, phylo)
    binary_stack <- clean_stack_and_phylo[[1]]
    pruned_tree <- clean_stack_and_phylo[[2]]
  }
  
  # Aggregate SpatialPolygons for regions
  all_sp_polygons <- do.call(raster::bind, list_sp_regions)   
  #plot(all_sp_polygons)
  
  # Check CRS matching
  if (!raster::compareCRS(all_sp_polygons, proba_stack))
  {
    all_sp_polygons <- spTransform(x = all_sp_polygons, CRSobj = proba_stack@crs)
  }
  
  # Extract community matrix per regions
  binary_mat_per_regions <- raster::extract(binary_stack, all_sp_polygons)
  
  ### Generate the map of lists for pmap function crossing all argument combination.
  map_cross <- purrr::cross(list(regional_data = binary_mat_per_regions, diversity_type = diversity_type, index_family = index_family, beta_part = beta_part, aggreg_type = aggreg_type))
  
  # Function to revert list structure
  revert_list_str <- function(ls) 
  { 
    # Get sub-elements in same order
    x <- lapply(ls, `[`, names(ls[[1]]))
    # Stack and reslice
    apply(do.call(rbind, x), 2, as.list) 
  }
  
  # Revert list structure
  map_cross_revert <- revert_list_str(map_cross)
  # str(map_cross_revert)
  
  ### Compute betadiversity values for all regions, all diversity types, all index families, all partitions, and all aggregation types
  region_values <- purrr::pmap(.l = map_cross_revert, .f = compute_betadiversity, phylo = pruned_tree, quiet = quiet)
  
  ### Format output
  
  # Build final ggplot df
  within_regions_ggplot_df <- data.frame(regions = rep(x = region_names, length.out = length(map_cross)), diversity_type = unlist(map_cross_revert$diversity_type), index_family = unlist(map_cross_revert$index_family), beta_part = unlist(map_cross_revert$beta_part), aggreg_type = unlist(map_cross_revert$aggreg_type), index_value = unlist(region_values))
  
  # Convert into Regions x Indices matrix
  within_regions_df <- within_regions_ggplot_df %>% 
    tidyr::pivot_wider(data = ., names_from = c(diversity_type, index_family, beta_part, aggreg_type), names_sep = "_", values_from = index_value)
  
  # Print output in the requested format
  if (output_type == "region_df")
  {
    cat(paste0("Indices within regions\n"))
    print(within_regions_df)
  }
  
  if (output_type == "ggplot_df")
  {
    cat(paste0("Indices within regions\n"))
    print(within_regions_ggplot_df)
  }
  
  ### Produce map(s) of regions with within regions betadiversity values if requested
  
  if(regions_map)
  {
    # Generate a mask for terrestrial areas to use as background
    continental_mask <- (calc(proba_stack, fun = sum) >= 0) - 1
    
    # Loop per index
    all_indices_stack <- stack()
    for (i in 2:ncol(within_regions_df))
    {
      # i <- 4
      
      # Loop per region
      index_stack <- stack()
      for (j in 1:nrow(within_regions_df))
      {
        # j <- 1
        
        index_stack <- addLayer(index_stack, rasterize(x = list_sp_regions[[j]], #### mudei aqui!
                                                       y = continental_mask, # Provide the grid to fill with CRS, bbox and resolution
                                                       field = as.numeric(within_regions_df[j,i]), # How to fill non empty cells. With the value of a variable in the df of the sp_obj, or directly with a fixed value ?
                                                       background = NA)) # Value to use to fill empty cells)
      }
      # plot(index_stack) 
      
      # Aggregate all regions in one layer
      final_index <- calc(x = index_stack, fun = median, na.rm = T)
      # Add null values for continental borders
      final_index_raster <- continental_mask
      final_index_raster@data@values[!is.na(final_index[])] <- final_index[!is.na(final_index[])]
      
      # plot(final_index_raster)
      
      # Add the index raster Layer to the final stack
      all_indices_stack <- addLayer(all_indices_stack, final_index_raster)
    }
    
    # Add index name to each layer
    # index_names <- stringr::str_split(string = names(within_regions_df)[-1], pattern = "__")
    # index_names <- lapply(X = index_names, FUN = function (x) {paste0("Diversity type: ", x[1], " ; Family: ", x[2], "\nPartition: ", x[3], " ; Aggregation: ", x[4])})
    names(all_indices_stack) <- names(within_regions_df)[-1]
    
    # Plot Raster Stack with all indices
    plot(all_indices_stack)
    
    # Add plot to the list output (return)
    final_output <- list(all_indices_stack)
    
  }
  
  # Export output
  if (output_type == "region_df")
  {
    final_output <- append(list(within_regions_df), final_output)
    return(final_output)
  }
  
  if (output_type == "ggplot_df")
  {
    final_output <- append(list(within_regions_ggplot_df), final_output)
    return(final_output)
  }
  
  if (output_type == "raster_only")
  {
    return(unlist(final_output))
  }
  
  if (!quiet) { cat(paste0(Sys.time(), " - Within regional Betadiversity Computation - Done\n")) }
}

#source(file = "./functions/compute_betadiversity.R")

### 12.1.3/ Compute betadiversity within regions ####

# Compute the indices and generate maps
within_regional_betadiversity <- compute_within_regional_betadiversity(sp_proba_stack, list_sp_regions, region_names,
                                                                       subsample_size = 10000,
                                                                       diversity_type = c("taxo", "phylo"),
                                                                       #diversity_type = c("taxo"),
                                                                       phylo = phylogeny_object, index_family = "sorensen",
                                                                       beta_part = c("turnover"),
                                                                       aggreg_type = c("mean_pairwise"),
                                                                       regions_map = T, ############## NÃO consegui plotar no mapa!!!!!!!!!!!!!!!!!!!!!!!
                                                                       #output_type = "ggplot_df",
                                                                       output_type = "region_df",
                                                                       quiet = F)

View(within_regional_betadiversity[[1]])
plot(within_regional_betadiversity[[2]])


# Save
save(within_regional_betadiversity, file = paste0("./outputs/Indices_maps/within_regional_betadiversity.RData"))

# Compute differences per regions
load(file = paste0("./outputs/Indices_maps/within_regional_betadiversity.RData"))

diff_PBD_TBD_within_regions <- within_regional_betadiversity[[2]][[2]] - within_regional_betadiversity[[2]][[1]]
plot(diff_PBD_TBD_within_regions, zlim = c(-0.3, 0.3), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
