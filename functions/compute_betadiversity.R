##### Functions to compute BetaDiversity indices from raster stack of SDM continuous outputs #####

# Author: Maël Doré
# Contact: mael.dore@gmail.com

### Contents ###

# 1/ Internal functions to compute betadiversity
# 2/ Betadiversity in reference to a focal point
# 3/ Betadiversity with a moving window
# 4/ Betadiversity within regions
# 5/ Betadiversity between regions
# 6/ Network graph for wihtin/between betadiversity of regions
# 7/ Pairwise betadiversity in RGB space

# Most functions allow to chose different types of betadiversity

  # Chose between "taxo" and "phylo" for the type of betadiversity indices with the "diversity_type" argument
  # Chose between "sorensen" and jaccard for the family of betadiversity indices with the "index_family" argument (See Baselga et al., 2012)
  # Chose between "total", "nestedness" and "turnover" for the type of partition of betadiversity with the "beta_part" argument (See Baselga et al., 2010, 2012)
  # Chose between "multi, mean_pairwise", "all_pairwise" and "focal_pairwise" for the way to compute a unique index from multiple sites with the "aggreg_type" argument 

library("betapart")

### 1/ Internal functions to compute betadiversity ####

### 1.1/ Function to binarize SDM continuous outputs following SDM scores ranking ####

# Input = Raster Stack of species/OMU probability of presence

# Output = Raster Stack of species/OMU binary presence/absence

# Do not apply a threshold to define presence/absence. 
# Use a ranking approach. For each community, the N species with the higher value from continuous SDM outputs are considered present.
# N is the expected local species richness such as the sum of species SDM continuous outputs.

### WARNING: make the assumption that SDM continuous outputs are comparable between species and can be ranked properly.
# Ideally, SDM continuous outputs should represent the probability of presence of each species.
# Strictly, it should only be applied in the case of presence/true absences models where the probability of detection of species
# is homogeneous across species, or the detection process is explicitely modelled in the SDM (See Guillera-Arroita et al. 2015)


binarize_output_with_ranks <- function (proba_stack)
{
  cat(paste0(Sys.time(), " - Binarization - Start\n"))
  
  # Extract community matrix
  proba_mat <- raster::getValues(proba_stack)
  
  # Extract community richness
  com_richness <- raster::getValues(calc(x = proba_stack, fun = sum))
  # Round it to get only integer as number of species expected
  rounded_richness <- round(x = com_richness, digits = 0)
  
  # Rank species probability of presence in each community
  com_rank_mat <- t(apply(X = proba_mat, MARGIN = 1, FUN = rank, na.last = "keep", ties.method = "random"))
  
  # Create 3D array to store binary data
  binary_array <- array(dim = c(proba_stack@nrows, proba_stack@ncols, ncol(proba_mat)),
                        dimnames = list(NULL, NULL, colnames(proba_mat)))
  
  # Fill the binary array, community by community
  for (i in 1:nrow(proba_mat))
  {
    # Extract species rank in the community
    com_rank <- com_rank_mat[i, ]
    # Binarize according to expected species richness
    binary_com <- com_rank > (max(com_rank) - rounded_richness[i])
    # Retrieve community indices in the 3D array
    RowCol <- raster::rowColFromCell(object = proba_stack, cell = i)
    # Fill data in the array
    binary_array[RowCol[1], RowCol[2], ] <- as.numeric(binary_com)
    
    if (i %% 1000 == 0) { cat(paste0(Sys.time(), " - ", i," on ",nrow(proba_mat),"\n")) }
  }
  
  # Write binary data in a Raster Brick
  binary_brick <- raster::brick(x = binary_array,
                                xmn = proba_stack@extent[1],
                                xmx = proba_stack@extent[2],
                                ymn = proba_stack@extent[3],
                                ymx = proba_stack@extent[4],
                                crs = proba_stack@crs)
  
  # Repair issue with min & max value
  binary_brick@data@min <- as.numeric(apply(X = binary_array, MARGIN = 3, FUN = min, na.rm = T))
  binary_brick@data@max <- as.numeric(apply(X = binary_array, MARGIN = 3, FUN = max, na.rm = T))
  
  # Convert to Raster Stack
  binary_stack <- raster::stack(binary_brick)
  
  return(binary_stack)
  
  cat(paste0(Sys.time(), " - Binarization - Done\n"))
}


### 1.2/ Function to compute any type of betadiversity from a sites x species matrix ####

# Input = sites x species matrix/df

# Output = Biodiversity index value


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

### 1.3/ Function to compute pairwise betadiversity matrix from a sites x species matrix ####

# Input = sites x species matrix/df

# Output = Matrix of pairwise betadiversity indices between sites


compute_pairwise_betadiversity <- function (regional_data, diversity_type, # "taxo" or "phylo"
                                            phylo, index_family, # "jaccard" or "sorensen"
                                            beta_part) # "total", "nestedness" or "turnover"
{
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
      # Compute pairwise taxonomic beta-diversity for all pairs of sites
      beta_results <- betapart::beta.pair(x = regional_data_clean, index.family = index_family)
    }
    
    # For phylogenetic beta-diversity based on Faith's Phylogenetic Diversity
    if (diversity_type == "phylo")
    {
      # Compute pairwise phylogenetic beta-diversity for all pairs of sites
      beta_results <- betapart::phylo.beta.pair(x = regional_data_clean, tree = phylo, index.family = index_family)
    }
    
    ### Extract results
    
    # Depending on the partition of beta-diversity
    if (beta_part == "turnover") {beta_mat <- as.matrix(beta_results[[1]])}
    if (beta_part == "nestedness") {beta_mat <- as.matrix(beta_results[[2]])}
    if (beta_part == "total") {beta_mat <- as.matrix(beta_results[[3]])}
    
  }  
  return(beta_mat)
}


### 1.4/ Function that matches species list in stack and phylogeny and clean them ####

# Inputs = Raster Stack with species names as layer names
#          Phylogeny with tips matching species names format

# Output = Raster Stack with only species found in the phylogeny
#          Phylogeny with only species found in the Raster Stack

# Prints = List of species removed from the Raster Stack and the phylogeny


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



### 2/ Function to compute betadiversity in reference to a focal point ####


# Inputs = Raster Stack of species/OMU probability of presence
#          Phylogeny of the clade (for PhyloBetaDiversity indices)
#          Coordinates of the focal site

# Output = Raster Layer of betadiversity in reference to the the focal site


### To improve
  # Allow to provide a sp polygon or sf object to define a focal region instead of a focal point
  # Aggregate info in a single line to use as the focal (duplicate it for each pixel so they appear with value = 0)
  # Add the option to compute for combination of parameters, all at the same time...

compute_focal_betadiversity <- function (proba_stack, diversity_type = "taxo",
                                         phylo = NULL, index_family = "sorensen",
                                         beta_part = "turnover", focal_site_coords)
{
  ### Prepare data
  
  # Binarize output following probability ranks if needed
  if (any(!is.element(unique(as.vector(getValues(proba_stack))), c(NA, 0, 1))))
  {
    binary_stack <- binarize_output_with_ranks(proba_stack)
  } else {
    binary_stack <- proba_stack
  }
  
  # If using the phylogenetic beta-diversity, need to match raster Stack and phylogeny species lists
  if (diversity_type == "phylo")
  {
    # Match raster Stack and Phylogeny species lists
    clean_stack_and_phylo <- match_stack_and_phylo(binary_stack, phylo)
    binary_stack <- clean_stack_and_phylo[[1]]
    pruned_tree <- clean_stack_and_phylo[[2]]
  }
  
  # Extract community matrix
  binary_mat <- raster::getValues(binary_stack)
  
  # Remove communities with NA
  binary_mat_clean <- na.omit(binary_mat)
  removed_com_indices <- as.numeric(attr(binary_mat_clean, which = "na.action"))
  
  ### Compute betadiversity with Baselga partitioning between turnover and nestedness
  
  # For taxonomical beta-diversity based on species identity
  if (diversity_type == "taxo")
  {
    beta_results <- betapart::beta.pair(x = binary_mat_clean, index.family = index_family)
  }
  
  # For phylogenetic beta-diversity based on Faith's Phylogenetic Diversity
  if (diversity_type == "phylo")
  {
    beta_results <- betapart::phylo.beta.pair(x = binary_mat_clean, tree = pruned_tree, index.family = index_family)
  }
  
  # Depending on the partition of beta-diversity
  if (beta_part == "turnover") {beta_mat <- as.matrix(beta_results[[1]])}
  if (beta_part == "nestedness") {beta_mat <- as.matrix(beta_results[[2]])}
  if (beta_part == "total") {beta_mat <- as.matrix(beta_results[[3]])}
  
  # Rebuilt matrix including NA communities
  beta_mat_full <- matrix(data = NA, nrow = nrow(binary_mat), ncol = nrow(binary_mat))
  if (length(removed_com_indices) > 0)
  {
    beta_mat_full[-removed_com_indices, -removed_com_indices] <- beta_mat
  } else {
    beta_mat_full <- beta_mat
  }
  
  ### Extract only result in reference to the focal point
  
  # Convert coordinates into the proper CRS if not in Latitude and Longitude
  if (!stringr::str_detect(string = proba_stack@crs, pattern = "\\+proj=longlat"))
  {
    focal_site_coords <- data.frame(lon = focal_site_coords[1], lat = focal_site_coords[2])
    coordinates(focal_site_coords) <- c("lon", "lat") # Transform into Spatial Points object
    proj4string(focal_site_coords) <- CRS("+init=epsg:4326") # Assign WGS84 CRS
    focal_site_coords <- sp::spTransform(x = focal_site_coords, CRSobj = proba_stack@crs)
    focal_site_coords <- focal_site_coords@coords
  }
  
  # Get the index of the focal site
  focal_site_index <- raster::cellFromXY(object = binary_stack, xy = focal_site_coords)
  # Extract associated beta-diversity indices
  beta_focal <- beta_mat_full[focal_site_index,]
  
  ### Format the output
  
  # Write focal betadiversity in a raster Layer
  focal_betadiv <- proba_stack[[1]]
  focal_betadiv@data@values <- beta_focal
  
  # Repair issue with min/max values
  focal_betadiv@data@min <- min(focal_betadiv[], na.rm = T)
  focal_betadiv@data@max <- max(focal_betadiv[], na.rm = T)
  
  return(focal_betadiv)
  
}


### 3/ Function to compute betadiversity with a moving window ####

### Compute betadiversity across multiple sites around a focal point within a moving window

# Need to specify the window size (width and height) in number of pixels from the focal point to the border of the window.
# By default, the value equals a fifth of the study area

# Can specify the way to compute a unique index from the multiple sites in the window using the "aggreg_type" argument
   # "all_pairwise" provides the mean value between all pairwise betadiversity indices across all pairs of sites within the window 
   # "focal_pairwise" provides the mean value between all pairwise betadiversity indices of all sites in reference to the focal site.   
   # "multi" provides the multi-sites betadiversity index as in Baselga et al., 2012


# Inputs = Raster Stack of species/OMU probability of presence
#          Phylogeny of the clade (for PhyloBetaDiversity indices)

# Output = Raster Layer of continuous betadiversity

### To improve
  # Deal with border rules
  # Compute more efficiently the "focal_pairwise" option without coputing the full pairwise matrix
  # Add the option to compute for combination of parameters, all at the same time...
  # Add the option to provide a template for the window (such as kernel for geomorpho)


compute_moving_betadiversity <- function (proba_stack, diversity_type = "taxo",
                                          phylo = NULL, index_family = "sorensen",
                                          beta_part = "turnover", aggreg_type = "all_pairwise",
                                          window_width = "default", window_height = "default", # In pixels, extended from the focal point
                                          border_rule = NULL)  # Add an option on how to deal with the borders
{
  cat(paste0(Sys.time(), " - Betadiversity Computation - Start\n"))
  
  ### Adjust default window size to a fifth of raster width and height
  
  if (window_width == "default")
  {
    window_width <- round(proba_stack@ncols/10)
    if (window_width == 0) {window_width <- 1}
  }
  if (window_height == "default")
  {
    window_height <- round(proba_stack@nrows/10)
    if (window_height == 0) {window_height <- 1}
  }
  
  ### Prepare data
  
  # Binarize output following probability ranks if needed
  if (any(!is.element(unique(as.vector(getValues(proba_stack))), c(NA, 0, 1))))
  {
    binary_stack <- binarize_output_with_ranks(proba_stack)
  } else {
    binary_stack <- proba_stack
  }
  
  # If using the phylogenetic beta-diversity, need to match raster Stack and phylogeny species lists
  if (diversity_type == "phylo")
  {
    # Match raster Stack and Phylogeny species lists
    clean_stack_and_phylo <- match_stack_and_phylo(binary_stack, phylo)
    binary_stack <- clean_stack_and_phylo[[1]]
    pruned_tree <- clean_stack_and_phylo[[2]]
  }
  
  # Extract community matrix
  binary_mat <- raster::getValues(binary_stack)
  
  ### Loop to move from pixel to pixel
  all_beta_values <- NA
  for (i in 1:nrow(binary_mat))
  {
    # i <- 1
    
    # Initiate control for focal data and regional data
    no_focal_data <- F
    no_regional_data <- F
    
    # Extract cell location
    RowCol <- raster::rowColFromCell(object = binary_stack, cell = i)
    
    ### Warning: Window size is reduced for right (add) and bottom borders (cutoff), but not for left and upper borders (shift window)
    # Add an option on how to deal with the borders with border_rule
    
    # Extract regional data
    regional_data <- raster::getValuesBlock(x = binary_stack, 
                                            row = RowCol[1] - window_width,
                                            col = RowCol[2] - window_height,
                                            nrow = window_width*2 + 1,
                                            ncol = window_height*2 + 1)
    ### Prepare assemblage data
    
    if (aggreg_type == "focal_pairwise") # For focal pairwise, need to put the focal site on the first row
    {
      # Extract data from the focal site
      ### Warning: does not work when the window reach the bottom border because extra rows are cut-off...
      focal_data <- regional_data[ceiling(nrow(regional_data)/2), ]
      # Check if focal site has data
      if (any(is.na(focal_data)))
      {
        no_focal_data <- T
        beta_value <- NA
      }
      
      # Remove focal data from the regional data
      regional_data_clean <- regional_data[-ceiling(nrow(regional_data)/2), ]
      
      # Clean sites with NA
      regional_data_clean <- na.omit(regional_data_clean) 
      
      # Add the focal site on the first row
      regional_data_clean <- rbind(focal_data, regional_data_clean)
      
    } else { # For other aggregation type, only need to clean sites with NA
      # Clean sites with NA
      regional_data_clean <- na.omit(regional_data) 
    }
    
    ### Check there are at least two communities with data in the region
    if (nrow(regional_data_clean) < 2)
    {
      no_regional_data <- T
      beta_value <- NA
    }
    
    ### Compute index if needed
    
    if (!no_focal_data & !no_regional_data)
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
          beta_results <- betapart::phylo.beta.multi(x = regional_data_clean, tree = pruned_tree, index.family = index_family)
        } else { 
          # Compute pairwise phylogenetic beta-diversity for all pairs of sites
          beta_results <- betapart::phylo.beta.pair(x = regional_data_clean, tree = pruned_tree, index.family = index_family)
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
      if (aggreg_type == "all_pairwise")
      {
        beta_value <- round(mean(beta_dist, na.rm = T), 5)
      }
      
      # For focal pairwise index, pixel value is the mean pairwise beta-diversity of for all other regional sites in reference to the focal site
      if (aggreg_type == "focal_pairwise")
      {
        beta_value <- round(mean(as.matrix(beta_dist)[1, -1], na.rm = T), 5)
      }
    }
    
    # Store value for that pixel
    all_beta_values[i] <- beta_value
    
    # Show i every 100 iterations and save a back-up file
    if (i %% 100 == 0) 
    {
      cat(paste0(Sys.time(), " - ", i," on ",nrow(binary_mat),"\n"))
      # save(all_FP, file = "./outputs/Indices_maps/backup_betadiversity_moving.RData")
    }
  }
  
  # Generate a mask for terrestrial areas
  continental_mask <- (calc(proba_stack, fun = sum) >= 0) - 1
  
  # Write beta-diversity values in a raster without continental borders
  betadiversity_raster_no_borders <- proba_stack[[1]]
  betadiversity_raster_no_borders@data@values <- all_beta_values
  
  # Apply continental mask to set NA pixels to NA values to avoid issue with extending borders due to the window
  betadiversity_raster_no_borders[is.na(continental_mask[])] <- NA
  
  # Write beta-diversity values in a raster with null values as background for terrestrial areas to show terrestrial borders
  betadiversity_raster <- continental_mask
  betadiversity_raster@data@values[!(is.na(betadiversity_raster_no_borders@data@values) | is.nan(betadiversity_raster_no_borders@data@values))] <- betadiversity_raster_no_borders@data@values[!(is.na(betadiversity_raster_no_borders@data@values) | is.nan(betadiversity_raster_no_borders@data@values))]
  
  # Repair issue with min/max values
  betadiversity_raster@data@max <- min(betadiversity_raster[], na.rm = T)
  betadiversity_raster@data@max <- max(betadiversity_raster[], na.rm = T)
  
  return(betadiversity_raster)
  
  cat(paste0(Sys.time(), " - Betadiversity Computation - Done\n"))
}


### 4/ Function to compute betadiversity within regions ####

### Compute a betadiversity value for each region based on betadiversity across all sites of each region

# Can specify the way to compute a unique index from the multiple sites in the region using the "aggreg_type" argument
   # "mean_pairwise" provides the mean value between all pairwise betadiversity indices across all pairs of sites within the window 
   # "multi" provides the multi-sites betadiversity index as in Baselga et al., 2012

# Support multiple index type, family, partition and aggregation at the same time

# Can specify a subsample size to reduce computation time and limit issues with spatial autocorrelation between pixels within each region

# Inputs = Raster Stack of species/OMU probability of presence
#          Phylogeny of the clade (for PhyloBetaDiversity indices)
#          List of SpatialPolygons for each region
#          Names of regions

# Outputs = Raster Stack of all continuous betadiversity requested
#           Dataframe of betadiversity values per regions with regions in rows and indices in columns (Optional: output_type = "region_df")
#           Dataframe of betadiversity values per regions with regions and type of indices in columns as in ggplot2 (Optional: output_type = "ggplot_df")

# Plot = Maps of betadiversity indices within regions (Optional: "regions_map" = T)

### To improve
  # All types of inputs: spPolygons + sf object + raster Layer (or Stack) as a region input


compute_within_regional_betadiversity <- function (proba_stack, list_shp_regions, region_names,
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
  
  # Add region names to each SpPolygonsDF
  for (i in 1: length(list_shp_regions))
  {
    n_poly <- nrow(list_shp_regions[[i]]@data)
    list_shp_regions[[i]]@data <- cbind(list_shp_regions[[i]]@data, data.frame("region" = rep(region_names[i], times = n_poly)))
  }
  
  # Check CRS matching between shp files and Raster Stack and project if necessary
  for (i in 1:length(list_shp_regions))
  {
    CRS_check <- raster::compareCRS(list_shp_regions[[i]], proba_stack)
    if(!CRS_check)
    {
      list_shp_regions[[i]] <- spTransform(x = list_shp_regions[[i]], CRSobj = proba_stack@crs)
    }
  }
  
  # Aggregate SpatialPolygons in a single shape file
  all_shp_polygons <- do.call(raster::bind, list_shp_regions)
  
  # all_shp_polygons@data
  
  # Extract community matrix per polygons
  binary_mat_per_polygons <- raster::extract(binary_stack, all_shp_polygons)
  
  # Aggregate per regions (in case some regions encompass multiple polygons)
  binary_mat_per_regions <- list()
  for (i in 1:length(region_names))
  {
    # i <- 2
    region <- region_names[i]
    polygon_indices <- which(region == all_shp_polygons@data$region)
    
    binary_mat_per_regions[[i]] <- do.call(what = "rbind", args = binary_mat_per_polygons[polygon_indices])
  }
  
  # ### Version if all regions have a single polygon
  # all_shp_polygons <- do.call(raster::bind, list_shp_regions)
  # binary_mat_per_regions <- raster::extract(binary_stack, all_shp_polygons)
  
  ###### Need to check if no more issue with all_shp_polygons !!! #####
  ## Copy-paste the section with CRS checking and data aggregation per region
  
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
  
  # Convert into Regions x Indices df
  within_regions_df <- within_regions_ggplot_df %>% 
    tidyr::pivot_wider(data = ., names_from = c(diversity_type, index_family, beta_part, aggreg_type), names_sep = "_", values_from = index_value) %>% 
    as.data.frame(.)
  
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
        
        index_stack <- addLayer(index_stack, rasterize(x = list_shp_regions[[j]], 
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



### 5/ Function to compute betadiversity between regions ####

### Compute pairwise betadiversity between regions aggregated as unique sites and provide associated pairwise betadiversity matrices

# Can specify the way to compute a unique index from the multiple sites in the region using the "aggreg_type" argument
  # "mean_pairwise" provides the mean value between all pairwise betadiversity indices across all pairs of sites within the window 
  # "multi" provides the multi-sites betadiversity index as in Baselga et al., 2012

# Support multiple index type, family, partition and aggregation at the same time

# Inputs = Raster Stack of species/OMU probability of presence
#          Phylogeny of the clade (for PhyloBetaDiversity indices)
#          List of SpatialPolygons for each region
#          Names of regions

# Output = List of pairwise betadiversity matrices. One per index requested.


### To improve
  # All types of inputs: spPolygons + sf object + raster Layer (or Stack) as a region input


compute_pairwise_regional_betadiversity <- function (proba_stack, list_shp_regions, region_names,
                                                     diversity_type = "taxo",
                                                     phylo = NULL, index_family = "sorensen",
                                                     beta_part = "turnover")
{
  cat(paste0(Sys.time(), " - Betadiversity Computation - Start\n"))
  
  ### Prepare input data
  
  ##### Could save computation time by just binarizing within the regions of interest, and not all the study area
  
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
  
  # Add region names to each SpPolygonsDF
  for (i in 1: length(list_shp_regions))
  {
    n_poly <- nrow(list_shp_regions[[i]]@data)
    list_shp_regions[[i]]@data <- cbind(list_shp_regions[[i]]@data, data.frame("region" = rep(region_names[i], times = n_poly)))
  }
  
  # Check CRS matching between shp files and Raster Stack and project if necessary
  for (i in 1:length(list_shp_regions))
  {
    CRS_check <- raster::compareCRS(list_shp_regions[[i]], proba_stack)
    if(!CRS_check)
    {
      list_shp_regions[[i]] <- spTransform(x = list_shp_regions[[i]], CRSobj = proba_stack@crs)
    }
  }
  
  # Aggregate SpatialPolygons in a single shp file
  all_shp_polygons <- do.call(raster::bind, list_shp_regions)
  
  # all_shp_polygons@data
  
  # Extract community matrix per polygons
  binary_mat_per_polygons <- raster::extract(binary_stack, all_shp_polygons)
  
  # Aggregate per regions (in case some regions encompass multiple polygons)
  binary_mat_per_regions <- list()
  for (i in 1:length(region_names))
  {
    # i <- 2
    region <- region_names[i]
    polygon_indices <- which(region == all_shp_polygons@data$region)
    
    binary_mat_per_regions[[i]] <- do.call(what = "rbind", args = binary_mat_per_polygons[polygon_indices])
  }  
  
  # ### Version if all regions have a single polygon
  # all_shp_polygons <- do.call(raster::bind, list_shp_regions)
  # binary_mat_per_regions <- raster::extract(binary_stack, all_shp_polygons)
  
  # Aggregate regional diversity
  occ_per_regions <- purrr::map(.x = binary_mat_per_regions, .f = function (x) {apply(X = x, MARGIN = 2, FUN = sum, na.rm = T)}) 
  PA_per_regions <- purrr::map(.x = occ_per_regions, .f = function (x) {x > 0})
  binary_per_regions <- purrr::map(.x = PA_per_regions, .f = as.numeric)
  aggregated_binary_mat <- t(as.data.frame(binary_per_regions))
  row.names(aggregated_binary_mat) <- region_names
  colnames(aggregated_binary_mat) <- names(binary_stack)
  
  ### Generate the map of lists for pmap function crossing all argument combination.
  map_cross <- purrr::cross(list(diversity_type = diversity_type, index_family = index_family, beta_part = beta_part))
  
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
  
  ### Compute betadiversity values for all regions, all diversity type, all index families, all partitions, and all aggregation type
  pairwise_matrices <- purrr::pmap(.l = map_cross_revert, .f = compute_pairwise_betadiversity, regional_data = aggregated_binary_mat, phylo = pruned_tree)
  
  # Format matrix names
  matrices_names <- purrr::pmap(.l = map_cross_revert, .f = paste, sep = "_")
  matrices_names <- lapply(X = matrices_names, FUN = stringr::str_split, pattern = "_", simplify = T)
  matrices_names <- lapply(X = matrices_names, FUN = function (x) {paste0("Diversity type: ", x[,1], " ; Family: ", x[,2], " ; Partition: ", x[,3])})
  names(pairwise_matrices) <- matrices_names
  
  ### Print results
  cat(paste0("\nPairwise Betadiversity between regions\n\n"))
  print(pairwise_matrices)
  
  # Export
  return(pairwise_matrices)
  
}


### 6/ Function to map as an interaction network the within and between regions betadiversity ####

### Network of betadiversity between/within regions
  # Background = map of regions
  # Nod position = controid of regions
  # Nod size = within regions betadiversity
  # Edges width = pairwise diversity between regions

# Can specify a subsample size to reduce computation time and limit issues with spatial autocorrelation between pixels within each region (if indices need ot be computed)

# Inputs = Dataframe with betadiversity values within each region. 2 columns: the region names, the index values
#          Matrix of betadiversity values between regions
#          List of SpatialPolygons for each region
#          Names of regions
#          Labels for regions
#
#          Raster Stack of species/OMU probability of presence (to compute indices)
#          Phylogeny of the clade (to compute PhyloBetaDiversity indices)

# Output = igraph or ggraph object to plot and customize

# Plot = Network of betadiversity between/within regions


### To improve
  # All types of inputs: spPolygons + sf object + raster Layer (or Stack) as a region input
  # Make names of region optional because can use labels only instead
  # Check if the output is plotable and customizable


map_betadiversity_within_and_between_regions <- function (
  
  # If indices are already computed
  region_df_within = NULL,   # To provide directly the df with indices values within each region. 2 columns: the region names, the index values
  region_matrix_between = NULL,  # To provide directly the matrix of betadiversity values between regions
  list_shp_regions, # To provide the list of SpatialPolygons of regions to extract centroid coordinates
  region_names, # To provide the names of regions
  region_labels = NULL, # To provide labels to use instead of full region names
  
  # If indices needs to be computed
  proba_stack = NULL, # Stack of species outputs from SDM
  subsample_size = NA, # To provide number of pixels to subsample globally and regularly in order to save computing time, especially if phylobetadiversity need to be computed
  diversity_type = "taxo", # Choose between taxonomic and phylogenetic diversity
  phylo = NULL, # Provide phylogeny if phylogenetic indices need to be computed
  index_family = "sorensen", # Choose the family of indices between "jaccard" and "sorensen" 
  beta_part = "turnover", # Choose the partition of betadiversity
  aggreg_type = "mean_pairwise", # Choose the aggregation method for within betadiversity
  
  # Plot options
  plot_type = "ggplot", # To select for plot type between "ggplot" and "igraph"
  nod_cex = 40, # To adjust nod size
  edge_cex = 1) # To adjust edge width
{
  # Check that there is only one index computed at once
  if(any(!(lapply(X = list(index_family, beta_part, aggreg_type), FUN = length) == 1)))
  {
    stop("Only one value per argument for index choice")
  }
  
  # Check if indices within regions need to be computed 
  if (is.null(region_df_within)) 
  {
    cat(paste0(Sys.time(), " - Compute within regions Betadiversity"))
    
    # Compute within region betadiversity if not provided
    region_df_within <- compute_within_regional_betadiversity(proba_stack = proba_stack,
                                                              list_shp_regions = list_shp_regions,
                                                              region_names = region_names,
                                                              subsample_size = subsample_size,
                                                              diversity_type = diversity_type,
                                                              phylo = phylo, 
                                                              index_family = index_family,
                                                              beta_part = beta_part,
                                                              aggreg_type = aggreg_type,
                                                              regions_map = F,
                                                              output_type = "region_df",
                                                              quiet = T)
    region_df_within <- as.data.frame(region_df_within[[1]])
  }
  
  # Check if indices between regions need to be computed 
  if (is.null(region_matrix_between)) 
  {
    cat(paste0(Sys.time(), " - Compute between regions Betadiversity"))
    
    # Compute between regions betadiversity if not provided
    region_matrix_between <- compute_pairwise_regional_betadiversity(proba_stack = proba_stack,
                                                                     list_shp_regions = list_shp_regions,
                                                                     region_names = region_names,
                                                                     diversity_type = diversity_type,
                                                                     phylo = phylo,
                                                                     index_family = index_family,
                                                                     beta_part = beta_part)
    region_matrix_between <- region_matrix_between[[1]]
  }
  
  ### Format region_df_within
  region_df_within <- as.data.frame(region_df_within)
  row.names(region_df_within) <- region_names <- region_df_within$regions
  region_df_within <- region_df_within[, 2, drop = F]
  names(region_df_within) <- "BetaDiversity_within"
  
  ### Clean out regions with no data
  
  # Check for no data within regions
  check_NaN <- is.nan(region_df_within[, 1])
  
  removed_regions <- row.names(region_df_within)[check_NaN]
  # Print the name of removed regions
  if (length(removed_regions) == 1)
  {
    cat(paste0("\nRemove 1 region with no data:\n"))
    print(removed_regions)
  }
  if (length(removed_regions) > 1)
  {
    cat(paste0("\nRemove ",length(removed_regions)," regions with no data:\n"))
    print(removed_regions)
  }
  # Clean out regions with no data
  region_df_within <- region_df_within[!check_NaN, , drop = F]
  region_matrix_between <- region_matrix_between[!check_NaN, !check_NaN]
  list_shp_regions <- list_shp_regions[!check_NaN]
  region_names <- region_names[!check_NaN]
  region_labels <- region_labels[!check_NaN]
  
  ### Aggregate all polygons in a single shape file
  all_shp_polygons <- do.call(raster::bind, list_shp_regions)
  
  ### Extract coordinates of regions centroids
  regions_centroids <- lapply(X = list_shp_regions, FUN = rgeos::gCentroid, byid=FALSE)
  regions_centroids <- lapply(X = regions_centroids, FUN = function (x) { x@coords })
  regions_centroids <- do.call(what = "rbind", args = regions_centroids)
  row.names(regions_centroids) <- region_names
  
  # # Version when regions all have a single polygon
  # regions_centroids <- rgeos::gCentroid(spgeom = all_shp_polygons, byid=TRUE)@coords
  
  
  ### Build interaction graph
  
  # Check Connection map with geosphere (https://www.r-graph-gallery.com/connection-map.html)
  
  # Check igraph into ggraph to combine with geom_raster and others
  
  # See sf and network as MULTILINESTRING (https://r-spatial.org/r/2019/09/26/spatial-networks.html)
  
  # Build network from adjacency matrix
  network <- igraph::graph_from_adjacency_matrix(adjmatrix = region_matrix_between, mode = "undirected", weighted = T, diag = F, add.colnames = T)
  
  # Update label if needed
  if (is.null(region_labels)) { region_labels <- region_names }
  
  if (plot_type == "igraph")
  {
    # Build scale for edge color
    standardized_edges <- round(c(region_matrix_between)/ max(c(region_matrix_between)) * 200)
    
    # Build scale for nod color
    standardized_nods <- round(region_df_within[,1]/ max(region_df_within[,1]) * 200)
    
    # Basic plot with igraph
    network_plot <- plot(x = network,
                         layout = regions_centroids, # Coordinates of region centroids)
                         
                         # Nod/vertex style
                         vertex.color = colorRampPalette(c("white", "dodgerblue"))(200)[standardized_nods],
                         vertex.size = region_df_within[,1]*nod_cex,
                         vertex.label = region_labels,
                         vertex.label.font = 2,
                         vertex.label.color = "black",
                         vertex.label.cex = 1,
                         vertex.label.dist = 0, # Distance between the label and the vertex
                         
                         # Edge style
                         edge.color = colorRampPalette(c("white", "orange", "red"))(200)[standardized_edges],
                         edge.width = c(region_matrix_between)*edge_cex
    )
  }
  
  # igraph into ggraph to combine with geom_sf and others
  if (plot_type == "ggplot")
  {
    # Convert sp to sf
    all_sf_polygons <- sf::st_as_sf(x = all_shp_polygons)
    
    # # Extract network info
    # network_plot <- ggraph::ggraph(network, layout = regions_centroids)
    # network_info <- ggraph::get_edges("short")(network_plot)
    
    # GGplot
    network_plot <- ggraph::ggraph(network, layout = regions_centroids) + 
      
      # Plot background map
      geom_sf(data = all_sf_polygons) +
      
      # # Plot edges
      # geom_edge_link(aes(colour = network_info$weight,
      #                    width = network_info$weight*edge_cex),
      #                alpha = 0.8) +
      
      # Plot edges
      geom_edge_link(aes(colour = weight,
                         width = weight),
                     alpha = 0.8) +
      
      # Legends for nods/vertex
      scale_edge_colour_gradient(name = "Between BD",
                                 low = "cadetblue1",
                                 high = "cornflowerblue",
                                 na.value = "grey50") +
      
      scale_edge_width_continuous(range = c(0.2, 4*edge_cex)) +
      
      # Plot nods/vertex
      geom_node_point(aes(fill = region_df_within$BetaDiversity_within,
                          col = region_df_within$BetaDiversity_within),
                      shape = 21,
                      size = region_df_within$BetaDiversity_within * nod_cex,
                      alpha = 0.8) +
      # Plot nods/vertex labels
      geom_node_text(aes(label = region_labels), col = "black") +
      
      # scale_size_continuous(range = c(0.2, 8*nod_cex)) +
      
      # Legends for nods/vertex
      scale_fill_gradient(name = "Within BD",
                          low = "yellow",
                          high = "red",
                          na.value = "grey50") +
      
      scale_color_gradient(name = "Within BD",
                           low = "yellow",
                           high = "red",
                           na.value = "grey50") +
      
      # Select legend to display or not
      guides(edge_width = "none") +
      
      # Custom theme
      theme_graph(base_family = NA)
  }
  
  # Plot
  print(network_plot)
  
  # Export
  return(network_plot)
  
}




### 7/ Function to map the RGB visualization of pairwise betadiversity between sites ####

### Apply NMDS in 3D on pairwise betadiversity between sites and map it using the RGB scheme

# Can specify the minimum number of species in a community to be included (tend to distort the betadiversity space otherwise)

# Can specify a subsample size to reduce computation time


# Inputs = Raster Stack of species/OMU probability of presence
#          Phylogeny of the clade (to compute PhyloBetaDiversity indices)

# Output = Raster Stack of the visualization of pairwise community betadiversity in RGB color bands.

# Plots = Plot community coordinates in 2D and 3D scatterplots (Optional: "NMDS_plot" = T)
#         Plot rasters of pairwise community betadiversity in RGB color bands and a composite RGB layer (Optional: "RGB_plot" = T)


### To improve
  # Make a better plot with ggplot or tmap or other...
  # Implement 2D RGBY scheme

map_pairwise_betadiversity <- function (proba_stack, diversity_type = "taxo",
                                        phylo = NULL, index_family = "sorensen",
                                        beta_part = "turnover",
                                        min_sp = 3, # Minimum number of species to include a community
                                        subsample_size = NA, # Provide number of pixels to subsample globally and regularly in order to save computing time, especially if phylobetadiversity is computed
                                        color_space = "3D",
                                        NMDS_plot = T, # To plot community coordinates in 2D and 3D scatterplots 
                                        RGB_plot = T)  # To plot rasters of pairwise community betadiversity in RGB color bands and a composite RGB layer
{
  ### Prepare data
  
  # Subsample sites in a standardized fashion such as density of sampling remains the same for all regions
  if (!is.na(subsample_size))
  {
    proba_stack <- raster::sampleRegular(x = proba_stack, size = subsample_size, asRaster = T)
  }
  
  # If using the phylogenetic beta-diversity, need to match raster Stack and phylogeny species lists
  if (diversity_type == "phylo")
  {
    # Match raster Stack and Phylogeny species lists
    clean_stack_and_phylo <- match_stack_and_phylo(binary_stack, phylo)
    binary_stack <- clean_stack_and_phylo[[1]]
    pruned_tree <- clean_stack_and_phylo[[2]]
  } else {
    binary_stack <- proba_stack
  }
  
  # Extract community matrix
  binary_mat <- raster::getValues(binary_stack)
  
  # Remove communities with NA
  binary_mat_clean_NA <- na.omit(binary_mat)
  NA_com_indices <- as.numeric(attr(binary_mat_clean_NA, which = "na.action"))
  # Remove communities with no or less than the minimum species threshold
  empty_com_indices <- which(!(apply(X = binary_mat, MARGIN = 1, FUN = sum) > min_sp))
  # Remove both NA and empty communities
  removed_com_indices <- c(NA_com_indices, empty_com_indices)
  binary_mat_clean <- binary_mat[-removed_com_indices, ]
  
  ### Compute betadiversity with Baselga partitioning between turnover and nestedness
  
  # For taxonomical beta-diversity based on species identity
  if (diversity_type == "taxo")
  {
    beta_results <- betapart::beta.pair(x = binary_mat_clean, index.family = index_family)
  }
  
  # For phylogenetic beta-diversity based on Faith's Phylogenetic Diversity
  if (diversity_type == "phylo")
  {
    beta_results <- betapart::phylo.beta.pair(x = binary_mat_clean, tree = pruned_tree, index.family = index_family)
  }
  
  # Depending on the partition of beta-diversity
  if (beta_part == "turnover") {beta_mat <- as.matrix(beta_results[[1]])}
  if (beta_part == "nestedness") {beta_mat <- as.matrix(beta_results[[2]])}
  if (beta_part == "total") {beta_mat <- as.matrix(beta_results[[3]])}
  
  ### Convert pairwise beta-diversity matrix into 2D or 3D NMDS
  
  # ?metaMDS # NMDS
  # ?raster::plotRGB
  
  if (color_space == "2D") 
  {
    NMDS <- vegan::metaMDS(comm = beta_mat, k = 3)
    
    # ?recluster.col        # To get coordinates in 2D RGBY color space
    # ?recluster.plot.col   # To plot in 2D RGBY space
    
    
    
  }
  
  
  if (color_space == "3D") 
  {
    NMDS <- vegan::metaMDS(comm = beta_mat, k = 3)
    
    # Scale NMDS output into 0 to 255 range
    color_data <- apply(X = NMDS$points, MARGIN = 2, FUN = scales::rescale, to = c(0, 255))
    color_data <- apply(X = color_data, MARGIN = 2, FUN = round)
    
    # Add empty sites
    color_raster_data <- matrix(data = NA, nrow = nrow(binary_mat), ncol = 3)
    color_raster_data[-removed_com_indices, ] <- color_data
    colnames(color_raster_data) <- c("Red", "Green", "Blue")
    
    # Create template raster Brick
    RGB_brick <- brick(stack(proba_stack[[1]], proba_stack[[1]], proba_stack[[1]]))
    names(RGB_brick) <- c("Red", "Green", "Blue")
    # Fill with data
    RGB_brick@data@values <- color_raster_data
    # Convert to raster Stack
    RGB_stack <- stack(RGB_brick)
    
    if(NMDS_plot)
    {
      par(mfrow = c(2,2))
      
      plot(color_data[, 1:2], type = "n", xlab = "Red Band (NMDS1)", ylab = "Green Band (NMDS2)")
      points(x = color_data[, 1:2], pch = 16, 
             col = rgb(red = color_data[, 1], green = color_data[, 2], blue = 0, maxColorValue = 255))
      
      plot(color_data[, c(1,3)], type = "n", xlab = "Red Band (NMDS1)", ylab = "Blue Band (NMDS3)")
      points(x = color_data[, c(1,3)], pch = 16, 
             col = rgb(red = color_data[, 1], green = 0, blue = color_data[, 3], maxColorValue = 255))
      
      plot(color_data[, c(3,2)], type = "n", xlab = "Blue Band (NMDS1)", ylab = "Green Band (NMDS3)")
      points(x = color_data[, c(3,2)], pch = 16, 
             col = rgb(red = 0, green = color_data[, 2], blue = color_data[, 3], maxColorValue = 255))
      
      ### 3D plot in RGB
      
      plot3D::scatter3D(x = color_data[, 1], y =  color_data[, 2], z =  color_data[, 3],
                        bty = "b2", colkey = FALSE, theta = 45, phi = 25,
                        ticktype = "detailed",
                        pch = 16, alpha = 0.7,
                        colvar = rep(NA, nrow(color_data)),
                        NAcol = rgb(red = color_data[, 1], green = color_data[, 2], blue = color_data[, 3], maxColorValue = 255),
                        main = "RGB plot from NMDS", xlab = "\nRed Band",
                        ylab = "\nGreen Band", zlab = "\nBlue Band")
      par(mfrow = c(1,1))
    }
    
  }
  
  ### Add the plot of the NMDS output with the chosen color scheme
  if (RGB_plot)
  {
    red_palette <- colorRampPalette(c("white","red"))(256)
    green_palette <- colorRampPalette(c("white","green"))(256)
    blue_palette <- colorRampPalette(c("white","blue"))(256)
    
    par(mfrow = c(2,2))
    
    image(RGB_stack[[1]], col = red_palette)
    title(main = "Red band (NMDS1)", cex.main = 1.5, line = 1.5)
    
    image(RGB_stack[[2]], col = green_palette, colNA = "aliceblue")
    title(main = "Green band (NMDS2)", cex.main = 1.5, line = 1.5)
    
    image(RGB_stack[[3]], col = blue_palette, colNA = "aliceblue")
    title(main = "Blue band (NMDS3)", cex.main = 1.5, line = 1.5)
    
    plotRGB(RGB_stack, axes = T, main = "RGB Plot", cex.main = 1.5)
    # title(main = "RGB Plot", cex.main = 1.5, line = 1.5)
    
    par(mfrow = c(1,1))
    
  }
  
  # Export
  return(RGB_stack)
}

