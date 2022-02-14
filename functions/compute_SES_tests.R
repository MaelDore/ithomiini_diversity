##### Functions to compute Standardized Effect-Size (SES) tests and map values and significance levels from Biodiversity index maps #####

# Author: Maël Doré
# Contact: mael.dore@gmail.com

### Contents ###

# 1/ SES-PhyloDiversity (PD) according to Species Richness (SR)
# 2/ SES-PhyoBetaDiversity (PBD) according to Taxonomic BetaDiversity (TBD)
# 3/ SES-PhyoBetaDiversity (PBD) according to Taxonomic BetaDiversity (TBD) within regions 


### SES-tests to account for deviance (Leprieur et al., 2012 ; Rosauer & Jetz, 2014). 
  # Test the significance of a phylogenetic index according to its taxonomic equivalent based only on species composition.
  # Null model = randomize names in the phylogeny to assume no effect of phylogeny on the index value
  # SES-value = (observed value - mean(null distribution)) / sd(null distribution)
  # Significance can be based on SES-value assuming Normal distribution of the null distribution or using quantiles

# Number of randomization is chosen with the "randomizations" argument. Default = 999.
# Can chose the alpha level of significance with the "alpha" argument. Default = 0.05.

# Can specify the coordinates of the community to use as example to plot the null distribution with the "example_coords" argument.


##### 1/ Function to map SES-PD according to SR #####

# Inputs = Raster Stack of species/OMU probability of presence
#          Phylogeny of the clade

# Outputs = Raster Stack with four maps
#               (a) SES-PD
#               (b) Significant areas following SES-PBD assuming Normal distribution of the null distribution
#               (c) p-values from randomized tests using quantiles
#               (d) Significant areas following p-values from randomized tests
#           Raster Stack of index values for randomized communities (optional: "export_null_rasters" = T)

# Prints = Four maps of (a) SES-PD, (b) Significant areas following SES-PD, (c) p-values from randomized tests using quantiles, (d) Significant areas following p-values from randomized tests
#          Example of null distribution for a community (Optional: "plot_distri" = T)



compute_SES_Faith_PD_vs_SR <- function (proba_stack, phylo, seed = 1, # Set seed to ensure reproducibility
                                        randomizations = 999, # To select the number of randomization
                                        alpha = 0.05,  # Threshold for significance for a two-sided test
                                        export_null_rasters = F, # Include the raster Stack of PD obtained from randomization in the output
                                        plot_distri = T, # To plot the null distribution for an "average" community as an example
                                        example_coords = NULL) # To provide the coordinates of the community to use as example, otherwise a community with median species richness is used.
{
  # Set the random seed for reproducibility
  set.seed(seed = seed)
  
  # Match raster Stack and Phylogeny species lists
  clean_stack_and_phylo <- match_stack_and_phylo(proba_stack, phylo)
  proba_stack_for_phylo <- clean_stack_and_phylo[[1]]
  pruned_tree <- clean_stack_and_phylo[[2]]
  
  # Compute observed Faith_PD
  cat(paste0(Sys.time(), " - Compute observed Faith's PD\n"))
  Faith_PD_obs <- compute_Faith_PD(proba_stack = proba_stack_for_phylo, phylo = pruned_tree, quiet = T)
  names(Faith_PD_obs) <- "Faith_PD_obs"
  
  # Loop for each randomization
  PD_stack <- stack(Faith_PD_obs)
  for (i in 1:randomizations)
  {
    cat(paste0(Sys.time(), " - Randomization - ", i," on ",randomizations,"\n"))
    
    # Randomize species in the phylogeny
    random_phylo <- pruned_tree
    random_phylo$tip.label <- sample(random_phylo$tip.label)
    
    # Compute Faith_PD for this randomized phylogeny
    random_PD <- compute_Faith_PD(proba_stack_for_phylo, random_phylo, quiet = T)
    names(random_PD) <- paste0("Faith_PD_random_", i)
    
    # Add to the final Raster Stack
    PD_stack <- addLayer(PD_stack, random_PD)
  }
  
  if(plot_distri)
  {
    # If no example community is provided
    if (is.null(example_coords))
    {
      # Sample an "average" community
      com_ranks <- rank(Faith_PD_obs[], ties.method = "random", na.last = "keep")
      sample_com_index <- which(round(median(com_ranks, na.rm = T)) == com_ranks)
      
      # Get its coordinates
      sample_com_coords <- round(raster::xyFromCell(object = Faith_PD_obs, cell = sample_com_index),3)
      
    } else { # If the example community coordinates are provided
      
      sample_com_coords <- example_coords
      
      # Use the provided coordinates to retrieve its index
      sample_com_index <- raster::cellFromXY(object = Faith_PD_obs, xy = sample_com_coords)
    }
    
    # Extract statistics
    com_null_distri <- getValues(PD_stack)[sample_com_index, ]
    com_PD_obs <- com_null_distri[1]
    mean_null_value <- mean(com_null_distri, na.rm = T)
    Q2.5_value <- quantile(com_null_distri, 0.025, na.rm = T)
    Q97.5_value <- quantile(com_null_distri, 0.975, na.rm = T)
    p_value <- round(ecdf(x = com_null_distri)(com_PD_obs), 3)
    SES_value <- round((com_PD_obs - mean_null_value) / sd(com_null_distri), 3)
    
    # Plot null distribution for an "average" community
    y_counts <- hist(com_null_distri, breaks = 30, plot = F)$counts
    hist(com_null_distri,
         breaks = 30, freq = TRUE, col = "gray", 
         main = paste0("Null distribution of the PD in a community\nLatitude: ", sample_com_coords[2], "; Longitude: ", sample_com_coords[1]), 
         xlab = "Faith's PD  [My]",
         cex.axis = 1.5, cex.lab = 1.3, cex.main = 1.6, lwd = 2)
    arrows(x0 = com_PD_obs, y0 = max(y_counts)*0.5, x1 = com_PD_obs, y1 = max(y_counts)*0.1, length = 0.1, lwd = 3)
    abline(v = mean(com_null_distri), lty = 2, lwd = 2)
    abline(v = Q2.5_value, lty = 2, lwd = 2, col = "red")
    abline(v = Q97.5_value, lty = 2, lwd = 2, col = "red")
    
    # Insert mean and quantiles legend
    legend(legend = c(paste0("Mean = ", format(round(mean_null_value, 1), nsmall = 1)), 
                      paste0("CI 2.5% = ", format(round(Q2.5_value, 1), nsmall = 1)),
                      paste0("CI 97.5% = ", format(round(Q97.5_value, 1), nsmall = 1))), 
           x = "topright", inset = c(0.02, 0.02), y.intersp = 1.2, lty = 2 , lwd = 2,
           col = c("black", "red"), cex = 1.2, bty = "o", bg = "white", box.col = NA)
    
    # Insert PD obs, SES and p-value legend
    legend(legend = c(paste0("PD obs = ", round(com_PD_obs, 2)),
                      paste0("SES-PD = ", format(SES_value, nsmall = 3)),
                      paste0("p = ", format(p_value, nsmall = 3))),
           x = "topleft", inset = c(0.02, 0.02),
           cex = 1.2, xjust = 1, bty = "o", bg = "white", box.col = NA)
  }
  
  ### Compute SES-values
  
  # Compute the mean of the null distribution
  mean_null_raster <- raster::calc(x = PD_stack, fun = mean, na.rm = T)
  # Compute the sd of the null distribution
  sd_null_raster <- raster::calc(x = PD_stack, fun = sd, na.rm = T)
  # Compute the SES-values
  SES_raster <- round(((Faith_PD_obs - mean_null_raster) / sd_null_raster), 3)
  
  ### Apply threshold to display only significant areas
  
  # Compute thresholds from alpha level
  threshold_lower <- qnorm(p = alpha/2, mean = 0, sd = 1, lower.tail = T)
  threshold_upper <- qnorm(p = alpha/2, mean = 0, sd = 1, lower.tail = F)
  
  # Create background raster
  continental_mask <- (SES_raster < Inf) - 1
  
  # Create masks
  SES_signif_lower_mask <- SES_raster < threshold_lower
  SES_signif_lower_mask[SES_signif_lower_mask[] == 0] <- NA
  SES_signif_upper_mask <- SES_raster > threshold_upper
  SES_signif_upper_mask[SES_signif_upper_mask[] == 0] <- NA
  
  # Apply masks
  SES_signif_lower <- raster::mask(x = SES_raster, mask = SES_signif_lower_mask)
  SES_signif_upper <- raster::mask(x = SES_raster, mask = SES_signif_upper_mask)
  
  # Aggregate all-in-one
  SES_signif <- mask(x = calc(stack(SES_signif_lower, SES_signif_upper), fun = sum, na.rm = T), continental_mask)
  
  ### Compute p-value
  
  # Loop per community
  PD_stack_values <- raster::getValues(PD_stack)
  all_p_values <- NA
  for (i in 1:nrow(PD_stack_values))
  {
    # Extract community data
    com_data <- PD_stack_values[i, ]
    
    # Check presence of data
    if(!any(is.na(com_data)) & (com_data[1] > 0))
    {
      # Compute p-value for this community
      com_p_value <- round(ecdf(x = com_data)(com_data[1]), 3)
      all_p_values[i] <- com_p_value
    } else {
      # Provide NA if no data
      all_p_values[i] <- NA
    }
  }
  
  # Put in raster
  p_value_raster <- continental_mask
  p_value_raster@data@values[!is.na(all_p_values)] <- all_p_values[!is.na(all_p_values)]
  
  ### Apply threshold to display only significant areas
  
  # Create masks
  p_value_signif_lower_mask <- p_value_raster < (alpha/2)
  p_value_signif_lower_mask[p_value_signif_lower_mask[] == 0] <- NA
  p_value_signif_upper_mask <- p_value_raster > (1-(alpha/2))
  p_value_signif_upper_mask[p_value_signif_upper_mask[] == 0] <- NA
  
  # Apply masks
  p_value_signif_lower <- raster::mask(x = p_value_raster, mask = p_value_signif_lower_mask)
  p_value_signif_upper <- raster::mask(x = p_value_raster, mask = p_value_signif_upper_mask)
  
  # Aggregate all-in-one
  # p_value_signif <- mask(x = calc(stack(p_value_signif_lower, p_value_signif_lower), fun = sum, na.rm = T), continental_mask)
  # Fix non-significant value to 0.5 to improve plot
  p_value_signif <- raster::cover(raster::cover(p_value_signif_lower, p_value_signif_upper), (continental_mask + 0.5))
  
  ### Plot the four raster Layers
  par(mfrow = c(2,2))
  
  plot(SES_raster, main = "SES-PD", col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
  plot(SES_signif, main = paste0("SES-PD significant areas\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
  
  plot(p_value_raster, zlim = c(0, 1), main = paste0("p-values from randomization tests\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
  plot(p_value_signif, zlim = c(0, 1), main = paste0("Significant areas from randomization tests\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
  
  par(mfrow = c(1,1))
  
  ### Export
  
  final_output <- list(SES_raster, SES_signif, p_value_raster, p_value_signif)
  # Include the Stack of randomized community if requested
  if (export_null_rasters) { final_output <- append(final_output, list(PD_stack)) }
  
  return(final_output)
  
}


##### 2/ Function to map SES-PBD according to TBD ######

### PhyloBetaDiversity based on moving window approach (See compute_moving_betadiversity)

# Need to specify the window size (width and height) in number of pixels from the focal point to the border of the window.
# By default, the value equals a fifth of the study area

# Can specify the way to compute a unique index from the multiple sites in the window using the "aggreg_type" argument
  # "all_pairwise" provides the mean value between all pairwise betadiversity indices across all pairs of sites within the window 
  # "focal_pairwise" provides the mean value between all pairwise betadiversity indices of all sites in reference to the focal site.   
  # "multi" provides the multi-sites betadiversity index as in Baselga et al., 2012

# Can specify the coordinates of the community to use as example to plot the null distribution with the "example_coords" argument.

# Inputs = Raster Stack of species/OMU probability of presence
#          Phylogeny of the clade

# Outputs = Raster Stack with four maps
#               (a) SES-PBD
#               (b) Significant areas following SES-PBD assuming Normal distribution of the null distribution
#               (c) p-values from randomized tests using quantiles
#               (d) Significant areas following p-values from randomized tests
#           Raster Stack of index values for randomized communities (optional: "export_null_rasters" = T)

# Prints = Four maps of (a) SES-PBD, (b) Significant areas following SES-PBD, (c) p-values from randomized tests using quantiles, (d) Significant areas following p-values from randomized tests
#          Example of null distribution for a community (Optional: "plot_distri" = T)


compute_SES_PBD_vs_TBD <- function (proba_stack, phylo, 
                                    index_family = "sorensen", # Chose the family of Beta-diversity indices
                                    beta_part = "turnover",    # Chose the partition
                                    aggreg_type = "all_pairwise",     # Chose the aggregation method
                                    window_width = "default", window_height = "default", # In pixels, extended from the focal point
                                    seed = 1, # Set seed to ensure reproducibility
                                    randomizations = 999, # To select the number of randomization
                                    alpha = 0.05,  # Threshold for significance for a two-sided test
                                    export_null_rasters = F, # Include the raster Stack of PBD obtained from randomization in the output
                                    plot_distri = T, # To plot the null distribution for an "average" community as an example
                                    example_coords = NULL) # To provide the coordinates of the community to use as example, otherwise a community with median species richness is used.
{
  # Set the random seed for reproducibility
  set.seed(seed = seed)
  
  # Match raster Stack and Phylogeny species lists
  clean_stack_and_phylo <- match_stack_and_phylo(proba_stack, phylo)
  proba_stack_for_phylo <- clean_stack_and_phylo[[1]]
  pruned_tree <- clean_stack_and_phylo[[2]]
  
  # Binarize output following probability ranks if needed
  if (any(!is.element(unique(as.vector(getValues(proba_stack_for_phylo))), c(NA, 0, 1))))
  {
    binary_stack <- binarize_output_with_ranks(proba_stack_for_phylo)
  } else {
    binary_stack <- proba_stack_for_phylo
  } 
  
  # Compute observed PBD based on moving window
  cat(paste0(Sys.time(), " - Compute observed PhyloBetaDiversity in a moving window\n"))
  PBD_obs <- compute_moving_betadiversity(binary_stack, diversity_type = "phylo",
                                          phylo = pruned_tree, index_family = index_family,
                                          beta_part = beta_part, aggreg_type = aggreg_type,
                                          window_width = window_width, window_height = window_height) # In pixels, extended from the focal point
  names(PBD_obs) <- "PBD_obs"
  
  # Loop for each randomization
  PBD_stack <- stack(PBD_obs)
  for (i in 1:randomizations)
  {
    cat(paste0(Sys.time(), " - Randomization - ", i," on ",randomizations,"\n"))
    
    # Randomize species in the phylogeny
    random_phylo <- pruned_tree
    random_phylo$tip.label <- sample(random_phylo$tip.label)
    
    # Compute PBD for this randomized phylogeny
    random_PBD <- compute_moving_betadiversity(binary_stack, diversity_type = "phylo",
                                               phylo = random_phylo, index_family = index_family,
                                               beta_part = beta_part, aggreg_type = aggreg_type,
                                               window_width = window_width, window_height = window_height)
    names(random_PBD) <- paste0("PBD_random_", i)
    
    # Add to the final Raster Stack
    PBD_stack <- addLayer(PBD_stack, random_PBD)
  }
  
  ### Plot the null distribution of randomized PD for an example community
  if(plot_distri)
  {
    # If no example community is provided
    if (is.null(example_coords))
    {
      # Sample an "average" community
      com_ranks <- rank(PBD_obs[], ties.method = "random", na.last = "keep")
      sample_com_index <- which(round(median(com_ranks, na.rm = T)) == com_ranks)
      
      # Get its coordinates
      sample_com_coords <- round(raster::xyFromCell(object = PBD_obs, cell = sample_com_index),3)
      
    } else { # If the example community coordinates are provided
      
      sample_com_coords <- example_coords
      
      # Use the provided coordinates to retrieve its index
      sample_com_index <- raster::cellFromXY(object = PBD_obs, xy = sample_com_coords)
    }
    
    # Extract statistics
    com_null_distri <- getValues(PBD_stack)[sample_com_index, ]
    com_PBD_obs <- com_null_distri[1]
    mean_null_value <- mean(com_null_distri, na.rm = T)
    Q2.5_value <- quantile(com_null_distri, 0.025, na.rm = T)
    Q97.5_value <- quantile(com_null_distri, 0.975, na.rm = T)
    p_value <- round(ecdf(x = com_null_distri)(com_PBD_obs), 3)
    SES_value <- round((com_PBD_obs - mean_null_value) / sd(com_null_distri), 3)
    
    # Plot null distribution for an "average" community
    y_counts <- hist(com_null_distri, breaks = 30, plot = F)$counts
    hist(com_null_distri,
         breaks = 30, freq = TRUE, col = "gray", 
         main = paste0("Null distribution of the PBD in a community\nLatitude: ", sample_com_coords[2], "; Longitude: ", sample_com_coords[1]), 
         xlab = "PhyloBetaDiversity (PBD) in a moving window",
         cex.axis = 1.5, cex.lab = 1.3, cex.main = 1.6, lwd = 2)
    arrows(x0 = com_PBD_obs, y0 = max(y_counts)*0.5, x1 = com_PBD_obs, y1 = max(y_counts)*0.1, length = 0.1, lwd = 3)
    abline(v = mean(com_null_distri), lty = 2, lwd = 2)
    abline(v = Q2.5_value, lty = 2, lwd = 2, col = "red")
    abline(v = Q97.5_value, lty = 2, lwd = 2, col = "red")
    
    # Insert mean and quantiles legend
    legend(legend = c(paste0("Mean = ", format(round(mean_null_value, 3), nsmall = 3)), 
                      paste0("CI 2.5% = ", format(round(Q2.5_value, 3), nsmall = 3)),
                      paste0("CI 97.5% = ", format(round(Q97.5_value, 3), nsmall = 3))), 
           x = "topright", inset = c(0.02, 0.02), y.intersp = 1.2, lty = 2 , lwd = 2,
           col = c("black", "red"), cex = 1.2, bty = "o", bg = "white", box.col = NA)
    
    # Insert PBD obs, SES and p-value legend
    legend(legend = c(paste0("PBD obs = ", format(round(com_PBD_obs, 3), nsmall = 3)),
                      paste0("SES-PBD = ", format(SES_value, nsmall = 3)),
                      paste0("p = ", format(p_value, nsmall = 3))),
           x = "topleft", inset = c(0.02, 0.02),
           cex = 1.2, xjust = 1, bty = "o", bg = "white", box.col = NA)
  }
  
  ### Compute SES-values
  
  # Compute the mean of the null distribution
  mean_null_raster <- raster::calc(x = PBD_stack, fun = mean, na.rm = T)
  # Compute the sd of the null distribution
  sd_null_raster <- raster::calc(x = PBD_stack, fun = sd, na.rm = T)
  # Compute the SES-values
  SES_raster <- round(((PBD_obs - mean_null_raster) / sd_null_raster), 3)
  
  ### Apply threshold to display only significant areas
  
  # Compute thresholds from alpha level
  threshold_lower <- qnorm(p = alpha/2, mean = 0, sd = 1, lower.tail = T)
  threshold_upper <- qnorm(p = alpha/2, mean = 0, sd = 1, lower.tail = F)
  
  # Create background raster
  continental_mask <- (SES_raster < Inf) - 1
  
  # Create masks
  SES_signif_lower_mask <- SES_raster < threshold_lower
  SES_signif_lower_mask[SES_signif_lower_mask[] == 0] <- NA
  SES_signif_upper_mask <- SES_raster > threshold_upper
  SES_signif_upper_mask[SES_signif_upper_mask[] == 0] <- NA
  
  # Apply masks
  SES_signif_lower <- raster::mask(x = SES_raster, mask = SES_signif_lower_mask)
  SES_signif_upper <- raster::mask(x = SES_raster, mask = SES_signif_upper_mask)
  
  # Aggregate all-in-one
  SES_signif <- mask(x = calc(stack(SES_signif_lower, SES_signif_upper), fun = sum, na.rm = T), continental_mask)
  
  ### Compute p-value
  
  # Loop per community
  PBD_stack_values <- raster::getValues(PBD_stack)
  all_p_values <- NA
  for (i in 1:nrow(PBD_stack_values))
  {
    # Extract community data
    com_data <- PBD_stack_values[i, ]
    
    # Check presence of data
    if(!any(is.na(com_data)) & (com_data[1] > 0))
    {
      # Compute p-value for this community
      com_p_value <- round(ecdf(x = com_data)(com_data[1]), 3)
      all_p_values[i] <- com_p_value
    } else {
      # Provide NA if no data
      all_p_values[i] <- NA
    }
  }
  
  # Put in raster
  p_value_raster <- continental_mask
  p_value_raster@data@values[!is.na(all_p_values)] <- all_p_values[!is.na(all_p_values)]
  
  ### Apply threshold to display only significant areas
  
  # Create masks
  p_value_signif_lower_mask <- p_value_raster < (alpha/2)
  p_value_signif_lower_mask[p_value_signif_lower_mask[] == 0] <- NA
  p_value_signif_upper_mask <- p_value_raster > (1-(alpha/2))
  p_value_signif_upper_mask[p_value_signif_upper_mask[] == 0] <- NA
  
  # Apply masks
  p_value_signif_lower <- raster::mask(x = p_value_raster, mask = p_value_signif_lower_mask)
  p_value_signif_upper <- raster::mask(x = p_value_raster, mask = p_value_signif_upper_mask)
  
  # Aggregate all-in-one
  # p_value_signif <- mask(x = calc(stack(p_value_signif_lower, p_value_signif_lower), fun = sum, na.rm = T), continental_mask)
  # Fix non-significant value to 0.5 to improve plot
  p_value_signif <- raster::cover(raster::cover(p_value_signif_lower, p_value_signif_upper), (continental_mask + 0.5))
  
  ### Plot the four raster Layers
  par(mfrow = c(2,2))
  
  plot(SES_raster, main = "SES-PBD", col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
  plot(SES_signif, main = paste0("SES-PBD significant areas\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
  
  plot(p_value_raster, zlim = c(0, 1), main = paste0("p-values from randomization tests\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
  plot(p_value_signif, zlim = c(0, 1), main = paste0("Significant areas from randomization tests\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
  
  par(mfrow = c(1,1))
  
  ### Export
  
  final_output <- list(SES_raster, SES_signif, p_value_raster, p_value_signif)
  # Include the Stack of randomized community if requested
  if (export_null_rasters) { final_output <- append(final_output, list(PBD_stack)) }
  
  return(final_output)
  
}


##### 3/ Function to compute SES-PBD according to TBD within regions ####

### PhyloBetaDiversity computed within each region (See compute_within_regional_betadiversity)

# Compute a betadiversity value for each region based on betadiversity across all sites of each region

# Can specify the way to compute a unique index from the multiple sites in the region using the "aggreg_type" argument
   # "mean_pairwise" provides the mean value between all pairwise betadiversity indices across all pairs of sites within the window 
   # "multi" provides the multi-sites betadiversity index as in Baselga et al., 2012

# Can specify a subsample size to reduce computation time and limit issues with spatial autocorrelation between pixels within each region

# Can provide the name of a region to use as example to plot the null distribution with the "example_region_name" argument.


# Inputs = Raster Stack of species/OMU probability of presence
#          Phylogeny of the clade
#          List of SpatialPolygons for each region
#          Names of regions

# Outputs = Raster Stack with four maps
#               (a) SES-PBD
#               (b) Significant areas following SES-PBD assuming Normal distribution of the null distribution
#               (c) p-values from randomized tests using quantiles
#               (d) Significant areas following p-values from randomized tests
#           Dataframe with indices values from all randomization (optional: "export_null_df" = T)
#           Dataframe with statistics from the SES-tests for each region  (optional: "export_region_stats_df" = T)

# Prints = Four maps of (a) SES-PBD, (b) Significant areas following SES-PBD, (c) p-values from randomized tests using quantiles, (d) Significant areas following p-values from randomized tests
#          Example of null distribution for a community (Optional: "plot_distri" = T)



compute_SES_PBD_vs_TBD_within_regions <- function (proba_stack, list_sp_regions, region_names,
                                                   subsample_size = NA, # Provide number of pixels to subsample globally and regularly in order to save computing time
                                                   phylo, # Provide the phylogeny
                                                   index_family = "sorensen", # Chose the family of Beta-diversity indices
                                                   beta_part = "turnover",    # Chose the partition
                                                   aggreg_type = "mean_pairwise",     # Chose the aggregation method
                                                   seed = 1, # Set seed to ensure reproducibility
                                                   randomizations = 999, # To select the number of randomization
                                                   alpha = 0.05,  # Threshold for significance for a two-sided test
                                                   export_null_df = F, # To include the df with indices values from all randomization in the output
                                                   export_region_stats_df = T, # To include the df with statistics from the SES-tests for each region
                                                   plot_distri = T, # To plot the null distribution for an "average" community as an example
                                                   example_region_name = NULL) # To provide the name of the region to use as example. Otherwise it is randomly chosen.
{
  
  # Check there is only one index computed at once
  if(any(!(lapply(X = list(index_family, beta_part, aggreg_type), FUN = length) == 1)))
  {
    stop("Only one value per argument for index choice")
  }
  
  # Set the random seed for reproducibility
  set.seed(seed = seed)
  
  # Match raster Stack and Phylogeny species lists
  clean_stack_and_phylo <- match_stack_and_phylo(proba_stack, phylo)
  proba_stack_for_phylo <- clean_stack_and_phylo[[1]]
  pruned_tree <- clean_stack_and_phylo[[2]]
  
  # Binarize output following probability ranks if needed
  if (any(!is.element(unique(as.vector(getValues(proba_stack_for_phylo))), c(NA, 0, 1))))
  {
    binary_stack <- binarize_output_with_ranks(proba_stack_for_phylo)
  } else {
    binary_stack <- proba_stack_for_phylo
  } 
  
  # Compute observed PBD within regions
  cat(paste0(Sys.time(), " - Compute observed PhyloBetaDiversity in a moving window within regions\n"))
  
  PBD_obs_results <- compute_within_regional_betadiversity(binary_stack, list_sp_regions = list_sp_regions, region_names = region_names,
                                                           subsample_size = subsample_size,
                                                           diversity_type = "phylo",
                                                           phylo = pruned_tree, index_family = index_family,
                                                           beta_part = beta_part, 
                                                           aggreg_type = aggreg_type,
                                                           regions_map = F, # To create and display a map of regions with "within" regional betadiversity values
                                                           output_type = "region_df", # To choose the type of df as output. "region_df" for a region x indices df. "ggplot_df" for a row per value. "Or"raster_only" to get just the raster Stack with all indices.
                                                           quiet = T) # To display progress at each iteration of the purrr:pmap()
  
  # Format df from the list output
  region_indices_df <- as.data.frame(PBD_obs_results[[1]])
  row.names(region_indices_df) <- region_indices_df$regions
  region_indices_df <- region_indices_df[, 2, drop = F]
  names(region_indices_df) <- "PBD_obs"
  
  # Loop for each randomization
  for (i in 1:randomizations)
  {
    cat(paste0(Sys.time(), " - Randomization - ", i," on ",randomizations,"\n"))
    
    # Randomize species in the phylogeny
    random_phylo <- pruned_tree
    random_phylo$tip.label <- sample(random_phylo$tip.label)
    
    # Compute PBD for this randomized phylogeny
    random_region_df <- compute_within_regional_betadiversity(binary_stack, list_sp_regions = list_sp_regions, region_names = region_names,
                                                              subsample_size = subsample_size,
                                                              diversity_type = "phylo",
                                                              phylo = random_phylo, index_family = index_family,
                                                              beta_part = beta_part, 
                                                              aggreg_type = aggreg_type,
                                                              regions_map = F, # To create and display a map of regions with "within" regional betadiversity values
                                                              output_type = "region_df", # To choose the type of df as output. "region_df" for a region x indices df. "ggplot_df" for a row per value. "Or"raster_only" to get just the raster Stack with all indices.
                                                              quiet = T) # To display progress at each iteration of the purrr:pmap()
    
    region_indices_df <- cbind(region_indices_df, random_region_df[[1]][, 2] )
  }
  
  
  ### Format output of null distribution
  
  # Rename columns of null distri df
  names(region_indices_df) <- c("PBD_obs", paste0("PBD_random_", 1:randomizations))
  
  # Check for regions with no data
  check_NaN <- apply(X = region_indices_df, MARGIN = 1, FUN = is.nan)
  check_NaN <- apply(X = check_NaN, MARGIN = 2, FUN = any)
  removed_regions <- row.names(region_indices_df)[check_NaN]
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
  region_indices_df <- region_indices_df[!check_NaN, ]
  list_sp_regions <- list_sp_regions[!check_NaN]
  region_names <- region_names[!check_NaN]
  
  
  ### To plot an histogram of null distribution for a region as an example
  if(plot_distri)
  {
    ### Selecting a region for the example
    
    # If no region is provided as example
    if (is.null(example_region_name))
    {
      # Sample a random region
      region_index <- sample(x = 1:nrow(region_indices_df), size = 1)
      region_null_distri <- as.numeric(region_indices_df[region_index, ])
      
      # Get its name
      example_region_name <- row.names(region_indices_df)[region_index]
      
    } else { # If the example region name is provided
      
      # Extract the target region null distribution
      region_index <- which(row.names(region_indices_df) == example_region_name)
      region_null_distri <- as.numeric(region_indices_df[region_index, ])
    }
    
    # Extract statistics
    region_PBD_obs <- region_null_distri[1]
    mean_null_value <- mean(region_null_distri, na.rm = T)
    Q2.5_value <- quantile(region_null_distri, 0.025)
    Q97.5_value <- quantile(region_null_distri, 0.975)
    p_value <- round(ecdf(x = region_null_distri)(region_PBD_obs), 3)
    SES_value <- round((region_PBD_obs - mean_null_value) / sd(region_null_distri), 3)
    
    # Plot null distribution for an "average" community
    y_counts <- hist(region_null_distri, breaks = 30, plot = F)$counts
    hist(region_null_distri,
         breaks = 30, freq = TRUE, col = "gray", 
         main = paste0("Null distribution of the PBD in ", example_region_name), 
         xlab = "PhyloBetaDiversity (PBD)",
         cex.axis = 1.5, cex.lab = 1.3, cex.main = 1.6, lwd = 2)
    arrows(x0 = region_PBD_obs, y0 = max(y_counts)*0.5, x1 = region_PBD_obs, y1 = max(y_counts)*0.1, length = 0.1, lwd = 3)
    abline(v = mean(region_null_distri), lty = 2, lwd = 2)
    abline(v = Q2.5_value, lty = 2, lwd = 2, col = "red")
    abline(v = Q97.5_value, lty = 2, lwd = 2, col = "red")
    
    # Insert mean and quantiles legend
    legend(legend = c(paste0("Mean = ", format(round(mean_null_value, 3), nsmall = 3)), 
                      paste0("CI 2.5% = ", format(round(Q2.5_value, 3), nsmall = 3)),
                      paste0("CI 97.5% = ", format(round(Q97.5_value, 3), nsmall = 3))), 
           x = "topright", inset = c(0.02, 0.02), y.intersp = 1.2, lty = 2 , lwd = 2,
           col = c("black", "red"), cex = 1.2, bty = "o", bg = "white", box.col = NA)
    
    # Insert PBD obs, SES and p-value legend
    legend(legend = c(paste0("PBD obs = ", format(round(region_PBD_obs, 3), nsmall = 3)),
                      paste0("SES-PBD = ", format(SES_value, nsmall = 3)),
                      paste0("p = ", format(p_value, nsmall = 3))),
           x = "topleft", inset = c(0.02, 0.02),
           cex = 1.2, xjust = 1, bty = "o", bg = "white", box.col = NA)
  }
  
  ### Function to compute p-value from rank in a null distribution
  compute_p_value_from_null_distribution <- function (null_distribution, obs_value = NULL, obs_index = NULL)
  {
    # If observed value is not provided, extract it from the null distribution based on the provided index
    if (is.null(obs_value))
    {
      obs_value <- null_distribution[obs_index]
    }
    
    p_value <- ecdf(x = null_distribution)(obs_value)
    
  }
  
  ### Compute stat per region from df
  region_stat_df <- region_indices_df %>% 
    mutate(mean_null_PBD = apply(X = ., MARGIN = 1, FUN = mean)) %>% 
    mutate(Q2.5_PBD = apply(X = ., MARGIN = 1, FUN = quantile, probs = 0.025)) %>% 
    mutate(Q97.5_PBD = apply(X = ., MARGIN = 1, FUN = quantile, probs = 0.975)) %>% 
    mutate(sd_PBD = apply(X = ., MARGIN = 1, FUN = sd)) %>% 
    mutate(SES_PBD = (PBD_obs - mean_null_PBD) / sd_PBD) %>% 
    mutate(p_value = apply(X = ., MARGIN = 1, FUN = compute_p_value_from_null_distribution, obs_index = 1)) %>% 
    select(PBD_obs, mean_null_PBD, Q2.5_PBD, Q97.5_PBD, sd_PBD, SES_PBD, p_value) %>% 
    round(., 3)
  
  ### Generate rasters
  
  # Aggregate SpatialPolygons for regions
  if (sum(!check_NaN) == 1)
  { # Case with only one region remaining with data
    all_sp_polygons <- list_sp_regions
    
  } else {
    # Case with multiple regions remaining with data
    all_sp_polygons <- do.call(raster::bind, list_sp_regions)
  }
  
  # Add the region stats
  all_sp_polygons@data <- data.frame(all_sp_polygons@data, regions = row.names(region_stat_df), region_stat_df)
  
  # Rasterize using the SES-PBD value
  SES_raster <- raster::rasterize(x = all_sp_polygons, y = binary_stack, 
                                  field = "SES_PBD", # Variable to use to fix the pixel values within regions
                                  background = NA)  # Value to use to fill empty cells
  
  # Rasterize using the p-values from randomization tests
  p_value_raster <- raster::rasterize(x = all_sp_polygons, y = binary_stack, 
                                      field = "p_value", # Variable to use to fix the pixel values within regions
                                      background = NA)  # Value to use to fill empty cells
  
  ### Apply threshold to display only significant areas for SES-PBD
  
  # Compute thresholds from alpha level
  threshold_lower <- qnorm(p = alpha/2, mean = 0, sd = 1, lower.tail = T)
  threshold_upper <- qnorm(p = alpha/2, mean = 0, sd = 1, lower.tail = F)
  
  # Create background raster
  continental_mask <- (SES_raster < Inf) - 1
  
  # Create masks
  SES_signif_lower_mask <- SES_raster < threshold_lower
  SES_signif_lower_mask[SES_signif_lower_mask[] == 0] <- NA
  SES_signif_upper_mask <- SES_raster > threshold_upper
  SES_signif_upper_mask[SES_signif_upper_mask[] == 0] <- NA
  
  # Apply masks
  SES_signif_lower <- raster::mask(x = SES_raster, mask = SES_signif_lower_mask)
  SES_signif_upper <- raster::mask(x = SES_raster, mask = SES_signif_upper_mask)
  
  # Aggregate all-in-one
  SES_signif <- mask(x = calc(stack(SES_signif_lower, SES_signif_upper), fun = sum, na.rm = T), continental_mask)
  
  ### Apply threshold to display only significant areas for p-values from randomization tests
  
  # Create masks
  p_value_signif_lower_mask <- p_value_raster < (alpha/2)
  p_value_signif_lower_mask[p_value_signif_lower_mask[] == 0] <- NA
  p_value_signif_upper_mask <- p_value_raster > (1-(alpha/2))
  p_value_signif_upper_mask[p_value_signif_upper_mask[] == 0] <- NA
  
  # Apply masks
  p_value_signif_lower <- raster::mask(x = p_value_raster, mask = p_value_signif_lower_mask)
  p_value_signif_upper <- raster::mask(x = p_value_raster, mask = p_value_signif_upper_mask)
  
  # Aggregate all-in-one
  # p_value_signif <- mask(x = calc(stack(p_value_signif_lower, p_value_signif_lower), fun = sum, na.rm = T), continental_mask)
  # Fix non-significant value to 0.5 to improve plot
  p_value_signif <- raster::cover(raster::cover(p_value_signif_lower, p_value_signif_upper), (continental_mask + 0.5))
  
  ### Plot the four raster Layers
  par(mfrow = c(2,2))
  
  plot(SES_raster, main = "SES-PBD", col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
  plot(SES_signif, main = paste0("SES-PBD significant areas\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
  
  plot(p_value_raster, zlim = c(0, 1), main = paste0("p-values from randomization tests\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
  plot(p_value_signif, zlim = c(0, 1), main = paste0("Significant areas from randomization tests\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
  
  par(mfrow = c(1,1))
  
  ### Export
  
  # Build the final output with the four rasters: 
  # 1/ SES-PBD
  # 2/ Significant areas following SES-PBD
  # 3/ p-values from randomized tests
  # 4/ Significant areas following p-values from randomized tests
  final_output <- list(SES_raster, SES_signif, p_value_raster, p_value_signif)
  
  # Include the df of null dstribution of PBD values if requested
  if (export_null_df) { final_output <- append(final_output, list(region_stat_df)) }
  
  # Include the df of stats per regions if requested
  if (export_region_stats_df) { final_output <- append(final_output, list(region_indices_df)) }
  
  return(final_output)
  
}
