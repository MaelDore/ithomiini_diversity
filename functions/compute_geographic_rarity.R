### Function to compute Geographic rarity from raster stack of SDM continuous outputs ###

# Author: Maël Doré
# Contact: mael.dore@gmail.com

### 2 levels

# Units/species geographic rarity weights
# Total and mean Community geographic rarity 

### Different type of indices

# Choose between "total" and "mean" community geographic rarity with the "index" argument
  # Total = sum of geographic rarity
  # Mean = standardized by richness

# Choose between "continuous" and "binary" with the "type" argument for the "weight_type" of rarity weights
  # Continuous = Compute Rarity weights based on an inverse exponential function with inflection point calibrated as the rarity threshold as such that communities host 25% of rare units/species in average
  # Binary = Apply a rarity threshold such as that communities host 25% of rare units/species in average


##### 1/ Function to compute geographic rarity from raster stack of SDM continuous outputs ####

### Wrapped-up around functions from the Rarity package (Boris Leroy, 2016))
# library("Rarity")
# citation("Rarity")
# Associated publications: 
# DOI: 10.1111/j.1752-4598.2011.00148.x
# DOI: 10.1111/ddi.12040
###

# Input = Raster Stack of units/species probability of presence

# Output = Dataframe of geographic ranges and rarity binary and continuous weights

compute_geographic_rarity <- function(proba_stack)
{
  
  ### Compute geographic ranges
  prob_mat <- raster::getValues(proba_stack) # Get matrix of probabilities of presence per units/species (columns) per community (rows)
  pixel_range <- apply(X = prob_mat, MARGIN = 2, FUN = sum, na.rm = T) # Compute expected number of occupied pixels
  pixel_area <- mean(raster::getValues(raster::area(proba_stack)), nar.rm = T) # Compute the mean area of a pixel in the stack
  area_range <- pixel_range*pixel_area  # Convert nb of pixels into range in km²
  
  ### Compute rarity weights based on geographic ranges and assemblages
  
  # Format matrix of assemblages. Rows = units/species). Columns = communities.
  assemblages <- t(tidyr::drop_na(as.data.frame(prob_mat))) # Clean communities with NA
  assemblages <- assemblages[, colSums(assemblages) >= 1] # Remove communities with less than one unit/species expected
  
  # Compute Rarity weights based on an inverse exponential function with inflection point calibrated as the rarity threshold as such that communities host 25% of rare units/species in average
  Leroy_weights_df <- Rarity::rWeights(occData = area_range, 
                                       wMethods = "W", rCutoff = "Leroy",
                                       normalised = T, rounding = 5,
                                       assemblages = assemblages)
  
  # Format final output
  geographic_rarity_df <- cbind(row.names(Leroy_weights_df), pixel_range, Leroy_weights_df)
  geographic_rarity_df <- geographic_rarity_df[, -ncol(geographic_rarity_df)]
  row.names(geographic_rarity_df) <- 1:nrow(geographic_rarity_df)
  names(geographic_rarity_df) <- c("Unit", "Range_pixel", "Range_km2", "Rarity_binary", "Rarity_weights")
  
  return(geographic_rarity_df)
}

### 2/ Function to compute community mean and total geographic rarity from SDM continuous outputs ####

# Input = Raster Stack of units/species probability of presence

# Output = Raster Layer of community geographic rarity


# Choose between "total" and "mean" community geographic rarity with the "index" argument
# Choose between "continuous" and "binary" with the "weight_type" argument for the type of rarity weights

# If using the weighted sum of continuous geographic rarity = Range-size weighted richness
# If using the weighted mean of continuous geographic rarity = Geographic rarity standardized by richness

# If using the weighted sum of binary geographic rarity = Number of rare units/species
# If using the weighted mean of binary geographic rarity = Proportion of rare units/species


compute_community_geographic_rarity <- function(proba_stack, index = "mean", weight_type = "continuous")
{
  # Compute the geographic rarity weights of units/species
  geographic_rarity_df <- compute_geographic_rarity(proba_stack)
  
  # Choose the type of weights to use
  if(weight_type == "continuous") # For weights derived from an negative exponential transformation of geographic ranges
  {
    rarity_weights <- geographic_rarity_df$Rarity_weights
  } else {  # Binary weights based on a cut-off
    rarity_weights <- geographic_rarity_df$Rarity_binary
  }
  
  # Get matrix of probabilities of presence per units/species (columns) per community (rows)
  proba_mat <- raster::getValues(proba_stack) 
  
  # Make a loop per community to compute community geographic rarity
  all_geographic_rarity <- rep(NA, nrow(proba_mat))
  for (k in 1:nrow(proba_mat)) 
  {
    proba_com <- proba_mat[k,] # Extract the probabilities of presence for the kth community
    
    # Compute Geographic rarity only if no NA is present in the community
    if (any(is.na(proba_com))) 
    {
      all_geographic_rarity[k] <- NA # if NA present, Geographic rarity = NA
    } else {
      
      if(index == "total") # Compute Geographic rarity as the weighted sum of the geographic rarity weights of units/species, weighted by their probability of presence in the community
      {
        geographic_rarity <- sum(rarity_weights * proba_com) 
      } else { # Compute mean Geographic rarity as the weighted mean of the geographic rarity weights of units/species, weighted by their probability of presence in the community
        geographic_rarity <- sum(rarity_weights * proba_com) / sum(proba_com)
      }
      all_geographic_rarity[k] <- geographic_rarity
    }
    
    # Show k every 1000 iterations and save a back-up file
    if (k %% 1000 == 0) {
      cat(paste0(Sys.time(), " - ", k," on ",nrow(proba_mat),"\n"))
      # save(all_geographic_rarity, file = "./outputs/Indices_maps/backup_geographic_rarity.RData")
    }
  }
  
  # Generate a mask for terrestrial areas
  continental_mask <- (readAll(calc(proba_stack, fun = sum)) >= 0) - 1
  
  # Write geographic rarity in a raster with null values as background for terrestrial areas
  geographic_rarity_raster <- continental_mask
  geographic_rarity_raster@data@values[!is.na(all_geographic_rarity)] <- all_geographic_rarity[!is.na(all_geographic_rarity)]
  
  # For a raster with only communities with data
  # geographic_rarity_raster <- proba_stack[[1]]
  # geographic_rarity_raster@data@values <- all_geographic_rarity
  
  # Repair issue with max value
  geographic_rarity_raster@data@max <- max(geographic_rarity_raster[], na.rm = T)
  
  return(geographic_rarity_raster)
}
