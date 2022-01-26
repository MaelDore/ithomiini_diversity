### Function to compute ring size from stack of species richness (and probability of presence) of mimicry rings ###

# Author: Maël Doré
# Contact: mael.dore@gmail.com

### Compute mimicry ring size in community ###

# (a) Mean as the local number of species divided by the local number of mimicry rings
# (b) Weighted Mean as the mean size of mimicry rings weighted by their probability of presence
# (c) Quantiles as the ith ring size, excluding the one rounding to zero

# Choose between "mean" and "weighted_mean" and "quantile" with the "type" argument
# Quantile must be provided as a probability (between 0 and 1) in the "quantile" argument


# Input = Raster Stack of species richness of mimicry rings
#         Raster Stack of probability of presence of mimicry rings (not necessary for quantiles)

# Output = Raster Layer of mean/quantile for mimicry ring size


compute_ring_size <- function(ring_richness_stack, ring_proba_stack = NULL, type = "mean", quantile = NULL)
{
  # If simple mean
  if (type == "mean")
  {
    # Compute species richness
    species_richness <- readAll(calc(ring_richness_stack, fun = sum))
    # Compute ring richness
    ring_richness <- readAll(calc(ring_proba_stack, fun = sum))
    # Compute mean ring size as species richness / ring richness
    ring_size_raster <- species_richness/ring_richness
    all_ring_size <- getValues(ring_size_raster)
    
  } else {  # If not, need to loop by community/pixel
    
    # Extract community data
    ring_richness_data <- getValues(ring_richness_stack)
    if (type == "weighted_mean") {ring_proba_data <- getValues(ring_proba_stack)}
    
    # Initiate vector
    all_ring_size <- vector()
    # Loop between communities/pixels
    for (i in 1:nrow(ring_richness_data))  
    {
      # i <- 57833
      
      # Extract local ring richness
      local_ring_richness <- ring_richness_data[i,]   
      
      # If weighted mean ring size
      if (type == "weighted_mean")
      {
        # Extract local ring probabilities of presence as weights
        local_ring_proba <- ring_proba_data[i,]
        # Compute mean ring size weighted by probabilities of presence
        all_ring_size[i] <- weighted.mean(x = local_ring_richness, w = local_ring_proba)
      }
      
      # If provided a quantile
      if (type == "quantile")  
      {
        # Need to round ring richness to avoid deflation of value due to numerous rings with richness values close to 0 but not null
        local_ring_richness <- round(local_ring_richness, digits = 0) 
        # Remove null values
        local_ring_richness_clean <- local_ring_richness[local_ring_richness != 0]
        
        # Extract quantile
        all_ring_size[i] <- quantile(x = local_ring_richness_clean, probs = quantile, na.rm = T)
      }
      
      # Show i every 10,000 iterations and save a back-up file
      if (i %% 10000 == 0) 
      {
        cat(paste0(Sys.time(), " - ", i," on ",nrow(ring_richness_data),"\n"))
        # save(all_ring_size, file = "./outputs/Indices_Maps/ring_size_backup.RData")
      }
    }
  }
  
  # Generate a mask for terrestrial areas
  continental_mask <- (readAll(calc(ring_richness_stack, fun = sum)) >= 0) - 1
  
  # Add null values as background for terrestrial areas
  ring_size_raster <- continental_mask
  ring_size_raster@data@values[!is.na(all_ring_size)] <- all_ring_size[!is.na(all_ring_size)]
  
  # # For a raster with only communities with data
  # ring_size_raster <- ring_richness_stack[[1]]
  # ring_size_raster@data@values <- all_ring_size
  
  # Repair issue with max value
  ring_size_raster@data@max <- max(ring_size_raster[], na.rm = T)
  
  return(ring_size_raster)
}
