### Function to compute Pielou_evenness from raster stack of SDM continuous outputs ###

# Author: Maël Doré
# Contact: mael.dore@gmail.com


# Input = Raster Stack of Species probability of presence

# Output = Raster Layer of Shannon's diversity

### Warning! Based on the assumption than SDM output (probability of presence/habitat suitability) are suitable proxies for species abundance


compute_Pielou_evenness <- function(proba_stack)
{
  # Function to compute Pielou's evenness index from vector of abundances (or by default, proxies)
  evenness <- function(x, na.rm) 
  {
    x <- round(x*1000, digits = 0) # Round values to avoid J > 1
    
    if (sum(x, na.rm = T) < 1)  # Cannot be computed when sum of rounded "abundances" < 1
    { 
      J <- NA
    } else {
      x <- x[x > 0] # Remove all 0 values to avoid error for log(0)
      y <- x/sum(x) # Compute frequencies
      H <- -sum(y*log(y)) # Compute H'
      Hmax <- log(sum(x)) # Compute Hmax
      J <- round(min(H/Hmax, 1), digits = 3) # Compute J
    }
    return(J) # Output
  }
  
  # Compute index raster  
  Pielou_evenness_temp <- readAll(calc(proba_stack, fun = evenness))
  
  # Generate a mask for terrestrial areas
  continental_mask <- (readAll(calc(proba_stack, fun = sum)) >= 0) - 1
  
  # Add null values as background for terrestrial areas
  Pielou_evenness <- continental_mask
  Pielou_evenness@data@values[!is.na(Pielou_evenness_temp@data@values)] <- Pielou_evenness_temp@data@values[!is.na(Pielou_evenness_temp@data@values)]
  
  # Repair issue with max value
  Pielou_evenness@data@max <- max(Pielou_evenness[], na.rm = T)
  
  return(Pielou_evenness)
}

