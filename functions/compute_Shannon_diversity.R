### Functions to compute Shannon's diversity from raster stack of SDM continuous outputs ###

# Author: Maël Doré
# Contact: mael.dore@gmail.com


# Input = Raster Stack of Species probability of presence

# Output = Raster Layer of Shannon's diversity

### Warning! Based on the assumption that SDM output (probability of presence/habitat suitability) are suitable proxies for species abundance


##### 1/ Classic version #####

compute_Shannon_diversity <- function(proba_stack)
{
  # Function to compute Shannon's diversity index from vector of abundances (or by default, proxies)
  shannon = function(x, na.rm) {
    x <- x[x>0] # Remove 0 values to avoid error for log(0)
    y <- x/sum(x) # Compute frequencies of each OMU
    h <- -sum(y*log(y)) # Compute Shannon's H'
    return(h) # Output
  }
  
  Shannon_diversity <- readAll(calc(proba_stack, fun = shannon))
  return(Shannon_diversity)
}


##### 2/ Jost transformation version #####

# Limit number of species/OMU integrated in the computation, corresponding to the local estimated species richness
# Can use Jost's transformation to make indices comparable. Select a similar number of species/OMU per pixel to avoid getting transformed values higher than actual species richness

compute_Shannon_diversity_Jost <- function(proba_stack, transfo = T, na.rm)
{
  # Function to compute Shannon's diversity index from vector of abundances (or by default, proxies)
  shannon_compatible <- function(x, na.rm) 
  {
    rich <- sum(x) # Local estimated species richness
    
    if (is.na(rich)) 
    {
      h <- NA # Case with NA
    } else {
      
      if (round(rich, 0) == 0) 
      {
        h <- 0 # Case with no species/OMU. Skip computation to avoid error such as log(0) = Inf
      } else { # Regular case
        
        x <- x[order(x, decreasing = T)] # Order the probabilities by decreasing values
        x <- head(x, n = round(rich, 0)) # Extract only the N highest probabilities. N = Rounded estimated richness
        y <- x/sum(x) # Compute frequencies
        h <- -sum(y*log(y)) # Compute H'
      }
    }
    
    # Output
    return(h) 
  }
  
  Shannon_diversity_Jost <- readAll(calc(proba_stack, fun = shannon_compatible))
  
  # If Jost transformation is to be applied, use exponential
  if (transfo)
  {
    Shannon_diversity_Jost <- exp(Shannon_diversity_Jost)
  }
  
  return(Shannon_diversity_Jost)
}

