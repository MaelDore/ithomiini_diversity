### Function to compute richness from raster stack of SDM continuous outputs ###

# Author: Maël Doré
# Contact: mael.dore@gmail.com


# Input = Raster Stack of species/OMU/mimicry ring probability of presence

# Output = Raster Layer of species/OMU/mimicry ring richness

compute_richness <- function(proba_stack)
{
  richness <- readAll(calc(proba_stack, fun = sum))
  return(richness)
}
