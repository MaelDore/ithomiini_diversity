
##### Script 14d: Compute biodiversity indices maps, clean version with functions #####

# Author: Maël Doré
# Contact: mael.dore@gmail.com


##### List of indices computed

# 1/ Species richness
# 2/ Species Shannon's Diversity
# 3/ Species Shannon's Diversity with compatibility with species richness: (a) Raw and (b) with Jost's transformation (exponential)
# 4/ Species Evenness

# 5/ Mimicry richness
# 6/ Mimicry Shannon's Diversity
# 7/ Mimicry Shannon's Diversity with compatibility with mimicry richness: (a) Raw and (b) with Jost's transformation (exponential)
# 8/ Mimicry Evenness

# 9/  Continuous Species Rarity = Range-size Weighted Species Richness. (a) Total and (b) Mean
# 10/ Categorical Species Rarity (25% threshold) : (a) Total and (b) Proportion
# 11/ Continuous Mimicry rarity = Range size weighted mimicry richness. (a) Total and (b) Mean

# 12/ Community vulnerability (Ring mean vulnerability). Weighted or not.

### 13/ Phylogeny-based indices

# 14/ Mean pairwise Phylogenetic Distance
# 15/ Faith's Phylogenetic Diversity
# 16/ Fair-proportion : (a) Total and (b) Mean


# 17/ Ring size: (a) Mean, (b) Weighted Mean, (c) Quantiles



###

# Inputs:
  # Stack of species 'probability' of presence
  # Phylogeny

# Outputs:
  # Lots of biodiversity index raster

###


### 0/ Prepare stuff ####

# Clean environment
rm(list = ls())

# Load libraries
library(raster)

# Change temp folder for raster files
rasterOptions(tmpdir = paste0("./temp"))
# To clean regularly temp folder
unlink(list.files(path = rasterOptions()$tmpdir, all.files = T, recursive = T, include.dirs = T, full.names = T), force = T, recursive = T)


# Load summary tables with ring membership by species/OMU
# load(file = paste0("./input_data/list.sp.RData"))
# load(file = paste0("./input_data/list.models.RData"))

# Create a list of mimicry rings
# mimicry_list <- as.character(unique(list.models$Mimicry.model)) # 44 Mimicry rings

# Load raster mask for continent borders
# continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_15.rds"))

### Functions to manipulate raster data ####

# Contrasting raster function
contrasting_raster <- function(x, zmin, zmax)
{
  Community_mask <- readRDS(file = "./input_data/Map_stuff/Community_mask.rds")
  continent_mask <- readRDS(file = paste0("./input_data/Map_stuff/continent_mask.rds"))
  
  x[x[] <= zmin] <- zmin  # Fix low values
  x[x[] >= zmax] <- zmax  # Fix high values
  
  x <- mask(x, mask = Community_mask)  # Cut out values that are outside Ithomiini range
  
  y <- continent_mask + zmin  # Create final new raster from continental mask with baseline = zmin
  y[!is.na(x[])] <- x[!is.na(x[])]  # Add initial raster values
  
  return(y)
}


##### 1/ Species richness ####

### Load the complete stack of species SDM outputs
sp_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))

### 1.1/ Function to compute species richness

compute_richness <- function(proba_stack)
{
  species_richness <- readAll(calc(proba_stack, fun = sum))
  return(species_richness)
}

source(file = "./functions/compute_richness.R")

### 1.2/ Index computation ####
  
sp_richness <- compute_richness(sp_proba_stack)

plot(sp_richness)

# Save
save(sp_richness, file = paste0("./outputs/Indices_maps/sp_richness.RData"))


##### 2/ Species Shannon's diversity index #####

### Warning! Based on the assumption than SDM output (probability of presence/habitat suitability) are suitable proxies for species abundance

### Load the stack of species SDM outputs
sp_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))

### 2.1/ Function to compute Shannon's diversity index ####

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

source(file = "./functions/compute_Shannon_diversity.R")

### 2.2/ Index computation ####

sp_Shannon_diversity <- compute_Shannon_diversity(sp_proba_stack)

plot(sp_Shannon_diversity)

# Save
save(sp_Shannon_diversity, file = paste0("./outputs/Indices_maps/sp_Shannon_diversity.RData"))


##### 3/ Species Shannon's Diversity with compatibility with species richness #####

### Warning! Based on the assumption that SDM output (probability of presence/habitat suitability) are suitable proxies for species abundance

# Limit number of species/OMU integrated in the computation, corresponding to the local estimated species richness
# Can use Jost's transformation to make indices comparable. Select a similar number of species/OMU per pixel to avoid getting transformed values higher than actual species richness


### Load the complete stack of species SDM outputs
sp_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))

### 3.1/ Function to compute Shannon's diversity index with Jost's transformation ####

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

source(file = "./functions/compute_Shannon_diversity.R")

### 3.2/ Index computation ####

sp_Shannon_diversity_Jost <- compute_Shannon_diversity_Jost(sp_proba_stack, transfo = T)

plot(sp_Shannon_diversity_Jost)

# Save
save(sp_Shannon_diversity_Jost, file = paste0("./outputs/Indices_maps/sp_Shannon_diversity_Jost.RData"))


##### 4/ Species Evenness = Pielou's evenness/equitability #####

### Warning! Based on the assumption that SDM output (probability of presence/habitat suitability) are suitable proxies for species abundance

### Load the complete stack of species SDM outputs
sp_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))

### 4.1/ Function to compute Pielou's evenness/equitability

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

source(file = "./functions/compute_Pielou_evenness.R")

### 4.2/ Index computation ####

sp_Pielou_evenness <- compute_Pielou_evenness(sp_proba_stack)

plot(sp_Pielou_evenness)

# Save
save(sp_Pielou_evenness, file = paste0("./outputs/Indices_maps/sp_Pielou_evenness.RData"))


##### 5/ Mimicry richness #####


### 5.1/ Mimicry Stacks generation ####

### Generate stack of mimicry probability of presence per ring
### Generate stack of mimicry richness in species per ring

# Create a list of mimicry rings
# mimicry_list <- as.character(unique(list.models$Mimicry.model)) # 44 Mimicry rings

# Load the complete stack of species SDM outputs
sp_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))

# Load the list of mimicry ring membership per species
N_rings <- 5
mimicry_ring_membership <- mimicry_list[(round(runif(n = nlayers(sp_proba_stack), min = 1, max = N_rings), digits = 0))]

# Initiate stacks
ring_proba_stack <- stack()
ring_richness_stack <- stack()

# Loop by mimicry ring
for (i in 1:length(mimicry_list))
{ 
  # i <- 1
  
  # Load ring name
  ring <- mimicry_list[i]
  
  # Get indices of species for this ring
  species_indices <- which(mimicry_ring_membership == ring)
  
  # Extract species layers to get a stack of only species from this ring
  sp_prob_stack_per_ring <- readAll(subset(x = sp_proba_stack, subset = species_indices))
  
  ### Compute species richness for this ring
  sp_richness_per_ring <- readAll(calc(sp_prob_stack_per_ring, fun = sum))

  ### Compute probability of presence for this ring
  
  # Function to compute probability of presence of at least one species/OMUs
  aggreg_prob = function(x, na.rm) { 
    y <- 1 - prod(1 - x) # Probability of presence of ring = probability of presence of at least one species/OMU = opposite of probability of absence of all species/OMU
    return(y) # Output
  }
  
  proba_per_ring <- readAll(calc(sp_prob_stack_per_ring, fun = aggreg_prob))
  
  ### Add layers to the final stack with all rings
  ring_richness_stack <- addLayer(ring_richness_stack, sp_richness_per_ring)
  ring_proba_stack <- addLayer(ring_proba_stack, proba_per_ring)

  # Check run
  if (i %% 10 == 0) {print(i)}
}

# Name layers with ring names
names(ring_richness_stack) <- names(ring_proba_stack) <- mimicry_list[1:N_rings]

# nlayers(ring_richness_stack) # N rings in the final stack

plot(ring_richness_stack)
plot(ring_proba_stack)

# Save stacks
save(ring_richness_stack, file = paste0("./outputs/Indices_stacks/ring_richness_stack.RData"))
saveRDS(ring_richness_stack, file = paste0("./outputs/Indices_stacks/ring_richness_stack.rds"))
save(ring_proba_stack, file = paste0("./outputs/Indices_stacks/ring_richness_stack.RData"))
saveRDS(ring_proba_stack, file = paste0("./outputs/Indices_stacks/ring_richness_stack.rds"))


### 5.2/ Index computation ####

source(file = "./functions/compute_richness.R")

### Load the ring probability stack
ring_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/ring_proba_stack.rds"))

ring_richness <- compute_richness(ring_proba_stack)

plot(ring_richness)

# Save
save(ring_richness, file = paste0("./outputs/Indices_maps/ring_richness.RData"))


##### 6/ Mimicry Shannon's diversity #####

# Warning! Shannon's diversity index using number of species in each ring as abundances for mimicry rings

### Load the ring probability stack
ring_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/ring_proba_stack.rds"))

### 6.1/ Index computation #####

source(file = "./functions/compute_Shannon_diversity.R")

ring_Shannon_diversity <- compute_Shannon_diversity(ring_proba_stack)

plot(ring_Shannon_diversity)

# Save
save(ring_Shannon_diversity, file = paste0("./outputs/Indices_maps/ring_Shannon_diversity.RData"))


##### 7/ Mimicry Shannon's diversity with compatibility with mimicry richness #####

# Warning! Shannon's diversity index using number of species in each ring as abundances for mimicry rings

# Limit number of rings integrated in the computation, corresponding to the local estimated ring richness
# Can use Jost's transformation to make indices comparable. Select a similar number of mimicry rings per pixel to avoid getting transformed values higher than actual mimicry ring richness


### Load the ring probability stack
ring_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/ring_proba_stack.rds"))


### 7.1/ Index computation ####

source(file = "./functions/compute_Shannon_diversity.R")

ring_Shannon_diversity_Jost <- compute_Shannon_diversity_Jost(ring_proba_stack, transfo = T)

plot(ring_Shannon_diversity_Jost)

# Save
save(ring_Shannon_diversity_Jost, file = paste0("./outputs/Indices_maps/ring_Shannon_diversity_Jost.RData"))



##### 8/ Mimicry Pielou's Evenness/Equitability #####

# Warning! Pielou's Evenness/Equitability index using number of species in each ring as abundances for mimicry rings


### Load the ring probability stack
ring_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/ring_proba_stack.rds"))


### 8.1/ Index computation ####

source(file = "./functions/compute_Pielou_evenness.R")

ring_Pielou_evenness <- compute_Pielou_evenness(ring_proba_stack)

plot(ring_Pielou_evenness)

# Save
save(ring_Pielou_evenness, file = paste0("./outputs/Indices_maps/ring_Pielou_evenness.RData"))


##### 9/ Species continuous geographic rarity #####

# If using the weighted sum of species continuous geographic rarity = Range-size weighted species richness
# If using the weighted mean of species continuous geographic rarity = Geographic rarity standardized by species richness

### Load the stack of species SDM outputs
sp_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))


### 9.1/ Function to compute geographic rarity from raster stack of SDM continuous outputs ####

# Wrapped-up around functions from the Rarity package (Boris Leroy, 2016))
# library("Rarity")
# citation("Rarity")
# Associated publications: 
   # DOI: 10.1111/j.1752-4598.2011.00148.x
   # DOI: 10.1111/ddi.12040

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
  
source(file = "./functions/compute_geographic_rarity.R")

compute_geographic_rarity(proba_stack)

### 9.2/ Function to compute community mean and total geographic rarity from SDM continuous outputs ####

# Choose between "total" and "mean" community geographic rarity with the "index" argument
# Choose between "continuous" and "binary" with the "type" argument for the type of rarity weights

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

source(file = "./functions/compute_geographic_rarity.R")

### 9.3/ Index computation ####

sp_continuous_geographic_rarity <- compute_community_geographic_rarity(sp_proba_stack, index = "mean")

plot(sp_continuous_geographic_rarity)

# Save
save(sp_continuous_geographic_rarity, file = paste0("./outputs/Indices_maps/sp_continuous_geographic_rarity.RData"))


##### 10/ Species binary geographic rarity (25% threshold) #####

# Threshold such as communities host 25% of rare units/species in average

# If using the weighted sum of species binary geographic rarity = Number of rare species
# If using the weighted mean of species binary geographic rarity = Proportion of rare species

### Load the stack of species SDM outputs
sp_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))

### 10.1/ Index computation ####

source(file = "./functions/compute_geographic_rarity.R")

# Number of rare species
sp_binary_geographic_rarity <- compute_community_geographic_rarity(sp_proba_stack, index = "total", weight_type = "binary")
# Proportion of rare species
sp_binary_geographic_rarity <- compute_community_geographic_rarity(sp_proba_stack, index = "mean", weight_type = "binary")

plot(sp_binary_geographic_rarity)

# Save
save(sp_binary_geographic_rarity, file = paste0("./outputs/Indices_maps/sp_binary_geographic_rarity.RData"))


##### 11/ Mimicry continuous geographic rarity #####

# If using the weighted sum of mimicry ring continuous geographic rarity = Range-size weighted mimicry richness
# If using the weighted mean of mimicry ring continuous geographic rarity = Geographic rarity standardized by mimicry richness

### Load the ring probability stack
ring_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/ring_proba_stack.rds"))

### 11.1/ Index computation ####

source(file = "./functions/compute_geographic_rarity.R")

# Mean mimicry Continuous geographic rarity
ring_continuous_geographic_rarity <- compute_community_geographic_rarity(ring_proba_stack, index = "total", weight_type = "continuous")
# Total mimicry Continuous geographic rarity
ring_continuous_geographic_rarity <- compute_community_geographic_rarity(ring_proba_stack, index = "mean", weight_type = "continuous")

plot(ring_continuous_geographic_rarity)

# Save
save(ring_continuous_geographic_rarity, file = paste0("./outputs/Indices_maps/ring_continuous_geographic_rarity.RData"))



##### 12/ Community vulnerability ######

# Based on Chazot et al., 2016 (DOI: 10.1007/978-3-319-22461-9_17)

# Compute community vulnerability as the sum of rings vulnerability approximated as the inverse of the local richness of each mimicry ring
# A ring with only one species as a vulnerability of 1/1 = 1. A ring of 4 species as a vulnerability of 1/4 = 0.25.
# Final vulnerability is standardized by total ring richness such as a community mean vulnerability

### Warning! Based on the assumption that SDM output (probability of presence/habitat suitability) are suitable proxies for species abundance, 
### thus that the number of species in a ring is a suitable proxy of the abundance of this mimicry ring in the community.

### Load the stack of species richness of mimicry ring
load(file = paste0("./outputs/Indices_stacks/ring_richness_stack.RData"))

### 12.1/ Function to compute community vulnerability from stack of species richness of mimicry rings ####

compute_vulnerability <- function(ring_richness_stack)
{
  # Function to compute Vulnerability index from vector of mimicry ring abundances (or by default, proxies)
  vulnerability = function(x, na.rm) {
    x <- round(x, digits = 0) # Need to round ring richness to avoid inflation of value due to numerous rings with richness values close to 0 but not null
    if (sum(x, na.rm = T) > 0) { # Computed only if local mimicry richness >= 1 once rounded
      x <- x[x>0] # Remove all 0 values to avoid error with 1/0
      V <- sum(1/x, na.rm = T) # Compute non-standardized vulnerability = sum of mimicry ring vulnerability.
      V <- V/length(x) # Standardization by local mimicry ring richness
    }else{
      V <- NA
    }
    return(V) # Output
  }
  
  community_vulnerability <- readAll(calc(ring_richness_stack, fun = vulnerability))
  
  # Generate a mask for terrestrial areas
  continental_mask <- (readAll(calc(ring_richness_stack, fun = sum)) >= 0) - 1
  
  # Add null values as background for terrestrial areas
  vulnerability_raster <- continental_mask
  vulnerability_raster@data@values[!is.na(community_vulnerability@data@values)] <- community_vulnerability@data@values[!is.na(community_vulnerability@data@values)]
  
  # Repair issue with max value
  vulnerability_raster@data@max <- max(vulnerability_raster[], na.rm = T)

  return(vulnerability_raster)
}

source(file = "./functions/compute_vulnerability.R")

### 12.2/ Function to compute weighted vulnerability (weighted by probability of presence of each ring)

compute_weighted_vulnerability <- function(ring_richness_stack, ring_proba_stack)
{
  # Extract community data
  ring_proba_data <- getValues(ring_proba_stack)
  # Need to round ring richness to avoid inflation of value due to numerous rings with richness values close to 0 but not null
  ring_richness_data <- round(getValues(ring_richness_stack), digits = 0) 
  
  weighted_vulnerability <- vector()
  for (i in 1:nrow(ring_richness_data))  # Loop between communities/pixels
  {
    local_ring_richness <- ring_richness_data[i,] # Extract species richness of mimicry rings for this community
    if (sum(local_ring_richness, na.rm = T) > 0) # Computed only if local mimicry richness >= 1 once rounded
    { 
      local_ring_richness_clean <- local_ring_richness[local_ring_richness > 0] # Remove all 0 values to avoid error with 1/0
      ring_vulnerability <- 1/local_ring_richness_clean # Compute non-standardized vulnerabilities = inverse of mimicry ring richness
      
      local_ring_proba <- ring_proba_data[i,] # Extract the probabilities of presence of mimicry ring in this community
      local_ring_proba_clean <- local_ring_proba[local_ring_richness > 0] # Remove all ring with richness = 0
      
      # Compute mean vulnerability weighted by probabilities of presence
      weighted_vulnerability[i] <- weighted.mean(x = ring_vulnerability, w = local_ring_proba_clean)
      
    } else { # If not enough richness, set to NA.
      weighted_vulnerability[i] <- NA
    }
  }
  
  table(weighted_vulnerability)
  
  # Generate a mask for terrestrial areas
  continental_mask <- (readAll(calc(ring_richness_stack, fun = sum)) >= 0) - 1
  
  # Add null values as background for terrestrial areas
  vulnerability_raster <- continental_mask
  vulnerability_raster@data@values[!is.na(weighted_vulnerability)] <- weighted_vulnerability[!is.na(weighted_vulnerability)]
  
  # Repair issue with max value
  vulnerability_raster@data@max <- max(vulnerability_raster[], na.rm = T)
  
  return(vulnerability_raster)
}

source(file = "./functions/compute_vulnerability.R")

### 12.3/ Index computation ####

community_vulnerability <- compute_vulnerability(ring_richness_stack)
community_vulnerability <- compute_weighted_vulnerability(ring_richness_stack, ring_proba_stack)

plot(community_vulnerability)

# Save
save(community_vulnerability, file = paste0("./outputs/Indices_maps/community_vulnerability.RData"))



#################################### Phylogeny-based Indices ##############################################

##### 13/ Load common stuff for phylogeny based indices #####

library(ape)
library(picante)
library(geiger)
library(treespace)


### Create function to aggregate probabilities to higher hierachical level (aggregate pixel, or go up on a phylogenetic tree)
aggreg_prob = function(x, na.rm) { 
  y <- 1-prod(1-x) 
  return(y) # Output
}


### Load phylogeny
phylogeny_object <- readRDS(file = "./input_data/Phylogenies/Final_phylogeny.rds")

### Load the complete stack of species SDM outputs
sp_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))

### Extract only the species in the phylogeny

# Extract only the 339 species included in the phylogeny from the stacks of sp probas
sp_proba_stack_phylo <- sp_proba_stack[[phylogeny_object$tip.label]]

# Save stacks
save(sp_proba_stack_phylo, file = paste0("./outputs/Indices_stacks/sp_proba_stack_phylo.RData"))
saveRDS(sp_proba_stack_phylo, file = paste0("./outputs/Indices_stacks/sp_proba_stack_phylo.rds"))

### Load directly the sp proba stack with only species in the phylogeny
sp_proba_stack_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/sp_proba_stack_phylo.rds"))


##### 14/ Mean pairwise Phylogenetic Distance (MPD) ######

?mpd # Cannot be used because the weighting scheme does not apply to our probability of presence of pairs of species 

### 14.1/ Function to compute MPD from SDM continuous outputs

compute_MPD <- function(proba_stack, phylo)
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
  
  # Get matrix of probabilities of presence per species (columns) per community (rows)
  proba_mat <- getValues(proba_stack_for_phylo)

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
  
  # Final nb of taxonomic units (species)
  n_leaves <- length(pruned_tree$tip.label)

  # Compute the patristic phylogenetic distances
  phylo_dist_mat <- ape::cophenetic.phylo(x = pruned_tree)

  # Loop to compute indices per community
  all_MPD <- NA
  for (k in 1:nrow(proba_mat)) {
    proba_com <- proba_mat[k,] # Extract the probabilities for the k community
    
    # Compute MPD only if no NA is present in the community
    if (any(is.na(proba_com))) {
      
      all_MPD[k] <- NA # if NA present, MPD = NA
      
    } else {
      
      # Compute a matrix of probability of existence of pairs in the community
      proba_pairs_mat <- matrix(data = NA, ncol = n_leaves, nrow = n_leaves)
      for (i in 1:n_leaves) {
        for (j in 1:n_leaves) {
          pair_probas <- c(proba_com[i], proba_com[j])
          proba_pairs_mat[i,j] <- prod(pair_probas) # Probability of presence of a pair = product of probability of presence of each of the species of the pair
        }
      }
      
      # Multiply them to weight phylogenetic distance by probability of presence of the pair of species
      weighted_phylo_dist_mat <- proba_pairs_mat * phylo_dist_mat
      
      # Extract unique pair values and get the mean
      total_phylo_dist <- sum(weighted_phylo_dist_mat[lower.tri(weighted_phylo_dist_mat, diag = F)]) # Extract weighted pair values and sum them
      MPD <- total_phylo_dist/sum(proba_pairs_mat[lower.tri(proba_pairs_mat, diag = F)]) # Divide by the sum of the weights to compute the weighted mean
      all_MPD[k] <- MPD # Store the MDP of each community in the final vector
      
    }
    # Show k every 1000 iterations and save a back-up file
    if (k %% 1000 == 0) 
    {
      cat(paste0(Sys.time(), " - ", k," on ",nrow(proba_mat),"\n"))
      # save(all_MPD, file = "./outputs/Indices_Maps/MPD_backup.RData")
    }
  }
  
  # Generate a mask for terrestrial areas
  continental_mask <- (readAll(calc(proba_stack, fun = sum)) >= 0) - 1
  
  # Write community MPD in a raster with null values as background for terrestrial areas
  MPD_raster <- continental_mask
  MPD_raster@data@values[!is.na(all_MPD)] <- all_MPD[!is.na(all_MPD)]
  
  # For a raster with only communities with data
  # MPD_raster <- proba_stack[[1]]
  # MPD_raster@data@values <- all_MPD
  
  # Repair issue with max value
  MPD_raster@data@max <- max(MPD_raster[], na.rm = T)
  
  return(MPD_raster)
}

source(file = "./functions/compute_MPD.R")

### 14.2/ Index computation ####

MPD <- compute_MPD(proba_stack = proba_stack_test, phylo = phylogeny_object)

plot(MPD)

# Save
save(MPD, file = paste0("./outputs/Indices_maps/MPD.RData"))



##### 15/ Faith's Phylogenetic Diversity #####

?pd # To compute Faith's phylogenetic distance (but not working with probabilities)

### 15.1/ Function to compute Faith's PD from SDM continuous outputs ####

compute_Faith_PD <- function(proba_stack, phylo)
{
  # Create function to aggregate probabilities to higher hierarchical level (aggregate pixel, or go up on a phylogenetic tree)
  aggreg_prob = function(x, na.rm) { 
    y <- 1 - prod(1 - x) 
    return(y) # Output
  }
  
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
  
  # Get matrix of probabilities of presence per species (columns) per community (rows)
  proba_mat <- getValues(proba_stack_for_phylo)
  
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
  
  # Final nb of taxonomic units (species)
  n_leaves <- length(pruned_tree$tip.label)
  
  # Generate matrix to store info on edges/branches from the pruned phylogeny
  branches = as.data.frame(matrix(NA, nrow(pruned_tree$edge), ncol = 4)) 
  names(branches) <- c("Starting nod", "Ending nod", "Length", "Proba_presence")
  branches[, 1:2] = pruned_tree$edge # Retrieve starting and ending node 
  branches[, 3] = round(pruned_tree$edge.length, 4) # Retrieve edge length
  
  # Make a loop per community
  all_PD <- NA
  for (k in 1:nrow(proba_mat)) {
    proba_com <- proba_mat[k,] # Extract the proba for the k community
    
    # Compute PD only if no NA is present in the community
    if (any(is.na(proba_com))) {
      
      all_PD[k] <- NA # if NA present, PD = NA
      
    } else {
      
      # Compute probability of presence of each edge/branch (i.e., a species descending from this branch) in this community
      for (i in 1:nrow(branches)) {
        descending_species = geiger::tips(pruned_tree, branches[i, 2]) # Retrieve the set of species descending from branch i
        index <- which(colnames(proba_mat) %in% descending_species) # Find position of the species in the stack/matrix
        prob_edge <- aggreg_prob(proba_com[index]) # Compute probability of presence of the edge
        branches[i, 4] <- prob_edge # Store info
      }
      PD <- round(sum(branches[,3] * branches[,4]), 4) # Compute PD as the weighted sum of the edge length (weighted by the probability of presence of this edge in the community)
      all_PD[k] <- PD
    }
    
    # Show k every 1000 iterations and save a backup
    if (k %% 1000 == 0) {
      cat(paste0(Sys.time(), " - ", k," on ",nrow(proba_mat),"\n"))
      # save(all_PD, file = "./outputs/Indices_Maps/PD_backup.RData")
    }
  }
  
  # Write community PD in a raster
  PD_raster <- proba_stack[[1]]
  PD_raster@data@values <- all_PD
  
  # Repair issue with max value
  PD_raster@data@max <- max(PD_raster[], na.rm = T)
  
  return(PD_raster)
  
}

source(file = "./functions/compute_Faith_PD.R")

### 15.2/ Index computation ####

Faith_PD <- compute_Faith_PD(sp_proba_stack, phylo = phylogeny_object)

plot(Faith_PD)

# Save
save(Faith_PD, file = paste0("./outputs/Indices_maps/Faith_PD.RData"))


##### 16/ Fair-proportions #####

# Fair-Proportion : distribute evolutionary time of each edge/branchs equally among its descending species

# At species level = Sum of all modified branches length from the root to the species. The most original species have the highest scores
# At community level = Sum of FP values of each species present in the community. The community with the most original species has the highest scores
# Can be standardized by the number of species to get the mean species FP index

### 16.1/ Function to compute species FP from a phylogeny ####

compute_species_FP <- function(phylo)
{
  # Generate table that store information on branches
  phylo_branches_df = as.data.frame(matrix(NA, nrow(phylo$edge), ncol = 5)) 
  names(phylo_branches_df) <- c("Starting_node","Ending_node","Branch_Length","Nb_desc_sp","FP")
  
  # Retrieve starting and ending nodes for each branches 
  phylo_branches_df[, c("Starting_node","Ending_node")] = phylo$edge 
  # Retrieve edge length
  phylo_branches_df$Branch_Length = phylo$edge.length 
  # Retrieve number of descending species/leaves (terminal nods)
  for (i in 1:nrow(phylo$edge)) 
  { 
    phylo_branches_df$Nb_desc_sp[i] = length(geiger::tips(phylo, phylo_branches_df$Ending_node[i]))
  }
  # Divide branch length by nb of descendant species for FP
  phylo_branches_df$FP = phylo_branches_df$Branch_Length/phylo_branches_df$Nb_desc_sp 

  # Generate table to store species Fair-Proportions
  species_FP_df <- data.frame(Taxon = phylo$tip.label, FP = NA)
  
  # Compute species Fair-Proportions
  for (i in 1:length(phylo$tip.label))
  {
    # Retrieve species name
    sp <- phylo$tip.label[i] 
    
    # Get indices of nods on the path from the root to the species terminal nod
    sp_path <- ape::nodepath(phy = phylo, 
                             from = phylo_branches_df$Starting_node[which.max(phylo_branches_df$Nb_desc_sp)], # Root
                             to = which(phylo$tip.label == sp)) # Tip
    
    # Compute FP index for the species by summing FP values of the branches on the path to the species tip in the phylogeny 
    species_FP_df$FP[species_FP_df$Taxon == sp] <- sum(phylo_branches_df[which(phylo_branches_df$Ending_node %in% sp_path), "FP"]) 
  }
  
  return(species_FP_df)
}

source(file = "./functions/compute_Fair_Proportions.R")

compute_species_FP(phylo = phylogeny_object)

### 16.2/ Function to compute community mean and total FP from SDM continuous outputs ####

# Choose between "total" FP and "mean" FP with the "index" argument

compute_community_FP <- function(proba_stack, phylo, index = "total")
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
  
  # Get matrix of probabilities of presence per species (columns) per community (rows)
  proba_mat <- getValues(proba_stack_for_phylo)
  
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
  
  # Compute the species FP df
  species_FP_df <- compute_species_FP(phylo = phylo)

  # Make a loop per community to compute community FP
  all_FP <- rep(NA, nrow(proba_mat))
  for (k in 1:nrow(proba_mat)) 
  {
    proba_com <- proba_mat[k,] # Extract the probabilities of presence for the kth community
    
    # Compute FP only if no NA is present in the community
    if (any(is.na(proba_com))) 
    {
      all_FP[k] <- NA # if NA present, FP = NA
    } else {
      
      if(index == "total") # Compute FP as the weighted sum of the FP index of species, weighted by their probability of presence in the community
      {
        FP <- sum(species_FP_df$FP * proba_com) 
      } else { # Compute mean FP as the weighted mean of the FP index of species, weighted by their probability of presence in the community
        FP <- sum(species_FP_df$FP * proba_com) / sum(proba_com)
      }
      all_FP[k] <- FP
    }
    
    # Show k every 1000 iterations and save a back-up file
    if (k %% 1000 == 0) {
      cat(paste0(Sys.time(), " - ", k," on ",nrow(proba_mat),"\n"))
      # save(all_FP, file = "./outputs/Indices_maps/backup_FP.RData")
    }
  }
  
  # Write community FP in a raster
  FP_raster <- proba_stack[[1]]
  FP_raster@data@values <- all_FP
  
  # Repair issue with max value
  FP_raster@data@max <- max(FP_raster[], na.rm = T)
  
  return(FP_raster)
}

source(file = "./functions/compute_Fair_Proportions.R")

### 16.3/ Index computation ####

community_FP <- compute_community_FP(proba_stack[[1:10]], phylo = phylogeny_object)

plot(community_FP)

# Save
save(community_FP, file = paste0("./outputs/Indices_maps/community_FP.RData"))




##### 17/ Ring size: (a) Mean, (b) Weighted Mean, (c) Quantiles #####

### Load the stack of species richness of mimicry ring
load(file = paste0("./outputs/Indices_stacks/ring_richness_stack.RData"))

### 17.1/ Function to compute ring size from stack of species richness of mimicry rings ####

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
  
source(file = "./functions/compute_ring_size.R")

### 17.2/ Index computation ####

ring_size <- compute_ring_size(ring_richness_stack, ring_proba_stack, type = "mean")
ring_size <- compute_ring_size(ring_richness_stack, type = "quantile", quantile = 0.5)
ring_size <- compute_ring_size(ring_richness_stack, type = "quantile", quantile = 0.9)
ring_size <- compute_ring_size(ring_richness_stack, type = "quantile", quantile = 0.1)
ring_size <- compute_ring_size(ring_richness_stack, ring_proba_stack, type = "weighted_mean")

plot(ring_size)

# Save
save(ring_size, file = paste0("./outputs/Indices_maps/ring_size.RData"))


