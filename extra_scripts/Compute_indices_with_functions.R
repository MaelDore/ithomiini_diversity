
##### Compute biodiversity indices maps, clean version with functions #####

# Author: Maël Doré
# Contact: mael.dore@gmail.com


##### See my notes in the dedicated "Article" folder for stuff to improve ####


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

# 18/ Betadiversity
# 18.2/ Betadiversity in reference to a focal point
# 18.3/ Betadiversity with a moving window
# 18.4/ Betadiversity within regions
# 18.5/ Betadiversity between regions
# 18.6/ Network graph for wihtin/between betadiversity of regions
# 18.7/ Pairwise betadiversity in RGB space

# 19/ Residuals
# 19.1/ From continuous indices
# 19.3/ From regional indices

# 20/ Standardized Effect-Size indices (randomization tests based on phylogeny permutations)
# 20.1/ SES-PD
# 20.2/ SES-PBD in moving window
# 20.3/ SES-PBD within regions



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
  # Create raster with null values for terrestrial areas
  continental_mask <- (calc(x, fun = sum) >= 0) - 1
  # plot(continental_mask)
  
  # Create raster with NA values outside of the species range (assuming a sum of zero is outside the range)
  community_mask_PA <- (calc(x, fun = sum) > 0) - 1
  community_mask <- community_mask_PA
  community_mask[community_mask[] < 0] <- NA
  # plot(community_mask)
  
  x[x[] <= zmin] <- zmin  # Fix low values
  x[x[] >= zmax] <- zmax  # Fix high values
  
  x <- mask(x, mask = community_mask)  # Cut out values that are outside species range
  
  y <- continental_mask + zmin  # Create final new raster from continental mask with baseline = zmin
  y[!is.na(x[])] <- x[!is.na(x[])]  # Add initial raster values
  
  return(y)
}


##### 1/ Species richness ####

### Load the complete stack of species SDM outputs
sp_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))

### 1.1/ Function to compute species richness

compute_richness <- function(proba_stack)
{
  richness <- calc(proba_stack, fun = sum)
  return(richness)
}

source(file = "./functions/compute_richness.R")

### 1.2/ Index computation ####

sp_richness <- compute_richness(sp_proba_stack)

plot(sp_richness)

# Save
save(sp_richness, file = paste0("./outputs/Indices_maps/sp_richness.RData"))

# Load
load(file = paste0("./outputs/Indices_maps/sp_richness.RData"))

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

### 5.1.1/ Function to generate raster Stacks of mimicry ring richness and mimicry ring probability of presence

generate_mimicry_ring_stack <- function (proba_stack, # Species or OMU continuous (or binary) raster Stack of SDM outputs
                                         mimicry_ring_membership, # Character string vector of mimicry ring membership per species/OMU in the same order than the Stack
                                         stack_type = c("proba", "richness"), quiet = F)
{
  if (!quiet) { cat(paste0(Sys.time(), " - Generation of mimicry ring stacks - Start\n")) }
  
  ### Initiate stacks
  if ("proba" %in% stack_type) { ring_proba_stack <- stack() }
  if ("richness" %in% stack_type) { ring_richness_stack <- stack() }

  ### Create function to aggregate probabilities
  if ("proba" %in% stack_type) 
  {
    # Function to compute probability of presence of at least one species/OMUs
    aggreg_prob <- function(x, na.rm) 
    { 
      y <- 1 - prod(1 - x) # Probability of presence of ring = probability of presence of at least one species/OMU = opposite of probability of absence of all species/OMU
      return(y) # Output
    }
  }

  # Generate list of mimicry_rings
  mimicry_list <- unique(mimicry_ring_membership)
  
  ### Loop by mimicry ring
  for (i in 1:length(mimicry_list))
  { 
    # i <- 1
    
    # Load ring name
    ring <- mimicry_list[i]
    
    # Get indices of species/OMU for this ring
    unit_indices <- which(mimicry_ring_membership == ring)
    
    # Extract species/OMU layers to get a stack of only species from this ring
    prob_stack_per_ring <- subset(x = proba_stack, subset = unit_indices)
    
    ### For probability stack
    if ("proba" %in% stack_type) 
    { 
      # Compute probability of presence for this ring
      proba_per_ring <- calc(prob_stack_per_ring, fun = aggreg_prob) 
      
      # Add layers to the final stack with all rings
      ring_proba_stack <- addLayer(ring_proba_stack, proba_per_ring)
    }
    
    ### For richness stack
    if ("richness" %in% stack_type) 
    { 
      # Compute species/OMU richness for this ring
      richness_per_ring <- calc(prob_stack_per_ring, fun = sum)
      
      # Add layers to the final stack with all rings
      ring_richness_stack <- addLayer(ring_richness_stack, richness_per_ring)
    }
    
    # Check run
    if ((i %% 10 == 0) & (!quiet))
    {
      cat(paste0(Sys.time(), " - Mimicry ring ", i," on ",length(mimicry_list),"\n"))
    }
  }
  
  ### Name layers with ring names
  if ("proba" %in% stack_type) { names(ring_proba_stack) <- mimicry_list }
  if ("richness" %in% stack_type) { names(ring_richness_stack) <- mimicry_list }

  # nlayers(ring_richness_stack) # N rings in the final stack
  
  # plot(ring_richness_stack)
  # plot(ring_proba_stack)
  
  ### Export

  # If only one type of stack is requested, provide the output as Raster Stack
  if (length(stack_type) == 1)
  {
    if (stack_type == "proba") { final_output <- ring_proba_stack }
    if (stack_type == "richness") { final_output <- ring_richness_stack }
    
  } else { # If the two types of stack are required, provide a list as output
    final_output <- list(ring_proba_stack, ring_richness_stack)
    names(final_output) <- c("ring_proba_stack", "ring_richness_stack")
  }

  return(final_output)
}

### 5.1.2/ Generate the mimicry ring stacks

# Load OMU stack and mimicry membership info
OMU_proba_stack <- readRDS("./outputs/Indices_stacks/All_OMU_stack_Jaccard.80.rds")
OMU_df <- readRDS("./input_data/list.models.rds")

# Check names and order of OMUs are identical in stack and dataframe
identical(x = names(OMU_proba_stack), OMU_df$Tag.model)

# Generate stacks
mimicry_ring_stacks <- generate_mimicry_ring_stack (proba_stack = OMU_proba_stack, # Species or OMU continuous (or binary) raster Stack of SDM outputs
                                                    mimicry_ring_membership = OMU_df$Mimicry.model, # Character string vector of mimicry ring membership per species/OMU in the same order than the Stack
                                                    stack_type = c("proba", "richness"))

# Extract both stacks individually
ring_proba_stack <- mimicry_ring_stacks$ring_proba_stack
ring_richness_stack <- mimicry_ring_stacks$ring_richness_stack

# Plot
plot(ring_proba_stack)
plot(ring_richness_stack)

# Save stacks
save(ring_proba_stack, file = paste0("./outputs/Indices_stacks/ring_richness_stack.RData"))
saveRDS(ring_proba_stack, file = paste0("./outputs/Indices_stacks/ring_richness_stack.rds"))
save(ring_richness_stack, file = paste0("./outputs/Indices_stacks/ring_richness_stack.RData"))
saveRDS(ring_richness_stack, file = paste0("./outputs/Indices_stacks/ring_richness_stack.rds"))



### 5.2/ Mimicry richness computation ####

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
aggreg_prob <- function(x, na.rm) 
{ 
  y <- 1-prod(1-x) 
  return(y) # Output
}


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
  
  # Rearrange stack so the species are in the same order than in the phylogeny
  proba_stack_for_phylo <- proba_stack_for_phylo[[c(pruned_tree$tip.label)]]
  
  # Return cleaned stack and phylogeny in a list
  cleaned_stack_and_phylo <- list(proba_stack_for_phylo, pruned_tree)
  return(cleaned_stack_and_phylo)
}

test <- match_stack_and_phylo(proba_stack, phylogeny_object)
clean_stack <- test[[1]]
pruned_tree <- test[[2]]

# Get matrix of probabilities of presence per species (columns) per community (rows)
proba_mat <- getValues(proba_stack_for_phylo)

# Final nb of taxonomic units (species)
n_leaves <- length(pruned_tree$tip.label)


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

compute_MPD <- function(proba_stack, phylo, quiet = F)
{
  # Match raster Stack and Phylogeny species lists
  clean_stack_and_phylo <- match_stack_and_phylo(proba_stack, phylo)
  proba_stack_for_phylo <- clean_stack_and_phylo[[1]]
  pruned_tree <- clean_stack_and_phylo[[2]]
  
  # Get matrix of probabilities of presence per species (columns) per community (rows)
  proba_mat <- getValues(proba_stack_for_phylo)
  
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
    if ((k %% 1000 == 0) & !quiet) 
    {
      cat(paste0(Sys.time(), " - ", k," on ",nrow(proba_mat),"\n"))
      # save(all_MPD, file = "./outputs/Indices_Maps/MPD_backup.RData")
    }
  }
  
  # Generate a mask for terrestrial areas
  continental_mask <- (calc(proba_stack, fun = sum) >= 0) - 1
  
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

MPD <- compute_MPD(proba_stack = sp_proba_stack, phylo = Ithomiini_phylo)

plot(MPD)

# Save
save(MPD, file = paste0("./outputs/Indices_maps/MPD.RData"))


names(sp_proba_stack)

##### 15/ Faith's Phylogenetic Diversity #####

?pd # To compute Faith's phylogenetic distance (but not working with probabilities)

### 15.1/ Function to compute Faith's PD from SDM continuous outputs ####

compute_Faith_PD <- function(proba_stack, phylo, quiet = F)
{
  # Create function to aggregate probabilities to higher hierarchical level (aggregate pixel, or go up on a phylogenetic tree)
  aggreg_prob <- function(x, na.rm) 
  { 
    y <- 1 - prod(1 - x) 
    return(y) # Output
  }
  
  # Match raster Stack and Phylogeny species lists
  clean_stack_and_phylo <- match_stack_and_phylo(proba_stack, phylo)
  proba_stack_for_phylo <- clean_stack_and_phylo[[1]]
  pruned_tree <- clean_stack_and_phylo[[2]]
  
  # Get matrix of probabilities of presence per species (columns) per community (rows)
  proba_mat <- getValues(proba_stack_for_phylo)
  
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
    if ((k %% 1000 == 0) & !quiet) 
    {
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

# Choose between "total" FP and "mean" FP with the "index" argument to select which index to compute
# Choose between "full" and "pruned" with the "phylo_type" argument to select which phylogeny to use to compute species FP

compute_community_FP <- function(proba_stack, phylo, phylo_type = "full", index = "total")
{
  # Match raster Stack and Phylogeny species lists
  clean_stack_and_phylo <- match_stack_and_phylo(proba_stack, phylo)
  proba_stack_for_phylo <- clean_stack_and_phylo[[1]]
  
  # Get matrix of probabilities of presence per species (columns) per community (rows)
  proba_mat <- getValues(proba_stack_for_phylo)
  
  # By default, compote the species FP on the full phylogeny
  if (phylo_type == "full") 
  {
    # Compute the species FP based on the full phylogeny (because the pruned tree may exhibit sampling bias)
    species_FP_df <- compute_species_FP(phylo = phylo)
    # Keep only the species present in the stack
    species_FP_df <- species_FP_df[species_FP_df$Taxon %in% names(proba_stack_for_phylo), ]
  }
  
  # If requested, use the pruned phylogeny with only the species from the stack to compute the species FP
  if (phylo_type == "pruned")
  {
    pruned_tree <- clean_stack_and_phylo[[2]]
    species_FP_df <- compute_species_FP(phylo = pruned_tree)
  }
  
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


##### 18/ Beta diversity #####

library("betapart")

?betapart::phylo.beta.pair()
?betapart::phylo.beta.multi()

### 18.1/ Binarize SDM continuous outputs following ranks ####

### 18.1.1/ Function to binarize following SDM scores ranking ####

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

source(file = "./functions/compute_betadiversity.R")

### 18.1.2/ Binarize SDM outputs following ranks ####

sp_binary_stack <- binarize_output_with_ranks(sp_proba_stack)
# binary_stack <- binarize_output_with_ranks(proba_stack_test)

plot(sp_binary_stack)

# Save
save(sp_binary_stack, file = paste0("./outputs/Indices_maps/sp_binary_stack.RData"))

load(file = paste0("./outputs/Indices_maps/sp_binary_stack.RData"))

### 18.2/ Betadiversity in reference to a focal point ####

### 18.2.1/ Function to compute betadiversity in reference to a focal point ####

##### Need to be improved to 
# Allow to provide a sp polygon or sf object to define the focal region
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

source(file = "./functions/compute_betadiversity.R")

### 18.2.2/ Compute focal betadiversity maps ####

# Find coordinates of the richest pixel to use as focal point
richest_com_coords <- raster::xyFromCell(object = sp_proba_stack, cell = which.max(apply(X = sp_proba_stack[], MARGIN = 1, FUN = sum, na.rm = T)))

focal_betadiversity <- compute_focal_betadiversity(sp_binary_stack,
                                                   index_family = "sorensen", 
                                                   beta_part = "turnover", 
                                                   focal_site_coords = richest_com_coords)

plot(focal_betadiversity)

focal_betaphylodiversity <- compute_focal_betadiversity(sp_binary_stack,
                                                        diversity_type = "phylo",
                                                        phylo = phylogeny_object,
                                                        index_family = "jaccard", 
                                                        beta_part = "turnover", 
                                                        focal_site_coords = richest_com_coords)

plot(focal_betaphylodiversity)

# Save
save(focal_betadiversity, file = paste0("./outputs/Indices_maps/focal_betadiversity.RData"))
save(focal_betaphylodiversity, file = paste0("./outputs/Indices_maps/focal_betaphylodiversity.RData"))


### 18.3/ Betadiversity with a moving window ####

### 18.3.1/ Function to compute betadiversity with a moving window ####

# Size of the window, default as a proportion
# Multi vs. mean pairwise

### Can be improved by adding the option to compute for combination of parameters, all at the same time...

compute_moving_betadiversity <- function (proba_stack, diversity_type = "taxo",
                                          phylo = NULL, index_family = "sorensen",
                                          beta_part = "turnover", aggreg_type = "all_pairwise",
                                          resolution_factor = 1,  # Provide a factor to aggregate data if needed (could be necessary to reduce computation time for wide window size)
                                          window_width = "default", window_height = "default", # In pixels. Must be an odd value if provided. Number of pixels after aggregation (if required)
                                          border_rule = NULL)  # Add an option on how to deal with the borders
{
  cat(paste0(Sys.time(), " - Betadiversity Computation - Start\n"))
  
  ### Aggregate data if needed
  
  if (resolution_factor != 1)
  {
    # Function to compute probability of presence of at least one species/OMUs
    aggreg_prob <- function(x, na.rm) 
    { 
      y <- 1 - prod(1 - x) # Probability of presence of ring = probability of presence of at least one species/OMU = opposite of probability of absence of all species/OMU
      return(y) # Output
    }
    
    proba_stack <- raster::aggregate(x = proba_stack,
                                     fact = resolution_factor,
                                     fun = aggreg_prob)
  }

  ### Adjust default window size to a fifth of the (aggregated) raster width and height
  
  if (window_width == "default")
  {
    window_width <- round(proba_stack@ncols/10)
    windows_width <- windows_width * 2 + 1
    if (window_width == 0) {window_width <- 3}
  }
  if (window_height == "default")
  {
    window_height <- round(proba_stack@nrows/10)
    window_height <- windows_width * 2 + 1
    if (window_height == 0) {window_height <- 3}
  }
  
  ### If provided, check that the window size parameters are odd integers
  if (((window_width + 1) %% 2) != 0) 
  {
    stop("Window width must be an odd integer")
  }
  if (((window_height + 1) %% 2) != 0) 
  {
    stop("Window height must be an odd integer")
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
                                            row = RowCol[1] - ((window_width - 1) / 2),
                                            col = RowCol[2] - ((window_height - 1) / 2),
                                            nrow = window_width,
                                            ncol = window_height)
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

source(file = "./functions/compute_betadiversity.R")

### 18.3.2/ Compute betadiversity with a moving window ####

# Choose resolution
resolution <- "3"
resolution <- "5"
resolution <- "11"
resolution <- "17"

# Run function
moving_betadiversity <- compute_moving_betadiversity(sp_binary_stack, diversity_type = "taxo",
                                                     phylo = NULL, index_family = "sorensen",
                                                     beta_part = "turnover", aggreg_type = "all_pairwise",
                                                     resolution_factor = 1, 
                                                     window_width = as.numeric(resolution),
                                                     window_height = as.numeric(resolution), # In pixels, extended from the focal point
                                                     border_rule = NULL)  # Add an option on how to deal with the borders

moving_betadiversity_focal <- compute_moving_betadiversity(sp_binary_stack, diversity_type = "taxo",
                                                           phylo = NULL, index_family = "sorensen",
                                                           beta_part = "turnover", aggreg_type = "focal_pairwise",
                                                           resolution_factor = 1, 
                                                           window_width = as.numeric(resolution),
                                                           window_height = as.numeric(resolution), # In pixels, extended from the focal point
                                                           border_rule = NULL)  # Add an option on how to deal with the borders


moving_betaphylodiversity <- compute_moving_betadiversity(sp_binary_stack, diversity_type = "phylo",
                                                          phylo = phylogeny_object, index_family = "sorensen",
                                                          beta_part = "turnover", aggreg_type = "all_pairwise",
                                                          resolution_factor = 1, 
                                                          window_width = as.numeric(resolution),
                                                          window_height = as.numeric(resolution), # In pixels, extended from the focal point
                                                          border_rule = NULL)  # Add an option on how to deal with the borders

moving_betaphylodiversity_focal <- compute_moving_betadiversity(sp_binary_stack, diversity_type = "phylo",
                                                                phylo = phylogeny_object, index_family = "sorensen",
                                                                beta_part = "turnover", aggreg_type = "focal_pairwise",
                                                                resolution_factor = 1, 
                                                                window_width = as.numeric(resolution),
                                                                window_height = as.numeric(resolution), # In pixels, extended from the focal point
                                                                border_rule = NULL)  # Add an option on how to deal with the borders

# Plot raw maps
plot(moving_betadiversity)
plot(moving_betaphylodiversity)

plot(moving_betadiversity_focal)
plot(moving_betaphylodiversity_focal)

# Manual threshold
hist(moving_betadiversity[])
hist(moving_betaphylodiversity[])
hist(moving_betadiversity_focal[])
hist(moving_betaphylodiversity_focal[])

# Automatic threshold
threshold <- round(x = quantile(x = moving_betadiversity[], p = 0.999, na.rm = T) * 10, digits = 0)/10 ; threshold
threshold <- round(x = quantile(x = moving_betaphylodiversity[], p = 0.999, na.rm = T) * 10, digits = 0)/10  ; threshold
threshold <- round(x = quantile(x = moving_betadiversity_focal[], p = 0.999, na.rm = T) * 10, digits = 0)/10  ; threshold
threshold <- round(x = quantile(x = moving_betaphylodiversity_focal[], p = 0.999, na.rm = T) * 10, digits = 0)/10  ; threshold

# Apply threshold
moving_betadiversity_contrasted <- contrasting_raster(x = moving_betadiversity, zmin = 0, zmax = 0.6)
moving_betaphylodiversity_contrasted <- contrasting_raster(x = moving_betaphylodiversity, zmin = 0, zmax = 0.4)
moving_betadiversity_focal_contrasted <- contrasting_raster(x = moving_betadiversity_focal, zmin = 0, zmax = 0.6)
moving_betaphylodiversity_focal_contrasted <- contrasting_raster(x = moving_betaphylodiversity_focal, zmin = 0, zmax = 0.4)

# Plot contrasted maps
plot(moving_betadiversity_contrasted)
plot(moving_betaphylodiversity_contrasted)

plot(moving_betadiversity_focal_contrasted)
plot(moving_betaphylodiversity_focal_contrasted)

# Save any resolution
# save(moving_betadiversity, file = paste0("./outputs/Indices_maps/moving_betadiversity_",resolution,"x",resolution,".RData"))
# save(moving_betaphylodiversity, file = paste0("./outputs/Indices_maps/moving_betaphylodiversity_",resolution,"x",resolution,".RData"))
# save(moving_betadiversity_focal, file = paste0("./outputs/Indices_maps/moving_betadiversity_focal_",resolution,"x",resolution,".RData"))
# save(moving_betaphylodiversity_focal, file = paste0("./outputs/Indices_maps/moving_betaphylodiversity_focal_",resolution,"x",resolution,".RData"))

# Load any resolution
load(file = paste0("./outputs/Indices_maps/moving_betadiversity_",resolution,"x",resolution,".RData"))
load(file = paste0("./outputs/Indices_maps/moving_betaphylodiversity_",resolution,"x",resolution,".RData"))
load(file = paste0("./outputs/Indices_maps/moving_betadiversity_",resolution,"x",resolution,"_focal.RData"))
load(file = paste0("./outputs/Indices_maps/moving_betaphylodiversity_",resolution,"x",resolution,"_focal.RData"))


## Export to build GIF

# Export betadiversity all_pairwise
pdf(file = paste0("./maps/Indices_maps/Betadiversity/moving_betadiversity_",resolution,"x",resolution,".pdf"),
    width = 8, height = 8)
# Apply threshold
moving_betadiversity_contrasted <- contrasting_raster(x = moving_betadiversity, zmin = 0, zmax = 0.6)
# Plot contrasted maps
plot(moving_betadiversity_contrasted, main = paste0("Betadiversity \n(", resolution,"x",resolution,")"))
dev.off()

# Export betaphylodiversity all_pairwise
pdf(file = paste0("./maps/Indices_maps/Betadiversity/moving_betaphylodiversity_",resolution,"x",resolution,".pdf"),
    width = 8, height = 8)
# Apply threshold
moving_betaphylodiversity_contrasted <- contrasting_raster(x = moving_betaphylodiversity, zmin = 0, zmax = 0.4)
# Plot contrasted maps
plot(moving_betaphylodiversity_contrasted, main = paste0("Betaphlyodiversity \n(", resolution,"x",resolution,")"))
dev.off()

# Export betadiversity focal
pdf(file = paste0("./maps/Indices_maps/Betadiversity/moving_betadiversity_",resolution,"x",resolution,"_focal.pdf"),
    width = 8, height = 8)
# Apply threshold
moving_betadiversity_focal_contrasted <- contrasting_raster(x = moving_betadiversity_focal, zmin = 0, zmax = 0.6)
# Plot contrasted maps
plot(moving_betadiversity_focal_contrasted, main = paste0("Betadiversity focal \n(", resolution,"x",resolution,")"))
dev.off()

# Export betaphylodiversity focal
pdf(file = paste0("./maps/Indices_maps/Betadiversity/moving_betaphylodiversity_",resolution,"x",resolution,"_focal.pdf"),
    width = 8, height = 8)
# Apply threshold
moving_betaphylodiversity_focal_contrasted <- contrasting_raster(x = moving_betaphylodiversity_focal, zmin = 0, zmax = 0.4)
# Plot contrasted maps
plot(moving_betaphylodiversity_focal_contrasted, main = paste0("Betaphylodiversity focal \n(", resolution,"x",resolution,")"))
dev.off()



# #### Test for getValuesBlock border rules ####
# 
# test <- raster::crop(x = sp_proba_stack[[1:5]], y = extent(-80.5, -79.2, -1.2, 0))
# test_matrix <- matrix(data = 1:125, nrow = 25, ncol = 5)
# colnames(test_matrix) <- colnames(test@data@values)
# test@data@values <- test_matrix
#  
# test_stack <- stack(test)
# 
# plot(test_stack)
# 
# regional_data <- raster::getValuesBlock(x = test_stack, 
#                                         row = 2,
#                                         col = 2,
#                                         nrow = 2*2 + 1,
#                                         ncol = 1*2 + 1)
# 
# regional_data


### 18.4/ Betadiversity within regions ####

### To improve ?
# All types of inputs: spPolygons + sf object + raster Layer (or Stack) as a "mask"


### 18.4.1/ Function to compute any type of betadiversity from a sites x species matrix ####

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

source(file = "./functions/compute_betadiversity.R")


### 18.4.2/ Function to compute betadiversity within regions ####

# Compute betadiversity across all pixels/sites of each region, and eventually map them

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

source(file = "./functions/compute_betadiversity.R")


### 18.4.3/ Compute betadiversity within regions ####

# Load spPolygon shape files for bioregions as example
load(file = "./input_data/Map_stuff/Bioregions/All_bioregions_in_figure.RData")

# list_shp_regions <- list(Caatinga_shp, Guyana_Shield_shp)
# region_names <- c("Caatinga", "Guyana_Shield")

# list_shp_regions <- list(Central_Andes_shp, Western_Amazon_shp)
# region_names <- c("Central Andes", "Western Amazon")

plot(binary_stack)

list_shp_regions <- list(Caatinga_shp, Caribbean_Islands_shp, full_CA_shp, Central_Andes_shp, Cerrado_shp, Chacos_shp, Coastal_desert_shp, Guyana_Shield_shp, Llanos_shp, Lower_Amazon_shp, Mata_Atlantica_shp5, Northern_Andes_shp, Pampas_shp, Pantanal_shp, Western_Amazon_shp, Western_Lowlands_shp)
region_names <- c("Caatinga", "Caribbean_Islands", "Central_America", "Central_Andes", "Cerrado", "Chacos", "Coastal_desert", "Guyana_Shield", "Llanos", "Lower_Amazon", "Mata_Atlantica", "Northern_Andes", "Pampas", "Pantanal", "Western_Amazon", "Western_Lowlands")

# ### Morrone's regions (dominion)
# Morrone_dominions_shp_files <- readRDS(file = "./input_data/Map_stuff/Patricia_data/Morrone_dominions_shp_files.rds")
# Morrone_dominions_names <- readRDS(Patricia_region_names, file = "./input_data/Map_stuff/Patricia_data/Morrone_dominions_names.rds")
# ##

# Compute the indices and generate maps
within_regional_betadiversity <- compute_within_regional_betadiversity(sp_binary_stack, list_shp_regions, region_names,
                                                                       subsample_size = 10000,
                                                                       diversity_type = c("taxo", "phylo"),
                                                                       # diversity_type = c("taxo"),
                                                                       phylo = phylo, index_family = "sorensen",
                                                                       beta_part = c("turnover"),
                                                                       aggreg_type = c("mean_pairwise"),
                                                                       regions_map = T,
                                                                       output_type = "region_df",
                                                                       quiet = F)

Morrone_within_regional_betadiversity <- compute_within_regional_betadiversity(sp_binary_stack, Morrone_dominions_shp_files, Morrone_dominions_names,
                                                                       subsample_size = 10000,
                                                                       diversity_type = c("taxo", "phylo"),
                                                                       # diversity_type = c("taxo"),
                                                                       phylo = phylo, index_family = "sorensen",
                                                                       beta_part = c("turnover"),
                                                                       aggreg_type = c("mean_pairwise"),
                                                                       regions_map = T,
                                                                       output_type = "region_df",
                                                                       quiet = F)

View(within_regional_betadiversity[[1]])
plot(within_regional_betadiversity[[2]], zlim = c(0, 0.6))

# Save
save(within_regional_betadiversity, file = paste0("./outputs/Indices_maps/within_regional_betadiversity.RData"))
save(Morrone_within_regional_betadiversity, file = paste0("./outputs/Indices_maps/Morrone_within_regional_betadiversity.RData"))

# Compute differences per regions
load(file = paste0("./outputs/Indices_maps/within_regional_betadiversity.RData"))

diff_PBD_TBD_within_regions <- within_regional_betadiversity[[2]][[2]] - within_regional_betadiversity[[2]][[1]]
plot(diff_PBD_TBD_within_regions, zlim = c(-0.3, 0.3), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))

### 18.5/ Betadiversity between regions ####

# Compute pairwise betadiversity between regions aggregated as unique sites and provide associated pairwise betadiversity matrices

### 18.5.1/ Function to compute pairwise betadiversity matrix from a sites x species matrix ####

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

source(file = "./functions/compute_betadiversity.R")


### 18.5.2/ Function to compute betadiversity between regions #### 

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

source(file = "./functions/compute_betadiversity.R")

### 18.5.3/ Compute pairwise betadiversity between regions ####

pairwise_regional_betadiversity <- compute_pairwise_regional_betadiversity(sp_binary_stack, list_shp_regions, region_names,
                                                                           diversity_type = c("taxo", "phylo"),
                                                                           # diversity_type = c("taxo"),
                                                                           phylo = phylo, index_family = "sorensen",
                                                                           beta_part = c("turnover"))

pairwise_regional_betadiversity


# Save
save(pairwise_regional_betadiversity, file = paste0("./outputs/Indices_maps/pairwise_regional_betadiversity.RData"))

load(file = paste0("./outputs/Indices_maps/pairwise_regional_betadiversity.RData"))


### 18.6/ Map as an interaction graph the within and between regions diversities ####

##### Include a graph plot to visualize betadiversity between regions ####
# Background = map.
# Nods controid of regions (can use the within regions value to define size!)
# Links = pairwise diversity

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
    edge_weights <- igraph::E(network)$weight
    
    standardized_edges <- round(edge_weights/ max(edge_weights) * 200)
    
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
                         edge.width = standardized_edges / 200 * edge_cex
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

library(tidyverse)
library(ggraph)

nod_cex <- 30
edge_cex <- 1

# Load spPolygon shape files for bioregions as example
load(file = "./input_data/Map_stuff/Bioregions/All_bioregions_in_figure.RData")
list_shp_regions <- list(Caatinga_shp, Caribbean_Islands_shp, full_CA_shp, Central_Andes_shp, Cerrado_shp, Chacos_shp, Coastal_desert_shp, Guyana_Shield_shp, Llanos_shp, Lower_Amazon_shp, Mata_Atlantica_shp5, Northern_Andes_shp, Pampas_shp, Pantanal_shp, Western_Amazon_shp, Western_Lowlands_shp)
region_names <- c("Caatinga", "Caribbean_Islands", "Central_America", "Central_Andes", "Cerrado", "Chacos", "Coastal_desert", "Guyana_Shield", "Llanos", "Lower_Amazon", "Mata_Atlantica", "Northern_Andes", "Pampas", "Pantanal", "Western_Amazon", "Western_Lowlands")
region_labels <- c("Ca", "CIs", "CAm", "CAn", "Ce", "Ch", "CD", "GS", "Ll", "LA", "MA", "NAn", "Pam", "Pan", "UA", "WL")

# Load within betadiversity
load(file = paste0("./outputs/Indices_maps/within_regional_betadiversity.RData"))
region_df_within_taxo <- as.data.frame(within_regional_betadiversity[[1]])[, c(1,2)]
region_df_within_phylo <- as.data.frame(within_regional_betadiversity[[1]])[, c(1,3)]

# Load between betadiversity
load(file = paste0("./outputs/Indices_maps/pairwise_regional_betadiversity.RData"))
region_df_between_taxo <- pairwise_regional_betadiversity[[1]]
region_df_between_phylo <- pairwise_regional_betadiversity[[2]]

# Plot
map_betadiversity_within_and_between_regions (region_df_within = region_df_within_phylo,   # To provide directly the df with indices values within each region. 2 columns: the region names, the index values
                                              region_matrix_between = region_df_between_phylo,  # To provide directly the matrix of betadiversity values between regions
                                              list_shp_regions, # To provide the list of SpatialPolygons of regions to extract centroid coordinates
                                              region_labels = region_labels, # To provide labels to use instead of full region names
                                              nod_cex = 50,
                                              edge_cex = 2,
                                              plot_type = "ggplot")

map_betadiversity_within_and_between_regions (region_df_within = region_df_within_taxo,   # To provide directly the df with indices values within each region. 2 columns: the region names, the index values
                                              region_matrix_between = region_df_between_taxo,  # To provide directly the matrix of betadiversity values between regions
                                              list_shp_regions, # To provide the list of SpatialPolygons of regions to extract centroid coordinates
                                              region_labels = region_labels, # To provide labels to use instead of full region names
                                              nod_cex = 50,
                                              edge_cex = 2,
                                              plot_type = "ggplot")

map_betadiversity_within_and_between_regions (region_df_within = region_df_within_phylo,   # To provide directly the df with indices values within each region. 2 columns: the region names, the index values
                                              region_matrix_between = region_df_between_phylo,  # To provide directly the matrix of betadiversity values between regions
                                              list_shp_regions, # To provide the list of SpatialPolygons of regions to extract centroid coordinates
                                              region_labels = region_labels, # To provide labels to use instead of full region names
                                              nod_cex = 100,
                                              edge_cex = 20,
                                              plot_type = "igraph")


### 18.7/ Make the RGB visualization after NMDS in 3D for the full pairwise approach ####

### 18.7.1/ Function to map pairwise betadiversity with RGB scheme ####

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

source(file = "./functions/compute_betadiversity.R")

### 18.7.2/ Map pairwise betadiversity with RGB scheme ####

betadiversity_RGB_map <- map_pairwise_betadiversity(sp_binary_stack, diversity_type = "taxo",
                                                    phylo = NULL, index_family = "sorensen",
                                                    beta_part = "turnover",
                                                    color_space = "3D",
                                                    min_sp = 3,
                                                    subsample_size = 5000,
                                                    RGB_plot = T, NMDS_plot = T)

betaphylodiversity_RGB_map <- map_pairwise_betadiversity(sp_binary_stack, diversity_type = "phylo",
                                                         phylo = phylogeny_object, index_family = "sorensen",
                                                         beta_part = "turnover",
                                                         color_space = "3D",
                                                         min_sp = 3,
                                                         subsample_size = 5000,
                                                         RGB_plot = T, NMDS_plot = T)

plot(betadiversity_RGB_map)
plotRGB(betadiversity_RGB_map, axes = T, main = "RGB Plot", cex.main = 1.5, colNA = "aliceblue")

plot(betaphylodiversity_RGB_map)
plotRGB(betaphylodiversity_RGB_map, axes = T, main = "RGB Plot", cex.main = 1.5, colNA = "aliceblue")

# Save
save(betadiversity_RGB_map, file = paste0("./outputs/Indices_maps/betadiversity_RGB_map.RData"))
save(betaphylodiversity_RGB_map, file = paste0("./outputs/Indices_maps/betaphylodiversity_RGB_map.RData"))




### 19/ Ist and Pst ####

### 19.1/ Taxonomic and Phylogenetic partionning (Ist & Pst) ####

### 19.1.1/ Function to compute Ist and Pst with a moving window ###

# Size of the window, default as a proportion
# Multi vs. mean pairwise

### Can be improved by adding the option to compute for combination of parameters, all at the same time...


### See the functions I created for Erika's project!

##### Next step = mimicry Ist... but annoying because you have to recompute mimicry ring richness stack from OMU at each randomization...
### Mimicry Ist is based on mimicry richness stack which must be recomputed each time. 

compute_moving_diversity_partitioning <- function (proba_stack, diversity_type = "taxo",
                                                   phylo = NULL, aggreg_type = "multi",
                                                   resolution_factor = 1,  # Provide a factor to aggregate data if needed (could be necessary to reduce computation time for wide window size)
                                                   window_width = "default", window_height = "default", # In pixels, extended from the focal point
                                                   border_rule = NULL)  # Add an option on how to deal with the borders
{
  cat(paste0(Sys.time(), " - Index Computation - Start\n"))
  
  ### Aggregate data if needed
  
  if (resolution_factor != 1)
  {
    # Function to compute probability of presence of at least one species/OMUs
    aggreg_prob <- function(x, na.rm) 
    { 
      y <- 1 - prod(1 - x) # Probability of presence of ring = probability of presence of at least one species/OMU = opposite of probability of absence of all species/OMU
      return(y) # Output
    }
    
    proba_stack <- raster::aggregate(x = proba_stack,
                                     fact = resolution_factor,
                                     fun = aggreg_prob)
  }
  
  ### Adjust default window size to a fifth of the (aggregated) raster width and height
  
  if (window_width == "default")
  {
    window_width <- round(proba_stack@ncols/10)
    windows_width <- windows_width * 2 + 1
    if (window_width == 0) {window_width <- 3}
  }
  if (window_height == "default")
  {
    window_height <- round(proba_stack@nrows/10)
    window_height <- windows_width * 2 + 1
    if (window_height == 0) {window_height <- 3}
  }
  
  ### If provided, check that the window size parameters are odd integers
  if (((window_width + 1) %% 2) != 0) 
  {
    stop("Window width must be an odd integer")
  }
  if (((window_height + 1) %% 2) != 0) 
  {
    stop("Window height must be an odd integer")
  }  
  
  # If using the phylogenetic beta-diversity, need to match raster Stack and phylogeny species lists
  if (diversity_type == "phylo")
  {
    # Match raster Stack and Phylogeny species lists
    clean_stack_and_phylo <- match_stack_and_phylo(proba_stack, phylo)
    proba_stack <- clean_stack_and_phylo[[1]]
    pruned_tree <- clean_stack_and_phylo[[2]]
  }
  
  # Extract community matrix
  proba_mat <- raster::getValues(proba_stack)
  
  
  ###### IN PROGRESS #######
  
  ### Loop to move from pixel to pixel
  all_index_values <- NA
  for (i in 1:nrow(proba_stack))
  {
    # i <- 1
    
    # Initiate control for focal data and regional data
    no_focal_data <- F
    no_regional_data <- F
    
    # Extract cell location
    RowCol <- raster::rowColFromCell(object = proba_stack, cell = i)
    
    ### Warning: Window size is reduced for right (add) and bottom borders (cutoff), but not for left and upper borders (shift window)
    # Add an option on how to deal with the borders with border_rule
    
    # Extract regional data
    regional_data <- raster::getValuesBlock(x = binary_stack, 
                                            row = RowCol[1] - ((window_width - 1) / 2),
                                            col = RowCol[2] - ((window_height - 1) / 2),
                                            nrow = window_width,
                                            ncol = window_height)
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
    all_index_values[i] <- beta_value
    
    # Show i every 100 iterations and save a back-up file
    if (i %% 100 == 0) 
    {
      cat(paste0(Sys.time(), " - ", i," on ",nrow(proba_stack),"\n"))
      # save(all_FP, file = "./outputs/Indices_maps/backup_betadiversity_moving.RData")
    }
  }
  
  # Generate a mask for terrestrial areas
  continental_mask <- (calc(proba_stack, fun = sum) >= 0) - 1
  
  # Write beta-diversity values in a raster without continental borders
  betadiversity_raster_no_borders <- proba_stack[[1]]
  betadiversity_raster_no_borders@data@values <- all_index_values
  
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

source(file = "./functions/compute_betadiversity.R")

### 19.1.2/ Compute Ist and Pst with a moving window ###

### 19.2/ Mimicry turnover with mimicry Ist ####

### 19.2.1/ Function to compute mimicry Ist ###

### 19.2.2/ Compute mimicry Ist ###



##### 20/ Map residuals between two indices #####

# Input = two index rasters
# Output = residuals Y ~ X
# Need to perform some sort of models to get the equation for the relationship
# LM on every pixels (Gumbs et al., 2020)
# GLS with spatial covariance structure (Devictor et al., 2010)
# GAM on every pixels
# Extract R² and put it on the map, with the equation
# Need to subsample to avoid spatial autocorrelation issues: N = 1000 or 1/10th of data
# Random (multiple time) or regular ?

### 20.1/ Function to compute and map residuals between two indices ####

library(tidyverse)

map_residuals <- function (y_index, x_index, 
                           method = "GAM", # Either LM for linear relationship, LM_quadratic to add quadratic effects, or GLS to account for spatial autocorrelation, or GAM to allow non-linear relationship.
                           gamma = 1,  # Smoothing parameter for GAM models. Higher than 1 = increase smoothness of the fit
                           corSpatial = "corExp", # To choose the type of spatial autocorrelation structure as in nlme::corSpatial()
                           subsample_size = NA,  # Number of sites to subsample
                           subsample_runs = 10,  # Number of random sample runs to aggregate
                           include_model = F, # To include the fitted model in the output
                           scatterplot = T,  # To display the scatter plot with regression line and R²
                           plot_variogram_correlogram = F) # To display the variogram and correlogram in case of spatial GLS
  
{
  y_values <- getValues(y_index)
  x_values <- getValues(x_index)
  indices_df_full <- data.frame(y = y_values, x = x_values)
  indices_df <- indices_df_full[apply(X = indices_df_full, MARGIN = 1, FUN = sum, na.rm = T) > 0, ] # Keep only non-null communities
  
  # See what is in common between methods and move it after, outside if
  
  if (method == "LM") 
  {
    model <- lm(data = indices_df, formula = y ~ x)
    
    # Extract R²
    R2 <- round(summary(model)$adj.r.squared, 3)
    
    # Print output
    cat(paste0("\nLM for ", names(y_index), " ~ ", names(x_index), ":\n\n"))
    summary(model)
  }
  
  if (method == "LM_quadratic") 
  {
    model <- lm(data = indices_df, formula = y ~ x + I(x^2))
    
    # Extract R²
    R2 <- round(summary(model)$adj.r.squared, 3)
    
    # Print output
    cat(paste0("\nLM with quadratic effects for ", names(y_index), " ~ ", names(x_index), ":\n\n"))
    summary(model)
  }
  
  if (method == "GAM")
  {
    model <- mgcv::gam(data = indices_df, formula = y ~ s(x), gamma = gamma)
    
    # Extract R²
    R2 <- round(summary(model)$r.sq, 3)
    
    # Print output
    cat(paste0("\nGAM for ", names(y_index), " ~ ", names(x_index), ":\n\n"))
    summary(model)
  }
  
  if (method == "GLS")
  {
    # ?nlme::gls
    # ?nlme::corSpatial
    # 
    # ?nlme::corExp # Exponential spatial correlation structure
    # ?nlme::corLin # Linear spatial correlation structure
    # ?nlme::corGaus # Gaussian spatial correlation structure = Sinusoidal shape
    # ?nlme::corSpher # Spherical spatial correlation structure = Exponential structure with distance threshold for plateau
    # ?nlme::corRatio # Rational quadratic spatial correlation structure = Quadratic shape
    
    # Extract coordinates
    coords <- raster::coordinates(y_index)
    
    set.seed(seed = 1)
    
    ### Add a loop to do this N times and compute the mean coefficients to avois issue of random sampling ###
    
    ### Extract the sampled sites
    
    sample_indices <- sample(x = 1:nrow(indices_df), size = subsample_size)
    indices_df <- indices_df[sample_indices, ]
    sampled_coords <- coords[sample_indices, ]
    
    ### Compute correlation matrix depending on the type of Spatial correlation structure chosen
    if(corSpatial == "corExp") { corr_matrix <- nlme::corExp(form = ~ sampled_coords[, 1] + sampled_coords[, 2], nugget = TRUE) }
    if(corSpatial == "corLin") { corr_matrix <- nlme::corLin(form = ~ sampled_coords[, 1] + sampled_coords[, 2], nugget = TRUE) }
    if(corSpatial == "corGaus") { corr_matrix <- nlme::corGaus(form = ~ sampled_coords[, 1] + sampled_coords[, 2], nugget = TRUE) }
    if(corSpatial == "corSpher") { corr_matrix <- nlme::corGaus(form = ~ sampled_coords[, 1] + sampled_coords[, 2], nugget = TRUE) }
    if(corSpatial == "corRatio") { corr_matrix <- nlme::corRatio(form = ~ sampled_coords[, 1] + sampled_coords[, 2], nugget = TRUE) }
    
    ### Compute GLS models
    gls_null_model <- nlme::gls(data = indices_df, model = y ~ x + I(x^2))
    model <- nlme::gls(data = indices_df, model = y ~ x + I(x^2), correlation = corr_matrix)
    
    # Extract residuals for null model and true model
    resids_null <- residuals(gls_null_model)
    resids <- residuals(model)
    
    # Extract R²
    R2 <- round(performance::r2(model)[[1]], 3)
    
    # Print output
    cat(paste0("\nGLS for ", names(y_index), " ~ ", names(x_index), " with spatial autocorrelation structure:\n\n"))
    summary(model)
    
    ### Plot variogram and correlogram to check fit of the model (only one time)
    if (plot_variogram_correlogram)
    {
      ### Find a way to get the variogram on a plot grid! Recreate them from output ? ###
      
      par(mfrow = c(2,2))
      
      # Check the fit of the model to the semi-variogram for Pearson residuals
      gls_variogram <- nlme::Variogram(model, form = ~ sampled_coords[, 1] + sampled_coords[, 2], resType = "pearson")
      plot(gls_variogram, smooth = TRUE, ylim = c(0, 1.2), main = "Semi-variogram for Pearson residuals")
      
      # Check for the absence of trends on the semi-variogram for normalized residuals
      gls_variogram_norm <- nlme::Variogram(model, form = ~ sampled_coords[, 1] + sampled_coords[, 2], resType = "normalized")
      plot(gls_variogram_norm, smooth = TRUE, ylim = c(0, 1.2), main = "Semi-variogram for Normalized residuals")
      
      # Define bins size
      sp_dist <- geosphere::distm(sampled_coords, fun = geosphere::distHaversine)/1000
      bins <- round(x= 1 + 3.3*log(x = max(sp_dist), base = 2), digits = 0) # Struges's rule
      # bins <- 20
      increment <- round(max(sp_dist)/bins,0)
      
      # ?ncf::correlog
      
      # Generate correlogram of residuals for the GLS model
      gls_correlogram <- ncf::correlog(x = sampled_coords[, 1], y = sampled_coords[, 2], z = resids_null, latlon = T, increment = increment, resamp = 100, na.rm = T, quiet = T)
      # Plot correlogram for GLS with spatial autocorrelation structure
      plot(gls_correlogram$mean.of.class[1:10], gls_correlogram$correlation[1:10], type = "b", col = c("black","red")[(gls_correlogram$p > 0.05) + 1], pch = 16, cex = 1, lwd = 1.5, main = "Correlogram for GLS model", xlab = "Distance (km)", ylab = "Moran's I", cex.lab = 1.2, cex.axis = 1) ; abline(h = 0)
      
      # Generate correlogram of residuals for the null model
      gls_correlogram_null <- ncf::correlog(x = sampled_coords[, 1], y = sampled_coords[, 2], z = resids, latlon = T, increment = increment, resamp = 100, na.rm = T, quiet = T)
      # Plot correlogram for GLS without spatial autocorrelation structure
      plot(gls_correlogram_null$mean.of.class[1:10], gls_correlogram_null$correlation[1:10], type = "b", col = c("black","red")[(gls_correlogram_null$p > 0.05) + 1], pch = 16, cex = 1, lwd = 1.5, main = "Correlogram for Null model", xlab = "Distance (km)", ylab = "Moran's I", cex.lab = 1.2, cex.axis = 1) ; abline(h=0)
      
      par(mfrow = c(1,1))
    }
  }
  
  # Plot indices values with model predictions and GOF
  if (scatterplot)
  {
    # indices_df <- indices_df[sample(x = 1:nrow(indices_df), size = 100), ]
    
    # Extract residuals for the color scheme
    resids <- residuals(model)
    indices_df <- indices_df %>% 
      mutate(Residuals = resids)
    
    # Compute predict for the regression line
    predicts <- predict(object = model, newdata = data.frame(x = seq(from = min(x_values, na.rm = T), to = max(x_values, na.rm = T), length.out = 1000)))
    predicts_df <- data.frame(x = seq(from = min(x_values, na.rm = T), to = max(x_values, na.rm = T), length.out = 1000), y = predicts)
    
    # Extract AICc
    AICc <- round(MuMIn::AICc(model), 1)
    
    # Plot
    g <- ggplot(data = indices_df) +
      geom_point(data = indices_df, aes(x = x, y = y, color = Residuals), alpha = 0.7, show.legend = T) +
      geom_line(data = predicts_df, aes(x = x, y = y), size = 1.5, col = "red", linetype = "solid") +
      
      # Set limits to display origin
      xlim(0, max(x_values, na.rm = T)) +
      ylim(0, max(y_values, na.rm = T)) +
      
      # Display R²
      annotate("text", x = 0.12*max(x_values, na.rm = T), y = 0.95*max(y_values, na.rm = T), label = paste0("R² = ", R2), size = 5, fontface = 2) +
      
      # Display AICc
      annotate("text", x = 0.12*max(x_values, na.rm = T), y = 0.88*max(y_values, na.rm = T), label = paste0("AICc = ", AICc), size = 5, fontface = 2) +
      
      # Add legend for predict
      geom_segment(x = 0.00*max(x_values, na.rm = T), y = 0.81*max(y_values, na.rm = T), xend = 0.10*max(x_values, na.rm = T), yend = 0.81*max(y_values, na.rm = T), color = "#FF000080", alpha = 0.7, size = 1.5) +
      annotate("text", x = 0.25*max(x_values, na.rm = T), y = 0.81*max(y_values, na.rm = T), label = paste0(method, " predict"), size = 5, fontface = 2) +
      
      # Set axis labels
      ylab(label = names(y_index)) +
      xlab(label = names(x_index)) +
      
      # Set title
      ggtitle(paste0(names(y_index), " ~ ", names(x_index))) +
      
      # Arrange axes limits
      ggthemes::geom_rangeframe(data = indices_df, aes(x = x, y = y), size = 1.4) +
      
      # Set color gradient
      scale_color_gradient2(low = scales::muted("blue"), high = scales::muted("red")) +
      
      # Adjust 
      theme(panel.background = element_rect(fill = "white"),
            plot.title = element_text(face = 2, size = 18, hjust = 0.5, margin = margin(t = 5, b = 10)),
            
            legend.position = c(0.85, 0.3),
            legend.background = element_rect(fill = NA, colour = NA),
            legend.title = element_text(size = 14, vjust = 3, face = "bold"),
            legend.text = element_text(size = 12, face = "bold"),
            # # margin = margin(t = 0, unit = "pt")),
            # # legend.key.size = unit(1.5, 'lines'),
            # # legend.spacing.y = unit(1,"cm"),
            
            axis.line = element_line(size = 1.2),
            axis.ticks = element_line(size = 1.2),
            axis.ticks.length = unit(10, "pt"),
            axis.text = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            axis.title.x = element_text(margin = margin(t = 10, b = 5)),
            axis.title.y = element_text(margin = margin(l = 5, r = 10)))
    
    print(g)
  }
  
  ### Predict values for all cells
  predicts_values <- predict(object = model, newdata = indices_df_full)
  
  ### Put predicts in a Raster Layer
  predicts_raster <- y_index
  predicts_raster@data@values <- c(predicts_values)
  
  ### Compute residuals
  residuals_raster <- y_index - predicts_raster
  
  ### Export
  if (include_model) 
  { 
    return(list(residuals_raster, model)) 
  } else {
    return(residuals_raster)
  }
}

source(file = "./functions/compute_residuals.R")

### 20.2/ Compute residuals ####

### Should use SR only from species in the phylogeny to have the same species pool than PD
sp_richness_phylo_only <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_Jaccard.80_phylo.rds"))
sp_richness <- readRDS(file = paste0("./outputs/Indices_maps/sp_richness.rds"))
Faith_PD <- readRDS(file = paste0("./outputs/Indices_maps/Faith_PD.rds"))

plot(sp_richness)
plot(sp_richness_phylo_only)
plot(Faith_PD)

PD_SR_residuals <- map_residuals(y_index = Faith_PD, x_index = sp_richness_phylo_only, 
                                 method = "GAM", # Either LM for linear relationship, LM_quadratic to add quadratic effects, or GLS to account for spatial autocorrelation, or GAM to allow non-linear relationship.
                                 corSpatial = "corExp", # To choose the type of spatial autocorrelation structure as in nlme::corSpatial()
                                 subsample_size = NA,  # Number of sites to subsample
                                 subsample_runs = 10,  # Number of random sample runs to aggregate
                                 include_model = F, # To include the fitted model in the output
                                 scatterplot = T,  # To display the scatter plot with regression line and R²
                                 plot_variogram_correlogram = F) # To display the variogram and correlogram in case of spatial GLS


sp_richness_phylo_only <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_Jaccard.80_phylo.rds"))
sp_richness <- readRDS(file = paste0("./outputs/Indices_maps/sp_richness.rds"))
MPD <- readRDS(file = paste0("./outputs/Indices_maps/MPD.raster_Jaccard.80.rds"))

synchro_stack <- virtualspecies::synchroniseNA(stack(sp_richness_phylo_only, MPD))
sp_richness_phylo_only <- synchro_stack[[1]]
MPD <- synchro_stack[[2]]

names(sp_richness_phylo_only) <- "sp_richness"
names(MPD) <- "MPD"

MPD_SR_residuals <- map_residuals(y_index = MPD, x_index = sp_richness_phylo_only, 
                                  method = "GAM", # Either LM for linear relationship, LM_quadratic to add quadratic effects, or GLS to account for spatial autocorrelation, or GAM to allow non-linear relationship.
                                  gamma = 1, # Smoothing parameter for GAM models. Higher than 1 = increase smoothness of the fit
                                  corSpatial = "corExp", # To choose the type of spatial autocorrelation structure as in nlme::corSpatial()
                                  subsample_size = NA,  # Number of sites to subsample
                                  subsample_runs = 10,  # Number of random sample runs to aggregate
                                  include_model = F, # To include the fitted model in the output
                                  scatterplot = T,  # To display the scatter plot with regression line and R²
                                  plot_variogram_correlogram = F) # To display the variogram and correlogram in case of spatial GLS


sp_richness_phylo_only <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_Jaccard.80_phylo.rds"))
sp_richness <- readRDS(file = paste0("./outputs/Indices_maps/sp_richness.rds"))
load(file = paste0("./outputs/Indices_maps/ring_richness.RData"))

names(sp_richness_phylo_only) <- "sp_richness"
names(ring_richness) <- "ring_richness"

MR_SR_residuals <- map_residuals(y_index = ring_richness, x_index = sp_richness_phylo_only, 
                                 method = "GAM", # Either LM for linear relationship, LM_quadratic to add quadratic effects, or GLS to account for spatial autocorrelation, or GAM to allow non-linear relationship.
                                 corSpatial = "corExp", # To choose the type of spatial autocorrelation structure as in nlme::corSpatial()
                                 subsample_size = NA,  # Number of sites to subsample
                                 subsample_runs = 10,  # Number of random sample runs to aggregate
                                 include_model = F, # To include the fitted model in the output
                                 scatterplot = T,  # To display the scatter plot with regression line and R²
                                 plot_variogram_correlogram = F) # To display the variogram and correlogram in case of spatial GLS



load(file = paste0("./outputs/Indices_maps/moving_betadiversity.RData"))
load(file = paste0("./outputs/Indices_maps/moving_betaphylodiversity.RData"))
load(file = paste0("./outputs/Indices_maps/moving_betadiversity_5x5.RData"))
load(file = paste0("./outputs/Indices_maps/moving_betaphylodiversity_5x5.RData"))
load(file = paste0("./outputs/Indices_maps/moving_betadiversity_11x11.RData"))
load(file = paste0("./outputs/Indices_maps/moving_betaphylodiversity_11x11.RData"))
load(file = paste0("./outputs/Indices_maps/moving_betadiversity_17x17.RData"))
load(file = paste0("./outputs/Indices_maps/moving_betaphylodiversity_17x17.RData"))

TBD_PBD_residuals <- map_residuals(y_index = moving_betaphylodiversity, x_index = moving_betadiversity, 
                                   method = "GAM", # Either LM for linear relationship, LM_quadratic to add quadratic effects, or GLS to account for spatial autocorrelation, or GAM to allow non-linear relationship.
                                   corSpatial = "corExp", # To choose the type of spatial autocorrelation structure as in nlme::corSpatial()
                                   subsample_size = NA,  # Number of sites to subsample
                                   subsample_runs = 10,  # Number of random sample runs to aggregate
                                   include_model = F, # To include the fitted model in the output
                                   scatterplot = T,  # To display the scatter plot with regression line and R²
                                   plot_variogram_correlogram = F) # To display the variogram and correlogram in case of spatial GLS


# Plot
plot(PD_SR_residuals, main = "PD ~ SR residuals", col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
plot(MPD_SR_residuals, main = "MPD ~ SR residuals", col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200), zlim = c(-15, 15))
plot(MR_SR_residuals, main = "MR ~ SR residuals", col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
plot(TBD_PBD_residuals, main = "PBD ~ TBD residuals", col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200), zlim = c(-0.15, 0.15))


# Save
save(PD_SR_residuals, file = paste0("./outputs/Indices_maps/PD_SR_residuals.RData"))
save(MPD_SR_residuals, file = paste0("./outputs/Indices_maps/MPD_SR_residuals.RData"))
save(MR_SR_residuals, file = paste0("./outputs/Indices_maps/MR_SR_residuals.RData"))
save(TBD_PBD_residuals, file = paste0("./outputs/Indices_maps/TBD_PBD_residuals.RData"))

load(file = paste0("./outputs/Indices_maps/PD_SR_residuals.RData"))
load(file = paste0("./outputs/Indices_maps/MPD_SR_residuals.RData"))
load(file = paste0("./outputs/Indices_maps/MR_SR_residuals.RData"))
load(file = paste0("./outputs/Indices_maps/TBD_PBD_residuals.RData"))

# ### Check for the trends in the regression plot of TBD ~ PBD
# test_red <- (moving_betaphylodiversity > 0.25) & (moving_betadiversity > 0.45)
# plot(test_red)
# 
# test_blue <- (moving_betaphylodiversity < 0.25) & (moving_betadiversity > 0.45)
# plot(test_blue)

### 20.3/ Residuals within regions ####

# 20.3.1/ Function to compute residuals within regions

map_residuals_within_regions <- function (y_index, x_index, 
                                          list_shp_regions, # Provide region Spatial Polygon to define borders
                                          method = "GAM", # Either LM for linear relationship, LM_quadratic to add quadratic effects, or GLS to account for spatial autocorrelation, or GAM to allow non-linear relationship.
                                          corSpatial = "corExp", # To choose the type of spatial autocorrelation structure as in nlme::corSpatial() for GLS
                                          include_model = F, # To include the fitted model in the output
                                          scatterplot = T,  # To display the scatter plot with regression line and R²
                                          plot_variogram_correlogram = F) # To display the variogram and correlogram in case of spatial GLS
  
{
  y_values <- sapply(X = list_shp_regions, FUN = raster::extract, x = y_index)
  y_values <- lapply(X = y_values, FUN = unlist) # Step useful if multiple polygons for one region
  y_values_regions <- round(unlist(lapply(X = y_values, FUN = median, na.rm = T)),3)
  
  x_values <- sapply(X = list_shp_regions, FUN = raster::extract, x = x_index)
  x_values <- lapply(X = x_values, FUN = unlist) # Step useful if multiple polygons for one region
  x_values_regions <- round(unlist(lapply(X = x_values, FUN = median, na.rm = T)),3)
  
  indices_df <- data.frame(y = y_values_regions, x = x_values_regions)
  
  # See what is in common between methods and move it after, outside if
  
  if (method == "LM") 
  {
    model <- lm(data = indices_df, formula = y ~ x)
    
    # Extract R²
    R2 <- round(summary(model)$adj.r.squared, 3)
    
    # Print output
    cat(paste0("\nLM for ", names(y_index), " ~ ", names(x_index), ":\n\n"))
    summary(model)
  }
  
  if (method == "LM_quadratic") 
  {
    model <- lm(data = indices_df, formula = y ~ x + I(x^2))
    
    # Extract R²
    R2 <- round(summary(model)$adj.r.squared, 3)
    
    # Print output
    cat(paste0("\nLM with quadratic effects for ", names(y_index), " ~ ", names(x_index), ":\n\n"))
    summary(model)
  }
  
  if (method == "GAM")
  {
    model <- mgcv::gam(data = indices_df, formula = y ~ s(x))
    
    # Extract R²
    R2 <- round(summary(model)$r.sq, 3)
    
    # Print output
    cat(paste0("\nGAM for ", names(y_index), " ~ ", names(x_index), ":\n\n"))
    summary(model)
  }
  
  if (method == "GLS")
  {
    # ?nlme::gls
    # ?nlme::corSpatial
    # 
    # ?nlme::corExp # Exponential spatial correlation structure
    # ?nlme::corLin # Linear spatial correlation structure
    # ?nlme::corGaus # Gaussian spatial correlation structure = Sinusoidal shape
    # ?nlme::corSpher # Spherical spatial correlation structure = Exponential structure with distance threshold for plateau
    # ?nlme::corRatio # Rational quadratic spatial correlation structure = Quadratic shape
    
    ### Extract coordinates of regions centroids
    centroids_coords <- lapply(X = list_shp_regions, FUN = rgeos::gCentroid, byid=FALSE)
    centroids_coords <- lapply(X = centroids_coords, FUN = function (x) { x@coords })
    centroids_coords <- do.call(what = "rbind", args = centroids_coords)
    row.names(centroids_coords) <- region_names
    
    # # Version when regions all have a single polygon
    # all_shp_polygons <- do.call(raster::bind, list_shp_regions)
    # centroids_coords <- rgeos::gCentroid(spgeom = all_shp_polygons, byid=TRUE)@coords
    
    set.seed(seed = 1)
    
    ### Compute correlation matrix depending on the type of Spatial correlation structure chosen
    if(corSpatial == "corExp") { corr_matrix <- nlme::corExp(form = ~ centroids_coords[, 1] + centroids_coords[, 2], nugget = TRUE) }
    if(corSpatial == "corLin") { corr_matrix <- nlme::corLin(form = ~ centroids_coords[, 1] + centroids_coords[, 2], nugget = TRUE) }
    if(corSpatial == "corGaus") { corr_matrix <- nlme::corGaus(form = ~ centroids_coords[, 1] + centroids_coords[, 2], nugget = TRUE) }
    if(corSpatial == "corSpher") { corr_matrix <- nlme::corGaus(form = ~ centroids_coords[, 1] + centroids_coords[, 2], nugget = TRUE) }
    if(corSpatial == "corRatio") { corr_matrix <- nlme::corRatio(form = ~ centroids_coords[, 1] + centroids_coords[, 2], nugget = TRUE) }
    
    ### Compute GLS models
    gls_null_model <- nlme::gls(data = indices_df, model = y ~ x + I(x^2))
    model <- nlme::gls(data = indices_df, model = y ~ x + I(x^2), correlation = corr_matrix)
    
    # Extract residuals for null model and true model
    resids_null <- residuals(gls_null_model)
    resids <- residuals(model)
    
    # Extract R²
    R2 <- round(performance::r2(model)[[1]], 3)
    
    # Print output
    cat(paste0("\nGLS for ", names(y_index), " ~ ", names(x_index), " with spatial autocorrelation structure:\n\n"))
    summary(model)
    
    ### Plot variogram and correlogram to check fit of the model (only one time)
    if (plot_variogram_correlogram)
    {
      ### Find a way to get the variogram on a plot grid! Recreate them from output ? ###
      
      par(mfrow = c(2,2))
      
      # Check the fit of the model to the semi-variogram for Pearson residuals
      gls_variogram <- nlme::Variogram(model, form = ~ centroids_coords[, 1] + centroids_coords[, 2], resType = "pearson")
      plot(gls_variogram, smooth = TRUE, ylim = c(0, 1.2), main = "Semi-variogram for Pearson residuals")
      
      # Check for the absence of trends on the semi-variogram for normalized residuals
      gls_variogram_norm <- nlme::Variogram(model, form = ~ centroids_coords[, 1] + centroids_coords[, 2], resType = "normalized")
      plot(gls_variogram_norm, smooth = TRUE, ylim = c(0, 1.2), main = "Semi-variogram for Normalized residuals")
      
      # Define bins size
      sp_dist <- geosphere::distm(centroids_coords, fun = geosphere::distHaversine)/1000
      bins <- round(x= 1 + 3.3*log(x = max(sp_dist), base = 2), digits = 0) # Struges's rule
      # bins <- 20
      increment <- round(max(sp_dist)/bins,0)
      
      # ?ncf::correlog
      
      # Generate correlogram of residuals for the GLS model
      gls_correlogram <- ncf::correlog(x = centroids_coords[, 1], y = centroids_coords[, 2], z = resids_null, latlon = T, increment = increment, resamp = 100, na.rm = T, quiet = T)
      # Plot correlogram for GLS with spatial autocorrelation structure
      plot(gls_correlogram$mean.of.class[1:10], gls_correlogram$correlation[1:10], type = "b", col = c("black","red")[(gls_correlogram$p > 0.05) + 1], pch = 16, cex = 1, lwd = 1.5, main = "Correlogram for GLS model", xlab = "Distance (km)", ylab = "Moran's I", cex.lab = 1.2, cex.axis = 1) ; abline(h = 0)
      
      # Generate correlogram of residuals for the null model
      gls_correlogram_null <- ncf::correlog(x = centroids_coords[, 1], y = centroids_coords[, 2], z = resids, latlon = T, increment = increment, resamp = 100, na.rm = T, quiet = T)
      # Plot correlogram for GLS without spatial autocorrelation structure
      plot(gls_correlogram_null$mean.of.class[1:10], gls_correlogram_null$correlation[1:10], type = "b", col = c("black","red")[(gls_correlogram_null$p > 0.05) + 1], pch = 16, cex = 1, lwd = 1.5, main = "Correlogram for Null model", xlab = "Distance (km)", ylab = "Moran's I", cex.lab = 1.2, cex.axis = 1) ; abline(h=0)
      
      par(mfrow = c(1,1))
    }
  }
  
  # Plot indices values with model predictions and GOF
  if (scatterplot)
  {
    # indices_df <- indices_df[sample(x = 1:nrow(indices_df), size = 100), ]
    
    # Extract residuals for the color scheme
    resids <- residuals(model)
    indices_df <- indices_df %>% 
      mutate(Residuals = resids)
    
    # Compute predict for the regression line
    predicts <- predict(object = model, newdata = data.frame(x = seq(from = min(x_values_regions, na.rm = T), to = max(x_values_regions, na.rm = T), length.out = 1000)))
    predicts_df <- data.frame(x = seq(from = min(x_values_regions, na.rm = T), to = max(x_values_regions, na.rm = T), length.out = 1000), y = predicts)
    
    # Extract AICc
    AICc <- round(MuMIn::AICc(model), 1)
    
    # Plot
    g <- ggplot(data = indices_df) +
      geom_point(data = indices_df, aes(x = x, y = y, color = Residuals), alpha = 0.7, show.legend = T) +
      geom_line(data = predicts_df, aes(x = x, y = y), size = 1.5, col = "red", linetype = "solid") +
      
      # Set limits to display origin
      xlim(0, max(x_values_regions, na.rm = T)) +
      ylim(0, max(y_values_regions, na.rm = T)) +
      
      # Display R²
      annotate("text", x = 0.12*max(x_values_regions, na.rm = T), y = 0.95*max(y_values_regions, na.rm = T), label = paste0("R² = ", R2), size = 5, fontface = 2) +
      
      # Display AICc
      annotate("text", x = 0.12*max(x_values_regions, na.rm = T), y = 0.88*max(y_values_regions, na.rm = T), label = paste0("AICc = ", AICc), size = 5, fontface = 2) +
      
      # Add legend for predict
      geom_segment(x = 0.00*max(x_values_regions, na.rm = T), y = 0.81*max(y_values_regions, na.rm = T), xend = 0.10*max(x_values_regions, na.rm = T), yend = 0.81*max(y_values_regions, na.rm = T), color = "#FF000080", alpha = 0.7, size = 1.5) +
      annotate("text", x = 0.25*max(x_values_regions, na.rm = T), y = 0.81*max(y_values_regions, na.rm = T), label = paste0(method, " predict"), size = 5, fontface = 2) +
      
      # Set axis labels
      ylab(label = names(y_index)) +
      xlab(label = names(x_index)) +
      
      # Set title
      ggtitle(paste0(names(y_index), " ~ ", names(x_index))) +
      
      # Arrange axes limits
      ggthemes::geom_rangeframe(data = indices_df, aes(x = x, y = y), size = 1.4) +
      
      # Set color gradient
      scale_color_gradient2(low = scales::muted("blue"), high = scales::muted("red")) +
      
      # Adjust 
      theme(panel.background = element_rect(fill = "white"),
            plot.title = element_text(face = 2, size = 18, hjust = 0.5, margin = margin(t = 5, b = 10)),
            
            legend.position = c(0.85, 0.3),
            legend.background = element_rect(fill = NA, colour = NA),
            legend.title = element_text(size = 14, vjust = 3, face = "bold"),
            legend.text = element_text(size = 12, face = "bold"),
            # # margin = margin(t = 0, unit = "pt")),
            # # legend.key.size = unit(1.5, 'lines'),
            # # legend.spacing.y = unit(1,"cm"),
            
            axis.line = element_line(size = 1.2),
            axis.ticks = element_line(size = 1.2),
            axis.ticks.length = unit(10, "pt"),
            axis.text = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            axis.title.x = element_text(margin = margin(t = 10, b = 5)),
            axis.title.y = element_text(margin = margin(l = 5, r = 10)))
    
    print(g)
  }
  
  ### Predict values for all regions
  predicts_values <- predict(object = model, newdata = indices_df)
  
  ### Compute residuals per regions
  residuals_regions <- y_values_regions - predicts_values
  
  ### Turn residuals into a Raster Layer
  
  # Generate a mask for terrestrial areas to use as background
  continental_mask <- (y_index < Inf) - 1
  
  # Initiate stack per region
  region_stack <- stack()
  # Loop per region
  for (i in 1:length(list_shp_regions))
  {
    # i <- 1
    
    region_stack <- addLayer(region_stack, rasterize(x = list_shp_regions[[i]], 
                                                     y = continental_mask, # Provide the grid to fill with CRS, bbox and resolution
                                                     field = residuals_regions[i], # How to fill non empty cells. With the value of a variable in the df of the sp_obj, or directly with a fixed value ?
                                                     background = NA)) # Value to use to fill empty cells)
  }
  # plot(region_stack) 
  
  # Aggregate all regions in one layer
  all_regions_raster <- calc(x = region_stack, fun = median, na.rm = T)
  # Add null values for continental borders
  residuals_raster <- continental_mask
  residuals_raster@data@values[!is.na(all_regions_raster[])] <- all_regions_raster[!is.na(all_regions_raster[])]
  
  # plot(all_regions_raster)
  # plot(residuals_raster)
  
  ### Export
  if (include_model) 
  { 
    return(list(residuals_raster, model)) 
  } else {
    return(residuals_raster)
  }
}

source(file = "./functions/compute_residuals.R")

# 20.3.2/ Compute residuals within regions


# Load within betadiversity
load(file = paste0("./outputs/Indices_maps/within_regional_betadiversity.RData"))
TBD_within_regions <- within_regional_betadiversity[[2]][[1]]
PBD_within_regions <- within_regional_betadiversity[[2]][[2]]

### Should use TBD only from species in the phylogeny to have the same species pool than TBD

# Compute
PBD_TBD_regions_residuals <- map_residuals_within_regions (y_index = PBD_within_regions,
                                                           x_index = TBD_within_regions, 
                                                           list_shp_regions, # Provide region Spatial Polygon to define borders
                                                           method = "GAM", # Either LM for linear relationship, LM_quadratic to add quadratic effects, or GLS to account for spatial autocorrelation, or GAM to allow non-linear relationship.
                                                           corSpatial = "corExp", # To choose the type of spatial autocorrelation structure as in nlme::corSpatial() for GLS
                                                           include_model = F, # To include the fitted model in the output
                                                           scatterplot = T,  # To display the scatter plot with regression line and R²
                                                           plot_variogram_correlogram = F) # To display the variogram and correlogram in case of spatial GLS

# Plot
plot(PBD_TBD_regions_residuals, col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))

# Save
save(PBD_TBD_regions_residuals, file = paste0("./outputs/Indices_maps/PBD_TBD_regions_residuals.RData"))

load(file = paste0("./outputs/Indices_maps/PBD_TBD_regions_residuals.RData"))


##### 21/ Map SES-values #####

# SES-tests to account for deviance (Leprieur et al., 2012 ; Rosauer & Jetz, 2014). 
# Need to design a specific one for each case depending on the null model for randomization
  # Can select the type of null model:
    # "keep_phylo" = randomize tips on the tree. Keep the initial phylogeny
    # "Yule" = generate random phylogenies based on a Yule birth_only model parametrized on the phylogeny
    # "BD" = generate random phylogenies based on a Birth-Death model parametrized on the phylogeny
    # "AIC_selection" = chose between a Yule and a Birth-Death model based on AIC
# Can provide an additionnal parameter for subsampling in case you do not have data for all species in the studied clade
  # nsp_tot = total number of extant species in the studied clade (including the one absent from the phylogeny)

# For PD/MPD/TBD vs. SR => randomize names in the phylogeny OR generate phylogeny from null model before recomputing PD/MPD/TPD
# For mimetic richness vs. SR => randomize mimicry membership and recompute mimicry richness

# Two options for the output 
  # Report the SES value normalized according to the mean and sd of the null distribution
  # Report the quantile

### 21.1/ SES-PD according to SR ####

### 21.1.1/ Function to map SES-PD according to SR ####

compute_SES_Faith_PD_vs_SR <- function (proba_stack, phylo, seed = 1, # Set seed to ensure reproducibility
                                        null_model = "keep_phylo",  # "keep_phylo" , "Yule", "BD", "AIC_selection"
                                        nsp_tot = NULL, # To provide the total number of species in the clade for subsampling in case data do not encompasses all species in the studied clade
                                        randomizations = 999, # To select the number of randomization
                                        alpha = 0.05,  # Threshold for significance for a two-sided test
                                        export_null_rasters = F, # Include the raster Stack of PD obtained from randomization in the output
                                        export_null_model = F, # Include the best fitted model for simulating phylogeny
                                        plot_maps = T, # To plot the SES-maps
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
  
  ### Prepare null model for phylogeny generation
  
  # Compute the incomplete sampling fraction (i.e., proportion of species present in the tree vs. actual number of species in the clade)
  if (is.null(nsp_tot))  # If no total number of species is provided, we assume the phylogeny is complete
  { 
    nsp_tot <- length(phylo$tip.label)
  } 
  sampling_fraction <- length(phylo$tip.label) / nsp_tot 
  if (sampling_fraction > 1) { stop("'nsp_tot' must be higher than the number of tips in the phylogeny\nIf your phylogeny includes species outside of the studied clade, please remove them before running the function") }
  
  # Null model based on Yule Birth-only model to generate random phylogeny
  if (null_model %in% c("Yule", "AIC_selection"))
  {
    fit_Yule <- phytools::fit.yule(tree = phylo, rho = sampling_fraction)
    print(fit_Yule)
    
    if (null_model == "Yule")
    {
      best_fit <- fit_Yule
    }
  }
  
  # Null model based on Birth-Death model to generate random phylogeny
  if (null_model %in% c("BD", "AIC_selection"))
  {
    fit_BD <- phytools::fit.bd(tree = phylo, rho = sampling_fraction)
    print(fit_BD)
    
    if (null_model == "BD")
    {
      best_fit <- fit_BD
    }
  }
  
  # Chose the best model based on AIC
  if (null_model == "AIC_selection")
  {
    AIC_df <- AIC(fit_Yule, fit_BD) ; print(AIC_df)
    best_fit <- list(fit_Yule, fit_BD)[[which.min(AIC_df$AIC)]]
  }
  
  ### Loop for each randomization
  PD_stack <- stack(Faith_PD_obs)
  for (i in 1:randomizations)
  {
    cat(paste0("\n", Sys.time(), " - Randomization - ", i," on ",randomizations,"\n"))
    
    ### Generate phylogeny for the null model
    
    # Null model based on tip randomization
    if (null_model == "keep_phylo")
    {
      # Randomize species in the phylogeny
      random_phylo <- pruned_tree
      random_phylo$tip.label <- sample(random_phylo$tip.label)
    }
    
    # Null model based on Yule or Birth-Death model
    if (null_model %in% c("Yule", "BD", "AIC_selection"))
    {
      # Generate the phylogeny
      # ?phytools::pbtree
      
      full_tree <- phytools::pbtree(n = nsp_tot,
                                    b = best_fit$b,
                                    d = best_fit$d,
                                    nsim = 1,
                                    type = "continuous",
                                    extant.only = T)
      
      # Recale so the root has the same age than the initial phylogeny
      initial_phylo_depth <- max(ape::node.depth.edgelength(phylo))
      full_tree_rescaled <- geiger::rescale(x = full_tree, model = "depth", depth = initial_phylo_depth)
      
      # plot(random_phylo)
      
      # Trim the proper subsample of species
      random_phylo <- ape::keep.tip(full_tree_rescaled, sample(x = 1:nsp_tot, size = length(pruned_tree$tip.label)))
      random_phylo$tip.label <- sample(pruned_tree$tip.label)
    }
    
    ### Compute Faith_PD for this randomized phylogeny
    random_PD <- compute_Faith_PD(proba_stack_for_phylo, random_phylo, quiet = T)
    names(random_PD) <- paste0("Faith_PD_random_", i)
    
    ### Add to the final Raster Stack
    PD_stack <- addLayer(PD_stack, random_PD)
  }
  
  ### Plot null distribution for a community if requested
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
  
  if(plot_maps)
  {
    par(mfrow = c(2,2))
    
    plot(SES_raster, main = "SES-PD", col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    plot(SES_signif, main = paste0("SES-PD significant areas\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    
    plot(p_value_raster, zlim = c(0, 1), main = paste0("p-values from randomization tests\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    plot(p_value_signif, zlim = c(0, 1), main = paste0("Significant areas from randomization tests\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    
    par(mfrow = c(1,1))
  }
  
  ### Export
  
  final_output <- list(SES_raster, SES_signif, p_value_raster, p_value_signif)
  # Include the Stack of randomized community if requested
  if (export_null_rasters) { final_output <- append(final_output, list(PD_stack)) }
  # Include best null model for phylogeny simulation if requested
  if (null_model %in% c("Yule", "BD", "AIC_selection"))
  {
    if (export_null_model) { final_output <- append(final_output, list(best_fit)) }
  }

  return(final_output)
  
}

source(file = "./functions/compute_SES_tests.R")

### 21.1.2/ Map SES-PD according to SR ####

SES_Faith_PD_vs_SR <- compute_SES_Faith_PD_vs_SR(sp_proba_stack, phylo, seed = 1,
                                                 null_model = "keep_phylo",
                                                 randomizations = 99, # To select the number of randomization
                                                 alpha = 0.05,  # Thresold for significance for a two-sided test
                                                 export_null_rasters = F, # Include the raster Stack of PD obtained from randomization in the output
                                                 plot_distri = T) # To plot the null distribution for an "average" community as an example

SES_Faith_PD_vs_SR_model <- compute_SES_Faith_PD_vs_SR(sp_proba_stack, phylo, seed = 1,
                                                       null_model = "AIC_selection",  # "keep_phylo" , "Yule", "BD", "AIC_selection"
                                                       nsp_tot = 396, # To provide the total number of species in the clade for subsampling in case data do not encompasses all specie sin the studied clade
                                                       randomizations = 99, # To select the number of randomization
                                                       alpha = 0.05,  # Thresold for significance for a two-sided test
                                                       export_null_rasters = F, # Include the raster Stack of PD obtained from randomization in the output
                                                       plot_distri = T) # To plot the null distribution for an "average" community as an example

# Save
save(SES_Faith_PD_vs_SR, file = paste0("./outputs/Indices_maps/SES_Faith_PD_vs_SR.RData"))
save(SES_Faith_PD_vs_SR_model, file = paste0("./outputs/Indices_maps/SES_Faith_PD_vs_SR_model.RData"))

# Load
load(file = paste0("./outputs/Indices_maps/SES_Faith_PD_vs_SR.RData"))
load(file = paste0("./outputs/Indices_maps/SES_Faith_PD_vs_SR_model.RData"))



### 21.2/ SES-MPD according to SR ####

### 21.2.1/ Function to map SES-MPD according to SR ####

compute_SES_MPD_vs_SR <- function (proba_stack, phylo, seed = 1, # Set seed to ensure reproducibility
                                   null_model = "keep_phylo",  # "keep_phylo" , "Yule", "BD", "AIC_selection"
                                   nsp_tot = NULL, # To provide the total number of species in the clade for subsampling in case data do not encompasses all species in the studied clade
                                   randomizations = 999, # To select the number of randomization
                                   quiet = T, # To display or not internal progress for each MPD computation
                                   alpha = 0.05,  # Threshold for significance for a two-sided test
                                   export_null_rasters = F, # Include the raster Stack of PD obtained from randomization in the output
                                   export_null_model = F, # Include the best fitted model for simulating phylogeny
                                   plot_maps = T, # To plot the SES-maps
                                   plot_distri = T, # To plot the null distribution for an "average" community as an example
                                   example_coords = NULL) # To provide the coordinates of the community to use as example, otherwise a community with median species richness is used.
{
  # Set the random seed for reproducibility
  set.seed(seed = seed)
  
  # Match raster Stack and Phylogeny species lists
  clean_stack_and_phylo <- match_stack_and_phylo(proba_stack, phylo)
  proba_stack_for_phylo <- clean_stack_and_phylo[[1]]
  pruned_tree <- clean_stack_and_phylo[[2]]
  
  # Compute observed MPD
  cat(paste0(Sys.time(), " - Compute observed MPD\n"))
  MPD_obs <- compute_MPD(proba_stack = proba_stack_for_phylo, phylo = pruned_tree, quiet = quiet)
  names(MPD_obs) <- "MPD_obs"
  
  ### Prepare null model for phylogeny generation
  
  # Compute the incomplete sampling fraction (i.e., proportion of species present in the tree vs. actual number of species in the clade)
  if (is.null(nsp_tot))  # If no total number of species is provided, we assume the phylogeny is complete
  { 
    nsp_tot <- length(phylo$tip.label)
  } 
  sampling_fraction <- length(phylo$tip.label) / nsp_tot 
  if (sampling_fraction > 1) { stop("'nsp_tot' must be higher than the number of tips in the phylogeny\nIf your phylogeny includes species outside of the studied clade, please remove them before running the function") }
  
  # Null model based on Yule Birth-only model to generate random phylogeny
  if (null_model %in% c("Yule", "AIC_selection"))
  {
    fit_Yule <- phytools::fit.yule(tree = phylo, rho = sampling_fraction)
    print(fit_Yule)
    
    if (null_model == "Yule")
    {
      best_fit <- fit_Yule
    }
  }
  
  # Null model based on Birth-Death model to generate random phylogeny
  if (null_model %in% c("BD", "AIC_selection"))
  {
    fit_BD <- phytools::fit.bd(tree = phylo, rho = sampling_fraction)
    print(fit_BD)
    
    if (null_model == "BD")
    {
      best_fit <- fit_BD
    }
  }
  
  # Chose the best model based on AIC
  if (null_model == "AIC_selection")
  {
    AIC_df <- AIC(fit_Yule, fit_BD) ; print(AIC_df)
    best_fit <- list(fit_Yule, fit_BD)[[which.min(AIC_df$AIC)]]
  }
  
  
  ### Loop for each randomization
  MPD_stack <- stack(MPD_obs)
  for (i in 1:randomizations)
  {
    cat(paste0(Sys.time(), " - Randomization - ", i," on ",randomizations,"\n"))
    
    ## Generate phylogeny for the null model
    
    # Null model based on tip randomization
    if (null_model == "keep_phylo")
    {
      # Randomize species in the phylogeny
      random_phylo <- pruned_tree
      random_phylo$tip.label <- sample(random_phylo$tip.label)
    }
    
    # Null model based on Yule or Birth-Death model
    if (null_model %in% c("Yule", "BD", "AIC_selection"))
    {
      # Generate the phylogeny
      # ?phytools::pbtree
      
      full_tree <- phytools::pbtree(n = nsp_tot,
                                    b = best_fit$b,
                                    d = best_fit$d,
                                    nsim = 1,
                                    type = "continuous",
                                    extant.only = T)
      
      # Recale so the root has the same age than the initial phylogeny
      initial_phylo_depth <- max(ape::node.depth.edgelength(phylo))
      full_tree_rescaled <- geiger::rescale(x = full_tree, model = "depth", depth = initial_phylo_depth)
      
      # plot(random_phylo)
      
      # Trim the proper subsample of species
      random_phylo <- ape::keep.tip(full_tree_rescaled, sample(x = 1:nsp_tot, size = length(pruned_tree$tip.label)))
      random_phylo$tip.label <- sample(pruned_tree$tip.label)
    }
    
    ## Compute MPD for this randomized phylogeny
    random_MPD <- compute_MPD(proba_stack_for_phylo, random_phylo, quiet = quiet)
    names(random_MPD) <- paste0("MPD_random_", i)
    
    ## Add to the final Raster Stack
    MPD_stack <- addLayer(MPD_stack, random_MPD)
  }
  
  ### Plot null distribution for a community if requested
  if(plot_distri)
  {
    # If no example community is provided
    if (is.null(example_coords))
    {
      # Sample an "average" community
      com_ranks <- rank(MPD_obs[], ties.method = "random", na.last = "keep")
      sample_com_index <- which(round(median(com_ranks, na.rm = T)) == com_ranks)
      
      # Get its coordinates
      sample_com_coords <- round(raster::xyFromCell(object = MPD_obs, cell = sample_com_index),3)
      
    } else { # If the example community coordinates are provided
      
      sample_com_coords <- example_coords
      
      # Use the provided coordinates to retrieve its index
      sample_com_index <- raster::cellFromXY(object = MPD_obs, xy = sample_com_coords)
    }
    
    # Extract statistics
    com_null_distri <- getValues(MPD_stack)[sample_com_index, ]
    com_MPD_obs <- com_null_distri[1]
    mean_null_value <- mean(com_null_distri, na.rm = T)
    Q2.5_value <- quantile(com_null_distri, 0.025, na.rm = T)
    Q97.5_value <- quantile(com_null_distri, 0.975, na.rm = T)
    p_value <- round(ecdf(x = com_null_distri)(com_MPD_obs), 3)
    SES_value <- round((com_MPD_obs - mean_null_value) / sd(com_null_distri), 3)
    
    # Plot null distribution for an "average" community
    y_counts <- hist(com_null_distri, breaks = 30, plot = F)$counts
    hist(com_null_distri,
         breaks = 30, freq = TRUE, col = "gray", 
         main = paste0("Null distribution of the MPD in a community\nLatitude: ", sample_com_coords[2], "; Longitude: ", sample_com_coords[1]), 
         xlab = "Mean pairwise Phylogenetic Distance  [My]",
         cex.axis = 1.5, cex.lab = 1.3, cex.main = 1.6, lwd = 2)
    arrows(x0 = com_MPD_obs, y0 = max(y_counts)*0.5, x1 = com_MPD_obs, y1 = max(y_counts)*0.1, length = 0.1, lwd = 3)
    abline(v = mean(com_null_distri), lty = 2, lwd = 2)
    abline(v = Q2.5_value, lty = 2, lwd = 2, col = "red")
    abline(v = Q97.5_value, lty = 2, lwd = 2, col = "red")
    
    # Insert mean and quantiles legend
    legend(legend = c(paste0("Mean = ", format(round(mean_null_value, 1), nsmall = 1)), 
                      paste0("CI 2.5% = ", format(round(Q2.5_value, 1), nsmall = 1)),
                      paste0("CI 97.5% = ", format(round(Q97.5_value, 1), nsmall = 1))), 
           x = "topright", inset = c(0.02, 0.02), y.intersp = 1.2, lty = 2 , lwd = 2,
           col = c("black", "red"), cex = 1.2, bty = "o", bg = "white", box.col = NA)
    
    # Insert MPD obs, SES and p-value legend
    legend(legend = c(paste0("MPD obs = ", round(com_MPD_obs, 2)),
                      paste0("SES-MPD = ", format(SES_value, nsmall = 3)),
                      paste0("p = ", format(p_value, nsmall = 3))),
           x = "topleft", inset = c(0.02, 0.02),
           cex = 1.2, xjust = 1, bty = "o", bg = "white", box.col = NA)
  }
  
  ### Compute SES-values
  
  # Compute the mean of the null distribution
  mean_null_raster <- raster::calc(x = MPD_stack, fun = mean, na.rm = T)
  # Compute the sd of the null distribution
  sd_null_raster <- raster::calc(x = MPD_stack, fun = sd, na.rm = T)
  # Compute the SES-values
  SES_raster <- round(((MPD_obs - mean_null_raster) / sd_null_raster), 3)
  
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
  MPD_stack_values <- raster::getValues(MPD_stack)
  all_p_values <- NA
  for (i in 1:nrow(MPD_stack_values))
  {
    # Extract community data
    com_data <- MPD_stack_values[i, ]
    
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
  
  if (plot_maps)
  {
    par(mfrow = c(2,2))
    
    plot(SES_raster, main = "SES-MPD", col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    plot(SES_signif, main = paste0("SES-MPD significant areas\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    
    plot(p_value_raster, zlim = c(0, 1), main = paste0("p-values from randomization tests\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    plot(p_value_signif, zlim = c(0, 1), main = paste0("Significant areas from randomization tests\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    
    par(mfrow = c(1,1))
  }

  ### Export
  
  final_output <- list(SES_raster, SES_signif, p_value_raster, p_value_signif)
  # Include the Stack of randomized community if requested
  if (export_null_rasters) { final_output <- append(final_output, list(PD_stack)) }
  # Include best null model for phylogeny simulation if requested
  if (null_model %in% c("Yule", "BD", "AIC_selection"))
  {
    if (export_null_model) { final_output <- append(final_output, list(best_fit)) }
  }

  return(final_output)
  
}

source(file = "./functions/compute_SES_tests.R")

### 21.1.2/ Map SES-MPD according to SR ####

SES_MPD_vs_SR <- compute_SES_MPD_vs_SR(proba_stack = sp_proba_stack,
                                       phylo = Ithomiini_phylo, seed = 1,
                                       randomizations = 9, # To select the number of randomization
                                       alpha = 0.05,  # Threshold for significance for a two-sided test
                                       export_null_rasters = F, # Include the raster Stack of PD obtained from randomization in the output
                                       plot_distri = T) # To plot the null distribution for an "average" community as an example


# Save
save(SES_MPD_vs_SR, file = paste0("./outputs/Indices_maps/SES_MPD_vs_SR.RData"))



### 21.3/ SES-PBD according to TBD ####

### 21.3.1/ Function to compute SES-PBD according to TBD ####

compute_SES_PBD_vs_TBD <- function (proba_stack, phylo, 
                                    index_family = "sorensen", # Chose the family of Beta-diversity indices
                                    beta_part = "turnover",    # Chose the partition
                                    aggreg_type = "all_pairwise",     # Chose the aggregation method
                                    window_width = "default", window_height = "default", # In pixels, extended from the focal point
                                    seed = 1, # Set seed to ensure reproducibility
                                    randomizations = 999, # To select the number of randomization
                                    alpha = 0.05,  # Threshold for significance for a two-sided test
                                    export_null_rasters = F, # Include the raster Stack of PBD obtained from randomization in the output
                                    plot_maps = T, # To plot the SES-maps
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
  
  if (plot_maps)
  {
    par(mfrow = c(2,2))
    
    plot(SES_raster, main = "SES-PBD", col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    plot(SES_signif, main = paste0("SES-PBD significant areas\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    
    plot(p_value_raster, zlim = c(0, 1), main = paste0("p-values from randomization tests\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    plot(p_value_signif, zlim = c(0, 1), main = paste0("Significant areas from randomization tests\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    
    par(mfrow = c(1,1))
  }

  
  ### Export
  
  final_output <- list(SES_raster, SES_signif, p_value_raster, p_value_signif)
  # Include the Stack of randomized community if requested
  if (export_null_rasters) { final_output <- append(final_output, list(PBD_stack)) }
  
  return(final_output)
  
}

source(file = "./functions/compute_SES_tests.R")

### 21.3.2/ Map SES-PBD according to TBD ####

SES_PBD_vs_TBD <- compute_SES_PBD_vs_TBD(binary_stack, phylo, 
                                         index_family = "sorensen", # Chose the family of Beta-diversity indices
                                         beta_part = "turnover",    # Chose the partition
                                         aggreg_type = "all_pairwise",     # Chose the aggregation method
                                         window_width = "default", window_height = "default", # In pixels, extended from the focal point
                                         seed = 1, randomizations = 99, # To select the number of randomization
                                         alpha = 0.05,  # Thresold for significance for a two-sided test
                                         export_null_rasters = F, # Include the raster Stack of PBD obtained from randomization in the output
                                         plot_distri = T, # To plot the null distribution for an "average" community as an example)
                                         example_coords = NULL) # Do not provide coordinates for the example community to display the null distribution

# Save
save(SES_PBD_vs_TBD, file = paste0("./outputs/Indices_maps/SES_PBD_vs_TBD.RData"))
# save(SES_PBD_vs_TBD, file = paste0("./outputs/Indices_maps/SES_PBD_vs_TBD_3x3.RData"))
# save(SES_PBD_vs_TBD, file = paste0("./outputs/Indices_maps/SES_PBD_vs_TBD_5x5.RData"))


### 21.4/ SES-PBD according to TBD within regions ####

### 21.4.1/ Function to compute SES-PBD according to TBD within regions ####


compute_SES_PBD_vs_TBD_within_regions <- function (proba_stack, list_shp_regions, region_names,
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
                                                   plot_maps = T, # To plot the SES-maps
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
  
  PBD_obs_results <- compute_within_regional_betadiversity(binary_stack, list_shp_regions = list_shp_regions, region_names = region_names,
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
    random_region_df <- compute_within_regional_betadiversity(binary_stack, list_shp_regions = list_shp_regions, region_names = region_names,
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
  list_shp_regions <- list_shp_regions[!check_NaN]
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
  
  # Add region names to each SpPolygonsDF
  for (i in 1: length(list_shp_regions))
  {
    n_poly <- nrow(list_shp_regions[[i]]@data)
    list_shp_regions[[i]]@data <- cbind(list_shp_regions[[i]]@data, data.frame("region" = rep(region_names[i], times = n_poly)))
  }
  
  # Aggregate SpatialPolygons for regions
  if (sum(!check_NaN) == 1)
  { # Case with only one region remaining with data
    all_shp_polygons <- list_shp_regions
    
  } else {
    # Case with multiple regions remaining with data
    all_shp_polygons <- do.call(raster::bind, list_shp_regions)
  }
  
  # Add the region stats
  region_stat_df_with_region_col <- data.frame(region = row.names(region_stat_df), region_stat_df)
  all_shp_polygons <- sp::merge(x = all_shp_polygons, y = region_stat_df_with_region_col, by = "region")
  
  ### Case with only one polygon per region
  # all_shp_polygons@data <- data.frame(all_shp_polygons@data, regions = row.names(region_stat_df), region_stat_df)
  
  # Rasterize using the SES-PBD value
  SES_raster <- raster::rasterize(x = all_shp_polygons, y = binary_stack, 
                                  field = "SES_PBD", # Variable to use to fix the pixel values within regions
                                  background = NA)  # Value to use to fill empty cells
  
  # Rasterize using the p-values from randomization tests
  p_value_raster <- raster::rasterize(x = all_shp_polygons, y = binary_stack, 
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
  
  if (plot_maps)
  {
    par(mfrow = c(2,2))
    
    plot(SES_raster, main = "SES-PBD", col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    plot(SES_signif, main = paste0("SES-PBD significant areas\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    
    plot(p_value_raster, zlim = c(0, 1), main = paste0("p-values from randomization tests\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    plot(p_value_signif, zlim = c(0, 1), main = paste0("Significant areas from randomization tests\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    
    par(mfrow = c(1,1))
  }

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

source(file = "./functions/compute_SES_tests.R")

### 21.4.2/ Map SES-PBD according to TBD within regions ####

# Load spPolygon shape files for bioregions as example
load(file = "./input_data/Map_stuff/Bioregions/All_bioregions_in_figure.RData")

# list_shp_regions <- list(Caatinga_shp, Guyana_Shield_shp)
# region_names <- c("Caatinga", "Guyana_Shield")

# list_shp_regions <- list(Central_Andes_shp, Western_Amazon_shp)
# region_names <- c("Central Andes", "Western Amazon")

# plot(binary_stack)

list_shp_regions <- list(Caatinga_shp, Caribbean_Islands_shp, full_CA_shp, Central_Andes_shp, Cerrado_shp, Chacos_shp, Coastal_desert_shp, Guyana_Shield_shp, Llanos_shp, Lower_Amazon_shp, Mata_Atlantica_shp5, Northern_Andes_shp, Pampas_shp, Pantanal_shp, Western_Amazon_shp, Western_Lowlands_shp)
region_names <- c("Caatinga", "Caribbean_Islands", "Central_America", "Central_Andes", "Cerrado", "Chacos", "Coastal_desert", "Guyana_Shield", "Llanos", "Lower_Amazon", "Mata_Atlantica", "Northern_Andes", "Pampas", "Pantanal", "Western_Amazon", "Western_Lowlands")


SES_PBD_vs_TBD_within_regions <- compute_SES_PBD_vs_TBD_within_regions(proba_stack = binary_stack, 
                                                                       list_shp_regions, region_names,
                                                                       subsample_size = 1000, # Provide number of pixels to subsample globally and regularly in order to save computing time
                                                                       # subsample_size = NA, # Provide number of pixels to subsample globally and regularly in order to save computing time
                                                                       phylo, # Provide the phylogeny
                                                                       index_family = "sorensen", # Chose the family of Beta-diversity indices
                                                                       beta_part = "turnover",    # Chose the partition
                                                                       aggreg_type = "mean_pairwise",     # Chose the aggregation method
                                                                       seed = 1, randomizations = 999, # To select the number of randomization
                                                                       alpha = 0.05,  # Threshold for significance for a two-sided test
                                                                       export_null_df = F, # To include the df with indices values from all randomization in the output
                                                                       export_region_stats_df = T, # To include the df with statistics from the tests for each region
                                                                       plot_distri = T, # To plot the null distribution for an "average" community as an example
                                                                       example_region_name = "Central Andes") # To provide the name of the region to use as example. Otherwise it is randomly chosen.


# Save
save(SES_PBD_vs_TBD_within_regions, file = paste0("./outputs/Indices_maps/SES_PBD_vs_TBD_within_regions.RData"))
# save(SES_PBD_vs_TBD_within_regions, file = paste0("./outputs/Indices_maps/SES_PBD_vs_TBD_within_regions_1000.RData"))



### 21.5/ SES-MR according to SR/OMU richness ####

##### Can be adapted for mimicry richness vs. SR => randomize mimicry membership and recompute mimicry richness #####

### But need to find a way to keep SR fixed while mimicry membership is randomized since it will be randomized across OMU and not species...
# Or simply use "OMU-richness" and randomize across OMU before to recompute mimicry richness without caring for species richness


### 21.5.1/ Function to compute SES-MR according to taxonomic richness ####

### Change for taxonomic richness but specify can be used at species level, and OMU level if polymorphism ####

compute_SES_MR_vs_taxo_richness <- function (proba_stack, # Species or OMU continuous (or binary) Raster Stack of SDM outputs
                                             mimicry_ring_membership, # Character string vector of mimicry ring membership per species/OMU in the same order than the Stack
                                             seed = 1, # Set seed to ensure reproducibility
                                             randomizations = 999, # To select the number of randomization
                                             alpha = 0.05,  # Threshold for significance for a two-sided test
                                             export_null_rasters = F, # Include the raster Stack of MR obtained from randomization in the output
                                             plot_maps = T, # To plot the SES-maps
                                             plot_distri = T, # To plot the null distribution for an "average" community as an example
                                             example_coords = NULL, # To provide the coordinates of the community to use as example, otherwise a community with median species richness is used.
                                             quiet = F) # Should each generation of mimicry ring stack be quiet or not
{
  
  # Set the random seed for reproducibility
  set.seed(seed = seed)
  
  cat(paste0(Sys.time(), " - Compute observed Mimicry richness\n"))
  
  # Generate observed mimicry proba stack
  ring_proba_stack_obs <- generate_mimicry_ring_stack(proba_stack = proba_stack, # Species or OMU continuous (or binary) Raster Stack of SDM outputs
                                                      mimicry_ring_membership = mimicry_ring_membership, # Character string vector of mimicry ring membership per species/OMU in the same order than the Stack
                                                      stack_type = "proba", quiet = quiet)
  
  # Compute observed mimicry richness (MR)
  MR_obs <- compute_richness(ring_proba_stack_obs)
  names(MR_obs) <- "MR_obs"
  
  
  ### Loop for each randomization
  MR_stack <- stack(MR_obs)
  for (i in 1:randomizations)
  {
    cat(paste0(Sys.time(), " - Randomization - ", i," on ",randomizations,"\n"))
    
    # Randomize mimicry ring membership
    random_mimicry_ring_membership <- as.character(mimicry_ring_membership)
    random_mimicry_ring_membership <- sample(random_mimicry_ring_membership)
    
    # Compute mimicry ring proba stack for this random mimicry ring membership scheme
    random_MR_stack <- generate_mimicry_ring_stack(proba_stack = proba_stack, # Species or OMU continuous (or binary) Raster Stack of SDM outputs
                                                   mimicry_ring_membership = random_mimicry_ring_membership, # Character string vector of mimicry ring membership per species/OMU in the same order than the Stack
                                                   stack_type = "proba", quiet = quiet)
    
    # Compute mimicry richness for this random mimicry ring membership scheme
    random_MR <- compute_richness(random_MR_stack)
    names(random_MR) <- paste0("MR_random_", i)
    
    # Add to the final Raster Stack
    MR_stack <- addLayer(MR_stack, random_MR)
  }
  
  if(plot_distri)
  {
    # If no example community is provided
    if (is.null(example_coords))
    {
      # Sample an "average" community
      com_ranks <- rank(MR_obs[], ties.method = "random", na.last = "keep")
      sample_com_index <- which(round(median(com_ranks, na.rm = T)) == com_ranks)
      
      # Get its coordinates
      sample_com_coords <- round(raster::xyFromCell(object = MR_obs, cell = sample_com_index),3)
      
    } else { # If the example community coordinates are provided
      
      sample_com_coords <- example_coords
      
      # Use the provided coordinates to retrieve its index
      sample_com_index <- raster::cellFromXY(object = MR_obs, xy = sample_com_coords)
    }
    
    # Extract statistics
    com_null_distri <- getValues(MR_stack)[sample_com_index, ]
    com_MR_obs <- com_null_distri[1]
    mean_null_value <- mean(com_null_distri, na.rm = T)
    Q2.5_value <- quantile(com_null_distri, 0.025, na.rm = T)
    Q97.5_value <- quantile(com_null_distri, 0.975, na.rm = T)
    p_value <- round(ecdf(x = com_null_distri)(com_MR_obs), 3)
    SES_value <- round((com_MR_obs - mean_null_value) / sd(com_null_distri), 3)
    
    # Plot null distribution for an "average" community
    y_counts <- hist(com_null_distri, breaks = 30, plot = F)$counts
    hist(com_null_distri,
         breaks = 30, freq = TRUE, col = "gray", 
         main = paste0("Null distribution of the mimicry richness in a community\nLatitude: ", sample_com_coords[2], "; Longitude: ", sample_com_coords[1]), 
         xlab = "Mimicry richness  [rings]",
         cex.axis = 1.5, cex.lab = 1.3, cex.main = 1.6, lwd = 2)
    arrows(x0 = com_MR_obs, y0 = max(y_counts)*0.5, x1 = com_MR_obs, y1 = max(y_counts)*0.1, length = 0.1, lwd = 3)
    abline(v = mean(com_null_distri), lty = 2, lwd = 2)
    abline(v = Q2.5_value, lty = 2, lwd = 2, col = "red")
    abline(v = Q97.5_value, lty = 2, lwd = 2, col = "red")
    
    # Insert mean and quantiles legend
    legend(legend = c(paste0("Mean = ", format(round(mean_null_value, 1), nsmall = 1)), 
                      paste0("CI 2.5% = ", format(round(Q2.5_value, 1), nsmall = 1)),
                      paste0("CI 97.5% = ", format(round(Q97.5_value, 1), nsmall = 1))), 
           x = "topright", inset = c(0.02, 0.02), y.intersp = 1.2, lty = 2 , lwd = 2,
           col = c("black", "red"), cex = 1.2, bty = "o", bg = "white", box.col = NA)
    
    # Insert MR obs, SES and p-value legend
    legend(legend = c(paste0("MR obs = ", round(com_MR_obs, 2)),
                      paste0("SES-MR = ", format(SES_value, nsmall = 3)),
                      paste0("p = ", format(p_value, nsmall = 3))),
           x = "topleft", inset = c(0.02, 0.02),
           cex = 1.2, xjust = 1, bty = "o", bg = "white", box.col = NA)
  }
  
  ### Compute SES-values
  
  # Compute the mean of the null distribution
  mean_null_raster <- raster::calc(x = MR_stack, fun = mean, na.rm = T)
  # Compute the sd of the null distribution
  sd_null_raster <- raster::calc(x = MR_stack, fun = sd, na.rm = T)
  # Compute the SES-values
  SES_raster <- round(((MR_obs - mean_null_raster) / sd_null_raster), 3)
  
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
  
  ### Compute p-values from randomization tests
  
  # Loop per community
  MR_stack_values <- raster::getValues(MR_stack)
  all_p_values <- NA
  for (i in 1:nrow(MR_stack_values))
  {
    # Extract community data
    com_data <- MR_stack_values[i, ]
    
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
  
  if (plot_maps)
  {
    par(mfrow = c(2,2))
    
    plot(SES_raster, main = "SES-MR", col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    plot(SES_signif, main = paste0("SES-MR significant areas\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    
    plot(p_value_raster, zlim = c(0, 1), main = paste0("p-values from randomization tests\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    plot(p_value_signif, zlim = c(0, 1), main = paste0("Significant areas from randomization tests\nalpha = ", alpha), col = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(200))
    
    par(mfrow = c(1,1))
  }

  ### Export
  
  final_output <- list(SES_raster, SES_signif, p_value_raster, p_value_signif)
  # Include the Stack of randomized community if requested
  if (export_null_rasters) { final_output <- append(final_output, list(MR_stack)) }
  
  return(final_output)
  
}

source(file = "./functions/compute_SES_tests.R")

### 21.5.2/ Map SES-MR according to SR ####

# Load OMU stack and mimicry membership info
OMU_proba_stack <- readRDS("./outputs/Indices_stacks/All_OMU_stack_Jaccard.80.rds")
OMU_df <- readRDS("./input_data/list.models.rds")

# Check names and order of OMUs are identical in stack and dataframe
identical(x = names(OMU_proba_stack), OMU_df$Tag.model)


SES_MR_vs_taxo_richness <- compute_SES_MR_vs_taxo_richness(proba_stack = OMU_proba_stack, # Species or OMU continuous (or binary) Raster Stack of SDM outputs
                                                           mimicry_ring_membership = OMU_df$Mimicry.model, # Character string vector of mimicry ring membership per species/OMU in the same order than the Stack
                                                           seed = 1, # Set seed to ensure reproducibility
                                                           randomizations = 999, # To select the number of randomization
                                                           alpha = 0.05,  # Threshold for significance for a two-sided test
                                                           export_null_rasters = F, # Include the raster Stack of MR obtained from randomization in the output
                                                           plot_distri = T, # To plot the null distribution for an "average" community as an example
                                                           example_coords = NULL, # To provide the coordinates of the community to use as example, otherwise a community with median species richness is used.
                                                           quiet = T)

# Save
save(SES_MR_vs_taxo_richness, file = paste0("./outputs/Indices_maps/SES_MR_vs_taxo_richness.RData"))





# ### Generate smaller test sets  ###
# 
# proba_stack_test <- stack(raster::crop(x = sp_proba_stack[[20:29]], y = extent(-81.5, -77.2, -2.8, 0)))
# proba_stack_test <- stack(raster::crop(x = sp_proba_stack[[20:29]], y = extent(-106, -80, 0, 22)))
# proba_stack_test <- stack(raster::crop(x = sp_proba_stack[[20:29]], y = extent(-80, -60, -20, 0)))
# plot(proba_stack_test)
# phylo <- phylogeny_object
#
###




###### Improvements for every functions ######

# Change function names to better distinguish map (output = raster) and compute (output = matrix or table) functions

# Check with calc

# If not work use the getValues and loop method

# Check purrr or apply instead of loop ?

# Computation by blocks instead of loading everything, such as in raster package...

# Parallelize heavy loops