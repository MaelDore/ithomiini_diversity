
###### Script 26: Elevation profile of indices ######

# Only for Jaccard.80

### Use only 4 indices

# 1/ Species richness
# 2/ Mean Continuous Species Rarity = Range-size Weighted Species Richness
# 3/ Mean pairwise Phylogenetic Distance
# 4/ Mimicry richness

###

###

# Inputs: 
   # Indices maps (4 main indices)
   # Elevation raster

# Outputs:
   # Scatter plots of communities diversity indices ~ elevation

###


# Effacer l'environnement
rm(list = ls())

library(raster)

# Load indices map
sp.richness <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
sp.rarity <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_Jaccard.80.rds")
MPD <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80.rds")
ring.richness <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))

# Create a mask of communities with presence of Ithomiini
Community_mask <- sp.richness
Community_mask[!(Community_mask > 0)] <- NA

plot(Community_mask)

# Load elevation raster
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_15.rds"))
DEM <- envData[["Elevation"]]

# Mask all rasters to constrain to Ithomiini communities
DEM_masked <- mask(DEM, Community_mask)
sp.richness_masked <- mask(sp.richness, Community_mask)
sp.rarity_masked <- mask(sp.rarity, Community_mask)
MPD_masked <- mask(MPD, Community_mask)
ring.richness_masked <- mask(ring.richness, Community_mask)

###

# Plot scatter plots

pairs(stack(DEM_masked, sp.richness_masked))
pairs(stack(DEM_masked, sp.rarity_masked))
pairs(stack(DEM_masked, MPD_masked))
pairs(stack(DEM_masked, ring.richness_masked))


?pairs
