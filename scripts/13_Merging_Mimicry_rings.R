
##### Script 13: Merging Mimicry rings #####

# Create map at mimicry ring level by computing "probability" of presence of any of the OMU of each mimicry ring

# Generate a continuous map for each 4 options between Jaccard/TSS and Buffer.80/95


##### 

# Inputs:
   # Stack of continuous EM and binary maps for each OMU, per mimicry ring from script 11

# Outputs:
   # Mimicry ring probability continuous map for each 4 option between Jaccard/TSS and Buffer.80/95 => aggregated with aggreg_proba function
   # Multiple pages PDF with all mimicry ring probability maps
   # Mimicry ring richness map for each 4 option between Jaccard/TSS and Buffer.80/95 => aggregated with sum
   # Multiple pages PDF with all mimicry ring richness maps
   # Generate a single pdf with all Mimicry ring richness for Jaccard.80

#####


### Prepare stuff

# Effacer l'environnement
rm(list = ls())

library(raster)

# Choose the resolution
# res <-  "5"
res <- "15"

# Load environmental stack to use as reference for extent and CRS and generate a mask for continental borders
library(raster)
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))
continent_mask <- envData[[1]]
continent_mask <- calc(continent_mask, fun = function(x, ...) {x*0}, na.rm = F)
# save(continent_mask, file = paste0("./input_data/Env_data/continent_mask_", res, ".RData"), version = "2")
# saveRDS(continent_mask, file = paste0("./input_data/Env_data/continent_mask_", res, ".rds"), version = "2")


# Load Summary table for OMU/unit and for Species
load(file = paste0("./input_data/list.sp.RData"))
load(file = paste0("./input_data/list.models.RData"))

# Mimicry list
mimicry.list <- as.character(unique(list.models$Mimicry.model)) # 44 Mimicry rings

# Make a list of modeled and non-modeled units (occ.unit = with only occurrences points)
modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]
rasterized_OMU <- list.models[is.na(list.models$Model_ID), ]

# Function to compute probability of presence with multiple OMUs
aggreg_prob = function(x, na.rm) { 
  y <- 1-prod(1-x) # Probability of presence of species = probability of presence of at least one OMU = opposite of probability of absence of all OMU
  return(y) # Output
}


### Loop by mimicry ring
for (i in 1:length(mimicry.list)) 
{
  # i <- 10
  
  # Load a random layer to initiate the stacks (to remove later)
  ring.stack_Jaccard.80 <- ring.stack_Jaccard.95 <- ring.stack_TSS.80 <- ring.stack_TSS.95 <- stack(continent_mask) 
  
  # Load Mimicry ring name
  ring <- as.character(mimicry.list[i])
  
  ### 1/ Make stacks of OMU maps ####
  
  OMU_list <- modeled_OMU[modeled_OMU$Mimicry.model == ring, ]
  if (nrow(OMU_list) > 0) {      # If at least one OMU was modeled
    
    for (j in 1:nrow(OMU_list))  # For each OMU
    { 
      unit <-  as.character(OMU_list$Tag.model[j]) # Load the unit name                  
      
      # Jaccard.80: Load the continuous maps and stack them all
      Ensemble_Jaccard_median_cropped_80 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_Jaccard_median_cropped_80.rds")) 
      ring.stack_Jaccard.80 <- addLayer(ring.stack_Jaccard.80, Ensemble_Jaccard_median_cropped_80)   
      
      # Jaccard.95: Load the continuous maps and stack them all
      Ensemble_Jaccard_median_cropped_95 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_Jaccard_median_cropped_95.rds")) 
      ring.stack_Jaccard.95 <- addLayer(ring.stack_Jaccard.95, Ensemble_Jaccard_median_cropped_95)
      
      # TSS.80: Load the continuous maps and stack them all
      Ensemble_TSS_median_cropped_80 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_TSS_median_cropped_80.rds")) 
      ring.stack_TSS.80 <- addLayer(ring.stack_TSS.80, Ensemble_TSS_median_cropped_80)   
      
      # TSS.95: Load the continuous maps and stack them all
      Ensemble_TSS_median_cropped_95 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_TSS_median_cropped_95.rds")) 
      ring.stack_TSS.95 <- addLayer(ring.stack_TSS.95, Ensemble_TSS_median_cropped_95)
    }
  }
  
  Binaries <- rasterized_OMU[rasterized_OMU$Mimicry.model == ring, ]
  if (nrow(Binaries) > 0) {  # If at least one OMU was rasterized
    
    for (j in 1:nrow(Binaries)) {  # For each  OMU
      
      unit <-  as.character(Binaries$Tag.model[j]) # Load the unit name 
      
      rasterized_map <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/rasterized_map_",unit,".rds"))
      
      # Jaccard.80: Load the binary maps and stack them all
      ring.stack_Jaccard.80 <- addLayer(ring.stack_Jaccard.80, rasterized_map)   
      
      # Jaccard.95: Load the binary maps and stack them all
      ring.stack_Jaccard.95 <- addLayer(ring.stack_Jaccard.95, rasterized_map)
      
      # TSS.80: Load the binary maps and stack them all
      ring.stack_TSS.80 <- addLayer(ring.stack_TSS.80, rasterized_map)   
      
      # TSS.95: Load the binary maps and stack them all
      ring.stack_TSS.95 <- addLayer(ring.stack_TSS.95, rasterized_map)
      
    }
  }
  
  # Drop the useless first layer used to initiate the stack
  ring.stack_Jaccard.80 <- dropLayer(ring.stack_Jaccard.80, i = 1)
  ring.stack_Jaccard.95 <- dropLayer(ring.stack_Jaccard.95, i = 1) 
  ring.stack_TSS.80 <- dropLayer(ring.stack_TSS.80, i = 1)
  ring.stack_TSS.95 <- dropLayer(ring.stack_TSS.95, i = 1)
  
  ### 2/ Compute probability map per ring ####
  
  ring.cont_Jaccard.80 <- calc(ring.stack_Jaccard.80, fun = aggreg_prob)
  ring.cont_Jaccard.95 <- calc(ring.stack_Jaccard.95, fun = aggreg_prob)
  ring.cont_TSS.80 <- calc(ring.stack_TSS.80, fun = aggreg_prob)
  ring.cont_TSS.95 <- calc(ring.stack_TSS.95, fun = aggreg_prob)
  
  plot(ring.cont_Jaccard.80, main = paste0(ring,"\nJaccard.80"))
  
  # Make sure to save the data, not just the link to the temp file!
  ring.cont_Jaccard.80 <- readAll(ring.cont_Jaccard.80)
  ring.cont_Jaccard.95 <- readAll(ring.cont_Jaccard.95)
  ring.cont_TSS.80 <- readAll(ring.cont_TSS.80)
  ring.cont_TSS.95 <- readAll(ring.cont_TSS.95)
  
  save(ring.cont_Jaccard.80, file = paste0("./outputs/Mimicry_rings_proba/Ring_proba_Jaccard.80_",ring,".RData"), version = "2")
  saveRDS(ring.cont_Jaccard.80, file = paste0("./outputs/Mimicry_rings_proba/Ring_proba_Jaccard.80_",ring,".rds"), version = "2")
  save(ring.cont_Jaccard.95, file = paste0("./outputs/Mimicry_rings_proba/Ring_proba_Jaccard.95_",ring,".RData"), version = "2")
  saveRDS(ring.cont_Jaccard.95, file = paste0("./outputs/Mimicry_rings_proba/Ring_proba_Jaccard.95_",ring,".rds"), version = "2")
  save(ring.cont_TSS.80, file = paste0("./outputs/Mimicry_rings_proba/Ring_proba_TSS.80_",ring,".RData"), version = "2")
  saveRDS(ring.cont_TSS.80, file = paste0("./outputs/Mimicry_rings_proba/Ring_proba_TSS.80_",ring,".rds"), version = "2")
  save(ring.cont_TSS.95, file = paste0("./outputs/Mimicry_rings_proba/Ring_proba_TSS.95_",ring,".RData"), version = "2")
  saveRDS(ring.cont_TSS.95, file = paste0("./outputs/Mimicry_rings_proba/Ring_proba_TSS.95_",ring,".rds"), version = "2")
  
  # Multiple pages pdf with all Ring proba maps
  pdf(file = paste0("./maps/Mimicry_rings_proba/Ring_proba_maps_",ring,".pdf"), height = 6, width = 7)
  plot(ring.cont_Jaccard.80, main = paste0(ring,"\nJaccard.80"))
  plot(ring.cont_Jaccard.95, main = paste0(ring,"\nJaccard.95"))
  plot(ring.cont_TSS.80, main = paste0(ring,"\nTSS.80"))
  plot(ring.cont_TSS.95, main = paste0(ring,"\nTSS.95"))
  dev.off()
  
  ### 3/ Compute richness map per ring (nb of OMU) ####
  
  ring.rich_Jaccard.80 <- calc(ring.stack_Jaccard.80, fun = sum)
  ring.rich_Jaccard.95 <- calc(ring.stack_Jaccard.95, fun = sum)
  ring.rich_TSS.80 <- calc(ring.stack_TSS.80, fun = sum)
  ring.rich_TSS.95 <- calc(ring.stack_TSS.95, fun = sum)
  
  plot(ring.rich_Jaccard.95, main = paste0(ring,"\nJaccard.95"))
  
  # Make sure to save the data, not just the link to the temp file!
  ring.rich_Jaccard.80 <- readAll(ring.rich_Jaccard.80)
  ring.rich_Jaccard.95 <- readAll(ring.rich_Jaccard.95)
  ring.rich_TSS.80 <- readAll(ring.rich_TSS.80)
  ring.rich_TSS.95 <- readAll(ring.rich_TSS.95)
  
  save(ring.rich_Jaccard.80, file = paste0("./outputs/Mimicry_ring_richness/Ring_rich_Jaccard.80_",ring,".RData"), version = "2")
  saveRDS(ring.rich_Jaccard.80, file = paste0("./outputs/Mimicry_ring_richness/Ring_rich_Jaccard.80_",ring,".rds"), version = "2")
  save(ring.rich_Jaccard.95, file = paste0("./outputs/Mimicry_ring_richness/Ring_rich_Jaccard.95_",ring,".RData"), version = "2")
  saveRDS(ring.rich_Jaccard.95, file = paste0("./outputs/Mimicry_ring_richness/Ring_rich_Jaccard.95_",ring,".rds"), version = "2")
  save(ring.rich_TSS.80, file = paste0("./outputs/Mimicry_ring_richness/Ring_rich_TSS.80_",ring,".RData"), version = "2")
  saveRDS(ring.rich_TSS.80, file = paste0("./outputs/Mimicry_ring_richness/Ring_rich_TSS.80_",ring,".rds"), version = "2")
  save(ring.rich_TSS.95, file = paste0("./outputs/Mimicry_ring_richness/Ring_rich_TSS.95_",ring,".RData"), version = "2")
  saveRDS(ring.rich_TSS.95, file = paste0("./outputs/Mimicry_ring_richness/Ring_rich_TSS.95_",ring,".rds"), version = "2")
  
  # Multiple pages pdf with all Ring richness maps ####
  pdf(file = paste0("./maps/Mimicry_ring_richness/Ring_rich_maps_",ring,".pdf"), height = 6, width = 7)
  plot(ring.rich_Jaccard.80, main = paste0(ring,"\nJaccard.80"))
  plot(ring.rich_Jaccard.95, main = paste0(ring,"\nJaccard.95"))
  plot(ring.rich_TSS.80, main = paste0(ring,"\nTSS.80"))
  plot(ring.rich_TSS.95, main = paste0(ring,"\nTSS.95"))
  dev.off()
  
  print(i)
  
}

# 4/ Generate a single pdf with all Mimicry ring richness for Jaccard.80 ####

# 4 rings per pages = 11 pages. Alphabetic order

# New color palette
pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")

# Load map stuff
crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".rds"))
bg_mask <- readRDS(file = "./input_data/Map_stuff/bg_mask.rds")


# Plot
pdf(file = paste0("./maps/Mimicry_ring_richness/All_ring_richness_Jaccard.80.pdf"), height = 13, width = 13)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
par(mfrow = c(2, 2))

for (i in mimicry.list)
{
  
  ring.rich_Jaccard.80 <- readRDS(file = paste0("./outputs/Mimicry_ring_richness/Ring_rich_Jaccard.80_",i,".rds"))
  
  # Extract info for legend ticks
  max <- ring.rich_Jaccard.80@data@max
  min <- ring.rich_Jaccard.80@data@min
  
  image(ring.rich_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mimicry ring richness\n",i), 
        cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
        ylab = "", xlab = "",
        legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
        legend  = F)
  plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
  # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
  # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
  plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
  scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
  prettymapr::addnortharrow(scale = 0.7, padin = c(0.2, 0.2), text.col = "#00000000")
  addRasterLegend(ring.rich_Jaccard.80, locs = axisTicks(usr = c(min, max), log = F, nint = 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
  graphics::text(x = -111, y = 4, font = 2, cex = 1.3, label = "Species")
  
}

par(mar = internal_margins)
par(mfrow = c(1, 1))

dev.off()



#### Copy-paste the section to compute stack of mimicry rings maps in Script 14a, Sections 5.1 & 6.1

# More logical to compute this here
