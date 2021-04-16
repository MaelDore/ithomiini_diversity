
##### Script 14a: Compute indices maps #####

## Compute 22 biodiversity index maps

# 1/ Species richness
# 2/ Species Shannon's Diversity
# 3/ Species Shannon's Diversity with compatibility with species richness: (a) Raw and (b) with Jost's transformation (exponential)
# 4/ Species Evenness

# 5/ Mimicry richness
# 6/ Mimicry Shannon's Diversity
# 7/ Mimicry Shannon's Diversity with compatibility with mimicry richness: (a) Raw and (b) with Jost's transformation (exponential)
# 8/ Mimicry Evenness

# 9/  Continuous Species Rarity = Range-size Weighted Species Richness. (a) Full and (b) mean
# 10/ Categorical Species Rarity (25% threshold) : (a) Total and (b) Proportion
# 11/ Continuous Mimicry rarity = Range size weighted mimicry richness. (a) Full and (b) mean

# 12/ Community vulnerability (Ring mean vulnerability)

### 13/ Phylogeny-based indices

# 14/ Mean pairwise Phylogenetic Distance: (a) Full and (b) Andes
# 15/ Faith's Phylogenetic Diversity
# 16/ Equal-Splits & Fair-proportion : (a) Sum and (b) mean (x2)


###

# Inputs:
   # Stack of species 'probability' of presence from script 12
   # Stack of mimicry ring 'probability' of presence from script 13
   # Stack of mimicry ring richness of presence from script 13
   # Phylogeny

# Outputs:
   # Lots of biodiversity index maps. One for each 4 options between Jaccard/TSS and Buffer.80/95
   # Put multiple pdf in main folder and individual maps in subfolder such as /TSS_Buffer.80/

###


### Prepare stuff

# Effacer l'environnement
rm(list = ls())

library(raster)
library(rangeBuilder)
library(gplots)

# Change temp folder for raster files
rasterOptions(tmpdir = paste0("./temp"))
# To clean regularly temp folder
unlink(list.files(path = rasterOptions()$tmpdir, all.files = T, recursive = T, include.dirs = T, full.names = T), force = T, recursive = T) # Clean raster store in temp


# Load Summary table for OMU/unit and for Species
load(file = paste0("./input_data/list.sp.RData"))
load(file = paste0("./input_data/list.models.RData"))

# Mimicry list
mimicry.list <- as.character(unique(list.models$Mimicry.model)) # 44 Mimicry rings

# Make a list of modeled and non-modeled units (occ.unit = with only occurrences points)
modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]
rasterized_OMU <- list.models[is.na(list.models$Model_ID), ]

# Color palette for plot
# cool = rainbow(49, s = 1, v = 1, start=rgb2hsv(col2rgb('yellow'))[1], end=rgb2hsv(col2rgb('blue'))[1])
# warm = rainbow(50, s = 1, v= 1, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
# pal_bl_red  = c(rev(cool), rev(warm))
# pal_bl_red <- c(gplots::col2hex("grey93"), pal_bl_red)

# New color palette
pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")


# Choose resolution
res <- "15"

# Load mask for continent borders
continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_", res, ".rds"))
# crop_mask_shp <- rasterToPolygons(x = continent_mask, dissolve = T, digits = 4)
# save(crop_mask_shp, file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".RData"), version = "2")
# saveRDS(crop_mask_shp, file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".rds"), version = "2")
crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".rds"))

# Merra_clim_5_neotrop <- readRDS(file = "./input_data/MERRA_Clim/MERRAClim_5m_cropped.rds")
# continent_mask_5min <- Merra_clim_5_neotrop[[1]]
# continent_mask_5min <- calc(continent_mask_5min, fun = function(x, ...) {x*0}, na.rm = F)
# crop_mask_shp_5min <- rasterToPolygons(x = continent_mask_5min, dissolve = T, digits = 4)

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


###### 1/ Species richness ######

### 1.1/ Species Stack generation ####
All_sp_proba_stack_Jaccard.80 <- All_sp_proba_stack_Jaccard.95 <- All_sp_proba_stack_TSS.80 <- All_sp_proba_stack_TSS.95  <- stack(continent_mask) # 1e temp layer used to initiate stack, to remove afterwards
# i <- 1

for (i in 1:nrow(list.sp)) # By sp
{ 
  # Load sp name
  sp <- as.character(list.sp$Sp_full[i]) 
  
  # Load Sp continuous probability map 
  sp.cont_Jaccard.80 <- readRDS(file = paste0("./outputs/By_species/",sp,"/cont_Jaccard.80_",sp,".rds"))
  sp.cont_Jaccard.95 <- readRDS(file = paste0("./outputs/By_species/",sp,"/cont_Jaccard.95_",sp,".rds"))
  sp.cont_TSS.80 <- readRDS(file = paste0("./outputs/By_species/",sp,"/cont_TSS.80_",sp,".rds"))
  sp.cont_TSS.95 <- readRDS(file = paste0("./outputs/By_species/",sp,"/cont_TSS.95_",sp,".rds"))
  
  # Build stack
  All_sp_proba_stack_Jaccard.80 <- addLayer(All_sp_proba_stack_Jaccard.80, sp.cont_Jaccard.80)
  All_sp_proba_stack_Jaccard.95 <- addLayer(All_sp_proba_stack_Jaccard.95, sp.cont_Jaccard.95)
  All_sp_proba_stack_TSS.80 <- addLayer(All_sp_proba_stack_TSS.80, sp.cont_TSS.80)
  All_sp_proba_stack_TSS.95 <- addLayer(All_sp_proba_stack_TSS.95, sp.cont_TSS.95)
  
  if (i %% 10 == 0) {print(i)}
}

# Drop the useless first layer used to initiate the stack
All_sp_proba_stack_Jaccard.80 <- dropLayer(All_sp_proba_stack_Jaccard.80, i = 1)
All_sp_proba_stack_Jaccard.95 <- dropLayer(All_sp_proba_stack_Jaccard.95, i = 1)
All_sp_proba_stack_TSS.80 <- dropLayer(All_sp_proba_stack_TSS.80, i = 1)
All_sp_proba_stack_TSS.95 <- dropLayer(All_sp_proba_stack_TSS.95, i = 1)

# Name layers with sp names
names(All_sp_proba_stack_Jaccard.80) <- names(All_sp_proba_stack_Jaccard.95) <- names(All_sp_proba_stack_TSS.80) <- names(All_sp_proba_stack_TSS.95)<- list.sp$Sp_full

# nlayers(All_sp_proba_stack_Jaccard.95) # 388 species in the final stack

# Save stacks
save(All_sp_proba_stack_Jaccard.80, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.RData"), version = "2")
saveRDS(All_sp_proba_stack_Jaccard.80, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"), version = "2")
save(All_sp_proba_stack_Jaccard.95, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.95.RData"), version = "2")
saveRDS(All_sp_proba_stack_Jaccard.95, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.95.rds"), version = "2")
save(All_sp_proba_stack_Jaccard.80, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80.RData"), version = "2")
saveRDS(All_sp_proba_stack_Jaccard.80, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80.rds"), version = "2")
save(All_sp_proba_stack_TSS.95, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.95.RData"), version = "2")
saveRDS(All_sp_proba_stack_TSS.95, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.95.rds"), version = "2")

### Load directly the complete stack
All_sp_proba_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))
All_sp_proba_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.95.rds"))
All_sp_proba_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80.rds"))
All_sp_proba_stack_TSS.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.95.rds"))


### 1.2/ Index computation ####
tot.sp.richness_Jaccard.80 <- readAll(calc(All_sp_proba_stack_Jaccard.80, fun = sum))
tot.sp.richness_Jaccard.95 <- readAll(calc(All_sp_proba_stack_Jaccard.95, fun = sum))
tot.sp.richness_TSS.80 <- readAll(calc(All_sp_proba_stack_TSS.80, fun = sum))
tot.sp.richness_TSS.95 <- readAll(calc(All_sp_proba_stack_TSS.95, fun = sum))


# Save
save(tot.sp.richness_Jaccard.80, file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.RData"), version = "2")
saveRDS(tot.sp.richness_Jaccard.80, file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"), version = "2")
save(tot.sp.richness_Jaccard.95, file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.95.RData"), version = "2")
saveRDS(tot.sp.richness_Jaccard.95, file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.95.rds"), version = "2")
save(tot.sp.richness_TSS.80, file = paste0("./outputs/Indices_maps/tot.sp.richness_TSS.80.RData"), version = "2")
saveRDS(tot.sp.richness_TSS.80, file = paste0("./outputs/Indices_maps/tot.sp.richness_TSS.80.rds"), version = "2")
save(tot.sp.richness_TSS.95, file = paste0("./outputs/Indices_maps/tot.sp.richness_TSS.95.RData"), version = "2")
saveRDS(tot.sp.richness_TSS.95, file = paste0("./outputs/Indices_maps/tot.sp.richness_TSS.95.rds"), version = "2")


### Load directly the final Species richness layer
tot.sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
tot.sp.richness_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.95.rds"))
tot.sp.richness_TSS.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_TSS.80.rds"))
tot.sp.richness_TSS.95 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_TSS.95.rds"))

### 1.3/ Plot Species richness ####

### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/tot.sp.richness_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.sp.richness_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Species richness \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_Jaccard.80, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/tot.sp.richness_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.sp.richness_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Species richness \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_Jaccard.95, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/tot.sp.richness_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.sp.richness_TSS.80, col = pal_bl_red_Mannion, main = paste0("Species richness \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_TSS.80, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/tot.sp.richness_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.sp.richness_TSS.95, col = pal_bl_red_Mannion, main = paste0("Species richness \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_TSS.95, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(tot.sp.richness_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/tot.sp.richness_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(tot.sp.richness_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Species richness \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_Jaccard.80, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")

image(tot.sp.richness_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Species richness \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_Jaccard.95, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")

image(tot.sp.richness_TSS.80, col = pal_bl_red_Mannion, main = paste0("Species richness \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_TSS.80, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")

image(tot.sp.richness_TSS.95, col = pal_bl_red_Mannion, main = paste0("Species richness \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_TSS.95, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")

par(mar = internal_margins)

dev.off()


# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])



###### 2/ Shannon's diversity index ######

# Function to compute Shannon's diversity index using habitat suitability as proxy of abundance (because we don't have a better option)
shannon = function(x, na.rm) {
  x <- x[x>0] # Remove 0 values to avoid error for log(0)
  y <- x/sum(x) # Compute frequencies of each OMU
  h <- -sum(y*log(y)) # Compute Shannon's H'
  return(h) # Output
}

### Load directly the complete stacks of species probabilities/suitability maps
All_sp_proba_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))
All_sp_proba_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.95.rds"))
All_sp_proba_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80.rds"))
All_sp_proba_stack_TSS.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.95.rds"))

### 2.1/ Index computation ####
sp.diversity_Jaccard.80 <- calc(All_sp_proba_stack_Jaccard.80, fun = shannon)*1
sp.diversity_Jaccard.95 <- calc(All_sp_proba_stack_Jaccard.95, fun = shannon)*1
sp.diversity_TSS.80 <- calc(All_sp_proba_stack_TSS.80, fun = shannon)*1
sp.diversity_TSS.95 <- calc(All_sp_proba_stack_TSS.95, fun = shannon)*1

# Save
save(sp.diversity_Jaccard.80, file = paste0("./outputs/Indices_maps/sp.diversity_Jaccard.80.RData"), version = "2")
saveRDS(sp.diversity_Jaccard.80, file = paste0("./outputs/Indices_maps/sp.diversity_Jaccard.80.rds"), version = "2")
save(sp.diversity_Jaccard.95, file = paste0("./outputs/Indices_maps/sp.diversity_Jaccard.95.RData"), version = "2")
saveRDS(sp.diversity_Jaccard.95, file = paste0("./outputs/Indices_maps/sp.diversity_Jaccard.95.rds"), version = "2")
save(sp.diversity_TSS.80, file = paste0("./outputs/Indices_maps/sp.diversity_TSS.80.RData"), version = "2")
saveRDS(sp.diversity_TSS.80, file = paste0("./outputs/Indices_maps/sp.diversity_TSS.80.rds"), version = "2")
save(sp.diversity_TSS.95, file = paste0("./outputs/Indices_maps/sp.diversity_TSS.95.RData"), version = "2")
saveRDS(sp.diversity_TSS.95, file = paste0("./outputs/Indices_maps/sp.diversity_TSS.95.rds"), version = "2")


### Load directly the final Species diversity layers
sp.diversity_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/sp.diversity_Jaccard.80.rds"))
sp.diversity_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_maps/sp.diversity_Jaccard.95.rds"))
sp.diversity_TSS.80 <- readRDS(file = paste0("./outputs/Indices_maps/sp.diversity_TSS.80.rds"))
sp.diversity_TSS.95 <- readRDS(file = paste0("./outputs/Indices_maps/sp.diversity_TSS.95.rds"))

### 2.2/ Plot Species Diversity raw version ####

### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/sp_diversity_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Shannon's species diversity \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_Jaccard.80, locs = seq(0, 5, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/sp_diversity_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Shannon's species diversity \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_Jaccard.95, locs = seq(0, 5, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/sp_diversity_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity_TSS.80, col = pal_bl_red_Mannion, main = paste0("Shannon's species diversity \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_TSS.80, locs = seq(0, 5, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/sp_diversity_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity_TSS.95, col = pal_bl_red_Mannion, main = paste0("Shannon's species diversity \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_TSS.95, locs = seq(0, 5, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# ## Tmap version
# 
# library(tmap)
# 
# tm_shape(sp.diversity_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")


### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/sp.diversity_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.diversity_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Shannon's species diversity \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_Jaccard.80, locs = seq(0, 5, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(sp.diversity_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Shannon's species diversity \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_Jaccard.95, locs = seq(0, 5, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(sp.diversity_TSS.80, col = pal_bl_red_Mannion, main = paste0("Shannon's species diversity \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_TSS.80, locs = seq(0, 5, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(sp.diversity_TSS.95, col = pal_bl_red_Mannion, main = paste0("Shannon's species diversity \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_TSS.95, locs = seq(0, 5, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

par(mar = internal_margins)

dev.off()


# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])


###### 3/ Species Shannon's Diversity with compatibility with species richness #####

# Use Jost's transformation to make indices comparable. Select a similar number of OMU per pixel to avoid getting transformed values higher than actual species richness

# Corrected Shannon's diversity index function
# Limit number of OMU integrated in the computation, corresponding to the local estimated species richness
shannon_compatible = function(x, na.rm) { # Déclaration des arguments en input
  rich <- sum(x) # Local estimated species richness
  
  if (is.na(rich)) {
    h <- NA # Case with NA
  } else {
    
    if (round(rich, 0) == 0) {
      h <- 0 # Case with no OMU. Skip computation to avoid error such as log(0) = Inf
      
    } else { # Regular case
      
      x <- x[order(x, decreasing = T)] # Order the proba. by decreasing values
      x <- head(x, n = round(rich, 0)) # Extract only the N highest proba. N = Rounded estimated richness
      y <- x/sum(x) # Calcul des pi => passage en fréquences
      h <- -sum(y*log(y)) # Calcul du H'
    }
  }
  return(h) # Output
}

### Load directly the complete stacks of species probabilities/suitability maps
All_sp_proba_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))
All_sp_proba_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.95.rds"))
All_sp_proba_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80.rds"))
All_sp_proba_stack_TSS.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.95.rds"))


### 3.1/ Index computation ####
sp.diversity.compatible_Jaccard.80 <- calc(All_sp_proba_stack_Jaccard.80, fun = shannon_compatible)*1 # Raw version
sp.diversity.compatible_Jost_Jaccard.80 <- exp(sp.diversity.compatible_Jaccard.80) - 1 # Jost's version

sp.diversity.compatible_Jaccard.95 <- calc(All_sp_proba_stack_Jaccard.95, fun = shannon_compatible)*1 # Raw version
sp.diversity.compatible_Jost_Jaccard.95 <- exp(sp.diversity.compatible_Jaccard.95) - 1 # Jost's version

sp.diversity.compatible_TSS.80 <- calc(All_sp_proba_stack_TSS.80, fun = shannon_compatible)*1 # Raw version
sp.diversity.compatible_Jost_TSS.80 <- exp(sp.diversity.compatible_TSS.80) - 1 # Jost's version

sp.diversity.compatible_TSS.95 <- calc(All_sp_proba_stack_TSS.95, fun = shannon_compatible)*1 # Raw version
sp.diversity.compatible_Jost_TSS.95 <- exp(sp.diversity.compatible_TSS.95) - 1 # Jost's version

# Save
save(sp.diversity.compatible_Jaccard.80, file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jaccard.80.RData"), version = "2")
saveRDS(sp.diversity.compatible_Jaccard.80, file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jaccard.80.rds"), version = "2")
save(sp.diversity.compatible_Jaccard.95, file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jaccard.95.RData"), version = "2")
saveRDS(sp.diversity.compatible_Jaccard.95, file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jaccard.95.rds"), version = "2")
save(sp.diversity.compatible_TSS.80, file = paste0("./outputs/Indices_maps/sp.diversity.compatible_TSS.80.RData"), version = "2")
saveRDS(sp.diversity.compatible_TSS.80, file = paste0("./outputs/Indices_maps/sp.diversity.compatible_TSS.80.rds"), version = "2")
save(sp.diversity.compatible_TSS.95, file = paste0("./outputs/Indices_maps/sp.diversity.compatible_TSS.95.RData"), version = "2")
saveRDS(sp.diversity.compatible_TSS.95, file = paste0("./outputs/Indices_maps/sp.diversity.compatible_TSS.95.rds"), version = "2")

save(sp.diversity.compatible_Jost_Jaccard.80, file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jost_Jaccard.80.RData"), version = "2")
saveRDS(sp.diversity.compatible_Jost_Jaccard.80, file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jost_Jaccard.80.rds"), version = "2")
save(sp.diversity.compatible_Jost_Jaccard.95, file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jost_Jaccard.95.RData"), version = "2")
saveRDS(sp.diversity.compatible_Jost_Jaccard.95, file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jost_Jaccard.95.rds"), version = "2")
save(sp.diversity.compatible_Jost_TSS.80, file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jost_TSS.80.RData"), version = "2")
saveRDS(sp.diversity.compatible_Jost_TSS.80, file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jost_TSS.80.rds"), version = "2")
save(sp.diversity.compatible_Jost_TSS.95, file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jost_TSS.95.RData"), version = "2")
saveRDS(sp.diversity.compatible_Jost_TSS.95, file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jost_TSS.95.rds"), version = "2")


### Load directly the final Species corrected diversity layers
sp.diversity.compatible_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jaccard.80.rds"))
sp.diversity.compatible_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jaccard.95.rds"))
sp.diversity.compatible_TSS.80 <- readRDS(file = paste0("./outputs/Indices_maps/sp.diversity.compatible_TSS.80.rds"))
sp.diversity.compatible_TSS.95 <- readRDS(file = paste0("./outputs/Indices_maps/sp.diversity.compatible_TSS.95.rds"))

sp.diversity.compatible_Jost_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jost_Jaccard.80.rds"))
sp.diversity.compatible_Jost_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jost_Jaccard.95.rds"))
sp.diversity.compatible_Jost_TSS.80 <- readRDS(file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jost_TSS.80.rds"))
sp.diversity.compatible_Jost_TSS.95 <- readRDS(file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jost_TSS.95.rds"))


### 3.2/ Plot Species Diversity Corrected version ####

### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/sp_diversity_compatible_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Shannon's species diversity \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jaccard.80, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/sp_diversity_compatible_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Shannon's species diversity \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jaccard.95, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/sp_diversity_compatible_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_TSS.80, col = pal_bl_red_Mannion, main = paste0("Shannon's species diversity \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_TSS.80, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/sp_diversity_compatible_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_TSS.95, col = pal_bl_red_Mannion, main = paste0("Shannon's species diversity \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_TSS.95, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# ## Tmap version
# 
# library(tmap)
# 
# tm_shape(sp.diversity.compatible_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")


### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/sp.diversity.compatible_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.diversity.compatible_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Shannon's species diversity \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jaccard.80, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(sp.diversity.compatible_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Shannon's species diversity \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jaccard.95, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(sp.diversity.compatible_TSS.80, col = pal_bl_red_Mannion, main = paste0("Shannon's species diversity \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_TSS.80, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(sp.diversity.compatible_TSS.95, col = pal_bl_red_Mannion, main = paste0("Shannon's species diversity \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_TSS.95, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

par(mar = internal_margins)

dev.off()

### 3.3/ Plot Species Diversity Jost version ####

### Individual plots

str(sp.diversity.compatible_Jost_Jaccard.80)

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/sp_diversity_compatible_Jost_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
image(sp.diversity.compatible_Jost_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Species diversity (Jost's effective species richness) \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_Jaccard.80, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/sp_diversity_compatible_Jost_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_Jost_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Species diversity (Jost's effective species richness) \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_Jaccard.95, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/sp_diversity_compatible_Jost_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_Jost_TSS.80, col = pal_bl_red_Mannion, main = paste0("Species diversity (Jost's effective species richness) \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_TSS.80, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/sp_diversity_compatible_Jost_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_Jost_TSS.95, col = pal_bl_red_Mannion, main = paste0("Species diversity (Jost's effective species richness) \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_TSS.95, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# ## Tmap version
# 
# library(tmap)
# 
# tm_shape(sp.diversity.compatible_Jost_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")


### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/sp.diversity.compatible_Jost_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.diversity.compatible_Jost_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Species diversity (Jost's effective species richness) \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_Jaccard.80, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1, label = "Species")

image(sp.diversity.compatible_Jost_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Species diversity (Jost's effective species richness) \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_Jaccard.95, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1, label = "Species")

image(sp.diversity.compatible_Jost_TSS.80, col = pal_bl_red_Mannion, main = paste0("Species diversity (Jost's effective species richness) \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_TSS.80, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1, label = "Species")

image(sp.diversity.compatible_Jost_TSS.95, col = pal_bl_red_Mannion, main = paste0("Species diversity (Jost's effective species richness) \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_TSS.95, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1, label = "Species")

par(mar = internal_margins)

dev.off()


# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])



###### 4/ Species Evenness #####

# Function to compute Pielou's evenness index using probabilities of presence as proxy for "abundances"
evenness = function(x, na.rm) { 
  x <- round(x*1000, digits = 0) # Round values to avoid J > 1
  if (sum(x, na.rm = T)<1) { # Cannot be computed when sum of rounded "abundances" < 1
    J <- NA
  } else {
    x <- x[x>0] # Remove all 0 values to avoid error for log(0)
    y <- x/sum(x) # Compute frequencies
    H <- -sum(y*log(y)) # Compute H'
    Hmax <- log(sum(x)) # Compute Hmax
    J <- round(min(H/Hmax,1), digits = 3) # Compute J
  }
  return(J) # Output
}

### Load directly the complete stack of species probabilities of presence
All_sp_proba_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))
All_sp_proba_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.95.rds"))
All_sp_proba_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80.rds"))
All_sp_proba_stack_TSS.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.95.rds"))

### 4.1/ Index computation ####
sp.evenness_Jaccard.80 <- calc(All_sp_proba_stack_Jaccard.80, fun = evenness)*1
sp.evenness_Jaccard.95 <- calc(All_sp_proba_stack_Jaccard.95, fun = evenness)*1
sp.evenness_TSS.80 <- calc(All_sp_proba_stack_TSS.80, fun = evenness)*1
sp.evenness_TSS.95 <- calc(All_sp_proba_stack_TSS.95, fun = evenness)*1

# Save
save(sp.evenness_Jaccard.80, file = paste0("./outputs/Indices_maps/sp.evenness_Jaccard.80.RData"), version = "2")
saveRDS(sp.evenness_Jaccard.80, file = paste0("./outputs/Indices_maps/sp.evenness_Jaccard.80.rds"), version = "2")
save(sp.evenness_Jaccard.95, file = paste0("./outputs/Indices_maps/sp.evenness_Jaccard.95.RData"), version = "2")
saveRDS(sp.evenness_Jaccard.95, file = paste0("./outputs/Indices_maps/sp.evenness_Jaccard.95.rds"), version = "2")
save(sp.evenness_TSS.80, file = paste0("./outputs/Indices_maps/sp.evenness_TSS.80.RData"), version = "2")
saveRDS(sp.evenness_TSS.80, file = paste0("./outputs/Indices_maps/sp.evenness_TSS.80.rds"), version = "2")
save(sp.evenness_TSS.95, file = paste0("./outputs/Indices_maps/sp.evenness_TSS.95.RData"), version = "2")
saveRDS(sp.evenness_TSS.95, file = paste0("./outputs/Indices_maps/sp.evenness_TSS.95.rds"), version = "2")

### 4.2/ Contrast maps ####

# Load directly the final Species evenness layer
sp.evenness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/sp.evenness_Jaccard.80.rds"))
sp.evenness_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_maps/sp.evenness_Jaccard.95.rds"))
sp.evenness_TSS.80 <- readRDS(file = paste0("./outputs/Indices_maps/sp.evenness_TSS.80.rds"))
sp.evenness_TSS.95 <- readRDS(file = paste0("./outputs/Indices_maps/sp.evenness_TSS.95.rds"))

# Contrast values by merging low values below 0.25
sp.evenness_Jaccard.80_contrasted <- sp.evenness_Jaccard.80
sp.evenness_Jaccard.80_contrasted[sp.evenness_Jaccard.80 < 0.25] <- 0.255
sp.evenness_Jaccard.95_contrasted <- sp.evenness_Jaccard.95
sp.evenness_Jaccard.95_contrasted[sp.evenness_Jaccard.95 < 0.25] <- 0.255
sp.evenness_TSS.80_contrasted <- sp.evenness_TSS.80
sp.evenness_TSS.80_contrasted[sp.evenness_TSS.80 < 0.25] <- 0.255
sp.evenness_TSS.95_contrasted <- sp.evenness_TSS.95
sp.evenness_TSS.95_contrasted[sp.evenness_TSS.95 < 0.25] <- 0.255

# Add 0 values for empty continental pixels (transformed into min values as 0.25)
temp <- continent_mask
temp[!is.na(sp.evenness_Jaccard.80_contrasted[])] <- sp.evenness_Jaccard.80_contrasted[!is.na(sp.evenness_Jaccard.80_contrasted[])]
sp.evenness_Jaccard.80_contrasted <- temp
sp.evenness_Jaccard.80_contrasted[sp.evenness_Jaccard.80_contrasted == 0] <- 0.25

temp <- continent_mask
temp[!is.na(sp.evenness_Jaccard.95_contrasted[])] <- sp.evenness_Jaccard.95_contrasted[!is.na(sp.evenness_Jaccard.95_contrasted[])]
sp.evenness_Jaccard.95_contrasted <- temp
sp.evenness_Jaccard.95_contrasted[sp.evenness_Jaccard.95_contrasted == 0] <- 0.25

temp <- continent_mask
temp[!is.na(sp.evenness_TSS.80_contrasted[])] <- sp.evenness_TSS.80_contrasted[!is.na(sp.evenness_TSS.80_contrasted[])]
sp.evenness_TSS.80_contrasted <- temp
sp.evenness_TSS.80_contrasted[sp.evenness_TSS.80_contrasted == 0] <- 0.25

temp <- continent_mask
temp[!is.na(sp.evenness_TSS.95_contrasted[])] <- sp.evenness_TSS.95_contrasted[!is.na(sp.evenness_TSS.95_contrasted[])]
sp.evenness_TSS.95_contrasted <- temp
sp.evenness_TSS.95_contrasted[sp.evenness_TSS.95_contrasted == 0] <- 0.25

# hist(sample(sp.evenness_Jaccard.80, 1000))

### 4.3/ Plot with contrasted scale merging low values (0.25 = 0-0.25) ####

### Individual plots for contrasted scale

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/sp.evenness_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.evenness_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Species evenness \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_Jaccard.80_contrasted, locs = seq(0.25, 0.4, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/sp.evenness_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.evenness_Jaccard.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Species evenness \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="         Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_Jaccard.95_contrasted, locs = seq(0.25, 0.45, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/sp.evenness_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.evenness_TSS.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Species evenness \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="         Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_TSS.80_contrasted, locs = seq(0.25, 0.4, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.80/sp.evenness_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.evenness_TSS.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Species evenness \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_TSS.95_contrasted, locs = seq(0.25, 0.45, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(sp.evenness_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/sp.evenness_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.evenness_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Species evenness \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_Jaccard.80_contrasted, locs = seq(0.25, 0.4, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")

image(sp.evenness_Jaccard.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Species evenness \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="         Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_Jaccard.95_contrasted, locs = seq(0.25, 0.45, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")

image(sp.evenness_TSS.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Species evenness \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="         Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_TSS.80_contrasted, locs = seq(0.25, 0.4, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")

image(sp.evenness_TSS.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Species evenness \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_TSS.95_contrasted, locs = seq(0.25, 0.45, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")

par(mar = internal_margins)

dev.off()


# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])






###### 5/ Mimicry richness ######

# Mimicry list
mimicry.list <- as.character(unique(list.models$Mimicry.model)) # 44 Mimicry rings

### 5.1/ Mimicry Stack generation ####

All_ring_proba_stack_Jaccard.80 <- All_ring_proba_stack_Jaccard.95 <- All_ring_proba_stack_TSS.80 <- All_ring_proba_stack_TSS.95  <- stack(continent_mask) # 1e temp layer used to initiate stack, to remove afterwards
# i <- 1

for (i in 1:length(mimicry.list)) # By ring
{ 
  # i <- 10
  
  # Load ring name
  ring <- mimicry.list[i]
  
  # Load Sp continuous probability map 
  ring.cont_Jaccard.80 <- readRDS(file = paste0("./outputs/Mimicry_rings_proba/Ring_proba_Jaccard.80_",ring,".rds"))
  ring.cont_Jaccard.95 <- readRDS(file = paste0("./outputs/Mimicry_rings_proba/Ring_proba_Jaccard.95_",ring,".rds"))
  ring.cont_TSS.80 <- readRDS(file = paste0("./outputs/Mimicry_rings_proba/Ring_proba_TSS.80_",ring,".rds"))
  ring.cont_TSS.95 <- readRDS(file = paste0("./outputs/Mimicry_rings_proba/Ring_proba_TSS.95_",ring,".rds"))
  
  # Build stack
  All_ring_proba_stack_Jaccard.80 <- addLayer(All_ring_proba_stack_Jaccard.80, ring.cont_Jaccard.80)
  All_ring_proba_stack_Jaccard.95 <- addLayer(All_ring_proba_stack_Jaccard.95, ring.cont_Jaccard.95)
  All_ring_proba_stack_TSS.80 <- addLayer(All_ring_proba_stack_TSS.80, ring.cont_TSS.80)
  All_ring_proba_stack_TSS.95 <- addLayer(All_ring_proba_stack_TSS.95, ring.cont_TSS.95)
  
  if (i %% 10 == 0) {print(i)}
}

# Drop the useless first layer used to initiate the stack
All_ring_proba_stack_Jaccard.80 <- dropLayer(All_ring_proba_stack_Jaccard.80, i = 1)
All_ring_proba_stack_Jaccard.95 <- dropLayer(All_ring_proba_stack_Jaccard.95, i = 1)
All_ring_proba_stack_TSS.80 <- dropLayer(All_ring_proba_stack_TSS.80, i = 1)
All_ring_proba_stack_TSS.95 <- dropLayer(All_ring_proba_stack_TSS.95, i = 1)

# Name layers with sp names
names(All_ring_proba_stack_Jaccard.80) <- names(All_ring_proba_stack_Jaccard.95) <- names(All_ring_proba_stack_TSS.80) <- names(All_ring_proba_stack_TSS.95)<- mimicry.list

# nlayers(All_ring_proba_stack_Jaccard.95) # 44 rings in the final stack

# Save stacks
save(All_ring_proba_stack_Jaccard.80, file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_Jaccard.80.RData"), version = "2")
saveRDS(All_ring_proba_stack_Jaccard.80, file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_Jaccard.80.RData"), version = "2")
save(All_ring_proba_stack_Jaccard.95, file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_Jaccard.95.RData"), version = "2")
saveRDS(All_ring_proba_stack_Jaccard.95, file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_Jaccard.95.RData"), version = "2")
save(All_ring_proba_stack_Jaccard.80, file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_TSS.80.RData"), version = "2")
saveRDS(All_ring_proba_stack_Jaccard.80, file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_TSS.80.RData"), version = "2")
save(All_ring_proba_stack_TSS.95, file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_TSS.95.RData"), version = "2")
saveRDS(All_ring_proba_stack_TSS.95, file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_TSS.95.RData"), version = "2")

### Load directly the complete stack
All_ring_proba_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_Jaccard.80.RData"))
All_ring_proba_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_Jaccard.95.RData"))
All_ring_proba_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_TSS.80.RData"))
All_ring_proba_stack_TSS.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_TSS.95.RData"))


### 5.2/ Index computation ####

ring.richness_Jaccard.80 <- calc(All_ring_proba_stack_Jaccard.80, fun = sum)
ring.richness_Jaccard.95 <- calc(All_ring_proba_stack_Jaccard.95, fun = sum)
ring.richness_TSS.80 <- calc(All_ring_proba_stack_TSS.80, fun = sum)
ring.richness_TSS.95 <- calc(All_ring_proba_stack_TSS.95, fun = sum)

# Save
save(ring.richness_Jaccard.80, file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.RData"), version = "2")
saveRDS(ring.richness_Jaccard.80, file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"), version = "2")
save(ring.richness_Jaccard.95, file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.95.RData"), version = "2")
saveRDS(ring.richness_Jaccard.95, file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.95.rds"), version = "2")
save(ring.richness_TSS.80, file = paste0("./outputs/Indices_maps/ring.richness_TSS.80.RData"), version = "2")
saveRDS(ring.richness_TSS.80, file = paste0("./outputs/Indices_maps/ring.richness_TSS.80.rds"), version = "2")
save(ring.richness_TSS.95, file = paste0("./outputs/Indices_maps/ring.richness_TSS.95.RData"), version = "2")
saveRDS(ring.richness_TSS.95, file = paste0("./outputs/Indices_maps/ring.richness_TSS.95.rds"), version = "2")


### Load directly the final Ring richness layer
ring.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))
ring.richness_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.95.rds"))
ring.richness_TSS.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_TSS.80.rds"))
ring.richness_TSS.95 <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_TSS.95.rds"))

### 5.3/ Plot Mimicry richness ####


### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/ring.richness_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.richness_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mimicry richness \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mimicry\n         rings", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.richness_Jaccard.80, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Mimicry\nrings")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/ring.richness_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.richness_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Mimicry richness \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mimicry\n         rings", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.richness_Jaccard.95, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Mimicry\nrings")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/ring.richness_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.richness_TSS.80, col = pal_bl_red_Mannion, main = paste0("Mimicry richness \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mimicry\n         rings", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.richness_TSS.80, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Mimicry\nrings")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/ring.richness_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.richness_TSS.95, col = pal_bl_red_Mannion, main = paste0("Mimicry richness \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mimicry\n         rings", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.richness_TSS.95, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Mimicry\nrings")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(ring.richness_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/ring.richness_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(ring.richness_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mimicry richness \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mimicry\n         rings", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.richness_Jaccard.80, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Mimicry\nrings")

image(ring.richness_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Mimicry richness \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mimicry\n         rings", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.richness_Jaccard.95, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Mimicry\nrings")

image(ring.richness_TSS.80, col = pal_bl_red_Mannion, main = paste0("Mimicry richness \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mimicry\n         rings", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.richness_TSS.80, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Mimicry\nrings")

image(ring.richness_TSS.95, col = pal_bl_red_Mannion, main = paste0("Mimicry richness \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mimicry\n         rings", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.richness_TSS.95, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Mimicry\nrings")

par(mar = internal_margins)

dev.off()


# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])


###### 6/ Mimicry diversity #####

# Shannon's diversity index using nb of species in each ring as abundances of rings

### 6.1/ Generate Mimicry rings richness stack #####

All_ring_rich_stack_Jaccard.80 <- All_ring_rich_stack_Jaccard.95 <- All_ring_rich_stack_TSS.80 <- All_ring_rich_stack_TSS.95  <- stack(continent_mask) # 1e temp layer used to initiate stack, to remove afterwards

# i <- 1
for (i in 1:length(mimicry.list)) { # By ring
  
  ring <- mimicry.list[i] # Load ring name
  
  # Load ring richness map 
  ring.rich_Jaccard.80 <- readRDS(file = paste0("./outputs/Mimicry_ring_richness/Ring_rich_Jaccard.80_",ring,".rds"))
  ring.rich_Jaccard.95 <- readRDS(file = paste0("./outputs/Mimicry_ring_richness/Ring_rich_Jaccard.95_",ring,".rds"))
  ring.rich_TSS.80 <- readRDS(file = paste0("./outputs/Mimicry_ring_richness/Ring_rich_TSS.80_",ring,".rds"))
  ring.rich_TSS.95 <- readRDS(file = paste0("./outputs/Mimicry_ring_richness/Ring_rich_TSS.95_",ring,".rds"))
  
  # Build stack
  All_ring_rich_stack_Jaccard.80 <- addLayer(All_ring_rich_stack_Jaccard.80, ring.rich_Jaccard.80)
  All_ring_rich_stack_Jaccard.95 <- addLayer(All_ring_rich_stack_Jaccard.95, ring.rich_Jaccard.95)
  All_ring_rich_stack_TSS.80 <- addLayer(All_ring_rich_stack_TSS.80, ring.rich_TSS.80)
  All_ring_rich_stack_TSS.95 <- addLayer(All_ring_rich_stack_TSS.95, ring.rich_TSS.95)
  
  if (i %% 10 == 0) {print(i)}
}

# Drop the useless first layer used to initiate the stack
All_ring_rich_stack_Jaccard.80 <- dropLayer(All_ring_rich_stack_Jaccard.80, i = 1)
All_ring_rich_stack_Jaccard.95 <- dropLayer(All_ring_rich_stack_Jaccard.95, i = 1)
All_ring_rich_stack_TSS.80 <- dropLayer(All_ring_rich_stack_TSS.80, i = 1)
All_ring_rich_stack_TSS.95 <- dropLayer(All_ring_rich_stack_TSS.95, i = 1)

# Name layers with sp names
names(All_ring_rich_stack_Jaccard.80) <- names(All_ring_rich_stack_Jaccard.95) <- names(All_ring_rich_stack_TSS.80) <- names(All_ring_rich_stack_TSS.95)<- mimicry.list

# nlayers(All_ring_rich_stack_Jaccard.95) # 44 rings in the final stack

# Save stacks
save(All_ring_rich_stack_Jaccard.80, file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_Jaccard.80.RData"), version = "2")
saveRDS(All_ring_rich_stack_Jaccard.80, file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_Jaccard.80.RData"), version = "2")
save(All_ring_rich_stack_Jaccard.95, file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_Jaccard.95.RData"), version = "2")
saveRDS(All_ring_rich_stack_Jaccard.95, file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_Jaccard.95.RData"), version = "2")
save(All_ring_rich_stack_Jaccard.80, file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_TSS.80.RData"), version = "2")
saveRDS(All_ring_rich_stack_Jaccard.80, file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_TSS.80.RData"), version = "2")
save(All_ring_rich_stack_TSS.95, file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_TSS.95.RData"), version = "2")
saveRDS(All_ring_rich_stack_TSS.95, file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_TSS.95.RData"), version = "2")

### Load directly the complete stacks of mimicry ring richness
All_ring_rich_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_Jaccard.80.RData"))
All_ring_rich_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_Jaccard.95.RData"))
All_ring_rich_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_TSS.80.RData"))
All_ring_rich_stack_TSS.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_TSS.95.RData"))

### 6.2/ Index computation #####

# Function to compute Shannon's diversity index using habitat suitability as proxy of abundance (because we don't have a better option)
shannon = function(x, na.rm) {
  x <- x[x>0] # Remove 0 values to avoid error for log(0)
  y <- x/sum(x) # Compute frequencies of each OMU
  h <- -sum(y*log(y)) # Compute Shannon's H'
  return(h) # Output
}

ring.diversity_Jaccard.80 <- calc(All_ring_rich_stack_Jaccard.80, fun = shannon)*1
ring.diversity_Jaccard.95 <- calc(All_ring_rich_stack_Jaccard.95, fun = shannon)*1
ring.diversity_TSS.80 <- calc(All_ring_rich_stack_TSS.80, fun = shannon)*1
ring.diversity_TSS.95 <- calc(All_ring_rich_stack_TSS.95, fun = shannon)*1

# Save
save(ring.diversity_Jaccard.80, file = paste0("./outputs/Indices_maps/ring.diversity_Jaccard.80.RData"), version = "2")
saveRDS(ring.diversity_Jaccard.80, file = paste0("./outputs/Indices_maps/ring.diversity_Jaccard.80.rds"), version = "2")
save(ring.diversity_Jaccard.95, file = paste0("./outputs/Indices_maps/ring.diversity_Jaccard.95.RData"), version = "2")
saveRDS(ring.diversity_Jaccard.95, file = paste0("./outputs/Indices_maps/ring.diversity_Jaccard.95.rds"), version = "2")
save(ring.diversity_TSS.80, file = paste0("./outputs/Indices_maps/ring.diversity_TSS.80.RData"), version = "2")
saveRDS(ring.diversity_TSS.80, file = paste0("./outputs/Indices_maps/ring.diversity_TSS.80.rds"), version = "2")
save(ring.diversity_TSS.95, file = paste0("./outputs/Indices_maps/ring.diversity_TSS.95.RData"), version = "2")
saveRDS(ring.diversity_TSS.95, file = paste0("./outputs/Indices_maps/ring.diversity_TSS.95.rds"), version = "2")


### Load directly the final Species richness layer
ring.diversity_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.diversity_Jaccard.80.rds"))
ring.diversity_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_maps/ring.diversity_Jaccard.95.rds"))
ring.diversity_TSS.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.diversity_TSS.80.rds"))
ring.diversity_TSS.95 <- readRDS(file = paste0("./outputs/Indices_maps/ring.diversity_TSS.95.rds"))

### 6.3/ Plot mimicry diversity - Raw version ####

### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/ring_diversity_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.diversity_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          H' index", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity_Jaccard.80, locs = seq(0, 3, 0.5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/ring_diversity_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.diversity_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          H' index", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity_Jaccard.95, locs = seq(0, 3, 0.5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/ring_diversity_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.diversity_TSS.80, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          H' index", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity_TSS.80, locs = seq(0, 3, 0.5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/ring_diversity_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.diversity_TSS.95, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          H' index", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity_TSS.95, locs = seq(0, 3, 0.5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# ## Tmap version
# 
# library(tmap)
# 
# tm_shape(ring.diversity_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")


### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/ring.diversity_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(ring.diversity_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          H' index", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity_Jaccard.80, locs = seq(0, 3, 0.5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(ring.diversity_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          H' index", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity_Jaccard.95, locs = seq(0, 3, 0.5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(ring.diversity_TSS.80, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          H' index", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity_TSS.80, locs = seq(0, 3, 0.5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(ring.diversity_TSS.95, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          H' index", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity_TSS.95, locs = seq(0, 3, 0.5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

par(mar = internal_margins)

dev.off()


# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])


###### 7/ Mimicry Shannon's Diversity with compatibility with mimicry richness #####

# Use Jost's transformation to make indices comparable. Select a similar number of OMU per pixel to avoid getting transformed values higher than actual species richness

# Corrected Shannon's diversity index function
# Limit number of OMU integrated in the computation, corresponding to the local estimated mimicry ring richness
# Modified version where local richness is provided by another raster
shannon_compatible_rich = function(x, rich, na.rm) {
  if (is.na(rich)) {
    h <- NA # Case with NA
  } else {
    if (round(rich, 0) == 0) {
      h <- 0 # Case with no OMU. Skip computation to avoid error such as log(0) = Inf
      
    } else { # Regular case
      
      x <- x[order(x, decreasing = T)] # Order the proba. by decreasing values
      x <- head(x, n = round(rich, 0)) # Extract only the N highest proba. N = Rounded estimated richness
      y <- x/sum(x) # Compute frequencies
      h <- -sum(y*log(y)) # Compute H'
    }
    
  }
  return(h) # Output
}


### Load directly the complete stacks of mimicry ring richness
All_ring_rich_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_Jaccard.80.RData"))
All_ring_rich_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_Jaccard.95.RData"))
All_ring_rich_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_TSS.80.RData"))
All_ring_rich_stack_TSS.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_TSS.95.RData"))

### Load directly the final Ring richness layer
ring.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))
ring.richness_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.95.rds"))
ring.richness_TSS.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_TSS.80.rds"))
ring.richness_TSS.95 <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_TSS.95.rds"))

plot(ring.richness_Jaccard.80)

# Index computation via overlay (does not work)
# ring.diversity.compatible_Jaccard.80 <- overlay(All_ring_rich_stack_Jaccard.80, ring.richness_Jaccard.80, fun = shannon_compatible_rich)*1 # Raw version
# ring.diversity.compatible_Jost_Jaccard.80 <- exp(ring.diversity.compatible_Jaccard.80) - 1 # Jost's version
# 
# ring.diversity.compatible_Jaccard.95 <- overlay(All_ring_rich_stack_Jaccard.95, fun = shannon_compatible_rich)*1 # Raw version
# ring.diversity.compatible_Jost_Jaccard.95 <- exp(ring.diversity.compatible_Jaccard.95) - 1 # Jost's version
# 
# ring.diversity.compatible_TSS.80 <- overlay(All_ring_rich_stack_TSS.80, fun = shannon_compatible_rich)*1 # Raw version
# ring.diversity.compatible_Jost_TSS.80 <- exp(ring.diversity.compatible_TSS.80) - 1 # Jost's version
# 
# ring.diversity.compatible_TSS.95 <- overlay(All_ring_rich_stack_TSS.95, fun = shannon_compatible_rich)*1 # Raw version
# ring.diversity.compatible_Jost_TSS.95 <- exp(ring.diversity.compatible_TSS.95) - 1 # Jost's version

### 7.1/ Index computation ####

# Index computation via a loop through each pixel ! For Jaccard.80
ring.diversity.compatible_Jaccard.80_values <- NA
for (i in 1:length(ring.richness_Jaccard.80)) {
  ring.diversity.compatible_Jaccard.80_values[i] <- shannon_compatible_rich(x = All_ring_rich_stack_Jaccard.80[i], rich = ring.richness_Jaccard.80[i])
  if ((i %% 1000) == 0) {
    print(i)
  }
}
ring.diversity.compatible_Jaccard.80 <- ring.richness_Jaccard.80
ring.diversity.compatible_Jaccard.80@data@values <- ring.diversity.compatible_Jaccard.80_values
ring.diversity.compatible_Jost_Jaccard.80 <- exp(ring.diversity.compatible_Jaccard.80) - 1

# Index computation via a loop through each pixel ! For Jaccard.95
ring.diversity.compatible_Jaccard.95_values <- NA
for (i in 1:length(ring.richness_Jaccard.95)) {
  ring.diversity.compatible_Jaccard.95_values[i] <- shannon_compatible_rich(x = All_ring_rich_stack_Jaccard.95[i], rich = ring.richness_Jaccard.95[i])
  if ((i %% 1000) == 0) {
    print(i)
  }
}
ring.diversity.compatible_Jaccard.95 <- ring.richness_Jaccard.95
ring.diversity.compatible_Jaccard.95@data@values <- ring.diversity.compatible_Jaccard.95_values
ring.diversity.compatible_Jost_Jaccard.95 <- exp(ring.diversity.compatible_Jaccard.95) - 1

# Index computation via a loop through each pixel ! For TSS.80
ring.diversity.compatible_TSS.80_values <- NA
for (i in 1:length(ring.richness_TSS.80)) {
  ring.diversity.compatible_TSS.80_values[i] <- shannon_compatible_rich(x = All_ring_rich_stack_TSS.80[i], rich = ring.richness_TSS.80[i])
  if ((i %% 1000) == 0) {
    print(i)
  }
}
ring.diversity.compatible_TSS.80 <- ring.richness_TSS.80
ring.diversity.compatible_TSS.80@data@values <- ring.diversity.compatible_TSS.80_values
ring.diversity.compatible_Jost_TSS.80 <- exp(ring.diversity.compatible_TSS.80) - 1

# Index computation via a loop through each pixel ! For TSS.95
ring.diversity.compatible_TSS.95_values <- NA
for (i in 1:length(ring.richness_TSS.95)) {
  ring.diversity.compatible_TSS.95_values[i] <- shannon_compatible_rich(x = All_ring_rich_stack_TSS.95[i], rich = ring.richness_TSS.95[i])
  if ((i %% 1000) == 0) {
    print(i)
  }
}
ring.diversity.compatible_TSS.95 <- ring.richness_TSS.95
ring.diversity.compatible_TSS.95@data@values <- ring.diversity.compatible_TSS.95_values
ring.diversity.compatible_Jost_TSS.95 <- exp(ring.diversity.compatible_TSS.95) - 1

# Save
save(ring.diversity.compatible_Jaccard.80, file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jaccard.80.RData"), version = "2")
saveRDS(ring.diversity.compatible_Jaccard.80, file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jaccard.80.rds"), version = "2")
save(ring.diversity.compatible_Jaccard.95, file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jaccard.95.RData"), version = "2")
saveRDS(ring.diversity.compatible_Jaccard.95, file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jaccard.95.rds"), version = "2")
save(ring.diversity.compatible_TSS.80, file = paste0("./outputs/Indices_maps/ring.diversity.compatible_TSS.80.RData"), version = "2")
saveRDS(ring.diversity.compatible_TSS.80, file = paste0("./outputs/Indices_maps/ring.diversity.compatible_TSS.80.rds"), version = "2")
save(ring.diversity.compatible_TSS.95, file = paste0("./outputs/Indices_maps/ring.diversity.compatible_TSS.95.RData"), version = "2")
saveRDS(ring.diversity.compatible_TSS.95, file = paste0("./outputs/Indices_maps/ring.diversity.compatible_TSS.95.rds"), version = "2")

save(ring.diversity.compatible_Jost_Jaccard.80, file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jost_Jaccard.80.RData"), version = "2")
saveRDS(ring.diversity.compatible_Jost_Jaccard.80, file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jost_Jaccard.80.rds"), version = "2")
save(ring.diversity.compatible_Jost_Jaccard.95, file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jost_Jaccard.95.RData"), version = "2")
saveRDS(ring.diversity.compatible_Jost_Jaccard.95, file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jost_Jaccard.95.rds"), version = "2")
save(ring.diversity.compatible_Jost_TSS.80, file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jost_TSS.80.RData"), version = "2")
saveRDS(ring.diversity.compatible_Jost_TSS.80, file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jost_TSS.80.rds"), version = "2")
save(ring.diversity.compatible_Jost_TSS.95, file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jost_TSS.95.RData"), version = "2")
saveRDS(ring.diversity.compatible_Jost_TSS.95, file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jost_TSS.95.rds"), version = "2")


### Load directly the final Species corrected diversity layers
ring.diversity.compatible_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jaccard.80.rds"))
ring.diversity.compatible_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jaccard.95.rds"))
ring.diversity.compatible_TSS.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.diversity.compatible_TSS.80.rds"))
ring.diversity.compatible_TSS.95 <- readRDS(file = paste0("./outputs/Indices_maps/ring.diversity.compatible_TSS.95.rds"))

ring.diversity.compatible_Jost_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jost_Jaccard.80.rds"))
ring.diversity.compatible_Jost_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jost_Jaccard.95.rds"))
ring.diversity.compatible_Jost_TSS.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jost_TSS.80.rds"))
ring.diversity.compatible_Jost_TSS.95 <- readRDS(file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jost_TSS.95.rds"))


### 7.2/ Plot Mimicry Diversity - Raw version ####

### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/ring_diversity_compatible_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.diversity.compatible_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          H' index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity.compatible_Jaccard.80, locs = seq(0, 3, 0.5), minmax = c(0, max(ring.diversity.compatible_Jaccard.80[], na.rm = T)), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/ring_diversity_compatible_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.diversity.compatible_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          H' index", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity.compatible_Jaccard.95, locs = seq(0, 3, 0.5), minmax = c(0, max(ring.diversity.compatible_Jaccard.95[], na.rm = T)), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/ring_diversity_compatible_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.diversity.compatible_TSS.80, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          H' index", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity.compatible_Jaccard.95, locs = seq(0, 3, 0.5), minmax = c(0, max(ring.diversity.compatible_TSS.80[], na.rm = T)), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/ring_diversity_compatible_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.diversity.compatible_TSS.95, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          H' index", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity.compatible_Jaccard.95, locs = seq(0, 3, 0.5), minmax = c(0, max(ring.diversity.compatible_TSS.95[], na.rm = T)), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# ## Tmap version
# 
# library(tmap)
# 
# tm_shape(ring.diversity.compatible_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")


### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/ring.diversity.compatible_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(ring.diversity.compatible_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          H' index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity.compatible_Jaccard.80, locs = seq(0, 3, 0.5), minmax = c(0, max(ring.diversity.compatible_Jaccard.80[], na.rm = T)), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(ring.diversity.compatible_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          H' index", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity.compatible_Jaccard.95, locs = seq(0, 3, 0.5), minmax = c(0, max(ring.diversity.compatible_Jaccard.95[], na.rm = T)), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(ring.diversity.compatible_TSS.80, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          H' index", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity.compatible_Jaccard.95, locs = seq(0, 3, 0.5), minmax = c(0, max(ring.diversity.compatible_TSS.80[], na.rm = T)), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(ring.diversity.compatible_TSS.95, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          H' index", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity.compatible_Jaccard.95, locs = seq(0, 3, 0.5), minmax = c(0, max(ring.diversity.compatible_TSS.95[], na.rm = T)), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

par(mar = internal_margins)

dev.off()

### 7.3/ Plot Mimicry Diversity - Jost version ####

### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/ring_diversity_compatible_Jost_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
image(ring.diversity.compatible_Jost_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity (Jost's effective mimicry richness) \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity.compatible_Jost_Jaccard.80, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Mimicry\nrings")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/ring_diversity_compatible_Jost_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.diversity.compatible_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          H' index", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity.compatible_Jaccard.95, locs = seq(0, 20, 5), minmax = c(0, max(ring.diversity.compatible_Jaccard.95[], na.rm = T)), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Mimicry\nrings")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/ring_diversity_compatible_Jost_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.diversity.compatible_TSS.80, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          H' index", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity.compatible_Jaccard.95, locs = seq(0, 20, 5), minmax = c(0, max(ring.diversity.compatible_TSS.80[], na.rm = T)), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Mimicry\nrings")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/ring_diversity_compatible_Jost_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.diversity.compatible_Jost_TSS.95, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity (Jost's effective mimicry richness) \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity.compatible_Jost_TSS.95, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Mimicry\nrings")
par(mar = internal_margins)
dev.off()

# ## Tmap version
# 
# library(tmap)
# 
# tm_shape(ring.diversity.compatible_Jost_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")


### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/ring.diversity.compatible_Jost_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(ring.diversity.compatible_Jost_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity (Jost's effective mimicry richness) \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity.compatible_Jost_Jaccard.80, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Mimicry\nrings")

image(ring.diversity.compatible_Jost_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity (Jost's effective mimicry richness) \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity.compatible_Jost_Jaccard.95, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Mimicry\nrings")

image(ring.diversity.compatible_Jost_TSS.80, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity (Jost's effective mimicry richness) \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity.compatible_Jost_TSS.80, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Mimicry\nrings")

image(ring.diversity.compatible_Jost_TSS.95, col = pal_bl_red_Mannion, main = paste0("Mimicry diversity (Jost's effective mimicry richness) \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.diversity.compatible_Jost_TSS.95, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Mimicry\nrings")

par(mar = internal_margins)

dev.off()


# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])



###### 8/ Mimicry Evenness #####

# Function to compute Pielou's evenness index using probabilities of presence as proxy for "abundances"
evenness = function(x, na.rm) { 
  x <- round(x*1000, digits = 0) # Round values to avoid J > 1
  if (sum(x, na.rm = T)<1) { # Cannot be computed when sum of rounded "abundances" < 1
    J <- NA
  } else {
    x <- x[x>0] # Remove all 0 values to avoid error for log(0)
    y <- x/sum(x) # Compute frequencies
    H <- -sum(y*log(y)) # Compute H'
    Hmax <- log(sum(x)) # Compute Hmax
    J <- round(min(H/Hmax,1), digits = 3) # Compute J
  }
  return(J) # Output
}

### Load directly the complete stacks of mimicry ring richness
All_ring_rich_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_Jaccard.80.RData"))
All_ring_rich_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_Jaccard.95.RData"))
All_ring_rich_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_TSS.80.RData"))
All_ring_rich_stack_TSS.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_TSS.95.RData"))

### 8.1/ Index computation ####
ring.evenness_Jaccard.80 <- calc(All_ring_rich_stack_Jaccard.80, fun = evenness)*1
ring.evenness_Jaccard.95 <- calc(All_ring_rich_stack_Jaccard.95, fun = evenness)*1
ring.evenness_TSS.80 <- calc(All_ring_rich_stack_TSS.80, fun = evenness)*1
ring.evenness_TSS.95 <- calc(All_ring_rich_stack_TSS.95, fun = evenness)*1

# Save
save(ring.evenness_Jaccard.80, file = paste0("./outputs/Indices_maps/ring.evenness_Jaccard.80.RData"), version = "2")
saveRDS(ring.evenness_Jaccard.80, file = paste0("./outputs/Indices_maps/ring.evenness_Jaccard.80.rds"), version = "2")
save(ring.evenness_Jaccard.95, file = paste0("./outputs/Indices_maps/ring.evenness_Jaccard.95.RData"), version = "2")
saveRDS(ring.evenness_Jaccard.95, file = paste0("./outputs/Indices_maps/ring.evenness_Jaccard.95.rds"), version = "2")
save(ring.evenness_TSS.80, file = paste0("./outputs/Indices_maps/ring.evenness_TSS.80.RData"), version = "2")
saveRDS(ring.evenness_TSS.80, file = paste0("./outputs/Indices_maps/ring.evenness_TSS.80.rds"), version = "2")
save(ring.evenness_TSS.95, file = paste0("./outputs/Indices_maps/ring.evenness_TSS.95.RData"), version = "2")
saveRDS(ring.evenness_TSS.95, file = paste0("./outputs/Indices_maps/ring.evenness_TSS.95.rds"), version = "2")

### 8.2/ Contrasting maps ####

### Load directly the final Species evenness layer
ring.evenness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.evenness_Jaccard.80.rds"))
ring.evenness_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_maps/ring.evenness_Jaccard.95.rds"))
ring.evenness_TSS.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.evenness_TSS.80.rds"))
ring.evenness_TSS.95 <- readRDS(file = paste0("./outputs/Indices_maps/ring.evenness_TSS.95.rds"))

# hist(ring.evenness_Jaccard.80)

# Contrast values by merging low values below 0.10
ring.evenness_Jaccard.80_contrasted <- ring.evenness_Jaccard.80
ring.evenness_Jaccard.80_contrasted[ring.evenness_Jaccard.80 < 0.10] <- 0.105
ring.evenness_Jaccard.95_contrasted <- ring.evenness_Jaccard.95
ring.evenness_Jaccard.95_contrasted[ring.evenness_Jaccard.95 < 0.10] <- 0.105
ring.evenness_TSS.80_contrasted <- ring.evenness_TSS.80
ring.evenness_TSS.80_contrasted[ring.evenness_TSS.80 < 0.10] <- 0.105
ring.evenness_TSS.95_contrasted <- ring.evenness_TSS.95
ring.evenness_TSS.95_contrasted[ring.evenness_TSS.95 < 0.10] <- 0.105


# Add 0 values for empty continental pixels (transformed into min values as 0.10)
temp <- continent_mask
temp[!is.na(ring.evenness_Jaccard.80_contrasted[])] <- ring.evenness_Jaccard.80_contrasted[!is.na(ring.evenness_Jaccard.80_contrasted[])]
ring.evenness_Jaccard.80_contrasted <- temp
ring.evenness_Jaccard.80_contrasted[ring.evenness_Jaccard.80_contrasted == 0] <- 0.10

temp <- continent_mask
temp[!is.na(ring.evenness_Jaccard.95_contrasted[])] <- ring.evenness_Jaccard.95_contrasted[!is.na(ring.evenness_Jaccard.95_contrasted[])]
ring.evenness_Jaccard.95_contrasted <- temp
ring.evenness_Jaccard.95_contrasted[ring.evenness_Jaccard.95_contrasted == 0] <- 0.10

temp <- continent_mask
temp[!is.na(ring.evenness_TSS.80_contrasted[])] <- ring.evenness_TSS.80_contrasted[!is.na(ring.evenness_TSS.80_contrasted[])]
ring.evenness_TSS.80_contrasted <- temp
ring.evenness_TSS.80_contrasted[ring.evenness_TSS.80_contrasted == 0] <- 0.10

temp <- continent_mask
temp[!is.na(ring.evenness_TSS.95_contrasted[])] <- ring.evenness_TSS.95_contrasted[!is.na(ring.evenness_TSS.95_contrasted[])]
ring.evenness_TSS.95_contrasted <- temp
ring.evenness_TSS.95_contrasted[ring.evenness_TSS.95_contrasted == 0] <- 0.10

### 8.3/ Plot with contrasted scale merging low values (0.10 = 0-0.10) ####

### Individual plots for contrasted scale

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/ring.evenness_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.evenness_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Mimicry evenness \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.evenness_Jaccard.80_contrasted, locs = seq(0.10, 0.25, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/ring.evenness_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.evenness_Jaccard.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Mimicry evenness \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="         Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.evenness_Jaccard.95_contrasted, locs = seq(0.10, 0.25, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/ring.evenness_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.evenness_TSS.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Mimicry evenness \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="         Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.evenness_TSS.80_contrasted, locs = seq(0.10, 0.25, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/ring.evenness_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.evenness_TSS.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Mimicry evenness \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.evenness_TSS.95_contrasted, locs = seq(0.10, 0.25, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(ring.evenness_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/ring.evenness_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(ring.evenness_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Mimicry evenness \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.evenness_Jaccard.80_contrasted, locs = seq(0.10, 0.25, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")

image(ring.evenness_Jaccard.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Mimicry evenness \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="         Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.evenness_Jaccard.95_contrasted, locs = seq(0.10, 0.25, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")

image(ring.evenness_TSS.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Mimicry evenness \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="         Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.evenness_TSS.80_contrasted, locs = seq(0.10, 0.25, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")

image(ring.evenness_TSS.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Mimicry evenness \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ring.evenness_TSS.95_contrasted, locs = seq(0.10, 0.25, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")

par(mar = internal_margins)

dev.off()


# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])






###### 9/ Continuous Species Rarity = Range-size Weighted Species Richness ######

### Load directly the stack of sp probas use to compute rarity based on geographic range
All_sp_proba_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))
All_sp_proba_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.95.rds"))
All_sp_proba_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80.rds"))
All_sp_proba_stack_TSS.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.95.rds"))

# Load maps of Species richness to use to compute mean rarity index among local species
sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
sp.richness_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.95.rds"))
sp.richness_TSS.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_TSS.80.rds"))
sp.richness_TSS.95 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_TSS.95.rds"))

### 9.1/ Index computation ####

## For Jaccard.80

# Compute weights based on geographic ranges
sp.range_size_Jaccard.80 <- NA
for (i in 1:nlayers(All_sp_proba_stack_Jaccard.80)) {
  # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
  # Compute the estimated number of pixels occupied by the species as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
  sp.range_size_Jaccard.80[i] <- sum(All_sp_proba_stack_Jaccard.80[[i]]@data@values, na.rm = T)*774.51/1000
}
sp.rarity_indices_Jaccard.80 <- 1-(sp.range_size_Jaccard.80/max(sp.range_size_Jaccard.80)) # Weights by max species extent = rarity indices

# Apply weights
sp.weighted.stack_Jaccard.80 <- All_sp_proba_stack_Jaccard.80 # Generate the stack to fill
for (i in 1:nlayers(sp.weighted.stack_Jaccard.80)){
  # Multiply species probability of presence by species rarity indices
  sp.weighted.stack_Jaccard.80@layers[[i]]@data@values <- All_sp_proba_stack_Jaccard.80@layers[[i]]@data@values * sp.rarity_indices_Jaccard.80[i]
}

# Sum = Richness weighted by rarity based on Range-size
sp.rarity_Jaccard.80 <- readAll(calc(sp.weighted.stack_Jaccard.80, fun = sum)) 
# Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each species)
sp.mean.rarity_Jaccard.80 <- sp.rarity_Jaccard.80/sp.richness_Jaccard.80

# Add 0 values for empty continental pixels (transformed into min values as 0.10)
temp <- continent_mask
temp[!is.na(sp.mean.rarity_Jaccard.80[])] <- sp.mean.rarity_Jaccard.80[!is.na(sp.mean.rarity_Jaccard.80[])]
sp.mean.rarity_Jaccard.80 <- temp

# Save species ranges and species continous rarity indices
save(sp.range_size_Jaccard.80, file = "./outputs/Indices_maps/sp.range_size_Jaccard.80.RData", version = "2")
saveRDS(sp.range_size_Jaccard.80, file = "./outputs/Indices_maps/sp.range_size_Jaccard.80.rds", version = "2")
save(sp.rarity_indices_Jaccard.80, file = "./outputs/Indices_maps/sp.rarity_indices_Jaccard.80.RData", version = "2")
saveRDS(sp.rarity_indices_Jaccard.80, file = "./outputs/Indices_maps/sp.rarity_indices_Jaccard.80.rds", version = "2")

# sp.rarity_indices_Jaccard.80 <- readRDS(file = "./outputs/Indices_maps/sp.rarity_indices_Jaccard.80.rds")
# hist(sp.rarity_indices_Jaccard.80)

# Save full and mean rarity maps
save(sp.rarity_Jaccard.80, file = "./outputs/Indices_maps/sp.rarity_Jaccard.80.RData", version = "2")
saveRDS(sp.rarity_Jaccard.80, file = "./outputs/Indices_maps/sp.rarity_Jaccard.80.rds", version = "2")
save(sp.mean.rarity_Jaccard.80, file = "./outputs/Indices_maps/sp.mean.rarity_Jaccard.80.RData", version = "2")
saveRDS(sp.mean.rarity_Jaccard.80, file = "./outputs/Indices_maps/sp.mean.rarity_Jaccard.80.rds", version = "2")

# Plot distribution of ranges and rarity indices
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/sp_range_distri_Jaccard.80.pdf"), height = 5.3, width = 6.5)
hist(sp.range_size_Jaccard.80, main = "Distribution of species geographic ranges\nJaccard.80", breaks = 20, xlab = "Geographic ranges in [10^3] km²")
range_threshold_Jaccard.80 <- round(quantile(sp.range_size_Jaccard.80, p=0.25), 1) # Range threshold value for qualitative rarity
abline(v = range_threshold_Jaccard.80, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",range_threshold_Jaccard.80*1000," km²"), col = "red", lty = 2, lwd = 2, x = "topright", cex = 1, bty ="n")
dev.off()

pdf(file = paste0("./maps/Indices_maps/Jaccard.80/sp_rarity_indices_distri_Jaccard.80.pdf"), height = 5.3, width = 6.5)
hist(sp.rarity_indices_Jaccard.80, main = "Distribution of species rarity indices\nJaccard.80", breaks = 20, xlab = "Rarity indices")
indice_threshold_Jaccard.80 <- round(quantile(sp.rarity_indices_Jaccard.80, p=0.75), 3) # Threshold value for qualitative rarity
abline(v = indice_threshold_Jaccard.80, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",indice_threshold_Jaccard.80), col = "red", lty = 2, lwd = 2, x = "topleft", cex = 1, bty ="n")
dev.off()


## Test other weighting scheme

# sp.range_size_Jaccard.80 <- readRDS(file = "./outputs/Indices_maps/sp.range_size_Jaccard.80.rds")
# hist(sp.range_size_Jaccard.80)
# 
# inv_weights <- max(sp.range_size_Jaccard.80)/sp.range_size_Jaccard.80
# inv_weights <- (weights-min(weights))/(max(weights)-min(weights)) # Normalized
# 
# hist(inv_weights)
# 
# plot(inv_weights ~ sp.range_size_Jaccard.80)

# Load stuff to compute new rarity index
sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
All_sp_proba_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))
load(file = "./outputs/Indices_maps/sp.range_size_Jaccard.80.RData")


library(Rarity)
?rWeights

# # Rarity weights based on inverse normalized range
# inv_weights <- rWeights(occData = sp.range_size_Jaccard.80, wMethods = "invQ", normalised = T, rounding = 5)$invQ
# hist(inv_weights)
# plot(norm_weights ~ sp.range_size_Jaccard.80)

sp_proba_brick_Jaccard.80 <- readAll(brick(All_sp_proba_stack_Jaccard.80))
sp_assemblage <- t(tidyr::drop_na(as.data.frame(sp_proba_brick_Jaccard.80@data@values)))
sp_assemblage <- sp_assemblage[, colSums(sp_assemblage) >= 1]

# # Rarity weights based on an inverse exponential function with inflexion point calibrated as the rarity threshold as such that 25% of species are rare
# Gaston_weights <- rWeights(occData = sp.range_size_Jaccard.80, wMethods = "W", rCutoff = "Gaston",
#                           normalised = T, rounding = 5, assemblages = sp_assemblage)$W
# hist(Gaston_weights)
# plot(Gaston_weights ~ sp.range_size_Jaccard.80)

# Rarity weights based on an inverse exponential function with inflexion point calibrated as the rarity threshold as such that communities host 25% of rare species in average
Leroy_weights_df <- rWeights(occData = sp.range_size_Jaccard.80, wMethods = "W", rCutoff = "Leroy",
                             normalised = T, rounding = 5, assemblages = sp_assemblage)
Leroy_weights <- Leroy_weights_df$W

hist(Leroy_weights)
plot(Leroy_weights ~ sp.range_size_Jaccard.80)

# tested_weights <- inv_weights
# tested_weights <- Gaston_weights
tested_weights <- Leroy_weights

# Apply weights
tested_sp.weighted.stack_Jaccard.80 <- All_sp_proba_stack_Jaccard.80 # Generate the stack to fill
for (i in 1:nlayers(tested_sp.weighted.stack_Jaccard.80)){
  # Multiply species probability of presence by species rarity indices
  tested_sp.weighted.stack_Jaccard.80@layers[[i]]@data@values <- All_sp_proba_stack_Jaccard.80@layers[[i]]@data@values * tested_weights[i]
}

# Sum = Richness weighted by rarity based on Range-size
tested_sp.rarity_Jaccard.80 <- readAll(calc(tested_sp.weighted.stack_Jaccard.80, fun = sum)) 
# Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each species)
tested_sp.mean.rarity_Jaccard.80 <- tested_sp.rarity_Jaccard.80/sp.richness_Jaccard.80


tested_sp.mean.rarity_Jaccard.80 <- tested_sp.mean.rarity_Jaccard.80 + 0.005
# Add 0 values for empty continental pixels (transformed into min values as 0.10)
temp <- continent_mask
temp[!is.na(tested_sp.mean.rarity_Jaccard.80[])] <- tested_sp.mean.rarity_Jaccard.80[!is.na(tested_sp.mean.rarity_Jaccard.80[])]
tested_sp.mean.rarity_Jaccard.80 <- temp

tested_sp.mean.rarity_Jaccard.80_contrasted <- contrasting_raster(x = tested_sp.mean.rarity_Jaccard.80, zmin = 0, zmax = 0.5)

plot(tested_sp.mean.rarity_Jaccard.80, col = pal_bl_red_Mannion)
plot(tested_sp.mean.rarity_Jaccard.80_contrasted, col = pal_bl_red_Mannion)
plot(sp.mean.rarity_Jaccard.80, col = pal_bl_red_Mannion)

pdf(file = paste0("./maps/Indices_maps/Jaccard.80/Comparison_sp_rarity_weights.pdf"), height = 5.3, width = 5.3)
plot(tested_sp.mean.rarity_Jaccard.80[] ~ sp.mean.rarity_Jaccard.80[], ylab = "Leroy's weighted rarity", xlab = "Linearly weighted rarity", main = "Comparison of weighting schemes for rarity")
dev.off()
cor.test(x = tested_sp.mean.rarity_Jaccard.80[], y = sp.mean.rarity_Jaccard.80[], method = "spearman", alternative = "two.sided", na.rm = T)

plot(tested_sp.mean.rarity_Jaccard.80[] ~ sp.richness_Jaccard.80[])
cor.test(x = tested_sp.mean.rarity_Jaccard.80[], y = sp.richness_Jaccard.80[], method = "spearman", alternative = "two.sided", na.rm = T)

threshold_75 <- quantile(x = tested_sp.mean.rarity_Jaccard.80[], p = 0.75, na.rm = T)
threshold_95 <- quantile(x = tested_sp.mean.rarity_Jaccard.80[], p = 0.95, na.rm = T)

plot(tested_sp.mean.rarity_Jaccard.80 > threshold_75)
plot(tested_sp.mean.rarity_Jaccard.80 > threshold_95)

# Plot all stuff to compare the two weighting scheme
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/Comparison_sp_rarity_indices_all.pdf"), height = 7.5, width = 25)
par(mfrow = c(2,7))

plot(tested_weights ~ sp.range_size_Jaccard.80, xlab = "Sp. range", ylab = "Sp. rarity", main = "Leroy's weighting scheme")
hist(tested_weights, xlab = "Sp. rarity weights", main = "Histo of weights")
hist(tested_sp.mean.rarity_Jaccard.80, xlab = 'Community mean geographic rarity', main = "Histo of community scores")
plot(tested_sp.mean.rarity_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = "Sp. geographic rarity")
plot(tested_sp.mean.rarity_Jaccard.80 > quantile(x = tested_sp.mean.rarity_Jaccard.80[], p = 0.75, na.rm = T), main = "Top 25%")
plot(tested_sp.mean.rarity_Jaccard.80 > quantile(x = tested_sp.mean.rarity_Jaccard.80[], p = 0.95, na.rm = T), main = "Top 5%")
plot(tested_sp.mean.rarity_Jaccard.80[] ~ sp.richness_Jaccard.80[], xlab = "Sp. richness", ylab = "Sp. geographic rarity", main = "Correlation with sp. richness")

plot(sp.rarity_indices_Jaccard.80 ~ sp.range_size_Jaccard.80, xlab = "Sp. range", ylab = "Sp. rarity", main = "Linear weighting scheme")
hist(sp.rarity_indices_Jaccard.80, xlab = "Sp. rarity weights", main = "Histo of weights")
hist(sp.mean.rarity_Jaccard.80, xlab = 'Community mean geographic rarity', main = "Histo of community scores")
plot(sp.mean.rarity_Jaccard.80, col = pal_bl_red_Mannion, main = "Sp. geographic rarity")
plot(sp.mean.rarity_Jaccard.80 > quantile(x = sp.mean.rarity_Jaccard.80[], p = 0.75, na.rm = T), main = "Top 25%")
plot(sp.mean.rarity_Jaccard.80 > quantile(x = sp.mean.rarity_Jaccard.80[], p = 0.95, na.rm = T), main = "Top 5%")
plot(sp.mean.rarity_Jaccard.80[] ~ sp.richness_Jaccard.80[], xlab = "Sp. richness", ylab = "Sp. geographic rarity", main = "Correlation with sp. richness")

par(mfrow = c(1,1))
dev.off()

# Plot the distribution of weights (species rarity indices) for Leroy's weighting
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/sp_rarity_indices_distri_Jaccard.80_Leroy_weights.pdf"), height = 5.3, width = 6.5)
hist(Leroy_weights, main = "Distribution of species rarity indices\nLeroy's weights", breaks = 20, xlab = "Rarity indices")
indice_threshold <- round(quantile(Leroy_weights, p = 0.75), 3) # Threshold value for Gaston's qualitative rarity
abline(v = 1-Leroy_weights_df$cut.off[1], col = "red", lty = 2, lwd = 2)
abline(v = indice_threshold, col = "red", lty = 3, lwd = 2)
legend(legend = c("Threshold for rarity", paste0("Leroy's scheme: ", round(1-Leroy_weights_df$cut.off[1], 3)), paste0("Gaston's 25%: ", indice_threshold)), col = "red", lty = c(NA, 2, 3), lwd = 2, x = "top", cex = 1, bty ="n")
dev.off()

# Save the stuff for the new rarity index based on Leroy (2013)'s weighting
sp.rarity_indices_Leroy_Jaccard.80 <- Leroy_weights
sp.rarity_Leroy_Jaccard.80 <- tested_sp.rarity_Jaccard.80
sp.mean.rarity_Leroy_Jaccard.80 <- tested_sp.mean.rarity_Jaccard.80
sp.mean.rarity_Leroy_Jaccard.80_contrasted <- tested_sp.mean.rarity_Jaccard.80_contrasted

# Save species ranges and species continous rarity indices
save(sp.rarity_indices_Leroy_Jaccard.80, file = "./outputs/Indices_maps/sp.rarity_indices_Leroy_Jaccard.80.RData", version = "2")
saveRDS(sp.rarity_indices_Leroy_Jaccard.80, file = "./outputs/Indices_maps/sp.rarity_indices_Leroy_Jaccard.80.rds", version = "2")

# Save full and mean rarity maps
save(sp.rarity_Leroy_Jaccard.80, file = "./outputs/Indices_maps/sp.rarity_Leroy_Jaccard.80.RData", version = "2")
saveRDS(sp.rarity_Leroy_Jaccard.80, file = "./outputs/Indices_maps/sp.rarity_Leroy_Jaccard.80.rds", version = "2")
save(sp.mean.rarity_Leroy_Jaccard.80, file = "./outputs/Indices_maps/sp.mean.rarity_Leroy_Jaccard.80.RData", version = "2")
saveRDS(sp.mean.rarity_Leroy_Jaccard.80, file = "./outputs/Indices_maps/sp.mean.rarity_Leroy_Jaccard.80.rds", version = "2")
save(sp.mean.rarity_Leroy_Jaccard.80_contrasted, file = "./outputs/Indices_maps/sp.mean.rarity_Leroy_Jaccard.80_contrasted.RData", version = "2")
saveRDS(sp.mean.rarity_Leroy_Jaccard.80_contrasted, file = "./outputs/Indices_maps/sp.mean.rarity_Leroy_Jaccard.80_contrasted.rds", version = "2")


## For Jaccard.95

# Compute weights based on geographic ranges
sp.range_size_Jaccard.95 <- NA
for (i in 1:nlayers(All_sp_proba_stack_Jaccard.95)) {
  # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
  # Compute the estimated number of pixels occupied by the species as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
  sp.range_size_Jaccard.95[i] <- sum(All_sp_proba_stack_Jaccard.95[[i]]@data@values, na.rm = T)*774.51/1000
}
sp.rarity_indices_Jaccard.95 <- 1-(sp.range_size_Jaccard.95/max(sp.range_size_Jaccard.95)) # Weights by max species extent = rarity indices

# Apply weights
sp.weighted.stack_Jaccard.95 <- All_sp_proba_stack_Jaccard.95 # Generate the stack to fill
for (i in 1:nlayers(sp.weighted.stack_Jaccard.95)){
  # Multiply species probability of presence by species rarity indices
  sp.weighted.stack_Jaccard.95@layers[[i]]@data@values <- All_sp_proba_stack_Jaccard.95@layers[[i]]@data@values * sp.rarity_indices_Jaccard.95[i]
}

# plot(sp.weighted.stack_Jaccard.95)

# Sum = Richness weighted by rarity based on Range-size
sp.rarity_Jaccard.95 <- readAll(calc(sp.weighted.stack_Jaccard.95, fun = sum)) 
# Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each species)
sp.mean.rarity_Jaccard.95 <- sp.rarity_Jaccard.95/sp.richness_Jaccard.95

# Add 0 values for empty continental pixels (transformed into min values as 0.10)
temp <- continent_mask
temp[!is.na(sp.mean.rarity_Jaccard.95[])] <- sp.mean.rarity_Jaccard.95[!is.na(sp.mean.rarity_Jaccard.95[])]
sp.mean.rarity_Jaccard.95 <- temp

# Save species ranges and species continuous rarity indices
save(sp.range_size_Jaccard.95, file = "./outputs/Indices_maps/sp.range_size_Jaccard.95.RData", version = "2")
saveRDS(sp.range_size_Jaccard.95, file = "./outputs/Indices_maps/sp.range_size_Jaccard.95.rds", version = "2")
save(sp.rarity_indices_Jaccard.95, file = "./outputs/Indices_maps/sp.rarity_indices_Jaccard.95.RData", version = "2")
saveRDS(sp.rarity_indices_Jaccard.95, file = "./outputs/Indices_maps/sp.rarity_indices_Jaccard.95.rds", version = "2")

# Save full and mean rarity maps
save(sp.rarity_Jaccard.95, file = "./outputs/Indices_maps/sp.rarity_Jaccard.95.RData", version = "2")
saveRDS(sp.rarity_Jaccard.95, file = "./outputs/Indices_maps/sp.rarity_Jaccard.95.rds", version = "2")
save(sp.mean.rarity_Jaccard.95, file = "./outputs/Indices_maps/sp.mean.rarity_Jaccard.95.RData", version = "2")
saveRDS(sp.mean.rarity_Jaccard.95, file = "./outputs/Indices_maps/sp.mean.rarity_Jaccard.95.rds", version = "2")

# Plot distribution of ranges and rarity indices
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/sp_range_distri_Jaccard.95.pdf"), height = 5.3, width = 6.5)
hist(sp.range_size_Jaccard.95, main = "Distribution of species geographic ranges\nJaccard.95", breaks = 20, xlab = "Geographic ranges in [10^3] km²")
range_threshold_Jaccard.95 <- round(quantile(sp.range_size_Jaccard.95, p=0.25), 1) # Range threshold value for qualitative rarity
abline(v = range_threshold_Jaccard.95, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",range_threshold_Jaccard.95*1000," km²"), col = "red", lty = 2, lwd = 2, x = "topright", cex = 1, bty ="n")
dev.off()

pdf(file = paste0("./maps/Indices_maps/Jaccard.95/sp_rarity_indices_distri_Jaccard.95.pdf"), height = 5.3, width = 6.5)
hist(sp.rarity_indices_Jaccard.95, main = "Distribution of species rarity indices\nJaccard.95", breaks = 20, xlab = "Rarity indices")
indice_threshold_Jaccard.95 <- round(quantile(sp.rarity_indices_Jaccard.95, p=0.75), 3) # Threshold value for qualitative rarity
abline(v = indice_threshold_Jaccard.95, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",indice_threshold_Jaccard.95), col = "red", lty = 2, lwd = 2, x = "topleft", cex = 1, bty ="n")
dev.off()


## For TSS.80

# Compute weights based on geographic ranges
sp.range_size_TSS.80 <- NA
for (i in 1:nlayers(All_sp_proba_stack_TSS.80)) {
  # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
  # Compute the estimated number of pixels occupied by the species as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
  sp.range_size_TSS.80[i] <- sum(All_sp_proba_stack_TSS.80[[i]]@data@values, na.rm = T)*774.51/1000
}
sp.rarity_indices_TSS.80 <- 1-(sp.range_size_TSS.80/max(sp.range_size_TSS.80)) # Weights by max species extent = rarity indices

# Apply weights
sp.weighted.stack_TSS.80 <- All_sp_proba_stack_TSS.80 # Generate the stack to fill
for (i in 1:nlayers(sp.weighted.stack_TSS.80)){
  # Multiply species probability of presence by species rarity indices
  sp.weighted.stack_TSS.80@layers[[i]]@data@values <- All_sp_proba_stack_TSS.80@layers[[i]]@data@values * sp.rarity_indices_TSS.80[i]
}

plot(sp.weighted.stack_TSS.80)

# Sum = Richness weighted by rarity based on Range-size
sp.rarity_TSS.80 <- readAll(calc(sp.weighted.stack_TSS.80, fun = sum)) 
# Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each species)
sp.mean.rarity_TSS.80 <- sp.rarity_TSS.80/sp.richness_TSS.80

# Add 0 values for empty continental pixels (transformed into min values as 0.10)
temp <- continent_mask
temp[!is.na(sp.mean.rarity_TSS.80[])] <- sp.mean.rarity_TSS.80[!is.na(sp.mean.rarity_TSS.80[])]
sp.mean.rarity_TSS.80 <- temp

# Save species ranges and species continous rarity indices
save(sp.range_size_TSS.80, file = "./outputs/Indices_maps/sp.range_size_TSS.80.RData", version = "2")
saveRDS(sp.range_size_TSS.80, file = "./outputs/Indices_maps/sp.range_size_TSS.80.rds", version = "2")
save(sp.rarity_indices_TSS.80, file = "./outputs/Indices_maps/sp.rarity_indices_TSS.80.RData", version = "2")
saveRDS(sp.rarity_indices_TSS.80, file = "./outputs/Indices_maps/sp.rarity_indices_TSS.80.rds", version = "2")

# Save full and mean rarity maps
save(sp.rarity_TSS.80, file = "./outputs/Indices_maps/sp.rarity_TSS.80.RData", version = "2")
saveRDS(sp.rarity_TSS.80, file = "./outputs/Indices_maps/sp.rarity_TSS.80.rds", version = "2")
save(sp.mean.rarity_TSS.80, file = "./outputs/Indices_maps/sp.mean.rarity_TSS.80.RData", version = "2")
saveRDS(sp.mean.rarity_TSS.80, file = "./outputs/Indices_maps/sp.mean.rarity_TSS.80.rds", version = "2")

# Plot distribution of ranges and rarity indices
pdf(file = paste0("./maps/Indices_maps/TSS.80/sp_range_distri_TSS.80.pdf"), height = 5.3, width = 6.5)
hist(sp.range_size_TSS.80, main = "Distribution of species geographic ranges\nTSS.80", breaks = 20, xlab = "Geographic ranges in [10^3] km²")
range_threshold_TSS.80 <- round(quantile(sp.range_size_TSS.80, p=0.25), 1) # Range threshold value for qualitative rarity
abline(v = range_threshold_TSS.80, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",range_threshold_TSS.80*1000," km²"), col = "red", lty = 2, lwd = 2, x = "topright", cex = 1, bty ="n")
dev.off()

pdf(file = paste0("./maps/Indices_maps/TSS.80/sp_rarity_indices_distri_TSS.80.pdf"), height = 5.3, width = 6.5)
hist(sp.rarity_indices_TSS.80, main = "Distribution of species rarity indices\nTSS.80", breaks = 20, xlab = "Rarity indices")
indice_threshold_TSS.80 <- round(quantile(sp.rarity_indices_TSS.80, p=0.75), 3) # Threshold value for qualitative rarity
abline(v = indice_threshold_TSS.80, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",indice_threshold_TSS.80), col = "red", lty = 2, lwd = 2, x = "topleft", cex = 1, bty ="n")
dev.off()


## For TSS.95

# Compute weights based on geographic ranges
sp.range_size_TSS.95 <- NA
for (i in 1:nlayers(All_sp_proba_stack_TSS.95)) {
  # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
  # Compute the estimated number of pixels occupied by the species as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
  sp.range_size_TSS.95[i] <- sum(All_sp_proba_stack_TSS.95[[i]]@data@values, na.rm = T)*774.51/1000
}
sp.rarity_indices_TSS.95 <- 1-(sp.range_size_TSS.95/max(sp.range_size_TSS.95)) # Weights by max species extent = rarity indices

# Apply weights
sp.weighted.stack_TSS.95 <- All_sp_proba_stack_TSS.95 # Generate the stack to fill
for (i in 1:nlayers(sp.weighted.stack_TSS.95)){
  # Multiply species probability of presence by species rarity indices
  sp.weighted.stack_TSS.95@layers[[i]]@data@values <- All_sp_proba_stack_TSS.95@layers[[i]]@data@values * sp.rarity_indices_TSS.95[i]
}

plot(sp.weighted.stack_TSS.95)

# Sum = Richness weighted by rarity based on Range-size
sp.rarity_TSS.95 <- readAll(calc(sp.weighted.stack_TSS.95, fun = sum)) 
# Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each species)
sp.mean.rarity_TSS.95 <- sp.rarity_TSS.95/sp.richness_TSS.95

# Add 0 values for empty continental pixels (transformed into min values as 0.10)
temp <- continent_mask
temp[!is.na(sp.mean.rarity_TSS.95[])] <- sp.mean.rarity_TSS.95[!is.na(sp.mean.rarity_TSS.95[])]
sp.mean.rarity_TSS.95 <- temp

# Save species ranges and species continuous rarity indices
save(sp.range_size_TSS.95, file = "./outputs/Indices_maps/sp.range_size_TSS.95.RData", version = "2")
saveRDS(sp.range_size_TSS.95, file = "./outputs/Indices_maps/sp.range_size_TSS.95.rds", version = "2")
save(sp.rarity_indices_TSS.95, file = "./outputs/Indices_maps/sp.rarity_indices_TSS.95.RData", version = "2")
saveRDS(sp.rarity_indices_TSS.95, file = "./outputs/Indices_maps/sp.rarity_indices_TSS.95.rds", version = "2")

# Save full and mean rarity maps
save(sp.rarity_TSS.95, file = "./outputs/Indices_maps/sp.rarity_TSS.95.RData", version = "2")
saveRDS(sp.rarity_TSS.95, file = "./outputs/Indices_maps/sp.rarity_TSS.95.rds", version = "2")
save(sp.mean.rarity_TSS.95, file = "./outputs/Indices_maps/sp.mean.rarity_TSS.95.RData", version = "2")
saveRDS(sp.mean.rarity_TSS.95, file = "./outputs/Indices_maps/sp.mean.rarity_TSS.95.rds", version = "2")

# Plot distribution of ranges and rarity indices
pdf(file = paste0("./maps/Indices_maps/TSS.95/sp_range_distri_TSS.95.pdf"), height = 5.3, width = 6.5)
hist(sp.range_size_TSS.95, main = "Distribution of species geographic ranges\nTSS.95", breaks = 20, xlab = "Geographic ranges in [10^3] km²")
range_threshold_TSS.95 <- round(quantile(sp.range_size_TSS.95, p=0.25), 1) # Range threshold value for qualitative rarity
abline(v = range_threshold_TSS.95, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",range_threshold_TSS.95*1000," km²"), col = "red", lty = 2, lwd = 2, x = "topright", cex = 1, bty ="n")
dev.off()

pdf(file = paste0("./maps/Indices_maps/TSS.95/sp_rarity_indices_distri_TSS.95.pdf"), height = 5.3, width = 6.5)
hist(sp.rarity_indices_TSS.95, main = "Distribution of species rarity indices\nTSS.95", breaks = 20, xlab = "Rarity indices")
indice_threshold_TSS.95 <- round(quantile(sp.rarity_indices_TSS.95, p=0.75), 3) # Threshold value for qualitative rarity
abline(v = indice_threshold_TSS.95, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",indice_threshold_TSS.95), col = "red", lty = 2, lwd = 2, x = "topleft", cex = 1, bty ="n")
dev.off()

### 9.2/ Plot community rarity maps ####

### Load directly the community rarity maps
sp.rarity_Jaccard.80 <- readRDS(file = "./outputs/Indices_maps/sp.rarity_Jaccard.80.rds")
sp.rarity_Jaccard.95 <- readRDS(file = "./outputs/Indices_maps/sp.rarity_Jaccard.95.rds")
sp.rarity_TSS.80 <- readRDS(file = "./outputs/Indices_maps/sp.rarity_TSS.80.rds")
sp.rarity_TSS.95 <- readRDS(file = "./outputs/Indices_maps/sp.rarity_TSS.95.rds")


### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/sp.rarity_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.rarity_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Species rarity \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_Jaccard.80, locs = seq(0, 100, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/sp.rarity_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.rarity_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Species rarity \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_Jaccard.95, locs = seq(0, 100, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")
par(mar = internal_margins)
dev.off()


# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/sp.rarity_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.rarity_TSS.80, col = pal_bl_red_Mannion, main = paste0("Species rarity \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_TSS.80, locs = seq(0, 100, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")
par(mar = internal_margins)
dev.off()


# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/sp.rarity_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.rarity_TSS.95, col = pal_bl_red_Mannion, main = paste0("Species rarity \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_TSS.95, locs = seq(0, 100, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(sp.rarity_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/sp.rarity_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.rarity_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Species rarity \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_Jaccard.80, locs = seq(0, 100, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")

image(sp.rarity_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Species rarity \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_Jaccard.95, locs = seq(0, 100, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")

image(sp.rarity_TSS.80, col = pal_bl_red_Mannion, main = paste0("Species rarity \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_TSS.80, locs = seq(0, 100, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")

image(sp.rarity_TSS.95, col = pal_bl_red_Mannion, main = paste0("Species rarity \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_TSS.95, locs = seq(0, 100, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")

par(mar = internal_margins)

dev.off()


### 9.3/ Plot mean species rarity maps ####

### Load directly the mean species rarity maps
sp.mean.rarity_Jaccard.80 <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_Jaccard.80.rds")
sp.mean.rarity_Jaccard.95 <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_Jaccard.95.rds")
sp.mean.rarity_TSS.80 <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_TSS.80.rds")
sp.mean.rarity_TSS.95 <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_TSS.95.rds")


### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/sp.mean.rarity_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.mean.rarity_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mean species rarity\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_Jaccard.80, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/sp.mean.rarity_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.mean.rarity_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Mean species rarity\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_Jaccard.95, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")
par(mar = internal_margins)
dev.off()


# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/sp.mean.rarity_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.mean.rarity_TSS.80, col = pal_bl_red_Mannion, main = paste0("Mean species rarity\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_TSS.80, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/sp.mean.rarity_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.mean.rarity_TSS.95, col = pal_bl_red_Mannion, main = paste0("Mean species rarity \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_TSS.95, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(sp.mean.rarity_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/sp.mean.rarity_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.mean.rarity_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mean species rarity \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_Jaccard.80, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")

image(sp.mean.rarity_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Mean species rarity\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_Jaccard.95, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")

image(sp.mean.rarity_TSS.80, col = pal_bl_red_Mannion, main = paste0("Mean species rarity\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_TSS.80, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")

image(sp.mean.rarity_TSS.95, col = pal_bl_red_Mannion, main = paste0("Mean species rarity\nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_TSS.95, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")

par(mar = internal_margins)

dev.off()


# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])


###### 10/ Categorical Species Rarity (25% threshold) ######

### Load directly the stack of sp probas use to compute rarity based on geographic range
All_sp_proba_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))
All_sp_proba_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.95.rds"))
All_sp_proba_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80.rds"))
All_sp_proba_stack_TSS.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.95.rds"))

### Load species ranges
sp.range_size_Jaccard.80 <- readRDS(file = "./outputs/Indices_maps/sp.range_size_Jaccard.80.rds")
sp.range_size_Jaccard.95 <- readRDS(file = "./outputs/Indices_maps/sp.range_size_Jaccard.95.rds")
sp.range_size_TSS.80 <- readRDS(file = "./outputs/Indices_maps/sp.range_size_TSS.80.rds")
sp.range_size_TSS.95 <- readRDS(file = "./outputs/Indices_maps/sp.range_size_TSS.95.rds")

### 10.1/ Extract rarity threshold ####
rarity_threshold_Jaccard.80 <- round(quantile(sp.range_size_Jaccard.80, probs = 0.25), 3) ; rarity_threshold_Jaccard.80
# Threshold surface to be considered a rare species = 95 086 km² 
rarity_threshold_Jaccard.95 <- quantile(sp.range_size_Jaccard.95, probs = 0.25) ; rarity_threshold_Jaccard.95
# Threshold surface to be considered a rare species = 105 990 km² 
rarity_threshold_TSS.80 <- quantile(sp.range_size_TSS.80, probs = 0.25) ; rarity_threshold_TSS.80
# Threshold surface to be considered a rare species = 95 086 km² 
rarity_threshold_TSS.95 <- quantile(sp.range_size_TSS.95, probs = 0.25) ; rarity_threshold_TSS.95
# Threshold surface to be considered a rare species = 105 990 km² 

### 10.2/ Plot species range distribution ####

# Plot species range distribution for Jaccard.80
ordered_sp.range_size_Jaccard.80 <- sp.range_size_Jaccard.80[order(sp.range_size_Jaccard.80)]

pdf(file = paste0("./maps/Indices_maps/Jaccard.80/sp_range_distri_barplot_Jaccard.80.pdf"), height = 5.3, width = 6.5)
barplot(ordered_sp.range_size_Jaccard.80, ylim = c(0,9000), ylab = "Surface in 10^3 km", xlab = "Species", main = "Distribution of Species range size")
abline(v = length(ordered_sp.range_size_Jaccard.80)/4, col = "red", lwd = 2, lty = 2)
legend(legend = c("Threshold 25 %", paste0("    ",round(rarity_threshold_Jaccard.80*1000, 0)," km²")), bty = "n", x = "top", text.col = "red", text.font = 2)
legend(legend = c("  RARE", "Species"), bty = "n", x = "left", text.col = "darkblue", text.font = 2)
dev.off()

# Plot species range distribution for Jaccard.95
ordered_sp.range_size_Jaccard.95 <- sp.range_size_Jaccard.95[order(sp.range_size_Jaccard.95)]

pdf(file = paste0("./maps/Indices_maps/Jaccard.95/sp_range_distri_barplot_Jaccard.95.pdf"), height = 5.3, width = 6.5)
barplot(ordered_sp.range_size_Jaccard.95, ylim = c(0,9000), ylab = "Surface in 10^3 km", xlab = "Species", main = "Distribution of Species range size")
abline(v = length(ordered_sp.range_size_Jaccard.95)/4, col = "red", lwd = 2, lty = 2)
legend(legend = c("Threshold 25 %", paste0("    ",round(rarity_threshold_Jaccard.95*1000, 0)," km²")), bty = "n", x = "top", text.col = "red", text.font = 2)
legend(legend = c("  RARE", "Species"), bty = "n", x = "left", text.col = "darkblue", text.font = 2)
dev.off()

# Plot species range distribution for TSS.80
ordered_sp.range_size_TSS.80 <- sp.range_size_TSS.80[order(sp.range_size_TSS.80)]

pdf(file = paste0("./maps/Indices_maps/TSS.80/sp_range_distri_barplot_TSS.80.pdf"), height = 5.3, width = 6.5)
barplot(ordered_sp.range_size_TSS.80, ylim = c(0,9000), ylab = "Surface in 10^3 km", xlab = "Species", main = "Distribution of Species range size")
abline(v = length(ordered_sp.range_size_TSS.80)/4, col = "red", lwd = 2, lty = 2)
legend(legend = c("Threshold 25 %", paste0("    ",round(rarity_threshold_TSS.80*1000, 0)," km²")), bty = "n", x = "top", text.col = "red", text.font = 2)
legend(legend = c("  RARE", "Species"), bty = "n", x = "left", text.col = "darkblue", text.font = 2)
dev.off()

# Plot species range distribution for TSS.95
ordered_sp.range_size_TSS.95 <- sp.range_size_TSS.95[order(sp.range_size_TSS.95)]

pdf(file = paste0("./maps/Indices_maps/TSS.95/sp_range_distri_barplot_TSS.95.pdf"), height = 5.3, width = 6.5)
barplot(ordered_sp.range_size_TSS.95, ylim = c(0,9000), ylab = "Surface in 10^3 km", xlab = "Species", main = "Distribution of Species range size")
abline(v = length(ordered_sp.range_size_TSS.95)/4, col = "red", lwd = 2, lty = 2)
legend(legend = c("Threshold 25 %", paste0("    ",round(rarity_threshold_TSS.95*1000, 0)," km²")), bty = "n", x = "top", text.col = "red", text.font = 2)
legend(legend = c("  RARE", "Species"), bty = "n", x = "left", text.col = "darkblue", text.font = 2)
dev.off()

### 10.3/ Extract only rare species to compute categorical rarity indices ####

# Detect rare species
rare.sp_Jaccard.80 <- sp.range_size_Jaccard.80 < rarity_threshold_Jaccard.80
rare.sp_Jaccard.95 <- sp.range_size_Jaccard.95 < rarity_threshold_Jaccard.95
rare.sp_TSS.80 <- sp.range_size_TSS.80 < rarity_threshold_TSS.80
rare.sp_TSS.95 <- sp.range_size_TSS.95 < rarity_threshold_TSS.95

sum(rare.sp_Jaccard.80) # 97 rare species on 388 (25%)

# Extract only rare species
Rare_sp_proba_stack_Jaccard.80 <- All_sp_proba_stack_Jaccard.80[[which(rare.sp_Jaccard.80)]]
Rare_sp_proba_stack_Jaccard.95 <- All_sp_proba_stack_Jaccard.95[[which(rare.sp_Jaccard.95)]]
Rare_sp_proba_stack_TSS.80 <- All_sp_proba_stack_TSS.80[[which(rare.sp_TSS.80)]]
Rare_sp_proba_stack_TSS.95 <- All_sp_proba_stack_TSS.95[[which(rare.sp_TSS.95)]]

### 10.4/ Index computation ####

# Load maps of Species richness to use to compute mean rarity index among local species
sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
sp.richness_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.95.rds"))
sp.richness_TSS.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_TSS.80.rds"))
sp.richness_TSS.95 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_TSS.95.rds"))

# Index computation for total of rare species and proportion of rare species
tot.rare.sp_Jaccard.80 <- readAll(calc(Rare_sp_proba_stack_Jaccard.80, fun = sum))
prop.rare.sp_Jaccard.80 <- tot.rare.sp_Jaccard.80/sp.richness_Jaccard.80*100
tot.rare.sp_Jaccard.95 <- readAll(calc(Rare_sp_proba_stack_Jaccard.95, fun = sum))
prop.rare.sp_Jaccard.95 <- tot.rare.sp_Jaccard.95/sp.richness_Jaccard.95*100
tot.rare.sp_TSS.80 <- readAll(calc(Rare_sp_proba_stack_TSS.80, fun = sum))
prop.rare.sp_TSS.80 <- tot.rare.sp_TSS.80/sp.richness_TSS.80*100
tot.rare.sp_TSS.95 <- readAll(calc(Rare_sp_proba_stack_TSS.95, fun = sum))
prop.rare.sp_TSS.95 <- tot.rare.sp_TSS.95/sp.richness_TSS.95*100

# Save maps
save(tot.rare.sp_Jaccard.80 , file = "./outputs/Indices_maps/tot.rare.sp_Jaccard.80.RData", version = "2")
saveRDS(tot.rare.sp_Jaccard.80 , file = "./outputs/Indices_maps/tot.rare.sp_Jaccard.80.rds", version = "2")
save(tot.rare.sp_Jaccard.95 , file = "./outputs/Indices_maps/tot.rare.sp_Jaccard.95.RData", version = "2")
saveRDS(tot.rare.sp_Jaccard.95 , file = "./outputs/Indices_maps/tot.rare.sp_Jaccard.95.rds", version = "2")
save(tot.rare.sp_TSS.80 , file = "./outputs/Indices_maps/tot.rare.sp_TSS.80.RData", version = "2")
saveRDS(tot.rare.sp_TSS.80 , file = "./outputs/Indices_maps/tot.rare.sp_TSS.80.rds", version = "2")
save(tot.rare.sp_TSS.95 , file = "./outputs/Indices_maps/tot.rare.sp_TSS.95.RData", version = "2")
saveRDS(tot.rare.sp_TSS.95 , file = "./outputs/Indices_maps/tot.rare.sp_TSS.95.rds", version = "2")

save(prop.rare.sp_Jaccard.80 , file = "./outputs/Indices_maps/prop.rare.sp_Jaccard.80.RData", version = "2")
saveRDS(prop.rare.sp_Jaccard.80 , file = "./outputs/Indices_maps/prop.rare.sp_Jaccard.80.rds", version = "2")
save(prop.rare.sp_Jaccard.95 , file = "./outputs/Indices_maps/prop.rare.sp_Jaccard.95.RData", version = "2")
saveRDS(prop.rare.sp_Jaccard.95 , file = "./outputs/Indices_maps/prop.rare.sp_Jaccard.95.rds", version = "2")
save(prop.rare.sp_TSS.80 , file = "./outputs/Indices_maps/prop.rare.sp_TSS.80.RData", version = "2")
saveRDS(prop.rare.sp_TSS.80 , file = "./outputs/Indices_maps/prop.rare.sp_TSS.80.rds", version = "2")
save(prop.rare.sp_TSS.95 , file = "./outputs/Indices_maps/prop.rare.sp_TSS.95.RData", version = "2")
saveRDS(prop.rare.sp_TSS.95 , file = "./outputs/Indices_maps/prop.rare.sp_TSS.95.rds", version = "2")

### 10.5 /Plots of rare species richness ####

### Load directly the final Rare Species richness layer
tot.rare.sp_Jaccard.80 <- readRDS(file = "./outputs/Indices_maps/tot.rare.sp_Jaccard.80.rds")
tot.rare.sp_Jaccard.95 <- readRDS(file = "./outputs/Indices_maps/tot.rare.sp_Jaccard.95.rds")
tot.rare.sp_TSS.80 <- readRDS(file = "./outputs/Indices_maps/tot.rare.sp_TSS.80.rds")
tot.rare.sp_TSS.95 <- readRDS(file = "./outputs/Indices_maps/tot.rare.sp_TSS.95.rds")


### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/tot.rare.sp_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.rare.sp_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Rare species richness \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_Jaccard.80, locs = seq(0, 15, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Rare\nspecies")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/tot.rare.sp_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.rare.sp_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Rare species richness \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_Jaccard.95, locs = seq(0, 15, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Rare\nspecies")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/tot.rare.sp_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.rare.sp_TSS.80, col = pal_bl_red_Mannion, main = paste0("Rare species richness \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_TSS.80, locs = seq(0, 15, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Rare\nspecies")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/tot.rare.sp_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.rare.sp_TSS.95, col = pal_bl_red_Mannion, main = paste0("Rare species richness \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_TSS.95, locs = seq(0, 15, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Rare\nspecies")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(tot.sp.richness_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/tot.rare.sp_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(tot.rare.sp_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Rare species richness \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_Jaccard.80, locs = seq(0, 15, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Rare\nspecies")

image(tot.rare.sp_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Rare species richness \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_Jaccard.95, locs = seq(0, 15, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Rare\nspecies")

image(tot.rare.sp_TSS.80, col = pal_bl_red_Mannion, main = paste0("Rare species richness \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_TSS.80, locs = seq(0, 15, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Rare\nspecies")

image(tot.rare.sp_TSS.95, col = pal_bl_red_Mannion, main = paste0("Rare species richness \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_TSS.95, locs = seq(0, 15, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Rare\nspecies")

par(mar = internal_margins)

dev.off()


### 10.5 /Plots of proportion of rare species ####

### Load directly the final Rare Species richness layer
prop.rare.sp_Jaccard.80 <- readRDS(file = "./outputs/Indices_maps/prop.rare.sp_Jaccard.80.rds")
prop.rare.sp_Jaccard.95 <- readRDS(file = "./outputs/Indices_maps/prop.rare.sp_Jaccard.95.rds")
prop.rare.sp_TSS.80 <- readRDS(file = "./outputs/Indices_maps/prop.rare.sp_TSS.80.rds")
prop.rare.sp_TSS.95 <- readRDS(file = "./outputs/Indices_maps/prop.rare.sp_TSS.95.rds")



# Add 0 values for empty continental pixels (transformed into min values as 0.25)
temp <- continent_mask
temp[!is.na(prop.rare.sp_Jaccard.80[])] <- prop.rare.sp_Jaccard.80[!is.na(prop.rare.sp_Jaccard.80)]
prop.rare.sp_Jaccard.80 <- temp

temp <- continent_mask
temp[!is.na(prop.rare.sp_Jaccard.95[])] <- prop.rare.sp_Jaccard.95[!is.na(prop.rare.sp_Jaccard.95)]
prop.rare.sp_Jaccard.95 <- temp

temp <- continent_mask
temp[!is.na(prop.rare.sp_TSS.80[])] <- prop.rare.sp_TSS.80[!is.na(prop.rare.sp_TSS.80)]
prop.rare.sp_TSS.80 <- temp

temp <- continent_mask
temp[!is.na(prop.rare.sp_TSS.95[])] <- prop.rare.sp_TSS.95[!is.na(prop.rare.sp_TSS.95)]
prop.rare.sp_TSS.95 <- temp

# Contrast by merging all values higher than 25 %
prop.rare.sp_Jaccard.80[prop.rare.sp_Jaccard.80 > 25] <-  25
prop.rare.sp_Jaccard.95[prop.rare.sp_Jaccard.95 > 25] <-  25
prop.rare.sp_TSS.80[prop.rare.sp_TSS.80 > 25] <-  25
prop.rare.sp_TSS.95[prop.rare.sp_TSS.95 > 25] <-  25


### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/prop.rare.sp_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(prop.rare.sp_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Proportion of rare species\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_Jaccard.80, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/prop.rare.sp_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(prop.rare.sp_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Proportion of rare species\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_Jaccard.95, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/prop.rare.sp_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(prop.rare.sp_TSS.80, col = pal_bl_red_Mannion, main = paste0("Proportion of rare species\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_TSS.80, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/prop.rare.sp_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(prop.rare.sp_TSS.95, col = pal_bl_red_Mannion, main = paste0("Proportion of rare species\nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_TSS.95, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(prop.sp.richness_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/prop.rare.sp_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(prop.rare.sp_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Proportion of rare species\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_Jaccard.80, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")

image(prop.rare.sp_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Proportion of rare species\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_Jaccard.95, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")

image(prop.rare.sp_TSS.80, col = pal_bl_red_Mannion, main = paste0("Proportion of rare species\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_TSS.80, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")

image(prop.rare.sp_TSS.95, col = pal_bl_red_Mannion, main = paste0("Proportion of rare species\nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_TSS.95, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")

par(mar = internal_margins)

dev.off()


# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])


###### 11/ Continuous Mimicry rarity = Range size weighted mimicry richness ######

### Load directly the stack of ring probas use to compute rarity based on geographic range
All_ring_proba_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_Jaccard.80.RData"))
All_ring_proba_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_Jaccard.95.RData"))
All_ring_proba_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_TSS.80.RData"))
All_ring_proba_stack_TSS.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_TSS.95.RData"))

# Load maps of Mimicry richness to use to compute mean rarity index among local rings
ring.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))
ring.richness_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.95.rds"))
ring.richness_TSS.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_TSS.80.rds"))
ring.richness_TSS.95 <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_TSS.95.rds"))

### 11.1/ Index computation and plot distribution of ranges and rarity indices ####

## For Jaccard.80

# Compute weights based on geographic ranges
ring.range_size_Jaccard.80 <- NA
for (i in 1:nlayers(All_ring_proba_stack_Jaccard.80)) {
  # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
  # Compute the estimated number of pixels occupied by the ring as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
  ring.range_size_Jaccard.80[i] <- sum(All_ring_proba_stack_Jaccard.80[[i]]@data@values, na.rm = T)*774.51/1000
}
ring.rarity_Jaccard.80 <- 1-(ring.range_size_Jaccard.80/max(ring.range_size_Jaccard.80)) # Weights by max ring extent = rarity indices

# Apply weights
ring.weighted.stack_Jaccard.80 <- All_ring_proba_stack_Jaccard.80 # Generate the stack to fill
for (i in 1:nlayers(ring.weighted.stack_Jaccard.80)){
  # Multiply ring probability of presence by ring rarity indices
  ring.weighted.stack_Jaccard.80@layers[[i]]@data@values <- All_ring_proba_stack_Jaccard.80@layers[[i]]@data@values * ring.rarity_Jaccard.80[i]
}

plot(ring.weighted.stack_Jaccard.80)

# Sum = Richness weighted by rarity based on Range-size
mimicry.rarity_Jaccard.80 <- calc(ring.weighted.stack_Jaccard.80, fun = sum) 
# Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each ring)
mimicry.mean.rarity_Jaccard.80 <- mimicry.rarity_Jaccard.80/ring.richness_Jaccard.80

# Add 0 values for empty continental pixels (transformed into min values as 0.10)
temp <- continent_mask
temp[!is.na(mimicry.mean.rarity_Jaccard.80[])] <- mimicry.mean.rarity_Jaccard.80[!is.na(mimicry.mean.rarity_Jaccard.80[])]
mimicry.mean.rarity_Jaccard.80 <- temp

# Save mimicry ring ranges and mimicry ring continous rarity indices
save(ring.range_size_Jaccard.80, file = "./outputs/Indices_maps/ring.range_size_Jaccard.80.RData", version = "2")
saveRDS(ring.range_size_Jaccard.80, file = "./outputs/Indices_maps/ring.range_size_Jaccard.80.rds", version = "2")
save(ring.rarity_Jaccard.80, file = "./outputs/Indices_maps/ring.rarity_Jaccard.80.RData", version = "2")
saveRDS(ring.rarity_Jaccard.80, file = "./outputs/Indices_maps/ring.rarity_Jaccard.80.rds", version = "2")

# Save full and mean rarity maps
save(mimicry.rarity_Jaccard.80, file = "./outputs/Indices_maps/mimicry.rarity_Jaccard.80.RData", version = "2")
saveRDS(mimicry.rarity_Jaccard.80, file = "./outputs/Indices_maps/mimicry.rarity_Jaccard.80.rds", version = "2")
save(mimicry.mean.rarity_Jaccard.80, file = "./outputs/Indices_maps/mimicry.mean.rarity_Jaccard.80.RData", version = "2")
saveRDS(mimicry.mean.rarity_Jaccard.80, file = "./outputs/Indices_maps/mimicry.mean.rarity_Jaccard.80.rds", version = "2")

# Plot distribution of ranges and rarity indices
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/ring_range_distri_Jaccard.80.pdf"), height = 5.3, width = 6.5)
hist(ring.range_size_Jaccard.80, main = "Distribution of mimicry ring geographic ranges\nJaccard.80", breaks = 20, xlab = "Geographic ranges in [10^3] km²")
range_threshold_Jaccard.80 <- round(quantile(ring.range_size_Jaccard.80, p=0.25), 1) # Range threshold value for qualitative rarity
abline(v = range_threshold_Jaccard.80, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",range_threshold_Jaccard.80," km²"), col = "red", lty = 2, lwd = 2, x = "topright", cex = 1, bty ="n")
dev.off()

pdf(file = paste0("./maps/Indices_maps/Jaccard.80/ring_rarity_indices_distri_Jaccard.80.pdf"), height = 5.3, width = 6.5)
hist(ring.rarity_Jaccard.80, main = "Distribution of mimicry ring rarity indices\nJaccard.80", breaks = 20, xlab = "Rarity indices")
indice_threshold_Jaccard.80 <- round(quantile(ring.rarity_Jaccard.80, p=0.75), 3) # Threshold value for qualitative rarity
abline(v = indice_threshold_Jaccard.80, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",indice_threshold_Jaccard.80), col = "red", lty = 2, lwd = 2, x = "topleft", cex = 1, bty ="n")
dev.off()


# Test other weighting scheme

# Load stuff to compute new rarity index
ring.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))
All_ring_proba_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_Jaccard.80.RData"))
load(file = "./outputs/Indices_maps/ring.range_size_Jaccard.80.RData")


library(Rarity)
?rWeights

# # Rarity weights based on inverse normalized range
# inv_weights <- rWeights(occData = ring.range_size_Jaccard.80, wMethods = "invQ", normalised = T, rounding = 5)$invQ
# hist(inv_weights)
# plot(inv_weights ~ ring.range_size_Jaccard.80)
# 
# plot(All_ring_proba_stack_Jaccard.80)

ring_proba_brick_Jaccard.80 <- brick(All_ring_proba_stack_Jaccard.80)
ring_assemblage <- t(tidyr::drop_na(as.data.frame(ring_proba_brick_Jaccard.80@data@values)))
ring_assemblage <- ring_assemblage[, colSums(ring_assemblage) >= 1]

# # Rarity weights based on an inverse exponential function with inflexion point calibrated as the rarity threshold as such that 25% of rings are rare
# Gaston_weights <- rWeights(occData = ring.range_size_Jaccard.80, wMethods = "W", rCutoff = "Gaston",
#                           normalised = T, rounding = 5, assemblages = ring_assemblage)$W
# hist(Gaston_weights)
# plot(Gaston_weights ~ ring.range_size_Jaccard.80)

# Rarity weights based on an inverse exponential function with inflexion point calibrated as the rarity threshold as such that communities host 25% of rare rings in average
Leroy_weights_df <- rWeights(occData = ring.range_size_Jaccard.80, wMethods = "W", rCutoff = "Leroy",
                             normalised = T, rounding = 5, assemblages = ring_assemblage)
Leroy_weights <- Leroy_weights_df$W

hist(Leroy_weights)
plot(Leroy_weights ~ ring.range_size_Jaccard.80)

# tested_weights <- inv_weights
# tested_weights <- Gaston_weights
tested_weights <- Leroy_weights

# Apply weights
tested_ring.weighted.stack_Jaccard.80 <- All_ring_proba_stack_Jaccard.80 # Generate the stack to fill
for (i in 1:nlayers(tested_ring.weighted.stack_Jaccard.80)){
  # Multiply ring probability of presence by ring rarity indices
  tested_ring.weighted.stack_Jaccard.80@layers[[i]]@data@values <- All_ring_proba_stack_Jaccard.80@layers[[i]]@data@values * tested_weights[i]
}

# Sum = Richness weighted by rarity based on Range-size
tested_ring.rarity_Jaccard.80 <- calc(tested_ring.weighted.stack_Jaccard.80, fun = sum)
# Divided by local richness = Mean ring rarity indices in the community (weighted by the probability of presence of each rings)
tested_ring.mean.rarity_Jaccard.80 <- tested_ring.rarity_Jaccard.80/ring.richness_Jaccard.80

tested_ring.mean.rarity_Jaccard.80 <- tested_ring.mean.rarity_Jaccard.80 + 0.005
# Add 0 values for empty continental pixels (transformed into min values as 0.10)
temp <- continent_mask
temp[!is.na(tested_ring.mean.rarity_Jaccard.80[])] <- tested_ring.mean.rarity_Jaccard.80[!is.na(tested_ring.mean.rarity_Jaccard.80[])]
tested_ring.mean.rarity_Jaccard.80 <- temp

tested_ring.mean.rarity_Jaccard.80_contrasted <- contrasting_raster(x = tested_ring.mean.rarity_Jaccard.80, zmin = 0, zmax = 0.75)

plot(tested_ring.mean.rarity_Jaccard.80, col = pal_bl_red_Mannion)
plot(tested_ring.mean.rarity_Jaccard.80_contrasted, col = pal_bl_red_Mannion)
plot(mimicry.mean.rarity_Jaccard.80, col = pal_bl_red_Mannion)

pdf(file = paste0("./maps/Indices_maps/Jaccard.80/Comparison_ring_rarity_weights.pdf"), height = 5.3, width = 5.3)
plot(tested_ring.mean.rarity_Jaccard.80[] ~ mimicry.mean.rarity_Jaccard.80[], ylab = "Leroy's weighted rarity", xlab = "Linearly weighted rarity", main = "Comparison of weighting schemes for rarity")
dev.off()
cor.test(x = tested_ring.mean.rarity_Jaccard.80[], y = mimicry.mean.rarity_Jaccard.80[], method = "spearman", alternative = "two.sided", na.rm = T)

plot(tested_ring.mean.rarity_Jaccard.80[] ~ ring.richness_Jaccard.80[])
cor.test(x = tested_ring.mean.rarity_Jaccard.80[], y = ring.richness_Jaccard.80[], method = "spearman", alternative = "two.sided", na.rm = T)

threshold_75 <- quantile(x = tested_ring.mean.rarity_Jaccard.80[], p = 0.75, na.rm = T)
threshold_95 <- quantile(x = tested_ring.mean.rarity_Jaccard.80[], p = 0.95, na.rm = T)

plot(tested_ring.mean.rarity_Jaccard.80 > threshold_75)
plot(tested_ring.mean.rarity_Jaccard.80 > threshold_95)

# Plot all stuff to compare the two weighting scheme
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/Comparison_ring_rarity_indices_all.pdf"), height = 7.5, width = 25)
par(mfrow = c(2,7))

plot(tested_weights ~ ring.range_size_Jaccard.80, xlab = "Ring range", ylab = "Ring rarity", main = "Leroy's weighting scheme")
hist(tested_weights, xlab = "Ring rarity weights", main = "Histo of weights")
hist(tested_ring.mean.rarity_Jaccard.80, xlab = 'Community mean geographic rarity', main = "Histo of community scores")
plot(tested_ring.mean.rarity_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = "Ring geographic rarity")
plot(tested_ring.mean.rarity_Jaccard.80 > quantile(x = tested_ring.mean.rarity_Jaccard.80[], p = 0.75, na.rm = T), main = "Top 25%")
plot(tested_ring.mean.rarity_Jaccard.80 > quantile(x = tested_ring.mean.rarity_Jaccard.80[], p = 0.95, na.rm = T), main = "Top 5%")
plot(tested_ring.mean.rarity_Jaccard.80[] ~ sp.richness_Jaccard.80[], xlab = "Sp. richness", ylab = "Ring geographic rarity", main = "Correlation with Sp. richness")

plot(ring.rarity_Jaccard.80 ~ ring.range_size_Jaccard.80, xlab = "Ring range", ylab = "Ring rarity", main = "Linear weighting scheme")
hist(ring.rarity_Jaccard.80, xlab = "Ring rarity weights", main = "Histo of weights")
hist(mimicry.mean.rarity_Jaccard.80, xlab = 'Community mean geographic rarity', main = "Histo of community scores")
plot(mimicry.mean.rarity_Jaccard.80, col = pal_bl_red_Mannion, main = "Ring geographic rarity")
plot(mimicry.mean.rarity_Jaccard.80 > quantile(x = mimicry.mean.rarity_Jaccard.80[], p = 0.75, na.rm = T), main = "Top 25%")
plot(mimicry.mean.rarity_Jaccard.80 > quantile(x = mimicry.mean.rarity_Jaccard.80[], p = 0.95, na.rm = T), main = "Top 5%")
plot(mimicry.mean.rarity_Jaccard.80[] ~ sp.richness_Jaccard.80[], xlab = "Sp. richness", ylab = "Ring geographic rarity", main = "Correlation with sp. richness")

par(mfrow = c(1,1))
dev.off()

# Plot the distribution of weights (ring rarity indices) for Leroy's weighting
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/ring_rarity_indices_distri_Jaccard.80_Leroy_weights.pdf"), height = 5.3, width = 6.5)
hist(Leroy_weights, main = "Distribution of ring rarity indices\nLeroy's weights", breaks = 20, xlab = "Rarity indices")
indice_threshold <- round(quantile(Leroy_weights, p = 0.75), 3) # Threshold value for Gaston's qualitative rarity
abline(v = 1-Leroy_weights_df$cut.off[1], col = "red", lty = 2, lwd = 2)
abline(v = indice_threshold, col = "red", lty = 3, lwd = 2)
legend(legend = c("Threshold for rarity", paste0("Leroy's scheme: ", round(1-Leroy_weights_df$cut.off[1], 3)), paste0("Gaston's 25%: ", indice_threshold)), col = "red", lty = c(NA, 2, 3), lwd = 2, x = "top", cex = 1, bty ="n")
dev.off()

# Save the stuff for the new rarity index based on Leroy (2013)'s weighting
ring.rarity_indices_Leroy_Jaccard.80 <- Leroy_weights
ring.rarity_Leroy_Jaccard.80 <- tested_ring.rarity_Jaccard.80
ring.mean.rarity_Leroy_Jaccard.80 <- tested_ring.mean.rarity_Jaccard.80
ring.mean.rarity_Leroy_Jaccard.80_contrasted <- tested_ring.mean.rarity_Jaccard.80_contrasted

# Save ring ranges and ring continuous rarity indices
save(ring.rarity_indices_Leroy_Jaccard.80, file = "./outputs/Indices_maps/ring.rarity_indices_Leroy_Jaccard.80.RData", version = "2")
saveRDS(ring.rarity_indices_Leroy_Jaccard.80, file = "./outputs/Indices_maps/ring.rarity_indices_Leroy_Jaccard.80.rds", version = "2")

# Save full and mean rarity maps
save(ring.rarity_Leroy_Jaccard.80, file = "./outputs/Indices_maps/ring.rarity_Leroy_Jaccard.80.RData", version = "2")
saveRDS(ring.rarity_Leroy_Jaccard.80, file = "./outputs/Indices_maps/ring.rarity_Leroy_Jaccard.80.rds", version = "2")
save(ring.mean.rarity_Leroy_Jaccard.80, file = "./outputs/Indices_maps/ring.mean.rarity_Leroy_Jaccard.80.RData", version = "2")
saveRDS(ring.mean.rarity_Leroy_Jaccard.80, file = "./outputs/Indices_maps/ring.mean.rarity_Leroy_Jaccard.80.rds", version = "2")
save(ring.mean.rarity_Leroy_Jaccard.80_contrasted, file = "./outputs/Indices_maps/ring.mean.rarity_Leroy_Jaccard.80_contrasted.RData", version = "2")
saveRDS(ring.mean.rarity_Leroy_Jaccard.80_contrasted, file = "./outputs/Indices_maps/ring.mean.rarity_Leroy_Jaccard.80_contrasted.rds", version = "2")


## For Jaccard.95

# Compute weights based on geographic ranges
ring.range_size_Jaccard.95 <- NA
for (i in 1:nlayers(All_ring_proba_stack_Jaccard.95)) {
  # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
  # Compute the estimated number of pixels occupied by the ring as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
  ring.range_size_Jaccard.95[i] <- sum(All_ring_proba_stack_Jaccard.95[[i]]@data@values, na.rm = T)*774.51/1000
}
ring.rarity_Jaccard.95 <- 1-(ring.range_size_Jaccard.95/max(ring.range_size_Jaccard.95)) # Weights by max ring extent = rarity indices

# Apply weights
ring.weighted.stack_Jaccard.95 <- All_ring_proba_stack_Jaccard.95 # Generate the stack to fill
for (i in 1:nlayers(ring.weighted.stack_Jaccard.95)){
  # Multiply ring probability of presence by ring rarity indices
  ring.weighted.stack_Jaccard.95@layers[[i]]@data@values <- All_ring_proba_stack_Jaccard.95@layers[[i]]@data@values * ring.rarity_Jaccard.95[i]
}

plot(ring.weighted.stack_Jaccard.95)

# Sum = Richness weighted by rarity based on Range-size
mimicry.rarity_Jaccard.95 <- calc(ring.weighted.stack_Jaccard.95, fun = sum) 
# Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each ring)
mimicry.mean.rarity_Jaccard.95 <- mimicry.rarity_Jaccard.95/ring.richness_Jaccard.95

# Add 0 values for empty continental pixels (transformed into min values as 0.10)
temp <- continent_mask
temp[!is.na(mimicry.mean.rarity_Jaccard.95[])] <- mimicry.mean.rarity_Jaccard.95[!is.na(mimicry.mean.rarity_Jaccard.95[])]
mimicry.mean.rarity_Jaccard.95 <- temp

# Save mimicry ring ranges and mimicry ring continous rarity indices
save(ring.range_size_Jaccard.95, file = "./outputs/Indices_maps/ring.range_size_Jaccard.95.RData", version = "2")
saveRDS(ring.range_size_Jaccard.95, file = "./outputs/Indices_maps/ring.range_size_Jaccard.95.rds", version = "2")
save(ring.rarity_Jaccard.95, file = "./outputs/Indices_maps/ring.rarity_Jaccard.95.RData", version = "2")
saveRDS(ring.rarity_Jaccard.95, file = "./outputs/Indices_maps/ring.rarity_Jaccard.95.rds", version = "2")

# Save full and mean rarity maps
save(mimicry.rarity_Jaccard.95, file = "./outputs/Indices_maps/mimicry.rarity_Jaccard.95.RData", version = "2")
saveRDS(mimicry.rarity_Jaccard.95, file = "./outputs/Indices_maps/mimicry.rarity_Jaccard.95.rds", version = "2")
save(mimicry.mean.rarity_Jaccard.95, file = "./outputs/Indices_maps/mimicry.mean.rarity_Jaccard.95.RData", version = "2")
saveRDS(mimicry.mean.rarity_Jaccard.95, file = "./outputs/Indices_maps/mimicry.mean.rarity_Jaccard.95.rds", version = "2")

# Plot distribution of ranges and rarity indices
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/ring_range_distri_Jaccard.95.pdf"), height = 5.3, width = 6.5)
hist(ring.range_size_Jaccard.95, main = "Distribution of mimicry ring geographic ranges\nJaccard.95", breaks = 20, xlab = "Geographic ranges in [10^3] km²")
range_threshold_Jaccard.95 <- round(quantile(ring.range_size_Jaccard.95, p=0.25), 1) # Range threshold value for qualitative rarity
abline(v = range_threshold_Jaccard.95, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",range_threshold_Jaccard.95," km²"), col = "red", lty = 2, lwd = 2, x = "topright", cex = 1, bty ="n")
dev.off()

pdf(file = paste0("./maps/Indices_maps/Jaccard.95/ring_rarity_indices_distri_Jaccard.95.pdf"), height = 5.3, width = 6.5)
hist(ring.rarity_Jaccard.95, main = "Distribution of mimicry ring rarity indices\nJaccard.95", breaks = 20, xlab = "Rarity indices")
indice_threshold_Jaccard.95 <- round(quantile(ring.rarity_Jaccard.95, p=0.75), 3) # Threshold value for qualitative rarity
abline(v = indice_threshold_Jaccard.95, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",indice_threshold_Jaccard.95), col = "red", lty = 2, lwd = 2, x = "topleft", cex = 1, bty ="n")
dev.off()


## For TSS.80

# Compute weights based on geographic ranges
ring.range_size_TSS.80 <- NA
for (i in 1:nlayers(All_ring_proba_stack_TSS.80)) {
  # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
  # Compute the estimated number of pixels occupied by the ring as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
  ring.range_size_TSS.80[i] <- sum(All_ring_proba_stack_TSS.80[[i]]@data@values, na.rm = T)*774.51/1000
}
ring.rarity_TSS.80 <- 1-(ring.range_size_TSS.80/max(ring.range_size_TSS.80)) # Weights by max ring extent = rarity indices

# Apply weights
ring.weighted.stack_TSS.80 <- All_ring_proba_stack_TSS.80 # Generate the stack to fill
for (i in 1:nlayers(ring.weighted.stack_TSS.80)){
  # Multiply ring probability of presence by ring rarity indices
  ring.weighted.stack_TSS.80@layers[[i]]@data@values <- All_ring_proba_stack_TSS.80@layers[[i]]@data@values * ring.rarity_TSS.80[i]
}

plot(ring.weighted.stack_TSS.80)

# Sum = Richness weighted by rarity based on Range-size
mimicry.rarity_TSS.80 <- calc(ring.weighted.stack_TSS.80, fun = sum) 
# Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each ring)
mimicry.mean.rarity_TSS.80 <- mimicry.rarity_TSS.80/ring.richness_TSS.80

# Add 0 values for empty continental pixels (transformed into min values as 0.10)
temp <- continent_mask
temp[!is.na(mimicry.mean.rarity_TSS.80[])] <- mimicry.mean.rarity_TSS.80[!is.na(mimicry.mean.rarity_TSS.80[])]
mimicry.mean.rarity_TSS.80 <- temp

# Save mimicry ring ranges and mimicry ring continous rarity indices
save(ring.range_size_TSS.80, file = "./outputs/Indices_maps/ring.range_size_TSS.80.RData", version = "2")
saveRDS(ring.range_size_TSS.80, file = "./outputs/Indices_maps/ring.range_size_TSS.80.rds", version = "2")
save(ring.rarity_TSS.80, file = "./outputs/Indices_maps/ring.rarity_TSS.80.RData", version = "2")
saveRDS(ring.rarity_TSS.80, file = "./outputs/Indices_maps/ring.rarity_TSS.80.rds", version = "2")

# Save full and mean rarity maps
save(mimicry.rarity_TSS.80, file = "./outputs/Indices_maps/mimicry.rarity_TSS.80.RData", version = "2")
saveRDS(mimicry.rarity_TSS.80, file = "./outputs/Indices_maps/mimicry.rarity_TSS.80.rds", version = "2")
save(mimicry.mean.rarity_TSS.80, file = "./outputs/Indices_maps/mimicry.mean.rarity_TSS.80.RData", version = "2")
saveRDS(mimicry.mean.rarity_TSS.80, file = "./outputs/Indices_maps/mimicry.mean.rarity_TSS.80.rds", version = "2")

# Plot distribution of ranges and rarity indices
pdf(file = paste0("./maps/Indices_maps/TSS.80/ring_range_distri_TSS.80.pdf"), height = 5.3, width = 6.5)
hist(ring.range_size_TSS.80, main = "Distribution of mimicry ring geographic ranges\nTSS.80", breaks = 20, xlab = "Geographic ranges in [10^3] km²")
range_threshold_TSS.80 <- round(quantile(ring.range_size_TSS.80, p=0.25), 1) # Range threshold value for qualitative rarity
abline(v = range_threshold_TSS.80, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",range_threshold_TSS.80," km²"), col = "red", lty = 2, lwd = 2, x = "topright", cex = 1, bty ="n")
dev.off()

pdf(file = paste0("./maps/Indices_maps/TSS.80/ring_rarity_indices_distri_TSS.80.pdf"), height = 5.3, width = 6.5)
hist(ring.rarity_TSS.80, main = "Distribution of mimicry ring rarity indices\nTSS.80", breaks = 20, xlab = "Rarity indices")
indice_threshold_TSS.80 <- round(quantile(ring.rarity_TSS.80, p=0.75), 3) # Threshold value for qualitative rarity
abline(v = indice_threshold_TSS.80, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",indice_threshold_TSS.80), col = "red", lty = 2, lwd = 2, x = "topleft", cex = 1, bty ="n")
dev.off()


## For TSS.95

# Compute weights based on geographic ranges
ring.range_size_TSS.95 <- NA
for (i in 1:nlayers(All_ring_proba_stack_TSS.95)) {
  # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
  # Compute the estimated number of pixels occupied by the ring as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
  ring.range_size_TSS.95[i] <- sum(All_ring_proba_stack_TSS.95[[i]]@data@values, na.rm = T)*774.51/1000
}
ring.rarity_TSS.95 <- 1-(ring.range_size_TSS.95/max(ring.range_size_TSS.95)) # Weights by max ring extent = rarity indices

# Apply weights
ring.weighted.stack_TSS.95 <- All_ring_proba_stack_TSS.95 # Generate the stack to fill
for (i in 1:nlayers(ring.weighted.stack_TSS.95)){
  # Multiply ring probability of presence by ring rarity indices
  ring.weighted.stack_TSS.95@layers[[i]]@data@values <- All_ring_proba_stack_TSS.95@layers[[i]]@data@values * ring.rarity_TSS.95[i]
}

plot(ring.weighted.stack_TSS.95)

# Sum = Richness weighted by rarity based on Range-size
mimicry.rarity_TSS.95 <- calc(ring.weighted.stack_TSS.95, fun = sum) 
# Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each ring)
mimicry.mean.rarity_TSS.95 <- mimicry.rarity_TSS.95/ring.richness_TSS.95

# Add 0 values for empty continental pixels (transformed into min values as 0.10)
temp <- continent_mask
temp[!is.na(mimicry.mean.rarity_TSS.95[])] <- mimicry.mean.rarity_TSS.95[!is.na(mimicry.mean.rarity_TSS.95[])]
mimicry.mean.rarity_TSS.95 <- temp

# Save mimicry ring ranges and mimicry ring continous rarity indices
save(ring.range_size_TSS.95, file = "./outputs/Indices_maps/ring.range_size_TSS.95.RData", version = "2")
saveRDS(ring.range_size_TSS.95, file = "./outputs/Indices_maps/ring.range_size_TSS.95.rds", version = "2")
save(ring.rarity_TSS.95, file = "./outputs/Indices_maps/ring.rarity_TSS.95.RData", version = "2")
saveRDS(ring.rarity_TSS.95, file = "./outputs/Indices_maps/ring.rarity_TSS.95.rds", version = "2")

# Save full and mean rarity maps
save(mimicry.rarity_TSS.95, file = "./outputs/Indices_maps/mimicry.rarity_TSS.95.RData", version = "2")
saveRDS(mimicry.rarity_TSS.95, file = "./outputs/Indices_maps/mimicry.rarity_TSS.95.rds", version = "2")
save(mimicry.mean.rarity_TSS.95, file = "./outputs/Indices_maps/mimicry.mean.rarity_TSS.95.RData", version = "2")
saveRDS(mimicry.mean.rarity_TSS.95, file = "./outputs/Indices_maps/mimicry.mean.rarity_TSS.95.rds", version = "2")

# Plot distribution of ranges and rarity indices
pdf(file = paste0("./maps/Indices_maps/TSS.95/ring_range_distri_TSS.95.pdf"), height = 5.3, width = 6.5)
hist(ring.range_size_TSS.95, main = "Distribution of mimicry ring geographic ranges\nTSS.95", breaks = 20, xlab = "Geographic ranges in [10^3] km²")
range_threshold_TSS.95 <- round(quantile(ring.range_size_TSS.95, p=0.25), 1) # Range threshold value for qualitative rarity
abline(v = range_threshold_TSS.95, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",range_threshold_TSS.95," km²"), col = "red", lty = 2, lwd = 2, x = "topright", cex = 1, bty ="n")
dev.off()

pdf(file = paste0("./maps/Indices_maps/TSS.95/ring_rarity_indices_distri_TSS.95.pdf"), height = 5.3, width = 6.5)
hist(ring.rarity_TSS.95, main = "Distribution of mimicry ring rarity indices\nTSS.95", breaks = 20, xlab = "Rarity indices")
indice_threshold_TSS.95 <- round(quantile(ring.rarity_TSS.95, p=0.75), 3) # Threshold value for qualitative rarity
abline(v = indice_threshold_TSS.95, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",indice_threshold_TSS.95), col = "red", lty = 2, lwd = 2, x = "topleft", cex = 1, bty ="n")
dev.off()

### 11.2/ Plot community rarity maps ####

### Load directly the community rarity maps
mimicry.rarity_Jaccard.80 <- readRDS(file = "./outputs/Indices_maps/mimicry.rarity_Jaccard.80.rds")
mimicry.rarity_Jaccard.95 <- readRDS(file = "./outputs/Indices_maps/mimicry.rarity_Jaccard.95.rds")
mimicry.rarity_TSS.80 <- readRDS(file = "./outputs/Indices_maps/mimicry.rarity_TSS.80.rds")
mimicry.rarity_TSS.95 <- readRDS(file = "./outputs/Indices_maps/mimicry.rarity_TSS.95.rds")


### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/ring.mimicry.rarity_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(mimicry.rarity_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mimicry rarity \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(mimicry.rarity_Jaccard.80, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/ring.mimicry.rarity_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(mimicry.rarity_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Mimicry rarity \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(mimicry.rarity_Jaccard.95, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")
par(mar = internal_margins)
dev.off()


# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/ring.mimicry.rarity_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(mimicry.rarity_TSS.80, col = pal_bl_red_Mannion, main = paste0("Mimicry rarity \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(mimicry.rarity_TSS.80, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")
par(mar = internal_margins)
dev.off()


# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/ring.mimicry.rarity_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(mimicry.rarity_TSS.95, col = pal_bl_red_Mannion, main = paste0("Mimicry rarity \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(mimicry.rarity_TSS.95, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(mimicry.rarity_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/ring.mimicry.rarity_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(mimicry.rarity_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mimicry rarity \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(mimicry.rarity_Jaccard.80, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")

image(mimicry.rarity_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Mimicry rarity \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(mimicry.rarity_Jaccard.95, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")

image(mimicry.rarity_TSS.80, col = pal_bl_red_Mannion, main = paste0("Mimicry rarity \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(mimicry.rarity_TSS.80, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")

image(mimicry.rarity_TSS.95, col = pal_bl_red_Mannion, main = paste0("Mimicry rarity \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(mimicry.rarity_TSS.95, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")

par(mar = internal_margins)

dev.off()


### 11.3/ Plot mean mimicry rarity maps ####

### Load directly the mean mimicry rarity maps
mimicry.mean.rarity_Jaccard.80 <- readRDS(file = "./outputs/Indices_maps/mimicry.mean.rarity_Jaccard.80.rds")
mimicry.mean.rarity_Jaccard.95 <- readRDS(file = "./outputs/Indices_maps/mimicry.mean.rarity_Jaccard.95.rds")
mimicry.mean.rarity_TSS.80 <- readRDS(file = "./outputs/Indices_maps/mimicry.mean.rarity_TSS.80.rds")
mimicry.mean.rarity_TSS.95 <- readRDS(file = "./outputs/Indices_maps/mimicry.mean.rarity_TSS.95.rds")


### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/ring.mimicry.mean.rarity_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(mimicry.mean.rarity_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mean mimicry rarity\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(mimicry.mean.rarity_Jaccard.80, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/ring.mimicry.mean.rarity_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(mimicry.mean.rarity_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Mean mimicry rarity\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(mimicry.mean.rarity_Jaccard.95, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")
par(mar = internal_margins)
dev.off()


# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/ring.mimicry.mean.rarity_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(mimicry.mean.rarity_TSS.80, col = pal_bl_red_Mannion, main = paste0("Mean mimicry rarity\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(mimicry.mean.rarity_TSS.80, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/ring.mimicry.mean.rarity_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(mimicry.mean.rarity_TSS.95, col = pal_bl_red_Mannion, main = paste0("Mean mimicry rarity\nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(mimicry.mean.rarity_TSS.95, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(mimicry.mean.rarity_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/ring.mimicry.mean.rarity_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(mimicry.mean.rarity_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mean mimicry rarity\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(mimicry.mean.rarity_Jaccard.80, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")

image(mimicry.mean.rarity_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Mean mimicry rarity\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(mimicry.mean.rarity_Jaccard.95, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")

image(mimicry.mean.rarity_TSS.80, col = pal_bl_red_Mannion, main = paste0("Mean mimicry rarity\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(mimicry.mean.rarity_TSS.80, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")

image(mimicry.mean.rarity_TSS.95, col = pal_bl_red_Mannion, main = paste0("Mean mimicry rarity\nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(mimicry.mean.rarity_TSS.95, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")

par(mar = internal_margins)

dev.off()


# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])



###### 12/ Community vulnerability ######

# Function to compute community vulnerability as the sum of rings vulnerability approximated as the inverse of the local richness of each mimicry ring
# A ring with only one species as a vulnerability of 1/1 = 1. A ring of 4 species as a vulnerability of 1/4 = 0.25.
# Final vulnerability is standardized by total ring richness such as a community mean vulnerability

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


### Load directly the complete stacks of mimicry ring richness
All_ring_rich_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_Jaccard.80.RData"))
All_ring_rich_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_Jaccard.95.RData"))
All_ring_rich_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_TSS.80.RData"))
All_ring_rich_stack_TSS.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_TSS.95.RData"))

### 12.1/ Compute indices ####

vulnerability_Jaccard.80 <- calc(x = All_ring_rich_stack_Jaccard.80, fun = vulnerability)*1
vulnerability_Jaccard.95 <- calc(x = All_ring_rich_stack_Jaccard.95, fun = vulnerability)*1
vulnerability_TSS.80 <- calc(x = All_ring_rich_stack_TSS.80, fun = vulnerability)*1
vulnerability_TSS.95 <- calc(x = All_ring_rich_stack_TSS.95, fun = vulnerability)*1

# Save
save(vulnerability_Jaccard.80, file = paste0("./outputs/Indices_maps/vulnerability_Jaccard.80.RData"), version = "2")
saveRDS(vulnerability_Jaccard.80, file = paste0("./outputs/Indices_maps/vulnerability_Jaccard.80.rds"), version = "2")
save(vulnerability_Jaccard.95, file = paste0("./outputs/Indices_maps/vulnerability_Jaccard.95.RData"), version = "2")
saveRDS(vulnerability_Jaccard.95, file = paste0("./outputs/Indices_maps/vulnerability_Jaccard.95.rds"), version = "2")
save(vulnerability_TSS.80, file = paste0("./outputs/Indices_maps/vulnerability_TSS.80.RData"), version = "2")
saveRDS(vulnerability_TSS.80, file = paste0("./outputs/Indices_maps/vulnerability_TSS.80.rds"), version = "2")
save(vulnerability_TSS.95, file = paste0("./outputs/Indices_maps/vulnerability_TSS.95.RData"), version = "2")
saveRDS(vulnerability_TSS.95, file = paste0("./outputs/Indices_maps/vulnerability_TSS.95.rds"), version = "2")

### 12.2/ Plot Vulnerability ####

# Load directly the final Ring richness layer
vulnerability_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/vulnerability_Jaccard.80.rds"))
vulnerability_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_maps/vulnerability_Jaccard.95.rds"))
vulnerability_TSS.80 <- readRDS(file = paste0("./outputs/Indices_maps/vulnerability_TSS.80.rds"))
vulnerability_TSS.95 <- readRDS(file = paste0("./outputs/Indices_maps/vulnerability_TSS.95.rds"))

hist(vulnerability_Jaccard.80)

# Contrast values by merging low values below 0.10
vulnerability_Jaccard.80_contrasted <- vulnerability_Jaccard.80
vulnerability_Jaccard.80_contrasted[vulnerability_Jaccard.80 < 0.20] <- 0.205
vulnerability_Jaccard.95_contrasted <- vulnerability_Jaccard.95
vulnerability_Jaccard.95_contrasted[vulnerability_Jaccard.95 < 0.20] <- 0.205
vulnerability_TSS.80_contrasted <- vulnerability_TSS.80
vulnerability_TSS.80_contrasted[vulnerability_TSS.80 < 0.20] <- 0.205
vulnerability_TSS.95_contrasted <- vulnerability_TSS.95
vulnerability_TSS.95_contrasted[vulnerability_TSS.95 < 0.20] <- 0.205

# Add 0 values for empty continental pixels (transformed into min values as 0.20)
temp <- continent_mask
temp[!is.na(vulnerability_Jaccard.80_contrasted[])] <- vulnerability_Jaccard.80_contrasted[!is.na(vulnerability_Jaccard.80_contrasted)]
vulnerability_Jaccard.80_contrasted <- temp
vulnerability_Jaccard.80_contrasted[vulnerability_Jaccard.80_contrasted == 0] <- 0.20

temp <- continent_mask
temp[!is.na(vulnerability_Jaccard.95_contrasted[])] <- vulnerability_Jaccard.95_contrasted[!is.na(vulnerability_Jaccard.95_contrasted)]
vulnerability_Jaccard.95_contrasted <- temp
vulnerability_Jaccard.95_contrasted[vulnerability_Jaccard.95_contrasted == 0] <- 0.20

temp <- continent_mask
temp[!is.na(vulnerability_TSS.80_contrasted[])] <- vulnerability_TSS.80_contrasted[!is.na(vulnerability_TSS.80_contrasted)]
vulnerability_TSS.80_contrasted <- temp
vulnerability_TSS.80_contrasted[vulnerability_TSS.80_contrasted == 0] <- 0.20

temp <- continent_mask
temp[!is.na(vulnerability_TSS.95_contrasted[])] <- vulnerability_TSS.95_contrasted[!is.na(vulnerability_TSS.95_contrasted)]
vulnerability_TSS.95_contrasted <- temp
vulnerability_TSS.95_contrasted[vulnerability_TSS.95_contrasted == 0] <- 0.20

hist(vulnerability_Jaccard.80_contrasted)

### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/vulnerability_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(vulnerability_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Community vulnerability\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Vulnerability\n         index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(vulnerability_Jaccard.80_contrasted, locs = seq(0.2, 1, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 7, font = 2, cex = 1.2, label = "Vulnerability\nindex")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/vulnerability_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(vulnerability_Jaccard.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Community vulnerability\nJaccard.95"), 
     cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
     ylab = "", xlab = "",
     legend.args=list(text="          Vulnerability\n         index", cex=1.2, line = 1, font = 2), 
     legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(vulnerability_Jaccard.95_contrasted, locs = seq(0.2, 1, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 7, font = 2, cex = 1.2, label = "Vulnerability\nindex")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/vulnerability_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(vulnerability_TSS.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Community vulnerability\nTSS.80"), 
     cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
     ylab = "", xlab = "",
     legend.args=list(text="          Vulnerability\n         index", cex=1.2, line = 1, font = 2), 
     legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(vulnerability_TSS.80_contrasted, locs = seq(0.2, 1, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 7, font = 2, cex = 1.2, label = "Vulnerability\nindex")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/vulnerability_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(vulnerability_TSS.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Community vulnerability\nTSS.95"), 
     cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
     ylab = "", xlab = "",
     legend.args=list(text="          Vulnerability\n         index", cex=1.2, line = 1, font = 2), 
     legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(vulnerability_TSS.95_contrasted, locs = seq(0.2, 1, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 7, font = 2, cex = 1.2, label = "Vulnerability\nindex")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(vulnerability_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/vulnerability_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(vulnerability_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Community vulnerability\nJaccard.80"), 
     cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
     ylab = "", xlab = "",
     legend.args=list(text="          Vulnerability\n         index", cex=1.2, line = 1, font = 2), 
     legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(vulnerability_Jaccard.80_contrasted, locs = seq(0.2, 1, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 7, font = 2, cex = 1.2, label = "Vulnerability\nindex")

image(vulnerability_Jaccard.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Community vulnerability\nJaccard.95"), 
     cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
     ylab = "", xlab = "",
     legend.args=list(text="          Vulnerability\n         index", cex=1.2, line = 1, font = 2), 
     legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(vulnerability_Jaccard.95_contrasted, locs = seq(0.2, 1, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 7, font = 2, cex = 1.2, label = "Vulnerability\nindex")

image(vulnerability_TSS.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Community vulnerability\nTSS.80"), 
     cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
     ylab = "", xlab = "",
     legend.args=list(text="          Vulnerability\n         index", cex=1.2, line = 1, font = 2), 
     legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(vulnerability_TSS.80_contrasted, locs = seq(0.2, 1, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 7, font = 2, cex = 1.2, label = "Vulnerability\nindex")

image(vulnerability_TSS.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Community vulnerability\nTSS.95"), 
     cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
     ylab = "", xlab = "",
     legend.args=list(text="          Vulnerability\n         index", cex=1.2, line = 1, font = 2), 
     legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(vulnerability_TSS.95_contrasted, locs = seq(0.2, 1, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 7, font = 2, cex = 1.2, label = "Vulnerability\nindex")

par(mar = internal_margins)

dev.off()


# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])




#################################### Phylogeny-based Indices ##############################################

###### 13/ Load common stuff for phylogeny based indices #####

# vignette("picante-intro")

library(ape)
library(picante)
library(geiger)
library(treespace)

# ?read.tree
# ?read.nexus

# ?cophenetic # To shift from tree to patristic distance matrix (sum of connecting branch lengths)
# ?vcv # to shift from tree to variance-covariance matrix (shared branch lengths)

### Prepare phylogeny

# # Load phylogeny
# phylo.Ithomiini <- read.nexus(file = "./input_data/Phylogenies/Ithomiini_MCC_annotated_outgroups.nex")
# 
# plot(phylo.Ithomiini)
# 
# # Load extra-group list
# EG_phylo <- read.csv2(file = "./input_data/Phylogenies/List_extra_groups_phylo.csv", header = F)
# EG_phylo <- as.character(EG_phylo$V1)
# # Prune extra-groups
# phylo.Ithomiini <- drop.tip(phy = phylo.Ithomiini, tip = EG_phylo)
# 
# # Rename taxa using format of the database (i.e., with dots between Genera and species)
# labels <- phylo.Ithomiini$tip.label
# sapply(X = labels, FUN = strsplit, split = "_")
# split.labels <- strsplit(x = labels, split = "_") 
# new.labels <- unlist(lapply(X = split.labels, FUN = paste, collapse ="."))
# 
# # Special cases to handle manually
# setdiff(new.labels, list.sp$Sp_full)
# # Renommer les taxons avec des "_"
# new.labels[new.labels == "Ithomia.terra.EAST"] <- "Ithomia.terraEAST"
# new.labels[new.labels == "Pseudoscada.timna.COSTARICA"] <- "Pseudoscada.timnaCOSTARICA"
# new.labels[new.labels == "Pseudoscada.timna.EASTERN"] <- "Pseudoscada.timnaEASTERN"
# new.labels[new.labels == "Pseudoscada.timna.WESTERN"] <- "Pseudoscada.timnaWESTERN"
# new.labels[new.labels == "Pteronymia.veia.WEST"] <- "Pteronymia.veiaWEST"
# new.labels[new.labels == "Pteronymia.veia.EAST"] <- "Pteronymia.veiaEAST"
# # Paste new labels on the phylogeny
# phylo.Ithomiini$tip.label <- new.labels
# # Remove the mistake Hypothyris.sp from the phylogeny
# phylo.Ithomiini <- drop.tip(phy = phylo.Ithomiini, tip = "Hypothyris.sp")
# 
# save(phylo.Ithomiini, file = paste0(internal.wd, "/Phylogenies/Final_phylogeny.RData"))

### Create function to aggregate probabilities to higher hierachical level (aggregate pixel, or go up on a phylogenetic tree)
aggreg_prob = function(x, na.rm) { 
  y <- 1-prod(1-x) 
  return(y) # Output
}


### Load directly the phylogeny with the good labels
load(file = "./input_data/Phylogenies/Final_phylogeny.RData")

phylo.Ithomiini <- readRDS(file = "./input_data/Phylogenies/Final_phylogeny.rds")


### Load directly the stacks of sp probas
All_sp_proba_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))
All_sp_proba_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.95.rds"))
All_sp_proba_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80.rds"))
All_sp_proba_stack_TSS.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.95.rds"))

### Need to extract only the species in the phylogeny

# Load Summary table for OMU/unit and for Species
load(file = paste0("./input_data/list.sp.RData"))

# Extract only the 339 species included in the phylogeny from list.sp and the sp.stack
list.sp_phylo <- list.sp[list.sp$Sp_full %in% phylo.Ithomiini$tip.label,]

# Extract only the 339 species included in the phylogeny from the stacks of sp probas
All_sp_proba_stack_Jaccard.80_phylo <- All_sp_proba_stack_Jaccard.80[[phylo.Ithomiini$tip.label]]
All_sp_proba_stack_Jaccard.95_phylo <- All_sp_proba_stack_Jaccard.95[[phylo.Ithomiini$tip.label]]
All_sp_proba_stack_TSS.80_phylo <- All_sp_proba_stack_TSS.80[[phylo.Ithomiini$tip.label]]
All_sp_proba_stack_TSS.95_phylo <- All_sp_proba_stack_TSS.95[[phylo.Ithomiini$tip.label]]

# Save stacks
save(All_sp_proba_stack_Jaccard.80_phylo, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80_phylo.RData"), version = "2")
saveRDS(All_sp_proba_stack_Jaccard.80_phylo, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80_phylo.rds"), version = "2")
save(All_sp_proba_stack_Jaccard.95_phylo, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.95_phylo.RData"), version = "2")
saveRDS(All_sp_proba_stack_Jaccard.95_phylo, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.95_phylo.rds"), version = "2")
save(All_sp_proba_stack_TSS.80_phylo, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80_phylo.RData"), version = "2")
saveRDS(All_sp_proba_stack_TSS.80_phylo, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80_phylo.rds"), version = "2")
save(All_sp_proba_stack_TSS.95_phylo, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.95_phylo.RData"), version = "2")
saveRDS(All_sp_proba_stack_TSS.95_phylo, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.95_phylo.rds"), version = "2")

# Load directly the sp proba stack with only species in the phylogeny
All_sp_proba_stack_Jaccard.80_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80_phylo.rds"))
All_sp_proba_stack_Jaccard.90_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.90_phylo.rds"))
All_sp_proba_stack_TSS.80_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80_phylo.rds"))
All_sp_proba_stack_TSS.90_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.90_phylo.rds"))


###### 14/ Mean pairwise Phylogenetic Distance (MPD) ######

?mpd # Cannot be used because the weighting scheme does not apply to our probability of presence of pairs of species 

# Get matrix of probas per species (columns) per community (lines)
sp.proba.mat_Jaccard.80 <- readAll(brick(All_sp_proba_stack_Jaccard.80_phylo))@data@values
sp.proba.mat_Jaccard.95 <- readAll(brick(All_sp_proba_stack_Jaccard.95_phylo))@data@values
sp.proba.mat_TSS.80 <- readAll(brick(All_sp_proba_stack_TSS.80_phylo))@data@values
sp.proba.mat_TSS.95 <- readAll(brick(All_sp_proba_stack_TSS.95_phylo))@data@values

# Save
save(sp.proba.mat_Jaccard.80, file = paste0("./outputs/Indices_stacks/sp.proba.mat_Jaccard.80.RData"), version = "2")
saveRDS(sp.proba.mat_Jaccard.80, file = paste0("./outputs/Indices_stacks/sp.proba.mat_Jaccard.80.rds"), version = "2")
save(sp.proba.mat_Jaccard.95, file = paste0("./outputs/Indices_stacks/sp.proba.mat_Jaccard.95.RData"), version = "2")
saveRDS(sp.proba.mat_Jaccard.95, file = paste0("./outputs/Indices_stacks/sp.proba.mat_Jaccard.95.rds"), version = "2")
save(sp.proba.mat_TSS.80, file = paste0("./outputs/Indices_stacks/sp.proba.mat_TSS.80.RData"), version = "2")
saveRDS(sp.proba.mat_TSS.80, file = paste0("./outputs/Indices_stacks/sp.proba.mat_TSS.80.rds"), version = "2")
save(sp.proba.mat_TSS.95, file = paste0("./outputs/Indices_stacks/sp.proba.mat_TSS.95.RData"), version = "2")
saveRDS(sp.proba.mat_TSS.95, file = paste0("./outputs/Indices_stacks/sp.proba.mat_TSS.95.rds"), version = "2")

# Compute the patristic phylogenetic distance
phylo.dist.mat <- cophenetic.phylo(x = phylo.Ithomiini)
# Get total number of species
n.leaves <- length(phylo.Ithomiini$tip.label)

### 14.1/ Indices computation ####

# Load matrix of species probas x communities
sp.proba.mat_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/sp.proba.mat_Jaccard.80.rds"))
sp.proba.mat_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/sp.proba.mat_Jaccard.95.rds"))
sp.proba.mat_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/sp.proba.mat_TSS.80.rds"))
sp.proba.mat_TSS.95 <- readRDS(file = paste0("./outputs/Indices_stacks/sp.proba.mat_TSS.95.rds"))


# Loop to compute index per community for Jaccard.80
all.MPD_Jaccard.80 <- NA
for (k in 1:nrow(sp.proba.mat_Jaccard.80)) {
  sp.proba.com <- sp.proba.mat_Jaccard.80[k,] # Extract the proba for the k community
  
  # Compute MPD only if no NA is present in the community
  if (any(is.na(sp.proba.com))) {
    
    all.MPD_Jaccard.80[k] <- NA # if NA present, MPD = NA
    
  } else {
    
    # Compute a matrix of proba of existence of pairs in the community
    proba.pairs.mat <- matrix(data = NA, ncol = n.leaves, nrow = n.leaves)
    for (i in 1:n.leaves) {
      for (j in 1:n.leaves) {
        pair.probas <- c(sp.proba.com[i], sp.proba.com[j])
        proba.pairs.mat[i,j] <- prod(pair.probas) # Probability of presence of a pair = product of proba. of presence of each of the species of the pair
      }
    }
    
    # Multiply them to ponderate phylogenetic distance by probability of presence of the pair of species
    weighted.pd.mat <- proba.pairs.mat*phylo.dist.mat
    
    # Extract unique pair values and get the mean
    total.phy.dist <- sum(weighted.pd.mat[lower.tri(weighted.pd.mat, diag = F)]) # Extract weighted pair values and sum them
    MPD <- total.phy.dist/sum(proba.pairs.mat[lower.tri(proba.pairs.mat, diag = F)]) # Divide by the sum of the weights to compute the weighted mean
    all.MPD_Jaccard.80[k] <- MPD # Store the MDP of each community in the final vector
    
  }
  # Show k every 1000 iterations and save a back-up file
  if (k%%1000==0) {
    cat(paste0(Sys.time(), " - ", k," on ",nrow(sp.proba.mat_Jaccard.80),"\n"))
    save(all.MPD_Jaccard.80, file = "./outputs/Indices_Maps/MPD_backup_Jaccard.80.RData", version = "2")
  }
}

# Write community MPD in a raster
MPD.raster_Jaccard.80 <- continent_mask
MPD.raster_Jaccard.80@data@values <-  all.MPD_Jaccard.80

# Save
save(MPD.raster_Jaccard.80, file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80.RData", version = "2")
saveRDS(MPD.raster_Jaccard.80, file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80.rds", version = "2")


# Loop to compute index per community for Jaccard.95
all.MPD_Jaccard.95 <- NA
for (k in 1:nrow(sp.proba.mat_Jaccard.95)) {
  sp.proba.com <- sp.proba.mat_Jaccard.95[k,] # Extract the proba for the k community
  
  # Compute MPD only if no NA is present in the community
  if (any(is.na(sp.proba.com))) {
    
    all.MPD_Jaccard.95[k] <- NA # if NA present, MPD = NA
    
  } else {
    
    # Compute a matrix of proba of existence of pairs in the community
    proba.pairs.mat <- matrix(data = NA, ncol = n.leaves, nrow = n.leaves)
    for (i in 1:n.leaves) {
      for (j in 1:n.leaves) {
        pair.probas <- c(sp.proba.com[i], sp.proba.com[j])
        proba.pairs.mat[i,j] <- prod(pair.probas) # Probability of presence of a pair = product of proba. of presence of each of the species of the pair
      }
    }
    
    # Multiply them to ponderate phylogenetic distance by probability of presence of the pair of species
    weighted.pd.mat <- proba.pairs.mat*phylo.dist.mat
    
    # Extract unique pair values and get the mean
    total.phy.dist <- sum(weighted.pd.mat[lower.tri(weighted.pd.mat, diag = F)]) # Extract weighted pair values and sum them
    MPD <- total.phy.dist/sum(proba.pairs.mat[lower.tri(proba.pairs.mat, diag = F)]) # Divide by the sum of the weights to compute the weighted mean
    all.MPD_Jaccard.95[k] <- MPD # Store the MDP of each community in the final vector
    
  }
  # Show k every 1000 iterations and save a back-up file
  if (k%%1000==0) {
    cat(paste0(Sys.time(), " - ", k," on ",nrow(sp.proba.mat_Jaccard.95),"\n"))
    save(all.MPD_Jaccard.95, file = "./outputs/Indices_Maps/MPD_backup_Jaccard.95.RData", version = "2")
  }
}

# Write community MPD in a raster
MPD.raster_Jaccard.95 <- continent_mask
MPD.raster_Jaccard.95@data@values <-  all.MPD_Jaccard.95
# plot(MPD.raster_Jaccard.95)

# Save
save(MPD.raster_Jaccard.95, file = "./outputs/Indices_Maps/MPD.raster_Jaccard.95.RData", version = "2")
saveRDS(MPD.raster_Jaccard.95, file = "./outputs/Indices_Maps/MPD.raster_Jaccard.95.rds", version = "2")


# Loop to compute index per community for TSS.80
all.MPD_TSS.80 <- NA
for (k in 1:nrow(sp.proba.mat_TSS.80)) {
  sp.proba.com <- sp.proba.mat_TSS.80[k,] # Extract the proba for the k community
  
  # Compute MPD only if no NA is present in the community
  if (any(is.na(sp.proba.com))) {
    
    all.MPD_TSS.80[k] <- NA # if NA present, MPD = NA
    
  } else {
    
    # Compute a matrix of proba of existence of pairs in the community
    proba.pairs.mat <- matrix(data = NA, ncol = n.leaves, nrow = n.leaves)
    for (i in 1:n.leaves) {
      for (j in 1:n.leaves) {
        pair.probas <- c(sp.proba.com[i], sp.proba.com[j])
        proba.pairs.mat[i,j] <- prod(pair.probas) # Probability of presence of a pair = product of proba. of presence of each of the species of the pair
      }
    }
    
    # Multiply them to ponderate phylogenetic distance by probability of presence of the pair of species
    weighted.pd.mat <- proba.pairs.mat*phylo.dist.mat
    
    # Extract unique pair values and get the mean
    total.phy.dist <- sum(weighted.pd.mat[lower.tri(weighted.pd.mat, diag = F)]) # Extract weighted pair values and sum them
    MPD <- total.phy.dist/sum(proba.pairs.mat[lower.tri(proba.pairs.mat, diag = F)]) # Divide by the sum of the weights to compute the weighted mean
    all.MPD_TSS.80[k] <- MPD # Store the MDP of each community in the final vector
    
  }
  # Show k every 1000 iterations and save a back-up file
  if (k%%1000==0) {
    cat(paste0(Sys.time(), " - ", k," on ",nrow(sp.proba.mat_TSS.80),"\n"))
    save(all.MPD_TSS.80, file = "./outputs/Indices_Maps/MPD_backup_TSS.80.RData", version = "2")
  }
}

# Write community MPD in a raster
MPD.raster_TSS.80 <- continent_mask
MPD.raster_TSS.80@data@values <-  all.MPD_TSS.80
# plot(MPD.raster_TSS.80)

# Save
save(MPD.raster_TSS.80, file = "./outputs/Indices_Maps/MPD.raster_TSS.80.RData", version = "2")
saveRDS(MPD.raster_TSS.80, file = "./outputs/Indices_Maps/MPD.raster_TSS.80.rds", version = "2")


# Loop to compute index per community for TSS.95
all.MPD_TSS.95 <- NA
for (k in 1:nrow(sp.proba.mat_TSS.95)) {
  sp.proba.com <- sp.proba.mat_TSS.95[k,] # Extract the proba for the k community
  
  # Compute MPD only if no NA is present in the community
  if (any(is.na(sp.proba.com))) {
    
    all.MPD_TSS.95[k] <- NA # if NA present, MPD = NA
    
  } else {
    
    # Compute a matrix of proba of existence of pairs in the community
    proba.pairs.mat <- matrix(data = NA, ncol = n.leaves, nrow = n.leaves)
    for (i in 1:n.leaves) {
      for (j in 1:n.leaves) {
        pair.probas <- c(sp.proba.com[i], sp.proba.com[j])
        proba.pairs.mat[i,j] <- prod(pair.probas) # Probability of presence of a pair = product of proba. of presence of each of the species of the pair
      }
    }
    
    # Multiply them to ponderate phylogenetic distance by probability of presence of the pair of species
    weighted.pd.mat <- proba.pairs.mat*phylo.dist.mat
    
    # Extract unique pair values and get the mean
    total.phy.dist <- sum(weighted.pd.mat[lower.tri(weighted.pd.mat, diag = F)]) # Extract weighted pair values and sum them
    MPD <- total.phy.dist/sum(proba.pairs.mat[lower.tri(proba.pairs.mat, diag = F)]) # Divide by the sum of the weights to compute the weighted mean
    all.MPD_TSS.95[k] <- MPD # Store the MDP of each community in the final vector
    
  }
  # Show k every 1000 iterations and save a back-up file
  if (k%%1000==0) {
    cat(paste0(Sys.time(), " - ", k," on ",nrow(sp.proba.mat_TSS.95),"\n"))
    save(all.MPD_TSS.95, file = "./outputs/Indices_Maps/MPD_backup_TSS.95.RData", version = "2")
  }
}

# Write community MPD in a raster
MPD.raster_TSS.95 <- continent_mask
MPD.raster_TSS.95@data@values <-  all.MPD_TSS.95
# plot(MPD.raster_TSS.95)

# Save
save(MPD.raster_TSS.95, file = "./outputs/Indices_Maps/MPD.raster_TSS.95.RData", version = "2")
saveRDS(MPD.raster_TSS.95, file = "./outputs/Indices_Maps/MPD.raster_TSS.95.rds", version = "2")


### 14.2/ Contrast rasters for plots ####

### Load directly the final MPD layers
MPD.raster_Jaccard.80 <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80.rds")
MPD.raster_Jaccard.95 <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.95.rds")
MPD.raster_TSS.80 <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_TSS.80.rds")
MPD.raster_TSS.95 <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_TSS.95.rds")

hist(MPD.raster_Jaccard.80)
quantile(MPD.raster_Jaccard.80[], na.rm = T, probs = seq(0, 1, 0.025))


# Contrast values by merging low values below 36.5 and above 45
MPD.raster_Jaccard.80_contrasted <- MPD.raster_Jaccard.80
MPD.raster_Jaccard.80_contrasted[MPD.raster_Jaccard.80 < 36.5] <- 36.7
MPD.raster_Jaccard.80_contrasted[MPD.raster_Jaccard.80 > 45] <- 45
MPD.raster_Jaccard.95_contrasted <- MPD.raster_Jaccard.95
MPD.raster_Jaccard.95_contrasted[MPD.raster_Jaccard.95 < 36.5] <- 36.7
MPD.raster_Jaccard.95_contrasted[MPD.raster_Jaccard.95 > 45] <- 45
MPD.raster_TSS.80_contrasted <- MPD.raster_TSS.80
MPD.raster_TSS.80_contrasted[MPD.raster_TSS.80 < 36.5] <- 36.7
MPD.raster_TSS.80_contrasted[MPD.raster_TSS.80 > 45] <- 45
MPD.raster_TSS.95_contrasted <- MPD.raster_TSS.95
MPD.raster_TSS.95_contrasted[MPD.raster_TSS.95 < 36.5] <- 36.7
MPD.raster_TSS.95_contrasted[MPD.raster_TSS.95 > 45] <- 45


# Add 0 values for empty continental pixels (transformed into min values as 36.5)
temp <- continent_mask
temp[!is.na(MPD.raster_Jaccard.80_contrasted[])] <- MPD.raster_Jaccard.80_contrasted[!is.na(MPD.raster_Jaccard.80_contrasted[])]
MPD.raster_Jaccard.80_contrasted <- temp
MPD.raster_Jaccard.80_contrasted[MPD.raster_Jaccard.80_contrasted == 0] <- 36.5

temp <- continent_mask
temp[!is.na(MPD.raster_Jaccard.95_contrasted[])] <- MPD.raster_Jaccard.95_contrasted[!is.na(MPD.raster_Jaccard.95_contrasted[])]
MPD.raster_Jaccard.95_contrasted <- temp
MPD.raster_Jaccard.95_contrasted[MPD.raster_Jaccard.95_contrasted == 0] <- 36.5

temp <- continent_mask
temp[!is.na(MPD.raster_TSS.80_contrasted[])] <- MPD.raster_TSS.80_contrasted[!is.na(MPD.raster_TSS.80_contrasted[])]
MPD.raster_TSS.80_contrasted <- temp
MPD.raster_TSS.80_contrasted[MPD.raster_TSS.80_contrasted == 0] <- 36.5

temp <- continent_mask
temp[!is.na(MPD.raster_TSS.95_contrasted[])] <- MPD.raster_TSS.95_contrasted[!is.na(MPD.raster_TSS.95_contrasted[])]
MPD.raster_TSS.95_contrasted <- temp
MPD.raster_TSS.95_contrasted[MPD.raster_TSS.95_contrasted == 0] <- 36.5

hist(MPD.raster_Jaccard.80_contrasted)

# Save contrasted rasters
save(MPD.raster_Jaccard.80_contrasted, file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80_contrasted.RData", version = "2")
saveRDS(MPD.raster_Jaccard.80_contrasted, file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80_contrasted.rds", version = "2")
save(MPD.raster_Jaccard.95_contrasted, file = "./outputs/Indices_Maps/MPD.raster_Jaccard.95_contrasted.RData", version = "2")
saveRDS(MPD.raster_Jaccard.95_contrasted, file = "./outputs/Indices_Maps/MPD.raster_Jaccard.95_contrasted.rds", version = "2")
save(MPD.raster_TSS.80_contrasted, file = "./outputs/Indices_Maps/MPD.raster_TSS.80_contrasted.RData", version = "2")
saveRDS(MPD.raster_TSS.80_contrasted, file = "./outputs/Indices_Maps/MPD.raster_TSS.80_contrasted.rds", version = "2")
save(MPD.raster_TSS.95_contrasted, file = "./outputs/Indices_Maps/MPD.raster_TSS.95_contrasted.RData", version = "2")
saveRDS(MPD.raster_TSS.95_contrasted, file = "./outputs/Indices_Maps/MPD.raster_TSS.95_contrasted.rds", version = "2")


### 14.3/ Plot full MPD maps ####

# Load directly the contrasted rasters to plot
MPD.raster_Jaccard.80_contrasted <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80_contrasted.rds")
MPD.raster_Jaccard.95_contrasted <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.95_contrasted.rds")
MPD.raster_TSS.80_contrasted <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_TSS.80_contrasted.rds")
MPD.raster_TSS.95_contrasted <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_TSS.95_contrasted.rds")

### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/MPD.raster_Jaccard.80_contrasted.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
image(MPD.raster_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean pairwise Phylogenetic Distance\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(MPD.raster_Jaccard.80_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/MPD.raster_Jaccard.95_contrasted.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
image(MPD.raster_Jaccard.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean pairwise Phylogenetic Distance\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(MPD.raster_Jaccard.95_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/MPD.raster_TSS.80_contrasted.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
image(MPD.raster_TSS.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean pairwise Phylogenetic Distance\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(MPD.raster_TSS.80_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/MPD.raster_TSS.95_contrasted.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
image(MPD.raster_TSS.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean pairwise Phylogenetic Distance\nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(MPD.raster_TSS.95_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(MPD.raster_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/MPD.raster_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(MPD.raster_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean pairwise Phylogenetic Distance\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(MPD.raster_Jaccard.80_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(MPD.raster_Jaccard.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean pairwise Phylogenetic Distance\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(MPD.raster_Jaccard.95_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(MPD.raster_TSS.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean pairwise Phylogenetic Distance\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(MPD.raster_TSS.80_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(MPD.raster_TSS.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean pairwise Phylogenetic Distance\nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(MPD.raster_TSS.95_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

par(mar = internal_margins)

dev.off()


### 14.4/ Plot MPD maps zoomed on the Andes ####

# Load directly the contrasted rasters to plot
MPD.raster_Jaccard.80_contrasted <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80_contrasted.rds")
MPD.raster_Jaccard.95_contrasted <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.95_contrasted.rds")
MPD.raster_TSS.80_contrasted <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_TSS.80_contrasted.rds")
MPD.raster_TSS.95_contrasted <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_TSS.95_contrasted.rds")

# Crop to Andes region
library(rgeos)
Andes_ext <- extent(c(-90, -59, -15, 14))
Andes.crop_mask_shp <- crop(crop_mask_shp, Andes_ext)

Andes.MPD.raster_Jaccard.80_contrasted <- crop(MPD.raster_Jaccard.80_contrasted, Andes_ext)
Andes.MPD.raster_Jaccard.95_contrasted <- crop(MPD.raster_Jaccard.95_contrasted, Andes_ext)
Andes.MPD.raster_TSS.80_contrasted <- crop(MPD.raster_TSS.80_contrasted, Andes_ext)
Andes.MPD.raster_TSS.95_contrasted <- crop(MPD.raster_TSS.95_contrasted, Andes_ext)


# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/Andes.MPD.raster_Jaccard.80_contrasted.pdf"), height = 5.3, width = 6)
internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
image(Andes.MPD.raster_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean pairwise Phylogenetic Distance\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(Andes.MPD.raster_Jaccard.80_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-88.22, -87.2, -13, 1))
graphics::text(x = -85.5, y = 4, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/Andes.MPD.raster_Jaccard.95_contrasted.pdf"), height = 5.3, width = 6)
internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
image(Andes.MPD.raster_Jaccard.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean pairwise Phylogenetic Distance\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(Andes.MPD.raster_Jaccard.95_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-88.22, -87.2, -13, 1))
graphics::text(x = -85.5, y = 4, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/Andes.MPD.raster_TSS.80_contrasted.pdf"), height = 5.3, width = 6)
internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
image(Andes.MPD.raster_TSS.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean pairwise Phylogenetic Distance\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(Andes.MPD.raster_TSS.80_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-88.22, -87.2, -13, 1))
graphics::text(x = -85.5, y = 4, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/Andes.MPD.raster_TSS.95_contrasted.pdf"), height = 5.3, width = 6)
internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
image(Andes.MPD.raster_TSS.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean pairwise Phylogenetic Distance\nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(Andes.MPD.raster_TSS.95_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-88.22, -87.2, -13, 1))
graphics::text(x = -85.5, y = 4, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(Andes.MPD.raster_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/Andes.MPD.raster_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(Andes.MPD.raster_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean pairwise Phylogenetic Distance\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(Andes.MPD.raster_Jaccard.80_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-88.22, -87.2, -13, 1))
graphics::text(x = -85.5, y = 4, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(Andes.MPD.raster_Jaccard.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean pairwise Phylogenetic Distance\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(Andes.MPD.raster_Jaccard.95_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-88.22, -87.2, -13, 1))
graphics::text(x = -85.5, y = 4, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(Andes.MPD.raster_TSS.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean pairwise Phylogenetic Distance\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(Andes.MPD.raster_TSS.80_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-88.22, -87.2, -13, 1))
graphics::text(x = -85.5, y = 4, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(Andes.MPD.raster_TSS.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean pairwise Phylogenetic Distance\nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(Andes.MPD.raster_TSS.95_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-88.22, -87.2, -13, 1))
graphics::text(x = -85.5, y = 4, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

par(mar = internal_margins)

dev.off()



# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])


###### 15/ Faith's Phylogenetic Diversity #####

?pd # To compute Faith's phylogenetic distance (but not working with probabilities)

# Generate matrix to store info on edges/branches
branches = matrix(NA, nrow(phylo.Ithomiini$edge), ncol = 4) 
names(branches) <- c("Starting nod", "Ending nod", "Length", "Proba_presence")
branches[, 1:2] = phylo.Ithomiini$edge # Retrieve starting and ending node 
branches[, 3] = round(phylo.Ithomiini$edge.length, 4) # Retrieve edge length

# Load matrix of species probas x communities
sp.proba.mat_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/sp.proba.mat_Jaccard.80.rds"))
sp.proba.mat_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/sp.proba.mat_Jaccard.95.rds"))
sp.proba.mat_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/sp.proba.mat_TSS.80.rds"))
sp.proba.mat_TSS.95 <- readRDS(file = paste0("./outputs/Indices_stacks/sp.proba.mat_TSS.95.rds"))

### 15.1/ Compute indices ####

# For each community, compute probability of presence of the edge as the probability for at least one descending species to be in the community

### Compute index for Jaccard.80

# Make a loop per community
all.PD_Jaccard.80 <- NA
for (k in 1:nrow(sp.proba.mat_Jaccard.80)) {
  sp.proba.com <- sp.proba.mat_Jaccard.80[k,] # Extract the proba for the k community
  
  # Compute PD only if no NA is present in the community
  if (any(is.na(sp.proba.com))) {
    
    all.PD_Jaccard.80[k] <- NA # if NA present, PD = NA
    
  } else {
    
    # Compute probability of presence of each edge/branch (i.e., a species descending from this branch) in this community
    for (i in 1:nrow(branches)) {
      leaves.node = tips(phylo.Ithomiini, branches[i, 2]) # Retrieve the set of species descending from branch i
      index <- which(colnames(sp.proba.mat_Jaccard.80) %in% leaves.node) # Find position of the species in the stack/matrix
      prob_edge <- aggreg_prob(sp.proba.com[index]) # Compute proba of presence of the edge
      branches[i, 4] <- prob_edge # Store info
    }
    PD <- round(sum(branches[,3]*branches[,4]),4) # Compute PD as the weighted sum of the edge length (weighted by the probabilty of presence of this edge in the community)
    all.PD_Jaccard.80[k]  <- PD
  }
  
  # Show k every 100 iterations and save a backup
  if (k %% 1000 == 0) {
    cat(paste0(Sys.time(), " - ", k," on ",nrow(sp.proba.mat_Jaccard.80),"\n"))
    save(all.PD_Jaccard.80, file = "./outputs/Indices_Maps/PD_backup_Jaccard.80.RData", version = "2")
  }
}

# Write PD values in a raster
PD.raster_Jaccard.80 <- continent_mask
PD.raster_Jaccard.80@data@values <-  all.PD_Jaccard.80
# plot(PD.raster_Jaccard.80)

# Repair issue with max value
PD.raster_Jaccard.80@data@max <- max(PD.raster_Jaccard.80[], na.rm = T)

# Save
save(PD.raster_Jaccard.80, file = "./outputs/Indices_Maps/PD.raster_Jaccard.80.RData", version = "2")
saveRDS(PD.raster_Jaccard.80, file = "./outputs/Indices_Maps/PD.raster_Jaccard.80.rds", version = "2")


### Compute index for Jaccard.95

# Make a loop per community
all.PD_Jaccard.95 <- NA
for (k in 1:nrow(sp.proba.mat_Jaccard.95)) {
  sp.proba.com <- sp.proba.mat_Jaccard.95[k,] # Extract the proba for the k community
  
  # Compute PD only if no NA is present in the community
  if (any(is.na(sp.proba.com))) {
    
    all.PD_Jaccard.95[k] <- NA # if NA present, PD = NA
    
  } else {
    
    # Compute probability of presence of each edge/branch (i.e., a species descending from this branch) in this community
    for (i in 1:nrow(branches)) {
      leaves.node = tips(phylo.Ithomiini, branches[i, 2]) # Retrieve the set of species descending from branch i
      index <- which(colnames(sp.proba.mat_Jaccard.95) %in% leaves.node) # Find position of the species in the stack/matrix
      prob_edge <- aggreg_prob(sp.proba.com[index]) # Compute proba of presence of the edge
      branches[i, 4] <- prob_edge # Store info
    }
    PD <- round(sum(branches[,3]*branches[,4]), 4) # Compute PD as the weighted sum of the edge length (weighted by the probabilty of presence of this edge in the community)
    all.PD_Jaccard.95[k]  <- PD
  }
  
  # Show k every 100 iterations and save a backup
  if (k %% 1000 == 0) {
    cat(paste0(Sys.time(), " - ", k," on ",nrow(sp.proba.mat_Jaccard.95),"\n"))
    save(all.PD_Jaccard.95, file = "./outputs/Indices_Maps/PD_backup_Jaccard.95.RData", version = "2")
  }
}

# Write PD values in a raster
PD.raster_Jaccard.95 <- continent_mask
PD.raster_Jaccard.95@data@values <-  all.PD_Jaccard.95
# plot(PD.raster_Jaccard.95)

# Repair issue with max value
PD.raster_Jaccard.95@data@max <- max(PD.raster_Jaccard.95[], na.rm = T)

# Save
save(PD.raster_Jaccard.95, file = "./outputs/Indices_Maps/PD.raster_Jaccard.95.RData", version = "2")
saveRDS(PD.raster_Jaccard.95, file = "./outputs/Indices_Maps/PD.raster_Jaccard.95.rds", version = "2")

### Compute index for TSS.80

# Make a loop per community
all.PD_TSS.80 <- NA
for (k in 1:nrow(sp.proba.mat_TSS.80)) {
  sp.proba.com <- sp.proba.mat_TSS.80[k,] # Extract the proba for the k community
  
  # Compute PD only if no NA is present in the community
  if (any(is.na(sp.proba.com))) {
    
    all.PD_TSS.80[k] <- NA # if NA present, PD = NA
    
  } else {
    
    # Compute probability of presence of each edge/branch (i.e., a species descending from this branch) in this community
    for (i in 1:nrow(branches)) {
      leaves.node = tips(phylo.Ithomiini, branches[i, 2]) # Retrieve the set of species descending from branch i
      index <- which(colnames(sp.proba.mat_TSS.80) %in% leaves.node) # Find position of the species in the stack/matrix
      prob_edge <- aggreg_prob(sp.proba.com[index]) # Compute proba of presence of the edge
      branches[i, 4] <- prob_edge # Store info
    }
    PD <- round(sum(branches[,3]*branches[,4]), 4) # Compute PD as the weighted sum of the edge length (weighted by the probabilty of presence of this edge in the community)
    all.PD_TSS.80[k]  <- PD
  }
  
  # Show k every 100 iterations and save a backup
  if (k %% 1000 == 0) {
    cat(paste0(Sys.time(), " - ", k," on ",nrow(sp.proba.mat_TSS.80),"\n"))
    save(all.PD_TSS.80, file = "./outputs/Indices_Maps/PD_backup_TSS.80.RData", version = "2")
  }
}

# Write PD values in a raster
PD.raster_TSS.80 <- continent_mask
PD.raster_TSS.80@data@values <-  all.PD_TSS.80
# plot(PD.raster_TSS.80)

# Repair issue with max value
PD.raster_TSS.80@data@max <- max(PD.raster_TSS.80[], na.rm = T)

# Save
save(PD.raster_TSS.80, file = "./outputs/Indices_Maps/PD.raster_TSS.80.RData", version = "2")
saveRDS(PD.raster_TSS.80, file = "./outputs/Indices_Maps/PD.raster_TSS.80.rds", version = "2")


### Compute index for TSS.95

# Make a loop per community
all.PD_TSS.95 <- NA
for (k in 1:nrow(sp.proba.mat_TSS.95)) {
  sp.proba.com <- sp.proba.mat_TSS.95[k,] # Extract the proba for the k community
  
  # Compute PD only if no NA is present in the community
  if (any(is.na(sp.proba.com))) {
    
    all.PD_TSS.95[k] <- NA # if NA present, PD = NA
    
  } else {
    
    # Compute probability of presence of each edge/branch (i.e., a species descending from this branch) in this community
    for (i in 1:nrow(branches)) {
      leaves.node = tips(phylo.Ithomiini, branches[i, 2]) # Retrieve the set of species descending from branch i
      index <- which(colnames(sp.proba.mat_TSS.95) %in% leaves.node) # Find position of the species in the stack/matrix
      prob_edge <- aggreg_prob(sp.proba.com[index]) # Compute proba of presence of the edge
      branches[i, 4] <- prob_edge # Store info
    }
    PD <- round(sum(branches[,3]*branches[,4]), 4) # Compute PD as the weighted sum of the edge length (weighted by the probabilty of presence of this edge in the community)
    all.PD_TSS.95[k]  <- PD
  }
  
  # Show k every 100 iterations and save a backup
  if (k %% 1000 == 0) {
    cat(paste0(Sys.time(), " - ", k," on ",nrow(sp.proba.mat_TSS.95),"\n"))
    save(all.PD_TSS.95, file = "./outputs/Indices_Maps/PD_backup_TSS.95.RData", version = "2")
  }
}

# Write PD values in a raster
PD.raster_TSS.95 <- continent_mask
PD.raster_TSS.95@data@values <-  all.PD_TSS.95
# plot(PD.raster_TSS.95)

# Repair issue with max value
PD.raster_TSS.95@data@max <- max(PD.raster_TSS.95[], na.rm = T)

# Save
save(PD.raster_TSS.95, file = "./outputs/Indices_Maps/PD.raster_TSS.95.RData", version = "2")
saveRDS(PD.raster_TSS.95, file = "./outputs/Indices_Maps/PD.raster_TSS.95.rds", version = "2")


### 15.2/ Plot PD maps ####

### Load directly the final PD layer
PD.raster_Jaccard.80 <- readRDS(file = "./outputs/Indices_Maps/PD.raster_Jaccard.80.rds")
PD.raster_Jaccard.95 <- readRDS(file = "./outputs/Indices_Maps/PD.raster_Jaccard.95.rds")
PD.raster_TSS.80 <- readRDS(file = "./outputs/Indices_Maps/PD.raster_TSS.80.rds")
PD.raster_TSS.95 <- readRDS(file = "./outputs/Indices_Maps/PD.raster_TSS.95.rds")

### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/PD.raster_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
image(PD.raster_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Faith's Phylogenetic Diversity\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(PD.raster_Jaccard.80, locs = seq(0, 800, 200), minmax = c(0, max(PD.raster_Jaccard.80[], na.rm = T)), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -107, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/PD.raster_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
image(PD.raster_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Faith's Phylogenetic Diversity\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(PD.raster_Jaccard.95, locs = seq(0, 800, 200), minmax = c(0, max(PD.raster_Jaccard.95[], na.rm = T)), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -107, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/PD.raster_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
image(PD.raster_TSS.80, col = pal_bl_red_Mannion, main = paste0("Faith's Phylogenetic Diversity\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(PD.raster_TSS.80, locs = seq(0, 800, 200), minmax = c(0, max(PD.raster_TSS.80[], na.rm = T)), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -107, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/PD.raster_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
image(PD.raster_TSS.95, col = pal_bl_red_Mannion, main = paste0("Faith's Phylogenetic Diversity\nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(PD.raster_TSS.95, locs = seq(0, 800, 200), minmax = c(0, max(PD.raster_TSS.95[], na.rm = T)), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -107, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(MPD.raster_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/PD.raster_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(PD.raster_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Faith's Phylogenetic Diversity\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(PD.raster_Jaccard.80, locs = seq(0, 800, 200), minmax = c(0, max(PD.raster_Jaccard.80[], na.rm = T)), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -107, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(PD.raster_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Faith's Phylogenetic Diversity\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(PD.raster_Jaccard.95, locs = seq(0, 800, 200), minmax = c(0, max(PD.raster_Jaccard.95[], na.rm = T)), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -107, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(PD.raster_TSS.80, col = pal_bl_red_Mannion, main = paste0("Faith's Phylogenetic Diversity\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(PD.raster_TSS.80, locs = seq(0, 800, 200), minmax = c(0, max(PD.raster_TSS.80[], na.rm = T)), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -107, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(PD.raster_TSS.95, col = pal_bl_red_Mannion, main = paste0("Faith's Phylogenetic Diversity\nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(PD.raster_TSS.95, locs = seq(0, 800, 200), minmax = c(0, max(PD.raster_TSS.95[], na.rm = T)), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -107, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

par(mar = internal_margins)

dev.off()


###### 16/ Equal-Splits & Fair-proportions ######

# Equal-Splits : distribute evolutionary time of each edge/branchs equally among its descending lineages
# Fair-Proportion : distribute evolutionary time of each edge/branchs equally among its descending species

# At species level = Sum of all modified branches length from the root to the species. The most original species have the highest scores
# At community level = Sum of ES/FP values of each species present in the community. The community with the most original species has the highest scores
# Can be standardized by the number of species to get the mean species ES/FP index

### 16.1/ Generate matrix to store info on edges ####

branches = matrix(NA, nrow(phylo.Ithomiini$edge), ncol = 7) 
branches[, 1:2] = phylo.Ithomiini$edge # Retrieve starting and ending node 
branches[, 3] = phylo.Ithomiini$edge.length # Retrieve edge length
for (i in 1:nrow(phylo.Ithomiini$edge)) { # Retrieve number of descending species/leaves (terminal nods)
  branches[i,4] = length(tips(phylo.Ithomiini, branches[i,2]))
}
for (i in 1:nrow(phylo.Ithomiini$edge)) { # Retrieve number of descending lineages (nb of branches branching out of the ending nod)
  branches[i,5] = sum(branches[,1]==branches[i,2])
}
branches[(branches[,5] == 0),5] <- 1 # Put value 1 for terminal branches

branches[,6] = branches[,3]/branches[,4] # Divide branch length by nb of descendant species for FP
branches[,7] = branches[,3]/branches[,5] # Divide branch length by nb of descendant lineages for ES

branches <- data.frame(branches)
names(branches) <- c("Starting_node","Ending_node","Branch_Length","Nb.desc.sp","Nb.desc.lineages","FP","ES")


### 16.2/ Generate df to store info on species FP and ES values ####

# (depend only on the phylogeny, not the community composition)

# Use species list with only species present in the phylogeny
Species_ES_FP_index_table <- data.frame(Sp.full = list.sp_phylo$Sp_full, FP = rep(NA,nrow(list.sp_phylo)), ES = rep(NA,nrow(list.sp_phylo)))
# Species_ES_FP_index_table <- data.frame(Sp.full = phylo.test$tip.label, FP = rep(NA,length(phylo.test$tip.label)), ES = rep(NA,length(phylo.test$tip.label)))
names.phylo <- phylo.Ithomiini$tip.label # Names of the leaves ordered by their index in the phylogeny

for (i in 1:length(names.phylo))
{
  # Retrieve sp name
  sp <- names.phylo[i] 
  
  # Get indices of nods on the path from the root to the species terminal nod
  sp.path <- nodepath(phy = phylo.Ithomiini, from = branches$Starting_node[which.max(branches$Nb.desc.sp)], to = which(names.phylo == sp))
  
  # Compute FP index for the species by summing FP values of the branches on the path to the species tip in the phylogeny 
  Species_ES_FP_index_table$FP[Species_ES_FP_index_table$Sp.full == sp] <- sum(branches[which(branches$Ending_node %in% sp.path), "FP"]) 
  
  # Compute ES index for the species by summing ES values of the branches on the path to the species tip in the phylogeny
  Species_ES_FP_index_table$ES[Species_ES_FP_index_table$Sp.full == sp] <- sum(branches[which(branches$Ending_node %in% sp.path), "ES"]) 
}

# Save the df for ES/FP species values
save(Species_ES_FP_index_table, file = "./outputs/Indices_Maps/Species_ES_FP_index_table.RData", version = "2")
# cor(Species_ES_FP_index_table$FP, Species_ES_FP_index_table$ES) # Correlation btw FP and ES indices = 0.888

### 16.3/ Indices computation ####

# For each community, compute FP and ES index by weighteing species index by their probability of presences

# Load directly the Species table of ES and FP index
load(file = "./outputs/Indices_Maps/Species_ES_FP_index_table.RData")

# Load matrix of species probas x communities
sp.proba.mat_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/sp.proba.mat_Jaccard.80.rds"))
sp.proba.mat_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/sp.proba.mat_Jaccard.95.rds"))
sp.proba.mat_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/sp.proba.mat_TSS.80.rds"))
sp.proba.mat_TSS.95 <- readRDS(file = paste0("./outputs/Indices_stacks/sp.proba.mat_TSS.95.rds"))


### Compute for Jaccard.80

# Make a loop per community
all.FP_Jaccard.80 <- all.mean.FP_Jaccard.80 <- all.ES_Jaccard.80 <- all.mean.ES_Jaccard.80 <- rep(NA, nrow(sp.proba.mat_Jaccard.80))
for (k in 1:nrow(sp.proba.mat_Jaccard.80)) 
{
  sp.proba.com <- sp.proba.mat_Jaccard.80[k,] # Extract the proba for the k community
  
  # Compute FP and ES only if no NA is present in the community
  if (any(is.na(sp.proba.com))) 
  {
    all.FP_Jaccard.80[k] <- NA # if NA present, FP = NA
    all.ES_Jaccard.80[k] <- NA # if NA present, ES = NA
    all.mean.FP_Jaccard.80[k] <- NA # if NA present, mean FP = NA
    all.mean.ES_Jaccard.80[k] <- NA # if NA present, mean ES = NA
    
  } else {
    FP <- sum(Species_ES_FP_index_table$FP*sp.proba.com) # Compute FP as the weighted sum of the FP index of species, weighted by their probability of presence in the community
    all.FP_Jaccard.80[k]  <- FP
    mean.FP <- FP/sum(sp.proba.com) # Compute mean FP as the weighted mean of the FP index of species, weighted by their probability of presence in the community
    all.mean.FP_Jaccard.80[k]  <- mean.FP
    ES <- sum(Species_ES_FP_index_table$ES*sp.proba.com) # Compute ES as the weighted sum of the ES index of species, weighted by their probability of presence in the community
    all.ES_Jaccard.80[k]  <- ES
    mean.ES <- ES/sum(sp.proba.com) # Compute mean ES as the weighted mean of the ES index of species, weighted by their probability of presence in the community
    all.mean.ES_Jaccard.80[k]  <- mean.ES
  }
  
  # Show k every 1000 iterations and save a back-up file
  if (k %% 1000 == 0) {
    cat(paste0(Sys.time(), " - ", k," on ",nrow(sp.proba.mat_Jaccard.80),"\n"))
    save(all.FP_Jaccard.80, all.mean.FP_Jaccard.80, all.ES_Jaccard.80, all.mean.ES_Jaccard.80, file = "./outputs/Indices_maps/backup.FP.ES_Jaccard.80.RData")
  }
  if(k == nrow(sp.proba.mat_Jaccard.80)) {
    cat(paste0(Sys.time(), " - ", k," on ",nrow(sp.proba.mat_Jaccard.80),"\n"))
    save(all.FP_Jaccard.80, all.mean.FP_Jaccard.80, all.ES_Jaccard.80, all.mean.ES_Jaccard.80, file = "./outputs/Indices_maps/backup.FP.ES_Jaccard.80.RData")
  }
}

# Write FP values in a raster
FP.raster_Jaccard.80 <- continent_mask
FP.raster_Jaccard.80@data@values <-  all.FP_Jaccard.80
FP.raster_Jaccard.80@data@max <- max(FP.raster_Jaccard.80[], na.rm = T) 
# plot(FP.raster_Jaccard.80)

# Write FP values in a raster
FP.mean.raster_Jaccard.80 <- continent_mask
FP.mean.raster_Jaccard.80@data@values <-  all.mean.FP_Jaccard.80
FP.mean.raster_Jaccard.80@data@max <- max(FP.mean.raster_Jaccard.80[], na.rm = T)
# plot(FP.mean.raster_Jaccard.80)

# Write ES values in a raster
ES.raster_Jaccard.80 <- continent_mask
ES.raster_Jaccard.80@data@values <-  all.ES_Jaccard.80
ES.raster_Jaccard.80@data@max <- max(ES.raster_Jaccard.80[], na.rm = T) 
# plot(ES.raster_Jaccard.80)

# Write ES values in a raster
ES.mean.raster_Jaccard.80 <- continent_mask
ES.mean.raster_Jaccard.80@data@values <-  all.mean.ES_Jaccard.80
ES.mean.raster_Jaccard.80@data@max <- max(ES.mean.raster_Jaccard.80[], na.rm = T)
# plot(ES.mean.raster_Jaccard.80)

# Save
save(FP.raster_Jaccard.80, file = "./outputs/Indices_Maps/FP.raster_Jaccard.80.RData", version = "2")
saveRDS(FP.raster_Jaccard.80, file = "./outputs/Indices_Maps/FP.raster_Jaccard.80.rds", version = "2")
save(FP.mean.raster_Jaccard.80, file = "./outputs/Indices_Maps/FP.mean.raster_Jaccard.80.RData", version = "2")
saveRDS(FP.mean.raster_Jaccard.80, file = "./outputs/Indices_Maps/FP.mean.raster_Jaccard.80.rds", version = "2")
save(ES.raster_Jaccard.80, file = "./outputs/Indices_Maps/ES.raster_Jaccard.80.RData", version = "2")
saveRDS(ES.raster_Jaccard.80, file = "./outputs/Indices_Maps/ES.raster_Jaccard.80.rds", version = "2")
save(ES.mean.raster_Jaccard.80, file = "./outputs/Indices_Maps/ES.mean.raster_Jaccard.80.RData", version = "2")
saveRDS(ES.mean.raster_Jaccard.80, file = "./outputs/Indices_Maps/ES.mean.raster_Jaccard.80.rds", version = "2")


### Compute for Jaccard.95

# Make a loop per community
all.FP_Jaccard.95 <- all.mean.FP_Jaccard.95 <- all.ES_Jaccard.95 <- all.mean.ES_Jaccard.95 <- rep(NA, nrow(sp.proba.mat_Jaccard.95))
for (k in 1:nrow(sp.proba.mat_Jaccard.95)) 
{
  sp.proba.com <- sp.proba.mat_Jaccard.95[k,] # Extract the proba for the k community
  
  # Compute FP and ES only if no NA is present in the community
  if (any(is.na(sp.proba.com))) 
  {
    all.FP_Jaccard.95[k] <- NA # if NA present, FP = NA
    all.ES_Jaccard.95[k] <- NA # if NA present, ES = NA
    all.mean.FP_Jaccard.95[k] <- NA # if NA present, mean FP = NA
    all.mean.ES_Jaccard.95[k] <- NA # if NA present, mean ES = NA
    
  } else {
    FP <- sum(Species_ES_FP_index_table$FP*sp.proba.com) # Compute FP as the weighted sum of the FP index of species, weighted by their probability of presence in the community
    all.FP_Jaccard.95[k]  <- FP
    mean.FP <- FP/sum(sp.proba.com) # Compute mean FP as the weighted mean of the FP index of species, weighted by their probability of presence in the community
    all.mean.FP_Jaccard.95[k]  <- mean.FP
    ES <- sum(Species_ES_FP_index_table$ES*sp.proba.com) # Compute ES as the weighted sum of the ES index of species, weighted by their probability of presence in the community
    all.ES_Jaccard.95[k]  <- ES
    mean.ES <- ES/sum(sp.proba.com) # Compute mean ES as the weighted mean of the ES index of species, weighted by their probability of presence in the community
    all.mean.ES_Jaccard.95[k]  <- mean.ES
  }
  
  # Show k every 1000 iterations and save a back-up file
  if (k %% 1000 == 0) {
    cat(paste0(Sys.time(), " - ", k," on ",nrow(sp.proba.mat_Jaccard.95),"\n"))
    save(all.FP_Jaccard.95, all.mean.FP_Jaccard.95, all.ES_Jaccard.95, all.mean.ES_Jaccard.95, file = "./outputs/Indices_maps/backup.FP.ES_Jaccard.95.RData")
  }
  if(k == nrow(sp.proba.mat_Jaccard.95)) {
    cat(paste0(Sys.time(), " - ", k," on ",nrow(sp.proba.mat_Jaccard.95),"\n"))
    save(all.FP_Jaccard.95, all.mean.FP_Jaccard.95, all.ES_Jaccard.95, all.mean.ES_Jaccard.95, file = "./outputs/Indices_maps/backup.FP.ES_Jaccard.95.RData")
  }
}

# Write FP values in a raster
FP.raster_Jaccard.95 <- continent_mask
FP.raster_Jaccard.95@data@values <-  all.FP_Jaccard.95
FP.raster_Jaccard.95@data@max <- max(FP.raster_Jaccard.95[], na.rm = T) 
# plot(FP.raster_Jaccard.95)

# Write FP values in a raster
FP.mean.raster_Jaccard.95 <- continent_mask
FP.mean.raster_Jaccard.95@data@values <-  all.mean.FP_Jaccard.95
FP.mean.raster_Jaccard.95@data@max <- max(FP.mean.raster_Jaccard.95[], na.rm = T)
# plot(FP.mean.raster_Jaccard.95)

# Write ES values in a raster
ES.raster_Jaccard.95 <- continent_mask
ES.raster_Jaccard.95@data@values <-  all.ES_Jaccard.95
ES.raster_Jaccard.95@data@max <- max(ES.raster_Jaccard.95[], na.rm = T) 
# plot(ES.raster_Jaccard.95)

# Write ES values in a raster
ES.mean.raster_Jaccard.95 <- continent_mask
ES.mean.raster_Jaccard.95@data@values <-  all.mean.ES_Jaccard.95
ES.mean.raster_Jaccard.95@data@max <- max(ES.mean.raster_Jaccard.95[], na.rm = T)
# plot(ES.mean.raster_Jaccard.95)

# Save
save(FP.raster_Jaccard.95, file = "./outputs/Indices_Maps/FP.raster_Jaccard.95.RData", version = "2")
saveRDS(FP.raster_Jaccard.95, file = "./outputs/Indices_Maps/FP.raster_Jaccard.95.rds", version = "2")
save(FP.mean.raster_Jaccard.95, file = "./outputs/Indices_Maps/FP.mean.raster_Jaccard.95.RData", version = "2")
saveRDS(FP.mean.raster_Jaccard.95, file = "./outputs/Indices_Maps/FP.mean.raster_Jaccard.95.rds", version = "2")
save(ES.raster_Jaccard.95, file = "./outputs/Indices_Maps/ES.raster_Jaccard.95.RData", version = "2")
saveRDS(ES.raster_Jaccard.95, file = "./outputs/Indices_Maps/ES.raster_Jaccard.95.rds", version = "2")
save(ES.mean.raster_Jaccard.95, file = "./outputs/Indices_Maps/ES.mean.raster_Jaccard.95.RData", version = "2")
saveRDS(ES.mean.raster_Jaccard.95, file = "./outputs/Indices_Maps/ES.mean.raster_Jaccard.95.rds", version = "2")


### Compute for TSS.80

# Make a loop per community
all.FP_TSS.80 <- all.mean.FP_TSS.80 <- all.ES_TSS.80 <- all.mean.ES_TSS.80 <- rep(NA, nrow(sp.proba.mat_TSS.80))
for (k in 1:nrow(sp.proba.mat_TSS.80)) 
{
  sp.proba.com <- sp.proba.mat_TSS.80[k,] # Extract the proba for the k community
  
  # Compute FP and ES only if no NA is present in the community
  if (any(is.na(sp.proba.com))) 
  {
    all.FP_TSS.80[k] <- NA # if NA present, FP = NA
    all.ES_TSS.80[k] <- NA # if NA present, ES = NA
    all.mean.FP_TSS.80[k] <- NA # if NA present, mean FP = NA
    all.mean.ES_TSS.80[k] <- NA # if NA present, mean ES = NA
    
  } else {
    FP <- sum(Species_ES_FP_index_table$FP*sp.proba.com) # Compute FP as the weighted sum of the FP index of species, weighted by their probability of presence in the community
    all.FP_TSS.80[k]  <- FP
    mean.FP <- FP/sum(sp.proba.com) # Compute mean FP as the weighted mean of the FP index of species, weighted by their probability of presence in the community
    all.mean.FP_TSS.80[k]  <- mean.FP
    ES <- sum(Species_ES_FP_index_table$ES*sp.proba.com) # Compute ES as the weighted sum of the ES index of species, weighted by their probability of presence in the community
    all.ES_TSS.80[k]  <- ES
    mean.ES <- ES/sum(sp.proba.com) # Compute mean ES as the weighted mean of the ES index of species, weighted by their probability of presence in the community
    all.mean.ES_TSS.80[k]  <- mean.ES
  }
  
  # Show k every 1000 iterations and save a back-up file
  if (k %% 1000 == 0) {
    cat(paste0(Sys.time(), " - ", k," on ",nrow(sp.proba.mat_TSS.80),"\n"))
    save(all.FP_TSS.80, all.mean.FP_TSS.80, all.ES_TSS.80, all.mean.ES_TSS.80, file = "./outputs/Indices_maps/backup.FP.ES_TSS.80.RData")
  }
  if(k == nrow(sp.proba.mat_TSS.80)) {
    cat(paste0(Sys.time(), " - ", k," on ",nrow(sp.proba.mat_TSS.80),"\n"))
    save(all.FP_TSS.80, all.mean.FP_TSS.80, all.ES_TSS.80, all.mean.ES_TSS.80, file = "./outputs/Indices_maps/backup.FP.ES_TSS.80.RData")
  }
}

# Write FP values in a raster
FP.raster_TSS.80 <- continent_mask
FP.raster_TSS.80@data@values <-  all.FP_TSS.80
FP.raster_TSS.80@data@max <- max(FP.raster_TSS.80[], na.rm = T) 
# plot(FP.raster_TSS.80)

# Write FP values in a raster
FP.mean.raster_TSS.80 <- continent_mask
FP.mean.raster_TSS.80@data@values <-  all.mean.FP_TSS.80
FP.mean.raster_TSS.80@data@max <- max(FP.mean.raster_TSS.80[], na.rm = T)
# plot(FP.mean.raster_TSS.80)

# Write ES values in a raster
ES.raster_TSS.80 <- continent_mask
ES.raster_TSS.80@data@values <-  all.ES_TSS.80
ES.raster_TSS.80@data@max <- max(ES.raster_TSS.80[], na.rm = T) 
# plot(ES.raster_TSS.80)

# Write ES values in a raster
ES.mean.raster_TSS.80 <- continent_mask
ES.mean.raster_TSS.80@data@values <-  all.mean.ES_TSS.80
ES.mean.raster_TSS.80@data@max <- max(ES.mean.raster_TSS.80[], na.rm = T)
# plot(ES.mean.raster_TSS.80)

# Save
save(FP.raster_TSS.80, file = "./outputs/Indices_Maps/FP.raster_TSS.80.RData", version = "2")
saveRDS(FP.raster_TSS.80, file = "./outputs/Indices_Maps/FP.raster_TSS.80.rds", version = "2")
save(FP.mean.raster_TSS.80, file = "./outputs/Indices_Maps/FP.mean.raster_TSS.80.RData", version = "2")
saveRDS(FP.mean.raster_TSS.80, file = "./outputs/Indices_Maps/FP.mean.raster_TSS.80.rds", version = "2")
save(ES.raster_TSS.80, file = "./outputs/Indices_Maps/ES.raster_TSS.80.RData", version = "2")
saveRDS(ES.raster_TSS.80, file = "./outputs/Indices_Maps/ES.raster_TSS.80.rds", version = "2")
save(ES.mean.raster_TSS.80, file = "./outputs/Indices_Maps/ES.mean.raster_TSS.80.RData", version = "2")
saveRDS(ES.mean.raster_TSS.80, file = "./outputs/Indices_Maps/ES.mean.raster_TSS.80.rds", version = "2")


### Compute for TSS.95

# Make a loop per community
all.FP_TSS.95 <- all.mean.FP_TSS.95 <- all.ES_TSS.95 <- all.mean.ES_TSS.95 <- rep(NA, nrow(sp.proba.mat_TSS.95))
for (k in 1:nrow(sp.proba.mat_TSS.95)) 
{
  sp.proba.com <- sp.proba.mat_TSS.95[k,] # Extract the proba for the k community
  
  # Compute FP and ES only if no NA is present in the community
  if (any(is.na(sp.proba.com))) 
  {
    all.FP_TSS.95[k] <- NA # if NA present, FP = NA
    all.ES_TSS.95[k] <- NA # if NA present, ES = NA
    all.mean.FP_TSS.95[k] <- NA # if NA present, mean FP = NA
    all.mean.ES_TSS.95[k] <- NA # if NA present, mean ES = NA
    
  } else {
    FP <- sum(Species_ES_FP_index_table$FP*sp.proba.com) # Compute FP as the weighted sum of the FP index of species, weighted by their probability of presence in the community
    all.FP_TSS.95[k]  <- FP
    mean.FP <- FP/sum(sp.proba.com) # Compute mean FP as the weighted mean of the FP index of species, weighted by their probability of presence in the community
    all.mean.FP_TSS.95[k]  <- mean.FP
    ES <- sum(Species_ES_FP_index_table$ES*sp.proba.com) # Compute ES as the weighted sum of the ES index of species, weighted by their probability of presence in the community
    all.ES_TSS.95[k]  <- ES
    mean.ES <- ES/sum(sp.proba.com) # Compute mean ES as the weighted mean of the ES index of species, weighted by their probability of presence in the community
    all.mean.ES_TSS.95[k]  <- mean.ES
  }
  
  # Show k every 1000 iterations and save a back-up file
  if (k %% 1000 == 0) {
    cat(paste0(Sys.time(), " - ", k," on ",nrow(sp.proba.mat_TSS.95),"\n"))
    save(all.FP_TSS.95, all.mean.FP_TSS.95, all.ES_TSS.95, all.mean.ES_TSS.95, file = "./outputs/Indices_maps/backup.FP.ES_TSS.95.RData")
  }
  if(k == nrow(sp.proba.mat_TSS.95)) {
    cat(paste0(Sys.time(), " - ", k," on ",nrow(sp.proba.mat_TSS.95),"\n"))
    save(all.FP_TSS.95, all.mean.FP_TSS.95, all.ES_TSS.95, all.mean.ES_TSS.95, file = "./outputs/Indices_maps/backup.FP.ES_TSS.95.RData")
  }
}

# Write FP values in a raster
FP.raster_TSS.95 <- continent_mask
FP.raster_TSS.95@data@values <-  all.FP_TSS.95
FP.raster_TSS.95@data@max <- max(FP.raster_TSS.95[], na.rm = T) 
# plot(FP.raster_TSS.95)

# Write FP values in a raster
FP.mean.raster_TSS.95 <- continent_mask
FP.mean.raster_TSS.95@data@values <-  all.mean.FP_TSS.95
FP.mean.raster_TSS.95@data@max <- max(FP.mean.raster_TSS.95[], na.rm = T)
# plot(FP.mean.raster_TSS.95)

# Write ES values in a raster
ES.raster_TSS.95 <- continent_mask
ES.raster_TSS.95@data@values <-  all.ES_TSS.95
ES.raster_TSS.95@data@max <- max(ES.raster_TSS.95[], na.rm = T) 
# plot(ES.raster_TSS.95)

# Write ES values in a raster
ES.mean.raster_TSS.95 <- continent_mask
ES.mean.raster_TSS.95@data@values <-  all.mean.ES_TSS.95
ES.mean.raster_TSS.95@data@max <- max(ES.mean.raster_TSS.95[], na.rm = T)
# plot(ES.mean.raster_TSS.95)

# Save
save(FP.raster_TSS.95, file = "./outputs/Indices_Maps/FP.raster_TSS.95.RData", version = "2")
saveRDS(FP.raster_TSS.95, file = "./outputs/Indices_Maps/FP.raster_TSS.95.rds", version = "2")
save(FP.mean.raster_TSS.95, file = "./outputs/Indices_Maps/FP.mean.raster_TSS.95.RData", version = "2")
saveRDS(FP.mean.raster_TSS.95, file = "./outputs/Indices_Maps/FP.mean.raster_TSS.95.rds", version = "2")
save(ES.raster_TSS.95, file = "./outputs/Indices_Maps/ES.raster_TSS.95.RData", version = "2")
saveRDS(ES.raster_TSS.95, file = "./outputs/Indices_Maps/ES.raster_TSS.95.rds", version = "2")
save(ES.mean.raster_TSS.95, file = "./outputs/Indices_Maps/ES.mean.raster_TSS.95.RData", version = "2")
saveRDS(ES.mean.raster_TSS.95, file = "./outputs/Indices_Maps/ES.mean.raster_TSS.95.rds", version = "2")

### 16.4/ Plot Sum of Fair-Proportions Maps ####

### Load directly the final Sum of FP layer
FP.raster_Jaccard.80 <- readRDS(file = "./outputs/Indices_Maps/FP.raster_Jaccard.80.rds")
FP.raster_Jaccard.95 <- readRDS(file = "./outputs/Indices_Maps/FP.raster_Jaccard.95.rds")
FP.raster_TSS.80 <- readRDS(file = "./outputs/Indices_Maps/FP.raster_TSS.80.rds")
FP.raster_TSS.95 <- readRDS(file = "./outputs/Indices_Maps/FP.raster_TSS.95.rds")

### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/FP.raster_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(FP.raster_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Sum of Fair-Proportions\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(FP.raster_Jaccard.80, locs = seq(0, 600, 100), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/FP.raster_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(FP.raster_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Sum of Fair-Proportions\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(FP.raster_Jaccard.95, locs = seq(0, 600, 100), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/FP.raster_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(FP.raster_TSS.80, col = pal_bl_red_Mannion, main = paste0("Sum of Fair-Proportions\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(FP.raster_TSS.80, locs = seq(0, 600, 100), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/FP.raster_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(FP.raster_TSS.95, col = pal_bl_red_Mannion, main = paste0("Sum of Fair-Proportions\nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(FP.raster_TSS.95, locs = seq(0, 600, 100), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(FP.raster_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/FP.raster_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(FP.raster_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Sum of Fair-Proportions\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(FP.raster_Jaccard.80, locs = seq(0, 600, 100), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(FP.raster_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Sum of Fair-Proportions\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(FP.raster_Jaccard.95, locs = seq(0, 600, 100), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(FP.raster_TSS.80, col = pal_bl_red_Mannion, main = paste0("Sum of Fair-Proportions\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(FP.raster_TSS.80, locs = seq(0, 600, 100), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(FP.raster_TSS.95, col = pal_bl_red_Mannion, main = paste0("Sum of Fair-Proportions\nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(FP.raster_TSS.95, locs = seq(0, 600, 100), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

par(mar = internal_margins)

dev.off()


### 16.5/ Contrast Mean Fair-Proportions Maps ####

### Load directly the final Mean FP layers
FP.mean.raster_Jaccard.80 <- readRDS(file = "./outputs/Indices_Maps/FP.mean.raster_Jaccard.80.rds")
FP.mean.raster_Jaccard.95 <- readRDS(file = "./outputs/Indices_Maps/FP.mean.raster_Jaccard.95.rds")
FP.mean.raster_TSS.80 <- readRDS(file = "./outputs/Indices_Maps/FP.mean.raster_TSS.80.rds")
FP.mean.raster_TSS.95 <- readRDS(file = "./outputs/Indices_Maps/FP.mean.raster_TSS.95.rds")

# hist(FP.mean.raster_Jaccard.80)
# 
# # Contrast values by merging low values below 3.5 and above 6.5
# FP.mean.raster_Jaccard.80_contrasted <- FP.mean.raster_Jaccard.80
# FP.mean.raster_Jaccard.80_contrasted[FP.mean.raster_Jaccard.80 < 3.5] <- 3.51
# FP.mean.raster_Jaccard.80_contrasted[FP.mean.raster_Jaccard.80 > 6.5] <- 6.5
# FP.mean.raster_Jaccard.95_contrasted <- FP.mean.raster_Jaccard.95
# FP.mean.raster_Jaccard.95_contrasted[FP.mean.raster_Jaccard.95 < 3.5] <- 3.51
# FP.mean.raster_Jaccard.95_contrasted[FP.mean.raster_Jaccard.95 > 6.5] <- 6.5
# FP.mean.raster_TSS.80_contrasted <- FP.mean.raster_TSS.80
# FP.mean.raster_TSS.80_contrasted[FP.mean.raster_TSS.80 < 3.5] <- 3.51
# FP.mean.raster_TSS.80_contrasted[FP.mean.raster_TSS.80 > 6.5] <- 6.5
# FP.mean.raster_TSS.95_contrasted <- FP.mean.raster_TSS.95
# FP.mean.raster_TSS.95_contrasted[FP.mean.raster_TSS.95 < 3.5] <- 3.51
# FP.mean.raster_TSS.80_contrasted[FP.mean.raster_TSS.80 > 6.5] <- 6.5
# 
# 
# # Add 0 values for empty continental pixels (transformed into min values as 3.5)
# temp <- continent_mask
# temp[!is.na(FP.mean.raster_Jaccard.80_contrasted[])] <- FP.mean.raster_Jaccard.80_contrasted[!is.na(FP.mean.raster_Jaccard.80_contrasted[])]
# FP.mean.raster_Jaccard.80_contrasted <- temp
# FP.mean.raster_Jaccard.80_contrasted[FP.mean.raster_Jaccard.80_contrasted == 0] <- 3.5
# 
# temp <- continent_mask
# temp[!is.na(FP.mean.raster_Jaccard.95_contrasted[])] <- FP.mean.raster_Jaccard.95_contrasted[!is.na(FP.mean.raster_Jaccard.95_contrasted[])]
# FP.mean.raster_Jaccard.95_contrasted <- temp
# FP.mean.raster_Jaccard.95_contrasted[FP.mean.raster_Jaccard.95_contrasted == 0] <- 3.5
# 
# temp <- continent_mask
# temp[!is.na(FP.mean.raster_TSS.80_contrasted[])] <- FP.mean.raster_TSS.80_contrasted[!is.na(FP.mean.raster_TSS.80_contrasted[])]
# FP.mean.raster_TSS.80_contrasted <- temp
# FP.mean.raster_TSS.80_contrasted[FP.mean.raster_TSS.80_contrasted == 0] <- 3.5
# 
# temp <- continent_mask
# temp[!is.na(FP.mean.raster_TSS.95_contrasted[])] <- FP.mean.raster_TSS.95_contrasted[!is.na(FP.mean.raster_TSS.95_contrasted[])]
# FP.mean.raster_TSS.95_contrasted <- temp
# FP.mean.raster_TSS.95_contrasted[FP.mean.raster_TSS.95_contrasted == 0] <- 3.5
# 
# hist(FP.mean.raster_Jaccard.80_contrasted)
# 
# # Save contrasted rasters
# save(FP.mean.raster_Jaccard.80_contrasted, file = "./outputs/Indices_Maps/FP.mean.raster_Jaccard.80_contrasted.RData", version = "2")
# saveRDS(FP.mean.raster_Jaccard.80_contrasted, file = "./outputs/Indices_Maps/FP.mean.raster_Jaccard.80_contrasted.rds", version = "2")
# save(FP.mean.raster_Jaccard.95_contrasted, file = "./outputs/Indices_Maps/FP.mean.raster_Jaccard.95_contrasted.RData", version = "2")
# saveRDS(FP.mean.raster_Jaccard.95_contrasted, file = "./outputs/Indices_Maps/FP.mean.raster_Jaccard.95_contrasted.rds", version = "2")
# save(FP.mean.raster_TSS.80_contrasted, file = "./outputs/Indices_Maps/FP.mean.raster_TSS.80_contrasted.RData", version = "2")
# saveRDS(FP.mean.raster_TSS.80_contrasted, file = "./outputs/Indices_Maps/FP.mean.raster_TSS.80_contrasted.rds", version = "2")
# save(FP.mean.raster_TSS.95_contrasted, file = "./outputs/Indices_Maps/FP.mean.raster_TSS.95_contrasted.RData", version = "2")
# saveRDS(FP.mean.raster_TSS.95_contrasted, file = "./outputs/Indices_Maps/FP.mean.raster_TSS.95_contrasted.rds", version = "2")

### 16.6/ Plot Mean Fair-Proportions Maps ####

# Contrasted plots are not better so we plot the raw ones

# Add 0 values for empty continental pixels (transformed into min values)
temp <- continent_mask
temp[!is.na(continent_mask[])] <- min(FP.mean.raster_Jaccard.80[], na.rm = T)
temp[!is.na(FP.mean.raster_Jaccard.80[])] <- FP.mean.raster_Jaccard.80[!is.na(FP.mean.raster_Jaccard.80[])]
FP.mean.raster_Jaccard.80 <- temp

temp <- continent_mask
temp[!is.na(continent_mask[])] <- min(FP.mean.raster_Jaccard.95[], na.rm = T)
temp[!is.na(FP.mean.raster_Jaccard.95[])] <- FP.mean.raster_Jaccard.95[!is.na(FP.mean.raster_Jaccard.95[])]
FP.mean.raster_Jaccard.95 <- temp

temp <- continent_mask
temp[!is.na(continent_mask[])] <- min(FP.mean.raster_TSS.80[], na.rm = T)
temp[!is.na(FP.mean.raster_TSS.80[])] <- FP.mean.raster_TSS.80[!is.na(FP.mean.raster_TSS.80[])]
FP.mean.raster_TSS.80 <- temp

temp <- continent_mask
temp[!is.na(continent_mask[])] <- min(FP.mean.raster_TSS.95[], na.rm = T)
temp[!is.na(FP.mean.raster_TSS.95[])] <- FP.mean.raster_TSS.95[!is.na(FP.mean.raster_TSS.95[])]
FP.mean.raster_TSS.95 <- temp


### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/FP.mean.raster_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(FP.mean.raster_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean Fair-Proportions\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(FP.mean.raster_Jaccard.80, locs = seq(3, 7, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/FP.mean.raster_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(FP.mean.raster_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Mean Fair-Proportions\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(FP.mean.raster_Jaccard.95, locs = seq(3, 7, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/FP.mean.raster_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(FP.mean.raster_TSS.80, col = pal_bl_red_Mannion, main = paste0("Mean Fair-Proportions\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(FP.mean.raster_TSS.80, locs = seq(3, 7, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/FP.mean.raster_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(FP.mean.raster_TSS.95, col = pal_bl_red_Mannion, main = paste0("Mean Fair-Proportions\nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(FP.mean.raster_TSS.95, locs = seq(3, 7, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

### Tmap version
#
# library(tmap)
# 
# tm_shape(FP.raster_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/FP.mean.raster_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(FP.mean.raster_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mean Fair-Proportions\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(FP.mean.raster_Jaccard.80, locs = seq(3, 7, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(FP.mean.raster_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Mean Fair-Proportions\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(FP.mean.raster_Jaccard.95, locs = seq(3, 7, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(FP.mean.raster_TSS.80, col = pal_bl_red_Mannion, main = paste0("Mean Fair-Proportions\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(FP.mean.raster_TSS.80, locs = seq(3, 7, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(FP.mean.raster_TSS.95, col = pal_bl_red_Mannion, main = paste0("Mean Fair-Proportions\nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(FP.mean.raster_TSS.95, locs = seq(3, 7, 1), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

par(mar = internal_margins)

dev.off()


### 16.7/ Plot Sum of Equal-Splits Maps ####

### Load directly the final Sum of ES layer
ES.raster_Jaccard.80 <- readRDS(file = "./outputs/Indices_Maps/ES.raster_Jaccard.80.rds")
ES.raster_Jaccard.95 <- readRDS(file = "./outputs/Indices_Maps/ES.raster_Jaccard.95.rds")
ES.raster_TSS.80 <- readRDS(file = "./outputs/Indices_Maps/ES.raster_TSS.80.rds")
ES.raster_TSS.95 <- readRDS(file = "./outputs/Indices_Maps/ES.raster_TSS.95.rds")

### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/ES.raster_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ES.raster_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Sum of Equal-Splits\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ES.raster_Jaccard.80, locs = seq(0, 1800, 300), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/ES.raster_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ES.raster_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Sum of Equal-Splits\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ES.raster_Jaccard.95, locs = seq(0, 1800, 300), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/ES.raster_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ES.raster_TSS.80, col = pal_bl_red_Mannion, main = paste0("Sum of Equal-Splits\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ES.raster_TSS.80, locs = seq(0, 1800, 300), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/ES.raster_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ES.raster_TSS.95, col = pal_bl_red_Mannion, main = paste0("Sum of Equal-Splits\nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ES.raster_TSS.95, locs = seq(0, 1800, 300), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(ES.raster_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/ES.raster_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(ES.raster_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Sum of Equal-Splits\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ES.raster_Jaccard.80, locs = seq(0, 1800, 300), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(ES.raster_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Sum of Equal-Splits\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ES.raster_Jaccard.95, locs = seq(0, 1800, 300), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(ES.raster_TSS.80, col = pal_bl_red_Mannion, main = paste0("Sum of Equal-Splits\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ES.raster_TSS.80, locs = seq(0, 1800, 300), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(ES.raster_TSS.95, col = pal_bl_red_Mannion, main = paste0("Sum of Equal-Splits\nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ES.raster_TSS.95, locs = seq(0, 1800, 300), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

par(mar = internal_margins)

dev.off()

### 16.8/ Contrast Mean Equal-Splits Maps ####

### Load directly the final Mean ES layers
ES.mean.raster_Jaccard.80 <- readRDS(file = "./outputs/Indices_Maps/ES.mean.raster_Jaccard.80.rds")
ES.mean.raster_Jaccard.95 <- readRDS(file = "./outputs/Indices_Maps/ES.mean.raster_Jaccard.95.rds")
ES.mean.raster_TSS.80 <- readRDS(file = "./outputs/Indices_Maps/ES.mean.raster_TSS.80.rds")
ES.mean.raster_TSS.95 <- readRDS(file = "./outputs/Indices_Maps/ES.mean.raster_TSS.95.rds")

plot(ES.mean.raster_Jaccard.80)
hist(ES.mean.raster_Jaccard.80)

# Contrast values by merging low values below 14.25 and above 15.25
ES.mean.raster_Jaccard.80_contrasted <- ES.mean.raster_Jaccard.80
ES.mean.raster_Jaccard.80_contrasted[ES.mean.raster_Jaccard.80 < 14.25] <- 14.27
ES.mean.raster_Jaccard.80_contrasted[ES.mean.raster_Jaccard.80 > 15.25] <- 15.25
ES.mean.raster_Jaccard.95_contrasted <- ES.mean.raster_Jaccard.95
ES.mean.raster_Jaccard.95_contrasted[ES.mean.raster_Jaccard.95 < 14.25] <- 14.27
ES.mean.raster_Jaccard.95_contrasted[ES.mean.raster_Jaccard.95 > 15.25] <- 15.25
ES.mean.raster_TSS.80_contrasted <- ES.mean.raster_TSS.80
ES.mean.raster_TSS.80_contrasted[ES.mean.raster_TSS.80 < 14.25] <- 14.27
ES.mean.raster_TSS.80_contrasted[ES.mean.raster_TSS.80 > 15.25] <- 15.25
ES.mean.raster_TSS.95_contrasted <- ES.mean.raster_TSS.95
ES.mean.raster_TSS.95_contrasted[ES.mean.raster_TSS.95 < 14.25] <- 14.27
ES.mean.raster_TSS.95_contrasted[ES.mean.raster_TSS.95 > 15.25] <- 15.25


# Add 0 values for empty continental pixels (transformed into min values as 14.25)
temp <- continent_mask
temp[!is.na(ES.mean.raster_Jaccard.80_contrasted[])] <- ES.mean.raster_Jaccard.80_contrasted[!is.na(ES.mean.raster_Jaccard.80_contrasted[])]
ES.mean.raster_Jaccard.80_contrasted <- temp
ES.mean.raster_Jaccard.80_contrasted[ES.mean.raster_Jaccard.80_contrasted == 0] <- 14.25

temp <- continent_mask
temp[!is.na(ES.mean.raster_Jaccard.95_contrasted[])] <- ES.mean.raster_Jaccard.95_contrasted[!is.na(ES.mean.raster_Jaccard.95_contrasted[])]
ES.mean.raster_Jaccard.95_contrasted <- temp
ES.mean.raster_Jaccard.95_contrasted[ES.mean.raster_Jaccard.95_contrasted == 0] <- 14.25

temp <- continent_mask
temp[!is.na(ES.mean.raster_TSS.80_contrasted[])] <- ES.mean.raster_TSS.80_contrasted[!is.na(ES.mean.raster_TSS.80_contrasted[])]
ES.mean.raster_TSS.80_contrasted <- temp
ES.mean.raster_TSS.80_contrasted[ES.mean.raster_TSS.80_contrasted == 0] <- 14.25

temp <- continent_mask
temp[!is.na(ES.mean.raster_TSS.95_contrasted[])] <- ES.mean.raster_TSS.95_contrasted[!is.na(ES.mean.raster_TSS.95_contrasted[])]
ES.mean.raster_TSS.95_contrasted <- temp
ES.mean.raster_TSS.95_contrasted[ES.mean.raster_TSS.95_contrasted == 0] <- 14.25

hist(ES.mean.raster_Jaccard.80_contrasted)

# Save contrasted rasters
save(ES.mean.raster_Jaccard.80_contrasted, file = "./outputs/Indices_Maps/ES.mean.raster_Jaccard.80_contrasted.RData", version = "2")
saveRDS(ES.mean.raster_Jaccard.80_contrasted, file = "./outputs/Indices_Maps/ES.mean.raster_Jaccard.80_contrasted.rds", version = "2")
save(ES.mean.raster_Jaccard.95_contrasted, file = "./outputs/Indices_Maps/ES.mean.raster_Jaccard.95_contrasted.RData", version = "2")
saveRDS(ES.mean.raster_Jaccard.95_contrasted, file = "./outputs/Indices_Maps/ES.mean.raster_Jaccard.95_contrasted.rds", version = "2")
save(ES.mean.raster_TSS.80_contrasted, file = "./outputs/Indices_Maps/ES.mean.raster_TSS.80_contrasted.RData", version = "2")
saveRDS(ES.mean.raster_TSS.80_contrasted, file = "./outputs/Indices_Maps/ES.mean.raster_TSS.80_contrasted.rds", version = "2")
save(ES.mean.raster_TSS.95_contrasted, file = "./outputs/Indices_Maps/ES.mean.raster_TSS.95_contrasted.RData", version = "2")
saveRDS(ES.mean.raster_TSS.95_contrasted, file = "./outputs/Indices_Maps/ES.mean.raster_TSS.95_contrasted.rds", version = "2")


### 16.9/ Plot Mean Equal-Splits Maps (contrasted) ####

### Load directly the final Mean ES layers
ES.mean.raster_Jaccard.80_contrasted <- readRDS(file = "./outputs/Indices_Maps/ES.mean.raster_Jaccard.80_contrasted.rds")
ES.mean.raster_Jaccard.95_contrasted <- readRDS(file = "./outputs/Indices_Maps/ES.mean.raster_Jaccard.95_contrasted.rds")
ES.mean.raster_TSS.80_contrasted <- readRDS(file = "./outputs/Indices_Maps/ES.mean.raster_TSS.80_contrasted.rds")
ES.mean.raster_TSS.95_contrasted <- readRDS(file = "./outputs/Indices_Maps/ES.mean.raster_TSS.95_contrasted.rds")


### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/ES.mean.raster_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ES.mean.raster_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean Equal-Splits\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ES.mean.raster_Jaccard.80_contrasted, locs = seq(14.3, 15.2, 0.15), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Indices_maps/Jaccard.95/ES.mean.raster_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ES.mean.raster_Jaccard.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean Equal-Splits\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ES.mean.raster_Jaccard.95_contrasted, locs = seq(14.3, 15.2, 0.15), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Indices_maps/TSS.80/ES.mean.raster_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ES.mean.raster_TSS.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean Equal-Splits\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ES.mean.raster_TSS.80_contrasted, locs = seq(14.3, 15.2, 0.15), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

# TSS.95
pdf(file = paste0("./maps/Indices_maps/TSS.95/ES.mean.raster_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ES.mean.raster_TSS.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean Equal-Splits\nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ES.mean.raster_TSS.95_contrasted, locs = seq(14.3, 15.2, 0.15), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()

### Tmap version
#
# library(tmap)
# 
# tm_shape(ES.raster_Jaccard.80) +
#   tm_raster(palette = pal_bl_red_Mannion) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/ES.mean.raster_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(ES.mean.raster_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean Equal-Splits\nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ES.mean.raster_Jaccard.80_contrasted, locs = seq(14.3, 15.2, 0.15), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(ES.mean.raster_Jaccard.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean Equal-Splits\nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ES.mean.raster_Jaccard.95_contrasted, locs = seq(14.3, 15.2, 0.15), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(ES.mean.raster_TSS.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean Equal-Splits\nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ES.mean.raster_TSS.80_contrasted, locs = seq(14.3, 15.2, 0.15), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

image(ES.mean.raster_TSS.95_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean Equal-Splits\nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(ES.mean.raster_TSS.95_contrasted, locs = seq(14.3, 15.2, 0.15), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

par(mar = internal_margins)

dev.off()



