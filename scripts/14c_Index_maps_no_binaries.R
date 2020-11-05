
##### Script 14b: Compute indices maps but only with species with OMU modeled with SDM #####

## Create 9 index maps

# 1/ Species richness
# 2/ Species Shannon's Diversity
# 3/ Species Shannon's Diversity with compatibility with species richness: (a) Raw and (b) with Jost's transformation (exponential)
# 4/ Species Evenness

# 9/ Continuous Species Rarity = Range-size Weighted Species Richness. (a) Full and (b) mean
# 10/ Categorical Species Rarity (25% threshold) : (a) Total and (b) Proportion


## Generate an indice map version for each 4 options between Jaccard/TSS and Buffer.80/95

## Put multiple pdf in main folder and individual maps in subfolder such as /TSS_Buffer.80/


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

# Color palet for plot
cool = rainbow(49, s = 1, v = 1, start=rgb2hsv(col2rgb('yellow'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm = rainbow(50, s = 1, v= 1, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
pal_bl_red  = c(rev(cool), rev(warm))
pal_bl_red <- c(gplots::col2hex("grey93"), pal_bl_red)

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


###### 0/ Generate the list of species with no Binary maps (all OMU have a niche model) to keep for these maps

# Extract only the 236 species with no Binary maps (all OMU have a niche model) from list.sp
list.sp_no_binary <- list.sp[!list.sp$one.Binaries, ]
nrow(list.sp_no_binary) # 236 species out of 388 (60.8%) have no binaries, only modeled OMU

# Extract their names
no_binary_sp <- list.sp_no_binary$Sp_full


###### 1/ Species richness ######

### 1.1/ Species Stack generation ####

### Load directly the complete stack of species proba
All_sp_proba_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))
All_sp_proba_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.95.rds"))
All_sp_proba_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80.rds"))
All_sp_proba_stack_TSS.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.95.rds"))

names(All_sp_proba_stack_Jaccard.80)

# Keep only species with no binary maps
All_sp_proba_stack_Jaccard.80_no_binary <- All_sp_proba_stack_Jaccard.80[[no_binary_sp]]
All_sp_proba_stack_Jaccard.95_no_binary <- All_sp_proba_stack_Jaccard.95[[no_binary_sp]]
All_sp_proba_stack_TSS.80_no_binary <- All_sp_proba_stack_TSS.80[[no_binary_sp]]
All_sp_proba_stack_TSS.95_no_binary <- All_sp_proba_stack_TSS.95[[no_binary_sp]]

# Save new stacks
save(All_sp_proba_stack_Jaccard.80_no_binary, file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_Jaccard.80_no_binary.RData"), version = "2")
saveRDS(All_sp_proba_stack_Jaccard.80_no_binary, file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_Jaccard.80_no_binary.rds"), version = "2")
save(All_sp_proba_stack_Jaccard.95_no_binary, file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_Jaccard.95_no_binary.RData"), version = "2")
saveRDS(All_sp_proba_stack_Jaccard.95_no_binary, file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_Jaccard.95_no_binary.rds"), version = "2")
save(All_sp_proba_stack_Jaccard.80_no_binary, file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_TSS.80_no_binary.RData"), version = "2")
saveRDS(All_sp_proba_stack_Jaccard.80_no_binary, file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_TSS.80_no_binary.rds"), version = "2")
save(All_sp_proba_stack_TSS.95_no_binary, file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_TSS.95_no_binary.RData"), version = "2")
saveRDS(All_sp_proba_stack_TSS.95_no_binary, file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_TSS.95_no_binary.rds"), version = "2")


### 1.2/ Index computation ####
tot.sp.richness_Jaccard.80_no_binary <- readAll(calc(All_sp_proba_stack_Jaccard.80_no_binary, fun = sum))
tot.sp.richness_Jaccard.95_no_binary <- readAll(calc(All_sp_proba_stack_Jaccard.95_no_binary, fun = sum))
tot.sp.richness_TSS.80_no_binary <- readAll(calc(All_sp_proba_stack_TSS.80_no_binary, fun = sum))
tot.sp.richness_TSS.95_no_binary <- readAll(calc(All_sp_proba_stack_TSS.95_no_binary, fun = sum))


# Save
save(tot.sp.richness_Jaccard.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_Jaccard.80_no_binary.RData"), version = "2")
saveRDS(tot.sp.richness_Jaccard.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_Jaccard.80_no_binary.rds"), version = "2")
save(tot.sp.richness_Jaccard.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_Jaccard.95_no_binary.RData"), version = "2")
saveRDS(tot.sp.richness_Jaccard.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_Jaccard.95_no_binary.rds"), version = "2")
save(tot.sp.richness_TSS.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_TSS.80_no_binary.RData"), version = "2")
saveRDS(tot.sp.richness_TSS.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_TSS.80_no_binary.rds"), version = "2")
save(tot.sp.richness_TSS.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_TSS.95_no_binary.RData"), version = "2")
saveRDS(tot.sp.richness_TSS.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_TSS.95_no_binary.rds"), version = "2")


### Load directly the final Species richness layer
tot.sp.richness_Jaccard.80_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_Jaccard.80_no_binary.rds"))
tot.sp.richness_Jaccard.95_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_Jaccard.95_no_binary.rds"))
tot.sp.richness_TSS.80_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_TSS.80_no_binary.rds"))
tot.sp.richness_TSS.95_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_TSS.95_no_binary.rds"))

### 1.3/ Plot Species richness ####

### Individual plots

# Jaccard.80_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.80_no_binary/tot.sp.richness_Jaccard.80_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.sp.richness_Jaccard.80_no_binary, col = pal_bl_red, main = paste0("Species richness \nJaccard.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_Jaccard.80_no_binary, locs = seq(0, 80, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# Jaccard.95_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.95_no_binary/tot.sp.richness_Jaccard.95_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.sp.richness_Jaccard.95_no_binary, col = pal_bl_red, main = paste0("Species richness \nJaccard.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_Jaccard.95_no_binary, locs = seq(0, 80, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# TSS.80_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.80_no_binary/tot.sp.richness_TSS.80_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.sp.richness_TSS.80_no_binary, col = pal_bl_red, main = paste0("Species richness \nTSS.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_TSS.80_no_binary, locs = seq(0, 80, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# TSS.95_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.95_no_binary/tot.sp.richness_TSS.95_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.sp.richness_TSS.95_no_binary, col = pal_bl_red, main = paste0("Species richness \nTSS.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_TSS.95_no_binary, locs = seq(0, 80, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(tot.sp.richness_Jaccard.80_no_binary) +
#   tm_raster(palette = pal_bl_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/no_binary/tot.sp.richness_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(tot.sp.richness_Jaccard.80_no_binary, col = pal_bl_red, main = paste0("Species richness \nJaccard.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_Jaccard.80_no_binary, locs = seq(0, 80, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1.2, label = "Species")

image(tot.sp.richness_Jaccard.95_no_binary, col = pal_bl_red, main = paste0("Species richness \nJaccard.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_Jaccard.95_no_binary, locs = seq(0, 80, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1.2, label = "Species")

image(tot.sp.richness_TSS.80_no_binary, col = pal_bl_red, main = paste0("Species richness \nTSS.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_TSS.80_no_binary, locs = seq(0, 80, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1.2, label = "Species")

image(tot.sp.richness_TSS.95_no_binary, col = pal_bl_red, main = paste0("Species richness \nTSS.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_TSS.95_no_binary, locs = seq(0, 80, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1.2, label = "Species")

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
All_sp_proba_stack_Jaccard.80_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_Jaccard.80_no_binary.rds"))
All_sp_proba_stack_Jaccard.95_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_Jaccard.95_no_binary.rds"))
All_sp_proba_stack_TSS.80_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_TSS.80_no_binary.rds"))
All_sp_proba_stack_TSS.95_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_TSS.95_no_binary.rds"))

### 2.1/ Index computation ####
sp.diversity_Jaccard.80_no_binary <- calc(All_sp_proba_stack_Jaccard.80_no_binary, fun = shannon)*1
sp.diversity_Jaccard.95_no_binary <- calc(All_sp_proba_stack_Jaccard.95_no_binary, fun = shannon)*1
sp.diversity_TSS.80_no_binary <- calc(All_sp_proba_stack_TSS.80_no_binary, fun = shannon)*1
sp.diversity_TSS.95_no_binary <- calc(All_sp_proba_stack_TSS.95_no_binary, fun = shannon)*1

# Save
save(sp.diversity_Jaccard.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity_Jaccard.80_no_binary.RData"), version = "2")
saveRDS(sp.diversity_Jaccard.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity_Jaccard.80_no_binary.rds"), version = "2")
save(sp.diversity_Jaccard.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity_Jaccard.95_no_binary.RData"), version = "2")
saveRDS(sp.diversity_Jaccard.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity_Jaccard.95_no_binary.rds"), version = "2")
save(sp.diversity_TSS.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity_TSS.80_no_binary.RData"), version = "2")
saveRDS(sp.diversity_TSS.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity_TSS.80_no_binary.rds"), version = "2")
save(sp.diversity_TSS.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity_TSS.95_no_binary.RData"), version = "2")
saveRDS(sp.diversity_TSS.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity_TSS.95_no_binary.rds"), version = "2")


### Load directly the final Species diversity layers
sp.diversity_Jaccard.80_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/sp.diversity_Jaccard.80_no_binary.rds"))
sp.diversity_Jaccard.95_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/sp.diversity_Jaccard.95_no_binary.rds"))
sp.diversity_TSS.80_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/sp.diversity_TSS.80_no_binary.rds"))
sp.diversity_TSS.95_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/sp.diversity_TSS.95_no_binary.rds"))

### 2.2/ Plot Species Diversity raw version ####

### Individual plots

# Jaccard.80_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.80_no_binary/sp_diversity_Jaccard.80_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity_Jaccard.80_no_binary, col = pal_bl_red, main = paste0("Shannon's species diversity \nJaccard.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_Jaccard.80_no_binary, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# Jaccard.95_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.95_no_binary/sp_diversity_Jaccard.95_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity_Jaccard.95_no_binary, col = pal_bl_red, main = paste0("Shannon's species diversity \nJaccard.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_Jaccard.95_no_binary, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# TSS.80_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.80_no_binary/sp_diversity_TSS.80_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity_TSS.80_no_binary, col = pal_bl_red, main = paste0("Shannon's species diversity \nTSS.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_TSS.80_no_binary, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# TSS.95_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.95_no_binary/sp_diversity_TSS.95_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity_TSS.95_no_binary, col = pal_bl_red, main = paste0("Shannon's species diversity \nTSS.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_TSS.95_no_binary, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# ## Tmap version
# 
# library(tmap)
# 
# tm_shape(sp.diversity_Jaccard.80_no_binary) +
#   tm_raster(palette = pal_bl_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")


### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/no_binary/sp.diversity_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.diversity_Jaccard.80_no_binary, col = pal_bl_red, main = paste0("Shannon's species diversity \nJaccard.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_Jaccard.80_no_binary, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(sp.diversity_Jaccard.95_no_binary, col = pal_bl_red, main = paste0("Shannon's species diversity \nJaccard.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_Jaccard.95_no_binary, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(sp.diversity_TSS.80_no_binary, col = pal_bl_red, main = paste0("Shannon's species diversity \nTSS.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_TSS.80_no_binary, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(sp.diversity_TSS.95_no_binary, col = pal_bl_red, main = paste0("Shannon's species diversity \nTSS.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_TSS.95_no_binary, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
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
All_sp_proba_stack_Jaccard.80_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_Jaccard.80_no_binary.rds"))
All_sp_proba_stack_Jaccard.95_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_Jaccard.95_no_binary.rds"))
All_sp_proba_stack_TSS.80_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_TSS.80_no_binary.rds"))
All_sp_proba_stack_TSS.95_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_TSS.95_no_binary.rds"))


### 3.1/ Index computation ####
sp.diversity.compatible_Jaccard.80_no_binary <- calc(All_sp_proba_stack_Jaccard.80_no_binary, fun = shannon_compatible)*1 # Raw version
sp.diversity.compatible_Jost_Jaccard.80_no_binary <- exp(sp.diversity.compatible_Jaccard.80_no_binary) - 1 # Jost's version

sp.diversity.compatible_Jaccard.95_no_binary <- calc(All_sp_proba_stack_Jaccard.95_no_binary, fun = shannon_compatible)*1 # Raw version
sp.diversity.compatible_Jost_Jaccard.95_no_binary <- exp(sp.diversity.compatible_Jaccard.95_no_binary) - 1 # Jost's version

sp.diversity.compatible_TSS.80_no_binary <- calc(All_sp_proba_stack_TSS.80_no_binary, fun = shannon_compatible)*1 # Raw version
sp.diversity.compatible_Jost_TSS.80_no_binary <- exp(sp.diversity.compatible_TSS.80_no_binary) - 1 # Jost's version

sp.diversity.compatible_TSS.95_no_binary <- calc(All_sp_proba_stack_TSS.95_no_binary, fun = shannon_compatible)*1 # Raw version
sp.diversity.compatible_Jost_TSS.95_no_binary <- exp(sp.diversity.compatible_TSS.95_no_binary) - 1 # Jost's version

# Save
save(sp.diversity.compatible_Jaccard.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_Jaccard.80_no_binary.RData"), version = "2")
saveRDS(sp.diversity.compatible_Jaccard.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_Jaccard.80_no_binary.rds"), version = "2")
save(sp.diversity.compatible_Jaccard.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_Jaccard.95_no_binary.RData"), version = "2")
saveRDS(sp.diversity.compatible_Jaccard.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_Jaccard.95_no_binary.rds"), version = "2")
save(sp.diversity.compatible_TSS.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_TSS.80_no_binary.RData"), version = "2")
saveRDS(sp.diversity.compatible_TSS.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_TSS.80_no_binary.rds"), version = "2")
save(sp.diversity.compatible_TSS.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_TSS.95_no_binary.RData"), version = "2")
saveRDS(sp.diversity.compatible_TSS.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_TSS.95_no_binary.rds"), version = "2")

save(sp.diversity.compatible_Jost_Jaccard.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_Jost_Jaccard.80_no_binary.RData"), version = "2")
saveRDS(sp.diversity.compatible_Jost_Jaccard.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_Jost_Jaccard.80_no_binary.rds"), version = "2")
save(sp.diversity.compatible_Jost_Jaccard.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_Jost_Jaccard.95_no_binary.RData"), version = "2")
saveRDS(sp.diversity.compatible_Jost_Jaccard.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_Jost_Jaccard.95_no_binary.rds"), version = "2")
save(sp.diversity.compatible_Jost_TSS.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_Jost_TSS.80_no_binary.RData"), version = "2")
saveRDS(sp.diversity.compatible_Jost_TSS.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_Jost_TSS.80_no_binary.rds"), version = "2")
save(sp.diversity.compatible_Jost_TSS.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_Jost_TSS.95_no_binary.RData"), version = "2")
saveRDS(sp.diversity.compatible_Jost_TSS.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_Jost_TSS.95_no_binary.rds"), version = "2")


### Load directly the final Species corrected diversity layers
sp.diversity.compatible_Jaccard.80_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_Jaccard.80_no_binary.rds"))
sp.diversity.compatible_Jaccard.95_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_Jaccard.95_no_binary.rds"))
sp.diversity.compatible_TSS.80_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_TSS.80_no_binary.rds"))
sp.diversity.compatible_TSS.95_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_TSS.95_no_binary.rds"))

sp.diversity.compatible_Jost_Jaccard.80_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_Jost_Jaccard.80_no_binary.rds"))
sp.diversity.compatible_Jost_Jaccard.95_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_Jost_Jaccard.95_no_binary.rds"))
sp.diversity.compatible_Jost_TSS.80_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_Jost_TSS.80_no_binary.rds"))
sp.diversity.compatible_Jost_TSS.95_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/sp.diversity.compatible_Jost_TSS.95_no_binary.rds"))


### 3.2/ Plot Species Diversity Corrected version ####

### Individual plots

# Jaccard.80_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.80_no_binary/sp_diversity_compatible_Jaccard.80_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_Jaccard.80_no_binary, col = pal_bl_red, main = paste0("Shannon's species diversity \nJaccard.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jaccard.80_no_binary, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# Jaccard.95_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.95_no_binary/sp_diversity_compatible_Jaccard.95_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_Jaccard.95_no_binary, col = pal_bl_red, main = paste0("Shannon's species diversity \nJaccard.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jaccard.95_no_binary, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# TSS.80_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.80_no_binary/sp_diversity_compatible_TSS.80_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_TSS.80_no_binary, col = pal_bl_red, main = paste0("Shannon's species diversity \nTSS.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_TSS.80_no_binary, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# TSS.95_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.95_no_binary/sp_diversity_compatible_TSS.95_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_TSS.95_no_binary, col = pal_bl_red, main = paste0("Shannon's species diversity \nTSS.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_TSS.95_no_binary, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# ## Tmap version
# 
# library(tmap)
# 
# tm_shape(sp.diversity.compatible_Jaccard.80_no_binary) +
#   tm_raster(palette = pal_bl_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")


### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/no_binary/sp.diversity.compatible_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.diversity.compatible_Jaccard.80_no_binary, col = pal_bl_red, main = paste0("Shannon's species diversity \nJaccard.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jaccard.80_no_binary, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(sp.diversity.compatible_Jaccard.95_no_binary, col = pal_bl_red, main = paste0("Shannon's species diversity \nJaccard.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jaccard.95_no_binary, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(sp.diversity.compatible_TSS.80_no_binary, col = pal_bl_red, main = paste0("Shannon's species diversity \nTSS.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_TSS.80_no_binary, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(sp.diversity.compatible_TSS.95_no_binary, col = pal_bl_red, main = paste0("Shannon's species diversity \nTSS.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_TSS.95_no_binary, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

par(mar = internal_margins)

dev.off()

### 3.3/ Plot Species Diversity Jost version ####

### Individual plots

# Jaccard.80_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.80_no_binary/sp_diversity_compatible_Jost_Jaccard.80_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
image(sp.diversity.compatible_Jost_Jaccard.80_no_binary, col = pal_bl_red, main = paste0("Species diversity (Jost's effective species richness) \nJaccard.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_Jaccard.80_no_binary, locs = seq(0, 80, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# Jaccard.95_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.95_no_binary/sp_diversity_compatible_Jost_Jaccard.95_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_Jost_Jaccard.95_no_binary, col = pal_bl_red, main = paste0("Species diversity (Jost's effective species richness) \nJaccard.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_Jaccard.95_no_binary, locs = seq(0, 80, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# TSS.80_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.80_no_binary/sp_diversity_compatible_Jost_TSS.80_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_Jost_TSS.80_no_binary, col = pal_bl_red, main = paste0("Species diversity (Jost's effective species richness) \nTSS.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_TSS.80_no_binary, locs = seq(0, 80, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# TSS.95_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.95_no_binary/sp_diversity_compatible_Jost_TSS.95_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_Jost_TSS.95_no_binary, col = pal_bl_red, main = paste0("Species diversity (Jost's effective species richness) \nTSS.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_TSS.95_no_binary, locs = seq(0, 80, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# ## Tmap version
# 
# library(tmap)
# 
# tm_shape(sp.diversity.compatible_Jost_Jaccard.80_no_binary) +
#   tm_raster(palette = pal_bl_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")


### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/no_binary/sp.diversity.compatible_Jost_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.diversity.compatible_Jost_Jaccard.80_no_binary, col = pal_bl_red, main = paste0("Species diversity (Jost's effective species richness) \nJaccard.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_Jaccard.80_no_binary, locs = seq(0, 80, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1, label = "Species")

image(sp.diversity.compatible_Jost_Jaccard.95_no_binary, col = pal_bl_red, main = paste0("Species diversity (Jost's effective species richness) \nJaccard.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_Jaccard.95_no_binary, locs = seq(0, 80, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1, label = "Species")

image(sp.diversity.compatible_Jost_TSS.80_no_binary, col = pal_bl_red, main = paste0("Species diversity (Jost's effective species richness) \nTSS.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_TSS.80_no_binary, locs = seq(0, 80, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1, label = "Species")

image(sp.diversity.compatible_Jost_TSS.95_no_binary, col = pal_bl_red, main = paste0("Species diversity (Jost's effective species richness) \nTSS.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_TSS.95_no_binary, locs = seq(0, 80, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1, label = "Species")

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
All_sp_proba_stack_Jaccard.80_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_Jaccard.80_no_binary.rds"))
All_sp_proba_stack_Jaccard.95_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_Jaccard.95_no_binary.rds"))
All_sp_proba_stack_TSS.80_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_TSS.80_no_binary.rds"))
All_sp_proba_stack_TSS.95_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_TSS.95_no_binary.rds"))

### 4.1/ Index computation ####
sp.evenness_Jaccard.80_no_binary <- calc(All_sp_proba_stack_Jaccard.80_no_binary, fun = evenness)*1
sp.evenness_Jaccard.95_no_binary <- calc(All_sp_proba_stack_Jaccard.95_no_binary, fun = evenness)*1
sp.evenness_TSS.80_no_binary <- calc(All_sp_proba_stack_TSS.80_no_binary, fun = evenness)*1
sp.evenness_TSS.95_no_binary <- calc(All_sp_proba_stack_TSS.95_no_binary, fun = evenness)*1

# Save
save(sp.evenness_Jaccard.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.evenness_Jaccard.80_no_binary.RData"), version = "2")
saveRDS(sp.evenness_Jaccard.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.evenness_Jaccard.80_no_binary.rds"), version = "2")
save(sp.evenness_Jaccard.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.evenness_Jaccard.95_no_binary.RData"), version = "2")
saveRDS(sp.evenness_Jaccard.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.evenness_Jaccard.95_no_binary.rds"), version = "2")
save(sp.evenness_TSS.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.evenness_TSS.80_no_binary.RData"), version = "2")
saveRDS(sp.evenness_TSS.80_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.evenness_TSS.80_no_binary.rds"), version = "2")
save(sp.evenness_TSS.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.evenness_TSS.95_no_binary.RData"), version = "2")
saveRDS(sp.evenness_TSS.95_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.evenness_TSS.95_no_binary.rds"), version = "2")

### 4.2/ Contrast maps ####

# Load directly the final Species evenness layer
sp.evenness_Jaccard.80_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/sp.evenness_Jaccard.80_no_binary.rds"))
sp.evenness_Jaccard.95_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/sp.evenness_Jaccard.95_no_binary.rds"))
sp.evenness_TSS.80_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/sp.evenness_TSS.80_no_binary.rds"))
sp.evenness_TSS.95_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/sp.evenness_TSS.95_no_binary.rds"))

hist(sp.evenness_Jaccard.80_no_binary)

# Contrast values by merging low values below 0.25
sp.evenness_Jaccard.80_no_binary_contrasted <- sp.evenness_Jaccard.80_no_binary
sp.evenness_Jaccard.80_no_binary_contrasted[sp.evenness_Jaccard.80_no_binary < 0.25] <- 0.255
sp.evenness_Jaccard.95_no_binary_contrasted <- sp.evenness_Jaccard.95_no_binary
sp.evenness_Jaccard.95_no_binary_contrasted[sp.evenness_Jaccard.95_no_binary < 0.25] <- 0.255
sp.evenness_TSS.80_no_binary_contrasted <- sp.evenness_TSS.80_no_binary
sp.evenness_TSS.80_no_binary_contrasted[sp.evenness_TSS.80_no_binary < 0.25] <- 0.255
sp.evenness_TSS.95_no_binary_contrasted <- sp.evenness_TSS.95_no_binary
sp.evenness_TSS.95_no_binary_contrasted[sp.evenness_TSS.95_no_binary < 0.25] <- 0.255

# Add 0 values for empty continental pixels (transformed into min values as 0.25)
temp <- continent_mask
temp[!is.na(sp.evenness_Jaccard.80_no_binary_contrasted[])] <- sp.evenness_Jaccard.80_no_binary_contrasted[!is.na(sp.evenness_Jaccard.80_no_binary_contrasted[])]
sp.evenness_Jaccard.80_no_binary_contrasted <- temp
sp.evenness_Jaccard.80_no_binary_contrasted[sp.evenness_Jaccard.80_no_binary_contrasted == 0] <- 0.25

temp <- continent_mask
temp[!is.na(sp.evenness_Jaccard.95_no_binary_contrasted[])] <- sp.evenness_Jaccard.95_no_binary_contrasted[!is.na(sp.evenness_Jaccard.95_no_binary_contrasted[])]
sp.evenness_Jaccard.95_no_binary_contrasted <- temp
sp.evenness_Jaccard.95_no_binary_contrasted[sp.evenness_Jaccard.95_no_binary_contrasted == 0] <- 0.25

temp <- continent_mask
temp[!is.na(sp.evenness_TSS.80_no_binary_contrasted[])] <- sp.evenness_TSS.80_no_binary_contrasted[!is.na(sp.evenness_TSS.80_no_binary_contrasted[])]
sp.evenness_TSS.80_no_binary_contrasted <- temp
sp.evenness_TSS.80_no_binary_contrasted[sp.evenness_TSS.80_no_binary_contrasted == 0] <- 0.25

temp <- continent_mask
temp[!is.na(sp.evenness_TSS.95_no_binary_contrasted[])] <- sp.evenness_TSS.95_no_binary_contrasted[!is.na(sp.evenness_TSS.95_no_binary_contrasted[])]
sp.evenness_TSS.95_no_binary_contrasted <- temp
sp.evenness_TSS.95_no_binary_contrasted[sp.evenness_TSS.95_no_binary_contrasted == 0] <- 0.25

# hist(sp.evenness_Jaccard.80_no_binary_contrasted)

### 4.3/ Plot with contrasted scale merging low values (0.25 = 0-0.25) ####

### Individual plots for contrasted scale

# Jaccard.80_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.80_no_binary/sp.evenness_Jaccard.80_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.evenness_Jaccard.80_no_binary_contrasted, col = pal_bl_red, main = paste0("Species evenness \nJaccard.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_Jaccard.80_no_binary_contrasted, locs = seq(0.25, 0.4, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")
par(mar = internal_margins)
dev.off()

# Jaccard.95_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.80_no_binary/sp.evenness_Jaccard.95_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.evenness_Jaccard.95_no_binary_contrasted, col = pal_bl_red, main = paste0("Species evenness \nJaccard.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="         Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_Jaccard.95_no_binary_contrasted, locs = seq(0.25, 0.40, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")
par(mar = internal_margins)
dev.off()

# TSS.80_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.80_no_binary/sp.evenness_TSS.80_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.evenness_TSS.80_no_binary_contrasted, col = pal_bl_red, main = paste0("Species evenness \nTSS.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="         Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_TSS.80_no_binary_contrasted, locs = seq(0.25, 0.4, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")
par(mar = internal_margins)
dev.off()

# TSS.95_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.80_no_binary/sp.evenness_TSS.95_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.evenness_TSS.95_no_binary_contrasted, col = pal_bl_red, main = paste0("Species evenness \nTSS.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_TSS.95_no_binary_contrasted, locs = seq(0.25, 0.40, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(sp.evenness_Jaccard.80_no_binary) +
#   tm_raster(palette = pal_bl_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/no_binary/sp.evenness_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.evenness_Jaccard.80_no_binary_contrasted, col = pal_bl_red, main = paste0("Species evenness \nJaccard.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_Jaccard.80_no_binary_contrasted, locs = seq(0.25, 0.4, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")

image(sp.evenness_Jaccard.95_no_binary_contrasted, col = pal_bl_red, main = paste0("Species evenness \nJaccard.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="         Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_Jaccard.95_no_binary_contrasted, locs = seq(0.25, 0.40, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")

image(sp.evenness_TSS.80_no_binary_contrasted, col = pal_bl_red, main = paste0("Species evenness \nTSS.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="         Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_TSS.80_no_binary_contrasted, locs = seq(0.25, 0.4, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")

image(sp.evenness_TSS.95_no_binary_contrasted, col = pal_bl_red, main = paste0("Species evenness \nTSS.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_TSS.95_no_binary_contrasted, locs = seq(0.25, 0.40, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")

par(mar = internal_margins)

dev.off()


# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])


###### 5/ to 8/ = Mimicry ring indices and maps. 

# Not so usefull and time consuming to rebuild because need to merge Mimicry rings without species not in the phylogeny. 
# Better to just focus on species-level indices and show there is no difference

###### 9/ Continuous Species Rarity = Range-size Weighted Species Richness ######

### Load directly the stack of sp probas use to compute rarity based on geographic range
All_sp_proba_stack_Jaccard.80_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_Jaccard.80_no_binary.rds"))
All_sp_proba_stack_Jaccard.95_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_Jaccard.95_no_binary.rds"))
All_sp_proba_stack_TSS.80_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_TSS.80_no_binary.rds"))
All_sp_proba_stack_TSS.95_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_TSS.95_no_binary.rds"))

# Load maps of Species richness to use to compute mean rarity index among local species
sp.richness_Jaccard.80_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_Jaccard.80_no_binary.rds"))
sp.richness_Jaccard.95_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_Jaccard.95_no_binary.rds"))
sp.richness_TSS.80_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_TSS.80_no_binary.rds"))
sp.richness_TSS.95_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_TSS.95_no_binary.rds"))

### 9.1/ Index computation ####

## For Jaccard.80_no_binary

# Compute weights based on geographic ranges
sp.range_size_Jaccard.80_no_binary <- NA
for (i in 1:nlayers(All_sp_proba_stack_Jaccard.80_no_binary)) {
  # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
  # Compute the estimated number of pixels occupied by the species as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
  sp.range_size_Jaccard.80_no_binary[i] <- sum(All_sp_proba_stack_Jaccard.80_no_binary[[i]]@data@values, na.rm = T)*774.51/1000
}
sp.rarity_indices_Jaccard.80_no_binary <- 1-(sp.range_size_Jaccard.80_no_binary/max(sp.range_size_Jaccard.80_no_binary)) # Weights by max species extent = rarity indices

# Apply weights
sp.weighted.stack_Jaccard.80_no_binary <- All_sp_proba_stack_Jaccard.80_no_binary # Generate the stack to fill
for (i in 1:nlayers(sp.weighted.stack_Jaccard.80_no_binary)){
  # Multiply species probability of presence by ring rarity indices
  sp.weighted.stack_Jaccard.80_no_binary@layers[[i]]@data@values <- All_sp_proba_stack_Jaccard.80_no_binary@layers[[i]]@data@values * sp.rarity_indices_Jaccard.80_no_binary[i]
}

# Sum = Richness weighted by rarity based on Range-size
sp.rarity_Jaccard.80_no_binary <- readAll(calc(sp.weighted.stack_Jaccard.80_no_binary, fun = sum)) 
# Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each species)
sp.mean.rarity_Jaccard.80_no_binary <- sp.rarity_Jaccard.80_no_binary/sp.richness_Jaccard.80_no_binary

# Add 0 values for empty continental pixels
temp <- continent_mask
temp[!is.na(sp.mean.rarity_Jaccard.80_no_binary[])] <- sp.mean.rarity_Jaccard.80_no_binary[!is.na(sp.mean.rarity_Jaccard.80_no_binary[])]
sp.mean.rarity_Jaccard.80_no_binary <- temp

# Save species ranges and species continous rarity indices
save(sp.range_size_Jaccard.80_no_binary, file = "./outputs/Indices_maps/no_binary/sp.range_size_Jaccard.80_no_binary.RData", version = "2")
saveRDS(sp.range_size_Jaccard.80_no_binary, file = "./outputs/Indices_maps/no_binary/sp.range_size_Jaccard.80_no_binary.rds", version = "2")
save(sp.rarity_indices_Jaccard.80_no_binary, file = "./outputs/Indices_maps/no_binary/sp.rarity_indices_Jaccard.80_no_binary.RData", version = "2")
saveRDS(sp.rarity_indices_Jaccard.80_no_binary, file = "./outputs/Indices_maps/no_binary/sp.rarity_indices_Jaccard.80_no_binary.rds", version = "2")

# Save full and mean rarity maps
save(sp.rarity_Jaccard.80_no_binary, file = "./outputs/Indices_maps/no_binary/sp.rarity_Jaccard.80_no_binary.RData", version = "2")
saveRDS(sp.rarity_Jaccard.80_no_binary, file = "./outputs/Indices_maps/no_binary/sp.rarity_Jaccard.80_no_binary.rds", version = "2")
save(sp.mean.rarity_Jaccard.80_no_binary, file = "./outputs/Indices_maps/no_binary/sp.mean.rarity_Jaccard.80_no_binary.RData", version = "2")
saveRDS(sp.mean.rarity_Jaccard.80_no_binary, file = "./outputs/Indices_maps/no_binary/sp.mean.rarity_Jaccard.80_no_binary.rds", version = "2")

# Plot distribution of ranges and rarity indices
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.80_no_binary/sp_range_distri_Jaccard.80_no_binary.pdf"), height = 5.3, width = 6.5)
hist(sp.range_size_Jaccard.80_no_binary, main = "Distribution of species geographic ranges\nJaccard.80_no_binary", breaks = 20, xlab = "Geographic ranges in [10^3] km²")
range_threshold_Jaccard.80_no_binary <- round(quantile(sp.range_size_Jaccard.80_no_binary, p=0.25), 1) # Range threshold value for qualitative rarity
abline(v = range_threshold_Jaccard.80_no_binary, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",range_threshold_Jaccard.80_no_binary," km²"), col = "red", lty = 2, lwd = 2, x = "topright", cex = 1, bty ="n")
dev.off()

pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.80_no_binary/sp_rarity_indices_distri_Jaccard.80_no_binary.pdf"), height = 5.3, width = 6.5)
hist(sp.rarity_indices_Jaccard.80_no_binary, main = "Distribution of species rarity indices\nJaccard.80_no_binary", breaks = 20, xlab = "Rarity indices")
indice_threshold_Jaccard.80_no_binary <- round(quantile(sp.rarity_indices_Jaccard.80_no_binary, p=0.75), 3) # Threshold value for qualitative rarity
abline(v = indice_threshold_Jaccard.80_no_binary, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",indice_threshold_Jaccard.80_no_binary), col = "red", lty = 2, lwd = 2, x = "topleft", cex = 1, bty ="n")
dev.off()


## For Jaccard.95_no_binary

# Compute weights based on geographic ranges
sp.range_size_Jaccard.95_no_binary <- NA
for (i in 1:nlayers(All_sp_proba_stack_Jaccard.95_no_binary)) {
  # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
  # Compute the estimated number of pixels occupied by the species as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
  sp.range_size_Jaccard.95_no_binary[i] <- sum(All_sp_proba_stack_Jaccard.95_no_binary[[i]]@data@values, na.rm = T)*774.51/1000
}
sp.rarity_indices_Jaccard.95_no_binary <- 1-(sp.range_size_Jaccard.95_no_binary/max(sp.range_size_Jaccard.95_no_binary)) # Weights by max species extent = rarity indices

# Apply weights
sp.weighted.stack_Jaccard.95_no_binary <- All_sp_proba_stack_Jaccard.95_no_binary # Generate the stack to fill
for (i in 1:nlayers(sp.weighted.stack_Jaccard.95_no_binary)){
  # Multiply species probability of presence by species rarity indices
  sp.weighted.stack_Jaccard.95_no_binary@layers[[i]]@data@values <- All_sp_proba_stack_Jaccard.95_no_binary@layers[[i]]@data@values * sp.rarity_indices_Jaccard.95_no_binary[i]
}

# plot(sp.weighted.stack_Jaccard.95_no_binary)

# Sum = Richness weighted by rarity based on Range-size
sp.rarity_Jaccard.95_no_binary <- readAll(calc(sp.weighted.stack_Jaccard.95_no_binary, fun = sum)) 
# Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each species)
sp.mean.rarity_Jaccard.95_no_binary <- sp.rarity_Jaccard.95_no_binary/sp.richness_Jaccard.95_no_binary

# Add 0 values for empty continental pixels (transformed into min values as 0.10)
temp <- continent_mask
temp[!is.na(sp.mean.rarity_Jaccard.95_no_binary[])] <- sp.mean.rarity_Jaccard.95_no_binary[!is.na(sp.mean.rarity_Jaccard.95_no_binary[])]
sp.mean.rarity_Jaccard.95_no_binary <- temp

# Save species ranges and species continuous rarity indices
save(sp.range_size_Jaccard.95_no_binary, file = "./outputs/Indices_maps/no_binary/sp.range_size_Jaccard.95_no_binary.RData", version = "2")
saveRDS(sp.range_size_Jaccard.95_no_binary, file = "./outputs/Indices_maps/no_binary/sp.range_size_Jaccard.95_no_binary.rds", version = "2")
save(sp.rarity_indices_Jaccard.95_no_binary, file = "./outputs/Indices_maps/no_binary/sp.rarity_indices_Jaccard.95_no_binary.RData", version = "2")
saveRDS(sp.rarity_indices_Jaccard.95_no_binary, file = "./outputs/Indices_maps/no_binary/sp.rarity_indices_Jaccard.95_no_binary.rds", version = "2")

# Save full and mean rarity maps
save(sp.rarity_Jaccard.95_no_binary, file = "./outputs/Indices_maps/no_binary/sp.rarity_Jaccard.95_no_binary.RData", version = "2")
saveRDS(sp.rarity_Jaccard.95_no_binary, file = "./outputs/Indices_maps/no_binary/sp.rarity_Jaccard.95_no_binary.rds", version = "2")
save(sp.mean.rarity_Jaccard.95_no_binary, file = "./outputs/Indices_maps/no_binary/sp.mean.rarity_Jaccard.95_no_binary.RData", version = "2")
saveRDS(sp.mean.rarity_Jaccard.95_no_binary, file = "./outputs/Indices_maps/no_binary/sp.mean.rarity_Jaccard.95_no_binary.rds", version = "2")

# Plot distribution of ranges and rarity indices
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.95_no_binary/sp_range_distri_Jaccard.95_no_binary.pdf"), height = 5.3, width = 6.5)
hist(sp.range_size_Jaccard.95_no_binary, main = "Distribution of species geographic ranges\nJaccard.95_no_binary", breaks = 20, xlab = "Geographic ranges in [10^3] km²")
range_threshold_Jaccard.95_no_binary <- round(quantile(sp.range_size_Jaccard.95_no_binary, p=0.25), 1) # Range threshold value for qualitative rarity
abline(v = range_threshold_Jaccard.95_no_binary, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",range_threshold_Jaccard.95_no_binary," km²"), col = "red", lty = 2, lwd = 2, x = "topright", cex = 1, bty ="n")
dev.off()

pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.95_no_binary/sp_rarity_indices_distri_Jaccard.95_no_binary.pdf"), height = 5.3, width = 6.5)
hist(sp.rarity_indices_Jaccard.95_no_binary, main = "Distribution of species rarity indices\nJaccard.95_no_binary", breaks = 20, xlab = "Rarity indices")
indice_threshold_Jaccard.95_no_binary <- round(quantile(sp.rarity_indices_Jaccard.95_no_binary, p=0.75), 3) # Threshold value for qualitative rarity
abline(v = indice_threshold_Jaccard.95_no_binary, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",indice_threshold_Jaccard.95_no_binary), col = "red", lty = 2, lwd = 2, x = "topleft", cex = 1, bty ="n")
dev.off()


## For TSS.80_no_binary

# Compute weights based on geographic ranges
sp.range_size_TSS.80_no_binary <- NA
for (i in 1:nlayers(All_sp_proba_stack_TSS.80_no_binary)) {
  # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
  # Compute the estimated number of pixels occupied by the species as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
  sp.range_size_TSS.80_no_binary[i] <- sum(All_sp_proba_stack_TSS.80_no_binary[[i]]@data@values, na.rm = T)*774.51/1000
}
sp.rarity_indices_TSS.80_no_binary <- 1-(sp.range_size_TSS.80_no_binary/max(sp.range_size_TSS.80_no_binary)) # Weights by max species extent = rarity indices

# Apply weights
sp.weighted.stack_TSS.80_no_binary <- All_sp_proba_stack_TSS.80_no_binary # Generate the stack to fill
for (i in 1:nlayers(sp.weighted.stack_TSS.80_no_binary)){
  # Multiply species probability of presence by species rarity indices
  sp.weighted.stack_TSS.80_no_binary@layers[[i]]@data@values <- All_sp_proba_stack_TSS.80_no_binary@layers[[i]]@data@values * sp.rarity_indices_TSS.80_no_binary[i]
}

plot(sp.weighted.stack_TSS.80_no_binary)

# Sum = Richness weighted by rarity based on Range-size
sp.rarity_TSS.80_no_binary <- readAll(calc(sp.weighted.stack_TSS.80_no_binary, fun = sum)) 
# Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each species)
sp.mean.rarity_TSS.80_no_binary <- sp.rarity_TSS.80_no_binary/sp.richness_TSS.80_no_binary

# Add 0 values for empty continental pixels (transformed into min values as 0.10)
temp <- continent_mask
temp[!is.na(sp.mean.rarity_TSS.80_no_binary[])] <- sp.mean.rarity_TSS.80_no_binary[!is.na(sp.mean.rarity_TSS.80_no_binary[])]
sp.mean.rarity_TSS.80_no_binary <- temp

# Save species ranges and species continous rarity indices
save(sp.range_size_TSS.80_no_binary, file = "./outputs/Indices_maps/no_binary/sp.range_size_TSS.80_no_binary.RData", version = "2")
saveRDS(sp.range_size_TSS.80_no_binary, file = "./outputs/Indices_maps/no_binary/sp.range_size_TSS.80_no_binary.rds", version = "2")
save(sp.rarity_indices_TSS.80_no_binary, file = "./outputs/Indices_maps/no_binary/sp.rarity_indices_TSS.80_no_binary.RData", version = "2")
saveRDS(sp.rarity_indices_TSS.80_no_binary, file = "./outputs/Indices_maps/no_binary/sp.rarity_indices_TSS.80_no_binary.rds", version = "2")

# Save full and mean rarity maps
save(sp.rarity_TSS.80_no_binary, file = "./outputs/Indices_maps/no_binary/sp.rarity_TSS.80_no_binary.RData", version = "2")
saveRDS(sp.rarity_TSS.80_no_binary, file = "./outputs/Indices_maps/no_binary/sp.rarity_TSS.80_no_binary.rds", version = "2")
save(sp.mean.rarity_TSS.80_no_binary, file = "./outputs/Indices_maps/no_binary/sp.mean.rarity_TSS.80_no_binary.RData", version = "2")
saveRDS(sp.mean.rarity_TSS.80_no_binary, file = "./outputs/Indices_maps/no_binary/sp.mean.rarity_TSS.80_no_binary.rds", version = "2")

# Plot distribution of ranges and rarity indices
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.80_no_binary/sp_range_distri_TSS.80_no_binary.pdf"), height = 5.3, width = 6.5)
hist(sp.range_size_TSS.80_no_binary, main = "Distribution of species geographic ranges\nTSS.80_no_binary", breaks = 20, xlab = "Geographic ranges in [10^3] km²")
range_threshold_TSS.80_no_binary <- round(quantile(sp.range_size_TSS.80_no_binary, p=0.25), 1) # Range threshold value for qualitative rarity
abline(v = range_threshold_TSS.80_no_binary, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",range_threshold_TSS.80_no_binary," km²"), col = "red", lty = 2, lwd = 2, x = "topright", cex = 1, bty ="n")
dev.off()

pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.80_no_binary/sp_rarity_indices_distri_TSS.80_no_binary.pdf"), height = 5.3, width = 6.5)
hist(sp.rarity_indices_TSS.80_no_binary, main = "Distribution of species rarity indices\nTSS.80_no_binary", breaks = 20, xlab = "Rarity indices")
indice_threshold_TSS.80_no_binary <- round(quantile(sp.rarity_indices_TSS.80_no_binary, p=0.75), 3) # Threshold value for qualitative rarity
abline(v = indice_threshold_TSS.80_no_binary, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",indice_threshold_TSS.80_no_binary), col = "red", lty = 2, lwd = 2, x = "topleft", cex = 1, bty ="n")
dev.off()


## For TSS.95_no_binary

# Compute weights based on geographic ranges
sp.range_size_TSS.95_no_binary <- NA
for (i in 1:nlayers(All_sp_proba_stack_TSS.95_no_binary)) {
  # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
  # Compute the estimated number of pixels occupied by the species as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
  sp.range_size_TSS.95_no_binary[i] <- sum(All_sp_proba_stack_TSS.95_no_binary[[i]]@data@values, na.rm = T)*774.51/1000
}
sp.rarity_indices_TSS.95_no_binary <- 1-(sp.range_size_TSS.95_no_binary/max(sp.range_size_TSS.95_no_binary)) # Weights by max species extent = rarity indices

# Apply weights
sp.weighted.stack_TSS.95_no_binary <- All_sp_proba_stack_TSS.95_no_binary # Generate the stack to fill
for (i in 1:nlayers(sp.weighted.stack_TSS.95_no_binary)){
  # Multiply species probability of presence by species rarity indices
  sp.weighted.stack_TSS.95_no_binary@layers[[i]]@data@values <- All_sp_proba_stack_TSS.95_no_binary@layers[[i]]@data@values * sp.rarity_indices_TSS.95_no_binary[i]
}

plot(sp.weighted.stack_TSS.95_no_binary)

# Sum = Richness weighted by rarity based on Range-size
sp.rarity_TSS.95_no_binary <- readAll(calc(sp.weighted.stack_TSS.95_no_binary, fun = sum)) 
# Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each species)
sp.mean.rarity_TSS.95_no_binary <- sp.rarity_TSS.95_no_binary/sp.richness_TSS.95_no_binary

# Add 0 values for empty continental pixels (transformed into min values as 0.10)
temp <- continent_mask
temp[!is.na(sp.mean.rarity_TSS.95_no_binary[])] <- sp.mean.rarity_TSS.95_no_binary[!is.na(sp.mean.rarity_TSS.95_no_binary[])]
sp.mean.rarity_TSS.95_no_binary <- temp

# Save species ranges and species continuous rarity indices
save(sp.range_size_TSS.95_no_binary, file = "./outputs/Indices_maps/no_binary/sp.range_size_TSS.95_no_binary.RData", version = "2")
saveRDS(sp.range_size_TSS.95_no_binary, file = "./outputs/Indices_maps/no_binary/sp.range_size_TSS.95_no_binary.rds", version = "2")
save(sp.rarity_indices_TSS.95_no_binary, file = "./outputs/Indices_maps/no_binary/sp.rarity_indices_TSS.95_no_binary.RData", version = "2")
saveRDS(sp.rarity_indices_TSS.95_no_binary, file = "./outputs/Indices_maps/no_binary/sp.rarity_indices_TSS.95_no_binary.rds", version = "2")

# Save full and mean rarity maps
save(sp.rarity_TSS.95_no_binary, file = "./outputs/Indices_maps/no_binary/sp.rarity_TSS.95_no_binary.RData", version = "2")
saveRDS(sp.rarity_TSS.95_no_binary, file = "./outputs/Indices_maps/no_binary/sp.rarity_TSS.95_no_binary.rds", version = "2")
save(sp.mean.rarity_TSS.95_no_binary, file = "./outputs/Indices_maps/no_binary/sp.mean.rarity_TSS.95_no_binary.RData", version = "2")
saveRDS(sp.mean.rarity_TSS.95_no_binary, file = "./outputs/Indices_maps/no_binary/sp.mean.rarity_TSS.95_no_binary.rds", version = "2")

# Plot distribution of ranges and rarity indices
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.95_no_binary/sp_range_distri_TSS.95_no_binary.pdf"), height = 5.3, width = 6.5)
hist(sp.range_size_TSS.95_no_binary, main = "Distribution of species geographic ranges\nTSS.95_no_binary", breaks = 20, xlab = "Geographic ranges in [10^3] km²")
range_threshold_TSS.95_no_binary <- round(quantile(sp.range_size_TSS.95_no_binary, p=0.25), 1) # Range threshold value for qualitative rarity
abline(v = range_threshold_TSS.95_no_binary, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",range_threshold_TSS.95_no_binary," km²"), col = "red", lty = 2, lwd = 2, x = "topright", cex = 1, bty ="n")
dev.off()

pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.95_no_binary/sp_rarity_indices_distri_TSS.95_no_binary.pdf"), height = 5.3, width = 6.5)
hist(sp.rarity_indices_TSS.95_no_binary, main = "Distribution of species rarity indices\nTSS.95_no_binary", breaks = 20, xlab = "Rarity indices")
indice_threshold_TSS.95_no_binary <- round(quantile(sp.rarity_indices_TSS.95_no_binary, p=0.75), 3) # Threshold value for qualitative rarity
abline(v = indice_threshold_TSS.95_no_binary, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",indice_threshold_TSS.95_no_binary), col = "red", lty = 2, lwd = 2, x = "topleft", cex = 1, bty ="n")
dev.off()



### 9.2/ Plot community rarity maps ####

### Load directly the community rarity maps
sp.rarity_Jaccard.80_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/sp.rarity_Jaccard.80_no_binary.rds")
sp.rarity_Jaccard.95_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/sp.rarity_Jaccard.95_no_binary.rds")
sp.rarity_TSS.80_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/sp.rarity_TSS.80_no_binary.rds")
sp.rarity_TSS.95_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/sp.rarity_TSS.95_no_binary.rds")


### Individual plots

# Jaccard.80_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.80_no_binary/sp.rarity_Jaccard.80_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.rarity_Jaccard.80_no_binary, col = pal_bl_red, main = paste0("Species rarity \nJaccard.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_Jaccard.80_no_binary, locs = seq(0, 70, 10), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")
par(mar = internal_margins)
dev.off()

# Jaccard.95_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.95_no_binary/sp.rarity_Jaccard.95_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.rarity_Jaccard.95_no_binary, col = pal_bl_red, main = paste0("Species rarity \nJaccard.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_Jaccard.95_no_binary, locs = seq(0, 70, 10), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")
par(mar = internal_margins)
dev.off()


# TSS.80_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.80_no_binary/sp.rarity_TSS.80_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.rarity_TSS.80_no_binary, col = pal_bl_red, main = paste0("Species rarity \nTSS.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_TSS.80_no_binary, locs = seq(0, 70, 10), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")
par(mar = internal_margins)
dev.off()


# TSS.95_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.95_no_binary/sp.rarity_TSS.95_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.rarity_TSS.95_no_binary, col = pal_bl_red, main = paste0("Species rarity \nTSS.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_TSS.95_no_binary, locs = seq(0, 70, 10), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(sp.rarity_Jaccard.80_no_binary) +
#   tm_raster(palette = pal_bl_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/no_binary/sp.rarity_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.rarity_Jaccard.80_no_binary, col = pal_bl_red, main = paste0("Species rarity \nJaccard.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_Jaccard.80_no_binary, locs = seq(0, 70, 10), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")

image(sp.rarity_Jaccard.95_no_binary, col = pal_bl_red, main = paste0("Species rarity \nJaccard.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_Jaccard.95_no_binary, locs = seq(0, 70, 10), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")

image(sp.rarity_TSS.80_no_binary, col = pal_bl_red, main = paste0("Species rarity \nTSS.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_TSS.80_no_binary, locs = seq(0, 70, 10), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")

image(sp.rarity_TSS.95_no_binary, col = pal_bl_red, main = paste0("Species rarity \nTSS.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_TSS.95_no_binary, locs = seq(0, 70, 10), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")

par(mar = internal_margins)

dev.off()


### 9.3/ Plot mean species rarity maps ####

### Load directly the mean species rarity maps
sp.mean.rarity_Jaccard.80_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/sp.mean.rarity_Jaccard.80_no_binary.rds")
sp.mean.rarity_Jaccard.95_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/sp.mean.rarity_Jaccard.95_no_binary.rds")
sp.mean.rarity_TSS.80_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/sp.mean.rarity_TSS.80_no_binary.rds")
sp.mean.rarity_TSS.95_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/sp.mean.rarity_TSS.95_no_binary.rds")


### Individual plots

# Jaccard.80_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.80_no_binary/sp.mean.rarity_Jaccard.80_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.mean.rarity_Jaccard.80_no_binary, col = pal_bl_red, main = paste0("Mean species rarity\nJaccard.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_Jaccard.80_no_binary, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")
par(mar = internal_margins)
dev.off()

# Jaccard.95_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.95_no_binary/sp.mean.rarity_Jaccard.95_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.mean.rarity_Jaccard.95_no_binary, col = pal_bl_red, main = paste0("Mean species rarity\nJaccard.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_Jaccard.95_no_binary, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")
par(mar = internal_margins)
dev.off()


# TSS.80_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.80_no_binary/sp.mean.rarity_TSS.80_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.mean.rarity_TSS.80_no_binary, col = pal_bl_red, main = paste0("Mean species rarity\nTSS.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_TSS.80_no_binary, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")
par(mar = internal_margins)
dev.off()

# TSS.95_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.95_no_binary/sp.mean.rarity_TSS.95_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.mean.rarity_TSS.95_no_binary, col = pal_bl_red, main = paste0("Mean species rarity \nTSS.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_TSS.95_no_binary, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(sp.mean.rarity_Jaccard.80_no_binary) +
#   tm_raster(palette = pal_bl_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/no_binary/sp.mean.rarity_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.mean.rarity_Jaccard.80_no_binary, col = pal_bl_red, main = paste0("Mean species rarity \nJaccard.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_Jaccard.80_no_binary, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")

image(sp.mean.rarity_Jaccard.95_no_binary, col = pal_bl_red, main = paste0("Mean species rarity\nJaccard.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_Jaccard.95_no_binary, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")

image(sp.mean.rarity_TSS.80_no_binary, col = pal_bl_red, main = paste0("Mean species rarity\nTSS.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_TSS.80_no_binary, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")

image(sp.mean.rarity_TSS.95_no_binary, col = pal_bl_red, main = paste0("Mean species rarity\nTSS.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_TSS.95_no_binary, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")

par(mar = internal_margins)

dev.off()

# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])



###### 10/ Categorical Species Rarity (25% threshold) ######

### Load directly the stack of sp probas use to compute rarity based on geographic range
All_sp_proba_stack_Jaccard.80_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_Jaccard.80_no_binary.rds"))
All_sp_proba_stack_Jaccard.95_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_Jaccard.95_no_binary.rds"))
All_sp_proba_stack_TSS.80_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_TSS.80_no_binary.rds"))
All_sp_proba_stack_TSS.95_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_TSS.95_no_binary.rds"))

### Load species ranges
sp.range_size_Jaccard.80_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/sp.range_size_Jaccard.80_no_binary.rds")
sp.range_size_Jaccard.95_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/sp.range_size_Jaccard.95_no_binary.rds")
sp.range_size_TSS.80_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/sp.range_size_TSS.80_no_binary.rds")
sp.range_size_TSS.95_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/sp.range_size_TSS.95_no_binary.rds")

### 10.1/ Extract rarity threshold ####
rarity_threshold_Jaccard.80_no_binary <- round(quantile(sp.range_size_Jaccard.80_no_binary, probs = 0.25), 3) ; rarity_threshold_Jaccard.80_no_binary
# Threshold surface to be considered a rare species = 148 112 km² 
rarity_threshold_Jaccard.95_no_binary <- quantile(sp.range_size_Jaccard.95_no_binary, probs = 0.25) ; rarity_threshold_Jaccard.95_no_binary
# Threshold surface to be considered a rare species = 198 331 km² 
rarity_threshold_TSS.80_no_binary <- quantile(sp.range_size_TSS.80_no_binary, probs = 0.25) ; rarity_threshold_TSS.80_no_binary
# Threshold surface to be considered a rare species = 148 112 km² 
rarity_threshold_TSS.95_no_binary <- quantile(sp.range_size_TSS.95_no_binary, probs = 0.25) ; rarity_threshold_TSS.95_no_binary
# Threshold surface to be considered a rare species = 198 411 km² 

### 10.2/ Plot species range distribution ####

# Plot species range distribution for Jaccard.80_no_binary
ordered_sp.range_size_Jaccard.80_no_binary <- sp.range_size_Jaccard.80_no_binary[order(sp.range_size_Jaccard.80_no_binary)]

pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.80_no_binary/sp_range_distri_barplot_Jaccard.80_no_binary.pdf"), height = 5.3, width = 6.5)
barplot(ordered_sp.range_size_Jaccard.80_no_binary, ylim = c(0,9000), ylab = "Surface in 10^3 km", xlab = "Species", main = "Distribution of Species range size")
abline(v = length(ordered_sp.range_size_Jaccard.80_no_binary)/4, col = "red", lwd = 2, lty = 2)
legend(legend = c("Threshold 25 %", paste0("    ",rarity_threshold_Jaccard.80_no_binary*1000," km²")), bty = "n", x = "top", text.col = "red", text.font = 2)
legend(legend = c("  RARE", "Species"), bty = "n", x = "left", text.col = "darkblue", text.font = 2)
dev.off()

# Plot species range distribution for Jaccard.95_no_binary
ordered_sp.range_size_Jaccard.95_no_binary <- sp.range_size_Jaccard.95_no_binary[order(sp.range_size_Jaccard.95_no_binary)]

pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.95_no_binary/sp_range_distri_barplot_Jaccard.95_no_binary.pdf"), height = 5.3, width = 6.5)
barplot(ordered_sp.range_size_Jaccard.95_no_binary, ylim = c(0,9000), ylab = "Surface in 10^3 km", xlab = "Species", main = "Distribution of Species range size")
abline(v = length(ordered_sp.range_size_Jaccard.95_no_binary)/4, col = "red", lwd = 2, lty = 2)
legend(legend = c("Threshold 25 %", paste0("    ",rarity_threshold_Jaccard.95_no_binary*1000," km²")), bty = "n", x = "top", text.col = "red", text.font = 2)
legend(legend = c("  RARE", "Species"), bty = "n", x = "left", text.col = "darkblue", text.font = 2)
dev.off()

# Plot species range distribution for TSS.80_no_binary
ordered_sp.range_size_TSS.80_no_binary <- sp.range_size_TSS.80_no_binary[order(sp.range_size_TSS.80_no_binary)]

pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.80_no_binary/sp_range_distri_barplot_TSS.80_no_binary.pdf"), height = 5.3, width = 6.5)
barplot(ordered_sp.range_size_TSS.80_no_binary, ylim = c(0,9000), ylab = "Surface in 10^3 km", xlab = "Species", main = "Distribution of Species range size")
abline(v = length(ordered_sp.range_size_TSS.80_no_binary)/4, col = "red", lwd = 2, lty = 2)
legend(legend = c("Threshold 25 %", paste0("    ",rarity_threshold_TSS.80_no_binary*1000," km²")), bty = "n", x = "top", text.col = "red", text.font = 2)
legend(legend = c("  RARE", "Species"), bty = "n", x = "left", text.col = "darkblue", text.font = 2)
dev.off()

# Plot species range distribution for TSS.95_no_binary
ordered_sp.range_size_TSS.95_no_binary <- sp.range_size_TSS.95_no_binary[order(sp.range_size_TSS.95_no_binary)]

pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.95_no_binary/sp_range_distri_barplot_TSS.95_no_binary.pdf"), height = 5.3, width = 6.5)
barplot(ordered_sp.range_size_TSS.95_no_binary, ylim = c(0,9000), ylab = "Surface in 10^3 km", xlab = "Species", main = "Distribution of Species range size")
abline(v = length(ordered_sp.range_size_TSS.95_no_binary)/4, col = "red", lwd = 2, lty = 2)
legend(legend = c("Threshold 25 %", paste0("    ",rarity_threshold_TSS.95_no_binary*1000," km²")), bty = "n", x = "top", text.col = "red", text.font = 2)
legend(legend = c("  RARE", "Species"), bty = "n", x = "left", text.col = "darkblue", text.font = 2)
dev.off()

### 10.3/ Extract only rare species to compute categorical rarity indices ####

# Detect rare species
rare.sp_Jaccard.80_no_binary <- sp.range_size_Jaccard.80_no_binary < rarity_threshold_Jaccard.80_no_binary
rare.sp_Jaccard.95_no_binary <- sp.range_size_Jaccard.95_no_binary < rarity_threshold_Jaccard.95_no_binary
rare.sp_TSS.80_no_binary <- sp.range_size_TSS.80_no_binary < rarity_threshold_TSS.80_no_binary
rare.sp_TSS.95_no_binary <- sp.range_size_TSS.95_no_binary < rarity_threshold_TSS.95_no_binary

sum(rare.sp_Jaccard.80_no_binary) # 59 rare species on 236 (25%)

# Extract only rare species
Rare_sp_proba_stack_Jaccard.80_no_binary <- All_sp_proba_stack_Jaccard.80_no_binary[[which(rare.sp_Jaccard.80_no_binary)]]
Rare_sp_proba_stack_Jaccard.95_no_binary <- All_sp_proba_stack_Jaccard.95_no_binary[[which(rare.sp_Jaccard.95_no_binary)]]
Rare_sp_proba_stack_TSS.80_no_binary <- All_sp_proba_stack_TSS.80_no_binary[[which(rare.sp_TSS.80_no_binary)]]
Rare_sp_proba_stack_TSS.95_no_binary <- All_sp_proba_stack_TSS.95_no_binary[[which(rare.sp_TSS.95_no_binary)]]

### 10.4/ Index computation ####

# Load maps of Species richness to use to compute mean rarity index among local species
sp.richness_Jaccard.80_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_Jaccard.80_no_binary.rds"))
sp.richness_Jaccard.95_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_Jaccard.95_no_binary.rds"))
sp.richness_TSS.80_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_TSS.80_no_binary.rds"))
sp.richness_TSS.95_no_binary <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/tot.sp.richness_TSS.95_no_binary.rds"))

# Index computation for total of rare species and proportion of rare species
tot.rare.sp_Jaccard.80_no_binary <- readAll(calc(Rare_sp_proba_stack_Jaccard.80_no_binary, fun = sum))
prop.rare.sp_Jaccard.80_no_binary <- tot.rare.sp_Jaccard.80_no_binary/sp.richness_Jaccard.80_no_binary*100
tot.rare.sp_Jaccard.95_no_binary <- readAll(calc(Rare_sp_proba_stack_Jaccard.95_no_binary, fun = sum))
prop.rare.sp_Jaccard.95_no_binary <- tot.rare.sp_Jaccard.95_no_binary/sp.richness_Jaccard.95_no_binary*100
tot.rare.sp_TSS.80_no_binary <- readAll(calc(Rare_sp_proba_stack_TSS.80_no_binary, fun = sum))
prop.rare.sp_TSS.80_no_binary <- tot.rare.sp_TSS.80_no_binary/sp.richness_TSS.80_no_binary*100
tot.rare.sp_TSS.95_no_binary <- readAll(calc(Rare_sp_proba_stack_TSS.95_no_binary, fun = sum))
prop.rare.sp_TSS.95_no_binary <- tot.rare.sp_TSS.95_no_binary/sp.richness_TSS.95_no_binary*100

# Save maps
save(tot.rare.sp_Jaccard.80_no_binary , file = "./outputs/Indices_maps/no_binary/tot.rare.sp_Jaccard.80_no_binary.RData", version = "2")
saveRDS(tot.rare.sp_Jaccard.80_no_binary , file = "./outputs/Indices_maps/no_binary/tot.rare.sp_Jaccard.80_no_binary.rds", version = "2")
save(tot.rare.sp_Jaccard.95_no_binary , file = "./outputs/Indices_maps/no_binary/tot.rare.sp_Jaccard.95_no_binary.RData", version = "2")
saveRDS(tot.rare.sp_Jaccard.95_no_binary , file = "./outputs/Indices_maps/no_binary/tot.rare.sp_Jaccard.95_no_binary.rds", version = "2")
save(tot.rare.sp_TSS.80_no_binary , file = "./outputs/Indices_maps/no_binary/tot.rare.sp_TSS.80_no_binary.RData", version = "2")
saveRDS(tot.rare.sp_TSS.80_no_binary , file = "./outputs/Indices_maps/no_binary/tot.rare.sp_TSS.80_no_binary.rds", version = "2")
save(tot.rare.sp_TSS.95_no_binary , file = "./outputs/Indices_maps/no_binary/tot.rare.sp_TSS.95_no_binary.RData", version = "2")
saveRDS(tot.rare.sp_TSS.95_no_binary , file = "./outputs/Indices_maps/no_binary/tot.rare.sp_TSS.95_no_binary.rds", version = "2")

save(prop.rare.sp_Jaccard.80_no_binary , file = "./outputs/Indices_maps/no_binary/prop.rare.sp_Jaccard.80_no_binary.RData", version = "2")
saveRDS(prop.rare.sp_Jaccard.80_no_binary , file = "./outputs/Indices_maps/no_binary/prop.rare.sp_Jaccard.80_no_binary.rds", version = "2")
save(prop.rare.sp_Jaccard.95_no_binary , file = "./outputs/Indices_maps/no_binary/prop.rare.sp_Jaccard.95_no_binary.RData", version = "2")
saveRDS(prop.rare.sp_Jaccard.95_no_binary , file = "./outputs/Indices_maps/no_binary/prop.rare.sp_Jaccard.95_no_binary.rds", version = "2")
save(prop.rare.sp_TSS.80_no_binary , file = "./outputs/Indices_maps/no_binary/prop.rare.sp_TSS.80_no_binary.RData", version = "2")
saveRDS(prop.rare.sp_TSS.80_no_binary , file = "./outputs/Indices_maps/no_binary/prop.rare.sp_TSS.80_no_binary.rds", version = "2")
save(prop.rare.sp_TSS.95_no_binary , file = "./outputs/Indices_maps/no_binary/prop.rare.sp_TSS.95_no_binary.RData", version = "2")
saveRDS(prop.rare.sp_TSS.95_no_binary , file = "./outputs/Indices_maps/no_binary/prop.rare.sp_TSS.95_no_binary.rds", version = "2")

### 10.5/ Plots of rare species richness ####

### Load directly the final Rare Species richness layer
tot.rare.sp_Jaccard.80_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/tot.rare.sp_Jaccard.80_no_binary.rds")
tot.rare.sp_Jaccard.95_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/tot.rare.sp_Jaccard.95_no_binary.rds")
tot.rare.sp_TSS.80_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/tot.rare.sp_TSS.80_no_binary.rds")
tot.rare.sp_TSS.95_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/tot.rare.sp_TSS.95_no_binary.rds")


### Individual plots

# Jaccard.80_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.80_no_binary/tot.rare.sp_Jaccard.80_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.rare.sp_Jaccard.80_no_binary, col = pal_bl_red, main = paste0("Rare species richness \nJaccard.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_Jaccard.80_no_binary, locs = seq(0, 14, 2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 7, font = 2, cex = 1.2, label = "Rare\nspecies")
par(mar = internal_margins)
dev.off()

# Jaccard.95_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.95_no_binary/tot.rare.sp_Jaccard.95_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.rare.sp_Jaccard.95_no_binary, col = pal_bl_red, main = paste0("Rare species richness \nJaccard.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_Jaccard.95_no_binary, locs = seq(0, 12, 2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 7, font = 2, cex = 1.2, label = "Rare\nspecies")
par(mar = internal_margins)
dev.off()

# TSS.80_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.80_no_binary/tot.rare.sp_TSS.80_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.rare.sp_TSS.80_no_binary, col = pal_bl_red, main = paste0("Rare species richness \nTSS.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_TSS.80_no_binary, locs = seq(0, 14, 2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 7, font = 2, cex = 1.2, label = "Rare\nspecies")
par(mar = internal_margins)
dev.off()

# TSS.95_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.95_no_binary/tot.rare.sp_TSS.95_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.rare.sp_TSS.95_no_binary, col = pal_bl_red, main = paste0("Rare species richness \nTSS.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_TSS.95_no_binary, locs = seq(0, 12, 2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 7, font = 2, cex = 1.2, label = "Rare\nspecies")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(tot.sp.richness_Jaccard.80_no_binary) +
#   tm_raster(palette = pal_bl_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/no_binary/tot.rare.sp_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(tot.rare.sp_Jaccard.80_no_binary, col = pal_bl_red, main = paste0("Rare species richness \nJaccard.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_Jaccard.80_no_binary, locs = seq(0, 14, 2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 7, font = 2, cex = 1.2, label = "Rare\nspecies")

image(tot.rare.sp_Jaccard.95_no_binary, col = pal_bl_red, main = paste0("Rare species richness \nJaccard.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_Jaccard.95_no_binary, locs = seq(0, 12, 2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 7, font = 2, cex = 1.2, label = "Rare\nspecies")

image(tot.rare.sp_TSS.80_no_binary, col = pal_bl_red, main = paste0("Rare species richness \nTSS.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_TSS.80_no_binary, locs = seq(0, 14, 2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 7, font = 2, cex = 1.2, label = "Rare\nspecies")

image(tot.rare.sp_TSS.95_no_binary, col = pal_bl_red, main = paste0("Rare species richness \nTSS.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_TSS.95_no_binary, locs = seq(0, 12, 2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 7, font = 2, cex = 1.2, label = "Rare\nspecies")

par(mar = internal_margins)

dev.off()


### 10.5 /Plots of proportion of rare species ####

### Load directly the final Rare Species richness layer
prop.rare.sp_Jaccard.80_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/prop.rare.sp_Jaccard.80_no_binary.rds")
prop.rare.sp_Jaccard.95_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/prop.rare.sp_Jaccard.95_no_binary.rds")
prop.rare.sp_TSS.80_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/prop.rare.sp_TSS.80_no_binary.rds")
prop.rare.sp_TSS.95_no_binary <- readRDS(file = "./outputs/Indices_maps/no_binary/prop.rare.sp_TSS.95_no_binary.rds")

hist(prop.rare.sp_Jaccard.80_no_binary)

# Add 0 values for empty continental pixels (transformed into min values as 0.25)
temp <- continent_mask
temp[!is.na(prop.rare.sp_Jaccard.80_no_binary[])] <- prop.rare.sp_Jaccard.80_no_binary[!is.na(prop.rare.sp_Jaccard.80_no_binary)]
prop.rare.sp_Jaccard.80_no_binary <- temp

temp <- continent_mask
temp[!is.na(prop.rare.sp_Jaccard.95_no_binary[])] <- prop.rare.sp_Jaccard.95_no_binary[!is.na(prop.rare.sp_Jaccard.95_no_binary)]
prop.rare.sp_Jaccard.95_no_binary <- temp

temp <- continent_mask
temp[!is.na(prop.rare.sp_TSS.80_no_binary[])] <- prop.rare.sp_TSS.80_no_binary[!is.na(prop.rare.sp_TSS.80_no_binary)]
prop.rare.sp_TSS.80_no_binary <- temp

temp <- continent_mask
temp[!is.na(prop.rare.sp_TSS.95_no_binary[])] <- prop.rare.sp_TSS.95_no_binary[!is.na(prop.rare.sp_TSS.95_no_binary)]
prop.rare.sp_TSS.95_no_binary <- temp

# Contrast by merging all values higher than 25 %
prop.rare.sp_Jaccard.80_no_binary[prop.rare.sp_Jaccard.80_no_binary > 25] <-  25
prop.rare.sp_Jaccard.95_no_binary[prop.rare.sp_Jaccard.95_no_binary > 25] <-  25
prop.rare.sp_TSS.80_no_binary[prop.rare.sp_TSS.80_no_binary > 25] <-  25
prop.rare.sp_TSS.95_no_binary[prop.rare.sp_TSS.95_no_binary > 25] <-  25


### Individual plots

# Jaccard.80_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.80_no_binary/prop.rare.sp_Jaccard.80_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(prop.rare.sp_Jaccard.80_no_binary, col = pal_bl_red, main = paste0("Proportion of rare species\nJaccard.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_Jaccard.80_no_binary, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")
par(mar = internal_margins)
dev.off()

# Jaccard.95_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/Jaccard.95_no_binary/prop.rare.sp_Jaccard.95_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(prop.rare.sp_Jaccard.95_no_binary, col = pal_bl_red, main = paste0("Proportion of rare species\nJaccard.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_Jaccard.95_no_binary, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")
par(mar = internal_margins)
dev.off()

# TSS.80_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.80_no_binary/prop.rare.sp_TSS.80_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(prop.rare.sp_TSS.80_no_binary, col = pal_bl_red, main = paste0("Proportion of rare species\nTSS.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_TSS.80_no_binary, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")
par(mar = internal_margins)
dev.off()

# TSS.95_no_binary
pdf(file = paste0("./maps/Indices_maps/no_binary/TSS.95_no_binary/prop.rare.sp_TSS.95_no_binary.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(prop.rare.sp_TSS.95_no_binary, col = pal_bl_red, main = paste0("Proportion of rare species\nTSS.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_TSS.95_no_binary, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(prop.sp.richness_Jaccard.80_no_binary) +
#   tm_raster(palette = pal_bl_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/no_binary/prop.rare.sp_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(prop.rare.sp_Jaccard.80_no_binary, col = pal_bl_red, main = paste0("Proportion of rare species\nJaccard.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_Jaccard.80_no_binary, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")

image(prop.rare.sp_Jaccard.95_no_binary, col = pal_bl_red, main = paste0("Proportion of rare species\nJaccard.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_Jaccard.95_no_binary, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")

image(prop.rare.sp_TSS.80_no_binary, col = pal_bl_red, main = paste0("Proportion of rare species\nTSS.80_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_TSS.80_no_binary, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")

image(prop.rare.sp_TSS.95_no_binary, col = pal_bl_red, main = paste0("Proportion of rare species\nTSS.95_no_binary"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_TSS.95_no_binary, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")

par(mar = internal_margins)

dev.off()

# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])



###### 11/ Continuous Mimicry rarity = Range size weighted mimicry richness ######

# Mimicry-level index. Not very useful to control for absence of species not in the phylogeny


###### 12/ Community vulnerability ######

# Computed on the basis of mimicry ring richness. Not very useful to control for absence of species not in the phylogeny


