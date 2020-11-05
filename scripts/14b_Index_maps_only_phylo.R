
##### Script 14b: Compute indices maps but only with species included in the phylogeny #####

## Create 9 index maps (no need to recompute phylogeny-based indices because they are already made only with species included in the phylogeny)

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


###### 0/ Generate the list of species included in the phylogeny to keep for these maps

### Load directly the phylogeny with the good labels
phylo.Ithomiini <- readRDS(file = "./input_data/Phylogenies/Final_phylogeny.rds")

# Name of the species to keep in all computation
phylo.Ithomiini$tip.label

# Extract only the 339 species included in the phylogeny from list.sp and the sp.stack
list.sp_phylo <- list.sp[list.sp$Sp_full %in% phylo.Ithomiini$tip.label,]


###### 1/ Species richness ######

### 1.1/ Species Stack generation ####

### Load directly the complete stack of species proba
All_sp_proba_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80_phylo.rds"))
All_sp_proba_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.95_phylo.rds"))
All_sp_proba_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80_phylo.rds"))
All_sp_proba_stack_TSS.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.95_phylo.rds"))

# Keep only species in the phylogeny
All_sp_proba_stack_Jaccard.80_phylo <- All_sp_proba_stack_Jaccard.80[[phylo.Ithomiini$tip.label]]
All_sp_proba_stack_Jaccard.95_phylo <- All_sp_proba_stack_Jaccard.95[[phylo.Ithomiini$tip.label]]
All_sp_proba_stack_TSS.80_phylo <- All_sp_proba_stack_TSS.80[[phylo.Ithomiini$tip.label]]
All_sp_proba_stack_TSS.95_phylo <- All_sp_proba_stack_TSS.95[[phylo.Ithomiini$tip.label]]

# Save new stacks
save(All_sp_proba_stack_Jaccard.80_phylo, file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_Jaccard.80_phylo.RData"), version = "2")
saveRDS(All_sp_proba_stack_Jaccard.80_phylo, file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_Jaccard.80_phylo.rds"), version = "2")
save(All_sp_proba_stack_Jaccard.95_phylo, file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_Jaccard.95_phylo.RData"), version = "2")
saveRDS(All_sp_proba_stack_Jaccard.95_phylo, file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_Jaccard.95_phylo.rds"), version = "2")
save(All_sp_proba_stack_Jaccard.80_phylo, file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_TSS.80_phylo.RData"), version = "2")
saveRDS(All_sp_proba_stack_Jaccard.80_phylo, file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_TSS.80_phylo.rds"), version = "2")
save(All_sp_proba_stack_TSS.95_phylo, file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_TSS.95_phylo.RData"), version = "2")
saveRDS(All_sp_proba_stack_TSS.95_phylo, file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_TSS.95_phylo.rds"), version = "2")


### 1.2/ Index computation ####
tot.sp.richness_Jaccard.80_phylo <- readAll(calc(All_sp_proba_stack_Jaccard.80_phylo, fun = sum))
tot.sp.richness_Jaccard.95_phylo <- readAll(calc(All_sp_proba_stack_Jaccard.95_phylo, fun = sum))
tot.sp.richness_TSS.80_phylo <- readAll(calc(All_sp_proba_stack_TSS.80_phylo, fun = sum))
tot.sp.richness_TSS.95_phylo <- readAll(calc(All_sp_proba_stack_TSS.95_phylo, fun = sum))


# Save
save(tot.sp.richness_Jaccard.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_Jaccard.80_phylo.RData"), version = "2")
saveRDS(tot.sp.richness_Jaccard.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_Jaccard.80_phylo.rds"), version = "2")
save(tot.sp.richness_Jaccard.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_Jaccard.95_phylo.RData"), version = "2")
saveRDS(tot.sp.richness_Jaccard.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_Jaccard.95_phylo.rds"), version = "2")
save(tot.sp.richness_TSS.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_TSS.80_phylo.RData"), version = "2")
saveRDS(tot.sp.richness_TSS.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_TSS.80_phylo.rds"), version = "2")
save(tot.sp.richness_TSS.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_TSS.95_phylo.RData"), version = "2")
saveRDS(tot.sp.richness_TSS.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_TSS.95_phylo.rds"), version = "2")


### Load directly the final Species richness layer
tot.sp.richness_Jaccard.80_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_Jaccard.80_phylo.rds"))
tot.sp.richness_Jaccard.95_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_Jaccard.95_phylo.rds"))
tot.sp.richness_TSS.80_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_TSS.80_phylo.rds"))
tot.sp.richness_TSS.95_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_TSS.95_phylo.rds"))

### 1.3/ Plot Species richness ####

### Individual plots

# Jaccard.80_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.80_phylo/tot.sp.richness_Jaccard.80_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.sp.richness_Jaccard.80_phylo, col = pal_bl_red, main = paste0("Species richness \nJaccard.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_Jaccard.80_phylo, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# Jaccard.95_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.95_phylo/tot.sp.richness_Jaccard.95_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.sp.richness_Jaccard.95_phylo, col = pal_bl_red, main = paste0("Species richness \nJaccard.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_Jaccard.95_phylo, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# TSS.80_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.80_phylo/tot.sp.richness_TSS.80_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.sp.richness_TSS.80_phylo, col = pal_bl_red, main = paste0("Species richness \nTSS.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_TSS.80_phylo, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# TSS.95_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.95_phylo/tot.sp.richness_TSS.95_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.sp.richness_TSS.95_phylo, col = pal_bl_red, main = paste0("Species richness \nTSS.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_TSS.95_phylo, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(tot.sp.richness_Jaccard.80_phylo) +
#   tm_raster(palette = pal_bl_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/phylo_only/tot.sp.richness_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(tot.sp.richness_Jaccard.80_phylo, col = pal_bl_red, main = paste0("Species richness \nJaccard.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_Jaccard.80_phylo, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")

image(tot.sp.richness_Jaccard.95_phylo, col = pal_bl_red, main = paste0("Species richness \nJaccard.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_Jaccard.95_phylo, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")

image(tot.sp.richness_TSS.80_phylo, col = pal_bl_red, main = paste0("Species richness \nTSS.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_TSS.80_phylo, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")

image(tot.sp.richness_TSS.95_phylo, col = pal_bl_red, main = paste0("Species richness \nTSS.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_TSS.95_phylo, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
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
All_sp_proba_stack_Jaccard.80_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_Jaccard.80_phylo.rds"))
All_sp_proba_stack_Jaccard.95_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_Jaccard.95_phylo.rds"))
All_sp_proba_stack_TSS.80_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_TSS.80_phylo.rds"))
All_sp_proba_stack_TSS.95_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_TSS.95_phylo.rds"))

### 2.1/ Index computation ####
sp.diversity_Jaccard.80_phylo <- calc(All_sp_proba_stack_Jaccard.80_phylo, fun = shannon)*1
sp.diversity_Jaccard.95_phylo <- calc(All_sp_proba_stack_Jaccard.95_phylo, fun = shannon)*1
sp.diversity_TSS.80_phylo <- calc(All_sp_proba_stack_TSS.80_phylo, fun = shannon)*1
sp.diversity_TSS.95_phylo <- calc(All_sp_proba_stack_TSS.95_phylo, fun = shannon)*1

# Save
save(sp.diversity_Jaccard.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity_Jaccard.80_phylo.RData"), version = "2")
saveRDS(sp.diversity_Jaccard.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity_Jaccard.80_phylo.rds"), version = "2")
save(sp.diversity_Jaccard.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity_Jaccard.95_phylo.RData"), version = "2")
saveRDS(sp.diversity_Jaccard.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity_Jaccard.95_phylo.rds"), version = "2")
save(sp.diversity_TSS.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity_TSS.80_phylo.RData"), version = "2")
saveRDS(sp.diversity_TSS.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity_TSS.80_phylo.rds"), version = "2")
save(sp.diversity_TSS.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity_TSS.95_phylo.RData"), version = "2")
saveRDS(sp.diversity_TSS.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity_TSS.95_phylo.rds"), version = "2")


### Load directly the final Species diversity layers
sp.diversity_Jaccard.80_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity_Jaccard.80_phylo.rds"))
sp.diversity_Jaccard.95_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity_Jaccard.95_phylo.rds"))
sp.diversity_TSS.80_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity_TSS.80_phylo.rds"))
sp.diversity_TSS.95_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity_TSS.95_phylo.rds"))

### 2.2/ Plot Species Diversity raw version ####

### Individual plots

# Jaccard.80_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.80_phylo/sp_diversity_Jaccard.80_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity_Jaccard.80_phylo, col = pal_bl_red, main = paste0("Shannon's species diversity \nJaccard.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_Jaccard.80_phylo, locs = seq(0, 5, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# Jaccard.95_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.95_phylo/sp_diversity_Jaccard.95_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity_Jaccard.95_phylo, col = pal_bl_red, main = paste0("Shannon's species diversity \nJaccard.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_Jaccard.95_phylo, locs = seq(0, 5, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# TSS.80_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.80_phylo/sp_diversity_TSS.80_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity_TSS.80_phylo, col = pal_bl_red, main = paste0("Shannon's species diversity \nTSS.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_TSS.80_phylo, locs = seq(0, 5, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# TSS.95_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.95_phylo/sp_diversity_TSS.95_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity_TSS.95_phylo, col = pal_bl_red, main = paste0("Shannon's species diversity \nTSS.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_TSS.95_phylo, locs = seq(0, 5, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# ## Tmap version
# 
# library(tmap)
# 
# tm_shape(sp.diversity_Jaccard.80_phylo) +
#   tm_raster(palette = pal_bl_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")


### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/phylo_only/sp.diversity_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.diversity_Jaccard.80_phylo, col = pal_bl_red, main = paste0("Shannon's species diversity \nJaccard.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_Jaccard.80_phylo, locs = seq(0, 5, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(sp.diversity_Jaccard.95_phylo, col = pal_bl_red, main = paste0("Shannon's species diversity \nJaccard.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_Jaccard.95_phylo, locs = seq(0, 5, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(sp.diversity_TSS.80_phylo, col = pal_bl_red, main = paste0("Shannon's species diversity \nTSS.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_TSS.80_phylo, locs = seq(0, 5, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(sp.diversity_TSS.95_phylo, col = pal_bl_red, main = paste0("Shannon's species diversity \nTSS.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity_TSS.95_phylo, locs = seq(0, 5, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
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
All_sp_proba_stack_Jaccard.80_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_Jaccard.80_phylo.rds"))
All_sp_proba_stack_Jaccard.95_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_Jaccard.95_phylo.rds"))
All_sp_proba_stack_TSS.80_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_TSS.80_phylo.rds"))
All_sp_proba_stack_TSS.95_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_TSS.95_phylo.rds"))


### 3.1/ Index computation ####
sp.diversity.compatible_Jaccard.80_phylo <- calc(All_sp_proba_stack_Jaccard.80_phylo, fun = shannon_compatible)*1 # Raw version
sp.diversity.compatible_Jost_Jaccard.80_phylo <- exp(sp.diversity.compatible_Jaccard.80_phylo) - 1 # Jost's version

sp.diversity.compatible_Jaccard.95_phylo <- calc(All_sp_proba_stack_Jaccard.95_phylo, fun = shannon_compatible)*1 # Raw version
sp.diversity.compatible_Jost_Jaccard.95_phylo <- exp(sp.diversity.compatible_Jaccard.95_phylo) - 1 # Jost's version

sp.diversity.compatible_TSS.80_phylo <- calc(All_sp_proba_stack_TSS.80_phylo, fun = shannon_compatible)*1 # Raw version
sp.diversity.compatible_Jost_TSS.80_phylo <- exp(sp.diversity.compatible_TSS.80_phylo) - 1 # Jost's version

sp.diversity.compatible_TSS.95_phylo <- calc(All_sp_proba_stack_TSS.95_phylo, fun = shannon_compatible)*1 # Raw version
sp.diversity.compatible_Jost_TSS.95_phylo <- exp(sp.diversity.compatible_TSS.95_phylo) - 1 # Jost's version

# Save
save(sp.diversity.compatible_Jaccard.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_Jaccard.80_phylo.RData"), version = "2")
saveRDS(sp.diversity.compatible_Jaccard.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_Jaccard.80_phylo.rds"), version = "2")
save(sp.diversity.compatible_Jaccard.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_Jaccard.95_phylo.RData"), version = "2")
saveRDS(sp.diversity.compatible_Jaccard.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_Jaccard.95_phylo.rds"), version = "2")
save(sp.diversity.compatible_TSS.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_TSS.80_phylo.RData"), version = "2")
saveRDS(sp.diversity.compatible_TSS.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_TSS.80_phylo.rds"), version = "2")
save(sp.diversity.compatible_TSS.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_TSS.95_phylo.RData"), version = "2")
saveRDS(sp.diversity.compatible_TSS.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_TSS.95_phylo.rds"), version = "2")

save(sp.diversity.compatible_Jost_Jaccard.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_Jost_Jaccard.80_phylo.RData"), version = "2")
saveRDS(sp.diversity.compatible_Jost_Jaccard.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_Jost_Jaccard.80_phylo.rds"), version = "2")
save(sp.diversity.compatible_Jost_Jaccard.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_Jost_Jaccard.95_phylo.RData"), version = "2")
saveRDS(sp.diversity.compatible_Jost_Jaccard.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_Jost_Jaccard.95_phylo.rds"), version = "2")
save(sp.diversity.compatible_Jost_TSS.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_Jost_TSS.80_phylo.RData"), version = "2")
saveRDS(sp.diversity.compatible_Jost_TSS.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_Jost_TSS.80_phylo.rds"), version = "2")
save(sp.diversity.compatible_Jost_TSS.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_Jost_TSS.95_phylo.RData"), version = "2")
saveRDS(sp.diversity.compatible_Jost_TSS.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_Jost_TSS.95_phylo.rds"), version = "2")


### Load directly the final Species corrected diversity layers
sp.diversity.compatible_Jaccard.80_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_Jaccard.80_phylo.rds"))
sp.diversity.compatible_Jaccard.95_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_Jaccard.95_phylo.rds"))
sp.diversity.compatible_TSS.80_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_TSS.80_phylo.rds"))
sp.diversity.compatible_TSS.95_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_TSS.95_phylo.rds"))

sp.diversity.compatible_Jost_Jaccard.80_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_Jost_Jaccard.80_phylo.rds"))
sp.diversity.compatible_Jost_Jaccard.95_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_Jost_Jaccard.95_phylo.rds"))
sp.diversity.compatible_Jost_TSS.80_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_Jost_TSS.80_phylo.rds"))
sp.diversity.compatible_Jost_TSS.95_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/sp.diversity.compatible_Jost_TSS.95_phylo.rds"))


### 3.2/ Plot Species Diversity Corrected version ####

### Individual plots

# Jaccard.80_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.80_phylo/sp_diversity_compatible_Jaccard.80_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_Jaccard.80_phylo, col = pal_bl_red, main = paste0("Shannon's species diversity \nJaccard.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jaccard.80_phylo, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# Jaccard.95_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.95_phylo/sp_diversity_compatible_Jaccard.95_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_Jaccard.95_phylo, col = pal_bl_red, main = paste0("Shannon's species diversity \nJaccard.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jaccard.95_phylo, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# TSS.80_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.80_phylo/sp_diversity_compatible_TSS.80_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_TSS.80_phylo, col = pal_bl_red, main = paste0("Shannon's species diversity \nTSS.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_TSS.80_phylo, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# TSS.95_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.95_phylo/sp_diversity_compatible_TSS.95_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_TSS.95_phylo, col = pal_bl_red, main = paste0("Shannon's species diversity \nTSS.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_TSS.95_phylo, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")
par(mar = internal_margins)
dev.off()

# ## Tmap version
# 
# library(tmap)
# 
# tm_shape(sp.diversity.compatible_Jaccard.80_phylo) +
#   tm_raster(palette = pal_bl_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")


### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/phylo_only/sp.diversity.compatible_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.diversity.compatible_Jaccard.80_phylo, col = pal_bl_red, main = paste0("Shannon's species diversity \nJaccard.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jaccard.80_phylo, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(sp.diversity.compatible_Jaccard.95_phylo, col = pal_bl_red, main = paste0("Shannon's species diversity \nJaccard.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jaccard.95_phylo, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(sp.diversity.compatible_TSS.80_phylo, col = pal_bl_red, main = paste0("Shannon's species diversity \nTSS.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_TSS.80_phylo, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

image(sp.diversity.compatible_TSS.95_phylo, col = pal_bl_red, main = paste0("Shannon's species diversity \nTSS.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Shannon's H'", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_TSS.95_phylo, locs = seq(0, 4, 1), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4, font = 2, cex = 1.2, label = "H' index")

par(mar = internal_margins)

dev.off()

### 3.3/ Plot Species Diversity Jost version ####

### Individual plots

# Jaccard.80_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.80_phylo/sp_diversity_compatible_Jost_Jaccard.80_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
image(sp.diversity.compatible_Jost_Jaccard.80_phylo, col = pal_bl_red, main = paste0("Species diversity (Jost's effective species richness) \nJaccard.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_Jaccard.80_phylo, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# Jaccard.95_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.95_phylo/sp_diversity_compatible_Jost_Jaccard.95_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_Jost_Jaccard.95_phylo, col = pal_bl_red, main = paste0("Species diversity (Jost's effective species richness) \nJaccard.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_Jaccard.95_phylo, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# TSS.80_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.80_phylo/sp_diversity_compatible_Jost_TSS.80_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_Jost_TSS.80_phylo, col = pal_bl_red, main = paste0("Species diversity (Jost's effective species richness) \nTSS.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_TSS.80_phylo, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# TSS.95_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.95_phylo/sp_diversity_compatible_Jost_TSS.95_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity.compatible_Jost_TSS.95_phylo, col = pal_bl_red, main = paste0("Species diversity (Jost's effective species richness) \nTSS.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_TSS.95_phylo, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

# ## Tmap version
# 
# library(tmap)
# 
# tm_shape(sp.diversity.compatible_Jost_Jaccard.80_phylo) +
#   tm_raster(palette = pal_bl_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")


### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/phylo_only/sp.diversity.compatible_Jost_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.diversity.compatible_Jost_Jaccard.80_phylo, col = pal_bl_red, main = paste0("Species diversity (Jost's effective species richness) \nJaccard.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_Jaccard.80_phylo, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1, label = "Species")

image(sp.diversity.compatible_Jost_Jaccard.95_phylo, col = pal_bl_red, main = paste0("Species diversity (Jost's effective species richness) \nJaccard.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_Jaccard.95_phylo, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1, label = "Species")

image(sp.diversity.compatible_Jost_TSS.80_phylo, col = pal_bl_red, main = paste0("Species diversity (Jost's effective species richness) \nTSS.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_TSS.80_phylo, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1, label = "Species")

image(sp.diversity.compatible_Jost_TSS.95_phylo, col = pal_bl_red, main = paste0("Species diversity (Jost's effective species richness) \nTSS.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.diversity.compatible_Jost_TSS.95_phylo, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
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
All_sp_proba_stack_Jaccard.80_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_Jaccard.80_phylo.rds"))
All_sp_proba_stack_Jaccard.95_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_Jaccard.95_phylo.rds"))
All_sp_proba_stack_TSS.80_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_TSS.80_phylo.rds"))
All_sp_proba_stack_TSS.95_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_TSS.95_phylo.rds"))

### 4.1/ Index computation ####
sp.evenness_Jaccard.80_phylo <- calc(All_sp_proba_stack_Jaccard.80_phylo, fun = evenness)*1
sp.evenness_Jaccard.95_phylo <- calc(All_sp_proba_stack_Jaccard.95_phylo, fun = evenness)*1
sp.evenness_TSS.80_phylo <- calc(All_sp_proba_stack_TSS.80_phylo, fun = evenness)*1
sp.evenness_TSS.95_phylo <- calc(All_sp_proba_stack_TSS.95_phylo, fun = evenness)*1

# Save
save(sp.evenness_Jaccard.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.evenness_Jaccard.80_phylo.RData"), version = "2")
saveRDS(sp.evenness_Jaccard.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.evenness_Jaccard.80_phylo.rds"), version = "2")
save(sp.evenness_Jaccard.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.evenness_Jaccard.95_phylo.RData"), version = "2")
saveRDS(sp.evenness_Jaccard.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.evenness_Jaccard.95_phylo.rds"), version = "2")
save(sp.evenness_TSS.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.evenness_TSS.80_phylo.RData"), version = "2")
saveRDS(sp.evenness_TSS.80_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.evenness_TSS.80_phylo.rds"), version = "2")
save(sp.evenness_TSS.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.evenness_TSS.95_phylo.RData"), version = "2")
saveRDS(sp.evenness_TSS.95_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.evenness_TSS.95_phylo.rds"), version = "2")

### 4.2/ Contrast maps ####

# Load directly the final Species evenness layer
sp.evenness_Jaccard.80_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/sp.evenness_Jaccard.80_phylo.rds"))
sp.evenness_Jaccard.95_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/sp.evenness_Jaccard.95_phylo.rds"))
sp.evenness_TSS.80_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/sp.evenness_TSS.80_phylo.rds"))
sp.evenness_TSS.95_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/sp.evenness_TSS.95_phylo.rds"))

hist(sp.evenness_Jaccard.80_phylo)

# Contrast values by merging low values below 0.25
sp.evenness_Jaccard.80_phylo_contrasted <- sp.evenness_Jaccard.80_phylo
sp.evenness_Jaccard.80_phylo_contrasted[sp.evenness_Jaccard.80_phylo < 0.25] <- 0.255
sp.evenness_Jaccard.95_phylo_contrasted <- sp.evenness_Jaccard.95_phylo
sp.evenness_Jaccard.95_phylo_contrasted[sp.evenness_Jaccard.95_phylo < 0.25] <- 0.255
sp.evenness_TSS.80_phylo_contrasted <- sp.evenness_TSS.80_phylo
sp.evenness_TSS.80_phylo_contrasted[sp.evenness_TSS.80_phylo < 0.25] <- 0.255
sp.evenness_TSS.95_phylo_contrasted <- sp.evenness_TSS.95_phylo
sp.evenness_TSS.95_phylo_contrasted[sp.evenness_TSS.95_phylo < 0.25] <- 0.255

# Add 0 values for empty continental pixels (transformed into min values as 0.25)
temp <- continent_mask
temp[!is.na(sp.evenness_Jaccard.80_phylo_contrasted[])] <- sp.evenness_Jaccard.80_phylo_contrasted[!is.na(sp.evenness_Jaccard.80_phylo_contrasted[])]
sp.evenness_Jaccard.80_phylo_contrasted <- temp
sp.evenness_Jaccard.80_phylo_contrasted[sp.evenness_Jaccard.80_phylo_contrasted == 0] <- 0.25

temp <- continent_mask
temp[!is.na(sp.evenness_Jaccard.95_phylo_contrasted[])] <- sp.evenness_Jaccard.95_phylo_contrasted[!is.na(sp.evenness_Jaccard.95_phylo_contrasted[])]
sp.evenness_Jaccard.95_phylo_contrasted <- temp
sp.evenness_Jaccard.95_phylo_contrasted[sp.evenness_Jaccard.95_phylo_contrasted == 0] <- 0.25

temp <- continent_mask
temp[!is.na(sp.evenness_TSS.80_phylo_contrasted[])] <- sp.evenness_TSS.80_phylo_contrasted[!is.na(sp.evenness_TSS.80_phylo_contrasted[])]
sp.evenness_TSS.80_phylo_contrasted <- temp
sp.evenness_TSS.80_phylo_contrasted[sp.evenness_TSS.80_phylo_contrasted == 0] <- 0.25

temp <- continent_mask
temp[!is.na(sp.evenness_TSS.95_phylo_contrasted[])] <- sp.evenness_TSS.95_phylo_contrasted[!is.na(sp.evenness_TSS.95_phylo_contrasted[])]
sp.evenness_TSS.95_phylo_contrasted <- temp
sp.evenness_TSS.95_phylo_contrasted[sp.evenness_TSS.95_phylo_contrasted == 0] <- 0.25

# hist(sp.evenness_Jaccard.80_phylo_contrasted)

### 4.3/ Plot with contrasted scale merging low values (0.25 = 0-0.25) ####

### Individual plots for contrasted scale

# Jaccard.80_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.80_phylo/sp.evenness_Jaccard.80_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.evenness_Jaccard.80_phylo_contrasted, col = pal_bl_red, main = paste0("Species evenness \nJaccard.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_Jaccard.80_phylo_contrasted, locs = seq(0.25, 0.4, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")
par(mar = internal_margins)
dev.off()

# Jaccard.95_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.80_phylo/sp.evenness_Jaccard.95_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.evenness_Jaccard.95_phylo_contrasted, col = pal_bl_red, main = paste0("Species evenness \nJaccard.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="         Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_Jaccard.95_phylo_contrasted, locs = seq(0.25, 0.45, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")
par(mar = internal_margins)
dev.off()

# TSS.80_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.80_phylo/sp.evenness_TSS.80_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.evenness_TSS.80_phylo_contrasted, col = pal_bl_red, main = paste0("Species evenness \nTSS.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="         Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_TSS.80_phylo_contrasted, locs = seq(0.25, 0.4, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")
par(mar = internal_margins)
dev.off()

# TSS.95_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.80_phylo/sp.evenness_TSS.95_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.evenness_TSS.95_phylo_contrasted, col = pal_bl_red, main = paste0("Species evenness \nTSS.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_TSS.95_phylo_contrasted, locs = seq(0.25, 0.45, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(sp.evenness_Jaccard.80_phylo) +
#   tm_raster(palette = pal_bl_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/phylo_only/sp.evenness_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.evenness_Jaccard.80_phylo_contrasted, col = pal_bl_red, main = paste0("Species evenness \nJaccard.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_Jaccard.80_phylo_contrasted, locs = seq(0.25, 0.4, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")

image(sp.evenness_Jaccard.95_phylo_contrasted, col = pal_bl_red, main = paste0("Species evenness \nJaccard.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="         Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_Jaccard.95_phylo_contrasted, locs = seq(0.25, 0.45, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")

image(sp.evenness_TSS.80_phylo_contrasted, col = pal_bl_red, main = paste0("Species evenness \nTSS.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="         Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_TSS.80_phylo_contrasted, locs = seq(0.25, 0.4, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Evenness")

image(sp.evenness_TSS.95_phylo_contrasted, col = pal_bl_red, main = paste0("Species evenness \nTSS.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evenness", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.evenness_TSS.95_phylo_contrasted, locs = seq(0.25, 0.45, 0.05), digits = 2, cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
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
All_sp_proba_stack_Jaccard.80_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_Jaccard.80_phylo.rds"))
All_sp_proba_stack_Jaccard.95_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_Jaccard.95_phylo.rds"))
All_sp_proba_stack_TSS.80_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_TSS.80_phylo.rds"))
All_sp_proba_stack_TSS.95_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_TSS.95_phylo.rds"))

# Load maps of Species richness to use to compute mean rarity index among local species
sp.richness_Jaccard.80_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_Jaccard.80_phylo.rds"))
sp.richness_Jaccard.95_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_Jaccard.95_phylo.rds"))
sp.richness_TSS.80_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_TSS.80_phylo.rds"))
sp.richness_TSS.95_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_TSS.95_phylo.rds"))

### 9.1/ Index computation ####

## For Jaccard.80_phylo

# Compute weights based on geographic ranges
sp.range_size_Jaccard.80_phylo <- NA
for (i in 1:nlayers(All_sp_proba_stack_Jaccard.80_phylo)) {
  # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
  # Compute the estimated number of pixels occupied by the species as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
  sp.range_size_Jaccard.80_phylo[i] <- sum(All_sp_proba_stack_Jaccard.80_phylo[[i]]@data@values, na.rm = T)*774.51/1000
}
sp.rarity_indices_Jaccard.80_phylo <- 1-(sp.range_size_Jaccard.80_phylo/max(sp.range_size_Jaccard.80_phylo)) # Weights by max species extent = rarity indices

# Apply weights
sp.weighted.stack_Jaccard.80_phylo <- All_sp_proba_stack_Jaccard.80_phylo # Generate the stack to fill
for (i in 1:nlayers(sp.weighted.stack_Jaccard.80_phylo)){
  # Multiply species probability of presence by ring rarity indices
  sp.weighted.stack_Jaccard.80_phylo@layers[[i]]@data@values <- All_sp_proba_stack_Jaccard.80_phylo@layers[[i]]@data@values * sp.rarity_indices_Jaccard.80_phylo[i]
}

# Sum = Richness weighted by rarity based on Range-size
sp.rarity_Jaccard.80_phylo <- readAll(calc(sp.weighted.stack_Jaccard.80_phylo, fun = sum)) 
# Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each species)
sp.mean.rarity_Jaccard.80_phylo <- sp.rarity_Jaccard.80_phylo/sp.richness_Jaccard.80_phylo

# Add 0 values for empty continental pixels (transformed into min values as 0.10)
temp <- continent_mask
temp[!is.na(sp.mean.rarity_Jaccard.80_phylo[])] <- sp.mean.rarity_Jaccard.80_phylo[!is.na(sp.mean.rarity_Jaccard.80_phylo[])]
sp.mean.rarity_Jaccard.80_phylo <- temp

# Save species ranges and species continous rarity indices
save(sp.range_size_Jaccard.80_phylo, file = "./outputs/Indices_maps/phylo_only/sp.range_size_Jaccard.80_phylo.RData", version = "2")
saveRDS(sp.range_size_Jaccard.80_phylo, file = "./outputs/Indices_maps/phylo_only/sp.range_size_Jaccard.80_phylo.rds", version = "2")
save(sp.rarity_indices_Jaccard.80_phylo, file = "./outputs/Indices_maps/phylo_only/sp.rarity_indices_Jaccard.80_phylo.RData", version = "2")
saveRDS(sp.rarity_indices_Jaccard.80_phylo, file = "./outputs/Indices_maps/phylo_only/sp.rarity_indices_Jaccard.80_phylo.rds", version = "2")

# Save full and mean rarity maps
save(sp.rarity_Jaccard.80_phylo, file = "./outputs/Indices_maps/phylo_only/sp.rarity_Jaccard.80_phylo.RData", version = "2")
saveRDS(sp.rarity_Jaccard.80_phylo, file = "./outputs/Indices_maps/phylo_only/sp.rarity_Jaccard.80_phylo.rds", version = "2")
save(sp.mean.rarity_Jaccard.80_phylo, file = "./outputs/Indices_maps/phylo_only/sp.mean.rarity_Jaccard.80_phylo.RData", version = "2")
saveRDS(sp.mean.rarity_Jaccard.80_phylo, file = "./outputs/Indices_maps/phylo_only/sp.mean.rarity_Jaccard.80_phylo.rds", version = "2")

# Plot distribution of ranges and rarity indices
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.80_phylo/sp_range_distri_Jaccard.80_phylo.pdf"), height = 5.3, width = 6.5)
hist(sp.range_size_Jaccard.80_phylo, main = "Distribution of species geographic ranges\nJaccard.80_phylo", breaks = 20, xlab = "Geographic ranges in [10^3] km²")
range_threshold_Jaccard.80_phylo <- round(quantile(sp.range_size_Jaccard.80_phylo, p=0.25), 1) # Range threshold value for qualitative rarity
abline(v = range_threshold_Jaccard.80_phylo, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",range_threshold_Jaccard.80_phylo," km²"), col = "red", lty = 2, lwd = 2, x = "topright", cex = 1, bty ="n")
dev.off()

pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.80_phylo/sp_rarity_indices_distri_Jaccard.80_phylo.pdf"), height = 5.3, width = 6.5)
hist(sp.rarity_indices_Jaccard.80_phylo, main = "Distribution of species rarity indices\nJaccard.80_phylo", breaks = 20, xlab = "Rarity indices")
indice_threshold_Jaccard.80_phylo <- round(quantile(sp.rarity_indices_Jaccard.80_phylo, p=0.75), 3) # Threshold value for qualitative rarity
abline(v = indice_threshold_Jaccard.80_phylo, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",indice_threshold_Jaccard.80_phylo), col = "red", lty = 2, lwd = 2, x = "topleft", cex = 1, bty ="n")
dev.off()


## For Jaccard.95_phylo

# Compute weights based on geographic ranges
sp.range_size_Jaccard.95_phylo <- NA
for (i in 1:nlayers(All_sp_proba_stack_Jaccard.95_phylo)) {
  # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
  # Compute the estimated number of pixels occupied by the species as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
  sp.range_size_Jaccard.95_phylo[i] <- sum(All_sp_proba_stack_Jaccard.95_phylo[[i]]@data@values, na.rm = T)*774.51/1000
}
sp.rarity_indices_Jaccard.95_phylo <- 1-(sp.range_size_Jaccard.95_phylo/max(sp.range_size_Jaccard.95_phylo)) # Weights by max species extent = rarity indices

# Apply weights
sp.weighted.stack_Jaccard.95_phylo <- All_sp_proba_stack_Jaccard.95_phylo # Generate the stack to fill
for (i in 1:nlayers(sp.weighted.stack_Jaccard.95_phylo)){
  # Multiply species probability of presence by species rarity indices
  sp.weighted.stack_Jaccard.95_phylo@layers[[i]]@data@values <- All_sp_proba_stack_Jaccard.95_phylo@layers[[i]]@data@values * sp.rarity_indices_Jaccard.95_phylo[i]
}

# plot(sp.weighted.stack_Jaccard.95_phylo)

# Sum = Richness weighted by rarity based on Range-size
sp.rarity_Jaccard.95_phylo <- readAll(calc(sp.weighted.stack_Jaccard.95_phylo, fun = sum)) 
# Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each species)
sp.mean.rarity_Jaccard.95_phylo <- sp.rarity_Jaccard.95_phylo/sp.richness_Jaccard.95_phylo

# Add 0 values for empty continental pixels (transformed into min values as 0.10)
temp <- continent_mask
temp[!is.na(sp.mean.rarity_Jaccard.95_phylo[])] <- sp.mean.rarity_Jaccard.95_phylo[!is.na(sp.mean.rarity_Jaccard.95_phylo[])]
sp.mean.rarity_Jaccard.95_phylo <- temp

# Save species ranges and species continuous rarity indices
save(sp.range_size_Jaccard.95_phylo, file = "./outputs/Indices_maps/phylo_only/sp.range_size_Jaccard.95_phylo.RData", version = "2")
saveRDS(sp.range_size_Jaccard.95_phylo, file = "./outputs/Indices_maps/phylo_only/sp.range_size_Jaccard.95_phylo.rds", version = "2")
save(sp.rarity_indices_Jaccard.95_phylo, file = "./outputs/Indices_maps/phylo_only/sp.rarity_indices_Jaccard.95_phylo.RData", version = "2")
saveRDS(sp.rarity_indices_Jaccard.95_phylo, file = "./outputs/Indices_maps/phylo_only/sp.rarity_indices_Jaccard.95_phylo.rds", version = "2")

# Save full and mean rarity maps
save(sp.rarity_Jaccard.95_phylo, file = "./outputs/Indices_maps/phylo_only/sp.rarity_Jaccard.95_phylo.RData", version = "2")
saveRDS(sp.rarity_Jaccard.95_phylo, file = "./outputs/Indices_maps/phylo_only/sp.rarity_Jaccard.95_phylo.rds", version = "2")
save(sp.mean.rarity_Jaccard.95_phylo, file = "./outputs/Indices_maps/phylo_only/sp.mean.rarity_Jaccard.95_phylo.RData", version = "2")
saveRDS(sp.mean.rarity_Jaccard.95_phylo, file = "./outputs/Indices_maps/phylo_only/sp.mean.rarity_Jaccard.95_phylo.rds", version = "2")

# Plot distribution of ranges and rarity indices
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.95_phylo/sp_range_distri_Jaccard.95_phylo.pdf"), height = 5.3, width = 6.5)
hist(sp.range_size_Jaccard.95_phylo, main = "Distribution of species geographic ranges\nJaccard.95_phylo", breaks = 20, xlab = "Geographic ranges in [10^3] km²")
range_threshold_Jaccard.95_phylo <- round(quantile(sp.range_size_Jaccard.95_phylo, p=0.25), 1) # Range threshold value for qualitative rarity
abline(v = range_threshold_Jaccard.95_phylo, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",range_threshold_Jaccard.95_phylo," km²"), col = "red", lty = 2, lwd = 2, x = "topright", cex = 1, bty ="n")
dev.off()

pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.95_phylo/sp_rarity_indices_distri_Jaccard.95_phylo.pdf"), height = 5.3, width = 6.5)
hist(sp.rarity_indices_Jaccard.95_phylo, main = "Distribution of species rarity indices\nJaccard.95_phylo", breaks = 20, xlab = "Rarity indices")
indice_threshold_Jaccard.95_phylo <- round(quantile(sp.rarity_indices_Jaccard.95_phylo, p=0.75), 3) # Threshold value for qualitative rarity
abline(v = indice_threshold_Jaccard.95_phylo, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",indice_threshold_Jaccard.95_phylo), col = "red", lty = 2, lwd = 2, x = "topleft", cex = 1, bty ="n")
dev.off()


## For TSS.80_phylo

# Compute weights based on geographic ranges
sp.range_size_TSS.80_phylo <- NA
for (i in 1:nlayers(All_sp_proba_stack_TSS.80_phylo)) {
  # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
  # Compute the estimated number of pixels occupied by the species as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
  sp.range_size_TSS.80_phylo[i] <- sum(All_sp_proba_stack_TSS.80_phylo[[i]]@data@values, na.rm = T)*774.51/1000
}
sp.rarity_indices_TSS.80_phylo <- 1-(sp.range_size_TSS.80_phylo/max(sp.range_size_TSS.80_phylo)) # Weights by max species extent = rarity indices

# Apply weights
sp.weighted.stack_TSS.80_phylo <- All_sp_proba_stack_TSS.80_phylo # Generate the stack to fill
for (i in 1:nlayers(sp.weighted.stack_TSS.80_phylo)){
  # Multiply species probability of presence by species rarity indices
  sp.weighted.stack_TSS.80_phylo@layers[[i]]@data@values <- All_sp_proba_stack_TSS.80_phylo@layers[[i]]@data@values * sp.rarity_indices_TSS.80_phylo[i]
}

plot(sp.weighted.stack_TSS.80_phylo)

# Sum = Richness weighted by rarity based on Range-size
sp.rarity_TSS.80_phylo <- readAll(calc(sp.weighted.stack_TSS.80_phylo, fun = sum)) 
# Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each species)
sp.mean.rarity_TSS.80_phylo <- sp.rarity_TSS.80_phylo/sp.richness_TSS.80_phylo

# Add 0 values for empty continental pixels (transformed into min values as 0.10)
temp <- continent_mask
temp[!is.na(sp.mean.rarity_TSS.80_phylo[])] <- sp.mean.rarity_TSS.80_phylo[!is.na(sp.mean.rarity_TSS.80_phylo[])]
sp.mean.rarity_TSS.80_phylo <- temp

# Save species ranges and species continous rarity indices
save(sp.range_size_TSS.80_phylo, file = "./outputs/Indices_maps/phylo_only/sp.range_size_TSS.80_phylo.RData", version = "2")
saveRDS(sp.range_size_TSS.80_phylo, file = "./outputs/Indices_maps/phylo_only/sp.range_size_TSS.80_phylo.rds", version = "2")
save(sp.rarity_indices_TSS.80_phylo, file = "./outputs/Indices_maps/phylo_only/sp.rarity_indices_TSS.80_phylo.RData", version = "2")
saveRDS(sp.rarity_indices_TSS.80_phylo, file = "./outputs/Indices_maps/phylo_only/sp.rarity_indices_TSS.80_phylo.rds", version = "2")

# Save full and mean rarity maps
save(sp.rarity_TSS.80_phylo, file = "./outputs/Indices_maps/phylo_only/sp.rarity_TSS.80_phylo.RData", version = "2")
saveRDS(sp.rarity_TSS.80_phylo, file = "./outputs/Indices_maps/phylo_only/sp.rarity_TSS.80_phylo.rds", version = "2")
save(sp.mean.rarity_TSS.80_phylo, file = "./outputs/Indices_maps/phylo_only/sp.mean.rarity_TSS.80_phylo.RData", version = "2")
saveRDS(sp.mean.rarity_TSS.80_phylo, file = "./outputs/Indices_maps/phylo_only/sp.mean.rarity_TSS.80_phylo.rds", version = "2")

# Plot distribution of ranges and rarity indices
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.80_phylo/sp_range_distri_TSS.80_phylo.pdf"), height = 5.3, width = 6.5)
hist(sp.range_size_TSS.80_phylo, main = "Distribution of species geographic ranges\nTSS.80_phylo", breaks = 20, xlab = "Geographic ranges in [10^3] km²")
range_threshold_TSS.80_phylo <- round(quantile(sp.range_size_TSS.80_phylo, p=0.25), 1) # Range threshold value for qualitative rarity
abline(v = range_threshold_TSS.80_phylo, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",range_threshold_TSS.80_phylo," km²"), col = "red", lty = 2, lwd = 2, x = "topright", cex = 1, bty ="n")
dev.off()

pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.80_phylo/sp_rarity_indices_distri_TSS.80_phylo.pdf"), height = 5.3, width = 6.5)
hist(sp.rarity_indices_TSS.80_phylo, main = "Distribution of species rarity indices\nTSS.80_phylo", breaks = 20, xlab = "Rarity indices")
indice_threshold_TSS.80_phylo <- round(quantile(sp.rarity_indices_TSS.80_phylo, p=0.75), 3) # Threshold value for qualitative rarity
abline(v = indice_threshold_TSS.80_phylo, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",indice_threshold_TSS.80_phylo), col = "red", lty = 2, lwd = 2, x = "topleft", cex = 1, bty ="n")
dev.off()


## For TSS.95_phylo

# Compute weights based on geographic ranges
sp.range_size_TSS.95_phylo <- NA
for (i in 1:nlayers(All_sp_proba_stack_TSS.95_phylo)) {
  # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
  # Compute the estimated number of pixels occupied by the species as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
  sp.range_size_TSS.95_phylo[i] <- sum(All_sp_proba_stack_TSS.95_phylo[[i]]@data@values, na.rm = T)*774.51/1000
}
sp.rarity_indices_TSS.95_phylo <- 1-(sp.range_size_TSS.95_phylo/max(sp.range_size_TSS.95_phylo)) # Weights by max species extent = rarity indices

# Apply weights
sp.weighted.stack_TSS.95_phylo <- All_sp_proba_stack_TSS.95_phylo # Generate the stack to fill
for (i in 1:nlayers(sp.weighted.stack_TSS.95_phylo)){
  # Multiply species probability of presence by species rarity indices
  sp.weighted.stack_TSS.95_phylo@layers[[i]]@data@values <- All_sp_proba_stack_TSS.95_phylo@layers[[i]]@data@values * sp.rarity_indices_TSS.95_phylo[i]
}

plot(sp.weighted.stack_TSS.95_phylo)

# Sum = Richness weighted by rarity based on Range-size
sp.rarity_TSS.95_phylo <- readAll(calc(sp.weighted.stack_TSS.95_phylo, fun = sum)) 
# Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each species)
sp.mean.rarity_TSS.95_phylo <- sp.rarity_TSS.95_phylo/sp.richness_TSS.95_phylo

# Add 0 values for empty continental pixels (transformed into min values as 0.10)
temp <- continent_mask
temp[!is.na(sp.mean.rarity_TSS.95_phylo[])] <- sp.mean.rarity_TSS.95_phylo[!is.na(sp.mean.rarity_TSS.95_phylo[])]
sp.mean.rarity_TSS.95_phylo <- temp

# Save species ranges and species continuous rarity indices
save(sp.range_size_TSS.95_phylo, file = "./outputs/Indices_maps/phylo_only/sp.range_size_TSS.95_phylo.RData", version = "2")
saveRDS(sp.range_size_TSS.95_phylo, file = "./outputs/Indices_maps/phylo_only/sp.range_size_TSS.95_phylo.rds", version = "2")
save(sp.rarity_indices_TSS.95_phylo, file = "./outputs/Indices_maps/phylo_only/sp.rarity_indices_TSS.95_phylo.RData", version = "2")
saveRDS(sp.rarity_indices_TSS.95_phylo, file = "./outputs/Indices_maps/phylo_only/sp.rarity_indices_TSS.95_phylo.rds", version = "2")

# Save full and mean rarity maps
save(sp.rarity_TSS.95_phylo, file = "./outputs/Indices_maps/phylo_only/sp.rarity_TSS.95_phylo.RData", version = "2")
saveRDS(sp.rarity_TSS.95_phylo, file = "./outputs/Indices_maps/phylo_only/sp.rarity_TSS.95_phylo.rds", version = "2")
save(sp.mean.rarity_TSS.95_phylo, file = "./outputs/Indices_maps/phylo_only/sp.mean.rarity_TSS.95_phylo.RData", version = "2")
saveRDS(sp.mean.rarity_TSS.95_phylo, file = "./outputs/Indices_maps/phylo_only/sp.mean.rarity_TSS.95_phylo.rds", version = "2")

# Plot distribution of ranges and rarity indices
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.95_phylo/sp_range_distri_TSS.95_phylo.pdf"), height = 5.3, width = 6.5)
hist(sp.range_size_TSS.95_phylo, main = "Distribution of species geographic ranges\nTSS.95_phylo", breaks = 20, xlab = "Geographic ranges in [10^3] km²")
range_threshold_TSS.95_phylo <- round(quantile(sp.range_size_TSS.95_phylo, p=0.25), 1) # Range threshold value for qualitative rarity
abline(v = range_threshold_TSS.95_phylo, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",range_threshold_TSS.95_phylo," km²"), col = "red", lty = 2, lwd = 2, x = "topright", cex = 1, bty ="n")
dev.off()

pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.95_phylo/sp_rarity_indices_distri_TSS.95_phylo.pdf"), height = 5.3, width = 6.5)
hist(sp.rarity_indices_TSS.95_phylo, main = "Distribution of species rarity indices\nTSS.95_phylo", breaks = 20, xlab = "Rarity indices")
indice_threshold_TSS.95_phylo <- round(quantile(sp.rarity_indices_TSS.95_phylo, p=0.75), 3) # Threshold value for qualitative rarity
abline(v = indice_threshold_TSS.95_phylo, col = "red", lty = 2, lwd = 2)
legend(legend = paste0("Threshold 75% for rarity\n",indice_threshold_TSS.95_phylo), col = "red", lty = 2, lwd = 2, x = "topleft", cex = 1, bty ="n")
dev.off()

### 9.2/ Plot community rarity maps ####

### Load directly the community rarity maps
sp.rarity_Jaccard.80_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/sp.rarity_Jaccard.80_phylo.rds")
sp.rarity_Jaccard.95_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/sp.rarity_Jaccard.95_phylo.rds")
sp.rarity_TSS.80_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/sp.rarity_TSS.80_phylo.rds")
sp.rarity_TSS.95_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/sp.rarity_TSS.95_phylo.rds")


### Individual plots

# Jaccard.80_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.80_phylo/sp.rarity_Jaccard.80_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.rarity_Jaccard.80_phylo, col = pal_bl_red, main = paste0("Species rarity \nJaccard.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_Jaccard.80_phylo, locs = seq(0, 100, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")
par(mar = internal_margins)
dev.off()

# Jaccard.95_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.95_phylo/sp.rarity_Jaccard.95_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.rarity_Jaccard.95_phylo, col = pal_bl_red, main = paste0("Species rarity \nJaccard.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_Jaccard.95_phylo, locs = seq(0, 100, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")
par(mar = internal_margins)
dev.off()


# TSS.80_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.80_phylo/sp.rarity_TSS.80_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.rarity_TSS.80_phylo, col = pal_bl_red, main = paste0("Species rarity \nTSS.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_TSS.80_phylo, locs = seq(0, 100, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")
par(mar = internal_margins)
dev.off()


# TSS.95_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.95_phylo/sp.rarity_TSS.95_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.rarity_TSS.95_phylo, col = pal_bl_red, main = paste0("Species rarity \nTSS.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_TSS.95_phylo, locs = seq(0, 100, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(sp.rarity_Jaccard.80_phylo) +
#   tm_raster(palette = pal_bl_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/phylo_only/sp.rarity_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.rarity_Jaccard.80_phylo, col = pal_bl_red, main = paste0("Species rarity \nJaccard.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_Jaccard.80_phylo, locs = seq(0, 100, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")

image(sp.rarity_Jaccard.95_phylo, col = pal_bl_red, main = paste0("Species rarity \nJaccard.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_Jaccard.95_phylo, locs = seq(0, 100, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")

image(sp.rarity_TSS.80_phylo, col = pal_bl_red, main = paste0("Species rarity \nTSS.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_TSS.80_phylo, locs = seq(0, 100, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")

image(sp.rarity_TSS.95_phylo, col = pal_bl_red, main = paste0("Species rarity \nTSS.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Community\n         rarity", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.rarity_TSS.95_phylo, locs = seq(0, 100, 20), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 5.5, font = 2, cex = 1.2, label = "Community\nrarity")

par(mar = internal_margins)

dev.off()


### 9.3/ Plot mean species rarity maps ####

### Load directly the mean species rarity maps
sp.mean.rarity_Jaccard.80_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/sp.mean.rarity_Jaccard.80_phylo.rds")
sp.mean.rarity_Jaccard.95_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/sp.mean.rarity_Jaccard.95_phylo.rds")
sp.mean.rarity_TSS.80_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/sp.mean.rarity_TSS.80_phylo.rds")
sp.mean.rarity_TSS.95_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/sp.mean.rarity_TSS.95_phylo.rds")


### Individual plots

# Jaccard.80_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.80_phylo/sp.mean.rarity_Jaccard.80_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.mean.rarity_Jaccard.80_phylo, col = pal_bl_red, main = paste0("Mean species rarity\nJaccard.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_Jaccard.80_phylo, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")
par(mar = internal_margins)
dev.off()

# Jaccard.95_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.95_phylo/sp.mean.rarity_Jaccard.95_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.mean.rarity_Jaccard.95_phylo, col = pal_bl_red, main = paste0("Mean species rarity\nJaccard.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_Jaccard.95_phylo, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")
par(mar = internal_margins)
dev.off()


# TSS.80_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.80_phylo/sp.mean.rarity_TSS.80_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.mean.rarity_TSS.80_phylo, col = pal_bl_red, main = paste0("Mean species rarity\nTSS.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_TSS.80_phylo, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")
par(mar = internal_margins)
dev.off()

# TSS.95_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.95_phylo/sp.mean.rarity_TSS.95_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.mean.rarity_TSS.95_phylo, col = pal_bl_red, main = paste0("Mean species rarity \nTSS.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_TSS.95_phylo, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(sp.mean.rarity_Jaccard.80_phylo) +
#   tm_raster(palette = pal_bl_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/phylo_only/sp.mean.rarity_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.mean.rarity_Jaccard.80_phylo, col = pal_bl_red, main = paste0("Mean species rarity \nJaccard.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_Jaccard.80_phylo, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")

image(sp.mean.rarity_Jaccard.95_phylo, col = pal_bl_red, main = paste0("Mean species rarity\nJaccard.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_Jaccard.95_phylo, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")

image(sp.mean.rarity_TSS.80_phylo, col = pal_bl_red, main = paste0("Mean species rarity\nTSS.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_TSS.80_phylo, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")

image(sp.mean.rarity_TSS.95_phylo, col = pal_bl_red, main = paste0("Mean species rarity\nTSS.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(sp.mean.rarity_TSS.95_phylo, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.2, label = "Rarity\nindex")

par(mar = internal_margins)

dev.off()

# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])



###### 10/ Categorical Species Rarity (25% threshold) ######

### Load directly the stack of sp probas use to compute rarity based on geographic range
All_sp_proba_stack_Jaccard.80_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_Jaccard.80_phylo.rds"))
All_sp_proba_stack_Jaccard.95_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_Jaccard.95_phylo.rds"))
All_sp_proba_stack_TSS.80_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_TSS.80_phylo.rds"))
All_sp_proba_stack_TSS.95_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_TSS.95_phylo.rds"))

### Load species ranges
sp.range_size_Jaccard.80_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/sp.range_size_Jaccard.80_phylo.rds")
sp.range_size_Jaccard.95_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/sp.range_size_Jaccard.95_phylo.rds")
sp.range_size_TSS.80_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/sp.range_size_TSS.80_phylo.rds")
sp.range_size_TSS.95_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/sp.range_size_TSS.95_phylo.rds")

### 10.1/ Extract rarity threshold ####
rarity_threshold_Jaccard.80_phylo <- round(quantile(sp.range_size_Jaccard.80_phylo, probs = 0.25), 3) ; rarity_threshold_Jaccard.80_phylo
# Threshold surface to be considered a rare species = 136 845 km² 
rarity_threshold_Jaccard.95_phylo <- quantile(sp.range_size_Jaccard.95_phylo, probs = 0.25) ; rarity_threshold_Jaccard.95_phylo
# Threshold surface to be considered a rare species = 155 610 km² 
rarity_threshold_TSS.80_phylo <- quantile(sp.range_size_TSS.80_phylo, probs = 0.25) ; rarity_threshold_TSS.80_phylo
# Threshold surface to be considered a rare species = 136 845 km² 
rarity_threshold_TSS.95_phylo <- quantile(sp.range_size_TSS.95_phylo, probs = 0.25) ; rarity_threshold_TSS.95_phylo
# Threshold surface to be considered a rare species = 156 065 km² 

### 10.2/ Plot species range distribution ####

# Plot species range distribution for Jaccard.80_phylo
ordered_sp.range_size_Jaccard.80_phylo <- sp.range_size_Jaccard.80_phylo[order(sp.range_size_Jaccard.80_phylo)]

pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.80_phylo/sp_range_distri_barplot_Jaccard.80_phylo.pdf"), height = 5.3, width = 6.5)
barplot(ordered_sp.range_size_Jaccard.80_phylo, ylim = c(0,9000), ylab = "Surface in 10^3 km", xlab = "Species", main = "Distribution of Species range size")
abline(v = length(ordered_sp.range_size_Jaccard.80_phylo)/4, col = "red", lwd = 2, lty = 2)
legend(legend = c("Threshold 25 %", paste0("    ",rarity_threshold_Jaccard.80_phylo*1000," km²")), bty = "n", x = "top", text.col = "red", text.font = 2)
legend(legend = c("  RARE", "Species"), bty = "n", x = "left", text.col = "darkblue", text.font = 2)
dev.off()

# Plot species range distribution for Jaccard.95_phylo
ordered_sp.range_size_Jaccard.95_phylo <- sp.range_size_Jaccard.95_phylo[order(sp.range_size_Jaccard.95_phylo)]

pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.95_phylo/sp_range_distri_barplot_Jaccard.95_phylo.pdf"), height = 5.3, width = 6.5)
barplot(ordered_sp.range_size_Jaccard.95_phylo, ylim = c(0,9000), ylab = "Surface in 10^3 km", xlab = "Species", main = "Distribution of Species range size")
abline(v = length(ordered_sp.range_size_Jaccard.95_phylo)/4, col = "red", lwd = 2, lty = 2)
legend(legend = c("Threshold 25 %", paste0("    ",rarity_threshold_Jaccard.95_phylo*1000," km²")), bty = "n", x = "top", text.col = "red", text.font = 2)
legend(legend = c("  RARE", "Species"), bty = "n", x = "left", text.col = "darkblue", text.font = 2)
dev.off()

# Plot species range distribution for TSS.80_phylo
ordered_sp.range_size_TSS.80_phylo <- sp.range_size_TSS.80_phylo[order(sp.range_size_TSS.80_phylo)]

pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.80_phylo/sp_range_distri_barplot_TSS.80_phylo.pdf"), height = 5.3, width = 6.5)
barplot(ordered_sp.range_size_TSS.80_phylo, ylim = c(0,9000), ylab = "Surface in 10^3 km", xlab = "Species", main = "Distribution of Species range size")
abline(v = length(ordered_sp.range_size_TSS.80_phylo)/4, col = "red", lwd = 2, lty = 2)
legend(legend = c("Threshold 25 %", paste0("    ",rarity_threshold_TSS.80_phylo*1000," km²")), bty = "n", x = "top", text.col = "red", text.font = 2)
legend(legend = c("  RARE", "Species"), bty = "n", x = "left", text.col = "darkblue", text.font = 2)
dev.off()

# Plot species range distribution for TSS.95_phylo
ordered_sp.range_size_TSS.95_phylo <- sp.range_size_TSS.95_phylo[order(sp.range_size_TSS.95_phylo)]

pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.95_phylo/sp_range_distri_barplot_TSS.95_phylo.pdf"), height = 5.3, width = 6.5)
barplot(ordered_sp.range_size_TSS.95_phylo, ylim = c(0,9000), ylab = "Surface in 10^3 km", xlab = "Species", main = "Distribution of Species range size")
abline(v = length(ordered_sp.range_size_TSS.95_phylo)/4, col = "red", lwd = 2, lty = 2)
legend(legend = c("Threshold 25 %", paste0("    ",rarity_threshold_TSS.95_phylo*1000," km²")), bty = "n", x = "top", text.col = "red", text.font = 2)
legend(legend = c("  RARE", "Species"), bty = "n", x = "left", text.col = "darkblue", text.font = 2)
dev.off()

### 10.3/ Extract only rare species to compute categorical rarity indices ####

# Detect rare species
rare.sp_Jaccard.80_phylo <- sp.range_size_Jaccard.80_phylo < rarity_threshold_Jaccard.80_phylo
rare.sp_Jaccard.95_phylo <- sp.range_size_Jaccard.95_phylo < rarity_threshold_Jaccard.95_phylo
rare.sp_TSS.80_phylo <- sp.range_size_TSS.80_phylo < rarity_threshold_TSS.80_phylo
rare.sp_TSS.95_phylo <- sp.range_size_TSS.95_phylo < rarity_threshold_TSS.95_phylo

sum(rare.sp_Jaccard.80_phylo) # 85 rare species on 339 (25%)

# Extract only rare species
Rare_sp_proba_stack_Jaccard.80_phylo <- All_sp_proba_stack_Jaccard.80_phylo[[which(rare.sp_Jaccard.80_phylo)]]
Rare_sp_proba_stack_Jaccard.95_phylo <- All_sp_proba_stack_Jaccard.95_phylo[[which(rare.sp_Jaccard.95_phylo)]]
Rare_sp_proba_stack_TSS.80_phylo <- All_sp_proba_stack_TSS.80_phylo[[which(rare.sp_TSS.80_phylo)]]
Rare_sp_proba_stack_TSS.95_phylo <- All_sp_proba_stack_TSS.95_phylo[[which(rare.sp_TSS.95_phylo)]]

### 10.4/ Index computation ####

# Load maps of Species richness to use to compute mean rarity index among local species
sp.richness_Jaccard.80_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_Jaccard.80_phylo.rds"))
sp.richness_Jaccard.95_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_Jaccard.95_phylo.rds"))
sp.richness_TSS.80_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_TSS.80_phylo.rds"))
sp.richness_TSS.95_phylo <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/tot.sp.richness_TSS.95_phylo.rds"))

# Index computation for total of rare species and proportion of rare species
tot.rare.sp_Jaccard.80_phylo <- readAll(calc(Rare_sp_proba_stack_Jaccard.80_phylo, fun = sum))
prop.rare.sp_Jaccard.80_phylo <- tot.rare.sp_Jaccard.80_phylo/sp.richness_Jaccard.80_phylo*100
tot.rare.sp_Jaccard.95_phylo <- readAll(calc(Rare_sp_proba_stack_Jaccard.95_phylo, fun = sum))
prop.rare.sp_Jaccard.95_phylo <- tot.rare.sp_Jaccard.95_phylo/sp.richness_Jaccard.95_phylo*100
tot.rare.sp_TSS.80_phylo <- readAll(calc(Rare_sp_proba_stack_TSS.80_phylo, fun = sum))
prop.rare.sp_TSS.80_phylo <- tot.rare.sp_TSS.80_phylo/sp.richness_TSS.80_phylo*100
tot.rare.sp_TSS.95_phylo <- readAll(calc(Rare_sp_proba_stack_TSS.95_phylo, fun = sum))
prop.rare.sp_TSS.95_phylo <- tot.rare.sp_TSS.95_phylo/sp.richness_TSS.95_phylo*100

# Save maps
save(tot.rare.sp_Jaccard.80_phylo , file = "./outputs/Indices_maps/phylo_only/tot.rare.sp_Jaccard.80_phylo.RData", version = "2")
saveRDS(tot.rare.sp_Jaccard.80_phylo , file = "./outputs/Indices_maps/phylo_only/tot.rare.sp_Jaccard.80_phylo.rds", version = "2")
save(tot.rare.sp_Jaccard.95_phylo , file = "./outputs/Indices_maps/phylo_only/tot.rare.sp_Jaccard.95_phylo.RData", version = "2")
saveRDS(tot.rare.sp_Jaccard.95_phylo , file = "./outputs/Indices_maps/phylo_only/tot.rare.sp_Jaccard.95_phylo.rds", version = "2")
save(tot.rare.sp_TSS.80_phylo , file = "./outputs/Indices_maps/phylo_only/tot.rare.sp_TSS.80_phylo.RData", version = "2")
saveRDS(tot.rare.sp_TSS.80_phylo , file = "./outputs/Indices_maps/phylo_only/tot.rare.sp_TSS.80_phylo.rds", version = "2")
save(tot.rare.sp_TSS.95_phylo , file = "./outputs/Indices_maps/phylo_only/tot.rare.sp_TSS.95_phylo.RData", version = "2")
saveRDS(tot.rare.sp_TSS.95_phylo , file = "./outputs/Indices_maps/phylo_only/tot.rare.sp_TSS.95_phylo.rds", version = "2")

save(prop.rare.sp_Jaccard.80_phylo , file = "./outputs/Indices_maps/phylo_only/prop.rare.sp_Jaccard.80_phylo.RData", version = "2")
saveRDS(prop.rare.sp_Jaccard.80_phylo , file = "./outputs/Indices_maps/phylo_only/prop.rare.sp_Jaccard.80_phylo.rds", version = "2")
save(prop.rare.sp_Jaccard.95_phylo , file = "./outputs/Indices_maps/phylo_only/prop.rare.sp_Jaccard.95_phylo.RData", version = "2")
saveRDS(prop.rare.sp_Jaccard.95_phylo , file = "./outputs/Indices_maps/phylo_only/prop.rare.sp_Jaccard.95_phylo.rds", version = "2")
save(prop.rare.sp_TSS.80_phylo , file = "./outputs/Indices_maps/phylo_only/prop.rare.sp_TSS.80_phylo.RData", version = "2")
saveRDS(prop.rare.sp_TSS.80_phylo , file = "./outputs/Indices_maps/phylo_only/prop.rare.sp_TSS.80_phylo.rds", version = "2")
save(prop.rare.sp_TSS.95_phylo , file = "./outputs/Indices_maps/phylo_only/prop.rare.sp_TSS.95_phylo.RData", version = "2")
saveRDS(prop.rare.sp_TSS.95_phylo , file = "./outputs/Indices_maps/phylo_only/prop.rare.sp_TSS.95_phylo.rds", version = "2")

### 10.5/ Plots of rare species richness ####

### Load directly the final Rare Species richness layer
tot.rare.sp_Jaccard.80_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/tot.rare.sp_Jaccard.80_phylo.rds")
tot.rare.sp_Jaccard.95_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/tot.rare.sp_Jaccard.95_phylo.rds")
tot.rare.sp_TSS.80_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/tot.rare.sp_TSS.80_phylo.rds")
tot.rare.sp_TSS.95_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/tot.rare.sp_TSS.95_phylo.rds")


### Individual plots

# Jaccard.80_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.80_phylo/tot.rare.sp_Jaccard.80_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.rare.sp_Jaccard.80_phylo, col = pal_bl_red, main = paste0("Rare species richness \nJaccard.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_Jaccard.80_phylo, locs = seq(0, 15, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Rare\nspecies")
par(mar = internal_margins)
dev.off()

# Jaccard.95_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.95_phylo/tot.rare.sp_Jaccard.95_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.rare.sp_Jaccard.95_phylo, col = pal_bl_red, main = paste0("Rare species richness \nJaccard.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_Jaccard.95_phylo, locs = seq(0, 15, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Rare\nspecies")
par(mar = internal_margins)
dev.off()

# TSS.80_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.80_phylo/tot.rare.sp_TSS.80_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.rare.sp_TSS.80_phylo, col = pal_bl_red, main = paste0("Rare species richness \nTSS.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_TSS.80_phylo, locs = seq(0, 15, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Rare\nspecies")
par(mar = internal_margins)
dev.off()

# TSS.95_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.95_phylo/tot.rare.sp_TSS.95_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.rare.sp_TSS.95_phylo, col = pal_bl_red, main = paste0("Rare species richness \nTSS.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_TSS.95_phylo, locs = seq(0, 15, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Rare\nspecies")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(tot.sp.richness_Jaccard.80_phylo) +
#   tm_raster(palette = pal_bl_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/phylo_only/tot.rare.sp_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(tot.rare.sp_Jaccard.80_phylo, col = pal_bl_red, main = paste0("Rare species richness \nJaccard.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_Jaccard.80_phylo, locs = seq(0, 15, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Rare\nspecies")

image(tot.rare.sp_Jaccard.95_phylo, col = pal_bl_red, main = paste0("Rare species richness \nJaccard.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_Jaccard.95_phylo, locs = seq(0, 15, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Rare\nspecies")

image(tot.rare.sp_TSS.80_phylo, col = pal_bl_red, main = paste0("Rare species richness \nTSS.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_TSS.80_phylo, locs = seq(0, 15, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Rare\nspecies")

image(tot.rare.sp_TSS.95_phylo, col = pal_bl_red, main = paste0("Rare species richness \nTSS.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.rare.sp_TSS.95_phylo, locs = seq(0, 15, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.2, label = "Rare\nspecies")

par(mar = internal_margins)

dev.off()


### 10.5 /Plots of proportion of rare species ####

### Load directly the final Rare Species richness layer
prop.rare.sp_Jaccard.80_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/prop.rare.sp_Jaccard.80_phylo.rds")
prop.rare.sp_Jaccard.95_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/prop.rare.sp_Jaccard.95_phylo.rds")
prop.rare.sp_TSS.80_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/prop.rare.sp_TSS.80_phylo.rds")
prop.rare.sp_TSS.95_phylo <- readRDS(file = "./outputs/Indices_maps/phylo_only/prop.rare.sp_TSS.95_phylo.rds")

hist(prop.rare.sp_Jaccard.80_phylo)

# Add 0 values for empty continental pixels (transformed into min values as 0.25)
temp <- continent_mask
temp[!is.na(prop.rare.sp_Jaccard.80_phylo[])] <- prop.rare.sp_Jaccard.80_phylo[!is.na(prop.rare.sp_Jaccard.80_phylo)]
prop.rare.sp_Jaccard.80_phylo <- temp

temp <- continent_mask
temp[!is.na(prop.rare.sp_Jaccard.95_phylo[])] <- prop.rare.sp_Jaccard.95_phylo[!is.na(prop.rare.sp_Jaccard.95_phylo)]
prop.rare.sp_Jaccard.95_phylo <- temp

temp <- continent_mask
temp[!is.na(prop.rare.sp_TSS.80_phylo[])] <- prop.rare.sp_TSS.80_phylo[!is.na(prop.rare.sp_TSS.80_phylo)]
prop.rare.sp_TSS.80_phylo <- temp

temp <- continent_mask
temp[!is.na(prop.rare.sp_TSS.95_phylo[])] <- prop.rare.sp_TSS.95_phylo[!is.na(prop.rare.sp_TSS.95_phylo)]
prop.rare.sp_TSS.95_phylo <- temp

# Contrast by merging all values higher than 25 %
prop.rare.sp_Jaccard.80_phylo[prop.rare.sp_Jaccard.80_phylo > 25] <-  25
prop.rare.sp_Jaccard.95_phylo[prop.rare.sp_Jaccard.95_phylo > 25] <-  25
prop.rare.sp_TSS.80_phylo[prop.rare.sp_TSS.80_phylo > 25] <-  25
prop.rare.sp_TSS.95_phylo[prop.rare.sp_TSS.95_phylo > 25] <-  25


### Individual plots

# Jaccard.80_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.80_phylo/prop.rare.sp_Jaccard.80_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(prop.rare.sp_Jaccard.80_phylo, col = pal_bl_red, main = paste0("Proportion of rare species\nJaccard.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_Jaccard.80_phylo, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")
par(mar = internal_margins)
dev.off()

# Jaccard.95_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/Jaccard.95_phylo/prop.rare.sp_Jaccard.95_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(prop.rare.sp_Jaccard.95_phylo, col = pal_bl_red, main = paste0("Proportion of rare species\nJaccard.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_Jaccard.95_phylo, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")
par(mar = internal_margins)
dev.off()

# TSS.80_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.80_phylo/prop.rare.sp_TSS.80_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(prop.rare.sp_TSS.80_phylo, col = pal_bl_red, main = paste0("Proportion of rare species\nTSS.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_TSS.80_phylo, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")
par(mar = internal_margins)
dev.off()

# TSS.95_phylo
pdf(file = paste0("./maps/Indices_maps/phylo_only/TSS.95_phylo/prop.rare.sp_TSS.95_phylo.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(prop.rare.sp_TSS.95_phylo, col = pal_bl_red, main = paste0("Proportion of rare species\nTSS.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_TSS.95_phylo, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(prop.sp.richness_Jaccard.80_phylo) +
#   tm_raster(palette = pal_bl_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Indices_maps/phylo_only/prop.rare.sp_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(prop.rare.sp_Jaccard.80_phylo, col = pal_bl_red, main = paste0("Proportion of rare species\nJaccard.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_Jaccard.80_phylo, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")

image(prop.rare.sp_Jaccard.95_phylo, col = pal_bl_red, main = paste0("Proportion of rare species\nJaccard.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_Jaccard.95_phylo, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")

image(prop.rare.sp_TSS.80_phylo, col = pal_bl_red, main = paste0("Proportion of rare species\nTSS.80_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_TSS.80_phylo, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")

image(prop.rare.sp_TSS.95_phylo, col = pal_bl_red, main = paste0("Proportion of rare species\nTSS.95_phylo"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="       Rare\n         species\n       %", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(prop.rare.sp_TSS.95_phylo, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 9, font = 2, cex = 1.2, label = "Rare\nspecies\n(%)")

par(mar = internal_margins)

dev.off()

# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])



###### 11/ Continuous Mimicry rarity = Range size weighted mimicry richness ######

# Mimicry-level index. Not very useful to control for absence of species not in the phylogeny


###### 12/ Community vulnerability ######

# Computed on the basis of mimicry ring richness. Not very useful to control for absence of species not in the phylogeny


