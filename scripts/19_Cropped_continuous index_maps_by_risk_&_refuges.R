
###### Script 19: Cropped index maps based on HF categories ######

# Crop continuous index maps using categories of threats based on HF levels

# Only for Jaccard.80

### Use only 4 indices

# 1/ Species richness
# 2/ Mean Continuous Species Rarity = Range-size Weighted Species Richness
# 3/ Mimicry richness
# 4/ Mean Continuous Mimicry rarity = Range-size Weighted Mimicry Richness
# Bonus/ Mean pairwise Phylogenetic Distance (removed from the final analyses)

###


# Inputs: 
   # Indices maps from script 14a
   # Threat map: HII (HF, Venter et al., 2016)

# Outputs:
   # Mask and shp files for HF high and low level values (areas by categories)
   # Maps of biodiversity indices cropped by risk and refuges based on HF
       # Per quantiles 5 & 25%
       # Per absolute values (1, 5 vs. 11, 23)

   # Boxplots of indices values per risk/refuges areas
      # Per quantiles 5 & 25%
      # Per absolute values (1, 5 vs. 11, 23)

###



### 1/ Load stuff ####

# Effacer l'environnement
rm(list = ls())

library(raster)
library(rgeos)
library(rangeBuilder)

res <- "15"

# Load mask for continent borders
continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_", res, ".rds"))
crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".rds"))
bg_mask <- readRDS(file = "./input_data/Map_stuff/bg_mask.rds")
bg_mask_pixel <- readRDS(file = "./input_data/Map_stuff/bg_mask_pixel.rds")

# Load indices map
sp.richness <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
# sp.rarity <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_Jaccard.80.rds")        # Linear weighting
sp.rarity <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_Leroy_Jaccard.80.rds")  # Leroy's weighting
MPD <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80.rds")
MPD_contrasted <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80_contrasted.rds")
ring.richness <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))
# ring.rarity <- readRDS(file = paste0("./outputs/Indices_maps/mimicry.mean.rarity_Jaccard.80.rds"))          # Linear weighting
ring.rarity <- readRDS(file = "./outputs/Indices_maps/ring.mean.rarity_Leroy_Jaccard.80.rds")        # Leroy's weighting


# Load shp files for maps
load(file = "./input_data/Map_stuff/bg_mask.RData") # Load bg shp
load(file = "./input_data/Map_stuff/country_borders.RData") # Load country borders

# Load color palette
library(tmaptools)
pal_grn_red_NA <- pal_grn_red <- rev(get_brewer_pal("RdYlGn", n = 200))
pal_grn_red_NA[1] <- "#EDEDED"
pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")

# Load HII
HII <- readRDS(file = "./input_data/HII/HII.rds")

# Create a mask of communities with presence of Ithomiini
Community_mask <- sp.richness
Community_mask[!(Community_mask > 0)] <- NA

# Extract number of community with Ithomiini
N.com <- sum(!is.na(Community_mask[]))

# Extract nb of communities per quantiles 5% and 25%
N.com_5 <- round(N.com *5/100, 0)
N.com_25 <- round(N.com *25/100, 0)

# Different nb of community for MPD because cannot be computed from community predicted with less than 1 species
N.com.MPD <- length(MPD[MPD > 0])


### 2/ Retrieve indices of most/least threatened communities and generate shp files of areas ####

### 2.1/ For quantiles ####

# Extract HII values for communities with Ithomiini
HII_masked <- mask(HII, Community_mask)

# Extract indices of communities with Ithomiini ordered by decreasing value of HII
HII_indices <- order(HII_masked[], decreasing = T)
HII_indices <- HII_indices[1:N.com]

HII_indices_5 <- HII_indices[1:N.com_5]
HII_indices_25 <- HII_indices[1:N.com_25]
HII_indices_75 <- HII_indices[(N.com-N.com_25):N.com]
HII_indices_95 <- HII_indices[(N.com-N.com_5):N.com]
HII_indices_25_75 <- c(HII_indices_25, HII_indices_75)

# Save indices
save(HII_indices_5, HII_indices_25, HII_indices_75, HII_indices_95, HII_indices_25_75, file = "./outputs/Threat_maps/HII_cropped/HII_rank_indices_files.rds", version = "2") 

# Generate masks
HII_mask_rank_5 <- continent_mask
HII_mask_rank_5[HII_indices_5] <- 1

plot(HII_mask_rank_5)

HII_mask_rank_25 <- continent_mask
HII_mask_rank_25[HII_indices_25] <- 1

plot(HII_mask_rank_25)

HII_mask_rank_75 <- continent_mask
HII_mask_rank_75[HII_indices_75] <- 1

plot(HII_mask_rank_75)

HII_mask_rank_95 <- continent_mask
HII_mask_rank_95[HII_indices_95] <- 1

plot(HII_mask_rank_95)

HII_mask_rank_25_75 <- continent_mask
HII_mask_rank_25_75[HII_indices_25_75] <- 1

plot(HII_mask_rank_25_75)

# Save masks files
save(HII_mask_rank_5, HII_mask_rank_25, HII_mask_rank_75, HII_mask_rank_95, HII_mask_rank_25_75, file = "./outputs/Threat_maps/HII_cropped/HII_masks_quantiles_files.rds", version = "2") 


# Generate shp files
HII_mask_rank_5_shp <- rasterToPolygons(x = HII_mask_rank_5, fun = function(x) {x == 1},
                                                dissolve = T)# To merge contiguous polygons of the same field value
HII_mask_rank_25_shp <- rasterToPolygons(x = HII_mask_rank_25, fun = function(x) {x == 1},
                                        dissolve = T)# To merge contiguous polygons of the same field value
HII_mask_rank_75_shp <- rasterToPolygons(x = HII_mask_rank_75, fun = function(x) {x == 1},
                                        dissolve = T)# To merge contiguous polygons of the same field value
HII_mask_rank_95_shp <- rasterToPolygons(x = HII_mask_rank_95, fun = function(x) {x == 1},
                                        dissolve = T)# To merge contiguous polygons of the same field value
HII_mask_rank_25_75_shp <- rasterToPolygons(x = HII_mask_rank_25_75, fun = function(x) {x == 1},
                                         dissolve = T)# To merge contiguous polygons of the same field value
# plot(HII_mask_rank_5_shp)


# sp.richness_mask_rank_5_shp_smooth <- gSimplify(sp.richness_mask_rank_5_shp, tol = 0.2, topologyPreserve = F)
# # plot(sp.richness_mask_rank_5_shp_smooth)
# sp.richness_mask_rank_25_shp_smooth <- gSimplify(sp.richness_mask_rank_25_shp, tol = 0.2, topologyPreserve = F)
# # plot(sp.richness_mask_rank_25_shp_smooth)

# Save shp files
save(HII_mask_rank_5_shp, HII_mask_rank_25_shp, HII_mask_rank_75_shp, HII_mask_rank_95_shp, HII_mask_rank_25_75_shp, file = "./outputs/Threat_maps/HII_cropped/HII_masked_rank_shp_files.rds", version = "2") 


### 2.1/ For absolute values ####

# Extract HII indices for communities with high/low HII

HII_indices_very_low <- which(HII_masked[] < 1)
HII_indices_low <- which(HII_masked[] < 5)
HII_indices_high <- which(HII_masked[] >= 11)
HII_indices_very_high <- which(HII_masked[] >= 23)
HII_indices_low_high <- c(HII_indices_low, HII_indices_high)
HII_indices_very_low_high <- c(HII_indices_very_low, HII_indices_high)

# Save indices
save(HII_indices_very_low, HII_indices_low, HII_indices_high, HII_indices_very_high, HII_indices_low_high, HII_indices_very_low_high, file = "./outputs/Threat_maps/HII_cropped/HII_indices_files.rds", version = "2") 


# Generate masks
HII_mask_very_low <- continent_mask
HII_mask_very_low[HII_indices_very_low] <- 1

plot(HII_mask_very_low)

HII_mask_low <- continent_mask
HII_mask_low[HII_indices_low] <- 1

plot(HII_mask_low)

HII_mask_high <- continent_mask
HII_mask_high[HII_indices_high] <- 1

plot(HII_mask_high)

HII_mask_very_high <- continent_mask
HII_mask_very_high[HII_indices_very_high] <- 1

plot(HII_mask_very_high)

HII_mask_low_high <- continent_mask
HII_mask_low_high[HII_indices_low_high] <- 1

plot(HII_mask_low_high)

HII_mask_very_low_high <- continent_mask
HII_mask_very_low_high[HII_indices_very_low_high] <- 1

plot(HII_mask_very_low_high)

# Save masks
save(HII_mask_very_low, HII_mask_low, HII_mask_high, HII_mask_very_high, HII_mask_low_high, HII_mask_very_low_high, file = "./outputs/Threat_maps/HII_cropped/HII_masks_bins_files.rds", version = "2") 


# Generate shp files
HII_mask_very_low_shp <- rasterToPolygons(x = HII_mask_very_low, fun = function(x) {x == 1},
                                        dissolve = T)# To merge contiguous polygons of the same field value
HII_mask_low_shp <- rasterToPolygons(x = HII_mask_low, fun = function(x) {x == 1},
                                         dissolve = T)# To merge contiguous polygons of the same field value
HII_mask_high_shp <- rasterToPolygons(x = HII_mask_high, fun = function(x) {x == 1},
                                         dissolve = T)# To merge contiguous polygons of the same field value
HII_mask_very_high_shp <- rasterToPolygons(x = HII_mask_very_high, fun = function(x) {x == 1},
                                         dissolve = T)# To merge contiguous polygons of the same field value
HII_mask_low_high_shp <- rasterToPolygons(x = HII_mask_low_high, fun = function(x) {x == 1},
                                           dissolve = T)# To merge contiguous polygons of the same field value
HII_mask_very_low_high_shp <- rasterToPolygons(x = HII_mask_very_low_high, fun = function(x) {x == 1},
                                          dissolve = T)# To merge contiguous polygons of the same field value

# plot(HII_mask_very_low_shp)


# sp.richness_mask_rank_5_shp_smooth <- gSimplify(sp.richness_mask_rank_5_shp, tol = 0.2, topologyPreserve = F)
# # plot(sp.richness_mask_rank_5_shp_smooth)
# sp.richness_mask_rank_25_shp_smooth <- gSimplify(sp.richness_mask_rank_25_shp, tol = 0.2, topologyPreserve = F)
# # plot(sp.richness_mask_rank_25_shp_smooth)

# Save shp files
save(HII_mask_very_low_shp, HII_mask_low_shp, HII_mask_high_shp, HII_mask_very_high_shp, HII_mask_low_high_shp, HII_mask_very_low_high_shp, file = "./outputs/Threat_maps/HII_cropped/HII_masked_values_shp_files.rds", version = "2") 



### Load mask rasters
load(file = "./outputs/Threat_maps/HII_cropped/HII_masks_quantiles_files.rds") 
load(file = "./outputs/Threat_maps/HII_cropped/HII_masks_bins_files.rds") 

### Load shp files
load(file = "./outputs/Threat_maps/HII_cropped/HII_masked_values_shp_files.rds")
load(file = "./outputs/Threat_maps/HII_cropped/HII_masked_rank_shp_files.rds") 


### 3/ Plot for species richness ####

# 3.1/ For quantiles ####

sp.richness_masked_HII_rank_25_75 <- sp.richness*HII_mask_rank_25_75

saveRDS(sp.richness_masked_HII_rank_25_75, file = "./outputs/Threat_maps/HII_cropped/sp.richness_masked_HII_rank_25_75")

plot(sp.richness_masked_HII_rank_25_75)

pdf(file = paste0("./maps/Threat_maps/HII_cropped/sp.richness_rank.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(sp.richness_masked_HII_rank_25_75, col = pal_bl_red_Mannion, main = paste0("Species richness: Threat and refuges\nQuantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(HII_mask_rank_25_shp, lwd = 1.0, border = "red", add = T)
plot(HII_mask_rank_5_shp, lwd = 0.8, border = "red4", col = NA, add = T)

plot(HII_mask_rank_75_shp, lwd = 1.0, border = "dodgerblue", add = T)
plot(HII_mask_rank_95_shp, lwd = 0.8, border = "navyblue", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(sp.richness, locs = seq(0, 120, 30), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -110, y = 1.0, font = 2, cex = 1.2, label = "Species")

legend(legend = "Top 5%", x = "bottomleft", bty = "n",
       pch = 22, col = "red4", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.20, 0.37))
legend(legend = "Top 25%", x = "bottomleft", bty = "n",
       pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.20, 0.30))
legend(legend = "Low 5%", x = "bottomleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.20, 0.23))
legend(legend = "Low 25%", x = "bottomleft", bty = "n",
       pch = 22, col = "dodgerblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.20, 0.16))

par(mar = internal_margins)
dev.off()


pdf(file = paste0("./maps/Threat_maps/HII_cropped/sp.richness_rank_light.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(sp.richness_masked_HII_rank_25_75, col = pal_bl_red_Mannion, main = paste0("Species richness: Threat and refuges\nQuantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(HII_mask_rank_25_shp, lwd = 1.0, border = "red", add = T)
plot(HII_mask_rank_75_shp, lwd = 1.0, border = "navyblue", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(sp.richness, locs = seq(0, 120, 30), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -110, y = 1.0, font = 2, cex = 1.2, label = "Species")

legend(legend = "Risk areas", x = "topleft", bty = "n",
       pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.19))
legend(legend = "Refuge areas", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.26))

par(mar = internal_margins)
dev.off()


# 3.2/ For absolute values ####

sp.richness_masked_HII_very_low_high <- sp.richness*HII_mask_very_low_high

saveRDS(sp.richness_masked_HII_very_low_high, file = "./outputs/Threat_maps/HII_cropped/sp.richness_masked_HII_very_low_high")


plot(sp.richness_masked_HII_very_low_high)

pdf(file = paste0("./maps/Threat_maps/HII_cropped/sp.richness_very_low_high.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(sp.richness_masked_HII_very_low_high, col = pal_bl_red_Mannion, main = paste0("Species richness: Threat and refuges\nAbsolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

# plot(HII_mask_low_shp, lwd = 1.0, border = "dodgerblue", add = T)
plot(HII_mask_very_low_shp, lwd = 1.0, border = "navyblue", col = NA, add = T)

plot(HII_mask_high_shp, lwd = 1.0, border = "red", add = T)
plot(HII_mask_very_high_shp, lwd = 0.8, border = "red4", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(sp.richness, locs = seq(0, 120, 30), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -110, y = 1.0, font = 2, cex = 1.2, label = "Species")

legend(legend = "Very high", x = "bottomleft", bty = "n",
       pch = 22, col = "red4", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.20, 0.30))
legend(legend = "High", x = "bottomleft", bty = "n",
       pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.20, 0.23))
legend(legend = "No threat", x = "bottomleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.20, 0.16))

par(mar = internal_margins)
dev.off()


pdf(file = paste0("./maps/Threat_maps/HII_cropped/sp.richness_very_low_high_light.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(sp.richness_masked_HII_very_low_high, col = pal_bl_red_Mannion, main = paste0("Species richness: Threat and refuges\nAbsolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

# plot(HII_mask_low_shp, lwd = 1.0, border = "dodgerblue", add = T)
plot(HII_mask_very_low_shp, lwd = 1.0, border = "navyblue", col = NA, add = T)

plot(HII_mask_high_shp, lwd = 1.0, border = "red", add = T)
# plot(HII_mask_very_high_shp, lwd = 0.8, border = "red4", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(sp.richness, locs = seq(0, 120, 30), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -110, y = 1.0, font = 2, cex = 1.2, label = "Species")

legend(legend = "Risk areas", x = "topleft", bty = "n",
       pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.19))
legend(legend = "Refuge areas", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.26))

par(mar = internal_margins)
dev.off()

### 4/ Plot for mean species rarity ####

# 4.1/ For quantiles ####

sp.rarity_masked_HII_rank_25_75 <- sp.rarity*HII_mask_rank_25_75

saveRDS(sp.rarity_masked_HII_rank_25_75, file = "./outputs/Threat_maps/HII_cropped/sp.rarity_masked_HII_rank_25_75")


plot(sp.rarity_masked_HII_rank_25_75)

pdf(file = paste0("./maps/Threat_maps/HII_cropped/sp.rarity_rank_light.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(sp.rarity_masked_HII_rank_25_75, col = pal_bl_red_Mannion, main = paste0("Mean species rarity: Threat and refuges\nQuantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(HII_mask_rank_25_shp, lwd = 1.0, border = "red", add = T)
plot(HII_mask_rank_75_shp, lwd = 1.0, border = "navyblue", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(sp.rarity, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -111, y = 2.5, font = 2, cex = 1.1, label = "Rarity\nindex")

legend(legend = "Risk areas", x = "topleft", bty = "n",
       pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.16))
legend(legend = "Refuge areas", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.23))

par(mar = internal_margins)
dev.off()


# 4.2/ For absolute values ####

sp.rarity_masked_HII_very_low_high <- sp.rarity*HII_mask_very_low_high

saveRDS(sp.rarity_masked_HII_very_low_high, file = "./outputs/Threat_maps/HII_cropped/sp.rarity_masked_HII_very_low_high")

plot(sp.rarity_masked_HII_very_low_high)

pdf(file = paste0("./maps/Threat_maps/HII_cropped/sp.rarity_very_low_high_light.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(sp.rarity_masked_HII_very_low_high, col = pal_bl_red_Mannion, main = paste0("Mean species rarity: Threat and refuges\nAbsolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

# plot(HII_mask_low_shp, lwd = 1.0, border = "dodgerblue", add = T)
plot(HII_mask_very_low_shp, lwd = 1.0, border = "navyblue", col = NA, add = T)

plot(HII_mask_high_shp, lwd = 1.0, border = "red", add = T)
# plot(HII_mask_very_high_shp, lwd = 0.8, border = "red4", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(sp.rarity, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -111, y = 2.5, font = 2, cex = 1.1, label = "Rarity\nindex")

legend(legend = "Risk areas", x = "topleft", bty = "n",
       pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.16))
legend(legend = "Refuge areas", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.23))

par(mar = internal_margins)
dev.off()

# ### 5/ Plot for MPD ####
# 
# # 5.1/ For quantiles ####
# 
# MPD_masked_HII_rank_25_75_contrasted <- MPD_contrasted*HII_mask_rank_25_75
# 
# saveRDS(MPD_masked_HII_rank_25_75_contrasted, file = "./outputs/Threat_maps/HII_cropped/MPD_masked_HII_rank_25_75_contrasted")
# 
# 
# plot(MPD_masked_HII_rank_25_75_contrasted)
# 
# pal_bl_red_Mannion_MDP <- c(rep("#EDEDED", 850), pal_bl_red_Mannion) # Special palette for contrasted MPD
# 
# pdf(file = paste0("./maps/Threat_maps/HII_cropped/MPD_rank_light.pdf"), height = 5.3, width = 6.5)
# 
# internal_margins <- par()$mar
# par(mar = c(3.1,3.5,3.5,2.1))
# 
# image(MPD_masked_HII_rank_25_75_contrasted, col = pal_bl_red_Mannion_MDP, main = paste0("MPD: Threat and refuges\nQuantiles"),
#       cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
#       legend  = T)
# plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)
# 
# plot(HII_mask_rank_25_shp, lwd = 1.0, border = "red", add = T)
# plot(HII_mask_rank_75_shp, lwd = 1.0, border = "navyblue", add = T)
# 
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
# addRasterLegend(MPD_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
# graphics::text(x = -110, y = 2.5, font = 2, cex = 1.1, label = "Evolutionary\nTime (My)")
# 
# legend(legend = "Risk areas", x = "topleft", bty = "n",
#        pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
#        text.font = 2, cex = 1.2, inset=c(0.02, 0.16))
# legend(legend = "Refuge areas", x = "topleft", bty = "n",
#        pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
#        text.font = 2, cex = 1.2, inset=c(0.02, 0.23))
# 
# par(mar = internal_margins)
# dev.off()
# 
# 
# # 5.2/ For absolute values ####
# 
# MPD_masked_HII_very_low_high_contrasted <- MPD_contrasted*HII_mask_very_low_high
# 
# saveRDS(MPD_masked_HII_very_low_high_contrasted, file = "./outputs/Threat_maps/HII_cropped/MPD_masked_HII_very_low_high_contrasted")
# 
# 
# plot(MPD_masked_HII_very_low_high_contrasted)
# 
# pdf(file = paste0("./maps/Threat_maps/HII_cropped/MPD_very_low_high_light.pdf"), height = 5.3, width = 6.5)
# 
# internal_margins <- par()$mar
# par(mar = c(3.1,3.5,3.5,2.1))
# 
# pal_bl_red_Mannion_MDP <- c(rep("#EDEDED", 850), pal_bl_red_Mannion) # Special palette for contrasted MPD
# 
# image(MPD_masked_HII_very_low_high_contrasted, col = pal_bl_red_Mannion_MDP, main = paste0("MPD: Threat and refuges\nAbsolute values"),
#       cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
#       legend  = T)
# plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)
# 
# # plot(HII_mask_low_shp, lwd = 1.0, border = "dodgerblue", add = T)
# plot(HII_mask_very_low_shp, lwd = 1.0, border = "navyblue", col = NA, add = T)
# 
# plot(HII_mask_high_shp, lwd = 1.0, border = "red", add = T)
# # plot(HII_mask_very_high_shp, lwd = 0.8, border = "red4", col = NA, add = T)
# 
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
# addRasterLegend(MPD_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
# graphics::text(x = -110, y = 2.5, font = 2, cex = 1.1, label = "Evolutionary\nTime (My)")
# 
# legend(legend = "Risk areas", x = "topleft", bty = "n",
#        pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
#        text.font = 2, cex = 1.2, inset=c(0.02, 0.16))
# legend(legend = "Refuge areas", x = "topleft", bty = "n",
#        pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
#        text.font = 2, cex = 1.2, inset=c(0.02, 0.23))
# 
# par(mar = internal_margins)
# dev.off()

### 6/ Plot for mimicry richness ####

# 6.1/ For quantiles ####

ring.richness_masked_HII_rank_25_75 <- ring.richness*HII_mask_rank_25_75

saveRDS(ring.richness_masked_HII_rank_25_75, file = "./outputs/Threat_maps/HII_cropped/ring.richness_masked_HII_rank_25_75")


plot(ring.richness_masked_HII_rank_25_75)

pdf(file = paste0("./maps/Threat_maps/HII_cropped/ring.richness_rank_light.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(ring.richness_masked_HII_rank_25_75, col = pal_bl_red_Mannion, main = paste0("Mimicry richness: Threat and refuges\nQuantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(HII_mask_rank_25_shp, lwd = 1.0, border = "red", add = T)
plot(HII_mask_rank_75_shp, lwd = 1.0, border = "navyblue", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(ring.richness, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -111, y = 2.5, font = 2, cex = 1.1, label = "Mimicry\nrings")

legend(legend = "Risk areas", x = "topleft", bty = "n",
       pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.16))
legend(legend = "Refuge areas", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.23))

par(mar = internal_margins)
dev.off()


# 6.2/ For absolute values ####

ring.richness_masked_HII_very_low_high <- ring.richness*HII_mask_very_low_high

saveRDS(ring.richness_masked_HII_very_low_high, file = "./outputs/Threat_maps/HII_cropped/ring.richness_masked_HII_very_low_high")


plot(ring.richness_masked_HII_very_low_high)

pdf(file = paste0("./maps/Threat_maps/HII_cropped/ring.richness_very_low_high_light.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(ring.richness_masked_HII_very_low_high, col = pal_bl_red_Mannion, main = paste0("Mimicry richness: Threat and refuges\nAbsolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

# plot(HII_mask_low_shp, lwd = 1.0, border = "dodgerblue", add = T)
plot(HII_mask_very_low_shp, lwd = 1.0, border = "navyblue", col = NA, add = T)

plot(HII_mask_high_shp, lwd = 1.0, border = "red", add = T)
# plot(HII_mask_very_high_shp, lwd = 0.8, border = "red4", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(ring.richness, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -111, y = 2.5, font = 2, cex = 1.1, label = "Mimicry\nrings")

legend(legend = "Risk areas", x = "topleft", bty = "n",
       pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.16))
legend(legend = "Refuge areas", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.23))

par(mar = internal_margins)
dev.off()

### 7/ Plot for mean mimicry rarity ####

# 7.1/ For quantiles ####

ring.rarity_masked_HII_rank_25_75 <- ring.rarity*HII_mask_rank_25_75

saveRDS(ring.rarity_masked_HII_rank_25_75, file = "./outputs/Threat_maps/HII_cropped/ring.rarity_masked_HII_rank_25_75")


plot(ring.rarity_masked_HII_rank_25_75)

pdf(file = paste0("./maps/Threat_maps/HII_cropped/ring.rarity_rank_light.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(ring.rarity_masked_HII_rank_25_75, col = pal_bl_red_Mannion, main = paste0("Mean mimicry rarity: Threat and refuges\nQuantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Mimicry\nrings", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(HII_mask_rank_25_shp, lwd = 1.0, border = "red", add = T)
plot(HII_mask_rank_75_shp, lwd = 1.0, border = "navyblue", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(ring.rarity, locs = seq(0, 0.6, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -111, y = 2.5, font = 2, cex = 1.1, label = "Rarity\nindex")

legend(legend = "Risk areas", x = "topleft", bty = "n",
       pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.16))
legend(legend = "Refuge areas", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.23))

par(mar = internal_margins)
dev.off()


# 7.2/ For absolute values ####

ring.rarity_masked_HII_very_low_high <- ring.rarity*HII_mask_very_low_high

saveRDS(ring.rarity_masked_HII_very_low_high, file = "./outputs/Threat_maps/HII_cropped/ring.rarity_masked_HII_very_low_high")

plot(ring.rarity_masked_HII_very_low_high)

pdf(file = paste0("./maps/Threat_maps/HII_cropped/ring.rarity_very_low_high_light.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(ring.rarity_masked_HII_very_low_high, col = pal_bl_red_Mannion, main = paste0("Mean mimicry rarity: Threat and refuges\nAbsolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Mimicry\nrings", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

# plot(HII_mask_low_shp, lwd = 1.0, border = "dodgerblue", add = T)
plot(HII_mask_very_low_shp, lwd = 1.0, border = "navyblue", col = NA, add = T)

plot(HII_mask_high_shp, lwd = 1.0, border = "red", add = T)
# plot(HII_mask_very_high_shp, lwd = 0.8, border = "red4", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(ring.rarity, locs = seq(0, 0.6, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -111, y = 2.5, font = 2, cex = 1.1, label = "Rarity\nindex")

legend(legend = "Risk areas", x = "topleft", bty = "n",
       pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.16))
legend(legend = "Refuge areas", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.23))

par(mar = internal_margins)
dev.off()



### 8/ All four indices ####

# 8.1/ For quantiles ####

sp.richness_masked_HII_rank_25_75 <- readRDS(file = "./outputs/Threat_maps/HII_cropped/sp.richness_masked_HII_rank_25_75")
sp.rarity_masked_HII_rank_25_75 <- readRDS(file = "./outputs/Threat_maps/HII_cropped/sp.rarity_masked_HII_rank_25_75")
# MPD_masked_HII_rank_25_75_contrasted <- readRDS(file = "./outputs/Threat_maps/HII_cropped/MPD_masked_HII_rank_25_75_contrasted")
ring.richness_masked_HII_rank_25_75 <- readRDS(file = "./outputs/Threat_maps/HII_cropped/ring.richness_masked_HII_rank_25_75")
ring.rarity_masked_HII_rank_25_75 <- readRDS(file = "./outputs/Threat_maps/HII_cropped/ring.rarity_masked_HII_rank_25_75")


pdf(file = paste0("./maps/Threat_maps/HII_cropped/All_indices_rank_light.pdf"), height = 10, width = 12)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
par(mfrow = c(2,2))


# A/ Species richness

image(sp.richness_masked_HII_rank_25_75, col = pal_bl_red_Mannion, main = paste0("Species richness: Threat and refuges\nQuantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(HII_mask_rank_25_shp, lwd = 1.0, border = "red", add = T)
plot(HII_mask_rank_75_shp, lwd = 1.0, border = "navyblue", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(sp.richness, locs = seq(0, 120, 30), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -110, y = 1.0, font = 2, cex = 1.2, label = "Species")

legend(legend = "Risk areas", x = "topleft", bty = "n",
       pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.19))
legend(legend = "Refuge areas", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.26))

legend(legend = "a", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0.0, 0.0))


# B/ Mean species rarity

image(sp.rarity_masked_HII_rank_25_75, col = pal_bl_red_Mannion, main = paste0("Mean species rarity: Threat and refuges\nQuantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(HII_mask_rank_25_shp, lwd = 1.0, border = "red", add = T)
plot(HII_mask_rank_75_shp, lwd = 1.0, border = "navyblue", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(sp.rarity, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -111, y = 2.5, font = 2, cex = 1.1, label = "Rarity\nindex")

legend(legend = "Risk areas", x = "topleft", bty = "n",
       pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.16))
legend(legend = "Refuge areas", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.23))

legend(legend = "b", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0.0, 0.0))


# # C/ Mean pairwise Phylogenetic Distance
# 
# pal_bl_red_Mannion_MDP <- c(rep("#EDEDED", 850), pal_bl_red_Mannion) # Special palette for contrasted MPD
# 
# image(MPD_masked_HII_rank_25_75_contrasted, col = pal_bl_red_Mannion_MDP, main = paste0("MPD: Threat and refuges\nQuantiles"),
#       cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
#       legend  = T)
# plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)
# 
# plot(HII_mask_rank_25_shp, lwd = 1.0, border = "red", add = T)
# plot(HII_mask_rank_75_shp, lwd = 1.0, border = "navyblue", add = T)
# 
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
# addRasterLegend(MPD_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
# graphics::text(x = -110, y = 2.5, font = 2, cex = 1.1, label = "Evolutionary\nTime (My)")
# 
# legend(legend = "Risk areas", x = "topleft", bty = "n",
#        pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
#        text.font = 2, cex = 1.2, inset=c(0.02, 0.16))
# legend(legend = "Refuge areas", x = "topleft", bty = "n",
#        pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
#        text.font = 2, cex = 1.2, inset=c(0.02, 0.23))
# 
# legend(legend = "c", x = "bottomright", bty = "n",
#        text.font = 2, cex = 1.8, inset=c(0.0, 0.0))


# C/ Mimicry richness

image(ring.richness_masked_HII_rank_25_75, col = pal_bl_red_Mannion, main = paste0("Mimicry richness: Threat and refuges\nQuantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(HII_mask_rank_25_shp, lwd = 1.0, border = "red", add = T)
plot(HII_mask_rank_75_shp, lwd = 1.0, border = "navyblue", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(ring.richness, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -111, y = 2.5, font = 2, cex = 1.1, label = "Mimicry\nrings")

legend(legend = "Risk areas", x = "topleft", bty = "n",
       pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.16))
legend(legend = "Refuge areas", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.23))

legend(legend = "c", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0.0, 0.0))

# D/ Mean mimicry rarity

image(ring.rarity_masked_HII_rank_25_75, col = pal_bl_red_Mannion, main = paste0("Mean mimicry rarity: Threat and refuges\nQuantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Mimicry\nrings", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(HII_mask_rank_25_shp, lwd = 1.0, border = "red", add = T)
plot(HII_mask_rank_75_shp, lwd = 1.0, border = "navyblue", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(ring.rarity, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -111, y = 2.5, font = 2, cex = 1.1, label = "Rarity\nindex")

legend(legend = "Risk areas", x = "topleft", bty = "n",
       pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.16))
legend(legend = "Refuge areas", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.23))

legend(legend = "d", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0.0, 0.0))

par(mar = internal_margins)
par(mfrow = c(1, 1))

dev.off()

# 8.2/ For absolute values ####

sp.richness_masked_HII_very_low_high <- readRDS(file = "./outputs/Threat_maps/HII_cropped/sp.richness_masked_HII_very_low_high")
sp.rarity_masked_HII_very_low_high <- readRDS(file = "./outputs/Threat_maps/HII_cropped/sp.rarity_masked_HII_very_low_high")
MPD_masked_HII_very_low_high_contrasted <- readRDS(file = "./outputs/Threat_maps/HII_cropped/MPD_masked_HII_very_low_high_contrasted")
ring.richness_masked_HII_very_low_high <- readRDS(file = "./outputs/Threat_maps/HII_cropped/ring.richness_masked_HII_very_low_high")
ring.rarity_masked_HII_very_low_high <- readRDS(file = "./outputs/Threat_maps/HII_cropped/ring.rarity_masked_HII_very_low_high")

pdf(file = paste0("./maps/Threat_maps/HII_cropped/All_indices_very_low_high_light.pdf"), height = 10, width = 12)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
par(mfrow = c(2,2))


# A/ Species richness

image(sp.richness_masked_HII_very_low_high, col = pal_bl_red_Mannion, main = paste0("Species richness: Threat and refuges\nAbsolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

# plot(HII_mask_low_shp, lwd = 1.0, border = "dodgerblue", add = T)
plot(HII_mask_very_low_shp, lwd = 1.0, border = "navyblue", col = NA, add = T)

plot(HII_mask_high_shp, lwd = 1.0, border = "red", add = T)
# plot(HII_mask_very_high_shp, lwd = 0.8, border = "red4", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(sp.richness, locs = seq(0, 120, 30), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -110, y = 1.0, font = 2, cex = 1.2, label = "Species")

legend(legend = "Risk areas", x = "topleft", bty = "n",
       pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.19))
legend(legend = "Refuge areas", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.26))

legend(legend = "a", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0.0, 0.0))


# B/ Mean species rarity

image(sp.rarity_masked_HII_very_low_high, col = pal_bl_red_Mannion, main = paste0("Mean species rarity: Threat and refuges\nAbsolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

# plot(HII_mask_low_shp, lwd = 1.0, border = "dodgerblue", add = T)
plot(HII_mask_very_low_shp, lwd = 1.0, border = "navyblue", col = NA, add = T)

plot(HII_mask_high_shp, lwd = 1.0, border = "red", add = T)
# plot(HII_mask_very_high_shp, lwd = 0.8, border = "red4", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(sp.rarity, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -111, y = 2.5, font = 2, cex = 1.1, label = "Rarity\nindex")

legend(legend = "Risk areas", x = "topleft", bty = "n",
       pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.16))
legend(legend = "Refuge areas", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.23))

legend(legend = "b", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0.0, 0.0))


# # C/ Mean pairwise Phylogenetic Distance
# 
# pal_bl_red_Mannion_MDP <- c(rep("#EDEDED", 850), pal_bl_red_Mannion) # Special palette for contrasted MPD
# 
# image(MPD_masked_HII_very_low_high_contrasted, col = pal_bl_red_Mannion_MDP, main = paste0("MPD: Threat and refuges\nAbsolute values"),
#       cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
#       legend  = T)
# plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)
# 
# # plot(HII_mask_low_shp, lwd = 1.0, border = "dodgerblue", add = T)
# plot(HII_mask_very_low_shp, lwd = 1.0, border = "navyblue", col = NA, add = T)
# 
# plot(HII_mask_high_shp, lwd = 1.0, border = "red", add = T)
# # plot(HII_mask_very_high_shp, lwd = 0.8, border = "red4", col = NA, add = T)
# 
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
# addRasterLegend(MPD_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
# graphics::text(x = -110, y = 2.5, font = 2, cex = 1.1, label = "Evolutionary\nTime (My)")
# 
# legend(legend = "Risk areas", x = "topleft", bty = "n",
#        pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
#        text.font = 2, cex = 1.2, inset=c(0.02, 0.16))
# legend(legend = "Refuge areas", x = "topleft", bty = "n",
#        pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
#        text.font = 2, cex = 1.2, inset=c(0.02, 0.23))
# 
# legend(legend = "c", x = "bottomright", bty = "n",
#        text.font = 2, cex = 1.8, inset=c(0.0, 0.0))
# 

# C/ Mimicry richness

image(ring.richness_masked_HII_very_low_high, col = pal_bl_red_Mannion, main = paste0("Mimicry richness: Threat and refuges\nAbsolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

# plot(HII_mask_low_shp, lwd = 1.0, border = "dodgerblue", add = T)
plot(HII_mask_very_low_shp, lwd = 1.0, border = "navyblue", col = NA, add = T)

plot(HII_mask_high_shp, lwd = 1.0, border = "red", add = T)
# plot(HII_mask_very_high_shp, lwd = 0.8, border = "red4", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(ring.richness, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -111, y = 2.5, font = 2, cex = 1.1, label = "Mimicry\nrings")

legend(legend = "Risk areas", x = "topleft", bty = "n",
       pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.16))
legend(legend = "Refuge areas", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.23))

legend(legend = "c", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0.0, 0.0))

# D/ Mean mimicry rarity

image(ring.rarity_masked_HII_very_low_high, col = pal_bl_red_Mannion, main = paste0("Mean mimicry rarity: Threat and refuges\nAbsolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Mimicry\nrings", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

# plot(HII_mask_low_shp, lwd = 1.0, border = "dodgerblue", add = T)
plot(HII_mask_very_low_shp, lwd = 1.0, border = "navyblue", col = NA, add = T)

plot(HII_mask_high_shp, lwd = 1.0, border = "red", add = T)
# plot(HII_mask_very_high_shp, lwd = 0.8, border = "red4", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(ring.rarity, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -111, y = 2.5, font = 2, cex = 1.1, label = "Rarity\nindex")

legend(legend = "Risk areas", x = "topleft", bty = "n",
       pch = 22, col = "red", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.16))
legend(legend = "Refuge areas", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.02, 0.23))

legend(legend = "d", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0.0, 0.0))


par(mar = internal_margins)
par(mfrow = c(1, 1))

dev.off()


### 9/ Boxplots for quantiles ####

load(file = "./outputs/Threat_maps/HII_cropped/HII_rank_indices_files.rds") 

# Load indices map
sp.richness <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
# sp.rarity <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_Jaccard.80.rds")        # Linear weighting
sp.rarity <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_Leroy_Jaccard.80.rds")  # Leroy's weighting
# MPD <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80.rds")
ring.richness <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))
# ring.rarity <- readRDS(file = paste0("./outputs/Indices_maps/mimicry.mean.rarity_Jaccard.80.rds"))   # Linear weighting
ring.rarity <- readRDS(file = "./outputs/Indices_maps/ring.mean.rarity_Leroy_Jaccard.80.rds")         # Leroy's weighting


# 9.1/ Extract values ####

library(raster)
library(tidyverse)

# all_indices_stack <- stack(sp.richness, sp.rarity, MPD, ring.richness)
all_indices_stack <- stack(sp.richness, sp.rarity, ring.richness, ring.rarity)
# names(all_indices_stack) <- c("sp.richness", "sp.rarity", "MPD", "ring.richness")
names(all_indices_stack) <- c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity")

values_5 <- data.frame(all_indices_stack[HII_indices_5])
values_25 <- data.frame(all_indices_stack[HII_indices_25])
values_75 <- data.frame(all_indices_stack[HII_indices_75])
values_95 <- data.frame(all_indices_stack[HII_indices_95])

values_5 <- values_5 %>% 
  pivot_longer(# cols = c("sp.richness", "sp.rarity", "MPD", "ring.richness"), # List of the name of the "columns" that are actually factor levels
               cols = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"), # List of the name of the "columns" that are actually factor levels
               names_to = "indices",   # name of the new column that is the real variable for these factor level
               values_to = "values", # name of the new column(s) that will get the values in the cells of the previous dataset
               values_drop_na = T) %>% # Drop rows for missing values ?  
  mutate(threshold = 5)

values_25 <- values_25 %>% 
  pivot_longer(# cols = c("sp.richness", "sp.rarity", "MPD", "ring.richness"), # List of the name of the "columns" that are actually factor levels
               cols = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"), # List of the name of the "columns" that are actually factor levels
               names_to = "indices",   # name of the new column that is the real variable for these factor level
               values_to = "values", # name of the new column(s) that will get the values in the cells of the previous dataset
               values_drop_na = T) %>% # Drop rows for missing values ?  
  mutate(threshold = 25)

values_75 <- values_75 %>% 
  pivot_longer(# cols = c("sp.richness", "sp.rarity", "MPD", "ring.richness"), # List of the name of the "columns" that are actually factor levels
               cols = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"), # List of the name of the "columns" that are actually factor levels
               names_to = "indices",   # name of the new column that is the real variable for these factor level
               values_to = "values", # name of the new column(s) that will get the values in the cells of the previous dataset
               values_drop_na = T) %>% # Drop rows for missing values ?  
  mutate(threshold = 75)

values_95 <- values_95 %>% 
  pivot_longer(# cols = c("sp.richness", "sp.rarity", "MPD", "ring.richness"), # List of the name of the "columns" that are actually factor levels
               cols = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"), # List of the name of the "columns" that are actually factor levels
               names_to = "indices",   # name of the new column that is the real variable for these factor level
               values_to = "values", # name of the new column(s) that will get the values in the cells of the previous dataset
               values_drop_na = T) %>% # Drop rows for missing values ?  
  mutate(threshold = 95)

HII_extract_for_quantiles_df <- rbind(values_5, values_25, values_75, values_95)
HII_extract_for_quantiles_df$threshold <- as.factor(HII_extract_for_quantiles_df$threshold)

saveRDS(HII_extract_for_quantiles_df, file = paste0("./outputs/Threat_maps/HII_cropped/HII_extract_for_quantiles_df.rds"))

# 9.2/ Plot boxplots ####

# 9.2.1/ Species richness

pdf(file = paste0("./graphs/Risks_refuges/Boxplot_HII_quantiles_sp.richness.pdf"), height = 5, width = 6)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

Boxplot_HII_quantiles_sp.richness <- HII_extract_for_quantiles_df %>% 
  filter(indices == "sp.richness") %>% 
  ggplot(aes(x = threshold, y = values, fill = threshold)) +
  # geom_violin(draw_quantiles = 0.5, lwd = 1.2) +
  geom_boxplot(lwd = 0.8, show.legend = F) +
  scale_fill_manual(values = c("5" = "red4", "25" = "red", "75" = "dodgerblue", "95" = "navyblue"),
                    labels=c("5" = "Top 5%", "25" = "Top 25%", "75" = "Low 25%", "95" = "Low 5%")) +
  scale_x_discrete(labels=c("5" = "Top 5%", "25" = "Top 25%", "75" = "Low 25%", "95" = "Low 5%")) +
  labs(title = "Species richness in risk areas and refuges", x = NULL, y = "Species richness", fill = "Regions") +
  theme_classic() +
  theme(title = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 5, 0)),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 15),
        axis.line = element_line(color = "black"))

print(Boxplot_HII_quantiles_sp.richness)

par(mar = internal_margins)
dev.off()

# 9.2.2/ Species rarity

pdf(file = paste0("./graphs/Risks_refuges/Boxplot_HII_quantiles_sp.rarity.pdf"), height = 5, width = 6)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

Boxplot_HII_quantiles_sp.rarity <- HII_extract_for_quantiles_df %>% 
  filter(indices == "sp.rarity") %>% 
  ggplot(aes(x = threshold, y = values, fill = threshold)) +
  # geom_violin(draw_quantiles = 0.5, lwd = 1.2) +
  geom_boxplot(lwd = 0.8, show.legend = F) +
  scale_fill_manual(values = c("5" = "red4", "25" = "red", "75" = "dodgerblue", "95" = "navyblue"),
                    labels=c("5" = "Top 5%", "25" = "Top 25%", "75" = "Low 25%", "95" = "Low 5%")) +
  scale_x_discrete(labels=c("5" = "Top 5%", "25" = "Top 25%", "75" = "Low 25%", "95" = "Low 5%")) +
  labs(title = "Mean species rarity in risk areas and refuges", x = NULL, y = "Mean species rarity", fill = "Regions") +
  theme_classic() +
  theme(title = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 5, 0)),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 15),
        axis.line = element_line(color = "black"))

print(Boxplot_HII_quantiles_sp.rarity)

par(mar = internal_margins)
dev.off()

# # 9.2.3/ MPD
# 
# pdf(file = paste0("./graphs/Risks_refuges/Boxplot_HII_quantiles_MPD.pdf"), height = 5, width = 6)
# 
# internal_margins <- par()$mar
# par(mar = c(3.1,3.5,3.5,2.1))
# 
# Boxplot_HII_quantiles_MPD <- HII_extract_for_quantiles_df %>% 
#   filter(indices == "MPD") %>% 
#   ggplot(aes(x = threshold, y = values, fill = threshold)) +
#   # geom_violin(draw_quantiles = 0.5, lwd = 1.2) +
#   geom_boxplot(lwd = 0.8, show.legend = F) +
#   scale_fill_manual(values = c("5" = "red4", "25" = "red", "75" = "dodgerblue", "95" = "navyblue"),
#                     labels=c("5" = "Top 5%", "25" = "Top 25%", "75" = "Low 25%", "95" = "Low 5%")) +
#   scale_x_discrete(labels=c("5" = "Top 5%", "25" = "Top 25%", "75" = "Low 25%", "95" = "Low 5%")) +
#   labs(title = "MPD in risk areas and refuges", x = NULL, y = "Mean pairwise Phylogenetic Distance", fill = "Regions") +
#   theme_classic() +
#   theme(title = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 5, 0)),
#         axis.text = element_text(size = 12, face = "bold", color = "black"),
#         axis.title = element_text(size = 15),
#         axis.line = element_line(color = "black"))
# 
# print(Boxplot_HII_quantiles_MPD)
# 
# par(mar = internal_margins)
# dev.off()

# 9.2.3/ Mimicry richness

pdf(file = paste0("./graphs/Risks_refuges/Boxplot_HII_quantiles_ring.richness.pdf"), height = 5, width = 6)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

Boxplot_HII_quantiles_ring.richness <- HII_extract_for_quantiles_df %>% 
  filter(indices == "ring.richness") %>% 
  ggplot(aes(x = threshold, y = values, fill = threshold)) +
  # geom_violin(draw_quantiles = 0.5, lwd = 1.2) +
  geom_boxplot(lwd = 0.8, show.legend = F) +
  scale_fill_manual(values = c("5" = "red4", "25" = "red", "75" = "dodgerblue", "95" = "navyblue"),
                    labels=c("5" = "Top 5%", "25" = "Top 25%", "75" = "Low 25%", "95" = "Low 5%")) +
  scale_x_discrete(labels=c("5" = "Top 5%", "25" = "Top 25%", "75" = "Low 25%", "95" = "Low 5%")) +
  labs(title = "Mimicry richness in risk areas and refuges", x = NULL, y = "Mimicry richness", fill = "Regions") +
  theme_classic() +
  theme(title = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 5, 0)),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 15),
        axis.line = element_line(color = "black"))

print(Boxplot_HII_quantiles_ring.richness)

par(mar = internal_margins)
dev.off()


# 9.2.4/ Mimicry rarity

pdf(file = paste0("./graphs/Risks_refuges/Boxplot_HII_quantiles_ring.rarity.pdf"), height = 5, width = 6)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

Boxplot_HII_quantiles_ring.rarity <- HII_extract_for_quantiles_df %>% 
  filter(indices == "ring.rarity") %>% 
  ggplot(aes(x = threshold, y = values, fill = threshold)) +
  # geom_violin(draw_quantiles = 0.5, lwd = 1.2) +
  geom_boxplot(lwd = 0.8, show.legend = F) +
  scale_fill_manual(values = c("5" = "red4", "25" = "red", "75" = "dodgerblue", "95" = "navyblue"),
                    labels=c("5" = "Top 5%", "25" = "Top 25%", "75" = "Low 25%", "95" = "Low 5%")) +
  scale_x_discrete(labels=c("5" = "Top 5%", "25" = "Top 25%", "75" = "Low 25%", "95" = "Low 5%")) +
  labs(title = "Mean mimicry rarity in risk areas and refuges", x = NULL, y = "Mean mimicry rarity", fill = "Regions") +
  theme_classic() +
  theme(title = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 5, 0)),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 15),
        axis.line = element_line(color = "black"))

print(Boxplot_HII_quantiles_ring.rarity)

par(mar = internal_margins)
dev.off()



# 9.2.5/ All four indices

library(gridExtra) # To plot gglots side-by-side

pdf(file = paste0("./graphs/Risks_refuges/Boxplot_HII_quantiles_all_indices.pdf"), height = 10, width = 12)

# gridExtra::grid.arrange(Boxplot_HII_quantiles_sp.richness, Boxplot_HII_quantiles_sp.rarity, Boxplot_HII_quantiles_MPD, Boxplot_HII_quantiles_ring.richness, ncol = 2)
gridExtra::grid.arrange(Boxplot_HII_quantiles_sp.richness, Boxplot_HII_quantiles_sp.rarity, Boxplot_HII_quantiles_ring.richness, Boxplot_HII_quantiles_ring.rarity, ncol = 2)

dev.off()


# 10/ Boxplots for absolute values ####

load(file = "./outputs/Threat_maps/HII_cropped/HII_indices_files.rds") 

# Load indices map
sp.richness <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
# sp.rarity <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_Jaccard.80.rds")        # Linear weighting
sp.rarity <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_Leroy_Jaccard.80.rds")  # Leroy's weighting
# MPD <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80.rds")
ring.richness <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))
# ring.rarity <- readRDS(file = paste0("./outputs/Indices_maps/mimicry.mean.rarity_Jaccard.80.rds"))   # Linear weighting
ring.rarity <- readRDS(file = "./outputs/Indices_maps/ring.mean.rarity_Leroy_Jaccard.80.rds")         # Leroy's weighting



# 10.1/ Extract values ####

library(raster)
library(tidyverse)

# all_indices_stack <- stack(sp.richness, sp.rarity, MPD, ring.richness)
all_indices_stack <- stack(sp.richness, sp.rarity, ring.richness, ring.rarity)
# names(all_indices_stack) <- c("sp.richness", "sp.rarity", "MPD", "ring.richness")
names(all_indices_stack) <- c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity")

values_very_high <- data.frame(all_indices_stack[HII_indices_very_high])
values_high <- data.frame(all_indices_stack[HII_indices_high])
values_low <- data.frame(all_indices_stack[HII_indices_low])
values_very_low <- data.frame(all_indices_stack[HII_indices_very_low])

values_very_high <- values_very_high %>% 
  pivot_longer(# cols = c("sp.richness", "sp.rarity", "MPD", "ring.richness"), # List of the name of the "columns" that are actually factor levels
               cols = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"), # List of the name of the "columns" that are actually factor levels
               names_to = "indices",   # name of the new column that is the real variable for these factor level
               values_to = "values", # name of the new column(s) that will get the values in the cells of the previous dataset
               values_drop_na = T) %>% # Drop rows for missing values ?  
  mutate(threshold = "very_high")

values_high <- values_high %>% 
  pivot_longer(# cols = c("sp.richness", "sp.rarity", "MPD", "ring.richness"), # List of the name of the "columns" that are actually factor levels
               cols = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"), # List of the name of the "columns" that are actually factor levels
               names_to = "indices",   # name of the new column that is the real variable for these factor level
               values_to = "values", # name of the new column(s) that will get the values in the cells of the previous dataset
               values_drop_na = T) %>% # Drop rows for missing values ?  
  mutate(threshold = "high")

values_low <- values_low %>% 
  pivot_longer(# cols = c("sp.richness", "sp.rarity", "MPD", "ring.richness"), # List of the name of the "columns" that are actually factor levels
               cols = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"), # List of the name of the "columns" that are actually factor levels
               names_to = "indices",   # name of the new column that is the real variable for these factor level
               values_to = "values", # name of the new column(s) that will get the values in the cells of the previous dataset
               values_drop_na = T) %>% # Drop rows for missing values ?  
  mutate(threshold = "low")

values_very_low <- values_very_low %>% 
  pivot_longer(# cols = c("sp.richness", "sp.rarity", "MPD", "ring.richness"), # List of the name of the "columns" that are actually factor levels
               cols = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"), # List of the name of the "columns" that are actually factor levels
               names_to = "indices",   # name of the new column that is the real variable for these factor level
               values_to = "values", # name of the new column(s) that will get the values in the cells of the previous dataset
               values_drop_na = T) %>% # Drop rows for missing values ?  
  mutate(threshold = "very_low")

HII_extract_for_bins_df <- rbind(values_very_high, values_high, values_low, values_very_low)
HII_extract_for_bins_df$threshold <- as.factor(HII_extract_for_bins_df$threshold)

saveRDS(HII_extract_for_bins_df, file = paste0("./outputs/Threat_maps/HII_cropped/HII_extract_for_bins_df.rds"))

# 10.2/ Plot boxplots ####

# 10.2.1/ Species richness

pdf(file = paste0("./graphs/Risks_refuges/Boxplot_HII_bins_sp.richness.pdf"), height = 5, width = 6)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,5.5,2.1))

Boxplot_HII_bins_sp.richness <- HII_extract_for_bins_df %>% 
  filter(indices == "sp.richness") %>% 
  ggplot(aes(x = threshold, y = values, fill = threshold)) +
  # geom_violin(draw_quantiles = 0.5, lwd = 1.2) +
  geom_boxplot(lwd = 0.8, show.legend = F) +
  scale_fill_manual(values = c("very_high" = "red4", "high" = "red", "low" = "dodgerblue", "very_low" = "navyblue"),
                    labels = c("very_high" = "Very high\n[23-100]", "high" = "High\n[11-100]", "low" = "Low\n[0-5[", "very_low" = "Very low\n[0-1[")) +
  scale_x_discrete(labels = c("very_high" = "Very high", "high" = "High",  "low" = "Low", "very_low" = "Very low"),
                   limits = c("very_high", "high", "low", "very_low")) +
  labs(title = "Species richness in risk areas and refuges", x = NULL, y = "Species richness", fill = "Regions") +
  theme_classic() +
  theme(title = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 5, 0)),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 15),
        axis.line = element_line(color = "black"))

print(Boxplot_HII_bins_sp.richness)

par(mar = internal_margins)
dev.off()

# 10.2.2/ Species rarity

pdf(file = paste0("./graphs/Risks_refuges/Boxplot_HII_bins_sp.rarity.pdf"), height = 5, width = 6)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

Boxplot_HII_bins_sp.rarity <- HII_extract_for_bins_df %>% 
  filter(indices == "sp.rarity") %>% 
  ggplot(aes(x = threshold, y = values, fill = threshold)) +
  # geom_violin(draw_quantiles = 0.5, lwd = 1.2) +
  geom_boxplot(lwd = 0.8, show.legend = F) +
  scale_fill_manual(values = c("very_high" = "red4", "high" = "red", "low" = "dodgerblue", "very_low" = "navyblue"),
                    labels=c("very_high" = "Very high\n[23-100]", "high" = "High\n[11-100]", "low" = "Low\n[0-5[", "very_low" = "Very low\n[0-1[")) +
  scale_x_discrete(labels = c("very_high" = "Very high", "high" = "High",  "low" = "Low", "very_low" = "Very low"),
                   limits = c("very_high", "high", "low", "very_low")) +
  labs(title = "Mean species rarity in risk areas and refuges", x = NULL, y = "Mean species rarity", fill = "Regions") +
  theme_classic() +
  theme(title = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 5, 0)),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 15),
        axis.line = element_line(color = "black"))

print(Boxplot_HII_bins_sp.rarity)

par(mar = internal_margins)
dev.off()

# # 10.2.3/ MPD
# 
# pdf(file = paste0("./graphs/Risks_refuges/Boxplot_HII_bins_MPD.pdf"), height = 5, width = 6)
# 
# internal_margins <- par()$mar
# par(mar = c(3.1,3.5,3.5,2.1))
# 
# Boxplot_HII_bins_MPD <- HII_extract_for_bins_df %>% 
#   filter(indices == "MPD") %>% 
#   ggplot(aes(x = threshold, y = values, fill = threshold)) +
#   # geom_violin(draw_quantiles = 0.5, lwd = 1.2) +
#   geom_boxplot(lwd = 0.8, show.legend = F) +
#   scale_fill_manual(values = c("very_high" = "red4", "high" = "red", "low" = "dodgerblue", "very_low" = "navyblue"),
#                     labels=c("very_high" = "Very high\n[23-100]", "high" = "High\n[11-100]", "low" = "Low\n[0-5[", "very_low" = "Very low\n[0-1[")) +
#   scale_x_discrete(labels = c("very_high" = "Very high", "high" = "High",  "low" = "Low", "very_low" = "Very low"),
#                    limits = c("very_high", "high", "low", "very_low")) +
#   labs(title = "MPD in risk areas and refuges", x = NULL, y = "Mean pairwise Phylogenetic Distance", fill = "Regions") +
#   theme_classic() +
#   theme(title = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 5, 0)),
#         axis.text = element_text(size = 12, face = "bold", color = "black"),
#         axis.title = element_text(size = 15),
#         axis.line = element_line(color = "black"))
# 
# print(Boxplot_HII_bins_MPD)
# 
# par(mar = internal_margins)
# dev.off()

# 10.2.3/ Mimicry richness

pdf(file = paste0("./graphs/Risks_refuges/Boxplot_HII_bins_ring.richness.pdf"), height = 5, width = 6)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

Boxplot_HII_bins_ring.richness <- HII_extract_for_bins_df %>% 
  filter(indices == "ring.richness") %>% 
  ggplot(aes(x = threshold, y = values, fill = threshold)) +
  # geom_violin(draw_quantiles = 0.5, lwd = 1.2) +
  geom_boxplot(lwd = 0.8, show.legend = F) +
  scale_fill_manual(values = c("very_high" = "red4", "high" = "red", "low" = "dodgerblue", "very_low" = "navyblue"),
                    labels=c("very_high" = "Very high\n[23-100]", "high" = "High\n[11-100]", "low" = "Low\n[0-5[", "very_low" = "Very low\n[0-1[")) +
  scale_x_discrete(labels = c("very_high" = "Very high", "high" = "High",  "low" = "Low", "very_low" = "Very low"),
                   limits = c("very_high", "high", "low", "very_low")) +
  labs(title = "Mimicry richness in risk areas and refuges", x = NULL, y = "Mimicry richness", fill = "Regions") +
  theme_classic() +
  theme(title = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 5, 0)),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 15),
        axis.line = element_line(color = "black"))

print(Boxplot_HII_bins_ring.richness)

par(mar = internal_margins)
dev.off()

# 10.2.4/ Mimicry rarity

pdf(file = paste0("./graphs/Risks_refuges/Boxplot_HII_bins_ring.rarity.pdf"), height = 5, width = 6)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

Boxplot_HII_bins_ring.rarity <- HII_extract_for_bins_df %>% 
  filter(indices == "ring.rarity") %>% 
  ggplot(aes(x = threshold, y = values, fill = threshold)) +
  # geom_violin(draw_quantiles = 0.5, lwd = 1.2) +
  geom_boxplot(lwd = 0.8, show.legend = F) +
  scale_fill_manual(values = c("very_high" = "red4", "high" = "red", "low" = "dodgerblue", "very_low" = "navyblue"),
                    labels=c("very_high" = "Very high\n[23-100]", "high" = "High\n[11-100]", "low" = "Low\n[0-5[", "very_low" = "Very low\n[0-1[")) +
  scale_x_discrete(labels = c("very_high" = "Very high", "high" = "High",  "low" = "Low", "very_low" = "Very low"),
                   limits = c("very_high", "high", "low", "very_low")) +
  labs(title = "Mean mimicry rarity in risk areas and refuges", x = NULL, y = "Mean mimicry rarity", fill = "Regions") +
  theme_classic() +
  theme(title = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 5, 0)),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 15),
        axis.line = element_line(color = "black"))

print(Boxplot_HII_bins_ring.rarity)

par(mar = internal_margins)
dev.off()

# 10.2.5/ All four indices

library(gridExtra) # To plot gglots side-by-side

pdf(file = paste0("./graphs/Risks_refuges/Boxplot_HII_bins_all_indices.pdf"), height = 10, width = 12)

# gridExtra::grid.arrange(Boxplot_HII_bins_sp.richness, Boxplot_HII_bins_sp.rarity, Boxplot_HII_bins_MPD, Boxplot_HII_bins_ring.richness, ncol = 2)
gridExtra::grid.arrange(Boxplot_HII_bins_sp.richness, Boxplot_HII_bins_sp.rarity, Boxplot_HII_bins_ring.richness, Boxplot_HII_bins_ring.rarity, ncol = 2)

dev.off()
