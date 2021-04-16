###### Script 18: Cropped Threat Maps ######

# Crop continuous HF maps keeping only proportions of most diverse/rare communities

# Only for Jaccard.80

### Use only 4 indices

# 1/ Species richness
# 2/ Mean Continuous Species Rarity = Range-size Weighted Species Richness
# 3/ Mimicry richness
# 4/ Mean Continuous Mimicry rarity = Range-size Weighted Mimicry Richness
# Bonus/ Mean pairwise Phylogenetic Distance (removed from the final analyses)

###


# Inputs: 
   # Index maps from script 14a
   # Threat map: HII (HF, Venter et al., 2016)

# Outputs:
   # Mask of the top 5%, 10%, 25% highest scoring communities for all the indices
   # Cropped HF maps for the top 5%, 10%, 25% highest scoring communities for all the indices

###

# Effacer l'environnement
rm(list = ls())

library(raster)
library(rgeos)
library(rangeBuilder)

# Load data for map masks
most_diverse_rank_global_map_indices <- load(file = "./outputs/Threat_curves/Most_diverse/most_diverse_rank_global_map_indices.rds")

res <- "15"

# Load mask for continent borders
continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_", res, ".rds"))
crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".rds"))
bg_mask <- readRDS(file = "./input_data/Map_stuff/bg_mask.rds")
bg_mask_pixel <- readRDS(file = "./input_data/Map_stuff/bg_mask_pixel.rds")

# Load HII
HII <- readRDS(file = "./input_data/HII/HII.rds")
HII_contrasted <- readRDS(file = "./input_data/HII/HII_contrasted.rds")

# Load shp files for maps
load(file = "./input_data/Map_stuff/bg_mask.RData") # Load bg shp
load(file = "./input_data/Map_stuff/country_borders.RData") # Load country borders

# Load color palette
library(tmaptools)
pal_grn_red_NA <- pal_grn_red <- rev(get_brewer_pal("RdYlGn", n = 200))
pal_grn_red_NA[1] <- "#EDEDED"


### 1/ Species richness ####

# Generate masks
sp.richness_mask_rank_5 <- continent_mask
sp.richness_mask_rank_5[sp.richness_map_indices[[1]]] <- 1

sp.richness_mask_rank_10 <- continent_mask
sp.richness_mask_rank_10[sp.richness_map_indices[[2]]] <- 1
# Need a version without continental border to differentiate null value inside the area of interest from null value for continent
sp.richness_mask_rank_10_NA <- sp.richness_mask_rank_10
sp.richness_mask_rank_10_NA[sp.richness_mask_rank_10 == 0] <- NA

plot(sp.richness_mask_rank_10)
plot(sp.richness_mask_rank_10_NA)

sp.richness_mask_rank_25 <- continent_mask
sp.richness_mask_rank_25[sp.richness_map_indices[[3]]] <- 1
# Need a version without continental border to differentiate null value inside the area of interest from null value for continent
sp.richness_mask_rank_25_NA <- sp.richness_mask_rank_25
sp.richness_mask_rank_25_NA[sp.richness_mask_rank_25 == 0] <- NA

plot(sp.richness_mask_rank_25)
plot(sp.richness_mask_rank_25_NA)

# Save masks
save(sp.richness_mask_rank_5, sp.richness_mask_rank_10, sp.richness_mask_rank_10_NA, sp.richness_mask_rank_25, sp.richness_mask_rank_25_NA, file = "./outputs/Threat_maps/Cropped/sp.richness_masks_rank_files.RData", version = "2") 

# Generate shp files
sp.richness_mask_rank_5_shp <- rasterToPolygons(x = sp.richness_mask_rank_5, fun = function(x) {x == 1},
                                                dissolve = T)# To merge contiguous polygons of the same field value
sp.richness_mask_rank_10_shp <- rasterToPolygons(x = sp.richness_mask_rank_10, fun = function(x) {x == 1},
                                                 dissolve = T)# To merge contiguous polygons of the same field value
sp.richness_mask_rank_25_shp <- rasterToPolygons(x = sp.richness_mask_rank_25, fun = function(x) {x == 1},
                                                 dissolve = T)# To merge contiguous polygons of the same field value
sp.richness_mask_rank_25_minus_10_shp <- gDifference(sp.richness_mask_rank_25_shp, sp.richness_mask_rank_10_shp)
sp.richness_mask_rank_10_minus_5_shp <- gDifference(sp.richness_mask_rank_10_shp, sp.richness_mask_rank_5_shp)

sp.richness_mask_rank_5_shp_smooth <- gSimplify(sp.richness_mask_rank_5_shp, tol = 0.2, topologyPreserve = F)
# plot(sp.richness_mask_rank_5_shp_smooth)
sp.richness_mask_rank_25_shp_smooth <- gSimplify(sp.richness_mask_rank_25_shp, tol = 0.2, topologyPreserve = F)
# plot(sp.richness_mask_rank_25_shp_smooth)

# Save shp files
save(sp.richness_mask_rank_5_shp, sp.richness_mask_rank_5_shp_smooth, sp.richness_mask_rank_10_shp, sp.richness_mask_rank_25_shp, sp.richness_mask_rank_25_shp_smooth, sp.richness_mask_rank_25_minus_10_shp, sp.richness_mask_rank_10_minus_5_shp, file = "./outputs/Threat_maps/Cropped/sp.richness_masked_rank_shp_files.rds", version = "2") 

# Contrast the null value inside the area of interest
sp.richness_masked_rank_5 <- sp.richness_mask_rank_5*HII_contrasted

sp.richness_masked_rank_10_NA <- sp.richness_mask_rank_10_NA*HII_contrasted
sp.richness_masked_rank_10_NA[sp.richness_masked_rank_10_NA[] < 0.2] <- 0.2
sp.richness_masked_rank_10 <- continent_mask
sp.richness_masked_rank_10 <- merge(sp.richness_masked_rank_10_NA, continent_mask) 

sp.richness_masked_rank_25_NA <- sp.richness_mask_rank_25_NA*HII_contrasted
sp.richness_masked_rank_25_NA[sp.richness_masked_rank_25_NA[] < 0.2] <- 0.2
sp.richness_masked_rank_25 <- continent_mask
sp.richness_masked_rank_25 <- merge(sp.richness_masked_rank_25_NA, continent_mask) 

# Save cropped map files
saveRDS(sp.richness_masked_rank_5, file = "./outputs/Threat_maps/Cropped/sp.richness_masked_rank_5.rds", version = "2")
saveRDS(sp.richness_masked_rank_10, file = "./outputs/Threat_maps/Cropped/sp.richness_masked_rank_10.rds", version = "2")
saveRDS(sp.richness_masked_rank_25, file = "./outputs/Threat_maps/Cropped/sp.richness_masked_rank_25.rds", version = "2")

rgb(t(col2rgb("red")), maxColorValue = 255)

## 1.1/ Plot all areas ####
pdf(file = paste0("./maps/Threat_maps/Cropped/sp.richness_rank_all_areas.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.richness_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/10/25% for Species Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey80", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

# plot(sp.richness_mask_rank_25_minus_10_shp, lwd = 0.4, border = "#33333390", col = "#00000060", add = T)
# plot(sp.richness_mask_rank_10_minus_5_shp, lwd = 0.4, border = NA, col = "#00000030", add = T)
# plot(sp.richness_mask_rank_5_shp, lwd = 0.4, border = "#33333390", add = T)

# plot(sp.richness_mask_rank_25_shp, lwd = 0.4, border = "#33333390", add = T)
# plot(sp.richness_mask_rank_10_minus_5_shp, lwd = 0.4, border = "#33333390", col = "#00000045", add = T)
# plot(sp.richness_mask_rank_5_shp, lwd = 0.4, border = "#33333390", col = "#00000090", add = T)

plot(sp.richness_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)
plot(sp.richness_mask_rank_5_shp_smooth, lwd = 2, border = "mediumblue", col = NA, add = T)

# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mar = internal_margins)
dev.off()


## 1.2/ Plot area 5% ####
pdf(file = paste0("./maps/Threat_maps/Cropped/sp.richness_rank_5.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.richness_masked_rank_5, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5% for Species Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(sp.richness_mask_rank_5_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mar = internal_margins)
dev.off()

## 1.3/ Plot area 10% ####
pdf(file = paste0("./maps/Threat_maps/Cropped/sp.richness_rank_10.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.richness_masked_rank_10, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 10% for Species Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(sp.richness_mask_rank_10_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mar = internal_margins)
dev.off()

## 1.4/ Plot area 25% ####
pdf(file = paste0("./maps/Threat_maps/Cropped/sp.richness_rank_25.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.richness_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 25% for Species Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(sp.richness_mask_rank_25_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mar = internal_margins)
dev.off()


## 1.5/ Multiple plot ####
# Make a multiple plot with a plot with all full 3 areas, and 1 facet for each

pdf(file = paste0("./maps/Threat_maps/Cropped/sp.richness_rank.pdf"), height = 11, width = 13)

par(mfrow = c(2, 2))
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.richness_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/10/25% for Species Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)

plot(sp.richness_mask_rank_25_minus_10_shp, lwd = 0.4, border = "#33333390", col = "#00000060", add = T)
plot(sp.richness_mask_rank_10_minus_5_shp, lwd = 0.4, border = NA, col = "#00000030", add = T)
plot(sp.richness_mask_rank_5_shp, lwd = 0.4, border = "#33333390", add = T)

# plot(sp.richness_mask_rank_25_shp, lwd = 0.4, border = "#33333390", add = T)
# plot(sp.richness_mask_rank_10_minus_5_shp, lwd = 0.4, border = "#33333390", col = "#00000045", add = T)
# plot(sp.richness_mask_rank_5_shp, lwd = 0.4, border = "#33333390", col = "#00000090", add = T)

# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(sp.richness_masked_rank_5, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5% for Species Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(sp.richness_mask_rank_5_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(sp.richness_masked_rank_10, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 10% for Species Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(sp.richness_mask_rank_10_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(sp.richness_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 25% for Species Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(sp.richness_mask_rank_25_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mfrow = c(1, 1))
par(mar = internal_margins)
dev.off()


### 2/ Species rarity ####

# Generate masks
sp.rarity_mask_rank_5 <- continent_mask
sp.rarity_mask_rank_5[sp.rarity_map_indices[[1]]] <- 1

sp.rarity_mask_rank_10 <- continent_mask
sp.rarity_mask_rank_10[sp.rarity_map_indices[[2]]] <- 1
# Need a version without continental border to differentiate null value inside the area of interest from null value for contient
sp.rarity_mask_rank_10_NA <- sp.rarity_mask_rank_10
sp.rarity_mask_rank_10_NA[sp.rarity_mask_rank_10 == 0] <- NA

plot(sp.rarity_mask_rank_10)
plot(sp.rarity_mask_rank_10_NA)

sp.rarity_mask_rank_25 <- continent_mask
sp.rarity_mask_rank_25[sp.rarity_map_indices[[3]]] <- 1
# Need a version without continental border to differentiate null value inside the area of interest from null value for contient
sp.rarity_mask_rank_25_NA <- sp.rarity_mask_rank_25
sp.rarity_mask_rank_25_NA[sp.rarity_mask_rank_25 == 0] <- NA

plot(sp.rarity_mask_rank_25)
plot(sp.rarity_mask_rank_25_NA)

# Save masks
save(sp.rarity_mask_rank_5, sp.rarity_mask_rank_10, sp.rarity_mask_rank_10_NA, sp.rarity_mask_rank_25, sp.rarity_mask_rank_25_NA, file = "./outputs/Threat_maps/Cropped/sp.rarity_masks_rank_files.RData", version = "2")

# Generate shp files
sp.rarity_mask_rank_5_shp <- rasterToPolygons(x = sp.rarity_mask_rank_5, fun = function(x) {x == 1},
                                                dissolve = T)# To merge contiguous polygons of the same field value
sp.rarity_mask_rank_10_shp <- rasterToPolygons(x = sp.rarity_mask_rank_10, fun = function(x) {x == 1},
                                                 dissolve = T)# To merge contiguous polygons of the same field value
sp.rarity_mask_rank_25_shp <- rasterToPolygons(x = sp.rarity_mask_rank_25, fun = function(x) {x == 1},
                                                 dissolve = T)# To merge contiguous polygons of the same field value
sp.rarity_mask_rank_25_minus_10_shp <- gDifference(sp.rarity_mask_rank_25_shp, sp.rarity_mask_rank_10_shp)
sp.rarity_mask_rank_10_minus_5_shp <- gDifference(sp.rarity_mask_rank_10_shp, sp.rarity_mask_rank_5_shp)

sp.rarity_mask_rank_5_shp_smooth <- gSimplify(sp.rarity_mask_rank_5_shp, tol = 0.2, topologyPreserve = F)
# plot(sp.rarity_mask_rank_5_shp)
# plot(sp.rarity_mask_rank_5_shp_smooth)
sp.rarity_mask_rank_25_shp_smooth <- gSimplify(sp.rarity_mask_rank_25_shp, tol = 0.2, topologyPreserve = F)
# plot(sp.rarity_mask_rank_25_shp_smooth)

# Save shp files
save(sp.rarity_mask_rank_5_shp, sp.rarity_mask_rank_5_shp_smooth, sp.rarity_mask_rank_10_shp, sp.rarity_mask_rank_25_shp, sp.rarity_mask_rank_25_shp_smooth, sp.rarity_mask_rank_25_minus_10_shp, sp.rarity_mask_rank_10_minus_5_shp, file = "./outputs/Threat_maps/Cropped/sp.rarity_masked_rank_shp_files.rds", version = "2") 


# Contrast the null value inside the area of interest
sp.rarity_masked_rank_5 <- sp.rarity_mask_rank_5*HII_contrasted

sp.rarity_masked_rank_10_NA <- sp.rarity_mask_rank_10_NA*HII_contrasted
sp.rarity_masked_rank_10_NA[sp.rarity_masked_rank_10_NA[] < 0.2] <- 0.2
sp.rarity_masked_rank_10 <- continent_mask
sp.rarity_masked_rank_10 <- merge(sp.rarity_masked_rank_10_NA, continent_mask) 

sp.rarity_masked_rank_25_NA <- sp.rarity_mask_rank_25_NA*HII_contrasted
sp.rarity_masked_rank_25_NA[sp.rarity_masked_rank_25_NA[] < 0.2] <- 0.2
sp.rarity_masked_rank_25 <- continent_mask
sp.rarity_masked_rank_25 <- merge(sp.rarity_masked_rank_25_NA, continent_mask) 

# Save cropped map files
saveRDS(sp.rarity_masked_rank_5, file = "./outputs/Threat_maps/Cropped/sp.rarity_masked_rank_5.rds", version = "2")
saveRDS(sp.rarity_masked_rank_10, file = "./outputs/Threat_maps/Cropped/sp.rarity_masked_rank_10.rds", version = "2")
saveRDS(sp.rarity_masked_rank_25, file = "./outputs/Threat_maps/Cropped/sp.rarity_masked_rank_25.rds", version = "2")



## 2.1/ Plot all areas ####
pdf(file = paste0("./maps/Threat_maps/Cropped/sp.rarity_rank_all_areas.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.rarity_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/10/25% for Mean Species Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

# plot(sp.rarity_mask_rank_25_minus_10_shp, lwd = 0.4, border = "#33333390", col = "#00000060", add = T)
# plot(sp.rarity_mask_rank_10_minus_5_shp, lwd = 0.4, border = NA, col = "#00000030", add = T)
# plot(sp.rarity_mask_rank_5_shp, lwd = 0.4, border = "#33333390", add = T)

# plot(sp.rarity_mask_rank_25_shp, lwd = 0.4, border = "#33333390", add = T)
# plot(sp.rarity_mask_rank_10_minus_5_shp, lwd = 0.4, border = "#33333390", col = "#00000045", add = T)
# plot(sp.rarity_mask_rank_5_shp, lwd = 0.4, border = "#33333390", col = "#00000090", add = T)

plot(sp.rarity_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)
plot(sp.rarity_mask_rank_5_shp_smooth, lwd = 2, border = "mediumblue", col = NA, add = T)



# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mar = internal_margins)
dev.off()

## 2.2/ Plot area 5% ####
pdf(file = paste0("./maps/Threat_maps/Cropped/sp.rarity_rank_5.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.rarity_masked_rank_5, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5% for Mean Species Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(sp.rarity_mask_rank_5_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mar = internal_margins)
dev.off()

## 2.3/ Plot area 10% ####
pdf(file = paste0("./maps/Threat_maps/Cropped/sp.rarity_rank_10.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.rarity_masked_rank_10, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 10% for Mean Species Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(sp.rarity_mask_rank_10_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mar = internal_margins)
dev.off()

## 2.4/ Plot area 25% ####
pdf(file = paste0("./maps/Threat_maps/Cropped/sp.rarity_rank_25.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.rarity_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 25% for Mean Species Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(sp.rarity_mask_rank_25_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mar = internal_margins)
dev.off()


## 2.5/ Multiple plot ####
# Make a multiple plot with a plot with all full 3 areas, and 1 facet for each

pdf(file = paste0("./maps/Threat_maps/Cropped/sp.rarity_rank.pdf"), height = 11, width = 13)

par(mfrow = c(2, 2))
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.rarity_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/10/25% for Mean Species Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)

plot(sp.rarity_mask_rank_25_minus_10_shp, lwd = 0.4, border = "#33333390", col = "#00000060", add = T)
plot(sp.rarity_mask_rank_10_minus_5_shp, lwd = 0.4, border = NA, col = "#00000030", add = T)
plot(sp.rarity_mask_rank_5_shp, lwd = 0.4, border = "#33333390", add = T)

# plot(sp.rarity_mask_rank_25_shp, lwd = 0.4, border = "#33333390", add = T)
# plot(sp.rarity_mask_rank_10_minus_5_shp, lwd = 0.4, border = "#33333390", col = "#00000045", add = T)
# plot(sp.rarity_mask_rank_5_shp, lwd = 0.4, border = "#33333390", col = "#00000090", add = T)

# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(sp.rarity_masked_rank_5, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5% for Mean Species Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(sp.rarity_mask_rank_5_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(sp.rarity_masked_rank_10, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 10% for Mean Species Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(sp.rarity_mask_rank_10_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(sp.rarity_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 25% for Mean Species Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(sp.rarity_mask_rank_25_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mfrow = c(1, 1))
par(mar = internal_margins)
dev.off()

# ### 3/ MPD ####
# 
# # Generate masks
# MPD_mask_rank_5 <- continent_mask
# MPD_mask_rank_5[MPD_map_indices[[1]]] <- 1
# # Need a version without continental border to differentiate null value inside the area of interest from null value for contient
# MPD_mask_rank_5_NA <- MPD_mask_rank_5
# MPD_mask_rank_5_NA[MPD_mask_rank_5 == 0] <- NA
# 
# MPD_mask_rank_10 <- continent_mask
# MPD_mask_rank_10[MPD_map_indices[[2]]] <- 1
# # Need a version without continental border to differentiate null value inside the area of interest from null value for contient
# MPD_mask_rank_10_NA <- MPD_mask_rank_10
# MPD_mask_rank_10_NA[MPD_mask_rank_10 == 0] <- NA
# 
# plot(MPD_mask_rank_10)
# plot(MPD_mask_rank_10_NA)
# 
# MPD_mask_rank_25 <- continent_mask
# MPD_mask_rank_25[MPD_map_indices[[3]]] <- 1
# # Need a version without continental border to differentiate null value inside the area of interest from null value for contient
# MPD_mask_rank_25_NA <- MPD_mask_rank_25
# MPD_mask_rank_25_NA[MPD_mask_rank_25 == 0] <- NA
# 
# plot(MPD_mask_rank_25)
# plot(MPD_mask_rank_25_NA)
# 
# # Save masks
# save(MPD_mask_rank_5, MPD_mask_rank_10, MPD_mask_rank_10_NA, MPD_mask_rank_25, MPD_mask_rank_25_NA, file = "./outputs/Threat_maps/Cropped/MPD_masks_rank_files.RData", version = "2")
# 
# 
# # Generate shp files
# MPD_mask_rank_5_shp <- rasterToPolygons(x = MPD_mask_rank_5, fun = function(x) {x == 1},
#                                               dissolve = T)# To merge contiguous polygons of the same field value
# MPD_mask_rank_10_shp <- rasterToPolygons(x = MPD_mask_rank_10, fun = function(x) {x == 1},
#                                                dissolve = T)# To merge contiguous polygons of the same field value
# MPD_mask_rank_25_shp <- rasterToPolygons(x = MPD_mask_rank_25, fun = function(x) {x == 1},
#                                                dissolve = T)# To merge contiguous polygons of the same field value
# MPD_mask_rank_25_minus_10_shp <- gDifference(MPD_mask_rank_25_shp, MPD_mask_rank_10_shp)
# MPD_mask_rank_10_minus_5_shp <- gDifference(MPD_mask_rank_10_shp, MPD_mask_rank_5_shp)
# 
# MPD_mask_rank_5_shp_smooth <- gSimplify(MPD_mask_rank_5_shp, tol = 0.2, topologyPreserve = F)
# plot(MPD_mask_rank_5_shp_smooth)
# MPD_mask_rank_25_shp_smooth <- gSimplify(MPD_mask_rank_25_shp, tol = 0.2, topologyPreserve = F)
# plot(MPD_mask_rank_25_shp_smooth)
# 
# # Save shp files
# save(MPD_mask_rank_5_shp, MPD_mask_rank_5_shp_smooth, MPD_mask_rank_10_shp, MPD_mask_rank_25_shp, MPD_mask_rank_25_shp_smooth, MPD_mask_rank_25_minus_10_shp, MPD_mask_rank_10_minus_5_shp, file = "./outputs/Threat_maps/Cropped/MPD_masked_rank_shp_files.rds", version = "2") 
# 
# 
# # Contrast the null value inside the area of interest
# MPD_masked_rank_5_NA <- MPD_mask_rank_5_NA*HII_contrasted
# MPD_masked_rank_5_NA[MPD_masked_rank_5_NA[] < 0.2] <- 0.2
# MPD_masked_rank_5 <- continent_mask
# MPD_masked_rank_5 <- merge(MPD_masked_rank_5_NA, continent_mask)
# 
# MPD_masked_rank_10_NA <- MPD_mask_rank_10_NA*HII_contrasted
# MPD_masked_rank_10_NA[MPD_masked_rank_10_NA[] < 0.2] <- 0.2
# MPD_masked_rank_10 <- continent_mask
# MPD_masked_rank_10 <- merge(MPD_masked_rank_10_NA, continent_mask) 
# 
# MPD_masked_rank_25_NA <- MPD_mask_rank_25_NA*HII_contrasted
# MPD_masked_rank_25_NA[MPD_masked_rank_25_NA[] < 0.2] <- 0.2
# MPD_masked_rank_25 <- continent_mask
# MPD_masked_rank_25 <- merge(MPD_masked_rank_25_NA, continent_mask) 
# 
# # Save cropped map files
# saveRDS(MPD_masked_rank_5, file = "./outputs/Threat_maps/Cropped/MPD_masked_rank_5.rds", version = "2")
# saveRDS(MPD_masked_rank_10, file = "./outputs/Threat_maps/Cropped/MPD_masked_rank_10.rds", version = "2")
# saveRDS(MPD_masked_rank_25, file = "./outputs/Threat_maps/Cropped/MPD_masked_rank_25.rds", version = "2")
# 
# ## 3.1/ Plot all areas ####
# pdf(file = paste0("./maps/Threat_maps/Cropped/MPD_rank_all_areas.pdf"), height = 5.3, width = 6.5)
# 
# internal_margins <- par()$mar
# par(mar = c(3.1, 3.5,3.5,2.1))
# 
# image(MPD_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/10/25% for MPD"),
#       cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
#       legend  = T)
# plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)
# 
# # plot(MPD_mask_rank_25_minus_10_shp, lwd = 0.4, border = "#33333390", col = "#00000060", add = T)
# # plot(MPD_mask_rank_10_minus_5_shp, lwd = 0.4, border = NA, col = "#00000030", add = T)
# # plot(MPD_mask_rank_5_shp, lwd = 0.4, border = "#33333390", add = T)
# 
# # plot(MPD_mask_rank_25_shp, lwd = 0.4, border = "#33333390", add = T)
# # plot(MPD_mask_rank_10_minus_5_shp, lwd = 0.4, border = "#33333390", col = "#00000045", add = T)
# # plot(MPD_mask_rank_5_shp, lwd = 0.4, border = "#33333390", col = "#00000090", add = T)
# 
# plot(MPD_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)
# plot(MPD_mask_rank_5_shp_smooth, lwd = 2, border = "mediumblue", col = NA, add = T)
# 
# # abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# # abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# # abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
# addRasterLegend(MPD_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
# graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")
# 
# par(mar = internal_margins)
# dev.off()
# 
# ## 3.2/ Plot area 5% ####
# pdf(file = paste0("./maps/Threat_maps/Cropped/MPD_rank_5.pdf"), height = 5.3, width = 6.5)
# 
# internal_margins <- par()$mar
# par(mar = c(3.1, 3.5,3.5,2.1))
# 
# image(MPD_masked_rank_5, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5% for MPD"),
#       cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
#       legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(MPD_mask_rank_5_shp, lwd = 0.4, border = "grey20", add = T)
# # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# # abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# # abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# # abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
# addRasterLegend(MPD_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
# graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")
# 
# par(mar = internal_margins)
# dev.off()
# 
# ## 3.3/ Plot area 10% ####
# pdf(file = paste0("./maps/Threat_maps/Cropped/MPD_rank_10.pdf"), height = 5.3, width = 6.5)
# 
# internal_margins <- par()$mar
# par(mar = c(3.1, 3.5,3.5,2.1))
# 
# image(MPD_masked_rank_10, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 10% for MPD"),
#       cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
#       legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(MPD_mask_rank_10_shp, lwd = 0.4, border = "grey20", add = T)
# # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# # abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# # abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# # abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
# addRasterLegend(MPD_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
# graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")
# 
# par(mar = internal_margins)
# dev.off()
# 
# ## 3.4/ Plot area 25% ####
# pdf(file = paste0("./maps/Threat_maps/Cropped/MPD_rank_25.pdf"), height = 5.3, width = 6.5)
# 
# internal_margins <- par()$mar
# par(mar = c(3.1, 3.5,3.5,2.1))
# 
# image(MPD_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 25% for MPD"),
#       cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
#       legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(MPD_mask_rank_25_shp, lwd = 0.4, border = "grey20", add = T)
# # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# # abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# # abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# # abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
# addRasterLegend(MPD_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
# graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")
# 
# par(mar = internal_margins)
# dev.off()
# 
# 
# ## 3.5/ Multiple plot ####
# # Make a multiple plot with a plot with all full 3 areas, and 1 facet for each
# 
# pdf(file = paste0("./maps/Threat_maps/Cropped/MPD_rank.pdf"), height = 11, width = 13)
# 
# par(mfrow = c(2, 2))
# internal_margins <- par()$mar
# par(mar = c(3.1, 3.5,3.5,2.1))
# 
# image(MPD_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/10/25% for MPD"),
#       cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
#       legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# 
# plot(MPD_mask_rank_25_minus_10_shp, lwd = 0.4, border = "#33333390", col = "#00000060", add = T)
# plot(MPD_mask_rank_10_minus_5_shp, lwd = 0.4, border = NA, col = "#00000030", add = T)
# plot(MPD_mask_rank_5_shp, lwd = 0.4, border = "#33333390", add = T)
# 
# # plot(MPD_mask_rank_25_shp, lwd = 0.4, border = "#33333390", add = T)
# # plot(MPD_mask_rank_10_minus_5_shp, lwd = 0.4, border = "#33333390", col = "#00000045", add = T)
# # plot(MPD_mask_rank_5_shp, lwd = 0.4, border = "#33333390", col = "#00000090", add = T)
# 
# # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# # abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# # abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# # abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
# addRasterLegend(MPD_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
# graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")
# 
# image(MPD_masked_rank_5, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5% for MPD"),
#       cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
#       legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(MPD_mask_rank_5_shp, lwd = 0.4, border = "grey20", add = T)
# # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# # abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# # abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# # abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
# addRasterLegend(MPD_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
# graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")
# 
# image(MPD_masked_rank_10, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 10% for MPD"),
#       cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
#       legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(MPD_mask_rank_10_shp, lwd = 0.4, border = "grey20", add = T)
# # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# # abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# # abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# # abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
# addRasterLegend(MPD_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
# graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")
# 
# image(MPD_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 25% for MPD"),
#       cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
#       legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(MPD_mask_rank_25_shp, lwd = 0.4, border = "grey20", add = T)
# # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# # abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# # abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# # abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
# addRasterLegend(MPD_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
# graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")
# 
# par(mfrow = c(1, 1))
# par(mar = internal_margins)
# dev.off()


### 4/ Mimicry richness ####

# Generate masks
ring.richness_mask_rank_5 <- continent_mask
ring.richness_mask_rank_5[ring.richness_map_indices[[1]]] <- 1
# Need a version without continental border to differentiate null value inside the area of interest from null value for contient
ring.richness_mask_rank_5_NA <- ring.richness_mask_rank_5
ring.richness_mask_rank_5_NA[ring.richness_mask_rank_5 == 0] <- NA

ring.richness_mask_rank_10 <- continent_mask
ring.richness_mask_rank_10[ring.richness_map_indices[[2]]] <- 1
# Need a version without continental border to differentiate null value inside the area of interest from null value for contient
ring.richness_mask_rank_10_NA <- ring.richness_mask_rank_10
ring.richness_mask_rank_10_NA[ring.richness_mask_rank_10 == 0] <- NA

plot(ring.richness_mask_rank_10)
plot(ring.richness_mask_rank_10_NA)

ring.richness_mask_rank_25 <- continent_mask
ring.richness_mask_rank_25[ring.richness_map_indices[[3]]] <- 1
# Need a version without continental border to differentiate null value inside the area of interest from null value for contient
ring.richness_mask_rank_25_NA <- ring.richness_mask_rank_25
ring.richness_mask_rank_25_NA[ring.richness_mask_rank_25 == 0] <- NA

plot(ring.richness_mask_rank_25)
plot(ring.richness_mask_rank_25_NA)

# Save masks
save(ring.richness_mask_rank_5, ring.richness_mask_rank_10, ring.richness_mask_rank_10_NA, ring.richness_mask_rank_25, ring.richness_mask_rank_25_NA, file = "./outputs/Threat_maps/Cropped/ring.richness_masks_rank_files.RData", version = "2")


# Generate shp files
ring.richness_mask_rank_5_shp <- rasterToPolygons(x = ring.richness_mask_rank_5, fun = function(x) {x == 1},
                                        dissolve = T)# To merge contiguous polygons of the same field value
ring.richness_mask_rank_10_shp <- rasterToPolygons(x = ring.richness_mask_rank_10, fun = function(x) {x == 1},
                                         dissolve = T)# To merge contiguous polygons of the same field value
ring.richness_mask_rank_25_shp <- rasterToPolygons(x = ring.richness_mask_rank_25, fun = function(x) {x == 1},
                                         dissolve = T)# To merge contiguous polygons of the same field value
ring.richness_mask_rank_25_minus_10_shp <- gDifference(ring.richness_mask_rank_25_shp, ring.richness_mask_rank_10_shp)
ring.richness_mask_rank_10_minus_5_shp <- gDifference(ring.richness_mask_rank_10_shp, ring.richness_mask_rank_5_shp)

ring.richness_mask_rank_5_shp_smooth <- gSimplify(ring.richness_mask_rank_5_shp, tol = 0.2, topologyPreserve = F)
# plot(ring.richness_mask_rank_5_shp_smooth)
ring.richness_mask_rank_25_shp_smooth <- gSimplify(ring.richness_mask_rank_25_shp, tol = 0.2, topologyPreserve = F)
# plot(ring.richness_mask_rank_25_shp_smooth)

# Save shp files
save(ring.richness_mask_rank_5_shp, ring.richness_mask_rank_5_shp_smooth, ring.richness_mask_rank_10_shp, ring.richness_mask_rank_25_shp, ring.richness_mask_rank_25_shp_smooth, ring.richness_mask_rank_25_minus_10_shp, ring.richness_mask_rank_10_minus_5_shp, file = "./outputs/Threat_maps/Cropped/ring.richness_masked_rank_shp_files.rds", version = "2") 


# Contrast the null value inside the area of interest
ring.richness_masked_rank_5_NA <- ring.richness_mask_rank_5_NA*HII_contrasted
ring.richness_masked_rank_5_NA[ring.richness_masked_rank_5_NA[] < 0.2] <- 0.2
ring.richness_masked_rank_5 <- continent_mask
ring.richness_masked_rank_5 <- merge(ring.richness_masked_rank_5_NA, continent_mask)

ring.richness_masked_rank_10_NA <- ring.richness_mask_rank_10_NA*HII_contrasted
ring.richness_masked_rank_10_NA[ring.richness_masked_rank_10_NA[] < 0.2] <- 0.2
ring.richness_masked_rank_10 <- continent_mask
ring.richness_masked_rank_10 <- merge(ring.richness_masked_rank_10_NA, continent_mask) 

ring.richness_masked_rank_25_NA <- ring.richness_mask_rank_25_NA*HII_contrasted
ring.richness_masked_rank_25_NA[ring.richness_masked_rank_25_NA[] < 0.2] <- 0.2
ring.richness_masked_rank_25 <- continent_mask
ring.richness_masked_rank_25 <- merge(ring.richness_masked_rank_25_NA, continent_mask) 

# Save cropped map files
saveRDS(ring.richness_masked_rank_5, file = "./outputs/Threat_maps/Cropped/ring.richness_masked_rank_5.rds", version = "2")
saveRDS(ring.richness_masked_rank_10, file = "./outputs/Threat_maps/Cropped/ring.richness_masked_rank_10.rds", version = "2")
saveRDS(ring.richness_masked_rank_25, file = "./outputs/Threat_maps/Cropped/ring.richness_masked_rank_25.rds", version = "2")

## 4.1/ Plot all areas ####
pdf(file = paste0("./maps/Threat_maps/Cropped/ring.richness_rank_all_areas.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(ring.richness_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/10/25% for Mimicry Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

# plot(ring.richness_mask_rank_25_minus_10_shp, lwd = 0.4, border = "#33333390", col = "#00000060", add = T)
# plot(ring.richness_mask_rank_10_minus_5_shp, lwd = 0.4, border = NA, col = "#00000030", add = T)
# plot(ring.richness_mask_rank_5_shp, lwd = 0.4, border = "#33333390", add = T)

# plot(ring.richness_mask_rank_25_shp, lwd = 0.4, border = "#33333390", add = T)
# plot(ring.richness_mask_rank_10_minus_5_shp, lwd = 0.4, border = "#33333390", col = "#00000045", add = T)
# plot(ring.richness_mask_rank_5_shp, lwd = 0.4, border = "#33333390", col = "#00000090", add = T)

plot(ring.richness_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)
plot(ring.richness_mask_rank_5_shp_smooth, lwd = 2, border = "mediumblue", col = NA, add = T)

# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mar = internal_margins)
dev.off()

## 4.2/ Plot area 5% ####
pdf(file = paste0("./maps/Threat_maps/Cropped/ring.richness_rank_5.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(ring.richness_masked_rank_5, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5% for Mimicry Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(ring.richness_mask_rank_5_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mar = internal_margins)
dev.off()

## 4.3/ Plot area 10% ####
pdf(file = paste0("./maps/Threat_maps/Cropped/ring.richness_rank_10.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(ring.richness_masked_rank_10, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 10% for Mimicry Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(ring.richness_mask_rank_10_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mar = internal_margins)
dev.off()

## 4.4/ Plot area 25% ####
pdf(file = paste0("./maps/Threat_maps/Cropped/ring.richness_rank_25.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(ring.richness_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 25% for Mimicry Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(ring.richness_mask_rank_25_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mar = internal_margins)
dev.off()


## 4.5/ Multiple plot ####
# Make a multiple plot with a plot with all full 3 areas, and 1 facet for each

pdf(file = paste0("./maps/Threat_maps/Cropped/ring.richness_rank.pdf"), height = 11, width = 13)

par(mfrow = c(2, 2))
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(ring.richness_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/10/25% for Mimicry Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)

plot(ring.richness_mask_rank_25_minus_10_shp, lwd = 0.4, border = "#33333390", col = "#00000060", add = T)
plot(ring.richness_mask_rank_10_minus_5_shp, lwd = 0.4, border = NA, col = "#00000030", add = T)
plot(ring.richness_mask_rank_5_shp, lwd = 0.4, border = "#33333390", add = T)

# plot(ring.richness_mask_rank_25_shp, lwd = 0.4, border = "#33333390", add = T)
# plot(ring.richness_mask_rank_10_minus_5_shp, lwd = 0.4, border = "#33333390", col = "#00000045", add = T)
# plot(ring.richness_mask_rank_5_shp, lwd = 0.4, border = "#33333390", col = "#00000090", add = T)

# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(ring.richness_masked_rank_5, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5% for Mimicry Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(ring.richness_mask_rank_5_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(ring.richness_masked_rank_10, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 10% for Mimicry Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(ring.richness_mask_rank_10_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(ring.richness_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 25% for Mimicry Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(ring.richness_mask_rank_25_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mfrow = c(1, 1))
par(mar = internal_margins)
dev.off()


### 5/ Mimicry rarity ####

# Generate masks
ring.rarity_mask_rank_5 <- continent_mask
ring.rarity_mask_rank_5[ring.rarity_map_indices[[1]]] <- 1

ring.rarity_mask_rank_10 <- continent_mask
ring.rarity_mask_rank_10[ring.rarity_map_indices[[2]]] <- 1
# Need a version without continental border to differentiate null value inside the area of interest from null value for contient
ring.rarity_mask_rank_10_NA <- ring.rarity_mask_rank_10
ring.rarity_mask_rank_10_NA[ring.rarity_mask_rank_10 == 0] <- NA

plot(ring.rarity_mask_rank_10)
plot(ring.rarity_mask_rank_10_NA)

ring.rarity_mask_rank_25 <- continent_mask
ring.rarity_mask_rank_25[ring.rarity_map_indices[[3]]] <- 1
# Need a version without continental border to differentiate null value inside the area of interest from null value for contient
ring.rarity_mask_rank_25_NA <- ring.rarity_mask_rank_25
ring.rarity_mask_rank_25_NA[ring.rarity_mask_rank_25 == 0] <- NA

plot(ring.rarity_mask_rank_25)
plot(ring.rarity_mask_rank_25_NA)

# Save masks
save(ring.rarity_mask_rank_5, ring.rarity_mask_rank_10, ring.rarity_mask_rank_10_NA, ring.rarity_mask_rank_25, ring.rarity_mask_rank_25_NA, file = "./outputs/Threat_maps/Cropped/ring.rarity_masks_rank_files.RData", version = "2")

# Generate shp files
ring.rarity_mask_rank_5_shp <- rasterToPolygons(x = ring.rarity_mask_rank_5, fun = function(x) {x == 1},
                                              dissolve = T)# To merge contiguous polygons of the same field value
ring.rarity_mask_rank_10_shp <- rasterToPolygons(x = ring.rarity_mask_rank_10, fun = function(x) {x == 1},
                                               dissolve = T)# To merge contiguous polygons of the same field value
ring.rarity_mask_rank_25_shp <- rasterToPolygons(x = ring.rarity_mask_rank_25, fun = function(x) {x == 1},
                                               dissolve = T)# To merge contiguous polygons of the same field value
ring.rarity_mask_rank_25_minus_10_shp <- gDifference(ring.rarity_mask_rank_25_shp, ring.rarity_mask_rank_10_shp)
ring.rarity_mask_rank_10_minus_5_shp <- gDifference(ring.rarity_mask_rank_10_shp, ring.rarity_mask_rank_5_shp)

ring.rarity_mask_rank_5_shp_smooth <- gSimplify(ring.rarity_mask_rank_5_shp, tol = 0.2, topologyPreserve = F)
# plot(ring.rarity_mask_rank_5_shp)
# plot(ring.rarity_mask_rank_5_shp_smooth)
ring.rarity_mask_rank_25_shp_smooth <- gSimplify(ring.rarity_mask_rank_25_shp, tol = 0.2, topologyPreserve = F)
# plot(ring.rarity_mask_rank_25_shp)
# plot(ring.rarity_mask_rank_25_shp_smooth)

# Save shp files
save(ring.rarity_mask_rank_5_shp, ring.rarity_mask_rank_5_shp_smooth, ring.rarity_mask_rank_10_shp, ring.rarity_mask_rank_25_shp, ring.rarity_mask_rank_25_shp_smooth, ring.rarity_mask_rank_25_minus_10_shp, ring.rarity_mask_rank_10_minus_5_shp, file = "./outputs/Threat_maps/Cropped/ring.rarity_masked_rank_shp_files.rds", version = "2") 


# Contrast the null value inside the area of interest
ring.rarity_masked_rank_5 <- ring.rarity_mask_rank_5*HII_contrasted

ring.rarity_masked_rank_10_NA <- ring.rarity_mask_rank_10_NA*HII_contrasted
ring.rarity_masked_rank_10_NA[ring.rarity_masked_rank_10_NA[] < 0.2] <- 0.2
ring.rarity_masked_rank_10 <- continent_mask
ring.rarity_masked_rank_10 <- merge(ring.rarity_masked_rank_10_NA, continent_mask) 

ring.rarity_masked_rank_25_NA <- ring.rarity_mask_rank_25_NA*HII_contrasted
ring.rarity_masked_rank_25_NA[ring.rarity_masked_rank_25_NA[] < 0.2] <- 0.2
ring.rarity_masked_rank_25 <- continent_mask
ring.rarity_masked_rank_25 <- merge(ring.rarity_masked_rank_25_NA, continent_mask) 

# Save cropped map files
saveRDS(ring.rarity_masked_rank_5, file = "./outputs/Threat_maps/Cropped/ring.rarity_masked_rank_5.rds", version = "2")
saveRDS(ring.rarity_masked_rank_10, file = "./outputs/Threat_maps/Cropped/ring.rarity_masked_rank_10.rds", version = "2")
saveRDS(ring.rarity_masked_rank_25, file = "./outputs/Threat_maps/Cropped/ring.rarity_masked_rank_25.rds", version = "2")



## 5.1/ Plot all areas ####
pdf(file = paste0("./maps/Threat_maps/Cropped/ring.rarity_rank_all_areas.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(ring.rarity_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/10/25% for Mean Mimicry Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

# plot(ring.rarity_mask_rank_25_minus_10_shp, lwd = 0.4, border = "#33333390", col = "#00000060", add = T)
# plot(ring.rarity_mask_rank_10_minus_5_shp, lwd = 0.4, border = NA, col = "#00000030", add = T)
# plot(ring.rarity_mask_rank_5_shp, lwd = 0.4, border = "#33333390", add = T)

# plot(ring.rarity_mask_rank_25_shp, lwd = 0.4, border = "#33333390", add = T)
# plot(ring.rarity_mask_rank_10_minus_5_shp, lwd = 0.4, border = "#33333390", col = "#00000045", add = T)
# plot(ring.rarity_mask_rank_5_shp, lwd = 0.4, border = "#33333390", col = "#00000090", add = T)

plot(ring.rarity_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)
plot(ring.rarity_mask_rank_5_shp_smooth, lwd = 2, border = "mediumblue", col = NA, add = T)



# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mar = internal_margins)
dev.off()

## 5.2/ Plot area 5% ####
pdf(file = paste0("./maps/Threat_maps/Cropped/ring.rarity_rank_5.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(ring.rarity_masked_rank_5, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5% for Mean Mimicry Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(ring.rarity_mask_rank_5_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mar = internal_margins)
dev.off()

## 5.3/ Plot area 10% ####
pdf(file = paste0("./maps/Threat_maps/Cropped/ring.rarity_rank_10.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(ring.rarity_masked_rank_10, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 10% for Mean Mimicry Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(ring.rarity_mask_rank_10_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mar = internal_margins)
dev.off()

## 5.4/ Plot area 25% ####
pdf(file = paste0("./maps/Threat_maps/Cropped/ring.rarity_rank_25.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(ring.rarity_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 25% for Mean Mimicry Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(ring.rarity_mask_rank_25_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mar = internal_margins)
dev.off()


## 5.5/ Multiple plot ####
# Make a multiple plot with a plot with all full 3 areas, and 1 facet for each

pdf(file = paste0("./maps/Threat_maps/Cropped/ring.rarity_rank.pdf"), height = 11, width = 13)

par(mfrow = c(2, 2))
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(ring.rarity_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/10/25% for Mean Mimicry Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)

plot(ring.rarity_mask_rank_25_minus_10_shp, lwd = 0.4, border = "#33333390", col = "#00000060", add = T)
plot(ring.rarity_mask_rank_10_minus_5_shp, lwd = 0.4, border = NA, col = "#00000030", add = T)
plot(ring.rarity_mask_rank_5_shp, lwd = 0.4, border = "#33333390", add = T)

# plot(ring.rarity_mask_rank_25_shp, lwd = 0.4, border = "#33333390", add = T)
# plot(ring.rarity_mask_rank_10_minus_5_shp, lwd = 0.4, border = "#33333390", col = "#00000045", add = T)
# plot(ring.rarity_mask_rank_5_shp, lwd = 0.4, border = "#33333390", col = "#00000090", add = T)

# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(ring.rarity_masked_rank_5, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5% for Mean Mimicry Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(ring.rarity_mask_rank_5_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(ring.rarity_masked_rank_10, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 10% for Mean Mimicry Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(ring.rarity_mask_rank_10_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(ring.rarity_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 25% for Mean Mimicry Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(ring.rarity_mask_rank_25_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mfrow = c(1, 1))
par(mar = internal_margins)
dev.off()



### 6/ Final plots with all 4 indices and 3 areas ####

# Load maps
sp.richness_masked_rank_5 <- readRDS(file = "./outputs/Threat_maps/Cropped/sp.richness_masked_rank_5.rds")
sp.richness_masked_rank_10 <- readRDS(file = "./outputs/Threat_maps/Cropped/sp.richness_masked_rank_10.rds")
sp.richness_masked_rank_25 <- readRDS(file = "./outputs/Threat_maps/Cropped/sp.richness_masked_rank_25.rds")

sp.rarity_masked_rank_5 <- readRDS(file = "./outputs/Threat_maps/Cropped/sp.rarity_masked_rank_5.rds")
sp.rarity_masked_rank_10 <- readRDS(file = "./outputs/Threat_maps/Cropped/sp.rarity_masked_rank_10.rds")
sp.rarity_masked_rank_25 <- readRDS(file = "./outputs/Threat_maps/Cropped/sp.rarity_masked_rank_25.rds")

# MPD_masked_rank_5 <- readRDS(file = "./outputs/Threat_maps/Cropped/MPD_masked_rank_5.rds")
# MPD_masked_rank_10 <- readRDS(file = "./outputs/Threat_maps/Cropped/MPD_masked_rank_10.rds")
# MPD_masked_rank_25 <- readRDS(file = "./outputs/Threat_maps/Cropped/MPD_masked_rank_25.rds")

ring.richness_masked_rank_5 <- readRDS(file = "./outputs/Threat_maps/Cropped/ring.richness_masked_rank_5.rds")
ring.richness_masked_rank_10 <- readRDS(file = "./outputs/Threat_maps/Cropped/ring.richness_masked_rank_10.rds")
ring.richness_masked_rank_25 <- readRDS(file = "./outputs/Threat_maps/Cropped/ring.richness_masked_rank_25.rds")

ring.rarity_masked_rank_5 <- readRDS(file = "./outputs/Threat_maps/Cropped/ring.rarity_masked_rank_5.rds")
ring.rarity_masked_rank_10 <- readRDS(file = "./outputs/Threat_maps/Cropped/ring.rarity_masked_rank_10.rds")
ring.rarity_masked_rank_25 <- readRDS(file = "./outputs/Threat_maps/Cropped/ring.rarity_masked_rank_25.rds")

# Load shp files
load(file = "./outputs/Threat_maps/Cropped/sp.richness_masked_rank_shp_files.rds")
load(file = "./outputs/Threat_maps/Cropped/sp.rarity_masked_rank_shp_files.rds")
# load(file = "./outputs/Threat_maps/Cropped/MPD_masked_rank_shp_files.rds")
load(file = "./outputs/Threat_maps/Cropped/ring.richness_masked_rank_shp_files.rds")
load(file = "./outputs/Threat_maps/Cropped/ring.rarity_masked_rank_shp_files.rds")

# Plot a multiple pdf


pdf(file = paste0("./maps/Threat_maps/Cropped/global_rank_all_indices_all_areas.pdf"), height = 11, width = 13)

# Page 1: all areas ####

par(mfrow = c(2, 2))
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(sp.richness_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/10/25% for Species Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)

plot(sp.richness_mask_rank_25_minus_10_shp, lwd = 0.4, border = "#33333390", col = "#00000060", add = T)
plot(sp.richness_mask_rank_10_minus_5_shp, lwd = 0.4, border = NA, col = "#00000030", add = T)
plot(sp.richness_mask_rank_5_shp, lwd = 0.4, border = "#33333390", add = T)

# plot(sp.richness_mask_rank_25_shp, lwd = 0.4, border = "#33333390", add = T)
# plot(sp.richness_mask_rank_10_minus_5_shp, lwd = 0.4, border = "#33333390", col = "#00000045", add = T)
# plot(sp.richness_mask_rank_5_shp, lwd = 0.4, border = "#33333390", col = "#00000090", add = T)

# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")


image(sp.rarity_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/10/25% for Mean Species Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)

plot(sp.rarity_mask_rank_25_minus_10_shp, lwd = 0.4, border = "#33333390", col = "#00000060", add = T)
plot(sp.rarity_mask_rank_10_minus_5_shp, lwd = 0.4, border = NA, col = "#00000030", add = T)
plot(sp.rarity_mask_rank_5_shp, lwd = 0.4, border = "#33333390", add = T)

# plot(sp.rarity_mask_rank_25_shp, lwd = 0.4, border = "#33333390", add = T)
# plot(sp.rarity_mask_rank_10_minus_5_shp, lwd = 0.4, border = "#33333390", col = "#00000045", add = T)
# plot(sp.rarity_mask_rank_5_shp, lwd = 0.4, border = "#33333390", col = "#00000090", add = T)

# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")


# image(MPD_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/10/25% for MPD"),
#       cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
#       legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# 
# plot(MPD_mask_rank_25_minus_10_shp, lwd = 0.4, border = "#33333390", col = "#00000060", add = T)
# plot(MPD_mask_rank_10_minus_5_shp, lwd = 0.4, border = NA, col = "#00000030", add = T)
# plot(MPD_mask_rank_5_shp, lwd = 0.4, border = "#33333390", add = T)
# 
# # plot(MPD_mask_rank_25_shp, lwd = 0.4, border = "#33333390", add = T)
# # plot(MPD_mask_rank_10_minus_5_shp, lwd = 0.4, border = "#33333390", col = "#00000045", add = T)
# # plot(MPD_mask_rank_5_shp, lwd = 0.4, border = "#33333390", col = "#00000090", add = T)
# 
# # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# # abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# # abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# # abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
# addRasterLegend(MPD_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
# graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")


image(ring.richness_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/10/25% for Mimicry Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)

plot(ring.richness_mask_rank_25_minus_10_shp, lwd = 0.4, border = "#33333390", col = "#00000060", add = T)
plot(ring.richness_mask_rank_10_minus_5_shp, lwd = 0.4, border = NA, col = "#00000030", add = T)
plot(ring.richness_mask_rank_5_shp, lwd = 0.4, border = "#33333390", add = T)

# plot(ring.richness_mask_rank_25_shp, lwd = 0.4, border = "#33333390", add = T)
# plot(ring.richness_mask_rank_10_minus_5_shp, lwd = 0.4, border = "#33333390", col = "#00000045", add = T)
# plot(ring.richness_mask_rank_5_shp, lwd = 0.4, border = "#33333390", col = "#00000090", add = T)

# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")


image(ring.rarity_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/10/25% for Mean Mimicry Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)

plot(ring.rarity_mask_rank_25_minus_10_shp, lwd = 0.4, border = "#33333390", col = "#00000060", add = T)
plot(ring.rarity_mask_rank_10_minus_5_shp, lwd = 0.4, border = NA, col = "#00000030", add = T)
plot(ring.rarity_mask_rank_5_shp, lwd = 0.4, border = "#33333390", add = T)

# plot(ring.rarity_mask_rank_25_shp, lwd = 0.4, border = "#33333390", add = T)
# plot(ring.rarity_mask_rank_10_minus_5_shp, lwd = 0.4, border = "#33333390", col = "#00000045", add = T)
# plot(ring.rarity_mask_rank_5_shp, lwd = 0.4, border = "#33333390", col = "#00000090", add = T)

# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")


# Page 2: Top 5% areas ####

image(sp.richness_masked_rank_5, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5% for Species Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(sp.richness_mask_rank_5_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(sp.rarity_masked_rank_5, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5% for Mean Species Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(sp.rarity_mask_rank_5_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

# image(MPD_masked_rank_5, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5% for MPD"),
#       cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
#       legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(MPD_mask_rank_5_shp, lwd = 0.4, border = "grey20", add = T)
# # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# # abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# # abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# # abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
# addRasterLegend(MPD_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
# graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(ring.richness_masked_rank_5, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5% for Mimicry Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(ring.richness_mask_rank_5_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(ring.rarity_masked_rank_5, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5% for Mean Mimicry Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(ring.rarity_mask_rank_5_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")


# Page 3: Top 10% areas ####

image(sp.richness_masked_rank_10, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 10% for Species Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(sp.richness_mask_rank_10_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(sp.rarity_masked_rank_10, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 10% for Mean Species Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(sp.rarity_mask_rank_10_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

# image(MPD_masked_rank_10, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 10% for MPD"),
#       cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
#       legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(MPD_mask_rank_10_shp, lwd = 0.4, border = "grey20", add = T)
# # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# # abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# # abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# # abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
# addRasterLegend(MPD_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
# graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(ring.richness_masked_rank_10, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 10% for Mimicry Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(ring.richness_mask_rank_10_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(ring.rarity_masked_rank_10, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 10% for Mean Mimicry Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(ring.rarity_mask_rank_10_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")


# Page 4: Top 25% areas ####

image(sp.richness_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 25% for Species Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(sp.richness_mask_rank_25_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(sp.rarity_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 25% for Mean Species Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(sp.rarity_mask_rank_25_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

# image(MPD_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 25% for MPD"),
#       cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
#       legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(MPD_mask_rank_25_shp, lwd = 0.4, border = "grey20", add = T)
# # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# # abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# # abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# # abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
# addRasterLegend(MPD_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
# graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(ring.richness_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 25% for Mimicry Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(ring.richness_mask_rank_25_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

image(ring.rarity_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 25% for Mean Mimicry Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(ring.rarity_mask_rank_25_shp, lwd = 0.4, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

par(mfrow = c(1, 1))
par(mar = internal_margins)
dev.off()



### 7/ Final plots with all 4 indices and 2 areas ####

# Load maps
sp.richness_masked_rank_5 <- readRDS(file = "./outputs/Threat_maps/Cropped/sp.richness_masked_rank_5.rds")
sp.richness_masked_rank_25 <- readRDS(file = "./outputs/Threat_maps/Cropped/sp.richness_masked_rank_25.rds")

sp.rarity_masked_rank_5 <- readRDS(file = "./outputs/Threat_maps/Cropped/sp.rarity_masked_rank_5.rds")
sp.rarity_masked_rank_25 <- readRDS(file = "./outputs/Threat_maps/Cropped/sp.rarity_masked_rank_25.rds")

# MPD_masked_rank_5 <- readRDS(file = "./outputs/Threat_maps/Cropped/MPD_masked_rank_5.rds")
# MPD_masked_rank_25 <- readRDS(file = "./outputs/Threat_maps/Cropped/MPD_masked_rank_25.rds")

ring.richness_masked_rank_5 <- readRDS(file = "./outputs/Threat_maps/Cropped/ring.richness_masked_rank_5.rds")
ring.richness_masked_rank_25 <- readRDS(file = "./outputs/Threat_maps/Cropped/ring.richness_masked_rank_25.rds")

ring.rarity_masked_rank_5 <- readRDS(file = "./outputs/Threat_maps/Cropped/ring.rarity_masked_rank_5.rds")
ring.rarity_masked_rank_25 <- readRDS(file = "./outputs/Threat_maps/Cropped/ring.rarity_masked_rank_25.rds")


# Load shp files
load(file = "./outputs/Threat_maps/Cropped/sp.richness_masked_rank_shp_files.rds")
load(file = "./outputs/Threat_maps/Cropped/sp.rarity_masked_rank_shp_files.rds")
# load(file = "./outputs/Threat_maps/Cropped/MPD_masked_rank_shp_files.rds")
load(file = "./outputs/Threat_maps/Cropped/ring.richness_masked_rank_shp_files.rds")
load(file = "./outputs/Threat_maps/Cropped/ring.rarity_masked_rank_shp_files.rds")

# Plot a multiple pdf


pdf(file = paste0("./maps/Threat_maps/Cropped/global_rank_all_indices_5_25_areas.pdf"), height = 11, width = 13)

par(mfrow = c(2, 2))
internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

# Species richness
image(sp.richness_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/25% for Species Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey80", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.richness_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)
plot(sp.richness_mask_rank_5_shp_smooth, lwd = 2, border = "mediumblue", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")


# Species rarity
image(sp.rarity_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/25% for Mean Species Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.rarity_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)
plot(sp.rarity_mask_rank_5_shp_smooth, lwd = 2, border = "mediumblue", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")


# # MPD
# image(MPD_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/25% for MPD"),
#       cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
#       legend  = T)
# plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)
# 
# plot(MPD_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)
# plot(MPD_mask_rank_5_shp_smooth, lwd = 2, border = "mediumblue", col = NA, add = T)
# 
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
# addRasterLegend(MPD_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
# graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")


# Mimicry richness
image(ring.richness_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/25% for Mimicry Richness"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.richness_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)
plot(ring.richness_mask_rank_5_shp_smooth, lwd = 2, border = "mediumblue", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.richness_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")

# Mimicry rarity
image(ring.rarity_masked_rank_25, col = pal_grn_red_NA, main = paste0("Human Influence Index\nTop 5/25% for Mean Mimicry Rarity"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.rarity_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)
plot(ring.rarity_mask_rank_5_shp_smooth, lwd = 2, border = "mediumblue", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.rarity_masked_rank_25, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")


par(mfrow = c(1, 1))
par(mar = internal_margins)
dev.off()


### 8/ Final plot for article = 4 indices + 2 boxplot (5/25%) ####

library(ggplot2)

# Load maps
sp.richness_masked_rank_5 <- readRDS(file = "./outputs/Threat_maps/Cropped/sp.richness_masked_rank_5.rds")
sp.richness_masked_rank_25 <- readRDS(file = "./outputs/Threat_maps/Cropped/sp.richness_masked_rank_25.rds")

sp.rarity_masked_rank_5 <- readRDS(file = "./outputs/Threat_maps/Cropped/sp.rarity_masked_rank_5.rds")
sp.rarity_masked_rank_25 <- readRDS(file = "./outputs/Threat_maps/Cropped/sp.rarity_masked_rank_25.rds")

# MPD_masked_rank_5 <- readRDS(file = "./outputs/Threat_maps/Cropped/MPD_masked_rank_5.rds")
# MPD_masked_rank_25 <- readRDS(file = "./outputs/Threat_maps/Cropped/MPD_masked_rank_25.rds")

ring.richness_masked_rank_5 <- readRDS(file = "./outputs/Threat_maps/Cropped/ring.richness_masked_rank_5.rds")
ring.richness_masked_rank_25 <- readRDS(file = "./outputs/Threat_maps/Cropped/ring.richness_masked_rank_25.rds")

ring.rarity_masked_rank_5 <- readRDS(file = "./outputs/Threat_maps/Cropped/ring.rarity_masked_rank_5.rds")
ring.rarity_masked_rank_25 <- readRDS(file = "./outputs/Threat_maps/Cropped/ring.rarity_masked_rank_25.rds")


# Load boxplots
most_diverse_rank_global_boxplot_5_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_rank_global_boxplot_5_df.rds")
most_diverse_rank_global_boxplot_25_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_rank_global_boxplot_25_df.rds")

# Load shp files
load(file = "./outputs/Threat_maps/Cropped/sp.richness_masked_rank_shp_files.rds")
load(file = "./outputs/Threat_maps/Cropped/sp.rarity_masked_rank_shp_files.rds")
# load(file = "./outputs/Threat_maps/Cropped/MPD_masked_rank_shp_files.rds")
load(file = "./outputs/Threat_maps/Cropped/ring.richness_masked_rank_shp_files.rds")
load(file = "./outputs/Threat_maps/Cropped/ring.rarity_masked_rank_shp_files.rds")

# Plot a multiple pdf

pdf(file = paste0("./for_article/threat_maps_2.pdf"), height = 11, width = 8)

par(mfrow = c(3, 2))
internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

# A/ Species richness
image(sp.richness_masked_rank_25, col = pal_grn_red_NA, main = paste0("Species richness"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.richness_mask_rank_25_shp_smooth, lwd = 1.2, border = "dodgerblue", add = T)
plot(sp.richness_mask_rank_5_shp_smooth, lwd = 1.5, border = "navyblue", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(sp.richness_masked_rank_25, locs = seq(0, 30, 10), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -110, y = 1.5, font = 2, cex = 1.2, label = "HII (%)")
text(x = - 100.7, y = -3, cex = 1.2, labels = "-100")


legend(legend = "A", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))
legend(legend = "Top 5%", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.04, 0.17))
legend(legend = "Top 25%", x = "topleft", bty = "n",
       pch = 22, col = "dodgerblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.04, 0.25))

# B/ Species rarity
image(sp.rarity_masked_rank_25, col = pal_grn_red_NA, main = paste0("Mean species rarity"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex = 1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.rarity_mask_rank_25_shp_smooth, lwd = 1.2, border = "dodgerblue", add = T)
plot(sp.rarity_mask_rank_5_shp_smooth, lwd = 1.5, border = "navyblue", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(sp.rarity_masked_rank_25, locs = seq(0, 30, 10), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -110, y = 1.5, font = 2, cex = 1.2, label = "HII (%)")
text(x = - 100.7, y = -3, cex = 1.2, labels = "-100")

legend(legend = "B", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))
legend(legend = "Top 5%", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.04, 0.17))
legend(legend = "Top 25%", x = "topleft", bty = "n",
       pch = 22, col = "dodgerblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.04, 0.25))

# # C/ MPD
# image(MPD_masked_rank_25, col = pal_grn_red_NA, main = paste0("Mean pairwise Phylogenetic Distance"),
#       cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
#       legend  = T)
# plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)
# 
# plot(MPD_mask_rank_25_shp_smooth, lwd = 1.2, border = "dodgerblue", add = T)
# plot(MPD_mask_rank_5_shp_smooth, lwd = 1.5, border = "navyblue", col = NA, add = T)
# 
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
# addRasterLegend(MPD_masked_rank_25, locs = seq(0, 30, 10), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
# graphics::text(x = -110, y = 1.5, font = 2, cex = 1.2, label = "HII (%)")
# text(x = - 100.7, y = -3, cex = 1.2, labels = "-100")
# 
# legend(legend = "C", x = "bottomright", bty = "n",
#        text.font = 2, cex = 1.8, inset=c(0, 0))
# legend(legend = "Top 5%", x = "topleft", bty = "n",
#        pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
#        text.font = 2, cex = 1.2, inset=c(0.04, 0.17))
# legend(legend = "Top 25%", x = "topleft", bty = "n",
#        pch = 22, col = "dodgerblue", pt.cex = 2, pt.lwd = 2.5,
#        text.font = 2, cex = 1.2, inset=c(0.04, 0.25))

# C/ Mimicry richness
image(ring.richness_masked_rank_25, col = pal_grn_red_NA, main = paste0("Mimicry richness"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.richness_mask_rank_25_shp_smooth, lwd = 1.2, border = "dodgerblue", add = T)
plot(ring.richness_mask_rank_5_shp_smooth, lwd = 1.5, border = "navyblue", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(ring.richness_masked_rank_25, locs = seq(0, 30, 10), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -110, y = 1.5, font = 2, cex = 1.2, label = "HII (%)")
text(x = - 100.7, y = -3, cex = 1.2, labels = "-100")

legend(legend = "C", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))
legend(legend = "Top 5%", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.04, 0.17))
legend(legend = "Top 25%", x = "topleft", bty = "n",
       pch = 22, col = "dodgerblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.04, 0.25))

# D/ Mimicry rarity
image(ring.rarity_masked_rank_25, col = pal_grn_red_NA, main = paste0("Mean mimicry rarity"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     HII", cex = 1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.rarity_mask_rank_25_shp_smooth, lwd = 1.2, border = "dodgerblue", add = T)
plot(ring.rarity_mask_rank_5_shp_smooth, lwd = 1.5, border = "navyblue", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))
addRasterLegend(ring.rarity_masked_rank_25, locs = seq(0, 30, 10), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, -3))
graphics::text(x = -110, y = 1.5, font = 2, cex = 1.2, label = "HII (%)")
text(x = - 100.7, y = -3, cex = 1.2, labels = "-100")

legend(legend = "D", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))
legend(legend = "Top 5%", x = "topleft", bty = "n",
       pch = 22, col = "navyblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.04, 0.17))
legend(legend = "Top 25%", x = "topleft", bty = "n",
       pch = 22, col = "dodgerblue", pt.cex = 2, pt.lwd = 2.5,
       text.font = 2, cex = 1.2, inset=c(0.04, 0.25))


# Boxplots graph par
par(mar = c(4.1,5.5,3.1,2.1))

# E/ Boxplot for Top 5%

# Draw the boxplot
boxplot(formula = most_diverse_rank_global_boxplot_5_df$HII ~ most_diverse_rank_global_boxplot_5_df$index,
        col = c("red", "dodgerblue", "darkorange", "limegreen"),
        ylab = "Human Influence Index  [%]",
        xlab = NA,
        main = "Top 5%", cex.main = 1.6,
        pch = 16,
        axes = F,
        cex.axis = 1.3,
        cex.lab = 1.5)
# Add the grid
grid(lty = 1) 
# Redraw the boxplot on the grid
boxplot(formula = most_diverse_rank_global_boxplot_5_df$HII ~ most_diverse_rank_global_boxplot_5_df$index,
        col = c("red", "dodgerblue", "darkorange", "limegreen"),
        ylab = "Human Influence Index  [%]",
        xlab = NA,
        main = "Top 5%", cex.main = 1.6,
        pch = 16,
        names = NA,
        cex.axis = 1.3,
        add = T)
# Draw x-labels with different size
# axis(1, cex.axis = 1.05, labels = c("Sp. richness", "Sp. rarity", "MPD", "Mim. richness"), font = 2, at = 1:4)
mtext(side = 1, cex = 0.73, 
      # text = c("Sp. richness", "Sp. rarity", "MPD", "Mim. richness"),
      text = c("Sp. rich.", "Sp. rarity", "Mim. rich.", "Mim. rarity"),
      font = 2, at = 1:4, line = 1)
# Add letter
legend(legend = "E", x = "topright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))

# F/ Boxplot for Top 25%

# Draw the boxplot
boxplot(formula = most_diverse_rank_global_boxplot_25_df$HII ~ most_diverse_rank_global_boxplot_25_df$index,
        col = c("red", "dodgerblue", "darkorange", "limegreen"),
        ylab = "Human Influence Index  [%]",
        xlab = NA,
        main = "Top 25%", cex.main = 1.6,
        pch = 16,
        axes = F,
        cex.axis = 1.3,
        cex.lab = 1.5)
# Add the grid
grid(lty = 1) 
# Redraw the boxplot on the grid
boxplot(formula = most_diverse_rank_global_boxplot_25_df$HII ~ most_diverse_rank_global_boxplot_25_df$index,
        col = c("red", "dodgerblue", "darkorange", "limegreen"),
        ylab = "Human Influence Index  [%]",
        xlab = NA,
        main = "Top 25%", cex.main = 1.6,
        pch = 16,
        names = NA,
        cex.axis = 1.3,
        add = T)
# Draw x-labels with different size
# axis(1, cex.axis = 1.05, labels = c("Sp. richness", "Sp. rarity", "MPD", "Mim. richness"), font = 2, at = 1:4)
mtext(side = 1, cex = 0.73, 
      # text = c("Sp. richness", "Sp. rarity", "MPD", "Mim. richness"),
      text = c("Sp. rich.", "Sp. rarity", "Mim. rich.", "Mim. rarity"),
      font = 2, at = 1:4, line = 1)

# Add letter
legend(legend = "F", x = "topright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))

par(mfrow = c(1, 1))
par(mar = internal_margins)
dev.off()

