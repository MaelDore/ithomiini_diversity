####################### Index maps for paper ############################

### Maps for Jaccard.80

### 9 Indices + 3 Zoom on Andes

# A/ Species richness
# B/ Species Shannon's Diversity with compatibility with species richness with Jost's transformation (exponential)
# C/ Mean Species rarity

# D/ Faith's Phylogenetic Diversity
# E/ Sum of Fair-proportion
# F/ Mean pairwise Phylogenetic Distance

# G/ Mimicry richness
# H/ Mimicry Shannon's Diversity with compatibility with mimicry richness with Jost's transformation (exponential)
# I/ Mean mimicry rarity

### Zoom on Andes

# J/ Species richness
# K/ Mean pairwise Phylogenetic Distance
# L/ Mean mimicry rarity


# Clean environment
rm(list = ls())

### 1/ Load stuff ####
library(raster)
library(rangeBuilder)
library(gplots)

# New color palette
pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")
pal_bl_light_red_Mannion <- pal_bl_red_Mannion[7:200]

# Choose resolution
res <- "15"
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))

# Load mask for continent borders
continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_", res, ".rds"))
crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".rds"))
bg_mask_pixel <- readRDS(file = "./input_data/Map_stuff/bg_mask_pixel.rds")
bg_mask <- readRDS(file = "./input_data/Map_stuff/bg_mask.rds")

# Load stuff for plot on Andes
Andes_ext <- extent(c(-90, -59, -15, 14))
Close_Andes_ext <- extent(c(-90, -67.2, -15, 14))
DEM <- envData[["Elevation"]]
Andes_DEM <- crop(DEM, Andes_ext)
Close_Andes_DEM <- crop(DEM, Close_Andes_ext)
Andes_bbox <- as(extent(Andes_DEM), 'SpatialPolygons') 


### 2/ Load maps ####
tot.sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
sp.diversity.compatible_Jost_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jost_Jaccard.80.rds"))
sp.mean.rarity_Jaccard.80 <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_Jaccard.80.rds")
PD.raster_Jaccard.80 <- readRDS(file = "./outputs/Indices_Maps/PD.raster_Jaccard.80.rds")
FP.raster_Jaccard.80 <- readRDS(file = "./outputs/Indices_Maps/FP.raster_Jaccard.80.rds")
MPD.raster_Jaccard.80_contrasted <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80_contrasted.rds")
ring.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))
ring.diversity.compatible_Jost_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jost_Jaccard.80.rds"))
mimicry.mean.rarity_Jaccard.80 <- readRDS(file = "./outputs/Indices_maps/mimicry.mean.rarity_Jaccard.80.rds")


### 3/ Generate Andes maps ####

# 3.1/ Species richness

Andes.tot.sp.richness_Jaccard.80 <- crop(tot.sp.richness_Jaccard.80, Andes_ext)

save(Andes.tot.sp.richness_Jaccard.80, file = "./outputs/Indices_maps/Andes/Andes.sp.richness_Jaccard.80.RData", version = "2")
saveRDS(Andes.tot.sp.richness_Jaccard.80, file = "./outputs/Indices_maps/Andes/Andes.sp.richness_Jaccard.80.rds", version = "2")

# 3.2/ MPD

Andes.MPD_Jaccard.80 <- crop(MPD.raster_Jaccard.80_contrasted, Andes_ext)

save(Andes.MPD_Jaccard.80, file = "./outputs/Indices_maps/Andes/Andes.MPD_Jaccard.80.RData", version = "2")
saveRDS(Andes.MPD_Jaccard.80, file = "./outputs/Indices_maps/Andes/Andes.MPD_Jaccard.80.rds", version = "2")

# 3.3/ Mimicry rarity

Andes.mimicry.mean.rarity_Jaccard.80 <- crop(mimicry.mean.rarity_Jaccard.80, Andes_ext)

save(Andes.mimicry.mean.rarity_Jaccard.80, file = "./outputs/Indices_maps/Andes/Andes.mimicry.mean.rarity_Jaccard.80.RData", version = "2")
saveRDS(Andes.mimicry.mean.rarity_Jaccard.80, file = "./outputs/Indices_maps/Andes/Andes.mimicry.mean.rarity_Jaccard.80.rds", version = "2")

### 4/ Plot final figure ####
pdf(file = paste0("./for_article/Index_maps.pdf"), height = 15.5, width = 12)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
par(mfrow = c(4, 3))

# A/ Species richness
image(tot.sp.richness_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Species richness"), 
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.55, padin = c(0.1, 0.1), text.col = "#00000000")
addRasterLegend(tot.sp.richness_Jaccard.80, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.1, label = "Species")
legend(legend = "A", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))

# B/ Species Shannon's diversity
image(sp.diversity.compatible_Jost_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Species Shannon's diversity"), 
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.55, padin = c(0.1, 0.1), text.col = "#00000000")
addRasterLegend(sp.diversity.compatible_Jost_Jaccard.80, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.1, label = "Species")
legend(legend = "B", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))

# C/ Mean species rarity
image(sp.mean.rarity_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mean species rarity"), 
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.55, padin = c(0.1, 0.1), text.col = "#00000000")
addRasterLegend(sp.mean.rarity_Jaccard.80, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.1, label = "Rarity\nindex")
legend(legend = "C", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))

# D/ Faith's PD
image(PD.raster_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Faith's Phylogenetic Diversity"), 
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.55, padin = c(0.1, 0.1), text.col = "#00000000")
addRasterLegend(PD.raster_Jaccard.80, locs = seq(0, 800, 200), minmax = c(0, max(PD.raster_Jaccard.80[], na.rm = T)), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -107, y = 6, font = 2, cex = 1.1, label = "Evolutionary\nTime (My)")
legend(legend = "D", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))

# E/ Sum of Fair-Proportions
image(FP.raster_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Sum of Fair-Proportions"), 
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.55, padin = c(0.1, 0.1), text.col = "#00000000")
addRasterLegend(FP.raster_Jaccard.80, locs = seq(0, 600, 100), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.1, label = "Evolutionary\nTime (My)")
legend(legend = "E", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))

# F/ Mean pairwise Phylogenetic Distance
image(MPD.raster_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Mean pairwise Phylogenetic Distance"), 
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.55, padin = c(0.1, 0.1), text.col = "#00000000")
addRasterLegend(MPD.raster_Jaccard.80_contrasted, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -108, y = 6, font = 2, cex = 1.1, label = "Evolutionary\nTime (My)")
legend(legend = "F", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))

# G/ Mimicry richness
image(ring.richness_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mimicry richness"), 
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mimicry\n         rings", cex=1.2, line = 1, font = 2), 
      legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.55, padin = c(0.1, 0.1), text.col = "#00000000")
addRasterLegend(ring.richness_Jaccard.80, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.1, label = "Mimicry\nrings")
legend(legend = "G", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))

# H/ Mimicry Shannon's Diversity
image(ring.diversity.compatible_Jost_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mimicry Shannon's diversity"), 
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.55, padin = c(0.1, 0.1), text.col = "#00000000")
addRasterLegend(ring.diversity.compatible_Jost_Jaccard.80, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6, font = 2, cex = 1.1, label = "Mimicry\nrings")
legend(legend = "H", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))

# I/ Mean mimicry rarity
image(mimicry.mean.rarity_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Mean mimicry rarity"), 
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.55, padin = c(0.1, 0.1), text.col = "#00000000")
addRasterLegend(mimicry.mean.rarity_Jaccard.80, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5.5, font = 2, cex = 1.1, label = "Rarity\nindex")
legend(legend = "I", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))

# J/ Species richness in Andes

Andes.tot.sp.richness_Jaccard.80@data@min <- 0

image(Andes.tot.sp.richness_Jaccard.80, col = pal_bl_light_red_Mannion, main = paste0("Species richness"), 
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2),  
      legend = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Andes_bbox, lwd = 1, border = "grey20", add = T)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1),
        drawlabels = F, col = "black", add = T)
addRasterLegend(Andes.tot.sp.richness_Jaccard.80, locs = seq(0, 120, 30), cex.axis = 1.1, ramp = pal_bl_light_red_Mannion, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.75, 4.25))
graphics::text(x = -86.7, y = 6, font = 2, cex = 1.0, label = "Species")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 0.9, text.font = 2)
graphics::text(x = -86.7, y = -10.8, font = 2, cex = 1.0, label = "Elevation")
scalebar(d = 500, type = "line", lwd = 4, divs = 4, xy = c(-82, -13.5), label = c("", "500 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.55, padin = c(0.1, 0.1), text.col = "#00000000")
legend(legend = "J", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0.0, 0.0))

# K/ MPD in Andes

image(Andes.MPD_Jaccard.80, col = pal_bl_light_red_Mannion, main = paste0("Mean pairwise Phylogenetic Distance"), 
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2),  
      legend = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Andes_bbox, lwd = 1, border = "grey20", add = T)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1),
        drawlabels = F, col = "black", add = T)
addRasterLegend(Andes.MPD_Jaccard.80, locs = seq(38, 44, 2), cex.axis = 1.1, ramp = pal_bl_light_red_Mannion, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.75, 4.25))
graphics::text(x = -86.2, y = 6, font = 2, cex = 1.0, label = "Evolutionary\nTime (My)")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 0.9, text.font = 2)
graphics::text(x = -86.7, y = -10.8, font = 2, cex = 1.0, label = "Elevation")
scalebar(d = 500, type = "line", lwd = 4, divs = 4, xy = c(-82, -13.5), label = c("", "500 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.55, padin = c(0.1, 0.1), text.col = "#00000000")
legend(legend = "K", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0.0, 0.0))

?addnortharrow

# L/ Mimicry rarity in Andes

Andes.mimicry.mean.rarity_Jaccard.80@data@min <- 0

image(Andes.mimicry.mean.rarity_Jaccard.80, col = pal_bl_red_Mannion[2:200], main = paste0("Mean mimicry rarity"), 
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2),  
      legend = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Andes_bbox, lwd = 1, border = "grey20", add = T)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1),
        drawlabels = F, col = "black", add = T)
addRasterLegend(Andes.mimicry.mean.rarity_Jaccard.80, locs = seq(0, 0.8, 0.2), cex.axis = 1.1, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.75, 4.25))
graphics::text(x = -86.7, y = 6, font = 2, cex = 1.0, label = "Rarity\nindex")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 0.9, text.font = 2)
graphics::text(x = -86.7, y = -10.8, font = 2, cex = 1.0, label = "Elevation")
scalebar(d = 500, type = "line", lwd = 4, divs = 4, xy = c(-82, -13.5), label = c("", "500 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.55, padin = c(0.1, 0.1), text.col = "#00000000")
legend(legend = "L", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0.0, 0.0))

par(mar = internal_margins)
par(mfrow = c(1, 1))

dev.off()




