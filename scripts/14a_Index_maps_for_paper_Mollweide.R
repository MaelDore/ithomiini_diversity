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
library(prettymapr)
library(rangeBuilder)

# New color palette
pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")

# Load mask for continent borders, plot border, and grid
grid_Mollweide_out <- readRDS(file = "./input_data/Map_stuff/grid_Mollweide_out.rds")
large_bg_mask_Mollweide <- readRDS(file = "./input_data/Map_stuff/large_bg_mask_Mollweide.rds")
bbox_sp_Mollweide <- readRDS(file = "./input_data/Map_stuff/bbox_sp_Mollweide.rds")
Andes_grid_Mollweide_out <- readRDS(file = "./input_data/Map_stuff/Andes_grid_Mollweide_out.rds")

# Load stuff for plot on Andes
Andes_ext <- extent(c(-90, -59, -15, 14))
Close_Andes_ext <- extent(c(-90, -67.2, -18, 14))

envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_15.rds"))
DEM <- envData[["Elevation"]]
Andes_DEM <- crop(DEM, Andes_ext)
Close_Andes_DEM <- crop(DEM, Close_Andes_ext)
Andes_bbox <- as(extent(Andes_DEM), 'SpatialPolygons')
Andes_bbox@proj4string@projargs <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" 

saveRDS(Andes_bbox, file = "./input_data/Map_stuff/Andes_bbox.rds", version = "2")

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

# 3.4/ Generate bbox for Andes with Mollweide projection

Andes_bbox_Mollweide <- SpatialPolygons(list(Polygons(list(Polygon(coords = matrix(data = c(-1654, 1879, # top-left
                                                                                            1760, 1879, # top-right
                                                                                            1760, -2014, # bottom-right
                                                                                            -1654, -2014, # bottom-left
                                                                                            -1654, 1879), # top-left
                                                                                   ncol = 2, byrow = T), hole = F))
                                                      , ID = 1))
                                        , proj4string = crs("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))



### 4/ Project all maps in Mollweide projection ####
tot.sp.richness_Jaccard.80_Mollweide <- projectRaster(from = tot.sp.richness_Jaccard.80, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                      crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                      alignOnly = F)
sp.diversity_Jaccard.80_Mollweide <- projectRaster(from = sp.diversity.compatible_Jost_Jaccard.80, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                   crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                   alignOnly = F)
sp.mean.rarity_Jaccard.80_Mollweide <- projectRaster(from = sp.mean.rarity_Jaccard.80, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                     crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                     alignOnly = F)
PD.raster_Jaccard.80_Mollweide <- projectRaster(from = PD.raster_Jaccard.80, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                alignOnly = F)
FP.raster_Jaccard.80_Mollweide <- projectRaster(from = FP.raster_Jaccard.80, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                alignOnly = F)
MPD.raster_Jaccard.80_Mollweide <- projectRaster(from = MPD.raster_Jaccard.80_contrasted, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                 crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                 alignOnly = F)
ring.richness_Jaccard.80_Mollweide <- projectRaster(from = ring.richness_Jaccard.80, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                    crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                    alignOnly = F)
ring.diversity_Jaccard.80_Mollweide <- projectRaster(from = ring.diversity.compatible_Jost_Jaccard.80, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                     crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                     alignOnly = F)
mimicry.mean.rarity_Jaccard.80_Mollweide <- projectRaster(from = mimicry.mean.rarity_Jaccard.80, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                          crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                          alignOnly = F)
mimicry.mean.rarity_Jaccard.80_Mollweide <- projectRaster(from = mimicry.mean.rarity_Jaccard.80, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                          crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                          alignOnly = F)
Andes.tot.sp.richness_Jaccard.80_Mollweide <- projectRaster(from = Andes.tot.sp.richness_Jaccard.80, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                                crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                                alignOnly = F)
Andes.MPD_Jaccard.80_Mollweide <- projectRaster(from = Andes.MPD_Jaccard.80, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                            crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                            alignOnly = F)
Andes.mimicry.mean.rarity_Jaccard.80_Mollweide <- projectRaster(from = Andes.mimicry.mean.rarity_Jaccard.80, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                                crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                                alignOnly = F)
Close_Andes_DEM_Mollweide <- projectRaster(from = Close_Andes_DEM, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                           crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                           alignOnly = F)

### 5/ Plot final figure ####

pdf(file = paste0("./for_article/Index_maps_Mollweide.pdf"), height = 15.5, width = 12)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5, 3.5, 2.1))
par(mar = c(3.1, 3.5, 2.7, 2.1))
par(mfrow = c(4, 3))

# A/ Species richness ####

# Plot without axis
image(tot.sp.richness_Jaccard.80_Mollweide, col = pal_bl_red_Mannion,
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "",
      legend.args = list(text = "          Species", cex = 1.2, line = 1, font = 2), 
      legend = F)
title("Species richness", cex.main = 1.4, line = 1)

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-2600, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.45, padin = c(0.15, 0.15), text.col = "#00000000")
rangeBuilder::addRasterLegend(tot.sp.richness_Jaccard.80, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
rangeBuilder::addRasterLegend(tot.sp.richness_Jaccard.80, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
graphics::text(x = -3550, y = 430, font = 2, cex = 1.1, label = "Species")
legend(legend = "a", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0, -0.008))

# B/ Species Shannon's diversity ####

# Plot without axis
image(sp.diversity_Jaccard.80_Mollweide, col = pal_bl_red_Mannion,
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "",
      legend.args = list(text = "          Species", cex = 1.2, line = 1, font = 2), 
      legend = F)
title("Species Shannon's diversity", cex.main = 1.4, line = 1)

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-2600, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.45, padin = c(0.15, 0.15), text.col = "#00000000")
rangeBuilder::addRasterLegend(sp.diversity_Jaccard.80_Mollweide, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
rangeBuilder::addRasterLegend(sp.diversity_Jaccard.80_Mollweide, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
graphics::text(x = -3550, y = 430, font = 2, cex = 1.1, label = "Species")
legend(legend = "b", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0, -0.008))


# C/ Mean species rarity ####

# Plot without axis
image(sp.mean.rarity_Jaccard.80_Mollweide, col = pal_bl_red_Mannion,
      # main = "Species richness\nMollweide projection", cex.main = 1.4,
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "",
      legend.args = list(text = "          Rarity\n           index", cex = 1.2, line = 1, font = 2), 
      legend = F)
title("Mean species rarity", cex.main = 1.4, line = 1)

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-2600, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.45, padin = c(0.15, 0.15), text.col = "#00000000")
rangeBuilder::addRasterLegend(sp.mean.rarity_Jaccard.80_Mollweide, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
rangeBuilder::addRasterLegend(sp.mean.rarity_Jaccard.80_Mollweide, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
graphics::text(x = -3650, y = 560, font = 2, cex = 1.1, label = "Rarity\nindex")
legend(legend = "c", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0, -0.008))


# D/ Faith's PD ####

# Plot without axis
image(PD.raster_Jaccard.80_Mollweide, col = pal_bl_red_Mannion,
      # main = "Species richness\nMollweide projection", cex.main = 1.4,
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "",
      legend.args = list(text = "          Evolutionary\n         Time (My)", cex = 1.2, line = 1, font = 2), 
      legend = F)
title("Faith's PD", cex.main = 1.4, line = 1)

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-2600, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.45, padin = c(0.15, 0.15), text.col = "#00000000")
rangeBuilder::addRasterLegend(PD.raster_Jaccard.80_Mollweide, locs = seq(0, 800, 200), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
rangeBuilder::addRasterLegend(PD.raster_Jaccard.80_Mollweide, locs = seq(0, 800, 200), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
graphics::text(x = -3350, y = 560, font = 2, cex = 1.1, label = "Evolutionary\nTime (My)")
legend(legend = "d", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0, -0.008))


# E/ Sum of Fair-Proportions ####

# Plot without axis
image(FP.raster_Jaccard.80_Mollweide, col = pal_bl_red_Mannion,
      # main = "Species richness\nMollweide projection", cex.main = 1.4,
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "",
      legend.args = list(text = "          Evolutionary\n         Time (My)", cex = 1.2, line = 1, font = 2), 
      legend = F)
title("Sum of Fair-Proportions", cex.main = 1.4, line = 1)

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-2600, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.45, padin = c(0.15, 0.15), text.col = "#00000000")
rangeBuilder::addRasterLegend(FP.raster_Jaccard.80_Mollweide, locs = seq(0, 500, 100), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
rangeBuilder::addRasterLegend(FP.raster_Jaccard.80_Mollweide, locs = seq(0, 500, 100), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
graphics::text(x = -3350, y = 560, font = 2, cex = 1.1, label = "Evolutionary\nTime (My)")
legend(legend = "e", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0, -0.008))


# F/ Mean pairwise Phylogenetic Distance ####

# Plot without axis
image(MPD.raster_Jaccard.80_Mollweide, col = pal_bl_red_Mannion,
      # main = "Species richness\nMollweide projection", cex.main = 1.4,
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "",
      legend.args = list(text = "          Evolutionary\n         Time (My)", cex = 1.2, line = 1, font = 2), 
      legend = F)
title("Mean pairwise Phylogenetic Distance", cex.main = 1.4, line = 1)

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-2600, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.45, padin = c(0.15, 0.15), text.col = "#00000000")
rangeBuilder::addRasterLegend(MPD.raster_Jaccard.80_Mollweide, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
rangeBuilder::addRasterLegend(MPD.raster_Jaccard.80_Mollweide, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
graphics::text(x = -3350, y = 560, font = 2, cex = 1.1, label = "Evolutionary\nTime (My)")
legend(legend = "f", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0, -0.008))


# G/ Mimicry richness ####

# Plot without axis
image(ring.richness_Jaccard.80_Mollweide, col = pal_bl_red_Mannion,
      # main = "Species richness\nMollweide projection", cex.main = 1.4,
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "",
      legend.args = list(text = "          Mimicry\n         rings", cex = 1.2, line = 1, font = 2), 
      legend = F)
title("Mimicry richness", cex.main = 1.4, line = 1)

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-2600, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.45, padin = c(0.15, 0.15), text.col = "#00000000")
rangeBuilder::addRasterLegend(ring.richness_Jaccard.80_Mollweide, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
rangeBuilder::addRasterLegend(ring.richness_Jaccard.80_Mollweide, locs = seq(0, 25, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
graphics::text(x = -3600, y = 560, font = 2, cex = 1.1, label = "Mimicry\nrings")
legend(legend = "g", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0, -0.008))


# H/ Mimicry Shannon's Diversity ####

# Plot without axis
image(ring.diversity_Jaccard.80_Mollweide, col = pal_bl_red_Mannion,
      # main = "Species richness\nMollweide projection", cex.main = 1.4,
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "",
      legend.args = list(text = "          Mimicry\n         rings", cex = 1.2, line = 1, font = 2), 
      legend = F)
title(" Mimicry Shannon's diversity", cex.main = 1.4, line = 1)

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-2600, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.45, padin = c(0.15, 0.15), text.col = "#00000000")
rangeBuilder::addRasterLegend(ring.diversity_Jaccard.80_Mollweide, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
rangeBuilder::addRasterLegend(ring.diversity_Jaccard.80_Mollweide, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
graphics::text(x = -3600, y = 620, font = 2, cex = 1.1, label = "Mimicry\nrings")
legend(legend = "h", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0, -0.008))


# I/ Mean mimicry rarity ####

# Plot without axis
image(mimicry.mean.rarity_Jaccard.80_Mollweide, col = pal_bl_red_Mannion,
      # main = "Species richness\nMollweide projection", cex.main = 1.4,
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "",
      legend.args = list(text = "          Rarity\n           index", cex = 1.2, line = 1, font = 2), 
      legend = F)
title("Mean mimicry rarity", cex.main = 1.4, line = 1)

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-2600, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.45, padin = c(0.15, 0.15), text.col = "#00000000")
rangeBuilder::addRasterLegend(mimicry.mean.rarity_Jaccard.80_Mollweide, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
rangeBuilder::addRasterLegend(mimicry.mean.rarity_Jaccard.80_Mollweide, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
graphics::text(x = -3650, y = 560, font = 2, cex = 1.1, label = "Rarity\nindex")
legend(legend = "i", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0, -0.008))


# J/ Species richness in Andes ####

# Plot without axis
image(tot.sp.richness_Jaccard.80_Mollweide, col = pal_bl_red_Mannion,
      # main = "Species richness\nMollweide projection", cex.main = 1.4,
      xlim = c(-1654, 1760), ylim = c(-2014, 1879), axes = F,
      xlab = "", ylab = "",
      legend.args = list(text = "          Species", cex = 1.2, line = 1, font = 2), 
      legend = F)
title("Species richness", cex.main = 1.4, line = 1)

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-1550, -420, 520, 1500), labels = c("90°E", "80°E", "70°E", "60°E"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-1850, -1220, -600, 0, 610, 1230), labels = c("15°S", "10°S", "5°S", "0°", "5°N", "10°N"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Andes_grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(Andes_bbox_Mollweide, lwd = 2, border = "black", col = NA, add = T)

contour(x = Close_Andes_DEM_Mollweide, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1),
        drawlabels = F, col = "black", add = T)
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 0.9, text.font = 2)
graphics::text(x = -1300, y = -1450, font = 2, cex = 1.0, label = "Elevation")

scalebar(d = 500, type = "line", lwd = 3.5, divs = 4, xy = c(-800, -1900), label = c("", "500 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.45, padin = c(0.15, 0.15), text.col = "#00000000")
rangeBuilder::addRasterLegend(tot.sp.richness_Jaccard.80_Mollweide, locs = seq(0, 120, 30), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-1460, -1350, -1250, 640))
rangeBuilder::addRasterLegend(tot.sp.richness_Jaccard.80_Mollweide, locs = seq(0, 120, 30), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-1460, -1350, -1250, 640))
graphics::text(x = -1300, y = 840, font = 2, cex = 1.1, label = "Species")

legend(legend = "j", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0, -0.008))


# K/ MPD in Andes ####

# Plot without axis
image(MPD.raster_Jaccard.80_Mollweide, col = pal_bl_red_Mannion,
      # main = "Species richness\nMollweide projection", cex.main = 1.4,
      xlim = c(-1654, 1760), ylim = c(-2014, 1879), axes = F,
      xlab = "", ylab = "",
      legend.args = list(text = "          Evolutionary\n         Time (My)", cex = 1.2, line = 1, font = 2), 
      legend = F)
title("Mean pairwise Phylogenetic Distance", cex.main = 1.4, line = 1)

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-1550, -420, 520, 1500), labels = c("90°E", "80°E", "70°E", "60°E"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-1850, -1220, -600, 0, 610, 1230), labels = c("15°S", "10°S", "5°S", "0°", "5°N", "10°N"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Andes_grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(Andes_bbox_Mollweide, lwd = 2, border = "black", col = NA, add = T)

contour(x = Close_Andes_DEM_Mollweide, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1),
        drawlabels = F, col = "black", add = T)
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 0.9, text.font = 2)
graphics::text(x = -1300, y = -1450, font = 2, cex = 1.0, label = "Elevation")

scalebar(d = 500, type = "line", lwd = 3.5, divs = 4, xy = c(-800, -1900), label = c("", "500 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.45, padin = c(0.15, 0.15), text.col = "#00000000")
rangeBuilder::addRasterLegend(MPD.raster_Jaccard.80_Mollweide, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-1460, -1350, -1250, 640))
rangeBuilder::addRasterLegend(MPD.raster_Jaccard.80_Mollweide, locs = seq(38, 44, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-1460, -1350, -1250, 640))
graphics::text(x = -1200, y = 940, font = 2, cex = 1.1, label = "Evolutionary\nTime (My)")

legend(legend = "k", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0, -0.008))


# L/ Mimicry rarity in Andes ####


# Plot without axis
image(mimicry.mean.rarity_Jaccard.80_Mollweide, col = pal_bl_red_Mannion,
      # main = "Species richness\nMollweide projection", cex.main = 1.4,
      xlim = c(-1654, 1760), ylim = c(-2014, 1879), axes = F,
      xlab = "", ylab = "",
      legend.args = list(text = "          Rarity\n           index", cex = 1.2, line = 1, font = 2), 
      legend = F)
title("Mean mimicry rarity", cex.main = 1.4, line = 1)

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-1550, -420, 520, 1500), labels = c("90°E", "80°E", "70°E", "60°E"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-1850, -1220, -600, 0, 610, 1230), labels = c("15°S", "10°S", "5°S", "0°", "5°N", "10°N"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Andes_grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(Andes_bbox_Mollweide, lwd = 2, border = "black", col = NA, add = T)

contour(x = Close_Andes_DEM_Mollweide, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1),
        drawlabels = F, col = "black", add = T)
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 0.9, text.font = 2)
graphics::text(x = -1300, y = -1450, font = 2, cex = 1.0, label = "Elevation")

scalebar(d = 500, type = "line", lwd = 3.5, divs = 4, xy = c(-800, -1900), label = c("", "500 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.45, padin = c(0.15, 0.15), text.col = "#00000000")
rangeBuilder::addRasterLegend(mimicry.mean.rarity_Jaccard.80_Mollweide, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-1460, -1350, -1250, 640))
rangeBuilder::addRasterLegend(mimicry.mean.rarity_Jaccard.80_Mollweide, locs = seq(0, 0.8, 0.2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-1460, -1350, -1250, 640))
graphics::text(x = -1320, y = 940, font = 2, cex = 1.1, label = "Rarity\nindex")

legend(legend = "l", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0, -0.008))


par(mar = internal_margins)
par(mfrow = c(1, 1))

dev.off()


