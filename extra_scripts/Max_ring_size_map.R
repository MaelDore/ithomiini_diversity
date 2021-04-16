#### Generate map of max ring size ####


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

### Load bioregions simplified shp ###
Central_America_shp2 <- readRDS(file = "./input_data/Map_stuff/Bioregions/Central_America_shp2.rds")
East_Ecuador_shp2 <- readRDS(file = "./input_data/Map_stuff/Bioregions/East_Ecuador_shp2.rds")
Peruvian_shp2 <- readRDS(file = "./input_data/Map_stuff/Bioregions/Peruvian_shp2.rds")
Mata_Atlantica_shp2 <- readRDS(file = "./input_data/Map_stuff/Bioregions/Mata_Atlantica_shp2.rds")

# Load mimicry ring richness stack
All_ring_rich_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_Jaccard.80.RData"))

# Load Mean mimicry ring size map for later plot
load(file = "./outputs/Correlation_tests/mean.ring.size.RData")


#### 2/ Compute index as max richness among all mimicry rings for each pixel #####

ring.max.size_Jaccard.80 <- calc(All_ring_rich_stack_Jaccard.80, fun = max)*1

#### 3/ Project into Mollweide's ####

ring.max.size_Jaccard.80_Mollweide <- projectRaster(from = ring.max.size_Jaccard.80, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                    crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                    alignOnly = F)

mean.ring.size_Jaccard.80_Mollweide <- projectRaster(from = mean.ring.size, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                     crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                     alignOnly = F)

Central_America_shp2_Mollweide <- spTransform(x = Central_America_shp2, CRSobj = crs("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
East_Ecuador_shp2_Mollweide <- spTransform(x = East_Ecuador_shp2, CRSobj = crs("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
Peruvian_shp2_Mollweide <- spTransform(x = Peruvian_shp2, CRSobj = crs("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
Mata_Atlantica_shp2_Mollweide <- spTransform(x = Mata_Atlantica_shp2, CRSobj = crs("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))


#### 4/ Plot alone ####

pdf(file = paste0("./maps/Indices_maps/ring.max.size.pdf"), height = 5.3, width = 6.5)

# Plot without axis
image(ring.max.size_Jaccard.80_Mollweide, col = pal_bl_red_Mannion,
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "",
      legend.args = list(text = "          Ring max size", cex = 1.2, line = 1, font = 2), 
      legend = F)
title("Mimicry ring max size", cex.main = 1.4, line = 1)

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-2600, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.45, padin = c(0.15, 0.15), text.col = "#00000000")
rangeBuilder::addRasterLegend(ring.max.size_Jaccard.80_Mollweide, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
rangeBuilder::addRasterLegend(ring.max.size_Jaccard.80_Mollweide, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
graphics::text(x = -3550, y = 730, font = 2, cex = 1.1, label = "Mimicry ring\n max size")
legend(legend = "b", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0, -0.008))

dev.off()


#### 5/ Plot with Mean ring size map ####



pdf(file = paste0("./maps/Indices_maps/mean_max_ring.size.pdf"), height = 4.8, width = 11)

internal_margins <- par()$mar
par(mar = c(2.6,2.5,2.6,1.0), mfrow = c(1,2))

# Set margin parameters: Bottom, left, top, right

# A/ Mean ring size

mean.ring.size_Jaccard.80_Mollweide@data@min <- 0

# Plot without axis
image(mean.ring.size_Jaccard.80_Mollweide, col = pal_bl_red_Mannion,
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "",
      # legend.args = list(text = "          Mean ring size", cex = 1.2, line = 1, font = 2), 
      legend = F)
title("Mimicry ring mean size", cex.main = 1.4, line = 1)

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)

plot(Central_America_shp2_Mollweide, col = NA, border = "red", lwd = 1.5, add = T)
plot(East_Ecuador_shp2_Mollweide, col = NA, border = "darkgreen", lwd = 1.5, add = T)
plot(Peruvian_shp2_Mollweide, col = NA, border = "blue", lwd = 1.5, add = T)
plot(Mata_Atlantica_shp2_Mollweide, col = NA, border = "purple", lwd = 1.5, add = T)

legend(x = "bottomleft", 
       # legend = c("Central America", "East Ecuador", "East Peru", "Brazilian Atlantic Forest"),
       legend = c("CA", "EE", "EP", "BAF"),
       border = c("red", "darkgreen", "blue", "purple") , fill = NA,
       inset = c(0.18,0.15),
       cex = 1.2,
       bty = "n") 

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-2600, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.45, padin = c(0.15, 0.15), text.col = "#00000000")
rangeBuilder::addRasterLegend(mean.ring.size_Jaccard.80_Mollweide, locs = seq(0, 6, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
rangeBuilder::addRasterLegend(mean.ring.size_Jaccard.80_Mollweide, locs = seq(0, 6, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
graphics::text(x = -3550, y = 750, font = 2, cex = 0.9, label = "Mimicry ring\n mean size")
legend(legend = "a", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0, -0.008))

# B/ Max ring size

# Plot without axis
image(ring.max.size_Jaccard.80_Mollweide, col = pal_bl_red_Mannion,
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "",
      # legend.args = list(text = "          Max ring size", cex = 1.2, line = 1, font = 2), 
      legend = F)
title("Mimicry ring max size", cex.main = 1.4, line = 1)

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)

plot(Central_America_shp2_Mollweide, col = NA, border = "red", lwd = 1.5, add = T)
plot(East_Ecuador_shp2_Mollweide, col = NA, border = "darkgreen", lwd = 1.5, add = T)
plot(Peruvian_shp2_Mollweide, col = NA, border = "blue", lwd = 1.5, add = T)
plot(Mata_Atlantica_shp2_Mollweide, col = NA, border = "purple", lwd = 1.5, add = T)

legend(x = "bottomleft", 
       # legend = c("Central America", "East Ecuador", "East Peru", "Brazilian Atlantic Forest"),
       legend = c("CA", "EE", "EP", "BAF"),
       border = c("red", "darkgreen", "blue", "purple") , fill = NA,
       inset = c(0.18,0.15),
       cex = 1.2,
       bty = "n")  

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-2600, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.45, padin = c(0.15, 0.15), text.col = "#00000000")
rangeBuilder::addRasterLegend(ring.max.size_Jaccard.80_Mollweide, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
rangeBuilder::addRasterLegend(ring.max.size_Jaccard.80_Mollweide, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))
graphics::text(x = -3550, y = 750, font = 2, cex = 0.9, label = "Mimicry ring\n max size")
legend(legend = "b", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0, -0.008))


par(mar = internal_margins, mfrow = c(1,1))
dev.off()

