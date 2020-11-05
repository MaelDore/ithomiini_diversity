

# Clean environment
rm(list = ls())

### 1/ Generate spatial object needed for the plots ####


### Generate the bbox rectangle
bbox_sp <- SpatialPolygons(list(Polygons(list(Polygon(coords = matrix(data = c(-120, 28, -30, 28, -30, -37, -120, -37, -120, 28),
                                                                      ncol = 2, byrow = T), hole = F))
                                         , ID = 1))
                           , proj4string = crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

save(bbox_sp, file = "./input_data/Map_stuff/bbox_sp.RData", version = "2")

### Generate the country borders

# Simpler world dataset
# library(spData)
# data("world")

# More advance world dataset
library("rnaturalearth")
library("rnaturalearthdata")

# Function to filter what you want from the global world dataset
?ne_countries

world_map <- ne_countries(scale = "medium",
                          continent = NULL,     # To extract only a specific continent
                          country = NULL,       # To extract a specific country
                          returnclass = "sf")   # To get an sp object (with geom in list, an attributes in separated table) or an sf object (with geom nested in a sfc) as output

# Extract only the geometry
world_map <- world_map$geometry

# Extract only the borders within the bbox
country_borders <- st_intersection(world_map, st_as_sf(bbox_sp))

# Save
save(country_borders, file = "./input_data/Map_stuff/country_borders.RData", version = "2")

### Get a background shp to color

# Difference operation work only on sp objects and not sf for no reason...

# Convert worldmap to sp object
world_map_sp <- as(world_map, Class = "Spatial")

# Get the difference
bg_mask <- gDifference(bbox_sp, world_map_sp)
# Get the intersect
continent_mask_smooth <- gIntersection(bbox_sp, world_map_sp)

# Save
save(bg_mask, file = "./input_data/Map_stuff/bg_mask.RData", version = "2")
save(continent_mask_smooth, file = "./input_data/Map_stuff/continent_mask_smooth.RData", version = "2")

### Get a background shp to color, with pixellated borders
library(rgeos)
continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_", res, ".rds"))
crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".rds"))
Neotrop_bbox <- as(extent(continent_mask), 'SpatialPolygons')
bg_mask_pixel <- gDifference(Neotrop_bbox, crop_mask_shp)

save(bg_mask_pixel, file = "./input_data/Map_stuff/bg_mask_pixel.RData", version = "2")
saveRDS(bg_mask_pixel, file = "./input_data/Map_stuff/bg_mask_pixel.rds", version = "2")

### 2/ Load Stuffs ####

load(file = "./input_data/Map_stuff/bg_mask.RData") # Load bg shp
load(file = "./input_data/Map_stuff/country_borders.RData") # Load country borders
rivers <- st_read("./input_data/Map_stuff/Major_rivers/MajorRivers.shp") # Load rivers

# New color palette
pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")

# Outer borders
crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".rds"))

# Load maps to plot
tot.sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))



### 3/ ggplot2 version ####

# It is shit because adding sf object on raster does not work, and NA value are not transparent so you cannot add a nice background color and a grid

library(rasterVis)
library(sf)
library(ggspatial)


tot.sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))

gplot(tot.sp.richness_Jaccard.80) + # Use wrapper function to add a raster to a ggplot2
  geom_tile(aes(fill = value)) + # Plot the z-value of the raster as background
  scale_fill_gradientn(colours = pal_bl_red_Mannion, na.value = rgb(0,0,0,alpha = 0)) + # Provide manually the color scale + use a transparent color for NA vales
  # geom_sf(data = world_map, col = "grey20") +
  # geom_sf(data = st_as_sf(crop_mask_shp), col = "grey20") +
  annotation_scale(location = "tr", width_hint = 0.3, # Add the scale bar
                   pad_x = unit(0.75, "cm"), pad_y = unit(0.5, "cm"),) +   
  annotation_north_arrow(location = "br", which_north = "true", 
                         pad_x = unit(0.75, "cm"), pad_y = unit(0.5, "cm"),  
                         style = north_arrow_fancy_orienteering) +  # Add the north arrow
  coord_sf(crs = tot.sp.richness_Jaccard.80@crs, xlim = c(-120, -30), ylim = c(-37, 28), expand = FALSE) +
  xlab(NULL) + ylab(NULL) +
  theme(panel.grid.major = element_line(color = gray(.5), linetype = 2, size = 0.5),  # To add a grid in the background (below filled polygons)
        panel.background = element_rect(fill = "aliceblue"))  # To add a background color

# Plotting sf objects works alone, but not with the raster called in the ggplot :-(
ggplot() +
  # geom_sf(data = world_map, col = "grey20") +
  geom_sf(data = st_as_sf(crop_mask_shp), col = "grey20") 


### 4/ Base plot version ####

### 4.1/ With image() ####

# See also package prettymapr to add scalebar and north arrow on baseplot/image

library(prettymapr)
library(rangeBuilder)

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/tot.sp.richness_Jaccard.80_3.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(tot.sp.richness_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Species richness \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      # asp = 1, colNA = "aliceblue",
      # xlim = tot.sp.richness_Jaccard.80@extent[1:2], ylim = tot.sp.richness_Jaccard.80@extent[3:4],
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(tot.sp.richness_Jaccard.80, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()



### 4.2/ With plot() ####

# Allow to have a colNA, but do a shitty job with the axis that do not adapt to the plot 
# when the aspect ratio is fixed (default = 1) (leave blank gaps on the sides), 
# or do not allow axis to be rescale to fit better in the pdf when asp = NA 
# (meaning the asp is the one needed to fit the plot, but axis will keep a strict 1:1 ratio)

# Jaccard.80
pdf(file = paste0("./maps/Indices_maps/Jaccard.80/tot.sp.richness_Jaccard.80_2.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
plot(tot.sp.richness_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Species richness \nJaccard.80"), 
     cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
     asp = 0.97, colNA = "aliceblue",
     xlim = tot.sp.richness_Jaccard.80@extent[1:2], ylim = tot.sp.richness_Jaccard.80@extent[3:4],
     ylab = "", xlab = "",
     legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
     legend  = F )
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(world, lwd = 1, border = "grey20", color = "", add = T)
plot(world_map, lwd = 1, border = "grey20", color = "", add = T)
addRasterLegend(tot.sp.richness_Jaccard.80, locs = seq(0, 120, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()

### 5/ Try reprojections ####

# Aitoff. Not pretty
tot.sp.richness_Jaccard.80_Aitoff <- projectRaster(from = tot.sp.richness_Jaccard.80, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                   crs = "+proj=aitoff +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                   alignOnly = F)

image(tot.sp.richness_Jaccard.80, col = pal_bl_red_Mannion)
image(tot.sp.richness_Jaccard.80_Aitoff, col = pal_bl_red_Mannion, main = "Aitoff projection")

# Mollweide. Quite okay
tot.sp.richness_Jaccard.80_Mollweide <- projectRaster(from = tot.sp.richness_Jaccard.80, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                      crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                      alignOnly = F)

image(tot.sp.richness_Jaccard.80, col = pal_bl_red_Mannion)
image(tot.sp.richness_Jaccard.80_Mollweide, col = pal_bl_red_Mannion, main = "Mollweide projection")

# Generate grid lines (graticule) only on oceans areas
?gridlines
large_bbox_sp <- SpatialPolygons(list(Polygons(list(Polygon(coords = matrix(data = c(-140, 40, -10, 40, -10, -50, -140, -50, -140, 40),
                                                                            ncol = 2, byrow = T), hole = F))
                                               , ID = 1))
                                 , proj4string = crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
grid_Mollweide <- gridlines(x = large_bbox_sp, easts = c(-140,-120,-100,-80,-60,-40,-20, 0), norths = c(-40,-30,-20,-10,0,10,20,30))
grid_Mollweide <- spTransform(grid_Mollweide, CRSobj = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
crop_mask_shp_Mollweide <- spTransform(crop_mask_shp, CRSobj = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
grid_Mollweide_out <- rgeos::gDifference(spgeom1 = grid_Mollweide, spgeom2 = crop_mask_shp_Mollweide)

saveRDS(grid_Mollweide_out, file = "./input_data/Map_stuff/grid_Mollweide_out.rds", version = "2")

# Generate grid lines (graticule) only on oceans areas, for Andes region
?gridlines
Andes_grid_Mollweide <- gridlines(x = large_bbox_sp, easts = c(-90,-80,-70,-60), norths = c(-15,-10,-5,0,5,10))
Andes_grid_Mollweide <- spTransform(Andes_grid_Mollweide, CRSobj = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
crop_mask_shp_Mollweide <- spTransform(crop_mask_shp, CRSobj = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
Andes_grid_Mollweide_out <- rgeos::gDifference(spgeom1 = Andes_grid_Mollweide, spgeom2 = crop_mask_shp_Mollweide)

saveRDS(Andes_grid_Mollweide_out, file = "./input_data/Map_stuff/Andes_grid_Mollweide_out.rds", version = "2")


# Generate bg mask for Mollweide projection
library("rnaturalearth")
world_map <- ne_countries(scale = "medium",
                          continent = NULL,     # To extract only a specific continent
                          country = NULL,       # To extract a specific country
                          returnclass = "sf")   # To get an sp object (with geom in list, an attributes in separated table) or an sf object (with geom nested in a sfc) as output
world_map <- world_map$geometry # Extract only the geometry
world_map_sp <- as(world_map, Class = "Spatial") # Convert worldmap to sp object
large_bg_mask <- rgeos::gDifference(large_bbox_sp, world_map_sp) # Get the difference
large_bg_mask_Mollweide <- spTransform(large_bg_mask, CRSobj = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")

saveRDS(large_bg_mask_Mollweide, file = "./input_data/Map_stuff/large_bg_mask_Mollweide.rds", version = "2")


# Generate bbox for Mollweide projection
bbox_sp_Mollweide <- SpatialPolygons(list(Polygons(list(Polygon(coords = matrix(data = c(-4600, 3400, 4600, 3400, 4600, -4450, -4600, -4450, -4600, 3400),
                                                                                ncol = 2, byrow = T), hole = F))
                                                   , ID = 1))
                                     , proj4string = crs("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))

saveRDS(bbox_sp_Mollweide, file = "./input_data/Map_stuff/bbox_sp_Mollweide.rds", version = "2")



### 5.1/ Plot Version Mollweide ####

library(prettymapr)
library(rangeBuilder)

pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")

tot.sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))

grid_Mollweide_out <- readRDS(file = "./input_data/Map_stuff/grid_Mollweide_out.rds")
large_bg_mask_Mollweide <- readRDS(file = "./input_data/Map_stuff/large_bg_mask_Mollweide.rds")
bbox_sp_Mollweide <- readRDS(file = "./input_data/Map_stuff/bbox_sp_Mollweide.rds")

# Reproject
tot.sp.richness_Jaccard.80_Mollweide <- projectRaster(from = tot.sp.richness_Jaccard.80, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                      crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                      alignOnly = F)

# Save final plot
pdf(file = paste0("./maps/Indices_maps/Test_Index_map.pdf"), height = 5.5, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5, 2.7, 2.1))

# Plot without axis
image(tot.sp.richness_Jaccard.80_Mollweide, col = pal_bl_red_Mannion,
      # main = "Species richness\nMollweide projection", cex.main = 1.4,
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

par(mar = internal_margins)

dev.off()



lindr
