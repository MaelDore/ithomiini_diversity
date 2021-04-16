
##### Script 22: Map occurrences #####

# A pretty map of all occurrences (Figure S1.1)

# Input: 
   # Occurrences data set

# Outputs:
   # Map of all occurrences in the study dataset (Figure S1.1)



# Remove environment
rm(list = ls())

### 1/ Load stuff ####

# Packages
library(raster)
library(prettymapr)
library(rangeBuilder)
library(sf)

# Load occurrence dataset
load(file = "./input_data/Ithomiini_final.RData")
Ithomiini_records <- readRDS(file = "./input_data/Ithomiini_records.rds")

# Load terrestrial mask as grid reference
continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_15.rds"))
plot(continent_mask, col = "red")

continent_mask_shp <- rasterToPolygons(x = continent_mask)

# Load mask for continent borders, plot border, and grid
grid_Mollweide_out <- readRDS(file = "./input_data/Map_stuff/grid_Mollweide_out.rds")
large_bg_mask_Mollweide <- readRDS(file = "./input_data/Map_stuff/large_bg_mask_Mollweide.rds")
bbox_sp_Mollweide <- readRDS(file = "./input_data/Map_stuff/bbox_sp_Mollweide.rds")

load(file = "./input_data/Map_stuff/country_borders.RData")
country_borders <- as(country_borders, "Spatial")
# plot(country_borders)

rivers <- st_read("./input_data/Map_stuff/Major_rivers/MajorRivers.shp") # Load rivers
rivers <- as(rivers, "Spatial")
xmin = -120 ; xmax = -30 ; ymin = -37 ; ymax = 28
e <- extent(c(xmin,xmax,ymin,ymax))
rivers <- crop(rivers, e)

### 2/ Generate spatial point object of occurrences

Occ_shp <- SpatialPoints(coords = Ithomiini_final[, c("Longitude", "Latitude")],
                         proj4string = continent_mask@crs)
# Occ_shp <- SpatialPoints(coords = Ithomiini_records[, c("Longitude", "Latitude")],
#                          proj4string = continent_mask@crs)

# Filter out occurrences falling in the sea
Occ_shp <- rgeos::gIntersection(spgeom1 = country_borders, spgeom2 = Occ_shp)

# plot(Occ_shp, pch = 16)

### 3/ Project all stuff in Mollweide projection ####

# Projection function
Mollweide_shp_projection <-  function(x) # Shp to project
{
  x_name <- deparse(substitute(x)) # Get the name of the initial shp as a character string
  
  new_shp <- spTransform(x, CRSobj = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
  
  # Generate new object with "_Mollweide" suffix in the global environment
  eval(call("<<-", as.name(paste0(x_name, "_Mollweide")), new_shp))
}

Mollweide_projection <- function(x) # Raster to project
{
  x_name <- deparse(substitute(x)) # Get the name of the initial raster as a character string
  
  # Project into Mollweide projection
  new_map <- projectRaster(from = x, 
                           method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                           crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                           alignOnly = F)
  
  # Generate new object with "_Mollweide" suffix in the global environment
  eval(call("<<-", as.name(paste0(x_name, "_Mollweide")), new_map))
}

# Apply projection
Mollweide_shp_projection(country_borders)
Mollweide_shp_projection(rivers)
Mollweide_shp_projection(Occ_shp)

Mollweide_projection(continent_mask)


### 4/ Map

pdf(file = paste0("./supplementaries/Occurrence_map.pdf"), height = 6, width = 6)

internal_margins <- par()$mar
par(mar = c(3.1, 3.1, 1.5, 1.6))

# Plot raster background without axis
image(x = continent_mask_Mollweide, col = "antiquewhite",
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.4, lwd = 0.2, lwd.ticks = 1, gap.axis = 0, padj = 0.5)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.4, lwd = 0.2, lwd.ticks = 1, gap.axis = 0)

# Add background, borders and graticules
plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
plot(country_borders_Mollweide, lwd = 1, border = "#FF000030", col = NA, add = T)
plot(rivers_Mollweide, lwd = 1, col = "#6495ED60", add = T)

# Add occurrence points
plot(Occ_shp_Mollweide, pch = 16, cex = 0.5, col = "#00000090", add = T)

# Add scale bar in legend
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4000, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1.1)
prettymapr::addnortharrow(scale = 0.45, padin = c(0.15, 0.15), text.col = "#00000000")
legend(legend = "Country borders", x = -4450, y = -2350, cex = 1.1, lty = 1, lwd = 2, col = "red", bty = "n")
legend(legend = "Country borders", x = -4450, y = -2350, cex = 1.1, lty = 1, lwd = 2, col = "red", bty = "n")
legend(legend = "Rivers", x = -4450, y = -2850, cex = 1.1, lty = 1, lwd = 2, col = "#6495ED", bty = "n")
legend(legend = "Rivers", x = -4450, y = -2850, cex = 1.1, lty = 1, lwd = 2, col = "#6495ED", bty = "n")

par(mar = internal_margins)

dev.off()
