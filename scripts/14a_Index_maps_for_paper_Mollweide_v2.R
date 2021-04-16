####################### Index maps for paper ############################

### Maps for Jaccard.80

### Version with 6 Indices + 3 Zoom on Andes

# A/ Species richness
# B/ Mean Species geographic rarity
# C/ Faith's Phylogenetic Diversity

# D/ Mimicry richness
# E/ Mean mimicry geographic rarity
# F/ Weighted mean ring size

### Zoom on Andes

# G/ Species richness
# H/ Mean mimicry geographic rarity
# I/ Weighted mean ring size (with bioregions ?)

# Clean environment
rm(list = ls())

### 1/ Load stuff ####

# Packages
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

# Load bioregions simplifiedshp
Central_America_shp2 <- readRDS(file = "./input_data/Map_stuff/Bioregions/Central_America_shp2.rds")
East_Ecuador_shp2 <- readRDS(file = "./input_data/Map_stuff/Bioregions/East_Ecuador_shp2.rds")
Peruvian_shp2 <- readRDS(file = "./input_data/Map_stuff/Bioregions/Peruvian_shp2.rds")
Mata_Atlantica_shp5 <- readRDS(file = "./input_data/Map_stuff/Bioregions/Mata_Atlantica_shp5.rds")

### Load maps 
sp_richness <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
# sp_mean_geo_rarity <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_Jaccard.80.rds")       # Linear weighting
sp_mean_geo_rarity <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_Leroy_Jaccard.80_contrasted.rds")   # Leroy's weighting
Faith_PD <- readRDS(file = "./outputs/Indices_Maps/PD.raster_Jaccard.80.rds")

ring_richness <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))
# ring_mean_geo_rarity <- readRDS(file = "./outputs/Indices_maps/mimicry.mean.rarity_Jaccard.80.rds")  # Linear weighting
ring_mean_geo_rarity <- readRDS(file = "./outputs/Indices_maps/ring.mean.rarity_Leroy_Jaccard.80_contrasted.rds")    # Leroy's weighting
weighted_mean_ring_size <- readRDS(file = "./outputs/Indices_maps/weighted_mean_ring_size.rds")


### 2/ Generate Andes maps ####

# 2.1/ Function to crop and save maps

crop_to_Andes <- function(x) # Raster to project
{
  x_name <- deparse(substitute(x)) # Get the name of the initial raster as a character string
  
  # Crop into Andes
  new_map <- crop(x, Andes_ext)
  
  # Generate new object with "Andes_" prefix in the global environment
  eval(call("<<-", as.name(paste0("Andes_", x_name)), new_map))
  
  # Save new object in outputs folder
  saveRDS(new_map, file = paste0("./outputs/Indices_maps/Andes/Andes_", x_name, ".rds"))
}

# 2.2/ Crop maps

crop_to_Andes(sp_richness)
crop_to_Andes(ring_mean_geo_rarity)
crop_to_Andes(weighted_mean_ring_size)

# 2.4/ Generate bbox for Andes with Mollweide projection

Andes_bbox_Mollweide <- SpatialPolygons(list(Polygons(list(Polygon(coords = matrix(data = c(-1654, 1879, # top-left
                                                                                            1760, 1879, # top-right
                                                                                            1760, -2014, # bottom-right
                                                                                            -1654, -2014, # bottom-left
                                                                                            -1654, 1879), # top-left
                                                                                   ncol = 2, byrow = T), hole = F))
                                                      , ID = 1))
                                        , proj4string = crs("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))

save(Andes_bbox_Mollweide, file = "./input_data/Map_stuff/Andes_bbox_Mollweide.RData")

### 3/ Project all maps in Mollweide projection ####

# 3.1/ Projection functions
Mollweide_projection <- function(x) # Raster to project
{
  new_map <- projectRaster(from = x, 
                           method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                           crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                           alignOnly = F)
  return(new_map)
}

Mollweide_projection_with_save <- function(x) # Raster to project
{
  x_name <- deparse(substitute(x)) # Get the name of the initial raster as a character string
  
  # Project into Mollweide projection
  new_map <- projectRaster(from = x, 
                           method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                           crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                           alignOnly = F)
  
  # Generate new object with "_Mollweide" suffix in the global environment
  eval(call("<<-", as.name(paste0(x_name, "_Mollweide")), new_map))
  
  # Save new object in outputs folder
  # save(eval(parse(text = paste0(x_name, "_Mollweide"))), file = paste0("./outputs/Indices_maps/mollweide_projections/",x_name, "_Mollweide.RData"))
  saveRDS(new_map, file = paste0("./outputs/Indices_maps/mollweide_projections/",x_name, "_Mollweide.rds"))
}

Mollweide_projection_with_save_and_name <- function(x, name) # Raster to project
{
  # Project into Mollweide projection
  new_map <- projectRaster(from = x, 
                           method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                           crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                           alignOnly = F)
  
  # Generate new object with "_Mollweide" suffix in the global environment
  eval(call("<<-", as.name(paste0(name, "_Mollweide")), new_map))
  
  # Save new object in outputs folder
  # save(eval(parse(text = paste0(name, "_Mollweide"))), file = paste0("./outputs/Indices_maps/mollweide_projections/",name, "_Mollweide.RData"))
  saveRDS(new_map, file = paste0("./outputs/Indices_maps/mollweide_projections/",name, "_Mollweide.rds"))
}

Mollweide_shp_projection <-  function(x) # Shp to project
{
  x_name <- deparse(substitute(x)) # Get the name of the initial shp as a character string
  
  new_shp <- spTransform(x, CRSobj = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
  
  # Generate new object with "_Mollweide" suffix in the global environment
  eval(call("<<-", as.name(paste0(x_name, "_Mollweide")), new_shp))
}

# 3.2/ Project all maps

# Put all maps in a list

list_all_maps <- list(sp_richness, sp_mean_geo_rarity, Faith_PD, ring_richness, ring_mean_geo_rarity, weighted_mean_ring_size, Andes_sp_richness, Andes_ring_mean_geo_rarity, Andes_weighted_mean_ring_size)
names(list_all_maps) <- c("sp_richness", "sp_mean_geo_rarity", "Faith_PD", "ring_richness", "ring_mean_geo_rarity", "weighted_mean_ring_size", "Andes_sp_richness", "Andes_ring_mean_geo_rarity", "Andes_weighted_mean_ring_size")

# Loop to generate projected maps
for (i in 1:length(list_all_maps))
{
  Mollweide_projection_with_save_and_name(x = list_all_maps[[i]], name = names(list_all_maps)[i])
}


# 3.3/ Project DEM
Close_Andes_DEM_Mollweide <- Mollweide_projection(Close_Andes_DEM)
save(Close_Andes_DEM_Mollweide, file = "./input_data/Map_stuff/Close_Andes_DEM_Mollweide.RData")

# 3.4/ Project bioregions shp
Mollweide_shp_projection(Central_America_shp2)
Mollweide_shp_projection(East_Ecuador_shp2)
Mollweide_shp_projection(Peruvian_shp2)
Mollweide_shp_projection(Mata_Atlantica_shp5)

### 4/ Plot final figure ####

# 4.1/ Plotting function ####

{
  map_indices_Mollweide <- function(x,                                    # Raster to map
                                  color_palette = pal_bl_red_Mannion,   # Color palette
                                  main_title,                           # Main title
                                  main_title_cex = 1.5,                 # Main title size
                                  
                                  xlim = c(-4600, 4600),   # Limit of plot on x-axis (Longitude)
                                  ylim = c(-4450, 3400),    # Limit of plot on y-axis (Latitude)
                                  axis_cex = 1.4,             # Axes size
                                  
                                  xlab = "",                # X-axis label
                                  ylab = "",                # Y-axis label
                                  x_axis_breaks = c(-3930, -2170, -420, 1500, 3050),            # X-axis tick breaks
                                  y_axis_breaks = c(-3650, -2450, -1220, 0, 1230, 2445),        # Y-axis tick breaks
                                  x_axis_labels = c("120°E", "100°E", "80°E", "60°E", "40°E"),      # X-axis tick labels
                                  y_axis_labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"),  # Y-axis tick labels
                                  
                                  legend_title,             # Legend title
                                  legend_title_cex = 1.4,   # Legend title size
                                  legend_title_x = -3550,   # Legend title x position
                                  legend_title_y = 430,     # Legend title y position
                                  legend_cex = 1.4,         # Legend size
                                  legend_breaks,            # Legend tick positions
                                  legend_location = c(-4100, -3800, -3950, 0),  # Legend position
                                  
                                  scale_bar_position = c(-2600, -4000),  # Scale bar position
                                  
                                  arrow_scale = 0.45,           # North arrow size
                                  arrow_padin = c(0.15, 0.15),  # North arrow position adjustement
                                  
                                  facet_letter = "",                  # Small case letter for facet
                                  facet_letter_col = "black",         # Color of case letter for facet
                                  facet_letter_cex = 2.2,             # Size of small case letter for facet
                                  facet_letter_inset = c(0, -0.008))  # Position adjustment of small case letter for facet

{
  # Plot raster background without axis
  image(x, col = color_palette,
        xlim = xlim, ylim = ylim, axes = F,
        xlab = xlab, ylab = ylab)
  title(main = main_title, cex.main = main_title_cex, line = 1)
  
  # Generate axes with manual positioning of ticks
  axis(side = 1, at = x_axis_breaks, labels = x_axis_labels, cex.axis = axis_cex, lwd = 0.2, lwd.ticks = 1, gap.axis = 0, padj = 0.5)
  axis(side = 2, at = y_axis_breaks, labels = y_axis_labels, cex.axis = axis_cex, lwd = 0.2, lwd.ticks = 1, gap.axis = 0)
  
  # Add background, borders and graticules
  plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
  plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
  plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
  
  # Add scale bar in legend
  scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = scale_bar_position, label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1.2)
  prettymapr::addnortharrow(scale = arrow_scale, padin = arrow_padin, text.col = "#00000000")
  rangeBuilder::addRasterLegend(x, locs = legend_breaks, cex.axis = legend_cex, ramp = color_palette, ncolors = 200, border = T, location = legend_location)
  rangeBuilder::addRasterLegend(x, locs = legend_breaks, cex.axis = legend_cex, ramp = color_palette, ncolors = 200, border = T, location = legend_location)
  graphics::text(x = legend_title_x, y = legend_title_y, font = 2, cex = legend_title_cex, label = legend_title)
  
  # Add facet letter
  legend(legend = facet_letter, x = "bottomright", bty = "n", text.col = facet_letter_col,
         text.font = 2, cex = facet_letter_cex, inset = facet_letter_inset)
  
  }
}

{
  map_indices_Andes_Mollweide <- function(x,                                    # Raster to map
                                          color_palette = pal_bl_red_Mannion,   # Color palette
                                          main_title,                           # Main title
                                          main_title_cex = 1.5,                 # Main title size
                                          
                                          xlim = c(-1654, 1760),    # Limit of plot on x-axis (Longitude)
                                          ylim = c(-2014, 1879),    # Limit of plot on y-axis (Latitude)
                                          axis_cex = 1.4,             # Axes size
                                          
                                          xlab = "",                # X-axis label
                                          ylab = "",                # Y-axis label
                                          x_axis_breaks = c(-1550, -420, 520, 1500),
                                          x_axis_labels = c("90°E", "80°E", "70°E", "60°E"),
                                          y_axis_breaks = c(-1850, -1220, -600, 0, 610, 1230), 
                                          y_axis_labels = c("15°S", "10°S", "5°S", "0°", "5°N", "10°N"),
                                          
                                          legend_title,             # Legend title
                                          legend_title_cex = 1.4,   # Legend title size
                                          legend_title_x = -1300,   # Legend title x position
                                          legend_title_y = 840,     # Legend title y position
                                          legend_cex = 1.4,         # Legend size
                                          legend_breaks,            # Legend tick positions
                                          legend_location = c(-1460, -1350, -1250, 640),  # Legend position
                                          
                                          elevation = TRUE, # Should elevation contour and assoicated legend added?
                                          
                                          scale_bar_position = c(-800, -1900),  # Scale bar position
                                          
                                          arrow_scale = 0.45,           # North arrow size
                                          arrow_padin = c(0.15, 0.15),  # North arrow position adjustement
                                          
                                          facet_letter = "",                  # Small case letter for facet
                                          facet_letter_col = "black",         # Color of case letter for facet
                                          facet_letter_cex = 2.2,             # Size of small case letter for facet
                                          facet_letter_inset = c(0, -0.008))  # Position adjustment of small case letter for facet
  
  {
    # Plot raster background without axis
    image(x, col = color_palette,
          xlim = xlim, ylim = ylim, axes = F,
          xlab = xlab, ylab = ylab)
    title(main = main_title, cex.main = main_title_cex, line = 1)
    
    # Generate axes with manual positioning of ticks
    axis(side = 1, at = x_axis_breaks, labels = x_axis_labels, cex.axis = axis_cex, lwd = 0.2, lwd.ticks = 1, gap.axis = 0, padj = 0.5)
    axis(side = 2, at = y_axis_breaks, labels = y_axis_labels, cex.axis = axis_cex, lwd = 0.2, lwd.ticks = 1, gap.axis = 0)
    
    # Add background, borders and graticules
    plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
    plot(Andes_grid_Mollweide_out, lty = 92, col = "grey80", add = T)
    plot(Andes_bbox_Mollweide, lwd = 2, border = "black", col = NA, add = T)
    
    if (elevation == TRUE) 
    {
      # Add elevation contour 
      contour(x = Close_Andes_DEM_Mollweide, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1),
              drawlabels = F, col = "black", add = T)
      legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 1.0, text.font = 2)
      graphics::text(x = -1300, y = -1500, font = 2, cex = 1.1, label = "Elevation")
    }
    
    # Add scale bar in legend
    scalebar(d = 500, type = "line", lwd = 3.5, divs = 4, xy = scale_bar_position, label = c("", "500 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1.2)
    prettymapr::addnortharrow(scale = arrow_scale, padin = arrow_padin, text.col = "#00000000")
    rangeBuilder::addRasterLegend(x, locs = legend_breaks, cex.axis = legend_cex, ramp = color_palette, ncolors = 200, border = T, location = legend_location)
    rangeBuilder::addRasterLegend(x, locs = legend_breaks, cex.axis = legend_cex, ramp = color_palette, ncolors = 200, border = T, location = legend_location)
    graphics::text(x = legend_title_x, y = legend_title_y, font = 2, cex = legend_title_cex, label = legend_title)
    
    # Add facet letter
    legend(legend = facet_letter, x = "bottomright", bty = "n", text.col = facet_letter_col,
           text.font = 2, cex = facet_letter_cex, inset = facet_letter_inset)
    
  }
}


{
  map_indices_Andes_Mollweide_down <- function(x,                                    # Raster to map
                                    color_palette = pal_bl_red_Mannion,   # Color palette
                                    main_title,                           # Main title
                                    main_title_cex = 1.5,                 # Main title size
                                    
                                    xlim = c(-1654, 1760),    # Limit of plot on x-axis (Longitude)
                                    ylim = c(-2014, 1879),    # Limit of plot on y-axis (Latitude)
                                    axis_cex = 1.4,             # Axes size
                                    
                                    xlab = "",                # X-axis label
                                    ylab = "",                # Y-axis label
                                    x_axis_breaks = c(-1550, -420, 520, 1500),
                                    x_axis_labels = c("90°E", "80°E", "70°E", "60°E"),
                                    y_axis_breaks = c(-1850, -1220, -600, 0, 610, 1230), 
                                    y_axis_labels = c("15°S", "10°S", "5°S", "0°", "5°N", "10°N"),
                                    
                                    legend_title,             # Legend title
                                    legend_title_cex = 1.4,   # Legend title size
                                    legend_title_x = -1300,   # Legend title x position
                                    legend_title_y = 290,     # Legend title y position
                                    legend_cex = 1.4,         # Legend size
                                    legend_breaks,            # Legend tick positions
                                    legend_location = c(-1460, -1350, -1800, 90),  # Legend position
                                    
                                    elevation = TRUE,            # Should elevation contour and associated legend added?
                                    y_adjust_lvl = 0,            # Adjustment for elevation level legend
                                    y_adjust_elevation = 0,      # Adjustment for elevation label
                                    elevation_title_cex = 1.3,   # Size of Elevation title
                                    elevation_legend_cex = 1.2,  # Size of Elevation legend
                                    
                                    scale_bar_position = c(-1000, -1800),  # Scale bar position
                                    
                                    arrow_scale = 0.45,           # North arrow size
                                    arrow_padin = c(0.15, 0.15),  # North arrow position adjustment
                                    
                                    facet_letter = "",                  # Small case letter for facet
                                    facet_letter_col = "black",         # Color of case letter for facet
                                    facet_letter_cex = 2.2,             # Size of small case letter for facet
                                    facet_letter_inset = c(0, -0.008))  # Position adjustment of small case letter for facet
  
  {
    # Plot raster background without axis
    image(x, col = color_palette,
          xlim = xlim, ylim = ylim, axes = F,
          xlab = xlab, ylab = ylab)
    title(main = main_title, cex.main = main_title_cex, line = 1)
    
    # Generate axes with manual positioning of ticks
    axis(side = 1, at = x_axis_breaks, labels = x_axis_labels, cex.axis = axis_cex, lwd = 0.2, lwd.ticks = 1, gap.axis = 0, padj = 0.5)
    axis(side = 2, at = y_axis_breaks, labels = y_axis_labels, cex.axis = axis_cex, lwd = 0.2, lwd.ticks = 1, gap.axis = 0)
    
    # Add background, borders and graticules
    plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
    plot(Andes_grid_Mollweide_out, lty = 92, col = "grey80", add = T)
    plot(Andes_bbox_Mollweide, lwd = 2, border = "black", col = NA, add = T)
    
    if (elevation == TRUE) 
    {
      # Add elevation contour 
      contour(x = Close_Andes_DEM_Mollweide, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1),
              drawlabels = F, col = "black", add = T)
      # Add legend for elevation
      legend(x = "bottomleft", inset = c(0, 0.63 + y_adjust_lvl), title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = elevation_legend_cex, text.font = 2)
      graphics::text(x = -1250, y = 1040 + y_adjust_elevation, font = 2, cex = elevation_title_cex, label = "Elevation")
    }

    # Add scale bar in legend
    scalebar(d = 500, type = "line", lwd = 3.5, divs = 4, xy = scale_bar_position, label = c("", "500 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1.2)
    prettymapr::addnortharrow(scale = arrow_scale, padin = arrow_padin, text.col = "#00000000")
    rangeBuilder::addRasterLegend(x, locs = legend_breaks, cex.axis = legend_cex, ramp = color_palette, ncolors = 200, border = T, location = legend_location)
    rangeBuilder::addRasterLegend(x, locs = legend_breaks, cex.axis = legend_cex, ramp = color_palette, ncolors = 200, border = T, location = legend_location)
    graphics::text(x = legend_title_x, y = legend_title_y, font = 2, cex = legend_title_cex, label = legend_title)
    
    # Add facet letter
    legend(legend = facet_letter, x = "bottomright", bty = "n", text.col = facet_letter_col,
           text.font = 2, cex = facet_letter_cex, inset = facet_letter_inset)
    
  }
}

?legend

pdf(file = paste0("./for_article/Index_maps_Mollweide_v2.pdf"), height = 12, width = 12)

internal_margins <- par()$mar
# par(mar = c(3.1, 3.5, 3.5, 2.1))
par(mar = c(3.1, 3.1, 2.7, 1.6))
par(mfrow = c(3, 3))

# A/ Species richness ####

map_indices_Mollweide(x = sp_richness_Mollweide,
                      main_title = "Species richness",
                      legend_title = "Species",
                      legend_title_x = -3650,
                      legend_title_y = 430,
                      legend_breaks = seq(0, 120, 20), 
                      facet_letter = "(a)")

# B/ Species mean geographic rarity ####

# Need to add a few values to lower values to avoid being colored as null outer-range values
sp_mean_geo_rarity_Mollweide[(sp_mean_geo_rarity_Mollweide[] > 0) & (sp_mean_geo_rarity_Mollweide[] <= 0.005)] <- 0.005

map_indices_Mollweide(x = sp_mean_geo_rarity_Mollweide,
                      main_title = "Species geographic rarity",
                      legend_title = "Rarity\nindex",
                      legend_title_x = -3650,
                      legend_title_y = 870,
                      legend_breaks = seq(0, 0.5, 0.1), 
                      facet_letter = "(b)")

# C/ Faith's PD ####

map_indices_Mollweide(x = Faith_PD_Mollweide,
                      main_title = "Phylogenetic diversity",
                      legend_title = "Evolutionary\nTime (My)",
                      legend_title_x = -3150,
                      legend_title_y = 670,
                      legend_breaks = seq(0, 800, 200), 
                      facet_letter = "(c)")

# D/ Mimicry richness ####

map_indices_Mollweide(x = ring_richness_Mollweide,
                      main_title = "Mimicry richness",
                      legend_title = "Mimicry\nrings",
                      legend_title_x = -3670,
                      legend_title_y = 670,
                      legend_breaks = seq(0, 25, 5), 
                      facet_letter = "(d)")

# E/ Mimicry ring mean geographic rarity ####

map_indices_Mollweide(x = ring_mean_geo_rarity_Mollweide,
                      main_title = "Mimicry geographic rarity",
                      legend_title = "Rarity\nindex",
                      legend_title_x = -3650,
                      legend_title_y = 670,
                      legend_breaks = seq(0, 0.6, 0.2), 
                      facet_letter = "(e)")

# F/ Mean ring size ####

map_indices_Mollweide(x = weighted_mean_ring_size_Mollweide,
                      main_title = "Mean ring size",
                      legend_title = "Ring\nsize",
                      legend_title_x = -3750,
                      legend_title_y = 670,
                      legend_breaks = seq(0, 6, 1), 
                      facet_letter = "(f)")

# G/ Species richness in Andes ####

map_indices_Andes_Mollweide_down(x = sp_richness_Mollweide,
                            main_title = "Species richness",
                            legend_title = "Species",
                            y_adjust_lvl = 0.0,            # Adjustment for elevation level legend
                            y_adjust_elevation = 0, 
                            legend_breaks = seq(0, 120, 20), 
                            facet_letter = "(g)",
                            facet_letter_col = "white")

# H/ Mean species geographic rarity in Andes ####

map_indices_Andes_Mollweide_down(x = sp_mean_geo_rarity_Mollweide,
                                 color_palette = pal_bl_red_Mannion[-1],
                                 main_title = "Species geographic rarity",
                                 legend_title = "Rarity index",
                                 legend_title_x = -1100,
                                 legend_title_y = 350,
                                 y_adjust_lvl = 0.0,
                                 y_adjust_elevation = 0,
                                 legend_breaks = seq(0, 0.5, 0.1), 
                                 facet_letter = "(h)",
                                 facet_letter_col = "white")

# I/ Mean ring size in Andes ####

map_indices_Andes_Mollweide_down(x = weighted_mean_ring_size_Mollweide,
                            main_title = "Mean ring size",
                            legend_title = "Ring\nsize",
                            legend_title_x = -1330,
                            legend_title_y = 400,
                            elevation = FALSE,
                            legend_breaks = seq(0, 6, 1), 
                            facet_letter = "(i)",
                            facet_letter_col = "white")

plot(Central_America_shp2_Mollweide, col = NA, border = "#F8766D", lwd = 4, add = T)
plot(East_Ecuador_shp2_Mollweide, col = NA, border = "#7CAE00", lwd = 4, add = T)
plot(Peruvian_shp2_Mollweide, col = NA, border = "#00BFC4", lwd = 4, add = T)

graphics::text(x = -570, y = 1500, font = 2, cex = 1.4, label = "CA", col = "#F8766D")
graphics::text(x = -780, y = 10, font = 2, cex = 1.4, label = "EE", col = "#7CAE00")
graphics::text(x = -640, y = -1350, font = 2, cex = 1.4, label = "EP", col = "#00BFC4")

segments(x0 = -610, y0 = 1350 , x1 = -700, y1 = 1100, col = "#F8766D", lwd = 4)
segments(x0 = -620, y0 = -20 , x1 = -400, y1 = -100, col = "#7CAE00", lwd = 4)
segments(x0 = -490, y0 = -1320 , x1 = -250, y1 = -1150, col = "#00BFC4", lwd = 4)

par(mar = internal_margins)
par(mfrow = c(1, 1))

dev.off()


