
##### Script 16: Build maps of forest cover changes to assess threats to Ithomiini diversity #####

# Build maps of absolute and relative forest cover changes under 2 SSP-RCPs scenarios: SSP2-RCP4.5  and SSP5-RCP8.5

# Data from Land Use Harmonization 2 (CMIP6 framework)

#####

# Inputs:
   # Land cover changes from Land Use Harmonization 2 (CMIP6 framework)

# Outputs:
   # Maps of relative and absolute forest cover changes under 2 SSP-RCPs scenarios: SSP2-RCP4.5  and SSP5-RCP8.5
   # Global maps and regional maps (Andes and Mata Atlantica)

#####


# Remove environment
rm(list = ls())

# Load libraries
library(raster)
library(ncdf4) 

# Choose resolution
res <- "15"
# Load environmental stack to use as template for resolution, extent, and terrestrial boundaries
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))



### 1/ Extract the baseline as the mean values from 1970-2000 such as for the climate baseline ####

filepath <- "./input_data/LUH2/states.nc"

# Get raster Brick from NetCDF file for only one variable
nc_raster_primf <- brick(filepath, var = "primf") # Primary forests
nc_raster_secdf <- brick(filepath, var = "secdf") # Secondary forests
nlayers(nc_raster_primf) # As many layers as time step from 850 to 2015

# Extract the time step from 1970 to 2000. File start in 850
nc_raster_primf_cover_70_2000 <- subset(nc_raster_primf, 1121:1151)
nc_raster_secdf_cover_70_2000 <- subset(nc_raster_secdf, 1121:1151)
names(nc_raster_primf_cover_70_2000) <- names(nc_raster_secdf_cover_70_2000) <- as.character(1970:2000)

# Remove rasterBrick per variable to save RAM space
rm(nc_raster_primf, nc_raster_secdf)

# Sum the two forest covers
forest_cover_70_2000 <- nc_raster_primf_cover_70_2000 + nc_raster_secdf_cover_70_2000

# plot(forest_cover_70_2000[[1]])

# Crop and resample to Neotropics
forest_cover_70_2000_cropped <- resample(x = forest_cover_70_2000, y = envData, method = "bilinear") # Ajust extent and resolution to match the one of our SDM maps

# Compute mean
forest_cover_baseline_cropped <- calc(forest_cover_70_2000_cropped, fun = mean, na.rm = T)
# Compute interannual variation over 20 years
forest_cover_baseline_cropped_sd <- calc(forest_cover_70_2000_cropped, fun = sd, na.rm = T)

plot(forest_cover_baseline_cropped)
plot(forest_cover_baseline_cropped_sd)

# Save rasters
saveRDS(forest_cover_baseline_cropped, file = paste0("./input_data/LUH2/Forest_cover_baseline.rds"), version = "2")
saveRDS(forest_cover_baseline_cropped_sd, file = paste0("./input_data/LUH2/Forest_cover_baseline_sd.rds"), version = "2")



### 2/ Extract the projections for 2061 to 2081 (such as the climate models for 2070) ####

# Set Time period and target year
period <- "2061-2080" ; target <- "2070"

# Per SSP-RCP
SSP_list <- c("ssp245", "ssp585")
# SSP <- SSP_list[2]

for (i in 1:length(SSP_list)) 
{
  # i <- 1
  
  SSP <- SSP_list[i]
  
  if (SSP == "ssp245") {
    filepath <- "./input_data/LUH2/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-MESSAGE-ssp245-2-1-f_gn_2015-2100.nc"
  } else {
    filepath <- "./input_data/LUH2/multiple-states_input4MIPs_landState_ScenarioMIP_UofMD-MAGPIE-ssp585-2-1-f_gn_2015-2100.nc"
  }
  
  # Get raster Brick from NetCDF file for only one variable
  nc_raster_primf <- brick(filepath, var = "primf") # Primary forests
  nc_raster_secdf <- brick(filepath, var = "secdf") # Secondary forests
  nlayers(nc_raster_primf) # As many layers as time step
  
  # Sum the two forest covers
  forest_cover <- nc_raster_primf + nc_raster_secdf
  
  # Remove rasterBrick per variable to save RAM space
  rm(nc_raster_primf, nc_raster_secdf)
  
  # Extract the time step from 2061 to 2080. File start in 2015
  forest_cover_61_80 <- subset(forest_cover, 46:65)
  names(forest_cover_61_80) <- as.character(2061:2080)
  
  plot(forest_cover_61_80[[1]])
  
  # Crop and resample to Neotropics
  forest_cover_61_80_cropped <- resample(x = forest_cover_61_80, y = envData, method = "bilinear") # Ajust extent and resolution to match the one of our SDM maps
  
  # Compute mean
  forest_cover_2070_cropped <- calc(forest_cover_61_80_cropped, fun = mean, na.rm = T)
  # Compute interannual variation over 20 years
  forest_cover_2070_cropped_sd <- calc(forest_cover_61_80_cropped, fun = sd, na.rm = T)
  
  plot(forest_cover_2070_cropped)
  plot(forest_cover_2070_cropped_sd)
  
  # Save rasters
  saveRDS(forest_cover_2070_cropped, file = paste0("./input_data/LUH2/Forest_cover_",target,"_",SSP,".rds"), version = "2")
  saveRDS(forest_cover_2070_cropped_sd, file = paste0("./input_data/LUH2/Forest_cover_",target,"_sd_",SSP,".rds"), version = "2")
}



### 3/ Compute changes ####

rm(list = ls())

# Set Time period and target year
period <- "2061-2080" ; target <- "2070"

# Load forest covers
forest_cover_baseline <- readRDS(file = paste0("./input_data/LUH2/Forest_cover_baseline.rds"))
forest_cover_2070_ssp245 <- readRDS(file = paste0("./input_data/LUH2/Forest_cover_",target,"_ssp245.rds"))
forest_cover_2070_ssp585 <- readRDS(file = paste0("./input_data/LUH2/Forest_cover_",target,"_ssp585.rds"))

# Compute changes
forest_cover_raw_changes_ssp245 <- (forest_cover_2070_ssp245 - forest_cover_baseline) * 100
names(forest_cover_raw_changes_ssp245) <- "Absolute changes - SSP2-RCP4.5"
forest_cover_raw_changes_ssp585 <- (forest_cover_2070_ssp585 - forest_cover_baseline) * 100
names(forest_cover_raw_changes_ssp585) <- "Absolute changes - SSP5-RCP8.5"

forest_cover_rel_changes_ssp245 <- forest_cover_raw_changes_ssp245 / forest_cover_baseline
names(forest_cover_rel_changes_ssp245) <- "Relative changes - SSP2-RCP4.5"
forest_cover_rel_changes_ssp585 <- forest_cover_raw_changes_ssp585 / forest_cover_baseline
names(forest_cover_rel_changes_ssp585) <- "Relative changes - SSP5-RCP8.5"

# Save
saveRDS(forest_cover_raw_changes_ssp245, file = paste0("./outputs/LUH2/Forest_cover_changes_",target,"_ssp245.rds"))
saveRDS(forest_cover_raw_changes_ssp585, file = paste0("./outputs/LUH2/Forest_cover_changes_",target,"_ssp585.rds"))
saveRDS(forest_cover_rel_changes_ssp245, file = paste0("./outputs/LUH2/Forest_cover_changes_",target,"_ssp245.rds"))
saveRDS(forest_cover_rel_changes_ssp585, file = paste0("./outputs/LUH2/Forest_cover_changes_",target,"_ssp585.rds"))

### 4/ Contrast maps ####

library(tmaptools)
pal_grn_red <- rev(get_brewer_pal("RdYlGn", n = 200))
pal_red_grn <- get_brewer_pal("RdYlGn", n = 200)

continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_", res, ".rds"))

# SSP245 - Raw
plot(forest_cover_raw_changes_ssp245)

hist(forest_cover_raw_changes_ssp245)

forest_cover_raw_changes_ssp245_contrasted <- calc(x = forest_cover_raw_changes_ssp245, 
                                                   fun = function(x){ x[x < -30] <- -30 ; x[x > 30] <- 30 ; return(x)} )

plot(forest_cover_raw_changes_ssp245_contrasted, col = pal_red_grn)

saveRDS(forest_cover_raw_changes_ssp245_contrasted, file = paste0("./outputs/LUH2/Forest_cover_raw_changes_",target,"_ssp245_contrasted.rds"), version = "2")


# SSP585 - Raw
plot(forest_cover_raw_changes_ssp585)

hist(forest_cover_raw_changes_ssp585)

forest_cover_raw_changes_ssp585_contrasted <- calc(x = forest_cover_raw_changes_ssp585, 
                                                   fun = function(x){ x[x < -30] <- -30 ; x[x > 30] <- 30 ; return(x)} )

plot(forest_cover_raw_changes_ssp585_contrasted, col = pal_red_grn)

saveRDS(forest_cover_raw_changes_ssp585_contrasted, file = paste0("./outputs/LUH2/Forest_cover_raw_changes_",target,"_ssp585_contrasted.rds"), version = "2")

# SSP245 - Relative
plot(forest_cover_rel_changes_ssp245)

hist(forest_cover_rel_changes_ssp245)

forest_cover_rel_changes_ssp245_contrasted <- calc(x = forest_cover_rel_changes_ssp245, 
                                                   fun = function(x){ x[x < -100] <- -100 ; x[x > 100] <- 100 ; return(x)} )

# Fill NA values for terrestrial pixel with 0 (case of pixel with 0% forest as baseline)
temp <- continent_mask
temp[!is.na(forest_cover_rel_changes_ssp245_contrasted[])] <- forest_cover_rel_changes_ssp245_contrasted[!is.na(forest_cover_rel_changes_ssp245_contrasted[])]
forest_cover_rel_changes_ssp245_contrasted <- temp

plot(forest_cover_rel_changes_ssp245_contrasted, col = pal_red_grn)

saveRDS(forest_cover_rel_changes_ssp245_contrasted, file = paste0("./outputs/LUH2/Forest_cover_rel_changes_",target,"_ssp245_contrasted.rds"), version = "2")


# SSP585 - Relative
plot(forest_cover_rel_changes_ssp585)

hist(forest_cover_rel_changes_ssp585)

forest_cover_rel_changes_ssp585_contrasted <- calc(x = forest_cover_rel_changes_ssp585, 
                                                   fun = function(x){ x[x < -100] <- -100 ; x[x > 100] <- 100 ; return(x)} )

# Fill NA values for terrestrial pixel with 0 (case of pixel with 0% forest as baseline)
temp <- continent_mask
temp[!is.na(forest_cover_rel_changes_ssp585_contrasted[])] <- forest_cover_rel_changes_ssp585_contrasted[!is.na(forest_cover_rel_changes_ssp585_contrasted[])]
forest_cover_rel_changes_ssp585_contrasted <- temp

plot(forest_cover_rel_changes_ssp585_contrasted, col = pal_red_grn)

saveRDS(forest_cover_rel_changes_ssp585_contrasted, file = paste0("./outputs/LUH2/Forest_cover_rel_changes_",target,"_ssp585_contrasted.rds"), version = "2")


### 5/ Plot full maps ####

library(rangeBuilder)

# Load contrasted maps
forest_cover_raw_changes_ssp245_contrasted <- readRDS(file = paste0("./outputs/LUH2/Forest_cover_raw_changes_",target,"_ssp245_contrasted.rds"))
forest_cover_rel_changes_ssp245_contrasted <- readRDS(file = paste0("./outputs/LUH2/Forest_cover_rel_changes_",target,"_ssp245_contrasted.rds"))
forest_cover_raw_changes_ssp585_contrasted <- readRDS(file = paste0("./outputs/LUH2/Forest_cover_raw_changes_",target,"_ssp585_contrasted.rds"))
forest_cover_rel_changes_ssp585_contrasted <- readRDS(file = paste0("./outputs/LUH2/Forest_cover_rel_changes_",target,"_ssp585_contrasted.rds"))

crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".rds"))

# 5.1/ Absolute and relative changes ####

pdf(file = paste0("./maps/LUH2/All_ssp_forest_cover_changes_",period,".pdf"), height = 10, width = 12)

par(mfrow=(c(2,2)))

image(forest_cover_raw_changes_ssp245_contrasted, col = pal_red_grn, legend = T,
      ylab = "", xlab = "", main = "Absolute - SSP2-RCP4.5",
      cex.axis = 1.4, cex.main = 1.6)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(forest_cover_raw_changes_ssp245_contrasted, locs = seq(-30, 30, 10), cex.axis = 1.1, ramp = pal_red_grn, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 5, font = 2, cex = 1.2, label = "Changes (%)")

image(forest_cover_raw_changes_ssp585_contrasted, col = pal_red_grn, legend = T,
      ylab = "", xlab = "", main = "Absolute - SSP5-RCP8.5",
      cex.axis = 1.4, cex.main = 1.6)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(forest_cover_raw_changes_ssp245_contrasted, locs = seq(-30, 30, 10), cex.axis = 1.1, ramp = pal_red_grn, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 5, font = 2, cex = 1.2, label = "Changes (%)")

image(forest_cover_rel_changes_ssp245_contrasted, col = pal_red_grn, legend = T,
      ylab = "", xlab = "", main = "Relative - SSP2-RCP4.5",
      cex.axis = 1.4, cex.main = 1.6)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(forest_cover_rel_changes_ssp245_contrasted, locs = seq(-100, 100, 50), cex.axis = 1.1, ramp = pal_red_grn, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 5, font = 2, cex = 1.2, label = "Changes (%)")

image(forest_cover_rel_changes_ssp585_contrasted, col = pal_red_grn, legend = T,
      ylab = "", xlab = "", main = "Relative - SSP5-RCP8.5",
      cex.axis = 1.4, cex.main = 1.6)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(forest_cover_rel_changes_ssp245_contrasted, locs = seq(-100, 100, 50), cex.axis = 1.1, ramp = pal_red_grn, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 5, font = 2, cex = 1.2, label = "Changes (%)")

dev.off()

par(mfrow=(c(1,1)))


# 5.2/ Combine with indices and climate anomalies for contrasted plots ####

library(tmaptools)

# Color palet for diversity indices
# cool = rainbow(49, s = 1, v = 1, start=rgb2hsv(col2rgb('yellow'))[1], end=rgb2hsv(col2rgb('blue'))[1])
# warm = rainbow(50, s = 1, v= 1, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
# pal_bl_red  = c(rev(cool), rev(warm))
# pal_bl_red <- c(gplots::col2hex("grey93"), pal_bl_red)

# New color palette
pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")

pal_grn_red <- rev(get_brewer_pal("RdYlGn", n = 200))
pal_red_grn <- get_brewer_pal("RdYlGn", n = 200)

load(file = "./input_data/Map_stuff/bg_mask.RData") # Load bg shp
load(file = "./input_data/Map_stuff/country_borders.RData") # Load country borders
rivers <- st_read("./input_data/Map_stuff/Major_rivers/MajorRivers.shp") # Load rivers

crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".rds"))

# Load Indices maps
tot.sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
MPD.raster_Jaccard.80_contrasted <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80_contrasted.rds")

# Load climate anomalies maps
SSP245_anomalies_base_std_all_median_contrasted <- readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_all_median_contrasted_",period,"_ssp245.rds"))
SSP585_anomalies_base_std_all_median_contrasted <- readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_all_median_contrasted_",period,"_ssp585.rds"))

# Load forest cover absolute changes maps
forest_cover_raw_changes_ssp245_contrasted <- readRDS(file = paste0("./outputs/LUH2/Forest_cover_raw_changes_",target,"_ssp245_contrasted.rds"))
forest_cover_raw_changes_ssp585_contrasted <- readRDS(file = paste0("./outputs/LUH2/Forest_cover_raw_changes_",target,"_ssp585_contrasted.rds"))

### 5.2.1/ Full extent ####

pdf(file = paste0("./maps/LUH2/All_ssp_forest_cover_and_climate_changes_",period,".pdf"), height = 13, width = 10)

par(mfrow=(c(3,2)))

image(tot.sp.richness_Jaccard.80, col = pal_bl_red_Mannion, legend = T,
      ylab = "", xlab = "", main = "Species richness",
      cex.axis = 1.4, cex.main = 1.6)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(tot.sp.richness_Jaccard.80, locs = seq(0, 120, 20), cex.axis = 1.1, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 5, font = 2, cex = 1.1, label = "Species")

image(MPD.raster_Jaccard.80_contrasted, col = pal_bl_red_Mannion, legend = T,
      ylab = "", xlab = "", main = "Mean pairwise Phylogenetic Distance",
      cex.axis = 1.4, cex.main = 1.6)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(MPD.raster_Jaccard.80_contrasted, locs = seq(38, 44, 2), cex.axis = 1.1, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6.5, font = 2, cex = 1.1, label = "Evolutionary\nTime (My)")

image(SSP245_anomalies_base_std_all_median_contrasted, col = pal_grn_red, legend = T,
      ylab = "", xlab = "", main = "Climate anomalies - SSP2-RCP4.5",
      cex.axis = 1.4, cex.main = 1.6)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(SSP245_anomalies_base_std_all_median_contrasted, locs = seq(0.45, 0.7, 0.05), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6.5, font = 2, cex = 1.1, label = "Standardized\nanomalies")

image(SSP585_anomalies_base_std_all_median_contrasted, col = pal_grn_red, legend = T,
      ylab = "", xlab = "", main = "Climate anomalies - SSP5-RCP8.5",
      cex.axis = 1.4, cex.main = 1.6)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(SSP585_anomalies_base_std_all_median_contrasted, locs = seq(0.6, 1.1, 0.1), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6.5, font = 2, cex = 1.1, label = "Standardized\nanomalies")

image(forest_cover_raw_changes_ssp245_contrasted, col = pal_red_grn, legend = T,
      ylab = "", xlab = "", main = "Forest cover - SSP2-RCP4.5",
      cex.axis = 1.4, cex.main = 1.6)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(forest_cover_raw_changes_ssp245_contrasted, locs = seq(-30, 30, 10), cex.axis = 1.1, ramp = pal_red_grn, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6.5, font = 2, cex = 1.1, label = "Absolute\nchanges (%)")

image(forest_cover_raw_changes_ssp585_contrasted, col = pal_red_grn, legend = T,
      ylab = "", xlab = "", main = "Forest cover - SSP5-RCP8.5",
      cex.axis = 1.4, cex.main = 1.6)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(forest_cover_raw_changes_ssp245_contrasted, locs = seq(-30, 30, 10), cex.axis = 1.1, ramp = pal_red_grn, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 6.5, font = 2, cex = 1.1, label = "Absolute\nchanges (%)")

dev.off()

par(mfrow=(c(1,1)))


### 5.2.2/ Plot with focus on Andes and comparison with diversity indices ####

### Load environmental stack to use as template for resolution and extent
res <- "15"
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))

# Get elevation contours
Andes_ext <- extent(c(-90, -59, -15, 14))
Close_Andes_ext <- extent(c(-90, -67.2, -15, 14))
DEM <- envData[["Elevation"]]
Andes_DEM <- crop(DEM, Andes_ext)
Close_Andes_DEM <- crop(DEM, Close_Andes_ext)

# Get Mata bbox rectangle
Andes_bbox_sp <- SpatialPolygons(list(Polygons(list(Polygon(coords = matrix(data = c(-90, 14, -59, 14, -59, -15, -90, -15, -90, 14),
                                                                           ncol = 2, byrow = T), hole = F))
                                              , ID = 1))
                                , proj4string = crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
# Andes_bg_mask <- gIntersection(Andes_bbox_sp, bg_mask)
Andes_bg_mask <- gDifference(Andes_bbox_sp, crop_mask_shp)

# Load boundaries
crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".rds"))

# Load Indices maps
tot.sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
MPD.raster_Jaccard.80_contrasted <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80_contrasted.rds")

# Load climate anomalies maps
SSP245_anomalies_base_std_all_median_contrasted <- readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_all_median_contrasted_",period,"_ssp245.rds"))
SSP585_anomalies_base_std_all_median_contrasted <- readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_all_median_contrasted_",period,"_ssp585.rds"))

# Load forest cover absolute changes maps
forest_cover_raw_changes_ssp245_contrasted <- readRDS(file = paste0("./outputs/LUH2/Forest_cover_raw_changes_",target,"_ssp245_contrasted.rds"))
forest_cover_raw_changes_ssp585_contrasted <- readRDS(file = paste0("./outputs/LUH2/Forest_cover_raw_changes_",target,"_ssp585_contrasted.rds"))


# Crop maps and remove pixels not computed for indices maps
Andes_SSP245_anomalies_base_std_all_median_contrasted <- crop(SSP245_anomalies_base_std_all_median_contrasted, Andes_ext)
Andes_SSP245_anomalies_base_std_all_median_contrasted <- mask(x = Andes_SSP245_anomalies_base_std_all_median_contrasted, mask = Andes_DEM)
Andes_SSP585_anomalies_base_std_all_median_contrasted <- crop(SSP585_anomalies_base_std_all_median_contrasted, Andes_ext)
Andes_SSP585_anomalies_base_std_all_median_contrasted <- mask(x = Andes_SSP585_anomalies_base_std_all_median_contrasted, mask = Andes_DEM)

Andes_forest_cover_raw_changes_ssp245_contrasted <- crop(forest_cover_raw_changes_ssp245_contrasted, Andes_ext)
Andes_forest_cover_raw_changes_ssp245_contrasted <- mask(x = Andes_forest_cover_raw_changes_ssp245_contrasted, mask = Andes_DEM)
Andes_forest_cover_raw_changes_ssp585_contrasted <- crop(forest_cover_raw_changes_ssp585_contrasted, Andes_ext)
Andes_forest_cover_raw_changes_ssp585_contrasted <- mask(x = Andes_forest_cover_raw_changes_ssp585_contrasted, mask = Andes_DEM)

Andes.tot.sp.richness_Jaccard.80 <- crop(tot.sp.richness_Jaccard.80, Andes_ext)
Andes.MPD.raster_Jaccard.80_contrasted <- crop(MPD.raster_Jaccard.80_contrasted, Andes_ext)


pdf(file = paste0("./maps/LUH2/All_ssp_forest_cover_and_climate_changes_",period,"_Andes_2.pdf"), height = 13, width = 10)

par(mfrow=(c(3,2)))

image(Andes.tot.sp.richness_Jaccard.80, col = pal_bl_red_Mannion, legend = T,
      ylab = "", xlab = "", main = "Species richness",
      cex.axis = 1.4, cex.main = 1.6)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
        drawlabels = F, col = "black", add = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(Andes_bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
scalebar(d = 750, type = "line", lwd = 4, divs = 4, xy = c(-67, -13.5), label = c("", "750 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.6, padin = c(0.2, 0.4), pos = "topleft")
addRasterLegend(Andes.tot.sp.richness_Jaccard.80, locs = seq(20, 120, 20), cex.axis = 1.1, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.5, 4.5))
graphics::text(x = -86.7, y = 6, font = 2, cex = 1.1, label = "Species")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 1.1, text.font = 2)
graphics::text(x = -86.7, y = -10.5, font = 2, cex = 1.1, label = "Elevation")

image(Andes.MPD.raster_Jaccard.80_contrasted, col = pal_bl_red_Mannion, legend = T,
      ylab = "", xlab = "", main = "Mean pairwise Phylogenetic Distance",
      cex.axis = 1.4, cex.main = 1.6)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
        drawlabels = F, col = "black", add = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(Andes.MPD.raster_Jaccard.80_contrasted, locs = seq(38, 44, 2), cex.axis = 1.1, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.5, 4.5))
graphics::text(x = -86.5, y = 6.5, font = 2, cex = 1.1, label = "Evolutionary\nTime (My)")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 1.1, text.font = 2)
graphics::text(x = -86.7, y = -10.5, font = 2, cex = 1.1, label = "Elevation")

image(Andes_SSP_anomalies_base_std_all_median_indep_contrasted[[1]], col = pal_grn_red, legend = T,
      ylab = "", xlab = "", main = "Climate anomalies - SSP2-RCP4.5",
      cex.axis = 1.4, cex.main = 1.6)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
        drawlabels = F, col = "black", add = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(Andes_SSP_anomalies_base_std_all_median_indep_contrasted[[1]], locs = seq(0.45, 0.7, 0.05), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.5, 4.5))
graphics::text(x = -86.5, y = 6.5, font = 2, cex = 1.1, label = "Standardized\nanomalies")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 1.1, text.font = 2)
graphics::text(x = -86.7, y = -10.5, font = 2, cex = 1.1, label = "Elevation")

image(Andes_SSP_anomalies_base_std_all_median_indep_contrasted[[2]], col = pal_grn_red, legend = T,
      ylab = "", xlab = "", main = "Climate anomalies - SSP5 - RCP8.5",
      cex.axis = 1.4, cex.main = 1.6)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
        drawlabels = F, col = "black", add = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(Andes_SSP_anomalies_base_std_all_median_indep_contrasted[[2]], locs = seq(0.6, 1.1, 0.1), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.5, 4.5))
graphics::text(x = -86.5, y = 6.5, font = 2, cex = 1.1, label = "Standardized\nanomalies")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 1.1, text.font = 2)
graphics::text(x = -86.7, y = -10.5, font = 2, cex = 1.1, label = "Elevation")

image(Andes_forest_cover_raw_changes_ssp245_contrasted, col = pal_red_grn, legend = T,
      ylab = "", xlab = "", main = "Forest cover - SSP2-RCP4.5",
      cex.axis = 1.4, cex.main = 1.6)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
        drawlabels = F, col = "black", add = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(Andes_forest_cover_raw_changes_ssp245_contrasted, locs = seq(-30, 30, 10), cex.axis = 1.1, ramp = pal_red_grn, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.5, 4.5))
graphics::text(x = -86.5, y = 6.5, font = 2, cex = 1.1, label = "Absolute\nchanges (%)")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 1.1, text.font = 2)
graphics::text(x = -86.7, y = -10.5, font = 2, cex = 1.1, label = "Elevation")

image(Andes_forest_cover_raw_changes_ssp585_contrasted, col = pal_red_grn, legend = T,
      ylab = "", xlab = "", main = "Forest cover - SSP5 - RCP8.5", zlim = c(-30, 30),
      cex.axis = 1.4, cex.main = 1.6)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
        drawlabels = F, col = "black", add = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(Andes_forest_cover_raw_changes_ssp245_contrasted, locs = seq(-30, 30, 10), cex.axis = 1.1, ramp = pal_red_grn, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.5, 4.5))
graphics::text(x = -86.5, y = 6.5, font = 2, cex = 1.1, label = "Absolute\nchanges (%)")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 1.1, text.font = 2)
graphics::text(x = -86.7, y = -10.5, font = 2, cex = 1.1, label = "Elevation")

dev.off()

par(mfrow=(c(1,1)))


### 5.2.3/ Plot with focus on Mata Atlantica and comparison with diversity indices ####

### Load environmental stack to use as template for resolution and extent
res <- "15"
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))

# Get extent
Mata_ext <- extent(c(-61, -31, -31, -4))

# Get elevation
DEM <- envData[["Elevation"]]
Mata_DEM <- crop(DEM, Mata_ext)

# Get Mata bbox rectangle
Mata_bbox_sp <- SpatialPolygons(list(Polygons(list(Polygon(coords = matrix(data = c(-61, -4, -31, -4, -31, -31, -61, -31, -61, -4),
                                                                      ncol = 2, byrow = T), hole = F))
                                         , ID = 1))
                           , proj4string = crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
Mata_bg_mask <- gIntersection(Mata_bbox_sp, bg_mask)
Mata_bg_mask <- gDifference(Mata_bbox_sp, crop_mask_shp)


# Load Mata Atlantica shp
# Mata_borders <- st_read("./input_data/Map_stuff/mata_atlantica_border/mata_atlantica_border.shp") # Load Mata borders
Mata_borders <- readOGR(dsn = "./input_data/Map_stuff/mata_atlantica_border", layer = "mata_atlantica_border")

# Load boundaries
crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".rds"))

# Load Indices maps
tot.sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
MPD.raster_Jaccard.80_contrasted <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80_contrasted.rds")

# Load climate anomalies maps
SSP245_anomalies_base_std_all_median_contrasted <- readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_all_median_contrasted_",period,"_ssp245.rds"))
SSP585_anomalies_base_std_all_median_contrasted <- readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_all_median_contrasted_",period,"_ssp585.rds"))

# Load forest cover absolute changes maps
forest_cover_raw_changes_ssp245_contrasted <- readRDS(file = paste0("./outputs/LUH2/Forest_cover_raw_changes_",target,"_ssp245_contrasted.rds"))
forest_cover_raw_changes_ssp585_contrasted <- readRDS(file = paste0("./outputs/LUH2/Forest_cover_raw_changes_",target,"_ssp585_contrasted.rds"))


# Crop maps and remove pixels not computed for indices maps
Mata_SSP245_anomalies_base_std_all_median_contrasted <- crop(SSP245_anomalies_base_std_all_median_contrasted, Mata_ext)
Mata_SSP245_anomalies_base_std_all_median_contrasted <- mask(x = Mata_SSP245_anomalies_base_std_all_median_contrasted, mask = Mata_DEM)
Mata_SSP585_anomalies_base_std_all_median_contrasted <- crop(SSP585_anomalies_base_std_all_median_contrasted, Mata_ext)
Mata_SSP585_anomalies_base_std_all_median_contrasted <- mask(x = Mata_SSP585_anomalies_base_std_all_median_contrasted, mask = Mata_DEM)

Mata_forest_cover_raw_changes_ssp245_contrasted <- crop(forest_cover_raw_changes_ssp245_contrasted, Mata_ext)
Mata_forest_cover_raw_changes_ssp245_contrasted <- mask(x = Mata_forest_cover_raw_changes_ssp245_contrasted, mask = Mata_DEM)
Mata_forest_cover_raw_changes_ssp585_contrasted <- crop(forest_cover_raw_changes_ssp585_contrasted, Mata_ext)
Mata_forest_cover_raw_changes_ssp585_contrasted <- mask(x = Mata_forest_cover_raw_changes_ssp585_contrasted, mask = Mata_DEM)

Mata.tot.sp.richness_Jaccard.80 <- crop(tot.sp.richness_Jaccard.80, Mata_ext)
Mata.MPD.raster_Jaccard.80_contrasted <- crop(MPD.raster_Jaccard.80_contrasted, Mata_ext)


pdf(file = paste0("./maps/LUH2/All_ssp_forest_cover_and_climate_changes_",period,"_Mata.pdf"), height = 13, width = 10)

par(mfrow=(c(3,2)))

image(Mata.tot.sp.richness_Jaccard.80, col = pal_bl_red_Mannion, legend = T,
      ylab = "", xlab = "", main = "Species richness",
      cex.axis = 1.4, cex.main = 1.6)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(Mata_bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Mata_borders, lwd = 2, border = "darkgreen", add = T)
plot(Mata_bbox_sp, lwd = 1, border = "grey20", add = T)
scalebar(d = 500, type = "line", lwd = 4, divs = 4, xy = c(-59, -7), label = c("", "500 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15))
addRasterLegend(Mata.tot.sp.richness_Jaccard.80, locs = seq(5, 35, 5), cex.axis = 1.1, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-33.22, -32.2, -28, -14))
graphics::text(x = -33.1, y = -12.5, font = 2, cex = 1.1, label = "Species")
legend(x = "bottomright", legend = c("Mata Atlantica"), lty = c(1) , col = "darkgreen", lwd = 2, bty = "n", cex = 1.1, text.font = 2)

image(Mata.MPD.raster_Jaccard.80_contrasted, col = pal_bl_red_Mannion, legend = T,
      ylab = "", xlab = "", main = "Mean pairwise Phylogenetic Distance",
      cex.axis = 1.4, cex.main = 1.6)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(Mata_bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Mata_borders, lwd = 2, border = "darkgreen", add = T)
plot(Mata_bbox_sp, lwd = 1, border = "grey20", add = T)
scalebar(d = 500, type = "line", lwd = 4, divs = 4, xy = c(-59, -7), label = c("", "500 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15))
addRasterLegend(Mata.MPD.raster_Jaccard.80_contrasted, locs = seq(38, 44, 2), cex.axis = 1.1, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-33.22, -32.2, -28, -14))
graphics::text(x = -34.1, y = -12, font = 2, cex = 1.1, label = "Evolutionary\nTime (My)")
legend(x = "bottomright", legend = c("Mata Atlantica"), lty = c(1) , col = "darkgreen", lwd = 2, bty = "n", cex = 1.1, text.font = 2)

image(Mata_SSP245_anomalies_base_std_all_median_contrasted, col = pal_grn_red, legend = T,
      ylab = "", xlab = "", main = "Climate anomalies - SSP2-RCP4.5",
      cex.axis = 1.4, cex.main = 1.6)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(Mata_bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Mata_borders, lwd = 2, border = "darkgreen", add = T)
plot(Mata_bbox_sp, lwd = 1, border = "grey20", add = T)
scalebar(d = 500, type = "line", lwd = 4, divs = 4, xy = c(-59, -7), label = c("", "500 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15))
addRasterLegend(Mata_SSP245_anomalies_base_std_all_median_contrasted, locs = seq(0.45, 0.65, 0.05), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-33.22, -32.2, -28, -14))
graphics::text(x = -34, y = -12, font = 2, cex = 1.1, label = "Standardized\nanomalies")
legend(x = "bottomright", legend = c("Mata Atlantica"), lty = c(1) , col = "darkgreen", lwd = 2, bty = "n", cex = 1.1, text.font = 2)


image(Mata_SSP585_anomalies_base_std_all_median_contrasted, col = pal_grn_red, legend = T,
      ylab = "", xlab = "", main = "Climate anomalies - SSP5-RCP8.5",
      cex.axis = 1.4, cex.main = 1.6)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(Mata_bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Mata_borders, lwd = 2, border = "darkgreen", add = T)
plot(Mata_bbox_sp, lwd = 1, border = "grey20", add = T)
scalebar(d = 500, type = "line", lwd = 4, divs = 4, xy = c(-59, -7), label = c("", "500 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15))
addRasterLegend(Mata_SSP585_anomalies_base_std_all_median_contrasted, locs = seq(0.6, 1.1, 0.1), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-33.22, -32.2, -28, -14))
graphics::text(x = -34, y = -12, font = 2, cex = 1.1, label = "Standardized\nanomalies")
legend(x = "bottomright", legend = c("Mata Atlantica"), lty = c(1) , col = "darkgreen", lwd = 2, bty = "n", cex = 1.1, text.font = 2)


image(Mata_forest_cover_raw_changes_ssp245_contrasted, col = pal_red_grn, legend = T,
      ylab = "", xlab = "", main = "Forest cover - SSP2-RCP4.5",
      cex.axis = 1.4, cex.main = 1.6)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(Mata_bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Mata_borders, lwd = 2, border = "darkgreen", add = T)
plot(Mata_bbox_sp, lwd = 1, border = "grey20", add = T)
scalebar(d = 500, type = "line", lwd = 4, divs = 4, xy = c(-59, -7), label = c("", "500 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15))
addRasterLegend(Mata_forest_cover_raw_changes_ssp245_contrasted, locs = seq(-30, 30, 10), cex.axis = 1.1, ramp = pal_red_grn, ncolors = 200, border = T, location = c(-33.22, -32.2, -28, -14))
graphics::text(x = -33.8, y = -11.5, font = 2, cex = 1.1, label = "Absolute\nchanges (%)")
legend(x = "bottomright", legend = c("Mata Atlantica"), lty = c(1) , col = "darkgreen", lwd = 2, bty = "n", cex = 1.1, text.font = 2)

image(Mata_forest_cover_raw_changes_ssp585_contrasted, col = pal_red_grn, legend = T,
      ylab = "", xlab = "", main = "Forest cover - SSP5-RCP8.5",
      cex.axis = 1.4, cex.main = 1.6)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(Mata_bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Mata_borders, lwd = 2, border = "darkgreen", add = T)
plot(Mata_bbox_sp, lwd = 1, border = "grey20", add = T)
scalebar(d = 500, type = "line", lwd = 4, divs = 4, xy = c(-59, -7), label = c("", "500 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15))
addRasterLegend(Mata_forest_cover_raw_changes_ssp585_contrasted, locs = seq(-30, 30, 10), cex.axis = 1.1, ramp = pal_red_grn, ncolors = 200, border = T, location = c(-33.22, -32.2, -28, -14))
graphics::text(x = -33.8, y = -11.5, font = 2, cex = 1.1, label = "Absolute\nchanges (%)")
legend(x = "bottomright", legend = c("Mata Atlantica"), lty = c(1) , col = "darkgreen", lwd = 2, bty = "n", cex = 1.1, text.font = 2)

dev.off()

par(mfrow=(c(1,1)))
