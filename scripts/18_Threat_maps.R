
###### Script 18: Threat maps for diversity ######

# Compute maps of threats from HII weighted by diversity indices

# Only for Jaccard.80

### Use only 9 indices

# 1/ Species richness
# 2/ Species Shannon's Diversity with compatibility with species richness, with Jost's transformation (exponential)
# 3/ Mean Continuous Species Rarity = Range-size Weighted Species Richness

# 4/ Faith's Phylogenetic Diversity
# 5/ Sum of Fair-proportions
# 6/ Mean pairwise Phylogenetic Distance

# 7/ Mimicry richness
# 8/ Mimicry Shannon's Diversity with compatibility with mimicry richness, with Jost's transformation (exponential)
# 9/ Mean Continuous Mimicry rarity = Range size weighted mimicry richness

###


# Inputs: 
   # Indices maps
   # Threat map: HII (HF, Venter et al., 2016)

# Outputs:
   # Weighted maps of threat for each indice, normalized version
   # Maps of indices, weighted by threats

###


#####

# Effacer l'environnement
rm(list = ls())

library(raster)
library(rangeBuilder)
library(gplots)

# Change temp folder for raster files
rasterOptions(tmpdir = paste0("./temp"))
# To clean regularly temp folder
unlink(list.files(path = rasterOptions()$tmpdir, all.files = T, recursive = T, include.dirs = T, full.names = T), force = T, recursive = T) # Clean raster store in temp

# Load shp files for maps
load(file = "./input_data/Map_stuff/bg_mask.RData") # Load bg shp
load(file = "./input_data/Map_stuff/country_borders.RData") # Load country borders

# New color palette
pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")

library(tmaptools)
pal_grn_red_NA <- pal_grn_red <- rev(get_brewer_pal("RdYlGn", n = 200))
pal_grn_red_NA[1] <- "#EDEDED"


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

# Load stuff for plot on Mata Atlantica
Mata_ext <- extent(c(-61, -31, -31, -4))
Mata_bbox <- as(Mata_ext, 'SpatialPolygons') 
Mata_borders <- readOGR(dsn = "./input_data/Map_stuff/mata_atlantica_border", layer = "mata_atlantica_border")


### 0/ Retrieve and plot HII map ####

# HII <-  raster(x = "./input_data/HII/wildareas-v3-2009-human-footprint-geotiff/WGS84_HumanFootprint_Venter_2009_int.tif")
# HII <- crop(x = HII, y = continent_mask) # Crop first to ease the computation process of resampling
# HII <- resample(x = HII, y = continent_mask, method = "bilinear") # Resample to modeling resolution (i.e., 15')
# HII <- mask(x = HII, mask = continent_mask) # Keep only communities with data
# HII <- round(HII*2, 1) # Put on a scale of 0-100
# 
# hist(HII)
# 
# save(HII, file = "./input_data/HII/HII.RData", version = "2")
# saveRDS(HII, file = "./input_data/HII/HII.rds", version = "2")

# Load HII map
HII <- readRDS(file = "./input_data/HII/HII.rds")
HII_norm <- HII/HII@data@max

# # Plot a contrasted HII ####
# HII_contrasted <- HII
# HII_contrasted[HII[] > 30] <- 30
# 
# save(HII_contrasted, file = "./input_data/HII/HII_contrasted.RData", version = "2")
# saveRDS(HII_contrasted, file = "./input_data/HII/HII_contrasted.rds", version = "2")
# 
# pdf(file = paste0("./maps/Threat_maps/HII.pdf"), height = 5.3, width = 6.5)
# internal_margins <- par()$mar
# par(mar = c(3.1, 3.5,3.5,2.1))
# image(HII_contrasted, col = pal_grn_red, main = paste0("Human Influence Index"),
#       cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
#       legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# # abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# # abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# # abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
# addRasterLegend(HII_contrasted, locs = seq(0, 30, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
# graphics::text(x = -112, y = 4.5, font = 2, cex = 1.2, label = "HII")
# par(mar = internal_margins)
# dev.off()
# 
# # Plot zoom on Andes
# Andes.HII_contrasted <- crop(HII_contrasted, Andes_ext)
# 
# save(Andes.HII_contrasted, file = "./input_data/HII/Andes.HII_contrasted.RData", version = "2")
# saveRDS(Andes.HII_contrasted, file = "./input_data/HII/Andes.HII_contrasted.rds", version = "2")
# 
# pdf(file = paste0("./maps/Threat_maps/Andes/Andes.HII.pdf"), height = 5.3, width = 6.5)
# internal_margins <- par()$mar
# par(mar = c(3.1, 3.5,3.5,2.1))
# image(Andes.HII_contrasted, col = pal_grn_red, main = paste0("Human Influence Index"),
#       cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
#       legend  = T)
# # plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# plot(Andes_bbox, lwd = 1, border = "grey20", add = T)
# contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
#         drawlabels = F, col = "black", add = T)
# addRasterLegend(Andes.HII_contrasted, locs = seq(0, 30, 5), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.75, 4.25))
# graphics::text(x = -87, y = 6, font = 2, cex = 0.9, label = "HII")
# legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 0.9, text.font = 2)
# graphics::text(x = -86.7, y = -10.8, font = 2, cex = 0.9, label = "Elevation")
# scalebar(d = 750, type = "line", lwd = 4, divs = 4, xy = c(-67, -13.5), label = c("", "750 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.6, padin = c(0.2, 0.4), pos = "topleft")
# par(mar = internal_margins)
# dev.off()
# 
# # Plot zoom on Mata Atlantica
# Mata.HII <- crop(HII, Mata_ext)
# 
# hist(Mata.HII)
# 
# Mata.HII_contrasted <- Mata.HII
# Mata.HII_contrasted[Mata.HII[] > 30] <- 30
# 
# save(Mata.HII_contrasted, file = "./input_data/HII/Mata.HII_contrasted.RData", version = "2")
# saveRDS(Mata.HII_contrasted, file = "./input_data/HII/Mata.HII_contrasted.rds", version = "2")
# 
# 
# pdf(file = paste0("./maps/Threat_maps/Mata/Mata.HII.pdf"), height = 5.3, width = 6.5)
# internal_margins <- par()$mar
# par(mar = c(3.1, 3.5,3.5,2.1))
# image(Mata.HII_contrasted, col = pal_grn_red, main = paste0("Human Influence Index"),
#       cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
#       ylab = "", xlab = "",
#       legend.args=list(text="     HII", cex=1.2, line = 1, font = 2),
#       legend  = T)
# # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
# plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# plot(Mata_bbox, lwd = 1, border = "grey20", add = T)
# plot(Mata_borders, lwd = 2, border = "black", add = T)
# addRasterLegend(Mata.HII_contrasted, locs = seq(0, 30, 5), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-33.22, -32.2, -27, -13))
# graphics::text(x = -32.8, y = -11, font = 2, cex = 1.1, label = "HII")
# legend(x = "bottomright", legend = c("Mata Atlantica"), lty = c(1) , col = "black", lwd = 2, bty = "n", cex = 1.1, text.font = 2)
# scalebar(d = 1000, type = "line", lwd = 4, divs = 4, xy = c(-59, -7), label = c("", "1000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15))
# par(mar = internal_margins)
# dev.off()


### 1/ Species richness ####

# Load indice map
tot.sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))

# Normalized indice on global scale between [0-1]
min <- as.numeric(names(table(tot.sp.richness_Jaccard.80[])[2])) # Do not include communities with null values in the scaling
max <- tot.sp.richness_Jaccard.80@data@max
sp.richness_norm <- (tot.sp.richness_Jaccard.80 - min) / (max - min)
sp.richness_norm[sp.richness_norm[] < 0.002] <- 0.002 # Add minimal value for communities were the index was null
sp.richness_norm[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

save(sp.richness_norm, file = "./outputs/Threat_maps/sp.richness_norm.RData", version = "2")
saveRDS(sp.richness_norm, file = "./outputs/Threat_maps/sp.richness_norm.rds", version = "2")

### 1.1/ Compute threat map weighted by species richness ####
sp.richness_threat <- sp.richness_norm * HII_norm
sp.richness_threat_raw <- tot.sp.richness_Jaccard.80 * HII_norm

# Save
save(sp.richness_threat, file = "./outputs/Threat_maps/sp.richness_threat.RData", version = "2")
saveRDS(sp.richness_threat, file = "./outputs/Threat_maps/sp.richness_threat.rds", version = "2")
save(sp.richness_threat_raw, file = "./outputs/Threat_maps/Weighted_indices/sp.richness_threat_raw.RData", version = "2")
saveRDS(sp.richness_threat_raw, file = "./outputs/Threat_maps/Weighted_indices/sp.richness_threat_raw.rds", version = "2")


### 1.2/ Plot global map ####

### Weighted HII map

# Load threat index
sp.richness_threat <- readRDS(file = "./outputs/Threat_maps/sp.richness_threat.rds")

hist(sp.richness_threat)

# Contrasted version
sp.richness_threat_contrasted <- sp.richness_threat
sp.richness_threat_contrasted[sp.richness_threat[] > 0.15] <- 0.15
sp.richness_threat_contrasted[sp.richness_threat[] < 0.001] <- 0.001 # Avoid null values to be colored as ones where no Ithomiini are present
sp.richness_threat_contrasted[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

save(sp.richness_threat_contrasted, file = "./outputs/Threat_maps/sp.richness_threat_contrasted.RData", version = "2")
saveRDS(sp.richness_threat_contrasted, file = "./outputs/Threat_maps/sp.richness_threat_contrasted.rds", version = "2")


pdf(file = paste0("./maps/Threat_maps/sp.richness_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.richness_threat_contrasted, col = pal_grn_red_NA, main = paste0("Species richness\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.richness_threat_contrasted, locs = seq(0, 0.15, 0.03), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")
par(mar = internal_margins)
dev.off()

### Weighted indice map
# Load threat index
sp.richness_threat_raw <- readRDS(file = "./outputs/Threat_maps/sp.richness_threat_raw.rds")

hist(sp.richness_threat_raw)

# Contrasted version
sp.richness_threat_raw_contrasted <- sp.richness_threat_raw
sp.richness_threat_raw_contrasted[sp.richness_threat_raw[] > 20] <- 20
sp.richness_threat_raw_contrasted[sp.richness_threat_raw[] < 0.1] <- 0.1 # Avoid null values to be colored as ones where no Ithomiini are present
sp.richness_threat_raw_contrasted[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

save(sp.richness_threat_raw_contrasted, file = "./outputs/Threat_maps/Weighted_indices/sp.richness_threat_raw_contrasted.RData", version = "2")
saveRDS(sp.richness_threat_raw_contrasted, file = "./outputs/Threat_maps/Weighted_indices/sp.richness_threat_contrasted.rds", version = "2")


pdf(file = paste0("./maps/Threat_maps/Weighted_indices/sp.richness_threat_raw.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.richness_threat_raw_contrasted, col = pal_grn_red_NA, main = paste0("Species richness\nThreat level"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.richness_threat_raw_contrasted, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 5, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()


### 1.3/ Plot zoom on Andes ####

# Load and crop Indices maps
sp.richness_threat <- readRDS(file = "./outputs/Threat_maps/sp.richness_threat.rds")
Andes.sp.richness_threat <- crop(sp.richness_threat, Andes_ext)

# Contrasted version
hist(Andes.sp.richness_threat)

Andes.sp.richness_threat_contrasted <- Andes.sp.richness_threat
Andes.sp.richness_threat_contrasted[Andes.sp.richness_threat[] > 0.2] <- 0.2

save(Andes.sp.richness_threa_contrasted, file = "./ouputs/Threat_maps/Andes/Andes.sp.richness_threat_contrasted.RData", version = "2")
saveRDS(Andes.sp.richness_threat_contrasted, file = "./outputs/Threat_maps/Andes/Andes.sp.richness_threat_contrasted.rds", version = "2")

# Export plot
pdf(file = paste0("./maps/Threat_maps/Andes/Andes.sp.richness_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Andes.sp.richness_threat_contrasted, col = pal_grn_red, main = paste0("Species richness\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Andes_bbox, lwd = 1, border = "grey20", add = T)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
        drawlabels = F, col = "black", add = T)
addRasterLegend(Andes.sp.richness_threat_contrasted, locs = seq(0, 0.2, 0.05), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.75, 4.25))
graphics::text(x = -86.7, y = 6, font = 2, cex = 0.9, label = "Threat index")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 0.9, text.font = 2)
graphics::text(x = -86.7, y = -10.8, font = 2, cex = 0.9, label = "Elevation")
scalebar(d = 750, type = "line", lwd = 4, divs = 4, xy = c(-67, -13.5), label = c("", "750 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.6, padin = c(0.2, 0.4), pos = "topleft")
par(mar = internal_margins)
dev.off()


### 1.4/ Plot zoom on Mata Atlantica ####

# Load and crop Indices maps
sp.richness_threat <- readRDS(file = "./outputs/Threat_maps/sp.richness_threat.rds")
Mata.sp.richness_threat <- crop(sp.richness_threat, Mata_ext)

# Contrasted version
hist(Mata.sp.richness_threat)

Mata.sp.richness_threat_contrasted <- Mata.sp.richness_threat
Mata.sp.richness_threat_contrasted[Mata.sp.richness_threat[] > 0.1] <- 0.1

save(Mata.sp.richness_threa_contrasted, file = "./ouputs/Threat_maps/Mata/Mata.sp.richness_threat_contrasted.RData", version = "2")
saveRDS(Mata.sp.richness_threat_contrasted, file = "./outputs/Threat_maps/Mata/Mata.sp.richness_threat_contrasted.rds", version = "2")

# Export plot
pdf(file = paste0("./maps/Threat_maps/Mata/Mata.sp.richness_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Mata.sp.richness_threat_contrasted, col = pal_grn_red, main = paste0("Species richness\nThreat index"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2),
      legend  = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Mata_bbox, lwd = 1, border = "grey20", add = T)
plot(Mata_borders, lwd = 2, border = "black", add = T)
addRasterLegend(Mata.sp.richness_threat_contrasted, locs = seq(0, 0.1, 0.02), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-33.22, -32.2, -27, -13))
graphics::text(x = -32.8, y = -10.5, font = 2, cex = 1.1, label = "Threat\nindex")
legend(x = "bottomright", legend = c("Mata Atlantica"), lty = c(1) , col = "black", lwd = 2, bty = "n", cex = 1.1, text.font = 2)
scalebar(d = 1000, type = "line", lwd = 4, divs = 4, xy = c(-59, -7), label = c("", "1000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15))
par(mar = internal_margins)
dev.off()


### 2/ Species diversity ####

# Load indice map
sp.diversity.compatible_Jost_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jost_Jaccard.80.rds"))

# Normalized indice on global scale between [0-1]
min <- as.numeric(names(table(sp.diversity.compatible_Jost_Jaccard.80[])[2])) # Do not include communities with null values in the scaling
max <- sp.diversity.compatible_Jost_Jaccard.80@data@max
sp.diversity_norm <- (sp.diversity.compatible_Jost_Jaccard.80 - min) / (max - min)
sp.diversity_norm[sp.diversity_norm[] < 0.002] <- 0.002 # Add minimal value for communities were the index was null
sp.diversity_norm[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

save(sp.diversity_norm, file = "./outputs/Threat_maps/sp.diversity_norm.RData", version = "2")
saveRDS(sp.diversity_norm, file = "./outputs/Threat_maps/sp.diversity_norm.rds", version = "2")

### 2.1/ Compute threat map weighted by species richness ####
sp.diversity_threat <- sp.diversity_norm * HII_norm
sp.diversity_threat_raw <- sp.diversity.compatible_Jost_Jaccard.80 * HII_norm

# Save
save(sp.diversity_threat, file = "./outputs/Threat_maps/sp.diversity_threat.RData", version = "2")
saveRDS(sp.diversity_threat, file = "./outputs/Threat_maps/sp.diversity_threat.rds", version = "2")
save(sp.diversity_threat_raw, file = "./outputs/Threat_maps/Weighted_indices/sp.diversity_threat_raw.RData", version = "2")
saveRDS(sp.diversity_threat_raw, file = "./outputs/Threat_maps/Weighted_indices/sp.diversity_threat_raw.rds", version = "2")



### 2.2/ Plot global map ####

### Weighted HII

# Load threat index
sp.diversity_threat <- readRDS(file = "./outputs/Threat_maps/sp.diversity_threat.rds")

hist(sp.diversity_threat)

# Contrasted version
sp.diversity_threat_contrasted <- sp.diversity_threat
sp.diversity_threat_contrasted[sp.diversity_threat[] > 0.15] <- 0.15
sp.diversity_threat_contrasted[sp.diversity_threat[] < 0.001] <- 0.001 # Avoid null values to be colored as ones where no Ithomiini are present
sp.diversity_threat_contrasted[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

save(sp.diversity_threat_contrasted, file = "./outputs/Threat_maps/sp.diversity_threat_contrasted.RData", version = "2")
saveRDS(sp.diversity_threat_contrasted, file = "./outputs/Threat_maps/sp.diversity_threat_contrasted.rds", version = "2")


pdf(file = paste0("./maps/Threat_maps/sp.diversity_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity_threat_contrasted, col = pal_grn_red_NA, main = paste0("Species diversity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.diversity_threat_contrasted, locs = seq(0, 0.15, 0.03), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")
par(mar = internal_margins)
dev.off()

### Weighted index

# Load threat index
sp.diversity_threat_raw <- readRDS(file = "./outputs/Threat_maps/Weighted_indices/sp.diversity_threat_raw.rds")

hist(sp.diversity_threat_raw)

# Contrasted version
sp.diversity_threat_raw_contrasted <- sp.diversity_threat_raw
sp.diversity_threat_raw_contrasted[sp.diversity_threat_raw[] > 20] <- 20
sp.diversity_threat_raw_contrasted[sp.diversity_threat_raw[] < 0.1] <- 0.1 # Avoid null values to be colored as ones where no Ithomiini are present
sp.diversity_threat_raw_contrasted[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

save(sp.diversity_threat_raw_contrasted, file = "./outputs/Threat_maps/Weighted_indices/sp.diversity_threat_raw_contrasted.RData", version = "2")
saveRDS(sp.diversity_threat_raw_contrasted, file = "./outputs/Threat_maps/Weighted_indices/sp.diversity_threat_raw_contrasted.rds", version = "2")


pdf(file = paste0("./maps/Threat_maps/Weighted_indices/sp.diversity_threat_raw.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.diversity_threat_raw_contrasted, col = pal_grn_red_NA, main = paste0("Species diversity\nThreat level"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.diversity_threat_raw_contrasted, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 5.5, font = 2, cex = 1.2, label = "Species")
par(mar = internal_margins)
dev.off()


### 2.3/ Plot zoom on Andes ####

# Load and crop Indices maps
sp.diversity_threat <- readRDS(file = "./outputs/Threat_maps/sp.diversity_threat.rds")
Andes.sp.diversity_threat <- crop(sp.diversity_threat, Andes_ext)

# Contrasted version
hist(Andes.sp.diversity_threat)

Andes.sp.diversity_threat_contrasted <- Andes.sp.diversity_threat
Andes.sp.diversity_threat_contrasted[Andes.sp.diversity_threat[] > 0.2] <- 0.2

save(Andes.sp.diversity_threat_contrasted, file = "./outputs/Threat_maps/Andes/Andes.sp.diversity_threat_contrasted.RData", version = "2")
saveRDS(Andes.sp.diversity_threat_contrasted, file = "./outputs/Threat_maps/Andes/Andes.sp.diversity_threat_contrasted.rds", version = "2")

# Export plot
pdf(file = paste0("./maps/Threat_maps/Andes/Andes.sp.diversity_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Andes.sp.diversity_threat_contrasted, col = pal_grn_red, main = paste0("Species diversity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Andes_bbox, lwd = 1, border = "grey20", add = T)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
        drawlabels = F, col = "black", add = T)
addRasterLegend(Andes.sp.diversity_threat_contrasted, locs = seq(0, 0.2, 0.05), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.75, 4.25))
graphics::text(x = -86.7, y = 6, font = 2, cex = 0.9, label = "Threat index")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 0.9, text.font = 2)
graphics::text(x = -86.7, y = -10.8, font = 2, cex = 0.9, label = "Elevation")
scalebar(d = 750, type = "line", lwd = 4, divs = 4, xy = c(-67, -13.5), label = c("", "750 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.6, padin = c(0.2, 0.4), pos = "topleft")
par(mar = internal_margins)
dev.off()


### 2.4/ Plot zoom on Mata Atlantica ####

# Load and crop Indices maps
sp.diversity_threat <- readRDS(file = "./outputs/Threat_maps/sp.diversity_threat.rds")
Mata.sp.diversity_threat <- crop(sp.diversity_threat, Mata_ext)

# Contrasted version
hist(Mata.sp.diversity_threat)

Mata.sp.diversity_threat_contrasted <- Mata.sp.diversity_threat
Mata.sp.diversity_threat_contrasted[Mata.sp.diversity_threat[] > 0.1] <- 0.1

save(Mata.sp.diversity_threat_contrasted, file = "./outputs/Threat_maps/Mata/Mata.sp.diversity_threat_contrasted.RData", version = "2")
saveRDS(Mata.sp.diversity_threat_contrasted, file = "./outputs/Threat_maps/Mata/Mata.sp.diversity_threat_contrasted.rds", version = "2")

# Export plot
pdf(file = paste0("./maps/Threat_maps/Mata/Mata.sp.diversity_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Mata.sp.diversity_threat_contrasted, col = pal_grn_red, main = paste0("Species diversity\nThreat index"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2),
      legend  = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Mata_bbox, lwd = 1, border = "grey20", add = T)
plot(Mata_borders, lwd = 2, border = "black", add = T)
addRasterLegend(Mata.sp.diversity_threat_contrasted, locs = seq(0, 0.1, 0.02), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-33.22, -32.2, -27, -13))
graphics::text(x = -32.8, y = -10.5, font = 2, cex = 1.1, label = "Threat\nindex")
legend(x = "bottomright", legend = c("Mata Atlantica"), lty = c(1) , col = "black", lwd = 2, bty = "n", cex = 1.1, text.font = 2)
scalebar(d = 1000, type = "line", lwd = 4, divs = 4, xy = c(-59, -7), label = c("", "1000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15))
par(mar = internal_margins)
dev.off()


### 3/ Species rarity ####

# Load indice map
sp.mean.rarity_Jaccard.80 <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_Jaccard.80.rds")

# Normalized indice on global scale between [0-1]
min <- as.numeric(names(table(sp.mean.rarity_Jaccard.80[])[2])) # Do not include communities with null values in the scaling
max <- sp.mean.rarity_Jaccard.80@data@max
sp.rarity_norm <- (sp.mean.rarity_Jaccard.80 - min) / (max - min)
sp.rarity_norm[sp.rarity_norm[] < 0.002] <- 0.002 # Add minimal value for communities were the index was null
sp.rarity_norm[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

save(sp.rarity_norm, file = "./outputs/Threat_maps/sp.rarity_norm.RData", version = "2")
saveRDS(sp.rarity_norm, file = "./outputs/Threat_maps/sp.rarity_norm.rds", version = "2")


### 3.1/ Compute threat map weighted by species richness ####
sp.rarity_threat <- sp.rarity_norm * HII_norm
sp.rarity_threat_raw <- sp.mean.rarity_Jaccard.80 * HII_norm

# Save
save(sp.rarity_threat, file = "./outputs/Threat_maps/sp.rarity_threat.RData", version = "2")
saveRDS(sp.rarity_threat, file = "./outputs/Threat_maps/sp.rarity_threat.rds", version = "2")
save(sp.rarity_threat_raw, file = "./outputs/Threat_maps/Weighted_indices/sp.rarity_threat_raw.RData", version = "2")
saveRDS(sp.rarity_threat_raw, file = "./outputs/Threat_maps/Weighted_indices/sp.rarity_threat_raw.rds", version = "2")



### 3.2/ Plot global map ####

### Weighted HII

# Load threat index
sp.rarity_threat <- readRDS(file = "./outputs/Threat_maps/sp.rarity_threat.rds")

hist(sp.rarity_threat)

# Contrasted version
sp.rarity_threat_contrasted <- sp.rarity_threat
sp.rarity_threat_contrasted[sp.rarity_threat[] > 0.3] <- 0.3
sp.rarity_threat_contrasted[sp.rarity_threat[] < 0.002] <- 0.002 # Avoid null values to be colored as ones where no Ithomiini are present
sp.rarity_threat_contrasted[sp.mean.rarity_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

save(sp.rarity_threat_contrasted, file = "./outputs/Threat_maps/sp.rarity_threat_contrasted.RData", version = "2")
saveRDS(sp.rarity_threat_contrasted, file = "./outputs/Threat_maps/sp.rarity_threat_contrasted.rds", version = "2")


pdf(file = paste0("./maps/Threat_maps/sp.rarity_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.rarity_threat_contrasted, col = pal_grn_red_NA, main = paste0("Mean species rarity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.rarity_threat_contrasted, locs = seq(0, 0.3, 0.05), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")
par(mar = internal_margins)
dev.off()


### Weighted Index

# Load threat index
sp.rarity_threat_raw <- readRDS(file = "./outputs/Threat_maps/Weighted_indices/sp.rarity_threat_raw.rds")

hist(sp.rarity_threat_raw)

# Contrasted version
sp.rarity_threat_raw_contrasted <- sp.rarity_threat_raw
sp.rarity_threat_raw_contrasted[sp.rarity_threat_raw[] > 0.35] <- 0.35
sp.rarity_threat_raw_contrasted[sp.rarity_threat_raw[] < 0.002] <- 0.002 # Avoid null values to be colored as ones where no Ithomiini are present
sp.rarity_threat_raw_contrasted[sp.mean.rarity_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

save(sp.rarity_threat_raw_contrasted, file = "./outputs/Threat_maps/Weighted_indices/sp.rarity_threat_raw_contrasted.RData", version = "2")
saveRDS(sp.rarity_threat_raw_contrasted, file = "./outputs/Threat_maps/Weighted_indices/sp.rarity_threat_raw_contrasted.rds", version = "2")


pdf(file = paste0("./maps/Threat_maps/Weighted_indices/sp.rarity_threat_raw.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(sp.rarity_threat_raw_contrasted, col = pal_grn_red_NA, main = paste0("Mean species rarity\nThreat level"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Rarity\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.rarity_threat_raw_contrasted, locs = seq(0, 0.35, 0.1), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Rarity\nindex")
par(mar = internal_margins)
dev.off()


### 3.3/ Plot zoom on Andes ####

# Load and crop Indices maps
sp.rarity_threat <- readRDS(file = "./outputs/Threat_maps/sp.rarity_threat.rds")
Andes.sp.rarity_threat <- crop(sp.rarity_threat, Andes_ext)

# Contrasted version
hist(Andes.sp.rarity_threat)

Andes.sp.rarity_threat_contrasted <- Andes.sp.rarity_threat
Andes.sp.rarity_threat_contrasted[Andes.sp.rarity_threat[] > 0.30] <- 0.30

save(Andes.sp.rarity_threat_contrasted, file = "./outputs/Threat_maps/Andes/Andes.sp.rarity_threat_contrasted.RData", version = "2")
saveRDS(Andes.sp.rarity_threat_contrasted, file = "./outputs/Threat_maps/Andes/Andes.sp.rarity_threat_contrasted.rds", version = "2")

# Export plot
pdf(file = paste0("./maps/Threat_maps/Andes/Andes.sp.rarity_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Andes.sp.rarity_threat_contrasted, col = pal_grn_red, main = paste0("Mean species rarity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Andes_bbox, lwd = 1, border = "grey20", add = T)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
        drawlabels = F, col = "black", add = T)
addRasterLegend(Andes.sp.rarity_threat_contrasted, locs = seq(0, 0.3, 0.05), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.75, 4.25))
graphics::text(x = -86.7, y = 6, font = 2, cex = 0.9, label = "Threat index")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 0.9, text.font = 2)
graphics::text(x = -86.7, y = -10.8, font = 2, cex = 0.9, label = "Elevation")
scalebar(d = 750, type = "line", lwd = 4, divs = 4, xy = c(-67, -13.5), label = c("", "750 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.6, padin = c(0.2, 0.4), pos = "topleft")
par(mar = internal_margins)
dev.off()


### 3.4/ Plot zoom on Mata Atlantica ####

# Load and crop Indices maps
sp.rarity_threat <- readRDS(file = "./outputs/Threat_maps/sp.rarity_threat.rds")
Mata.sp.rarity_threat <- crop(sp.rarity_threat, Mata_ext)

# Contrasted version
hist(Mata.sp.rarity_threat)

Mata.sp.rarity_threat_contrasted <- Mata.sp.rarity_threat
Mata.sp.rarity_threat_contrasted[Mata.sp.rarity_threat[] > 0.3] <- 0.3

save(Mata.sp.rarity_threat_contrasted, file = "./outputs/Threat_maps/Mata/Mata.sp.rarity_threat_contrasted.RData", version = "2")
saveRDS(Mata.sp.rarity_threat_contrasted, file = "./outputs/Threat_maps/Mata/Mata.sp.rarity_threat_contrasted.rds", version = "2")

# Export plot
pdf(file = paste0("./maps/Threat_maps/Mata/Mata.sp.rarity_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Mata.sp.rarity_threat_contrasted, col = pal_grn_red, main = paste0("Mean species rarity\nThreat index"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2),
      legend  = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Mata_bbox, lwd = 1, border = "grey20", add = T)
plot(Mata_borders, lwd = 2, border = "black", add = T)
addRasterLegend(Mata.sp.rarity_threat_contrasted, locs = seq(0, 0.3, 0.05), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-33.22, -32.2, -27, -13))
graphics::text(x = -32.8, y = -10.5, font = 2, cex = 1.1, label = "Threat\nindex")
legend(x = "bottomright", legend = c("Mata Atlantica"), lty = c(1) , col = "black", lwd = 2, bty = "n", cex = 1.1, text.font = 2)
scalebar(d = 1000, type = "line", lwd = 4, divs = 4, xy = c(-59, -7), label = c("", "1000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15))
par(mar = internal_margins)
dev.off()



### 4/ Faith's Phylogenetic Diversity ####

# Load indice map
PD.raster_Jaccard.80 <- readRDS(file = "./outputs/Indices_Maps/PD.raster_Jaccard.80.rds")

# Normalized indice on global scale
min <- as.numeric(names(table(PD.raster_Jaccard.80[])[2])) # Do not include communities with null values in the scaling
max <- PD.raster_Jaccard.80@data@max
PD_norm <- (PD.raster_Jaccard.80 - min) / (max - min)
PD_norm[PD_norm[] < 0.002] <- 0.002 # Add minimal value for communities were the index was null
PD_norm[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

# plot(PD_norm, col = pal_grn_red_NA)

### 4.1/ Compute threat map weighted by species richness ####
PD_threat <- PD_norm * HII_norm
PD_threat_raw <- PD.raster_Jaccard.80 * HII_norm

# Save
save(PD_threat, file = "./outputs/Threat_maps/PD_threat.RData", version = "2")
saveRDS(PD_threat, file = "./outputs/Threat_maps/PD_threat.rds", version = "2")
save(PD_threat_raw, file = "./outputs/Threat_maps/Weighted_indices/PD_threat_raw.RData", version = "2")
saveRDS(PD_threat_raw, file = "./outputs/Threat_maps/Weighted_indices/PD_threat_raw.rds", version = "2")


### 4.2/ Plot global map ####

### Weighted HII

# Load threat index
PD_threat <- readRDS(file = "./outputs/Threat_maps/PD_threat.rds")

hist(PD_threat)

# Contrasted version
PD_threat_contrasted <- PD_threat
PD_threat_contrasted[PD_threat[] > 0.2] <- 0.2
PD_threat_contrasted[PD_threat[] < 0.001] <- 0.001 # Avoid null values to be colored as ones where no Ithomiini are present
PD_threat_contrasted[sp.mean.rarity_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

save(PD_threat_contrasted, file = "./outputs/Threat_maps/PD_threat_contrasted.RData", version = "2")
saveRDS(PD_threat_contrasted, file = "./outputs/Threat_maps/PD_threat_contrasted.rds", version = "2")


pdf(file = paste0("./maps/Threat_maps/PD_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(PD_threat_contrasted, col = pal_grn_red_NA, main = paste0("Faith's Phylogenetic Diversity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(PD_threat_contrasted, locs = seq(0, 0.2, 0.05), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")
par(mar = internal_margins)
dev.off()

### Weighted Index

# Load threat index
PD_threat_raw <- readRDS(file = "./outputs/Threat_maps/Weighted_indices/PD_threat_raw.rds")

hist(PD_threat_raw)

# Contrasted version
PD_threat_raw_contrasted <- PD_threat_raw
PD_threat_raw_contrasted[PD_threat_raw[] > 180] <- 180
PD_threat_raw_contrasted[PD_threat_raw[] < 1] <- 1 # Avoid null values to be colored as ones where no Ithomiini are present
PD_threat_raw_contrasted[sp.mean.rarity_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

save(PD_threat_raw_contrasted, file = "./outputs/Threat_maps/Weighted_indices/PD_threat_raw_contrasted.RData", version = "2")
saveRDS(PD_threat_raw_contrasted, file = "./outputs/Threat_maps/Weighted_indices/PD_threat_raw_contrasted.rds", version = "2")


pdf(file = paste0("./maps/Threat_maps/Weighted_indices/PD_threat_raw.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(PD_threat_raw_contrasted, col = pal_grn_red_NA, main = paste0("Faith's Phylogenetic Diversity\nThreat level"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(PD_threat_raw_contrasted, locs = seq(0, 180, 40), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -109, y = 6.5, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()


### 4.3/ Plot zoom on Andes ####

# Load and crop Indices maps
PD_threat <- readRDS(file = "./outputs/Threat_maps/PD_threat.rds")
Andes.PD_threat <- crop(PD_threat, Andes_ext)

# Contrasted version
hist(Andes.PD_threat)

Andes.PD_threat_contrasted <- Andes.PD_threat
Andes.PD_threat_contrasted[Andes.PD_threat[] > 0.20] <- 0.20

save(Andes.PD_threat_contrasted, file = "./outputs/Threat_maps/Andes/Andes.PD_threat_contrasted.RData", version = "2")
saveRDS(Andes.PD_threat_contrasted, file = "./outputs/Threat_maps/Andes/Andes.PD_threat_contrasted.rds", version = "2")

# Export plot
pdf(file = paste0("./maps/Threat_maps/Andes/Andes.PD_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Andes.PD_threat_contrasted, col = pal_grn_red, main = paste0("Faith's Phylogenetic Diversity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Andes_bbox, lwd = 1, border = "grey20", add = T)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
        drawlabels = F, col = "black", add = T)
addRasterLegend(Andes.PD_threat_contrasted, locs = seq(0, 0.2, 0.05), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.75, 4.25))
graphics::text(x = -86.7, y = 6, font = 2, cex = 0.9, label = "Threat index")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 0.9, text.font = 2)
graphics::text(x = -86.7, y = -10.8, font = 2, cex = 0.9, label = "Elevation")
scalebar(d = 750, type = "line", lwd = 4, divs = 4, xy = c(-67, -13.5), label = c("", "750 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.6, padin = c(0.2, 0.4), pos = "topleft")
par(mar = internal_margins)
dev.off()


### 4.4/ Plot zoom on Mata Atlantica ####

# Load and crop Indices maps
PD_threat <- readRDS(file = "./outputs/Threat_maps/PD_threat.rds")
Mata.PD_threat <- crop(PD_threat, Mata_ext)

# Contrasted version
hist(Mata.PD_threat)

Mata.PD_threat_contrasted <- Mata.PD_threat
Mata.PD_threat_contrasted[Mata.PD_threat[] > 0.15] <- 0.15

save(Mata.PD_threat_contrasted, file = "./outputs/Threat_maps/Mata/Mata.PD_threat_contrasted.RData", version = "2")
saveRDS(Mata.PD_threat_contrasted, file = "./outputs/Threat_maps/Mata/Mata.PD_threat_contrasted.rds", version = "2")

# Export plot
pdf(file = paste0("./maps/Threat_maps/Mata/Mata.PD_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Mata.PD_threat_contrasted, col = pal_grn_red, main = paste0("Faith's Phylogenetic Diversity\nThreat index"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2),
      legend  = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Mata_bbox, lwd = 1, border = "grey20", add = T)
plot(Mata_borders, lwd = 2, border = "black", add = T)
addRasterLegend(Mata.PD_threat_contrasted, locs = seq(0, 0.15, 0.03), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-33.22, -32.2, -27, -13))
graphics::text(x = -32.8, y = -10.5, font = 2, cex = 1.1, label = "Threat\nindex")
legend(x = "bottomright", legend = c("Mata Atlantica"), lty = c(1) , col = "black", lwd = 2, bty = "n", cex = 1.1, text.font = 2)
scalebar(d = 1000, type = "line", lwd = 4, divs = 4, xy = c(-59, -7), label = c("", "1000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15))
par(mar = internal_margins)
dev.off()



### 5/ Sum of Fair-proportions ####

# Load indice map
FP.raster_Jaccard.80 <- readRDS(file = "./outputs/Indices_Maps/FP.raster_Jaccard.80.rds")

# Normalized indice on global scale
min <- as.numeric(names(table(FP.raster_Jaccard.80[])[2])) # Do not include communities with null values in the scaling
max <- FP.raster_Jaccard.80@data@max
FP_norm <- (FP.raster_Jaccard.80 - min) / (max - min)
FP_norm[FP_norm[] < 0.002] <- 0.002 # Add minimal value for communities were the index was null
FP_norm[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

# plot(FP_norm, col = pal_grn_red_NA)

### 5.1/ Compute threat map weighted by species richness ####
FP_threat <- FP_norm * HII_norm
FP_threat_raw <- FP.raster_Jaccard.80 * HII_norm

# Save
save(FP_threat, file = "./outputs/Threat_maps/FP_threat.RData", version = "2")
saveRDS(FP_threat, file = "./outputs/Threat_maps/FP_threat.rds", version = "2")
save(FP_threat_raw, file = "./outputs/Threat_maps/Weighted_indices/FP_threat_raw.RData", version = "2")
saveRDS(FP_threat_raw, file = "./outputs/Threat_maps/Weighted_indices/FP_threat_raw.rds", version = "2")

### 5.2/ Plot global map ####

### Weighted HII

# Load threat index
FP_threat <- readRDS(file = "./outputs/Threat_maps/FP_threat.rds")

hist(FP_threat)

# Contrasted version
FP_threat_contrasted <- FP_threat
FP_threat_contrasted[FP_threat[] > 0.15] <- 0.15
FP_threat_contrasted[FP_threat[] < 0.001] <- 0.001 # Avoid null values to be colored as ones where no Ithomiini are present
FP_threat_contrasted[sp.mean.rarity_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

save(FP_threat_contrasted, file = "./outputs/Threat_maps/FP_threat_contrasted.RData", version = "2")
saveRDS(FP_threat_contrasted, file = "./outputs/Threat_maps/FP_threat_contrasted.rds", version = "2")


pdf(file = paste0("./maps/Threat_maps/FP_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(FP_threat_contrasted, col = pal_grn_red_NA, main = paste0("Sum of Fair-proportions\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(FP_threat_contrasted, locs = seq(0, 0.15, 0.03), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")
par(mar = internal_margins)
dev.off()

### Weighted Index

# Load threat index
FP_threat_raw <- readRDS(file = "./outputs/Threat_maps/Weighted_indices/FP_threat_raw.rds")

hist(FP_threat_raw)

# Contrasted version
FP_threat_raw_contrasted <- FP_threat_raw
FP_threat_raw_contrasted[FP_threat_raw[] > 75] <- 75
FP_threat_raw_contrasted[FP_threat_raw[] < 0.5] <- 0.5 # Avoid null values to be colored as ones where no Ithomiini are present
FP_threat_raw_contrasted[sp.mean.rarity_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

save(FP_threat_raw_contrasted, file = "./outputs/Threat_maps/Weighted_indices/FP_threat_raw_contrasted.RData", version = "2")
saveRDS(FP_threat_raw_contrasted, file = "./outputs/Threat_maps/Weighted_indices/FP_threat_raw_contrasted.rds", version = "2")


pdf(file = paste0("./maps/Threat_maps/Weighted_indices/FP_threat_raw.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(FP_threat_raw_contrasted, col = pal_grn_red_NA, main = paste0("Sum of Fair-proportions\nThreat level"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(FP_threat_raw_contrasted, locs = seq(0, 75, 20), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -109, y = 6.5, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()


### 5.3/ Plot zoom on Andes ####

# Load and crop Indices maps
FP_threat <- readRDS(file = "./outputs/Threat_maps/FP_threat.rds")
Andes.FP_threat <- crop(FP_threat, Andes_ext)

# Contrasted version
hist(Andes.FP_threat)

Andes.FP_threat_contrasted <- Andes.FP_threat
Andes.FP_threat_contrasted[Andes.FP_threat[] > 0.15] <- 0.15

save(Andes.FP_threat_contrasted, file = "./outputs/Threat_maps/Andes/Andes.FP_threat_contrasted.RData", version = "2")
saveRDS(Andes.FP_threat_contrasted, file = "./outputs/Threat_maps/Andes/Andes.FP_threat_contrasted.rds", version = "2")

# Export plot
pdf(file = paste0("./maps/Threat_maps/Andes/Andes.FP_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Andes.FP_threat_contrasted, col = pal_grn_red, main = paste0("Sum of Fair-proportions\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Andes_bbox, lwd = 1, border = "grey20", add = T)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
        drawlabels = F, col = "black", add = T)
addRasterLegend(Andes.FP_threat_contrasted, locs = seq(0, 0.15, 0.03), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.75, 4.25))
graphics::text(x = -86.7, y = 6, font = 2, cex = 0.9, label = "Threat index")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 0.9, text.font = 2)
graphics::text(x = -86.7, y = -10.8, font = 2, cex = 0.9, label = "Elevation")
scalebar(d = 750, type = "line", lwd = 4, divs = 4, xy = c(-67, -13.5), label = c("", "750 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.6, padin = c(0.2, 0.4), pos = "topleft")
par(mar = internal_margins)
dev.off()


### 5.4/ Plot zoom on Mata Atlantica ####

# Load and crop Indices maps
FP_threat <- readRDS(file = "./outputs/Threat_maps/FP_threat.rds")
Mata.FP_threat <- crop(FP_threat, Mata_ext)

# Contrasted version
hist(Mata.FP_threat)

Mata.FP_threat_contrasted <- Mata.FP_threat
Mata.FP_threat_contrasted[Mata.FP_threat[] > 0.10] <- 0.10

save(Mata.FP_threat_contrasted, file = "./outputs/Threat_maps/Mata/Mata.FP_threat_contrasted.RData", version = "2")
saveRDS(Mata.FP_threat_contrasted, file = "./outputs/Threat_maps/Mata/Mata.FP_threat_contrasted.rds", version = "2")

# Export plot
pdf(file = paste0("./maps/Threat_maps/Mata/Mata.FP_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Mata.FP_threat_contrasted, col = pal_grn_red, main = paste0("Sum of Fair-proportions\nThreat index"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2),
      legend  = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Mata_bbox, lwd = 1, border = "grey20", add = T)
plot(Mata_borders, lwd = 2, border = "black", add = T)
addRasterLegend(Mata.FP_threat_contrasted, locs = seq(0, 0.10, 0.02), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-33.22, -32.2, -27, -13))
graphics::text(x = -32.8, y = -10.5, font = 2, cex = 1.1, label = "Threat\nindex")
legend(x = "bottomright", legend = c("Mata Atlantica"), lty = c(1) , col = "black", lwd = 2, bty = "n", cex = 1.1, text.font = 2)
scalebar(d = 1000, type = "line", lwd = 4, divs = 4, xy = c(-59, -7), label = c("", "1000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15))
par(mar = internal_margins)
dev.off()


### 6/ Mean pairwise Phylogenetic Distance ####

# Load indice map
MPD.raster_Jaccard.80 <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80.rds")

# Normalized indice on global scale
min <- as.numeric(names(table(MPD.raster_Jaccard.80[])[2])) # Do not include communities with null values in the scaling
max <- max(MPD.raster_Jaccard.80[], na.rm = T)
MPD_norm <- (MPD.raster_Jaccard.80 - min) / (max - min)
MPD_norm[MPD_norm[] < 0.002] <- 0.002 # Add minimal value for communities were the index was null
MPD_norm[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

# plot(MPD_norm, col = pal_grn_red_NA, main = "Normalized MPD")

save(MPD_norm, file = "./outputs/Threat_maps/MPD_norm.RData", version = "2")
saveRDS(MPD_norm, file = "./outputs/Threat_maps/MPD_norm.rds", version = "2")

### 6.1/ Compute threat map weighted by species richness ####
MPD_threat <- MPD_norm * HII_norm
MPD_threat <- merge(MPD_threat, continent_mask) # 0 values for community with no Ithomiini
MPD_threat[(tot.sp.richness_Jaccard.80[] > 0) & (MPD_threat[] < 0.002)] <- 0.002 # Avoid null values to be colored as ones where no Ithomiini are present

MPD_threat_raw <- MPD.raster_Jaccard.80 * HII_norm
MPD_threat_raw <- merge(MPD_threat_raw, continent_mask) # 0 values for community with no Ithomiini
MPD_threat_raw[(tot.sp.richness_Jaccard.80[] > 0) & (MPD_threat_raw[] < 0.5)] <- 0.5 # Avoid null values to be colored as ones where no Ithomiini are present

# plot(MPD_threat, col = pal_grn_red_NA)
# plot(MPD_threat_raw, col = pal_grn_red_NA)

# Save
save(MPD_threat, file = "./outputs/Threat_maps/MPD_threat.RData", version = "2")
saveRDS(MPD_threat, file = "./outputs/Threat_maps/MPD_threat.rds", version = "2")
save(MPD_threat_raw, file = "./outputs/Threat_maps/Weighted_indices/MPD_threat_raw.RData", version = "2")
saveRDS(MPD_threat_raw, file = "./outputs/Threat_maps/Weighted_indices/MPD_threat_raw.rds", version = "2")

### 6.2/ Plot global map ####

### Weighted HII

# Load threat index
MPD_threat <- readRDS(file = "./outputs/Threat_maps/MPD_threat.rds")

hist(MPD_threat)

# Contrasted version
MPD_threat_contrasted <- MPD_threat
MPD_threat_contrasted[MPD_threat[] > 0.4] <- 0.4

save(MPD_threat_contrasted, file = "./outputs/Threat_maps/MPD_threat_contrasted.RData", version = "2")
saveRDS(MPD_threat_contrasted, file = "./outputs/Threat_maps/MPD_threat_contrasted.rds", version = "2")


pdf(file = paste0("./maps/Threat_maps/MPD_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(MPD_threat_contrasted, col = pal_grn_red_NA, main = paste0("Mean pairwise Phylogenetic Distance\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(MPD_threat_contrasted, locs = seq(0, 0.4, 0.1), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")
par(mar = internal_margins)
dev.off()


### Weighted Index

# Load threat index
MPD_threat_raw <- readRDS(file = "./outputs/Threat_maps/Weighted_indices/MPD_threat_raw.rds")

hist(MPD_threat_raw)

# Contrasted version
MPD_threat_raw_contrasted <- MPD_threat_raw
MPD_threat_raw_contrasted[MPD_threat_raw[] > 20] <- 20

save(MPD_threat_raw_contrasted, file = "./outputs/Threat_maps/Weighted_indices/MPD_threat_raw_contrasted.RData", version = "2")
saveRDS(MPD_threat_raw_contrasted, file = "./outputs/Threat_maps/Weighted_indices/MPD_threat_raw_contrasted.rds", version = "2")


pdf(file = paste0("./maps/Threat_maps/Weighted_indices/MPD_threat_raw.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(MPD_threat_raw_contrasted, col = pal_grn_red_NA, main = paste0("Mean pairwise Phylogenetic Distance\nThreat level"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(MPD_threat_raw_contrasted, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -109, y = 6.5, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")
par(mar = internal_margins)
dev.off()


### 6.3/ Plot zoom on Andes ####

# Load and crop Indices maps
MPD_threat <- readRDS(file = "./outputs/Threat_maps/MPD_threat.rds")
Andes.MPD_threat <- crop(MPD_threat, Andes_ext)

# Contrasted version
hist(Andes.MPD_threat)

Andes.MPD_threat_contrasted <- Andes.MPD_threat
Andes.MPD_threat_contrasted[Andes.MPD_threat[] > 0.4] <- 0.4

save(Andes.MPD_threat_contrasted, file = "./outputs/Threat_maps/Andes/Andes.MPD_threat_contrasted.RData", version = "2")
saveRDS(Andes.MPD_threat_contrasted, file = "./outputs/Threat_maps/Andes/Andes.MPD_threat_contrasted.rds", version = "2")

# Export plot
pdf(file = paste0("./maps/Threat_maps/Andes/Andes.MPD_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Andes.MPD_threat_contrasted, col = pal_grn_red, main = paste0("Mean pairwise Phylogenetic Distance\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Andes_bbox, lwd = 1, border = "grey20", add = T)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
        drawlabels = F, col = "black", add = T)
addRasterLegend(Andes.MPD_threat_contrasted, locs = seq(0, 0.4, 0.1), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.75, 4.25))
graphics::text(x = -86.7, y = 6, font = 2, cex = 0.9, label = "Threat index")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 0.9, text.font = 2)
graphics::text(x = -86.7, y = -10.8, font = 2, cex = 0.9, label = "Elevation")
scalebar(d = 750, type = "line", lwd = 4, divs = 4, xy = c(-67, -13.5), label = c("", "750 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.6, padin = c(0.2, 0.4), pos = "topleft")
par(mar = internal_margins)
dev.off()


### 6.4/ Plot zoom on Mata Atlantica ####

# Load and crop Indices maps
MPD_threat <- readRDS(file = "./outputs/Threat_maps/MPD_threat.rds")
Mata.MPD_threat <- crop(MPD_threat, Mata_ext)

# Contrasted version
hist(Mata.MPD_threat)

Mata.MPD_threat_contrasted <- Mata.MPD_threat
Mata.MPD_threat_contrasted[Mata.MPD_threat[] > 0.4] <- 0.4

save(Mata.MPD_threat_contrasted, file = "./outputs/Threat_maps/Mata/Mata.MPD_threat_contrasted.RData", version = "2")
saveRDS(Mata.MPD_threat_contrasted, file = "./outputs/Threat_maps/Mata/Mata.MPD_threat_contrasted.rds", version = "2")

# Export plot
pdf(file = paste0("./maps/Threat_maps/Mata/Mata.MPD_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Mata.MPD_threat_contrasted, col = pal_grn_red, main = paste0("Mean pairwise Phylogenetic Distance\nThreat index"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2),
      legend  = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Mata_bbox, lwd = 1, border = "grey20", add = T)
plot(Mata_borders, lwd = 2, border = "black", add = T)
addRasterLegend(Mata.MPD_threat_contrasted, locs = seq(0, 0.4, 0.1), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-33.22, -32.2, -27, -13))
graphics::text(x = -32.8, y = -10.5, font = 2, cex = 1.1, label = "Threat\nindex")
legend(x = "bottomright", legend = c("Mata Atlantica"), lty = c(1) , col = "black", lwd = 2, bty = "n", cex = 1.1, text.font = 2)
scalebar(d = 1000, type = "line", lwd = 4, divs = 4, xy = c(-59, -7), label = c("", "1000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15))
par(mar = internal_margins)
dev.off()



### 7/ Mimicry richness ####

# Load indice map
ring.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))

# Normalized indice on global scale
min <- as.numeric(names(table(ring.richness_Jaccard.80[])[2])) # Do not include communities with null values in the scaling
max <- max(ring.richness_Jaccard.80[], na.rm = T)
ring.richness_norm <- (ring.richness_Jaccard.80 - min) / (max - min)
ring.richness_norm[ring.richness_norm[] < 0.002] <- 0.002 # Add minimal value for communities were the index was null
ring.richness_norm[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

# plot(ring.richness_norm, col = pal_grn_red_NA)

save(ring.richness_norm, file = "./outputs/Threat_maps/ring.richness_norm.RData", version = "2")
saveRDS(ring.richness_norm, file = "./outputs/Threat_maps/ring.richness_norm.rds", version = "2")


### 7.1/ Compute threat map weighted by species richness ####
ring.richness_threat <- ring.richness_norm * HII_norm
ring.richness_threat_raw <- ring.richness_Jaccard.80 * HII_norm

# Save
save(ring.richness_threat, file = "./outputs/Threat_maps/ring.richness_threat.RData", version = "2")
saveRDS(ring.richness_threat, file = "./outputs/Threat_maps/ring.richness_threat.rds", version = "2")
save(ring.richness_threat_raw, file = "./outputs/Threat_maps/Weighted_indices/ring.richness_threat_raw.RData", version = "2")
saveRDS(ring.richness_threat_raw, file = "./outputs/Threat_maps/Weighted_indices/ring.richness_threat_raw.rds", version = "2")


### 7.2/ Plot global map ####

### Weighted HII

# Load threat index
ring.richness_threat <- readRDS(file = "./outputs/Threat_maps/ring.richness_threat.rds")

hist(ring.richness_threat)

# Contrasted version
tot.sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
ring.richness_threat_contrasted <- ring.richness_threat
ring.richness_threat_contrasted[ring.richness_threat[] > 0.2] <- 0.2
ring.richness_threat_contrasted[ring.richness_threat[] < 0.001] <- 0.001 # Avoid null values to be colored as ones where no Ithomiini are present
ring.richness_threat_contrasted[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

save(ring.richness_threat_contrasted, file = "./outputs/Threat_maps/ring.richness_threat_contrasted.RData", version = "2")
saveRDS(ring.richness_threat_contrasted, file = "./outputs/Threat_maps/ring.richness_threat_contrasted.rds", version = "2")


pdf(file = paste0("./maps/Threat_maps/ring.richness_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.richness_threat_contrasted, col = pal_grn_red_NA, main = paste0("Mimicry richness\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.richness_threat_contrasted, locs = seq(0, 0.2, 0.05), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")
par(mar = internal_margins)
dev.off()

### Weighted Index

# Load threat index
ring.richness_threat_raw <- readRDS(file = "./outputs/Threat_maps/Weighted_indices/ring.richness_threat.rds")

hist(ring.richness_threat_raw)

# Contrasted version
tot.sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
ring.richness_threat_raw_contrasted <- ring.richness_threat_raw
ring.richness_threat_raw_contrasted[ring.richness_threat_raw[] > 6] <- 6
ring.richness_threat_raw_contrasted[ring.richness_threat_raw[] < 0.05] <- 0.05 # Avoid null values to be colored as ones where no Ithomiini are present
ring.richness_threat_raw_contrasted[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

save(ring.richness_threat_raw_contrasted, file = "./outputs/Threat_maps/Weighted_indices/ring.richness_threat_raw_contrasted.RData", version = "2")
saveRDS(ring.richness_threat_raw_contrasted, file = "./outputs/Threat_maps/Weighted_indices/ring.richness_threat_raw_contrasted.rds", version = "2")


pdf(file = paste0("./maps/Threat_maps/Weighted_indices/ring.richness_threat_raw.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.richness_threat_raw_contrasted, col = pal_grn_red_NA, main = paste0("Mimicry richness\nThreat level"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mimicry\n         rings", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.richness_threat_raw_contrasted, locs = seq(0, 6, 2), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6, font = 2, cex = 1.2, label = "Mimicry\nrings")
par(mar = internal_margins)
dev.off()


### 7.3/ Plot zoom on Andes ####

# Load and crop Indices maps
ring.richness_threat <- readRDS(file = "./outputs/Threat_maps/ring.richness_threat.rds")
Andes.ring.richness_threat <- crop(ring.richness_threat, Andes_ext)

# Contrasted version
hist(Andes.ring.richness_threat)

Andes.ring.richness_threat_contrasted <- Andes.ring.richness_threat
Andes.ring.richness_threat_contrasted[Andes.ring.richness_threat[] > 0.30] <- 0.30

save(Andes.ring.richness_threat_contrasted, file = "./outputs/Threat_maps/Andes/Andes.ring.richness_threat_contrasted.RData", version = "2")
saveRDS(Andes.ring.richness_threat_contrasted, file = "./outputs/Threat_maps/Andes/Andes.ring.richness_threat_contrasted.rds", version = "2")

# Export plot
pdf(file = paste0("./maps/Threat_maps/Andes/Andes.ring.richness_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Andes.ring.richness_threat_contrasted, col = pal_grn_red, main = paste0("Mimicry richness\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Andes_bbox, lwd = 1, border = "grey20", add = T)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
        drawlabels = F, col = "black", add = T)
addRasterLegend(Andes.ring.richness_threat_contrasted, locs = seq(0, 0.3, 0.05), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.75, 4.25))
graphics::text(x = -86.7, y = 6, font = 2, cex = 0.9, label = "Threat index")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 0.9, text.font = 2)
graphics::text(x = -86.7, y = -10.8, font = 2, cex = 0.9, label = "Elevation")
scalebar(d = 750, type = "line", lwd = 4, divs = 4, xy = c(-67, -13.5), label = c("", "750 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.6, padin = c(0.2, 0.4), pos = "topleft")
par(mar = internal_margins)
dev.off()


### 7.4/ Plot zoom on Mata Atlantica ####

# Load and crop Indices maps
ring.richness_threat <- readRDS(file = "./outputs/Threat_maps/ring.richness_threat.rds")
Mata.ring.richness_threat <- crop(ring.richness_threat, Mata_ext)

# Contrasted version
hist(Mata.ring.richness_threat)

Mata.ring.richness_threat_contrasted <- Mata.ring.richness_threat
Mata.ring.richness_threat_contrasted[Mata.ring.richness_threat[] > 0.1] <- 0.1

save(Mata.ring.richness_threat_contrasted, file = "./outputs/Threat_maps/Mata/Mata.ring.richness_threat_contrasted.RData", version = "2")
saveRDS(Mata.ring.richness_threat_contrasted, file = "./outputs/Threat_maps/Mata/Mata.ring.richness_threat_contrasted.rds", version = "2")

# Export plot
pdf(file = paste0("./maps/Threat_maps/Mata/Mata.ring.richness_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Mata.ring.richness_threat_contrasted, col = pal_grn_red, main = paste0("Mimicry richness\nThreat index"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2),
      legend  = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Mata_bbox, lwd = 1, border = "grey20", add = T)
plot(Mata_borders, lwd = 2, border = "black", add = T)
addRasterLegend(Mata.ring.richness_threat_contrasted, locs = seq(0, 0.1, 0.025), digits = 3, cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-33.22, -32.2, -27, -13))
graphics::text(x = -32.8, y = -10.5, font = 2, cex = 1.1, label = "Threat\nindex")
legend(x = "bottomright", legend = c("Mata Atlantica"), lty = c(1) , col = "black", lwd = 2, bty = "n", cex = 1.1, text.font = 2)
scalebar(d = 1000, type = "line", lwd = 4, divs = 4, xy = c(-59, -7), label = c("", "1000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15))
par(mar = internal_margins)
dev.off()


### 8/ Mimicry diversity ####

# Load indice map
ring.diversity.compatible_Jost_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jost_Jaccard.80.rds"))

# Normalized indice on global scale
min <- as.numeric(names(table(ring.diversity.compatible_Jost_Jaccard.80[])[2])) # Do not include communities with null values in the scaling
max <- max(ring.diversity.compatible_Jost_Jaccard.80[], na.rm = T)
ring.diversity_norm <- (ring.diversity.compatible_Jost_Jaccard.80 - min) / (max - min)
ring.diversity_norm[ring.diversity_norm[] < 0.002] <- 0.002 # Add minimal value for communities were the index was null
ring.diversity_norm[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

# plot(ring.diversity_norm, col = pal_grn_red_NA)

### 8.1/ Compute threat map weighted by species richness ####
ring.diversity_threat <- ring.diversity_norm * HII_norm
ring.diversity_threat_raw <- ring.diversity.compatible_Jost_Jaccard.80 * HII_norm

# Save
save(ring.diversity_threat, file = "./outputs/Threat_maps/ring.diversity_threat.RData", version = "2")
saveRDS(ring.diversity_threat, file = "./outputs/Threat_maps/ring.diversity_threat.rds", version = "2")
save(ring.diversity_threat_raw, file = "./outputs/Threat_maps/Weighted_indices/ring.diversity_threat_raw.RData", version = "2")
saveRDS(ring.diversity_threat_raw, file = "./outputs/Threat_maps/Weighted_indices/ring.diversity_threat_raw.rds", version = "2")


### 8.2/ Plot global map ####

### Weighted HII

# Load threat index
ring.diversity_threat <- readRDS(file = "./outputs/Threat_maps/ring.diversity_threat.rds")

hist(ring.diversity_threat)

# Contrasted version
tot.sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
ring.diversity_threat_contrasted <- ring.diversity_threat
ring.diversity_threat_contrasted[ring.diversity_threat[] > 0.2] <- 0.2
ring.diversity_threat_contrasted[ring.diversity_threat[] < 0.001] <- 0.001 # Avoid null values to be colored as ones where no Ithomiini are present
ring.diversity_threat_contrasted[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

save(ring.diversity_threat_contrasted, file = "./outputs/Threat_maps/ring.diversity_threat_contrasted.RData", version = "2")
saveRDS(ring.diversity_threat_contrasted, file = "./outputs/Threat_maps/ring.diversity_threat_contrasted.rds", version = "2")


pdf(file = paste0("./maps/Threat_maps/ring.diversity_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.diversity_threat_contrasted, col = pal_grn_red_NA, main = paste0("Mimicry diversity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.diversity_threat_contrasted, locs = seq(0, 0.2, 0.05), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")
par(mar = internal_margins)
dev.off()

### Weighted Index

# Load threat index
ring.diversity_threat_raw <- readRDS(file = "./outputs/Threat_maps/Weighted_indices/ring.diversity_threat.rds")

hist(ring.diversity_threat_raw)

# Contrasted version
tot.sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
ring.diversity_threat_raw_contrasted <- ring.diversity_threat_raw
ring.diversity_threat_raw_contrasted[ring.diversity_threat_raw[] > 5] <- 5
ring.diversity_threat_raw_contrasted[ring.diversity_threat_raw[] < 0.05] <- 0.05 # Avoid null values to be colored as ones where no Ithomiini are present
ring.diversity_threat_raw_contrasted[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

save(ring.diversity_threat_raw_contrasted, file = "./outputs/Threat_maps/Weighted_indices/ring.diversity_threat_raw_contrasted.RData", version = "2")
saveRDS(ring.diversity_threat_raw_contrasted, file = "./outputs/Threat_maps/Weighted_indices/ring.diversity_threat_raw_contrasted.rds", version = "2")


pdf(file = paste0("./maps/Threat_maps/Weighted_indices/ring.diversity_threat_raw.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.diversity_threat_raw_contrasted, col = pal_grn_red_NA, main = paste0("Mimicry diversity\nThreat level"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mimicry\n         rings", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.diversity_threat_raw_contrasted, locs = seq(0, 5, 1), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Mimicry\nrings")
par(mar = internal_margins)
dev.off()


### 8.3/ Plot zoom on Andes ####

# Load and crop Indices maps
ring.diversity_threat <- readRDS(file = "./outputs/Threat_maps/ring.diversity_threat.rds")
Andes.ring.diversity_threat <- crop(ring.diversity_threat, Andes_ext)

# Contrasted version
hist(Andes.ring.diversity_threat)

Andes.ring.diversity_threat_contrasted <- Andes.ring.diversity_threat
Andes.ring.diversity_threat_contrasted[Andes.ring.diversity_threat[] > 0.25] <- 0.25

save(Andes.ring.diversity_threat_contrasted, file = "./outputs/Threat_maps/Andes/Andes.ring.diversity_threat_contrasted.RData", version = "2")
saveRDS(Andes.ring.diversity_threat_contrasted, file = "./outputs/Threat_maps/Andes/Andes.ring.diversity_threat_contrasted.rds", version = "2")

# Export plot
pdf(file = paste0("./maps/Threat_maps/Andes/Andes.ring.diversity_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Andes.ring.diversity_threat_contrasted, col = pal_grn_red, main = paste0("Mimicry diversity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Andes_bbox, lwd = 1, border = "grey20", add = T)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
        drawlabels = F, col = "black", add = T)
addRasterLegend(Andes.ring.diversity_threat_contrasted, locs = seq(0, 0.25, 0.05), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.75, 4.25))
graphics::text(x = -86.7, y = 6, font = 2, cex = 0.9, label = "Threat index")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 0.9, text.font = 2)
graphics::text(x = -86.7, y = -10.8, font = 2, cex = 0.9, label = "Elevation")
scalebar(d = 750, type = "line", lwd = 4, divs = 4, xy = c(-67, -13.5), label = c("", "750 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.6, padin = c(0.2, 0.4), pos = "topleft")
par(mar = internal_margins)
dev.off()


### 8.4/ Plot zoom on Mata Atlantica ####

# Load and crop Indices maps
ring.diversity_threat <- readRDS(file = "./outputs/Threat_maps/ring.diversity_threat.rds")
Mata.ring.diversity_threat <- crop(ring.diversity_threat, Mata_ext)

# Contrasted version
hist(Mata.ring.diversity_threat)

Mata.ring.diversity_threat_contrasted <- Mata.ring.diversity_threat
Mata.ring.diversity_threat_contrasted[Mata.ring.diversity_threat[] > 0.1] <- 0.1

save(Mata.ring.diversity_threat_contrasted, file = "./outputs/Threat_maps/Mata/Mata.ring.diversity_threat_contrasted.RData", version = "2")
saveRDS(Mata.ring.diversity_threat_contrasted, file = "./outputs/Threat_maps/Mata/Mata.ring.diversity_threat_contrasted.rds", version = "2")

# Export plot
pdf(file = paste0("./maps/Threat_maps/Mata/Mata.ring.diversity_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Mata.ring.diversity_threat_contrasted, col = pal_grn_red, main = paste0("Mimicry diversity\nThreat index"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2),
      legend  = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Mata_bbox, lwd = 1, border = "grey20", add = T)
plot(Mata_borders, lwd = 2, border = "black", add = T)
addRasterLegend(Mata.ring.diversity_threat_contrasted, locs = seq(0, 0.1, 0.025), digits = 3, cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-33.22, -32.2, -27, -13))
graphics::text(x = -32.8, y = -10.5, font = 2, cex = 1.1, label = "Threat\nindex")
legend(x = "bottomright", legend = c("Mata Atlantica"), lty = c(1) , col = "black", lwd = 2, bty = "n", cex = 1.1, text.font = 2)
scalebar(d = 1000, type = "line", lwd = 4, divs = 4, xy = c(-59, -7), label = c("", "1000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15))
par(mar = internal_margins)
dev.off()


### 9/ Mimicry rarity ####

# Load indice map
mimicry.mean.rarity_Jaccard.80 <- readRDS(file = "./outputs/Indices_maps/mimicry.mean.rarity_Jaccard.80.rds")

# Normalized indice on global scale
min <- as.numeric(names(table(mimicry.mean.rarity_Jaccard.80[])[2])) # Do not include communities with null values in the scaling
max <- max(mimicry.mean.rarity_Jaccard.80[], na.rm = T)
ring.rarity_norm <- (mimicry.mean.rarity_Jaccard.80 - min) / (max - min)
ring.rarity_norm[ring.rarity_norm[] < 0.002] <- 0.002 # Add minimal value for communities were the index was null
ring.rarity_norm[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

# plot(ring.rarity_norm, col = pal_grn_red_NA)

save(ring.rarity_norm, file = "./outputs/Threat_maps/ring.rarity_norm.RData", version = "2")
saveRDS(ring.rarity_norm, file = "./outputs/Threat_maps/ring.rarity_norm.rds", version = "2")


### 9.1/ Compute threat map weighted by species richness ####
ring.rarity_threat <- ring.rarity_norm * HII_norm
ring.rarity_threat_raw <- mimicry.mean.rarity_Jaccard.80 * HII_norm

# Save
save(ring.rarity_threat, file = "./outputs/Threat_maps/ring.rarity_threat.RData", version = "2")
saveRDS(ring.rarity_threat, file = "./outputs/Threat_maps/ring.rarity_threat.rds", version = "2")
save(ring.rarity_threat_raw, file = "./outputs/Threat_maps/Weighted_indices/ring.rarity_threat_raw.RData", version = "2")
saveRDS(ring.rarity_threat_raw, file = "./outputs/Threat_maps/Weighted_indices/ring.rarity_threat_raw.rds", version = "2")

### 9.2/ Plot global map ####

### Weighted HII

# Load threat index
ring.rarity_threat <- readRDS(file = "./outputs/Threat_maps/ring.rarity_threat.rds")

hist(ring.rarity_threat)

# Contrasted version
tot.sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
ring.rarity_threat_contrasted <- ring.rarity_threat
ring.rarity_threat_contrasted[ring.rarity_threat[] > 0.25] <- 0.25
ring.rarity_threat_contrasted[ring.rarity_threat[] < 0.002] <- 0.002 # Avoid null values to be colored as ones where no Ithomiini are present
ring.rarity_threat_contrasted[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

save(ring.rarity_threat_contrasted, file = "./outputs/Threat_maps/ring.rarity_threat_contrasted.RData", version = "2")
saveRDS(ring.rarity_threat_contrasted, file = "./outputs/Threat_maps/ring.rarity_threat_contrasted.rds", version = "2")


pdf(file = paste0("./maps/Threat_maps/ring.rarity_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.rarity_threat_contrasted, col = pal_grn_red_NA, main = paste0("Mean mimicry rarity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.rarity_threat_contrasted, locs = seq(0, 0.25, 0.05), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")
par(mar = internal_margins)
dev.off()


### Weighted Index

# Load threat index
ring.rarity_threat_raw <- readRDS(file = "./outputs/Threat_maps/Weighted_indices/ring.rarity_threat_raw.rds")

hist(ring.rarity_threat_raw)

# Contrasted version
tot.sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
ring.rarity_threat_raw_contrasted <- ring.rarity_threat_raw
ring.rarity_threat_raw_contrasted[ring.rarity_threat_raw[] > 0.25] <- 0.25
ring.rarity_threat_raw_contrasted[ring.rarity_threat_raw[] < 0.002] <- 0.002 # Avoid null values to be colored as ones where no Ithomiini are present
ring.rarity_threat_raw_contrasted[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent

save(ring.rarity_threat_raw_contrasted, file = "./outputs/Threat_maps/Weighted_indices/ring.rarity_threat_raw_contrasted.RData", version = "2")
saveRDS(ring.rarity_threat_raw_contrasted, file = "./outputs/Threat_maps/Weighted_indices/ring.rarity_threat_raw_contrasted.rds", version = "2")


pdf(file = paste0("./maps/Threat_maps/Weighted_indices/ring.rarity_threat_raw.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(ring.rarity_threat_raw_contrasted, col = pal_grn_red_NA, main = paste0("Mean mimicry rarity\nThreat level"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.rarity_threat_raw_contrasted, locs = seq(0, 0.25, 0.05), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Rarity\nindex")
par(mar = internal_margins)
dev.off()



### 9.3/ Plot zoom on Andes ####

# Load and crop Indices maps
ring.rarity_threat <- readRDS(file = "./outputs/Threat_maps/ring.rarity_threat.rds")
Andes.ring.rarity_threat <- crop(ring.rarity_threat, Andes_ext)

# Contrasted version
hist(Andes.ring.rarity_threat)

Andes.ring.rarity_threat_contrasted <- Andes.ring.rarity_threat
Andes.ring.rarity_threat_contrasted[Andes.ring.rarity_threat[] > 0.25] <- 0.25

save(Andes.ring.rarity_threat_contrasted, file = "./outputs/Threat_maps/Andes/Andes.ring.rarity_threat_contrasted.RData", version = "2")
saveRDS(Andes.ring.rarity_threat_contrasted, file = "./outputs/Threat_maps/Andes/Andes.ring.rarity_threat_contrasted.rds", version = "2")

# Export plot
pdf(file = paste0("./maps/Threat_maps/Andes/Andes.ring.rarity_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Andes.ring.rarity_threat_contrasted, col = pal_grn_red, main = paste0("Mean mimicry rarity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
# plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Andes_bbox, lwd = 1, border = "grey20", add = T)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
        drawlabels = F, col = "black", add = T)
addRasterLegend(Andes.ring.rarity_threat_contrasted, locs = seq(0, 0.25, 0.05), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.75, 4.25))
graphics::text(x = -86.7, y = 6, font = 2, cex = 0.9, label = "Threat index")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 0.9, text.font = 2)
graphics::text(x = -86.7, y = -10.8, font = 2, cex = 0.9, label = "Elevation")
scalebar(d = 750, type = "line", lwd = 4, divs = 4, xy = c(-67, -13.5), label = c("", "750 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.6, padin = c(0.2, 0.4), pos = "topleft")
par(mar = internal_margins)
dev.off()


### 9.4/ Plot zoom on Mata Atlantica ####

# Load and crop Indices maps
ring.rarity_threat <- readRDS(file = "./outputs/Threat_maps/ring.rarity_threat.rds")
Mata.ring.rarity_threat <- crop(ring.rarity_threat, Mata_ext)

# Contrasted version
hist(Mata.ring.rarity_threat)

Mata.ring.rarity_threat_contrasted <- Mata.ring.rarity_threat
Mata.ring.rarity_threat_contrasted[Mata.ring.rarity_threat[] > 0.15] <- 0.15

save(Mata.ring.rarity_threat_contrasted, file = "./outputs/Threat_maps/Mata/Mata.ring.rarity_threat_contrasted.RData", version = "2")
saveRDS(Mata.ring.rarity_threat_contrasted, file = "./outputs/Threat_maps/Mata/Mata.ring.rarity_threat_contrasted.rds", version = "2")

# Export plot
pdf(file = paste0("./maps/Threat_maps/Mata/Mata.ring.rarity_threat.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Mata.ring.rarity_threat_contrasted, col = pal_grn_red, main = paste0("Mean mimicry rarity\nThreat index"),
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2),
      legend  = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask_pixel, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(Mata_bbox, lwd = 1, border = "grey20", add = T)
plot(Mata_borders, lwd = 2, border = "black", add = T)
addRasterLegend(Mata.ring.rarity_threat_contrasted, locs = seq(0, 0.15, 0.03), digits = 3, cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-33.22, -32.2, -27, -13))
graphics::text(x = -32.8, y = -10.5, font = 2, cex = 1.1, label = "Threat\nindex")
legend(x = "bottomright", legend = c("Mata Atlantica"), lty = c(1) , col = "black", lwd = 2, bty = "n", cex = 1.1, text.font = 2)
scalebar(d = 1000, type = "line", lwd = 4, divs = 4, xy = c(-59, -7), label = c("", "1000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15))
par(mar = internal_margins)
dev.off()

#####



### 10/ Multiple plots with the 9 indices ####

## 10.1/ Full maps ####

# 10.1.1/ Independant contrast ####

pdf(file = paste0("./maps/Threat_maps/All_indices_threat.pdf"), height = 15.9, width = 19.5)
par(mfrow = c(3, 3))

# 1: Species richness
sp.richness_threat_contrasted <- readRDS(file = "./outputs/Threat_maps/sp.richness_threat_contrasted.rds")
image(sp.richness_threat_contrasted, col = pal_grn_red_NA, main = paste0("Species richness\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.richness_threat_contrasted, locs = seq(0, 0.15, 0.03), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")

# 2: Species Diversity
sp.diversity_threat_contrasted <- readRDS(file = "./outputs/Threat_maps/sp.diversity_threat_contrasted.rds")
image(sp.diversity_threat_contrasted, col = pal_grn_red_NA, main = paste0("Species diversity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.diversity_threat_contrasted, locs = seq(0, 0.15, 0.03), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")

# 3: Species rarity
sp.rarity_threat_contrasted <- readRDS(file = "./outputs/Threat_maps/sp.rarity_threat_contrasted.rds")
image(sp.rarity_threat_contrasted, col = pal_grn_red_NA, main = paste0("Mean species rarity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.rarity_threat_contrasted, locs = seq(0, 0.3, 0.05), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")

# 4: Faith's Phylogenetic Diversity
PD_threat_contrasted <- readRDS(file = "./outputs/Threat_maps/PD_threat_contrasted.rds")
image(PD_threat_contrasted, col = pal_grn_red_NA, main = paste0("Faith's Phylogenetic Diversity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(PD_threat_contrasted, locs = seq(0, 0.2, 0.05), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")

# 5: Sum of Fair-Proportions
FP_threat_contrasted <- readRDS(file = "./outputs/Threat_maps/FP_threat_contrasted.rds")
image(FP_threat_contrasted, col = pal_grn_red_NA, main = paste0("Sum of Fair-proportions\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(FP_threat_contrasted, locs = seq(0, 0.15, 0.03), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")

# 6: MPD
MPD_threat_contrasted <- readRDS(file = "./outputs/Threat_maps/MPD_threat_contrasted.rds")
image(MPD_threat_contrasted, col = pal_grn_red_NA, main = paste0("Mean pairwise Phylogenetic Distance\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(MPD_threat_contrasted, locs = seq(0, 0.4, 0.1), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")

# 7: Mimicry richness
ring.richness_threat_contrasted <- readRDS(file = "./outputs/Threat_maps/ring.richness_threat_contrasted.rds")
image(ring.richness_threat_contrasted, col = pal_grn_red_NA, main = paste0("Mimicry richness\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.richness_threat_contrasted, locs = seq(0, 0.2, 0.05), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")

# 8: Mimicry diversity
ring.diversity_threat_contrasted <- readRDS(file = "./outputs/Threat_maps/ring.diversity_threat_contrasted.rds")
image(ring.diversity_threat_contrasted, col = pal_grn_red_NA, main = paste0("Mimicry diversity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.diversity_threat_contrasted, locs = seq(0, 0.2, 0.05), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")

# 9: Mimicry rarity
ring.rarity_threat_contrasted <- readRDS(file = "./outputs/Threat_maps/ring.rarity_threat_contrasted.rds")
image(ring.rarity_threat_contrasted, col = pal_grn_red_NA, main = paste0("Mean mimicry rarity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.rarity_threat_contrasted, locs = seq(0, 0.25, 0.05), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")

par(mfrow = c(1, 1))

dev.off()

# 10.1.2/ Unique contrast version ####

# Max. contrast = 0.2
max.contrast <- 0.2

tot.sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))

pdf(file = paste0("./maps/Threat_maps/All_indices_threat_contrast_0.2.pdf"), height = 15.9, width = 19.5)
par(mfrow = c(3, 3))

# 1: Species richness
sp.richness_threat <- readRDS(file = "./outputs/Threat_maps/sp.richness_threat.rds")
# Contrasted version
sp.richness_threat_contrasted <- sp.richness_threat
sp.richness_threat_contrasted[sp.richness_threat[] > max.contrast] <- max.contrast
sp.richness_threat_contrasted[sp.richness_threat[] < 0.002] <- 0.002 # Avoid null values to be colored as ones where no Ithomiini are present
sp.richness_threat_contrasted[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent
image(sp.richness_threat_contrasted, col = pal_grn_red_NA, main = paste0("Species richness\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.richness_threat_contrasted, locs = seq(0, max.contrast, floor(max.contrast/4*100)/100), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")

# 2: Species Diversity
sp.diversity_threat <- readRDS(file = "./outputs/Threat_maps/sp.diversity_threat.rds")
# Contrasted version
sp.diversity_threat_contrasted <- sp.diversity_threat
sp.diversity_threat_contrasted[sp.diversity_threat[] > max.contrast] <- max.contrast
sp.diversity_threat_contrasted[sp.diversity_threat[] < 0.001] <- 0.001 # Avoid null values to be colored as ones where no Ithomiini are present
sp.diversity_threat_contrasted[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent
image(sp.diversity_threat_contrasted, col = pal_grn_red_NA, main = paste0("Species diversity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.diversity_threat_contrasted, locs = seq(0, max.contrast, floor(max.contrast/4*100)/100), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")

# 3: Species rarity
sp.rarity_threat <- readRDS(file = "./outputs/Threat_maps/sp.rarity_threat.rds")
# Contrasted version
sp.rarity_threat_contrasted <- sp.rarity_threat
sp.rarity_threat_contrasted[sp.rarity_threat[] > max.contrast] <- max.contrast
sp.rarity_threat_contrasted[sp.rarity_threat[] < 0.002] <- 0.002 # Avoid null values to be colored as ones where no Ithomiini are present
sp.rarity_threat_contrasted[sp.mean.rarity_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent
image(sp.rarity_threat_contrasted, col = pal_grn_red_NA, main = paste0("Mean species rarity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.rarity_threat_contrasted, locs = seq(0, max.contrast, floor(max.contrast/4*100)/100), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")

# 4: Faith's Phylogenetic Diversity
PD_threat <- readRDS(file = "./outputs/Threat_maps/PD_threat.rds")
# Contrasted version
PD_threat_contrasted <- PD_threat
PD_threat_contrasted[PD_threat[] > max.contrast] <- max.contrast
PD_threat_contrasted[PD_threat[] < 0.001] <- 0.001 # Avoid null values to be colored as ones where no Ithomiini are present
PD_threat_contrasted[sp.mean.rarity_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent
image(PD_threat_contrasted, col = pal_grn_red_NA, main = paste0("Faith's Phylogenetic Diversity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(PD_threat_contrasted, locs = seq(0, max.contrast, floor(max.contrast/4*100)/100), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")

# 5: Sum of Fair-Proportions
FP_threat <- readRDS(file = "./outputs/Threat_maps/FP_threat.rds")
# Contrasted version
FP_threat_contrasted <- FP_threat
FP_threat_contrasted[FP_threat[] > max.contrast] <- max.contrast
FP_threat_contrasted[FP_threat[] < 0.001] <- 0.001 # Avoid null values to be colored as ones where no Ithomiini are present
FP_threat_contrasted[sp.mean.rarity_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent
image(FP_threat_contrasted, col = pal_grn_red_NA, main = paste0("Sum of Fair-proportions\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(FP_threat_contrasted, locs = seq(0, max.contrast, floor(max.contrast/4*100)/100), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")

# 6: MPD
MPD_threat <- readRDS(file = "./outputs/Threat_maps/MPD_threat.rds")
# Contrasted version
MPD_threat_contrasted <- MPD_threat
MPD_threat_contrasted[MPD_threat[] > max.contrast] <- max.contrast
image(MPD_threat_contrasted, col = pal_grn_red_NA, main = paste0("Mean pairwise Phylogenetic Distance\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(MPD_threat_contrasted, locs = seq(0, max.contrast, floor(max.contrast/4*100)/100), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")

# 7: Mimicry richness
ring.richness_threat <- readRDS(file = "./outputs/Threat_maps/ring.richness_threat.rds")
# Contrasted version
ring.richness_threat_contrasted <- ring.richness_threat
ring.richness_threat_contrasted[ring.richness_threat[] > max.contrast] <- max.contrast
ring.richness_threat_contrasted[ring.richness_threat[] < 0.001] <- 0.001 # Avoid null values to be colored as ones where no Ithomiini are present
ring.richness_threat_contrasted[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent
image(ring.richness_threat_contrasted, col = pal_grn_red_NA, main = paste0("Mimicry richness\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.richness_threat_contrasted, locs = seq(0, max.contrast, floor(max.contrast/4*100)/100), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")

# 8: Mimicry diversity
ring.diversity_threat <- readRDS(file = "./outputs/Threat_maps/ring.diversity_threat.rds")
# Contrasted version
ring.diversity_threat_contrasted <- ring.diversity_threat
ring.diversity_threat_contrasted[ring.diversity_threat[] > max.contrast] <- max.contrast
ring.diversity_threat_contrasted[ring.diversity_threat[] < 0.001] <- 0.001 # Avoid null values to be colored as ones where no Ithomiini are present
ring.diversity_threat_contrasted[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent
image(ring.diversity_threat_contrasted, col = pal_grn_red_NA, main = paste0("Mimicry diversity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.diversity_threat_contrasted, locs = seq(0, max.contrast, floor(max.contrast/4*100)/100), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")

# 9: Mimicry rarity
ring.rarity_threat <- readRDS(file = "./outputs/Threat_maps/ring.rarity_threat.rds")
# Contrasted version
ring.rarity_threat_contrasted <- ring.rarity_threat
ring.rarity_threat_contrasted[ring.rarity_threat[] > max.contrast] <- max.contrast
ring.rarity_threat_contrasted[ring.rarity_threat[] < 0.002] <- 0.002 # Avoid null values to be colored as ones where no Ithomiini are present
ring.rarity_threat_contrasted[tot.sp.richness_Jaccard.80[] == 0] <- 0 # Add null values where Ithomiini are absent
image(ring.rarity_threat_contrasted, col = pal_grn_red_NA, main = paste0("Mean mimicry rarity\nThreat index"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Threat\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.rarity_threat_contrasted, locs = seq(0, max.contrast, floor(max.contrast/4*100)/100), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Threat\nindex")

par(mfrow = c(1, 1))

dev.off()




# Andes

# Mata


### Multiple plot with the 4 main indices: 3 columns (Full, Andes, Mata)
# Same contrast
# Independant contrast


hist(MPD_norm[MPD_norm>0])

hist(sp.rarity_norm[sp.rarity_norm>0])

hist(FP_norm[FP_norm>0])


# 10.1.3/ Weighted indices ####

pdf(file = paste0("./maps/Threat_maps/Weighted_indices/All_indices_threat_raw.pdf"), height = 15.9, width = 19.5)
par(mfrow = c(3, 3))

# 1: Species richness
sp.richness_threat_raw_contrasted <- readRDS(file = "./outputs/Threat_maps/Weighted_indices/sp.richness_threat_raw_contrasted.rds")
image(sp.richness_threat_raw_contrasted, col = pal_grn_red_NA, main = paste0("Species richness\nThreat level"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.richness_threat_raw_contrasted, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 5, font = 2, cex = 1.2, label = "Species")

# 2: Species Diversity
sp.diversity_threat_raw_contrasted <- readRDS(file = "./outputs/Threat_maps/Weighted_indices/sp.diversity_threat_raw_contrasted.rds")
image(sp.diversity_threat_raw_contrasted, col = pal_grn_red_NA, main = paste0("Species diversity\nThreat level"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Species", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.diversity_threat_raw_contrasted, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 5.5, font = 2, cex = 1.2, label = "Species")

# 3: Species rarity
sp.rarity_threat_raw_contrasted <- readRDS(file = "./outputs/Threat_maps/Weighted_indices/sp.rarity_threat_raw_contrasted.rds")
image(sp.rarity_threat_raw_contrasted, col = pal_grn_red_NA, main = paste0("Mean species rarity\nThreat level"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="        Rarity\n         Index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(sp.rarity_threat_raw_contrasted, locs = seq(0, 0.35, 0.1), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Rarity\nindex")

# 4: Faith's Phylogenetic Diversity
PD_threat_raw_contrasted <- readRDS(file = "./outputs/Threat_maps/Weighted_indices/PD_threat_raw_contrasted.rds")
image(PD_threat_raw_contrasted, col = pal_grn_red_NA, main = paste0("Faith's Phylogenetic Diversity\nThreat level"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(PD_threat_raw_contrasted, locs = seq(0, 180, 40), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -109, y = 6.5, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

# 5: Sum of Fair-Proportions
FP_threat_raw_contrasted <- readRDS(file = "./outputs/Threat_maps/Weighted_indices/FP_threat_raw_contrasted.rds")
image(FP_threat_raw_contrasted, col = pal_grn_red_NA, main = paste0("Sum of Fair-proportions\nThreat level"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(FP_threat_raw_contrasted, locs = seq(0, 75, 20), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -109, y = 6.5, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

# 6: MPD
MPD_threat_raw_contrasted <- readRDS(file = "./outputs/Threat_maps/Weighted_indices/MPD_threat_raw_contrasted.rds")
image(MPD_threat_raw_contrasted, col = pal_grn_red_NA, main = paste0("Mean pairwise Phylogenetic Distance\nThreat level"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Evolutionary\n         Time (My)", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(MPD_threat_raw_contrasted, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -109, y = 6.5, font = 2, cex = 1.2, label = "Evolutionary\nTime (My)")

# 7: Mimicry richness
ring.richness_threat_raw_contrasted <- readRDS(file = "./outputs/Threat_maps/Weighted_indices/ring.richness_threat_raw_contrasted.rds")
image(ring.richness_threat_raw_contrasted, col = pal_grn_red_NA, main = paste0("Mimicry richness\nThreat level"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mimicry\n         rings", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.richness_threat_raw_contrasted, locs = seq(0, 6, 2), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6, font = 2, cex = 1.2, label = "Mimicry\nrings")

# 8: Mimicry diversity
ring.diversity_threat_raw_contrasted <- readRDS(file = "./outputs/Threat_maps/Weighted_indices/ring.diversity_threat_raw_contrasted.rds")
image(ring.diversity_threat_raw_contrasted, col = pal_grn_red_NA, main = paste0("Mimicry diversity\nThreat level"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mimicry\n         rings", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.diversity_threat_raw_contrasted, locs = seq(0, 5, 1), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Mimicry\nrings")

# 9: Mimicry rarity
ring.rarity_threat_raw_contrasted <- readRDS(file = "./outputs/Threat_maps/Weighted_indices/ring.rarity_threat_raw_contrasted.rds")
image(ring.rarity_threat_raw_contrasted, col = pal_grn_red_NA, main = paste0("Mean mimicry rarity\nThreat level"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Rarity\n           index", cex=1.2, line = 1, font = 2), 
      legend  = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# abline(h = 0, lty = "59", lwd = 0.2, col = "black")
# abline(h = 23.43665, lty = "29", lwd = 0.2, col = "black")
# abline(h = -23.43665, lty = "29", lwd = 0.2, col = "black")
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
addRasterLegend(ring.rarity_threat_raw_contrasted, locs = seq(0, 0.25, 0.05), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Rarity\nindex")

par(mfrow = c(1, 1))

dev.off()


### 11/ Boxplot of weighted HII ####

sp.richness_threat <- readRDS(file = "./outputs/Threat_maps/sp.richness_threat.rds")
sp.rarity_threat <- readRDS(file = "./outputs/Threat_maps/sp.rarity_threat.rds")
MPD_threat <- readRDS(file = "./outputs/Threat_maps/MPD_threat.rds")
ring.richness_threat <- readRDS(file = "./outputs/Threat_maps/ring.richness_threat.rds")

sp.richness_threat_data <- sp.richness_threat[sp.richness_threat > 0]
sp.rarity_threat_data <- sp.rarity_threat[sp.rarity_threat > 0]
MPD_threat_data <- MPD_threat[MPD_threat > 0]
ring.richness_threat_data <- ring.richness_threat[ring.richness_threat > 0]

sp.richness_threat_data <- tibble(HII = sp.richness_threat_data) %>% 
  mutate(index = "sp.richness")
sp.rarity_threat_data <- tibble(HII = sp.rarity_threat_data) %>% 
  mutate(index = "sp.rarity")
MPD_threat_data <- tibble(HII = MPD_threat_data) %>% 
  mutate(index = "MPD")
ring.richness_threat_data <- tibble(HII = ring.richness_threat_data) %>% 
  mutate(index = "ring.richness")

weighted_HII_data <- rbind(sp.richness_threat_data, sp.rarity_threat_data, MPD_threat_data, ring.richness_threat_data)
weighted_HII_data <- weighted_HII_data %>% 
  mutate(index = factor(index, levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness")))

# Compute quantiles 95 to add to the plot
quant_95 <- as.numeric(with(weighted_HII_data, by(HII, INDICES = index, FUN = quantile, probs = 0.95)))
quant_95 <- tibble(quant_value = quant_95, index = c("sp.richness", "sp.rarity", "MPD", "ring.richness"))
weighted_HII_data <- weighted_HII_data %>% 
  left_join(y = quant_95, by = "index") %>% 
  mutate(index = factor(index, levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness"))) # Change factor level order

saveRDS(weighted_HII_data, "./outputs/Threat_maps/weighted_HII_data.rds", version = "2")


# Plot the boxplot
pdf(file = paste0("./graphs/Threat_curves/weighted_HII_boxplot.pdf"), height = 5.3, width = 6.5)
g <- ggplot(weighted_HII_data, aes(x = index, y = HII), show.legend = F) +
  geom_boxplot(aes(fill = index)) +
  geom_point(aes(x = index, y = quant_value, col = index), size = 5) +
  labs(y = "weighted HII", x = "Index",
       title = "Full map ; Weighted HII") +
  theme_bw() +
  scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  guides(color = "none")
print(g)
dev.off()
