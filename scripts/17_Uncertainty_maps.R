
##### Script 17: Map uncertainties in the modeling process #####

# Map uncertainties in the modeling process


## Different options to evaluate uncertainties

# 1/ Mean OMU uncertainty (no error propagation ?)

# 2/ Random sampling approach (using only one random submodel per OMU)

# 3/ Full submodels approach (sd between the richness maps computed for each submodels type (i.e., algos))

## Generate a map version for each 4 options between Jaccard/TSS and Buffer.80/95

## Put multiple pdf in main folder and individual maps in subfolder such as /TSS_Buffer.80/




### Prepare stuff

# Effacer l'environnement
rm(list = ls())

library(sf)
library(raster)
library(rangeBuilder)
library(gplots)
library(tidyverse)

# Change temp folder for raster files
rasterOptions(tmpdir = paste0("./temp"))
# To clean regularly temp folder
unlink(list.files(path = rasterOptions()$tmpdir, all.files = T, recursive = T, include.dirs = T, full.names = T), force = T, recursive = T) # Clean raster store in temp


# Load Summary table for OMU/unit and for Species
load(file = paste0("./input_data/list.sp.RData"))
load(file = paste0("./input_data/list.models.RData"))

# Make a list of modeled and non-modeled units (occ.unit = with only occurrences points)
modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]
rasterized_OMU <- list.models[is.na(list.models$Model_ID), ]

# Load shp files for maps
load(file = "./input_data/Map_stuff/bg_mask.RData") # Load bg shp
load(file = "./input_data/Map_stuff/country_borders.RData") # Load country borders
rivers <- st_read("./input_data/Map_stuff/Major_rivers/MajorRivers.shp") # Load rivers

# New color palette
pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")

library(tmaptools)
pal_grn_red <- rev(get_brewer_pal("RdYlGn", n = 200))

# Choose resolution
res <- "15"

# Load mask for continent borders
continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_", res, ".rds"))
crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".rds"))


###### 1/ Mean uncertainty ######

### 1.1/ OMU model uncertainty Stack generation ####
All_OMU_sd_stack_Jaccard.80 <- All_OMU_sd_stack_Jaccard.95 <- All_OMU_sd_stack_TSS.80 <- All_OMU_sd_stack_TSS.95  <- stack(continent_mask) # 1e temp layer used to initiate stack, to remove afterwards
# i <- 1

for (i in 1:nrow(modeled_OMU)) # By modeled OMU
{ 
  # Load OMU name
  unit <- as.character(modeled_OMU$Tag.model[i]) 
  
  # Load uncertainty maps
  Ensemble_sd_Jaccard.80 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_Jaccard_sd_cropped_80.rds"))
  Ensemble_sd_Jaccard.95 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_Jaccard_sd_cropped_95.rds"))
  Ensemble_sd_TSS.80 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_TSS_sd_cropped_80.rds"))
  Ensemble_sd_TSS.95 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_TSS_sd_cropped_95.rds"))

  # Build stack
  All_OMU_sd_stack_Jaccard.80 <- addLayer(All_OMU_sd_stack_Jaccard.80, Ensemble_sd_Jaccard.80)
  All_OMU_sd_stack_Jaccard.95 <- addLayer(All_OMU_sd_stack_Jaccard.95, Ensemble_sd_Jaccard.95)
  All_OMU_sd_stack_TSS.80 <- addLayer(All_OMU_sd_stack_TSS.80, Ensemble_sd_TSS.80)
  All_OMU_sd_stack_TSS.95 <- addLayer(All_OMU_sd_stack_TSS.95, Ensemble_sd_TSS.95)
  
  if (i %% 10 == 0) {print(i)}
}

# Drop the useless first layer used to initiate the stack
All_OMU_sd_stack_Jaccard.80 <- dropLayer(All_OMU_sd_stack_Jaccard.80, i = 1)
All_OMU_sd_stack_Jaccard.95 <- dropLayer(All_OMU_sd_stack_Jaccard.95, i = 1)
All_OMU_sd_stack_TSS.80 <- dropLayer(All_OMU_sd_stack_TSS.80, i = 1)
All_OMU_sd_stack_TSS.95 <- dropLayer(All_OMU_sd_stack_TSS.95, i = 1)

# Name layers with OMU names
names(All_OMU_sd_stack_Jaccard.80) <- names(All_OMU_sd_stack_Jaccard.95) <- names(All_OMU_sd_stack_TSS.80) <- names(All_OMU_sd_stack_TSS.95)<- as.character(modeled_OMU$Tag.model)

# nlayers(All_OMU_sd_stack_Jaccard.80) # 563 modeled OMU in the final stacks

# Save stacks
save(All_OMU_sd_stack_Jaccard.80, file = paste0("./outputs/Indices_stacks/All_OMU_sd_stack_Jaccard.80.RData"), version = "2")
saveRDS(All_OMU_sd_stack_Jaccard.80, file = paste0("./outputs/Indices_stacks/All_OMU_sd_stack_Jaccard.80.rds"), version = "2")
save(All_OMU_sd_stack_Jaccard.95, file = paste0("./outputs/Indices_stacks/All_OMU_sd_stack_Jaccard.95.RData"), version = "2")
saveRDS(All_OMU_sd_stack_Jaccard.95, file = paste0("./outputs/Indices_stacks/All_OMU_sd_stack_Jaccard.95.rds"), version = "2")
save(All_OMU_sd_stack_Jaccard.80, file = paste0("./outputs/Indices_stacks/All_OMU_sd_stack_TSS.80.RData"), version = "2")
saveRDS(All_OMU_sd_stack_Jaccard.80, file = paste0("./outputs/Indices_stacks/All_OMU_sd_stack_TSS.80.rds"), version = "2")
save(All_OMU_sd_stack_TSS.95, file = paste0("./outputs/Indices_stacks/All_OMU_sd_stack_TSS.95.RData"), version = "2")
saveRDS(All_OMU_sd_stack_TSS.95, file = paste0("./outputs/Indices_stacks/All_OMU_sd_stack_TSS.95.rds"), version = "2")

### Load directly the complete stack
All_OMU_sd_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_OMU_sd_stack_Jaccard.80.rds"))
All_OMU_sd_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_OMU_sd_stack_Jaccard.95.rds"))
All_OMU_sd_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_OMU_sd_stack_TSS.80.rds"))
All_OMU_sd_stack_TSS.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_OMU_sd_stack_TSS.95.rds"))


### 1.2/ Mean uncertainty computation ####
mean_uncertainty_Jaccard.80 <- readAll(calc(All_OMU_sd_stack_Jaccard.80, fun = mean, na.rm = T))
mean_uncertainty_Jaccard.95 <- readAll(calc(All_OMU_sd_stack_Jaccard.95, fun = mean, na.rm = T))
mean_uncertainty_TSS.80 <- readAll(calc(All_OMU_sd_stack_TSS.80, fun = mean, na.rm = T))
mean_uncertainty_TSS.95 <- readAll(calc(All_OMU_sd_stack_TSS.95, fun = mean, na.rm = T))

# Save
save(mean_uncertainty_Jaccard.80, file = paste0("./outputs/Uncertainties_maps/mean_uncertainty_Jaccard.80.RData"), version = "2")
saveRDS(mean_uncertainty_Jaccard.80, file = paste0("./outputs/Uncertainties_maps/mean_uncertainty_Jaccard.80.rds"), version = "2")
save(mean_uncertainty_Jaccard.95, file = paste0("./outputs/Uncertainties_maps/mean_uncertainty_Jaccard.95.RData"), version = "2")
saveRDS(mean_uncertainty_Jaccard.95, file = paste0("./outputs/Uncertainties_maps/mean_uncertainty_Jaccard.95.rds"), version = "2")
save(mean_uncertainty_TSS.80, file = paste0("./outputs/Uncertainties_maps/mean_uncertainty_TSS.80.RData"), version = "2")
saveRDS(mean_uncertainty_TSS.80, file = paste0("./outputs/Uncertainties_maps/mean_uncertainty_TSS.80.rds"), version = "2")
save(mean_uncertainty_TSS.95, file = paste0("./outputs/Uncertainties_maps/mean_uncertainty_TSS.95.RData"), version = "2")
saveRDS(mean_uncertainty_TSS.95, file = paste0("./outputs/Uncertainties_maps/mean_uncertainty_TSS.95.rds"), version = "2")


### Load directly the final mean Uncertainty layer
mean_uncertainty_Jaccard.80 <- readRDS(file = paste0("./outputs/Uncertainties_maps/mean_uncertainty_Jaccard.80.rds"))
mean_uncertainty_Jaccard.95 <- readRDS(file = paste0("./outputs/Uncertainties_maps/mean_uncertainty_Jaccard.95.rds"))
mean_uncertainty_TSS.80 <- readRDS(file = paste0("./outputs/Uncertainties_maps/mean_uncertainty_TSS.80.rds"))
mean_uncertainty_TSS.95 <- readRDS(file = paste0("./outputs/Uncertainties_maps/mean_uncertainty_TSS.95.rds"))

### 1.3/ Plot mean uncertainty ####

### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Uncertainties_maps/mean_uncertainty_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(mean_uncertainty_Jaccard.80, col = pal_grn_red, main = paste0("Mean Uncertainty \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mean\n            Uncertainty", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(mean_uncertainty_Jaccard.80, locs = seq(0, 0.08, 0.02), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 5.5, font = 2, cex = 1.2, label = "Mean\nUncertainty")
par(mar = internal_margins)
dev.off()

str(mean_uncertainty_Jaccard.95)

# Jaccard.95
pdf(file = paste0("./maps/Uncertainties_maps/mean_uncertainty_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(mean_uncertainty_Jaccard.95, col = pal_grn_red, main = paste0("Mean Uncertainty \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mean\n            Uncertainty", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(mean_uncertainty_Jaccard.95, locs = seq(0, 0.09, 0.02), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 5.5, font = 2, cex = 1.2, label = "Mean\nUncertainty")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Uncertainties_maps/mean_uncertainty_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(mean_uncertainty_TSS.80, col = pal_grn_red, main = paste0("Mean Uncertainty \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mean\n            Uncertainty", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(mean_uncertainty_TSS.80, locs = seq(0, 0.08, 0.02), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 5.5, font = 2, cex = 1.2, label = "Mean\nUncertainty")
par(mar = internal_margins)
dev.off()

str(mean_uncertainty_TSS.95)

# TSS.95
pdf(file = paste0("./maps/Uncertainties_maps/mean_uncertainty_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(mean_uncertainty_TSS.95, col = pal_grn_red, main = paste0("Mean Uncertainty \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mean\n            Uncertainty", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(mean_uncertainty_TSS.95, locs = seq(0, 0.09, 0.02), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 5.5, font = 2, cex = 1.2, label = "Mean\nUncertainty")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(mean_uncertainty_Jaccard.80) +
#   tm_raster(palette = pal_grn_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Uncertainties_maps/mean_uncertainty_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(mean_uncertainty_Jaccard.80, col = pal_grn_red, main = paste0("Mean Uncertainty \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mean\n            Uncertainty", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(mean_uncertainty_Jaccard.80, locs = seq(0, 0.08, 0.02), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 5.5, font = 2, cex = 1.2, label = "Mean\nUncertainty")

image(mean_uncertainty_Jaccard.95, col = pal_grn_red, main = paste0("Mean Uncertainty \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mean\n            Uncertainty", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(mean_uncertainty_Jaccard.95, locs = seq(0, 0.09, 0.02), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 5.5, font = 2, cex = 1.2, label = "Mean\nUncertainty")

image(mean_uncertainty_TSS.80, col = pal_grn_red, main = paste0("Mean Uncertainty \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mean\n            Uncertainty", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(mean_uncertainty_TSS.80, locs = seq(0, 0.08, 0.02), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 5.5, font = 2, cex = 1.2, label = "Mean\nUncertainty")

image(mean_uncertainty_TSS.95, col = pal_grn_red, main = paste0("Mean Uncertainty \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Mean\n            Uncertainty", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(mean_uncertainty_TSS.95, locs = seq(0, 0.09, 0.02), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 5.5, font = 2, cex = 1.2, label = "Mean\nUncertainty")

par(mar = internal_margins)

dev.off()


# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])



###### 2/ Richness Uncertainty : random sampling approach ######

# Load environmental stack to use as reference for extent and CRS and generate a mask for continental borders
continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_", res, ".rds"))

# Function to compute probability of presence with multiple OMUs
aggreg_prob = function(x, na.rm) { 
  y <- 1-prod(1-x) # Probability of presence of species = probability of presence of at least one OMU = opposite of probability of absence of all OMU
  return(y) # Output
}

### 2.1/ Make stack of richness maps obtained from a single randomly drawn submodel for each OMU ####

set.seed(513123)

# Loop to do 100 maps to get a sample big enough to estimate sd with a minimal error

Random_richness_stack_Jaccard.80 <- Random_richness_stack_Jaccard.95 <- Random_richness_stack_TSS.80 <- Random_richness_stack_TSS.95 <- stack(continent_mask) # 1e temp layer used to initiate stack, to remove afterwards
# for (k in 1:100) 
for (k in 1:25) 
{
  # k <- 1
  
  cat(paste0("\n", Sys.time()," ------ Iteration ", k, "/100 - Start ------\n"))
  
  ### 2.1.1/ Randomly choose one submodel per modeled OMU
  
  Random_sub_stack_Jaccard.80 <- Random_sub_stack_Jaccard.95 <- Random_sub_stack_TSS.80 <- Random_sub_stack_TSS.95 <- stack(continent_mask) # 1e temp layer used to initiate stack, to remove afterwards
  # i <- 1
  
  for (i in 1:nrow(modeled_OMU)) # By modeled OMU
  { 
    # Load OMU name
    unit <- as.character(modeled_OMU$Tag.model[i]) 
    
    # Load stack of all submodels
    Select_Jaccard_stack_cont_cropped_80 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Select_Jaccard_stack_cont_cropped_80.rds"))
    Select_Jaccard_stack_cont_cropped_95 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Select_Jaccard_stack_cont_cropped_95.rds"))
    
    Select_TSS_stack_cont_cropped_80 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Select_TSS_stack_cont_cropped_80.rds"))
    Select_TSS_stack_cont_cropped_95 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Select_TSS_stack_cont_cropped_95.rds"))
    
    # plot(Select_Jaccard_stack_cont_cropped_80)
    
    # Randomly select one submodel
    N_sub_Jaccard <- nlayers(Select_Jaccard_stack_cont_cropped_80)
    ID_sub_Jaccard <- sample(1:N_sub_Jaccard, 1) 
    
    # Randomly select one submodel
    N_sub_TSS <- nlayers(Select_TSS_stack_cont_cropped_80)
    ID_sub_TSS <- sample(1:N_sub_TSS, 1) 
    
    # Extract the submodel
    Random_sub_Jaccard.80 <- subset(Select_Jaccard_stack_cont_cropped_80, subset = ID_sub_Jaccard)
    Random_sub_Jaccard.95 <- subset(Select_Jaccard_stack_cont_cropped_95, subset = ID_sub_Jaccard)
    Random_sub_TSS.80 <- subset(Select_TSS_stack_cont_cropped_80, subset = ID_sub_TSS)
    Random_sub_TSS.95 <- subset(Select_TSS_stack_cont_cropped_95, subset = ID_sub_TSS)
    
    # Build stack
    Random_sub_stack_Jaccard.80 <- addLayer(Random_sub_stack_Jaccard.80, Random_sub_Jaccard.80)
    Random_sub_stack_Jaccard.95 <- addLayer(Random_sub_stack_Jaccard.95, Random_sub_Jaccard.95)
    Random_sub_stack_TSS.80 <- addLayer(Random_sub_stack_TSS.80, Random_sub_TSS.80)
    Random_sub_stack_TSS.95 <- addLayer(Random_sub_stack_TSS.95, Random_sub_TSS.95)
    
    if (i %% 10 == 0) {cat(paste0(i," / 563 Modeled OMUs\n"))}
  }
  
  # Drop the useless first layer used to initiate the stack
  Random_sub_stack_Jaccard.80 <- dropLayer(Random_sub_stack_Jaccard.80, i = 1)
  Random_sub_stack_Jaccard.95 <- dropLayer(Random_sub_stack_Jaccard.95, i = 1)
  Random_sub_stack_TSS.80 <- dropLayer(Random_sub_stack_TSS.80, i = 1)
  Random_sub_stack_TSS.95 <- dropLayer(Random_sub_stack_TSS.95, i = 1)
  
  # Name layers with OMU names
  names(Random_sub_stack_Jaccard.80) <- names(Random_sub_stack_Jaccard.95) <- names(Random_sub_stack_TSS.80) <- names(Random_sub_stack_TSS.95) <-as.character(modeled_OMU$Tag.model)
  
  # nlayers(Random_sub_stack_Jaccard.80) # 563 modeled OMU in the final stacks
  
  ### 2.1.2/ Merge at species level and stack the proba maps
  
  # Create map at species level by computing "probability" of presence of any of the OMU of the species
  # Generate a continous map for each 4 options between Jaccard/TSS and Buffer.80/95
  
  
  All_sp_proba_stack_Jaccard.80 <- All_sp_proba_stack_Jaccard.95 <- All_sp_proba_stack_TSS.80 <- All_sp_proba_stack_TSS.95 <- stack(continent_mask) 
  ### Loop by species
  for (i in 1:nrow(list.sp)) 
  {
    # i <- 23
    
    # Load a random layer to initiate the stacks of OMU per species (to remove later)
    sp.stack_Jaccard.80 <- sp.stack_Jaccard.95 <- sp.stack_TSS.80 <- sp.stack_TSS.95 <- stack(continent_mask) 
    
    # Load Species name
    sp <- as.character(list.sp$Sp_full[i]) 
    
    # Modeled OMU
    OMU_list <- modeled_OMU[modeled_OMU$Sp_full == sp, "Tag.model"]
    if (length(OMU_list) > 0)  # If at least one OMU was modeled
    {      
      sp.stack_Jaccard.80 <- stack(sp.stack_Jaccard.80, Random_sub_stack_Jaccard.80[[OMU_list]])
      sp.stack_Jaccard.95 <- stack(sp.stack_Jaccard.95, Random_sub_stack_Jaccard.95[[OMU_list]])
      
      sp.stack_TSS.80 <- stack(sp.stack_TSS.80, Random_sub_stack_TSS.80[[OMU_list]])
      sp.stack_TSS.95 <- stack(sp.stack_TSS.95, Random_sub_stack_TSS.95[[OMU_list]])
    }  
    
    # Rasterized OMU  
    Binaries <- rasterized_OMU[rasterized_OMU$Sp_full == sp,]
    if (nrow(Binaries) > 0)  # If at least one OMU was rasterized
    {  
      for (j in 1:nrow(Binaries))  # For each  OMU
      {  
        unit <-  as.character(Binaries$Tag.model[j]) # Load the unit name 
        
        # Load the rasterized map, with continental borders
        rasterized_map <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/rasterized_map_",unit,".rds"))
        
        # Jaccard.80: Load the binary maps and stack them all
        sp.stack_Jaccard.80 <- addLayer(sp.stack_Jaccard.80, rasterized_map)
        # Jaccard.95: Load the binary maps and stack them all
        sp.stack_Jaccard.95 <- addLayer(sp.stack_Jaccard.95, rasterized_map)
        # TSS.80: Load the binary maps and stack them all
        sp.stack_TSS.80 <- addLayer(sp.stack_TSS.80, rasterized_map)   
        # TSS.95: Load the binary maps and stack them all
        sp.stack_TSS.95 <- addLayer(sp.stack_TSS.95, rasterized_map)
      }
    }
    
    # Drop the useless first layer used to initiate the stack
    sp.stack_Jaccard.80 <- dropLayer(sp.stack_Jaccard.80, i = 1)
    sp.stack_Jaccard.95 <- dropLayer(sp.stack_Jaccard.95, i = 1)
    sp.stack_TSS.80 <- dropLayer(sp.stack_TSS.80, i = 1)
    sp.stack_TSS.95 <- dropLayer(sp.stack_TSS.95, i = 1)
    
    # Compute probability map per Species
    sp.cont_Jaccard.80 <- calc(sp.stack_Jaccard.80, fun = aggreg_prob)
    sp.cont_Jaccard.95 <- calc(sp.stack_Jaccard.95, fun = aggreg_prob)
    sp.cont_TSS.80 <- calc(sp.stack_TSS.80, fun = aggreg_prob)
    sp.cont_TSS.95 <- calc(sp.stack_TSS.95, fun = aggreg_prob)
    
    # plot(sp.cont_Jaccard.80)
    
    # Add to the final species stack
    All_sp_proba_stack_Jaccard.80 <- stack(All_sp_proba_stack_Jaccard.80, sp.cont_Jaccard.80)
    All_sp_proba_stack_Jaccard.95 <- stack(All_sp_proba_stack_Jaccard.95, sp.cont_Jaccard.95)
    
    All_sp_proba_stack_TSS.80 <- stack(All_sp_proba_stack_TSS.80, sp.cont_TSS.80)
    All_sp_proba_stack_TSS.95 <- stack(All_sp_proba_stack_TSS.95, sp.cont_TSS.95)
    
    if (i %% 10 == 0) {cat(paste0(i," / 388 species\n"))}
    
  }
  
  # Drop the useless first layer used to initiate the stack
  All_sp_proba_stack_Jaccard.80 <- dropLayer(All_sp_proba_stack_Jaccard.80, i = 1)
  All_sp_proba_stack_Jaccard.95 <- dropLayer(All_sp_proba_stack_Jaccard.95, i = 1)
  All_sp_proba_stack_TSS.80 <- dropLayer(All_sp_proba_stack_TSS.80, i = 1)
  All_sp_proba_stack_TSS.95 <- dropLayer(All_sp_proba_stack_TSS.95, i = 1)
  
  # Renames layer with sp names
  names(All_sp_proba_stack_Jaccard.80) <- names(All_sp_proba_stack_Jaccard.95) <- names(All_sp_proba_stack_TSS.80) <- names(All_sp_proba_stack_TSS.95) <- as.character(list.sp$Sp_full) 
  
  ### 2.1.3/ Compute richness map
  tot.sp.richness_Jaccard.80 <- readAll(calc(All_sp_proba_stack_Jaccard.80, fun = sum))
  tot.sp.richness_Jaccard.95 <- readAll(calc(All_sp_proba_stack_Jaccard.95, fun = sum))
  tot.sp.richness_TSS.80 <- readAll(calc(All_sp_proba_stack_TSS.80, fun = sum))
  tot.sp.richness_TSS.95 <- readAll(calc(All_sp_proba_stack_TSS.95, fun = sum))
  
  # Save randomly drawn richness maps for this iteration
  saveRDS(tot.sp.richness_Jaccard.80, file = paste0("./outputs/Uncertainties_maps/richness_maps_iterations/tot.sp.richness_Jaccard.80_",k,".rds"), version = "2")
  saveRDS(tot.sp.richness_Jaccard.95, file = paste0("./outputs/Uncertainties_maps/richness_maps_iterations/tot.sp.richness_Jaccard.95_",k,".rds"), version = "2")
  saveRDS(tot.sp.richness_TSS.80, file = paste0("./outputs/Uncertainties_maps/richness_maps_iterations/tot.sp.richness_TSS.80_",k,".rds"), version = "2")
  saveRDS(tot.sp.richness_TSS.95, file = paste0("./outputs/Uncertainties_maps/richness_maps_iterations/tot.sp.richness_TSS.95_",k,".rds"), version = "2")
  
  ### 2.1.4/ Stack randomly drawn richness maps
  Random_richness_stack_Jaccard.80 <- addLayer(Random_richness_stack_Jaccard.80, tot.sp.richness_Jaccard.80)
  Random_richness_stack_Jaccard.95 <- addLayer(Random_richness_stack_Jaccard.95, tot.sp.richness_Jaccard.95)
  Random_richness_stack_TSS.80 <- addLayer(Random_richness_stack_TSS.80, tot.sp.richness_TSS.80)
  Random_richness_stack_TSS.95 <- addLayer(Random_richness_stack_TSS.95, tot.sp.richness_TSS.95)
  
  # Clean environnement from stacks and maps, except the final richness stack
  rm(list = setdiff(ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")], str_subset(string = ls(), pattern = "(Random_richness_stack)")))
  
  cat(paste0("\n", Sys.time()," ------ Iteration ", k, "/100 - Done ------\n"))
  
}

### To make stack from saved layers (when using parallelized computation)
Random_richness_stack_Jaccard.80 <- Random_richness_stack_Jaccard.95 <- Random_richness_stack_TSS.80 <- Random_richness_stack_TSS.95 <- stack(continent_mask) # 1e temp layer used to initiate stack, to remove afterwards
for (i in 1:100) {
# for (i in 26:100) {
  Random_richness_stack_Jaccard.80 <- addLayer(Random_richness_stack_Jaccard.80, readRDS(file = paste0("./outputs/Uncertainties_maps/richness_maps_iterations/tot.sp.richness_Jaccard.80_",i,".rds")))
  Random_richness_stack_Jaccard.95 <- addLayer(Random_richness_stack_Jaccard.95, readRDS(file = paste0("./outputs/Uncertainties_maps/richness_maps_iterations/tot.sp.richness_Jaccard.95_",i,".rds")))
  Random_richness_stack_TSS.80 <- addLayer(Random_richness_stack_TSS.80, readRDS(file = paste0("./outputs/Uncertainties_maps/richness_maps_iterations/tot.sp.richness_TSS.80_",i,".rds")))
  Random_richness_stack_TSS.95 <- addLayer(Random_richness_stack_TSS.95, readRDS(file = paste0("./outputs/Uncertainties_maps/richness_maps_iterations/tot.sp.richness_TSS.95_",i,".rds")))
}
###

# Drop the useless first layer used to initiate the stack
Random_richness_stack_Jaccard.80 <- dropLayer(Random_richness_stack_Jaccard.80, i = 1)
Random_richness_stack_Jaccard.95 <- dropLayer(Random_richness_stack_Jaccard.95, i = 1)
Random_richness_stack_TSS.80 <- dropLayer(Random_richness_stack_TSS.80, i = 1)
Random_richness_stack_TSS.95 <- dropLayer(Random_richness_stack_TSS.95, i = 1)

# Save randomly drawn richness stacks
saveRDS(Random_richness_stack_Jaccard.80, file = paste0("./outputs/Uncertainties_maps/Random_richness_stack_Jaccard.80.rds"), version = "2")
saveRDS(Random_richness_stack_Jaccard.95, file = paste0("./outputs/Uncertainties_maps/Random_richness_stack_Jaccard.95.rds"), version = "2")
saveRDS(Random_richness_stack_TSS.80, file = paste0("./outputs/Uncertainties_maps/Random_richness_stack_TSS.80.rds"), version = "2")
saveRDS(Random_richness_stack_TSS.95, file = paste0("./outputs/Uncertainties_maps/Random_richness_stack_TSS.95.rds"), version = "2")





### 2.2/ Compute richness uncertainties based on randomly drawn model outputs ####

### Load the final Species richness layer
tot.sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
tot.sp.richness_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.95.rds"))
tot.sp.richness_TSS.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_TSS.80.rds"))
tot.sp.richness_TSS.95 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_TSS.95.rds"))

# Load randomly drawn richness stacks
Random_richness_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Uncertainties_maps/Random_richness_stack_Jaccard.80.rds"))
Random_richness_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Uncertainties_maps/Random_richness_stack_Jaccard.95.rds"))
Random_richness_stack_TSS.80 <- readRDS(file = paste0("./outputs/Uncertainties_maps/Random_richness_stack_TSS.80.rds"))
Random_richness_stack_TSS.95 <- readRDS(file = paste0("./outputs/Uncertainties_maps/Random_richness_stack_TSS.95.rds"))

# Compute Uncertainties as sd
Random_richness_sd_Jaccard.80 <- readAll(calc(Random_richness_stack_Jaccard.80, fun = sd, na.rm = T))
Random_richness_sd_Jaccard.95 <- readAll(calc(Random_richness_stack_Jaccard.95, fun = sd, na.rm = T))
Random_richness_sd_TSS.80 <- readAll(calc(Random_richness_stack_TSS.80, fun = sd, na.rm = T))
Random_richness_sd_TSS.95 <- readAll(calc(Random_richness_stack_TSS.95, fun = sd, na.rm = T))

# Compute Uncertainties as sd from the richness obtained from the Ensembles
Random_richness_var_to_EM_Jaccard.80 <- calc((Random_richness_stack_Jaccard.80 - tot.sp.richness_Jaccard.80) * (Random_richness_stack_Jaccard.80 - tot.sp.richness_Jaccard.80), fun = sum, na.rm = T)/nlayers(Random_richness_stack_Jaccard.80) # Variance = mean of squared differences to the "mean"
Random_richness_sd_to_EM_Jaccard.80 <- sqrt(Random_richness_var_to_EM_Jaccard.80) # SD = sqrt(var)

Random_richness_var_to_EM_Jaccard.95 <- calc((Random_richness_stack_Jaccard.95 - tot.sp.richness_Jaccard.95) * (Random_richness_stack_Jaccard.95 - tot.sp.richness_Jaccard.95), fun = sum, na.rm = T)/nlayers(Random_richness_stack_Jaccard.95) # Variance = mean of squared differences to the "mean"
Random_richness_sd_to_EM_Jaccard.95 <- sqrt(Random_richness_var_to_EM_Jaccard.95) # SD = sqrt(var)

Random_richness_var_to_EM_TSS.80 <- calc((Random_richness_stack_TSS.80 - tot.sp.richness_TSS.80) * (Random_richness_stack_TSS.80 - tot.sp.richness_TSS.80), fun = sum, na.rm = T)/nlayers(Random_richness_stack_TSS.80) # Variance = mean of squared differences to the "mean"
Random_richness_sd_to_EM_TSS.80 <- sqrt(Random_richness_var_to_EM_TSS.80) # SD = sqrt(var)

Random_richness_var_to_EM_TSS.95 <- calc((Random_richness_stack_TSS.95 - tot.sp.richness_TSS.95) * (Random_richness_stack_TSS.95 - tot.sp.richness_TSS.95), fun = sum, na.rm = T)/nlayers(Random_richness_stack_TSS.95) # Variance = mean of squared differences to the "mean"
Random_richness_sd_to_EM_TSS.95 <- sqrt(Random_richness_var_to_EM_TSS.95) # SD = sqrt(var)

# Save Uncertainties
save(Random_richness_sd_Jaccard.80, file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_Jaccard.80.RData"), version = "2")
saveRDS(Random_richness_sd_Jaccard.80, file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_Jaccard.80.rds"), version = "2")
save(Random_richness_sd_Jaccard.95, file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_Jaccard.95.RData"), version = "2")
saveRDS(Random_richness_sd_Jaccard.95, file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_Jaccard.95.rds"), version = "2")
save(Random_richness_sd_TSS.80, file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_TSS.80.RData"), version = "2")
saveRDS(Random_richness_sd_TSS.80, file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_TSS.80.rds"), version = "2")
save(Random_richness_sd_TSS.95, file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_TSS.95.RData"), version = "2")
saveRDS(Random_richness_sd_TSS.95, file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_TSS.95.rds"), version = "2")

save(Random_richness_sd_to_EM_Jaccard.80, file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_to_EM_Jaccard.80.RData"), version = "2")
saveRDS(Random_richness_sd_to_EM_Jaccard.80, file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_to_EM_Jaccard.80.rds"), version = "2")
save(Random_richness_sd_to_EM_Jaccard.95, file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_to_EM_Jaccard.95.RData"), version = "2")
saveRDS(Random_richness_sd_to_EM_Jaccard.95, file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_to_EM_Jaccard.95.rds"), version = "2")
save(Random_richness_sd_to_EM_TSS.80, file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_to_EM_TSS.80.RData"), version = "2")
saveRDS(Random_richness_sd_to_EM_TSS.80, file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_to_EM_TSS.80.rds"), version = "2")
save(Random_richness_sd_to_EM_TSS.95, file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_to_EM_TSS.95.RData"), version = "2")
saveRDS(Random_richness_sd_to_EM_TSS.95, file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_to_EM_TSS.95.rds"), version = "2")


### 2.3/ Plot richness uncertainty based on randomly drawn model outputs ####

## 2.3.1/ SD compared to the mean of randomly drawn richness ####

# Load directly the richness uncertainty based on randomly drawn model outputs
Random_richness_sd_Jaccard.80 <- readRDS(file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_Jaccard.80.rds"))
Random_richness_sd_Jaccard.95 <- readRDS(file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_Jaccard.95.rds"))
Random_richness_sd_TSS.80 <- readRDS(file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_TSS.80.rds"))
Random_richness_sd_TSS.95 <- readRDS(file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_TSS.95.rds"))


### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Uncertainties_maps/Random_richness_sd_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Random_richness_sd_Jaccard.80, col = pal_grn_red, main = paste0("Richness Uncertainty \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Uncertainty\n            (sd)", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(Random_richness_sd_Jaccard.80, locs = seq(0, 3, 0.5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -106, y = 4, font = 2, cex = 1.2, label = "Uncertainty (sd)")
par(mar = internal_margins)
dev.off()

# Jaccard.95
pdf(file = paste0("./maps/Uncertainties_maps/Random_richness_sd_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Random_richness_sd_Jaccard.95, col = pal_grn_red, main = paste0("Richness Uncertainty \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Uncertainty\n            (sd)", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(Random_richness_sd_Jaccard.95, locs = seq(0, 3, 0.5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -106, y = 4, font = 2, cex = 1.2, label = "Uncertainty (sd)")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Uncertainties_maps/Random_richness_sd_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Random_richness_sd_TSS.80, col = pal_grn_red, main = paste0("Richness Uncertainty \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Uncertainty\n            (sd)", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(Random_richness_sd_TSS.80, locs = seq(0, 3, 0.5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -106, y = 4, font = 2, cex = 1.2, label = "Uncertainty (sd)")
par(mar = internal_margins)
dev.off()

str(Random_richness_sd_TSS.95)

# TSS.95
pdf(file = paste0("./maps/Uncertainties_maps/Random_richness_sd_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Random_richness_sd_TSS.95, col = pal_grn_red, main = paste0("Richness Uncertainty \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Uncertainty\n            (sd)", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(Random_richness_sd_TSS.95, locs = seq(0, 3, 0.5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -106, y = 4, font = 2, cex = 1.2, label = "Uncertainty (sd)")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(Random_richness_sd_Jaccard.80) +
#   tm_raster(palette = pal_grn_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Uncertainties_maps/Random_richness_sd_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(Random_richness_sd_Jaccard.80, col = pal_grn_red, main = paste0("Richness Uncertainty \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Uncertainty\n            (sd)", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(Random_richness_sd_Jaccard.80, locs = seq(0, 3, 0.5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -106, y = 4, font = 2, cex = 1.2, label = "Uncertainty (sd)")

image(Random_richness_sd_Jaccard.95, col = pal_grn_red, main = paste0("Richness Uncertainty \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Uncertainty\n            (sd)", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(Random_richness_sd_Jaccard.95, locs = seq(0, 3, 0.5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -106, y = 4, font = 2, cex = 1.2, label = "Uncertainty (sd)")

image(Random_richness_sd_TSS.80, col = pal_grn_red, main = paste0("Richness Uncertainty \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Uncertainty\n            (sd)", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(Random_richness_sd_TSS.80, locs = seq(0, 3, 0.5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -106, y = 4, font = 2, cex = 1.2, label = "Uncertainty (sd)")

image(Random_richness_sd_TSS.95, col = pal_grn_red, main = paste0("Richness Uncertainty \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Uncertainty\n            (sd)", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(Random_richness_sd_TSS.95, locs = seq(0, 3, 0.5), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -106, y = 4, font = 2, cex = 1.2, label = "Uncertainty (sd)")

par(mar = internal_margins)

dev.off()


# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])


## 2.3.2/ SD compared to the richness from the Ensembles ####

# Load directly the richness uncertainty based on randomly drawn model outputs
Random_richness_sd_to_EM_Jaccard.80 <- readRDS(file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_to_EM_Jaccard.80.rds"))
Random_richness_sd_to_EM_Jaccard.95 <- readRDS(file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_to_EM_Jaccard.95.rds"))
Random_richness_sd_to_EM_TSS.80 <- readRDS(file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_to_EM_TSS.80.rds"))
Random_richness_sd_to_EM_TSS.95 <- readRDS(file = paste0("./outputs/Uncertainties_maps/Random_richness_sd_to_EM_TSS.95.rds"))


### Individual plots

# Jaccard.80
pdf(file = paste0("./maps/Uncertainties_maps/Random_richness_sd_to_EM_Jaccard.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Random_richness_sd_to_EM_Jaccard.80, col = pal_grn_red, main = paste0("Richness Uncertainty \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Uncertainty\n            (sd)", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(Random_richness_sd_to_EM_Jaccard.80, locs = seq(0, 6, 1), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -106, y = 4, font = 2, cex = 1.2, label = "Uncertainty (sd)")
par(mar = internal_margins)
dev.off()


# Jaccard.95
pdf(file = paste0("./maps/Uncertainties_maps/Random_richness_sd_to_EM_Jaccard.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Random_richness_sd_to_EM_Jaccard.95, col = pal_grn_red, main = paste0("Richness Uncertainty \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Uncertainty\n            (sd)", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(Random_richness_sd_to_EM_Jaccard.95, locs = seq(0, 0.09, 0.02), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -106, y = 4, font = 2, cex = 1.2, label = "Uncertainty (sd)")
par(mar = internal_margins)
dev.off()

# TSS.80
pdf(file = paste0("./maps/Uncertainties_maps/Random_richness_sd_to_EM_TSS.80.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Random_richness_sd_to_EM_TSS.80, col = pal_grn_red, main = paste0("Richness Uncertainty \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Uncertainty\n            (sd)", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(Random_richness_sd_to_EM_TSS.80, locs = seq(0, 0.08, 0.02), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -106, y = 4, font = 2, cex = 1.2, label = "Uncertainty (sd)")
par(mar = internal_margins)
dev.off()


# TSS.95
pdf(file = paste0("./maps/Uncertainties_maps/Random_richness_sd_to_EM_TSS.95.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(Random_richness_sd_to_EM_TSS.95, col = pal_grn_red, main = paste0("Richness Uncertainty \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Uncertainty\n            (sd)", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(Random_richness_sd_to_EM_TSS.95, locs = seq(0, 0.09, 0.02), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -106, y = 4, font = 2, cex = 1.2, label = "Uncertainty (sd)")
par(mar = internal_margins)
dev.off()


### Tmap version
#
# library(tmap)
# 
# tm_shape(Random_richness_sd_to_EM_Jaccard.80) +
#   tm_raster(palette = pal_grn_red) +
#   tm_shape(crop_mask_shp) +
#   tm_borders(lwd = 1.2, col = "grey20")
# 
# tmap_mode("view")

### Multiple pages pdf

pdf(file = paste0("./maps/Uncertainties_maps/Random_richness_sd_to_EM_all_maps.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

image(Random_richness_sd_to_EM_Jaccard.80, col = pal_grn_red, main = paste0("Richness Uncertainty \nJaccard.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Uncertainty\n            (sd)", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(Random_richness_sd_to_EM_Jaccard.80, locs = seq(0, 0.08, 0.02), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -106, y = 4, font = 2, cex = 1.2, label = "Uncertainty (sd)")

image(Random_richness_sd_to_EM_Jaccard.95, col = pal_grn_red, main = paste0("Richness Uncertainty \nJaccard.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Uncertainty\n            (sd)", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(Random_richness_sd_to_EM_Jaccard.95, locs = seq(0, 0.09, 0.02), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -106, y = 4, font = 2, cex = 1.2, label = "Uncertainty (sd)")

image(Random_richness_sd_to_EM_TSS.80, col = pal_grn_red, main = paste0("Richness Uncertainty \nTSS.80"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Uncertainty\n            (sd)", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(Random_richness_sd_to_EM_TSS.80, locs = seq(0, 0.08, 0.02), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -106, y = 4, font = 2, cex = 1.2, label = "Uncertainty (sd)")

image(Random_richness_sd_to_EM_TSS.95, col = pal_grn_red, main = paste0("Richness Uncertainty \nTSS.95"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="          Uncertainty\n            (sd)", cex=1.2, line = 1, font = 2), 
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
addRasterLegend(Random_richness_sd_to_EM_TSS.95, locs = seq(0, 0.09, 0.02), cex.axis = 1.2, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -106, y = 4, font = 2, cex = 1.2, label = "Uncertainty (sd)")

par(mar = internal_margins)

dev.off()


# Clean environnement from stacks and maps
rm(list = ls()[grep(x = ls(), pattern = "(Jaccard)|(TSS)")])
