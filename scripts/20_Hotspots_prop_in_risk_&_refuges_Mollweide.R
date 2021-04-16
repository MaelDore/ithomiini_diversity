
###### Script 20: Categorical threats maps within hotspots ######

# Only for Jaccard.80

# Map risk and refuges areas from Human Footprint within biodiversity hotspots (Figure 5)
# Compute proportions of overlaps between categories (Figure 6)

### Use only 4 indices

# 1/ Species richness
# 2/ Mean Continuous Species Rarity = Range-size Weighted Species Richness
# 3/ Mimicry richness
# 4/ Mean Continuous Mimicry rarity = Range-size Weighted Mimicry Richness
# Bonus/ Mean pairwise Phylogenetic Distance (removed from the final analyses)

### Two levels of hotspots: 25% and 5%

### Two types of bins for Human Footprint: quantiles and absolute values

###

# Inputs: 
    # Indices maps from script 14a
    # Threat map: HII (HF, Venter et al., 2016)
    # Raster masks and shp for hotspots from script 18
    # Raster mask for HF categories from script 19

# Outputs:
   # Map of overlapping regions
      # Per quantiles 5 & 25% (Figure 5) and per absolute values (1, 5 vs. 11, 23)
      # For hotspots = 5% & 25%
   # Bar plots of proportion of hotspots per HII bins
      # Per quantiles 5 & 25% (Figure 6) and per absolute values (1, 5 vs. 11, 23)
      # For hotspots = 5% & 25% 

###


### 1/ Load stuff ####

# Clean environment
rm(list = ls())

library(raster)
library(prettymapr)
library(rangeBuilder)

# New color palette
pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")

# Load mask for continent borders, plot border, and grid
grid_Mollweide_out <- readRDS(file = "./input_data/Map_stuff/grid_Mollweide_out.rds")
large_bg_mask_Mollweide <- readRDS(file = "./input_data/Map_stuff/large_bg_mask_Mollweide.rds")
bbox_sp_Mollweide <- readRDS(file = "./input_data/Map_stuff/bbox_sp_Mollweide.rds")

# Load masks for HII quantiles
load(file = "./outputs/Threat_maps/HII_cropped/HII_masks_quantiles_files.rds") 

# Load masks for HII bins based on absolute values
load(file = "./outputs/Threat_maps/HII_cropped/HII_masks_bins_files.rds")

# Load masks for index hotspots
load(file = "./outputs/Threat_maps/Cropped/sp.richness_masks_rank_files.RData") 
load(file = "./outputs/Threat_maps/Cropped/sp.rarity_masks_rank_files.RData") 
load(file = "./outputs/Threat_maps/Cropped/MPD_masks_rank_files.RData")
load(file = "./outputs/Threat_maps/Cropped/ring.richness_masks_rank_files.RData")
load(file = "./outputs/Threat_maps/Cropped/ring.rarity_masks_rank_files.RData")

# Load shp files for index hotspots
load(file = "./outputs/Threat_maps/Cropped/sp.richness_masked_rank_shp_files.rds") 
load(file = "./outputs/Threat_maps/Cropped/sp.rarity_masked_rank_shp_files.rds") 
load(file = "./outputs/Threat_maps/Cropped/MPD_masked_rank_shp_files.rds") 
load(file = "./outputs/Threat_maps/Cropped/ring.richness_masked_rank_shp_files.rds") 
load(file = "./outputs/Threat_maps/Cropped/ring.rarity_masked_rank_shp_files.rds") 

bg_mask <- readRDS(file = "./input_data/Map_stuff/bg_mask.rds")
crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_15.rds"))


### 2/ Compute categorical maps of HII ####

# 2.1/ For quantiles

plot(HII_mask_rank_5)  # red4 = 1
plot(HII_mask_rank_25) # red = 2
plot(HII_mask_rank_75) # dodgerblue = 3
plot(HII_mask_rank_95) # navyblue = 4

# Compute categorical maps of HII bins
HII_cat_map_rank <- HII_mask_rank_5*(-1) + HII_mask_rank_25*2 + HII_mask_rank_75*3 + HII_mask_rank_95
plot(HII_cat_map_rank)
saveRDS(HII_cat_map_rank, file = "./outputs/Threat_maps/Overlaps/HII_cat_map_rank")

# 2.2/ For absolute values
plot(HII_mask_very_high)  # red4 = 1
plot(HII_mask_high) # red = 2
plot(HII_mask_low) # dodgerblue = 3
plot(HII_mask_very_low) # navyblue = 4

# Compute categorical maps of HII bins
HII_cat_map_bins <- HII_mask_very_high*(-1) + HII_mask_high*2 + HII_mask_low*3 + HII_mask_very_low
plot(HII_cat_map_bins)
saveRDS(HII_cat_map_bins, file = "./outputs/Threat_maps/Overlaps/HII_cat_map_bins")

# Only very low category
HII_cat_map_bins_3cat <- HII_mask_very_high*(-1) + HII_mask_high*2 + HII_mask_very_low*3
plot(HII_cat_map_bins_3cat)
saveRDS(HII_cat_map_bins, file = "./outputs/Threat_maps/Overlaps/HII_cat_map_bins_3cat")


##### 3/ Map overlaps ####

### Categories :

# 0 = No hotspot
# 1 = No refuge or risk area
# 2 = High risk area (5%)
# 3 = Risk area (25%)
# 4 = Refuge (25%)
# 5 = Top refuge (5%)

### 3.1/ Species richness ####

# 3.1.1/ For quantiles ####

plot(sp.richness_mask_rank_25)
plot(sp.richness_mask_rank_5)

# Crop with hotspots = 25% highest species richness
sp.richness_cat_map_rank_hotspot_25 <- HII_cat_map_rank*sp.richness_mask_rank_25 # Crop for hotspot only
sp.richness_cat_map_rank_hotspot_25 <- sp.richness_cat_map_rank_hotspot_25 + sp.richness_mask_rank_25 # Add all hotspots
plot(sp.richness_cat_map_rank_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"))
plot(sp.richness_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)
saveRDS(sp.richness_cat_map_rank_hotspot_25, file = "./outputs/Threat_maps/Overlaps/sp.richness_cat_map_rank_hotspot_25")

# Crop with hotspots = 5% highest species richness
sp.richness_cat_map_rank_hotspot_5 <- HII_cat_map_rank*sp.richness_mask_rank_5
sp.richness_cat_map_rank_hotspot_5 <- sp.richness_cat_map_rank_hotspot_5 + sp.richness_mask_rank_5 # Add all hotspots
plot(sp.richness_cat_map_rank_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"))
plot(sp.richness_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)
saveRDS(sp.richness_cat_map_rank_hotspot_5, file = "./outputs/Threat_maps/Overlaps/sp.richness_cat_map_rank_hotspot_5")


pdf(file = paste0("./maps/Threat_maps/Overlaps/sp.richness_cat_map_rank_hotspot_25.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(sp.richness_cat_map_rank_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), main = paste0("Species richness hotspots (25%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.richness_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()


pdf(file = paste0("./maps/Threat_maps/Overlaps/sp.richness_cat_map_rank_hotspot_5.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(sp.richness_cat_map_rank_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), main = paste0("Species richness hotspots (5%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.richness_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()

# 3.1.2/ For absolute values ####

plot(sp.richness_mask_rank_25)
plot(sp.richness_mask_rank_5)

# Crop with hotspots = 25% highest species richness
sp.richness_cat_map_bins_3cat_hotspot_25 <- HII_cat_map_bins_3cat*sp.richness_mask_rank_25
sp.richness_cat_map_bins_3cat_hotspot_25 <- sp.richness_cat_map_bins_3cat_hotspot_25 + sp.richness_mask_rank_25 # Add all hotspots
plot(sp.richness_cat_map_bins_3cat_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"))
saveRDS(sp.richness_cat_map_bins_3cat_hotspot_25, file = "./outputs/Threat_maps/Overlaps/sp.richness_cat_map_bins_hotspot_25")

# Crop with hotspots = 5% highest species richness
sp.richness_cat_map_bins_3cat_hotspot_5 <- HII_cat_map_bins_3cat*sp.richness_mask_rank_5
sp.richness_cat_map_bins_3cat_hotspot_5 <- sp.richness_cat_map_bins_3cat_hotspot_5 + sp.richness_mask_rank_5 # Add all hotspots
plot(sp.richness_cat_map_bins_3cat_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"))
saveRDS(sp.richness_cat_map_bins_3cat_hotspot_5, file = "./outputs/Threat_maps/Overlaps/sp.richness_cat_map_bins_hotspot_5")


pdf(file = paste0("./maps/Threat_maps/Overlaps/sp.richness_cat_map_bins_3cat_hotspot_25.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(sp.richness_cat_map_bins_3cat_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Species richness hotspots (25%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.richness_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()


pdf(file = paste0("./maps/Threat_maps/Overlaps/sp.richness_cat_map_bins_3cat_hotspot_5.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(sp.richness_cat_map_bins_3cat_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Species richness hotspots (5%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.richness_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()

# 3.1.3/ All 4 options ####


pdf(file = paste0("./maps/Threat_maps/Overlaps/sp.richness_cat_map_all_options.pdf"), height = 10, width = 12)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
par(mfrow = c(2,2))


# A/ Hotspots 25%, for quantiles

image(sp.richness_cat_map_rank_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), main = paste0("Species richness hotspots (25%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.richness_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

# B/ Hotspots 25%, for absolute values

image(sp.richness_cat_map_rank_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), main = paste0("Species richness hotspots (5%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.richness_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

# C/ Hotspots 5%, for quantiles

image(sp.richness_cat_map_bins_3cat_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Species richness hotspots (25%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.richness_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

# D/ Hotspots 5%, for absolute values

image(sp.richness_cat_map_bins_3cat_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Species richness hotspots (5%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.richness_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
par(mfrow = c(1, 1))

dev.off()



### 3.2/ Species rarity ####

# 3.2.1/ For quantiles ####

plot(sp.rarity_mask_rank_25)
plot(sp.rarity_mask_rank_5)

# Crop with hotspots = 25% highest Mean species rarity
sp.rarity_cat_map_rank_hotspot_25 <- HII_cat_map_rank*sp.rarity_mask_rank_25 # Crop for hotspot only
sp.rarity_cat_map_rank_hotspot_25 <- sp.rarity_cat_map_rank_hotspot_25 + sp.rarity_mask_rank_25 # Add all hotspots
plot(sp.rarity_cat_map_rank_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"))
plot(sp.rarity_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)
saveRDS(sp.rarity_cat_map_rank_hotspot_25, file = "./outputs/Threat_maps/Overlaps/sp.rarity_cat_map_rank_hotspot_25")

# Crop with hotspots = 5% highest Mean species rarity
sp.rarity_cat_map_rank_hotspot_5 <- HII_cat_map_rank*sp.rarity_mask_rank_5
sp.rarity_cat_map_rank_hotspot_5 <- sp.rarity_cat_map_rank_hotspot_5 + sp.rarity_mask_rank_5 # Add all hotspots
plot(sp.rarity_cat_map_rank_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"))
plot(sp.rarity_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)
saveRDS(sp.rarity_cat_map_rank_hotspot_5, file = "./outputs/Threat_maps/Overlaps/sp.rarity_cat_map_rank_hotspot_5")


pdf(file = paste0("./maps/Threat_maps/Overlaps/sp.rarity_cat_map_rank_hotspot_25.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(sp.rarity_cat_map_rank_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), main = paste0("Mean species rarity hotspots (25%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.rarity_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()


pdf(file = paste0("./maps/Threat_maps/Overlaps/sp.rarity_cat_map_rank_hotspot_5.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(sp.rarity_cat_map_rank_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Mean species rarity hotspots (5%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.rarity_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()


# 3.2.2/ For absolute values ####

plot(sp.rarity_mask_rank_25)
plot(sp.rarity_mask_rank_5)

# Crop with hotspots = 25% highest Mean species rarity
sp.rarity_cat_map_bins_3cat_hotspot_25 <- HII_cat_map_bins_3cat*sp.rarity_mask_rank_25
sp.rarity_cat_map_bins_3cat_hotspot_25 <- sp.rarity_cat_map_bins_3cat_hotspot_25 + sp.rarity_mask_rank_25 # Add all hotspots
plot(sp.rarity_cat_map_bins_3cat_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"))
saveRDS(sp.rarity_cat_map_bins_3cat_hotspot_25, file = "./outputs/Threat_maps/Overlaps/sp.rarity_cat_map_bins_hotspot_25")

# Crop with hotspots = 5% highest Mean species rarity
sp.rarity_cat_map_bins_3cat_hotspot_5 <- HII_cat_map_bins_3cat*sp.rarity_mask_rank_5
sp.rarity_cat_map_bins_3cat_hotspot_5 <- sp.rarity_cat_map_bins_3cat_hotspot_5 + sp.rarity_mask_rank_5 # Add all hotspots
plot(sp.rarity_cat_map_bins_3cat_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"))
saveRDS(sp.rarity_cat_map_bins_3cat_hotspot_5, file = "./outputs/Threat_maps/Overlaps/sp.rarity_cat_map_bins_hotspot_5")


pdf(file = paste0("./maps/Threat_maps/Overlaps/sp.rarity_cat_map_bins_3cat_hotspot_25.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(sp.rarity_cat_map_bins_3cat_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Mean species rarity hotspots (25%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.rarity_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()


pdf(file = paste0("./maps/Threat_maps/Overlaps/sp.rarity_cat_map_bins_3cat_hotspot_5.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(sp.rarity_cat_map_bins_3cat_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Mean species rarity hotspots (5%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.rarity_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()


# 3.2.3/ All 4 options ####

pdf(file = paste0("./maps/Threat_maps/Overlaps/sp.rarity_cat_map_all_options.pdf"), height = 10, width = 12)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
par(mfrow = c(2,2))


# A/ Hotspots 25%, for quantiles

image(sp.rarity_cat_map_rank_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), main = paste0("Mean species rarity hotspots (25%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.rarity_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

# B/ Hotspots 5%, for quantiles

image(sp.rarity_cat_map_rank_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Mean species rarity hotspots (5%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.rarity_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

# C/ Hotspots 25%, for absolute values

image(sp.rarity_cat_map_bins_3cat_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Mean species rarity hotspots (25%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.rarity_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

# D/ Hotspots 5%, for absolute values

image(sp.rarity_cat_map_bins_3cat_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Mean species rarity hotspots (5%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(sp.rarity_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
par(mfrow = c(1, 1))

dev.off()


### 3.3/ MPD ####

# 3.3.1/ For quantiles ####

plot(MPD_mask_rank_25)
plot(MPD_mask_rank_5)

# Crop with hotspots = 25% highest MPD
MPD_cat_map_rank_hotspot_25 <- HII_cat_map_rank*MPD_mask_rank_25 # Crop for hotspot only
MPD_cat_map_rank_hotspot_25 <- MPD_cat_map_rank_hotspot_25 + MPD_mask_rank_25 # Add all hotspots
plot(MPD_cat_map_rank_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"))
plot(MPD_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)
saveRDS(MPD_cat_map_rank_hotspot_25, file = "./outputs/Threat_maps/Overlaps/MPD_cat_map_rank_hotspot_25")

# Crop with hotspots = 5% highest MPD
MPD_cat_map_rank_hotspot_5 <- HII_cat_map_rank*MPD_mask_rank_5
MPD_cat_map_rank_hotspot_5 <- MPD_cat_map_rank_hotspot_5 + MPD_mask_rank_5 # Add all hotspots
plot(MPD_cat_map_rank_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"))
plot(MPD_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)
saveRDS(MPD_cat_map_rank_hotspot_5, file = "./outputs/Threat_maps/Overlaps/MPD_cat_map_rank_hotspot_5")

pdf(file = paste0("./maps/Threat_maps/Overlaps/MPD_cat_map_rank_hotspot_25.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(MPD_cat_map_rank_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), main = paste0("MPD hotspots (25%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(MPD_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()


pdf(file = paste0("./maps/Threat_maps/Overlaps/MPD_cat_map_rank_hotspot_5.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(MPD_cat_map_rank_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), main = paste0("MPD hotspots (5%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(MPD_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()


# 3.3.2/ For absolute values ####

plot(MPD_mask_rank_25)
plot(MPD_mask_rank_5)

# Crop with hotspots = 25% highest MPD
MPD_cat_map_bins_3cat_hotspot_25 <- HII_cat_map_bins_3cat*MPD_mask_rank_25
MPD_cat_map_bins_3cat_hotspot_25 <- MPD_cat_map_bins_3cat_hotspot_25 + MPD_mask_rank_25 # Add all hotspots
plot(MPD_cat_map_bins_3cat_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"))
saveRDS(MPD_cat_map_bins_3cat_hotspot_25, file = "./outputs/Threat_maps/Overlaps/MPD_cat_map_bins_hotspot_25")


# Crop with hotspots = 5% highest MPD
MPD_cat_map_bins_3cat_hotspot_5 <- HII_cat_map_bins_3cat*MPD_mask_rank_5
MPD_cat_map_bins_3cat_hotspot_5 <- MPD_cat_map_bins_3cat_hotspot_5 + MPD_mask_rank_5 # Add all hotspots
plot(MPD_cat_map_bins_3cat_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"))
saveRDS(MPD_cat_map_bins_3cat_hotspot_5, file = "./outputs/Threat_maps/Overlaps/MPD_cat_map_bins_hotspot_5")


pdf(file = paste0("./maps/Threat_maps/Overlaps/MPD_cat_map_bins_3cat_hotspot_25.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(MPD_cat_map_bins_3cat_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("MPD hotspots (25%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(MPD_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()


pdf(file = paste0("./maps/Threat_maps/Overlaps/MPD_cat_map_bins_3cat_hotspot_5.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(MPD_cat_map_bins_3cat_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("MPD hotspots (5%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(MPD_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()


# 3.3.3/ All 4 options ####

pdf(file = paste0("./maps/Threat_maps/Overlaps/MPD_cat_map_all_options.pdf"), height = 10, width = 12)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
par(mfrow = c(2,2))


# A/ Hotspots 25%, for quantiles

image(MPD_cat_map_rank_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), main = paste0("MPD hotspots (25%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(MPD_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

# B/ Hotspots 5%, for quantiles

image(MPD_cat_map_rank_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), main = paste0("MPD hotspots (5%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(MPD_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

# C/ Hotspots 25%, for absolute values

image(MPD_cat_map_bins_3cat_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("MPD hotspots (25%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(MPD_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

# D/ Hotspots 5%, for absolute values

image(MPD_cat_map_bins_3cat_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("MPD hotspots (5%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(MPD_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
par(mfrow = c(1, 1))

dev.off()


### 4.3/ Mimicry richness ####

# 4.3.1/ For quantiles ####

plot(ring.richness_mask_rank_25)
plot(ring.richness_mask_rank_5)

# Crop with hotspots = 25% highest Mimicry richness
ring.richness_cat_map_rank_hotspot_25 <- HII_cat_map_rank*ring.richness_mask_rank_25 # Crop for hotspot only
ring.richness_cat_map_rank_hotspot_25 <- ring.richness_cat_map_rank_hotspot_25 + ring.richness_mask_rank_25 # Add all hotspots
plot(ring.richness_cat_map_rank_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"))
plot(ring.richness_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)
saveRDS(ring.richness_cat_map_rank_hotspot_25, file = "./outputs/Threat_maps/Overlaps/ring.richness_cat_map_rank_hotspot_25")

# Crop with hotspots = 5% highest Mimicry richness
ring.richness_cat_map_rank_hotspot_5 <- HII_cat_map_rank*ring.richness_mask_rank_5
ring.richness_cat_map_rank_hotspot_5 <- ring.richness_cat_map_rank_hotspot_5 + ring.richness_mask_rank_5 # Add all hotspots
plot(ring.richness_cat_map_rank_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"))
plot(ring.richness_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)
saveRDS(ring.richness_cat_map_rank_hotspot_5, file = "./outputs/Threat_maps/Overlaps/ring.richness_cat_map_rank_hotspot_5")

pdf(file = paste0("./maps/Threat_maps/Overlaps/ring.richness_cat_map_rank_hotspot_25.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(ring.richness_cat_map_rank_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), main = paste0("Mimicry richness hotspots (25%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.richness_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()


pdf(file = paste0("./maps/Threat_maps/Overlaps/ring.richness_cat_map_rank_hotspot_5.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(ring.richness_cat_map_rank_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Mimicry richness hotspots (5%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.richness_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()


# 4.3.2/ For absolute values ####

plot(ring.richness_mask_rank_25)
plot(ring.richness_mask_rank_5)

# Crop with hotspots = 25% highest Mimicry richness
ring.richness_cat_map_bins_3cat_hotspot_25 <- HII_cat_map_bins_3cat*ring.richness_mask_rank_25
ring.richness_cat_map_bins_3cat_hotspot_25 <- ring.richness_cat_map_bins_3cat_hotspot_25 + ring.richness_mask_rank_25 # Add all hotspots
plot(ring.richness_cat_map_bins_3cat_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"))
saveRDS(ring.richness_cat_map_bins_3cat_hotspot_25, file = "./outputs/Threat_maps/Overlaps/ring.richness_cat_map_bins_hotspot_25")

# Crop with hotspots = 5% highest Mimicry richness
ring.richness_cat_map_bins_3cat_hotspot_5 <- HII_cat_map_bins_3cat*ring.richness_mask_rank_5
ring.richness_cat_map_bins_3cat_hotspot_5 <- ring.richness_cat_map_bins_3cat_hotspot_5 + ring.richness_mask_rank_5 # Add all hotspots
plot(ring.richness_cat_map_bins_3cat_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"))
saveRDS(ring.richness_cat_map_bins_3cat_hotspot_5, file = "./outputs/Threat_maps/Overlaps/ring.richness_cat_map_bins_hotspot_5")


pdf(file = paste0("./maps/Threat_maps/Overlaps/ring.richness_cat_map_bins_3cat_hotspot_25.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(ring.richness_cat_map_bins_3cat_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Mimicry richness hotspots (25%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.richness_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()


pdf(file = paste0("./maps/Threat_maps/Overlaps/ring.richness_cat_map_bins_3cat_hotspot_5.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(ring.richness_cat_map_bins_3cat_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Mimicry richness hotspots (5%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.richness_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()


# 4.3.3/ All 4 options ####

pdf(file = paste0("./maps/Threat_maps/Overlaps/ring.richness_cat_map_all_options.pdf"), height = 10, width = 12)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
par(mfrow = c(2,2))


# A/ Hotspots 25%, for quantiles

image(ring.richness_cat_map_rank_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), main = paste0("Mimicry richness hotspots (25%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.richness_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

# B/ Hotspots 5%, for quantiles

image(ring.richness_cat_map_rank_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Mimicry richness hotspots (5%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.richness_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

# C/ Hotspots 25%, for absolute values

image(ring.richness_cat_map_bins_3cat_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Mimicry richness hotspots (25%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.richness_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

# D/ Hotspots 5%, for absolute values

image(ring.richness_cat_map_bins_3cat_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Mimicry richness hotspots (5%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.richness_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
par(mfrow = c(1, 1))

dev.off()

### 5.2/ Mimicry rarity ####

# 5.2.1/ For quantiles ####

plot(ring.rarity_mask_rank_25)
plot(ring.rarity_mask_rank_5)

# Crop with hotspots = 25% highest Mean species rarity
ring.rarity_cat_map_rank_hotspot_25 <- HII_cat_map_rank*ring.rarity_mask_rank_25 # Crop for hotspot only
ring.rarity_cat_map_rank_hotspot_25 <- ring.rarity_cat_map_rank_hotspot_25 + ring.rarity_mask_rank_25 # Add all hotspots
plot(ring.rarity_cat_map_rank_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"))
plot(ring.rarity_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)
saveRDS(ring.rarity_cat_map_rank_hotspot_25, file = "./outputs/Threat_maps/Overlaps/ring.rarity_cat_map_rank_hotspot_25")

# Crop with hotspots = 5% highest Mean species rarity
ring.rarity_cat_map_rank_hotspot_5 <- HII_cat_map_rank*ring.rarity_mask_rank_5
ring.rarity_cat_map_rank_hotspot_5 <- ring.rarity_cat_map_rank_hotspot_5 + ring.rarity_mask_rank_5 # Add all hotspots
plot(ring.rarity_cat_map_rank_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"))
plot(ring.rarity_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)
saveRDS(ring.rarity_cat_map_rank_hotspot_5, file = "./outputs/Threat_maps/Overlaps/ring.rarity_cat_map_rank_hotspot_5")

pdf(file = paste0("./maps/Threat_maps/Overlaps/ring.rarity_cat_map_rank_hotspot_25.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(ring.rarity_cat_map_rank_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), main = paste0("Mean mimicry rarity hotspots (25%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.rarity_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()


pdf(file = paste0("./maps/Threat_maps/Overlaps/ring.rarity_cat_map_rank_hotspot_5.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(ring.rarity_cat_map_rank_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Mean mimicry rarity hotspots (5%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.rarity_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()


# 5.2.2/ For absolute values ####

plot(ring.rarity_mask_rank_25)
plot(ring.rarity_mask_rank_5)

# Crop with hotspots = 25% highest Mean species rarity
ring.rarity_cat_map_bins_3cat_hotspot_25 <- HII_cat_map_bins_3cat*ring.rarity_mask_rank_25
ring.rarity_cat_map_bins_3cat_hotspot_25 <- ring.rarity_cat_map_bins_3cat_hotspot_25 + ring.rarity_mask_rank_25 # Add all hotspots
plot(ring.rarity_cat_map_bins_3cat_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"))
saveRDS(ring.rarity_cat_map_bins_3cat_hotspot_25, file = "./outputs/Threat_maps/Overlaps/ring.rarity_cat_map_bins_hotspot_25")

# Crop with hotspots = 5% highest Mean species rarity
ring.rarity_cat_map_bins_3cat_hotspot_5 <- HII_cat_map_bins_3cat*ring.rarity_mask_rank_5
ring.rarity_cat_map_bins_3cat_hotspot_5 <- ring.rarity_cat_map_bins_3cat_hotspot_5 + ring.rarity_mask_rank_5 # Add all hotspots
plot(ring.rarity_cat_map_bins_3cat_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"))
saveRDS(ring.rarity_cat_map_bins_3cat_hotspot_5, file = "./outputs/Threat_maps/Overlaps/ring.rarity_cat_map_bins_hotspot_5")



pdf(file = paste0("./maps/Threat_maps/Overlaps/ring.rarity_cat_map_bins_3cat_hotspot_25.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(ring.rarity_cat_map_bins_3cat_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Mean mimicry rarity hotspots (25%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.rarity_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()


pdf(file = paste0("./maps/Threat_maps/Overlaps/ring.rarity_cat_map_bins_3cat_hotspot_5.pdf"), height = 6, width = 7.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(ring.rarity_cat_map_bins_3cat_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Mean mimicry rarity hotspots (5%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.rarity_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
dev.off()



# 5.2.3/ All 4 options ####

pdf(file = paste0("./maps/Threat_maps/Overlaps/ring.rarity_cat_map_all_options.pdf"), height = 10, width = 12)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
par(mfrow = c(2,2))


# A/ Hotspots 25%, for quantiles

image(ring.rarity_cat_map_rank_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), main = paste0("Mean mimicry rarity hotspots (25%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.rarity_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

# B/ Hotspots 5%, for quantiles

image(ring.rarity_cat_map_rank_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Mean mimicry rarity hotspots (5%)\n Threat and refuges for quantiles"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.rarity_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

# C/ Hotspots 25%, for absolute values

image(ring.rarity_cat_map_bins_3cat_hotspot_25, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Mean mimicry rarity hotspots (25%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.rarity_mask_rank_25_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

# D/ Hotspots 5%, for absolute values

image(ring.rarity_cat_map_bins_3cat_hotspot_5, col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), main = paste0("Mean mimicry rarity hotspots (5%)\n Threat and refuges for absolute values"),
      cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
      ylab = "", xlab = "",
      legend.args=list(text="     Species", cex=1.2, line = 1, font = 2),
      legend  = T)
plot(crop_mask_shp, lwd = 1, border = "grey20", add = T)
# plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
# plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(ring.rarity_mask_rank_5_shp, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.1, 0.1))

legend(legend = "Risk areas (HII > 23)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Risk areas (HII > 11)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (HII < 1)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

par(mar = internal_margins)
par(mfrow = c(1, 1))

dev.off()


### 6/ Reproject in Mollweide for final plot with all indices ####

# 6.1/ For Hotspots 25% ####

sp.richness_cat_map_rank_hotspot_25 <- readRDS(file = "./outputs/Threat_maps/Overlaps/sp.richness_cat_map_rank_hotspot_25")
sp.rarity_cat_map_rank_hotspot_25 <- readRDS(file = "./outputs/Threat_maps/Overlaps/sp.rarity_cat_map_rank_hotspot_25")
MPD_cat_map_rank_hotspot_25 <- readRDS(file = "./outputs/Threat_maps/Overlaps/MPD_cat_map_rank_hotspot_25")
ring.richness_cat_map_rank_hotspot_25 <- readRDS(file = "./outputs/Threat_maps/Overlaps/ring.richness_cat_map_rank_hotspot_25")
ring.rarity_cat_map_rank_hotspot_25 <- readRDS(file = "./outputs/Threat_maps/Overlaps/ring.rarity_cat_map_rank_hotspot_25")

sp.richness_cat_map_rank_hotspot_25_Mollweide <- projectRaster(from = sp.richness_cat_map_rank_hotspot_25, method = "ngb", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                               crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                               alignOnly = F)
sp.rarity_cat_map_rank_hotspot_25_Mollweide <- projectRaster(from = sp.rarity_cat_map_rank_hotspot_25, method = "ngb", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                             crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                             alignOnly = F)
MPD_cat_map_rank_hotspot_25_Mollweide <- projectRaster(from = MPD_cat_map_rank_hotspot_25, method = "ngb", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                       crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                       alignOnly = F)
ring.richness_cat_map_rank_hotspot_25_Mollweide <- projectRaster(from = ring.richness_cat_map_rank_hotspot_25, method = "ngb", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                                 crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                                 alignOnly = F)
ring.rarity_cat_map_rank_hotspot_25_Mollweide <- projectRaster(from = ring.rarity_cat_map_rank_hotspot_25, method = "ngb", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                             crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                             alignOnly = F)


sp.richness_mask_rank_25_shp_Mollweide <- spTransform(x = sp.richness_mask_rank_25_shp, CRSobj = crs("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
sp.rarity_mask_rank_25_shp_Mollweide <- spTransform(x = sp.rarity_mask_rank_25_shp, CRSobj = crs("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
MPD_mask_rank_25_shp_Mollweide <- spTransform(x = MPD_mask_rank_25_shp, CRSobj = crs("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
ring.richness_mask_rank_25_shp_Mollweide <- spTransform(x = ring.richness_mask_rank_25_shp, CRSobj = crs("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
ring.rarity_mask_rank_25_shp_Mollweide <- spTransform(x = ring.rarity_mask_rank_25_shp, CRSobj = crs("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
crop_mask_shp_Mollweide <- spTransform(x = crop_mask_shp, CRSobj = crs("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))

# 6.2/ For Hotspots 5% ####

sp.richness_cat_map_rank_hotspot_5 <- readRDS(file = "./outputs/Threat_maps/Overlaps/sp.richness_cat_map_rank_hotspot_5")
sp.rarity_cat_map_rank_hotspot_5 <- readRDS(file = "./outputs/Threat_maps/Overlaps/sp.rarity_cat_map_rank_hotspot_5")
MPD_cat_map_rank_hotspot_5 <- readRDS(file = "./outputs/Threat_maps/Overlaps/MPD_cat_map_rank_hotspot_5")
ring.richness_cat_map_rank_hotspot_5 <- readRDS(file = "./outputs/Threat_maps/Overlaps/ring.richness_cat_map_rank_hotspot_5")
ring.rarity_cat_map_rank_hotspot_5 <- readRDS(file = "./outputs/Threat_maps/Overlaps/ring.rarity_cat_map_rank_hotspot_5")


sp.richness_cat_map_rank_hotspot_5_Mollweide <- projectRaster(from = sp.richness_cat_map_rank_hotspot_5, method = "ngb", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                               crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                               alignOnly = F)
sp.rarity_cat_map_rank_hotspot_5_Mollweide <- projectRaster(from = sp.rarity_cat_map_rank_hotspot_5, method = "ngb", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                             crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                             alignOnly = F)
MPD_cat_map_rank_hotspot_5_Mollweide <- projectRaster(from = MPD_cat_map_rank_hotspot_5, method = "ngb", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                       crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                       alignOnly = F)
ring.richness_cat_map_rank_hotspot_5_Mollweide <- projectRaster(from = ring.richness_cat_map_rank_hotspot_5, method = "ngb", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                                 crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                                 alignOnly = F)
ring.rarity_cat_map_rank_hotspot_5_Mollweide <- projectRaster(from = ring.rarity_cat_map_rank_hotspot_5, method = "ngb", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                              crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                              alignOnly = F)


sp.richness_mask_rank_5_shp_Mollweide <- spTransform(x = sp.richness_mask_rank_5_shp, CRSobj = crs("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
sp.rarity_mask_rank_5_shp_Mollweide <- spTransform(x = sp.rarity_mask_rank_5_shp, CRSobj = crs("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
MPD_mask_rank_5_shp_Mollweide <- spTransform(x = MPD_mask_rank_5_shp, CRSobj = crs("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
ring.richness_mask_rank_5_shp_Mollweide <- spTransform(x = ring.richness_mask_rank_5_shp, CRSobj = crs("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
ring.rarity_mask_rank_5_shp_Mollweide <- spTransform(x = ring.rarity_mask_rank_5_shp, CRSobj = crs("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
crop_mask_shp_Mollweide <- spTransform(x = crop_mask_shp, CRSobj = crs("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))


### 7/ Plot all indices for hotspots 25% and quantiles ####


pdf(file = paste0("./maps/Threat_maps/Overlaps/all_indices_cat_map_rank_hotspot_25_Mollweide.pdf"), height = 10, width = 12)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,2.5,2.1))
par(mfrow = c(2,2))

# A/ Species richness

# Plot without axis
image(sp.richness_cat_map_rank_hotspot_25_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), 
      main = paste0("Species richness"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(sp.richness_mask_rank_25_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "a", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0.01, 0))

# B/ Species rarity

# Plot without axis
image(sp.rarity_cat_map_rank_hotspot_25_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), 
      main = paste0("Mean species rarity"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(sp.rarity_mask_rank_25_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "b", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0.01, 0))

# # C/ MPD
# 
# # Plot without axis
# image(MPD_cat_map_rank_hotspot_25_Mollweide, 
#       col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), 
#       main = paste0("Mean pairwise Phylogenetic Distance"),
#       cex.axis = 1.4, cex.main = 1.6, 
#       xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
#       xlab = "", ylab = "")
# 
# # Generate axes with manual positioning of ticks
# axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
# axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
# 
# plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
# plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# # plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)
# 
# plot(MPD_mask_rank_25_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)
# 
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")
# 
# legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
#        fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
# legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
#        fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
# legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
#        fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
# legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
#        fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
# legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
#        fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))
# 
# legend(legend = "c", x = "bottomright", bty = "n",
#        text.font = 2, cex = 1.6, inset = c(0.01, 0))

# C/ Mimicry richness

# Plot without axis
image(ring.richness_cat_map_rank_hotspot_25_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), 
      main = paste0("Mimicry richness"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(ring.richness_mask_rank_25_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "c", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0.01, 0))

# D/ Mimicry rarity

# Plot without axis
image(ring.rarity_cat_map_rank_hotspot_25_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), 
      main = paste0("Mean mimicry rarity"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(ring.rarity_mask_rank_25_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "d", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0.01, 0))


par(mar = internal_margins)
par(mfrow = c(1, 1))

dev.off()


### 8/ Plot all indices for hotspots 5% and quantiles ####


pdf(file = paste0("./maps/Threat_maps/Overlaps/all_indices_cat_map_rank_hotspot_5_Mollweide.pdf"), height = 10, width = 12)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,2.5,2.1))
par(mfrow = c(2,2))

# A/ Species richness

# Plot without axis
image(sp.richness_cat_map_rank_hotspot_5_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), 
      main = paste0("Species richness"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(sp.richness_mask_rank_5_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "a", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0.01, 0))

# B/ Species rarity

# Plot without axis
image(sp.rarity_cat_map_rank_hotspot_5_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), 
      main = paste0("Mean species rarity"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(sp.rarity_mask_rank_5_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "b", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0.01, 0))

# # C/ MPD
# 
# # Plot without axis
# image(MPD_cat_map_rank_hotspot_5_Mollweide, 
#       col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), 
#       main = paste0("Mean pairwise Phylogenetic Distance"),
#       cex.axis = 1.4, cex.main = 1.6, 
#       xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
#       xlab = "", ylab = "")
# 
# # Generate axes with manual positioning of ticks
# axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
# axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
# 
# plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
# plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# # plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)
# 
# plot(MPD_mask_rank_5_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)
# 
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")
# 
# legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
#        fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
# legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
#        fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
# legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
#        fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
# legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
#        fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
# legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
#        fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))
# 
# legend(legend = "c", x = "bottomright", bty = "n",
#        text.font = 2, cex = 1.6, inset = c(0.01, 0))

# C/ Mimicry richness

# Plot without axis
image(ring.richness_cat_map_rank_hotspot_5_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), 
      main = paste0("Mimicry richness"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(ring.richness_mask_rank_5_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "c", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0.01, 0))

# D/ Mimicry rarity

# Plot without axis
image(ring.rarity_cat_map_rank_hotspot_5_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), 
      main = paste0("Mean mimicry rarity"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(ring.rarity_mask_rank_5_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "d", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0.01, 0))


par(mar = internal_margins)
par(mfrow = c(1, 1))

dev.off()


### 9/ Plot all indices for hotspots 5% and 25% and quantiles on the same 4 mini-maps ####

pdf(file = paste0("./maps/Threat_maps/Overlaps/all_indices_cat_map_rank_hotspot_5_25_Mollweide.pdf"), height = 10, width = 12)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,2.5,2.1))
par(mfrow = c(2,2))

# A/ Species richness

# Plot without axis
image(sp.richness_cat_map_rank_hotspot_25_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), 
      main = paste0("Species richness"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(sp.richness_mask_rank_25_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

plot(sp.richness_mask_rank_5_shp_Mollweide, lwd = 1, border = "#333333BB", density = 15, angle = 45, add = T)
plot(sp.richness_mask_rank_5_shp_Mollweide, lwd = 1, border = "#333333BB", density = 15, angle = 135, add = T)
plot(sp.richness_mask_rank_5_shp_Mollweide, lwd = 1.2, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Top 5% hotspots", x = "bottomleft", bty = "n",
       density = 15, angle = 45, 
       fill = "black", text.font = 2, cex = 1.2, inset=c(0.02, 0.48))
legend(legend = "", x = "bottomleft", bty = "n",
       density = 15, angle = 135,
       fill = "black", text.font = 2, cex = 1.2, inset=c(0.02, 0.48))
legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "a", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0.01, 0))

# B/ Species rarity

# Plot without axis
image(sp.rarity_cat_map_rank_hotspot_25_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), 
      main = paste0("Mean species rarity"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(sp.rarity_mask_rank_25_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

plot(sp.rarity_mask_rank_5_shp_Mollweide, lwd = 1, border = "#333333BB", density = 15, angle = 45, add = T)
plot(sp.rarity_mask_rank_5_shp_Mollweide, lwd = 1, border = "#333333BB", density = 15, angle = 135, add = T)
plot(sp.rarity_mask_rank_5_shp_Mollweide, lwd = 1.2, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Top 5% hotspots", x = "bottomleft", bty = "n",
       density = 15, angle = 45, 
       fill = "black", text.font = 2, cex = 1.2, inset=c(0.02, 0.48))
legend(legend = "", x = "bottomleft", bty = "n",
       density = 15, angle = 135,
       fill = "black", text.font = 2, cex = 1.2, inset=c(0.02, 0.48))
legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "b", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0.01, 0))


# # C/ MPD
# 
# # Plot without axis
# image(MPD_cat_map_rank_hotspot_25_Mollweide, 
#       col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), 
#       main = paste0("Mean pairwise Phylogenetic Distance"),
#       cex.axis = 1.4, cex.main = 1.6, 
#       xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
#       xlab = "", ylab = "")
# 
# # Generate axes with manual positioning of ticks
# axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
# axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
# 
# plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
# plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# # plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)
# 
# plot(MPD_mask_rank_25_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)
# 
# plot(MPD_mask_rank_5_shp_Mollweide, lwd = 1, border = "#333333BB", density = 15, angle = 45, add = T)
# plot(MPD_mask_rank_5_shp_Mollweide, lwd = 1, border = "#333333BB", density = 15, angle = 135, add = T)
# plot(MPD_mask_rank_5_shp_Mollweide, lwd = 1.2, border = "#333333BB", add = T)
# 
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")
# 
# legend(legend = "Top 5% hotspots", x = "bottomleft", bty = "n",
#        density = 15, angle = 45, 
#        fill = "black", text.font = 2, cex = 1.2, inset=c(0.02, 0.48))
# legend(legend = "", x = "bottomleft", bty = "n",
#        density = 15, angle = 135,
#        fill = "black", text.font = 2, cex = 1.2, inset=c(0.02, 0.48))
# legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
#        fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
# legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
#        fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
# legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
#        fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
# legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
#        fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
# legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
#        fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))
# 
# legend(legend = "c", x = "bottomright", bty = "n",
#        text.font = 2, cex = 1.6, inset = c(0.01, 0))

# C/ Mimicry richness

# Plot without axis
image(ring.richness_cat_map_rank_hotspot_25_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), 
      main = paste0("Mimicry richness"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(ring.richness_mask_rank_25_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

plot(ring.richness_mask_rank_5_shp_Mollweide, lwd = 1, border = "#333333BB", density = 15, angle = 45, add = T)
plot(ring.richness_mask_rank_5_shp_Mollweide, lwd = 1, border = "#333333BB", density = 15, angle = 135, add = T)
plot(ring.richness_mask_rank_5_shp_Mollweide, lwd = 1.2, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Top 5% hotspots", x = "bottomleft", bty = "n",
       density = 15, angle = 45, 
       fill = "black", text.font = 2, cex = 1.2, inset=c(0.02, 0.48))
legend(legend = "", x = "bottomleft", bty = "n",
       density = 15, angle = 135,
       fill = "black", text.font = 2, cex = 1.2, inset=c(0.02, 0.48))
legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "c", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0.01, 0))

# D/ Species rarity

# Plot without axis
image(ring.rarity_cat_map_rank_hotspot_25_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), 
      main = paste0("Mean species rarity"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(ring.rarity_mask_rank_25_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

plot(ring.rarity_mask_rank_5_shp_Mollweide, lwd = 1, border = "#333333BB", density = 15, angle = 45, add = T)
plot(ring.rarity_mask_rank_5_shp_Mollweide, lwd = 1, border = "#333333BB", density = 15, angle = 135, add = T)
plot(ring.rarity_mask_rank_5_shp_Mollweide, lwd = 1.2, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Top 5% hotspots", x = "bottomleft", bty = "n",
       density = 15, angle = 45, 
       fill = "black", text.font = 2, cex = 1.2, inset=c(0.02, 0.48))
legend(legend = "", x = "bottomleft", bty = "n",
       density = 15, angle = 135,
       fill = "black", text.font = 2, cex = 1.2, inset=c(0.02, 0.48))
legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "d", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.6, inset = c(0.01, 0))


par(mar = internal_margins)
par(mfrow = c(1, 1))

dev.off()


### 10/ Plot all indices for hotspots 25% and 5% on 8 mini-maps and quantiles ####


pdf(file = paste0("./maps/Threat_maps/Overlaps/all_indices_cat_map_rank_hotspot_5_25_Mollweide_8_maps.pdf"), height = 16, width = 9)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,4.5,2.1))
par(mfrow = c(4,2))

# A/ Species richness for hotspot 25%

# Plot without axis
image(sp.richness_cat_map_rank_hotspot_25_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), 
      main = paste0("Species richness\nTop 25% hotspots"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(sp.richness_mask_rank_25_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "(a)", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset = c(0.01, 0))

# B/ Species richness for hotspot 5%

# Plot without axis
image(sp.richness_cat_map_rank_hotspot_5_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), 
      main = paste0("Species richness\nTop 5% hotspots"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(sp.richness_mask_rank_5_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "(b)", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset = c(0.01, 0))

# C/ Species rarity for hotspot 25%

# Plot without axis
image(sp.rarity_cat_map_rank_hotspot_25_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), 
      main = paste0("Species geographic rarity\nTop 25% hotspots"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(sp.rarity_mask_rank_25_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "(c)", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset = c(0.01, 0))

# D/ Species rarity for hotspot 5%

# Plot without axis
image(sp.rarity_cat_map_rank_hotspot_5_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), 
      main = paste0("Species geographic rarity\nTop 5% hotspots"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(sp.rarity_mask_rank_5_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "(d)", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset = c(0.01, 0))

# # E/ MPD for hotspot 25%
# 
# # Plot without axis
# image(MPD_cat_map_rank_hotspot_25_Mollweide, 
#       col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), 
#       main = paste0("Mean pairwise Phylogenetic Distance\nTop 25% hotspots"),
#       cex.axis = 1.4, cex.main = 1.6, 
#       xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
#       xlab = "", ylab = "")
# 
# # Generate axes with manual positioning of ticks
# axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
# axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
# 
# plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
# plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# # plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)
# 
# plot(MPD_mask_rank_25_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)
# 
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")
# 
# legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
#        fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
# legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
#        fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
# legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
#        fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
# legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
#        fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
# legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
#        fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))
# 
# legend(legend = "e", x = "bottomright", bty = "n",
#        text.font = 2, cex = 1.6, inset = c(0.01, 0))
# 
# # F/ MPD for hotspot 5%
# 
# # Plot without axis
# image(MPD_cat_map_rank_hotspot_5_Mollweide, 
#       col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), 
#       main = paste0("Mean pairwise Phylogenetic Distance\nTop 5% hotspots"),
#       cex.axis = 1.4, cex.main = 1.6, 
#       xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
#       xlab = "", ylab = "")
# 
# # Generate axes with manual positioning of ticks
# axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
# axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
# 
# plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
# plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
# plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# # plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)
# 
# plot(MPD_mask_rank_5_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)
# 
# scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
# prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")
# 
# legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
#        fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
# legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
#        fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
# legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
#        fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
# legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
#        fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
# legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
#        fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))
# 
# legend(legend = "f", x = "bottomright", bty = "n",
#        text.font = 2, cex = 1.6, inset = c(0.01, 0))

# E/ Mimicry richness for hotspot 25%

# Plot without axis
image(ring.richness_cat_map_rank_hotspot_25_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), 
      main = paste0("Mimicry richness\nTop 25% hotspots"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(ring.richness_mask_rank_25_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "(e)", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset = c(0.01, 0))

# F/ Mimicry richness for hotspot 5%

# Plot without axis
image(ring.richness_cat_map_rank_hotspot_5_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), 
      main = paste0("Mimicry richness\nTop 5% hotspots"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(ring.richness_mask_rank_5_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "(f)", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset = c(0.01, 0))

# G/ Mimicry rarity for hotspot 25%

# Plot without axis
image(ring.rarity_cat_map_rank_hotspot_25_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue", "navyblue"), 
      main = paste0("Mimicry geographic rarity\nTop 25% hotspots"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(ring.rarity_mask_rank_25_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "(g)", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset = c(0.01, 0))

# H/ Mimicry rarity for hotspot 5%

# Plot without axis
image(ring.rarity_cat_map_rank_hotspot_5_Mollweide, 
      col = c("#EDEDED", "papayawhip", "red4", "red", "dodgerblue"), 
      main = paste0("Mimicry geographic rarity\nTop 5% hotspots"),
      cex.axis = 1.4, cex.main = 1.6, 
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.3, lwd = 0.2, lwd.ticks = 1)

plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
# plot(crop_mask_shp_Mollweide, lwd = 1, border = "grey20", add = T)

plot(ring.rarity_mask_rank_5_shp_Mollweide, lwd = 0.4, border = "#333333BB", add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4200, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.15, 0.15), text.col = "#00000000")

legend(legend = "Risk areas (5%)", x = "bottomleft", bty = "n",
       fill = "red4", text.font = 2, cex = 1.2, inset=c(0.02, 0.41))
legend(legend = "Risk areas (25%)", x = "bottomleft", bty = "n",
       fill = "red", text.font = 2, cex = 1.2, inset=c(0.02, 0.34))
legend(legend = "Other hotspots", x = "bottomleft", bty = "n",
       fill = "papayawhip", text.font = 2, cex = 1.2, inset=c(0.02, 0.27))
legend(legend = "Refuge areas (25%)", x = "bottomleft", bty = "n",
       fill = "dodgerblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.20))
legend(legend = "Refuge areas (5%)", x = "bottomleft", bty = "n",
       fill = "navyblue", text.font = 2, cex = 1.2, inset=c(0.02, 0.13))

legend(legend = "(h)", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset = c(0.01, 0))


par(mar = internal_margins)
par(mfrow = c(1, 1))

dev.off()


### 11/ Bar plots for proportions of overlap with hotspots 25% ####

library(tidyverse)

### Categories :

# 0 = No hotspot
# 1 = No refuge or risk area
# 2 = High risk area (5%)
# 3 = Risk area (25%)
# 4 = Refuge (25%)
# 5 = Top refuge (5%)

# 11.1/ For quantiles and hotspots 25%

# Load categorical maps
sp.richness_cat_map_rank_hotspot_25 <- readRDS(file = "./outputs/Threat_maps/Overlaps/sp.richness_cat_map_rank_hotspot_25")
sp.rarity_cat_map_rank_hotspot_25 <- readRDS(file = "./outputs/Threat_maps/Overlaps/sp.rarity_cat_map_rank_hotspot_25")
# MPD_cat_map_rank_hotspot_25 <- readRDS(file = "./outputs/Threat_maps/Overlaps/MPD_cat_map_rank_hotspot_25")
ring.richness_cat_map_rank_hotspot_25 <- readRDS(file = "./outputs/Threat_maps/Overlaps/ring.richness_cat_map_rank_hotspot_25")
ring.rarity_cat_map_rank_hotspot_25 <- readRDS(file = "./outputs/Threat_maps/Overlaps/ring.rarity_cat_map_rank_hotspot_25")

# Extract data for species richness
sp.richness_N.hotspots <- sum(sp.richness_cat_map_rank_hotspot_25[] > 0, na.rm = T)
sp.richness_N.5 <- sum(sp.richness_cat_map_rank_hotspot_25[] == 2, na.rm = T)
sp.richness_N.25 <- sum(sp.richness_cat_map_rank_hotspot_25[] == 3, na.rm = T) + sp.richness_N.5
sp.richness_N.95 <- sum(sp.richness_cat_map_rank_hotspot_25[] == 5, na.rm = T)
sp.richness_N.75 <- sum(sp.richness_cat_map_rank_hotspot_25[] == 4, na.rm = T) + sp.richness_N.95

sp.richness_N.quantiles <- c(sp.richness_N.5, sp.richness_N.25, sp.richness_N.75, sp.richness_N.95)
sp.richness_N.quantiles_perc <- round(sp.richness_N.quantiles/sp.richness_N.hotspots*100, 1)
sp.richness_labels <-  paste0(sp.richness_N.quantiles, "\n(", sp.richness_N.quantiles_perc, "%)")
quantiles <- c("Risk areas (5%)", "Risk areas (25%)", "Refuges (25%)", "Refuges (5%)")
sp.richness_N.quantiles_df <- data.frame(cbind(sp.richness_N.quantiles, quantiles, sp.richness_labels))
sp.richness_N.quantiles_df$sp.richness_N.quantiles <- as.numeric(as.character(sp.richness_N.quantiles_df$sp.richness_N.quantiles))

saveRDS(sp.richness_N.quantiles_df, file = "./outputs/Threat_maps/Overlaps/sp.richness_N.quantiles_df.rds", version = "2")

# Extract data for species rarity

sp.rarity_N.hotspots <- sum(sp.rarity_cat_map_rank_hotspot_25[] > 0, na.rm = T)
sp.rarity_N.5 <- sum(sp.rarity_cat_map_rank_hotspot_25[] == 2, na.rm = T)
sp.rarity_N.25 <- sum(sp.rarity_cat_map_rank_hotspot_25[] == 3, na.rm = T) + sp.rarity_N.5
sp.rarity_N.95 <- sum(sp.rarity_cat_map_rank_hotspot_25[] == 5, na.rm = T)
sp.rarity_N.75 <- sum(sp.rarity_cat_map_rank_hotspot_25[] == 4, na.rm = T) + sp.rarity_N.95

sp.rarity_N.quantiles <- c(sp.rarity_N.5, sp.rarity_N.25, sp.rarity_N.75, sp.rarity_N.95)
sp.rarity_N.quantiles_perc <- round(sp.rarity_N.quantiles/sp.rarity_N.hotspots*100, 1)
sp.rarity_labels <-  paste0(sp.rarity_N.quantiles, "\n(", sp.rarity_N.quantiles_perc, "%)")
quantiles <- c("Risk areas (5%)", "Risk areas (25%)", "Refuges (25%)", "Refuges (5%)")
sp.rarity_N.quantiles_df <- data.frame(cbind(sp.rarity_N.quantiles, quantiles, sp.rarity_labels))
sp.rarity_N.quantiles_df$sp.rarity_N.quantiles <- as.numeric(as.character(sp.rarity_N.quantiles_df$sp.rarity_N.quantiles))

saveRDS(sp.rarity_N.quantiles_df, file = "./outputs/Threat_maps/Overlaps/sp.rarity_N.quantiles_df.rds", version = "2")


# Extract data for MPD

MPD_N.hotspots <- sum(MPD_cat_map_rank_hotspot_25[] > 0, na.rm = T)
MPD_N.5 <- sum(MPD_cat_map_rank_hotspot_25[] == 2, na.rm = T)
MPD_N.25 <- sum(MPD_cat_map_rank_hotspot_25[] == 3, na.rm = T) + MPD_N.5
MPD_N.95 <- sum(MPD_cat_map_rank_hotspot_25[] == 5, na.rm = T)
MPD_N.75 <- sum(MPD_cat_map_rank_hotspot_25[] == 4, na.rm = T) + MPD_N.95

MPD_N.quantiles <- c(MPD_N.5, MPD_N.25, MPD_N.75, MPD_N.95)
MPD_N.quantiles_perc <- round(MPD_N.quantiles/MPD_N.hotspots*100, 1)
MPD_labels <-  paste0(MPD_N.quantiles, "\n(", MPD_N.quantiles_perc, "%)")
quantiles <- c("Risk areas (5%)", "Risk areas (25%)", "Refuges (25%)", "Refuges (5%)")
MPD_N.quantiles_df <- data.frame(cbind(MPD_N.quantiles, quantiles, MPD_labels))
MPD_N.quantiles_df$MPD_N.quantiles <- as.numeric(as.character(MPD_N.quantiles_df$MPD_N.quantiles))

saveRDS(MPD_N.quantiles_df, file = "./outputs/Threat_maps/Overlaps/MPD_N.quantiles_df.rds", version = "2")


# Extract data for mimicry richness

ring.richness_N.hotspots <- sum(ring.richness_cat_map_rank_hotspot_25[] > 0, na.rm = T)
ring.richness_N.5 <- sum(ring.richness_cat_map_rank_hotspot_25[] == 2, na.rm = T)
ring.richness_N.25 <- sum(ring.richness_cat_map_rank_hotspot_25[] == 3, na.rm = T) + ring.richness_N.5
ring.richness_N.95 <- sum(ring.richness_cat_map_rank_hotspot_25[] == 5, na.rm = T)
ring.richness_N.75 <- sum(ring.richness_cat_map_rank_hotspot_25[] == 4, na.rm = T) + ring.richness_N.95

ring.richness_N.quantiles <- c(ring.richness_N.5, ring.richness_N.25, ring.richness_N.75, ring.richness_N.95)
ring.richness_N.quantiles_perc <- round(ring.richness_N.quantiles/ring.richness_N.hotspots*100, 1)
ring.richness_labels <-  paste0(ring.richness_N.quantiles, "\n(", ring.richness_N.quantiles_perc, "%)")
quantiles <- c("Risk areas (5%)", "Risk areas (25%)", "Refuges (25%)", "Refuges (5%)")
ring.richness_N.quantiles_df <- data.frame(cbind(ring.richness_N.quantiles, quantiles, ring.richness_labels))
ring.richness_N.quantiles_df$ring.richness_N.quantiles <- as.numeric(as.character(ring.richness_N.quantiles_df$ring.richness_N.quantiles))

saveRDS(ring.richness_N.quantiles_df, file = "./outputs/Threat_maps/Overlaps/ring.richness_N.quantiles_df.rds", version = "2")


# Extract data for mimicry rarity

ring.rarity_N.hotspots <- sum(ring.rarity_cat_map_rank_hotspot_25[] > 0, na.rm = T)
ring.rarity_N.5 <- sum(ring.rarity_cat_map_rank_hotspot_25[] == 2, na.rm = T)
ring.rarity_N.25 <- sum(ring.rarity_cat_map_rank_hotspot_25[] == 3, na.rm = T) + ring.rarity_N.5
ring.rarity_N.95 <- sum(ring.rarity_cat_map_rank_hotspot_25[] == 5, na.rm = T)
ring.rarity_N.75 <- sum(ring.rarity_cat_map_rank_hotspot_25[] == 4, na.rm = T) + ring.rarity_N.95

ring.rarity_N.quantiles <- c(ring.rarity_N.5, ring.rarity_N.25, ring.rarity_N.75, ring.rarity_N.95)
ring.rarity_N.quantiles_perc <- round(ring.rarity_N.quantiles/ring.rarity_N.hotspots*100, 1)
ring.rarity_labels <-  paste0(ring.rarity_N.quantiles, "\n(", ring.rarity_N.quantiles_perc, "%)")
quantiles <- c("Risk areas (5%)", "Risk areas (25%)", "Refuges (25%)", "Refuges (5%)")
ring.rarity_N.quantiles_df <- data.frame(cbind(ring.rarity_N.quantiles, quantiles, ring.rarity_labels))
ring.rarity_N.quantiles_df$ring.rarity_N.quantiles <- as.numeric(as.character(ring.rarity_N.quantiles_df$ring.rarity_N.quantiles))

saveRDS(ring.rarity_N.quantiles_df, file = "./outputs/Threat_maps/Overlaps/ring.rarity_N.quantiles_df.rds", version = "2")



# Build df with all data

names(sp.richness_N.quantiles_df) <- names(sp.rarity_N.quantiles_df) <- names(MPD_N.quantiles_df) <- names(ring.richness_N.quantiles_df) <- names(ring.rarity_N.quantiles_df) <- c("N.quantiles", "quantiles", "labels")
all_indices_N.quantiles_df <- rbind(sp.richness_N.quantiles_df,
                                    sp.rarity_N.quantiles_df,
                                    # MPD_N.quantiles_df,
                                    ring.richness_N.quantiles_df,
                                    ring.rarity_N.quantiles_df)
all_indices_N.quantiles_df$indices <- c(rep("sp.richness", 4),
                                        rep("sp.rarity", 4),
                                        # rep("MPD", 4),
                                        rep("ring.richness", 4),
                                        rep("ring.rarity", 4))
all_indices_N.quantiles_df$indices <- factor(all_indices_N.quantiles_df$indices,
                                             # levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness"))
                                             levels = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))

saveRDS(all_indices_N.quantiles_df, file = "./outputs/Threat_maps/Overlaps/all_indices_N.quantiles_df.rds", version = "2")

# Plot

all_indices_N.quantiles_df <- readRDS(file = "./outputs/Threat_maps/Overlaps/all_indices_N.quantiles_df.rds")

empty_rows <- data.frame(matrix(data = NA, nrow = 4, ncol = 4))
names(empty_rows) <- names(all_indices_N.quantiles_df)
# empty_rows[, 4] <- c("sp.richness", "sp.rarity", "MPD", "ring.richness")
empty_rows[, 4] <- c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity")
empty_rows[, 2] <- ""
empty_rows[, 1] <-  0
all_indices_N.quantiles_df <- rbind(all_indices_N.quantiles_df, empty_rows)

# Without parenthesis in labels
all_indices_N.quantiles_df$labels <- str_remove(string = str_remove(string = all_indices_N.quantiles_df$labels, pattern = "\\("), pattern = "\\)")

str_remove(string = all_indices_N.quantiles_df$labels, pattern = "\\(.+\\)") # Remove all % lines
str_remove(string = all_indices_N.quantiles_df$labels, pattern = "[\\)|\\(]") # Does work, but should

# With only %
all_indices_N.quantiles_df$labels <- str_remove(string = all_indices_N.quantiles_df$labels, pattern = ".+[\n]")

# With group names on two lines
all_indices_N.quantiles_df$quantiles <- as.factor(rep(c("Risk areas\n(5%)", "Risk areas\n(25%)", "Refuges\n(25%)", "Refuges\n(5%)"), 5))

# Remove labels overlapped by the 25% and 5% lines
all_indices_N.quantiles_df$labels[c(2, # 8,
                                    9,
                                    # 13,
                                    # 15,
                                    16
                                    )] <- ""

limits_df <- tibble(x = c(0, 150), y = c(0, 3500))

pdf(file = "./graphs/Overlaps/areas_overlap_quantiles_hotspot_25.pdf", height = 8, width = 12)

# levels(all_indices_N.quantiles_df$quantiles)

# Play with witdh and the extra empty level for N.quantiles to get a space inbetween group of bars

# internal_margins <- par()$mar
# par(mar = c(3.1, 3.5, 3.5, 2.1))
# external_margins <- par()$oma
# par(oma = c(0, 0, 3, 0))

# rgb2hsv(col2rgb("darkgoldenrod1"))
# hsv(h = 0.13, s = 0.7, v = 1)

g1 <- ggplot(all_indices_N.quantiles_df, aes(x = quantiles, y = N.quantiles, fill = indices)) +
  geom_col(position = position_dodge(width = 0.9), 
           # width = c(0.9,0.9,0.9,0.9,
           #           0.9,0.9,0.9,0.9,
           #           0.9,0.9,0.9,0.9,
           #           0.9,0.9,0.9,0.9, 
           #           0.1,0.1,0.1,0.1),
           show.legend = T) +
  labs(title = "", 
       y = "Communities", x = NULL) +
  # scale_x_discrete(limits = c("Risk areas (5%)", "Risk areas (25%)", "Refuges (25%)", "Refuges (5%)")) +
  scale_x_discrete(limits = c("Risk areas\n(25%)", "Refuges\n(25%)", "", "Risk areas\n(5%)", "Refuges\n(5%)")) +
  # scale_fill_manual(name = "Indices",
  #                   values = c("sp.richness" = "#FF4040", "sp.rarity" = "#80C2FF", "MPD" = "#FFD84D", "ring.richness" = "#7ACC29"),
  #                   labels = c("sp.richness" = "Species richness", "sp.rarity" = "Species rarity", "MPD" = "MPD", "ring.richness" = "Mimicry richness")) +
  scale_fill_manual(name = "Indices",
                    values = c("sp.richness" = "#FF4040", "sp.rarity" = "#80C2FF", "ring.richness" = "#FFD84D", "ring.rarity" = "#7ACC29"),
                    labels = c("sp.richness" = "Species richness", "sp.rarity" = "Species rarity", "ring.richness" = "Mimicry richness", "ring.rarity" = "Mimicry rarity")) +
   
  # theme_classic() +
  scale_y_continuous(limits = c(0, 4000), expand = c(0, 0), breaks = seq(0, 4000, 1000)) +
  # ylim = c(0, 3500) +
  # geom_rangeframe(# data = limits_df, aes(x = x, y = y),
  #                 size = 1.4, sides = "l") +
  theme(panel.background = element_rect(fill = NA),
        legend.position = c(0.82, 0.57),
        # legend.key = element_rect(fill = "white", colour = "white"),
        legend.title = element_text(size = 24, vjust = 3, face = "bold"),
        legend.text = element_text(size = 22, face = "plain"),
        # margin = margin(t = 0, unit = "pt")),
        legend.key.size = unit(2.5, 'lines'),
        # legend.spacing.y = unit(1,"cm"),
        # axis.line.x = element_line(color = NA),
        axis.line.y = element_line(color = "black", size = 1.3),
        axis.ticks = element_line(color = "black"),
        axis.ticks.x = element_line(size = NA),
        axis.ticks.y = element_line(size = 1.2),
        axis.ticks.length.y = unit(8, "pt"),
        axis.text = element_text(size = 22, face = "bold", color = "black"),
        axis.text.x = element_text(vjust = 0, size = 24, margin = margin(t = 8, b = 15)),
        axis.title = element_text(size = 28, face = "bold", color = "black"),
        # axis.title.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title.y = element_text(size = 28, margin = margin(l = 5, r = 10))) +
  geom_segment(aes(x = 0.5, y = 1338.5, xend = 2.5, yend = 1338.5), lwd = 1, color = "grey30", lty = "dashed") +
  
  # geom_segment(aes(x = 0.5, y = 1338.5, xend = 0.75, yend = 1338.5), lwd = 1, color = "grey30", lty = "dashed") +
  # geom_segment(aes(x = 1, y = 1338.5, xend = 1.55, yend = 1338.5), lwd = 1, color = "grey30", lty = "dashed") +
  # geom_segment(aes(x = 1.8, y = 1338.5, xend = 2.0, yend = 1338.5), lwd = 1, color = "grey30", lty = "dashed") +
  # geom_segment(aes(x = 2.48, y = 1338.5, xend = 2.5, yend = 1338.5), lwd = 1, color = "grey30", lty = "dashed") +
  
  geom_segment(aes(x = 3.5, y = 267.7, xend = 5.5, yend = 267.7), lwd = 1, color = "grey30", lty = "dashed") +
  geom_text(aes(label = labels), vjust=-0.8, size = 6, fontface = "bold", position = position_dodge(width = 0.9)) +
  
  annotate("text", x = 0.67, y = 1460, size = 6, fontface = "bold", label = "23.9%") +
  # annotate("text", x = 2.35, y = 1460, size = 6, fontface = "bold", label = "22.2%") +
  annotate("text", x = 4.11, y = 400, size = 6, fontface = "bold", label = "3%") +
  # annotate("text", x = 4.34, y = 400, size = 6, fontface = "bold", label = "3.0%") +
  # annotate("text", x = 4.90, y = 180, size = 6, fontface = "bold", label = "1.6%") +
  annotate("text", x = 5.35, y = 380, size = 6, fontface = "bold", label = "2.3%") +
   
  annotate("text", x = 3, y = 3875, size = 11, fontface = "bold", label = "Top 25% hotspots") +
  annotate("text", x = 0.7, y = 3875, size = 11, fontface = "bold", label = "(a)") +
  
  coord_cartesian(clip = "off") # To allow to draw label outside of the plot area

print(g1)

par(mar = internal_margins, oma = external_margins)

dev.off()



### 12/ Bar plots for proportions of overlap with hotspots 5% ####

library(tidyverse)

### Categories :

# 0 = No hotspot
# 1 = No refuge or risk area
# 2 = High risk area (5%)
# 3 = Risk area (25%)
# 4 = Refuge (25%)
# 5 = Top refuge (5%)

# 10.1/ For quantiles and hotspots 5%

# Load categorical maps
sp.richness_cat_map_rank_hotspot_5 <- readRDS(file = "./outputs/Threat_maps/Overlaps/sp.richness_cat_map_rank_hotspot_5")
sp.rarity_cat_map_rank_hotspot_5 <- readRDS(file = "./outputs/Threat_maps/Overlaps/sp.rarity_cat_map_rank_hotspot_5")
MPD_cat_map_rank_hotspot_5 <- readRDS(file = "./outputs/Threat_maps/Overlaps/MPD_cat_map_rank_hotspot_5")
ring.richness_cat_map_rank_hotspot_5 <- readRDS(file = "./outputs/Threat_maps/Overlaps/ring.richness_cat_map_rank_hotspot_5")
ring.rarity_cat_map_rank_hotspot_5 <- readRDS(file = "./outputs/Threat_maps/Overlaps/ring.rarity_cat_map_rank_hotspot_5")


# Extract data for species richness
sp.richness_N.hotspots <- sum(sp.richness_cat_map_rank_hotspot_5[] > 0, na.rm = T)
sp.richness_N.5 <- sum(sp.richness_cat_map_rank_hotspot_5[] == 2, na.rm = T)
sp.richness_N.25 <- sum(sp.richness_cat_map_rank_hotspot_5[] == 3, na.rm = T) + sp.richness_N.5
sp.richness_N.95 <- sum(sp.richness_cat_map_rank_hotspot_5[] == 5, na.rm = T)
sp.richness_N.75 <- sum(sp.richness_cat_map_rank_hotspot_5[] == 4, na.rm = T) + sp.richness_N.95

sp.richness_N.quantiles <- c(sp.richness_N.5, sp.richness_N.25, sp.richness_N.75, sp.richness_N.95)
sp.richness_N.quantiles_perc <- round(sp.richness_N.quantiles/sp.richness_N.hotspots*100, 1)
sp.richness_labels <-  paste0(sp.richness_N.quantiles, "\n(", sp.richness_N.quantiles_perc, "%)")
quantiles <- c("Risk areas (5%)", "Risk areas (25%)", "Refuges (25%)", "Refuges (5%)")
sp.richness_N.quantiles_df_5 <- data.frame(cbind(sp.richness_N.quantiles, quantiles, sp.richness_labels))
sp.richness_N.quantiles_df_5$sp.richness_N.quantiles <- as.numeric(as.character(sp.richness_N.quantiles_df_5$sp.richness_N.quantiles))

saveRDS(sp.richness_N.quantiles_df_5, file = "./outputs/Threat_maps/Overlaps/sp.richness_N.quantiles_df_5.rds", version = "2")

# Extract data for species rarity

sp.rarity_N.hotspots <- sum(sp.rarity_cat_map_rank_hotspot_5[] > 0, na.rm = T)
sp.rarity_N.5 <- sum(sp.rarity_cat_map_rank_hotspot_5[] == 2, na.rm = T)
sp.rarity_N.25 <- sum(sp.rarity_cat_map_rank_hotspot_5[] == 3, na.rm = T) + sp.rarity_N.5
sp.rarity_N.95 <- sum(sp.rarity_cat_map_rank_hotspot_5[] == 5, na.rm = T)
sp.rarity_N.75 <- sum(sp.rarity_cat_map_rank_hotspot_5[] == 4, na.rm = T) + sp.rarity_N.95

sp.rarity_N.quantiles <- c(sp.rarity_N.5, sp.rarity_N.25, sp.rarity_N.75, sp.rarity_N.95)
sp.rarity_N.quantiles_perc <- round(sp.rarity_N.quantiles/sp.rarity_N.hotspots*100, 1)
sp.rarity_labels <-  paste0(sp.rarity_N.quantiles, "\n(", sp.rarity_N.quantiles_perc, "%)")
quantiles <- c("Risk areas (5%)", "Risk areas (25%)", "Refuges (25%)", "Refuges (5%)")
sp.rarity_N.quantiles_df_5 <- data.frame(cbind(sp.rarity_N.quantiles, quantiles, sp.rarity_labels))
sp.rarity_N.quantiles_df_5$sp.rarity_N.quantiles <- as.numeric(as.character(sp.rarity_N.quantiles_df_5$sp.rarity_N.quantiles))

saveRDS(sp.rarity_N.quantiles_df_5, file = "./outputs/Threat_maps/Overlaps/sp.rarity_N.quantiles_df_5.rds", version = "2")


# Extract data for MPD

MPD_N.hotspots <- sum(MPD_cat_map_rank_hotspot_5[] > 0, na.rm = T)
MPD_N.5 <- sum(MPD_cat_map_rank_hotspot_5[] == 2, na.rm = T)
MPD_N.25 <- sum(MPD_cat_map_rank_hotspot_5[] == 3, na.rm = T) + MPD_N.5
MPD_N.95 <- sum(MPD_cat_map_rank_hotspot_5[] == 5, na.rm = T)
MPD_N.75 <- sum(MPD_cat_map_rank_hotspot_5[] == 4, na.rm = T) + MPD_N.95

MPD_N.quantiles <- c(MPD_N.5, MPD_N.25, MPD_N.75, MPD_N.95)
MPD_N.quantiles_perc <- round(MPD_N.quantiles/MPD_N.hotspots*100, 1)
MPD_labels <-  paste0(MPD_N.quantiles, "\n(", MPD_N.quantiles_perc, "%)")
quantiles <- c("Risk areas (5%)", "Risk areas (25%)", "Refuges (25%)", "Refuges (5%)")
MPD_N.quantiles_df_5 <- data.frame(cbind(MPD_N.quantiles, quantiles, MPD_labels))
MPD_N.quantiles_df_5$MPD_N.quantiles <- as.numeric(as.character(MPD_N.quantiles_df_5$MPD_N.quantiles))

saveRDS(MPD_N.quantiles_df_5, file = "./outputs/Threat_maps/Overlaps/MPD_N.quantiles_df_5.rds", version = "2")


# Extract data for mimicry richness

ring.richness_N.hotspots <- sum(ring.richness_cat_map_rank_hotspot_5[] > 0, na.rm = T)
ring.richness_N.5 <- sum(ring.richness_cat_map_rank_hotspot_5[] == 2, na.rm = T)
ring.richness_N.25 <- sum(ring.richness_cat_map_rank_hotspot_5[] == 3, na.rm = T) + ring.richness_N.5
ring.richness_N.95 <- sum(ring.richness_cat_map_rank_hotspot_5[] == 5, na.rm = T)
ring.richness_N.75 <- sum(ring.richness_cat_map_rank_hotspot_5[] == 4, na.rm = T) + ring.richness_N.95

ring.richness_N.quantiles <- c(ring.richness_N.5, ring.richness_N.25, ring.richness_N.75, ring.richness_N.95)
ring.richness_N.quantiles_perc <- round(ring.richness_N.quantiles/ring.richness_N.hotspots*100, 1)
ring.richness_labels <-  paste0(ring.richness_N.quantiles, "\n(", ring.richness_N.quantiles_perc, "%)")
quantiles <- c("Risk areas (5%)", "Risk areas (25%)", "Refuges (25%)", "Refuges (5%)")
ring.richness_N.quantiles_df_5 <- data.frame(cbind(ring.richness_N.quantiles, quantiles, ring.richness_labels))
ring.richness_N.quantiles_df_5$ring.richness_N.quantiles <- as.numeric(as.character(ring.richness_N.quantiles_df_5$ring.richness_N.quantiles))

saveRDS(ring.richness_N.quantiles_df_5, file = "./outputs/Threat_maps/Overlaps/ring.richness_N.quantiles_df_5.rds", version = "2")


# Extract data for mimicry rarity

ring.rarity_N.hotspots <- sum(ring.rarity_cat_map_rank_hotspot_5[] > 0, na.rm = T)
ring.rarity_N.5 <- sum(ring.rarity_cat_map_rank_hotspot_5[] == 2, na.rm = T)
ring.rarity_N.25 <- sum(ring.rarity_cat_map_rank_hotspot_5[] == 3, na.rm = T) + ring.rarity_N.5
ring.rarity_N.95 <- sum(ring.rarity_cat_map_rank_hotspot_5[] == 5, na.rm = T)
ring.rarity_N.75 <- sum(ring.rarity_cat_map_rank_hotspot_5[] == 4, na.rm = T) + ring.rarity_N.95

ring.rarity_N.quantiles <- c(ring.rarity_N.5, ring.rarity_N.25, ring.rarity_N.75, ring.rarity_N.95)
ring.rarity_N.quantiles_perc <- round(ring.rarity_N.quantiles/ring.rarity_N.hotspots*100, 1)
ring.rarity_labels <-  paste0(ring.rarity_N.quantiles, "\n(", ring.rarity_N.quantiles_perc, "%)")
quantiles <- c("Risk areas (5%)", "Risk areas (25%)", "Refuges (25%)", "Refuges (5%)")
ring.rarity_N.quantiles_df_5 <- data.frame(cbind(ring.rarity_N.quantiles, quantiles, ring.rarity_labels))
ring.rarity_N.quantiles_df_5$ring.rarity_N.quantiles <- as.numeric(as.character(ring.rarity_N.quantiles_df_5$ring.rarity_N.quantiles))

saveRDS(ring.rarity_N.quantiles_df_5, file = "./outputs/Threat_maps/Overlaps/ring.rarity_N.quantiles_df_5.rds", version = "2")


# Build df with all data

names(sp.richness_N.quantiles_df_5) <- names(sp.rarity_N.quantiles_df_5) <- names(MPD_N.quantiles_df_5) <- names(ring.richness_N.quantiles_df_5) <- names(ring.rarity_N.quantiles_df_5) <- c("N.quantiles", "quantiles", "labels")
all_indices_N.quantiles_df_5 <- rbind(sp.richness_N.quantiles_df_5,
                                      sp.rarity_N.quantiles_df_5,
                                      # MPD_N.quantiles_df_5,
                                      ring.richness_N.quantiles_df_5,
                                      ring.rarity_N.quantiles_df_5)
all_indices_N.quantiles_df_5$indices <- c(rep("sp.richness", 4),
                                          rep("sp.rarity", 4),
                                          # rep("MPD", 4),
                                          rep("ring.richness", 4),
                                          rep("ring.rarity", 4))
all_indices_N.quantiles_df_5$indices <- factor(all_indices_N.quantiles_df_5$indices, 
                                               # levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness"))
                                               levels = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))

saveRDS(all_indices_N.quantiles_df_5, file = "./outputs/Threat_maps/Overlaps/all_indices_N.quantiles_df_5.rds", version = "2")

# Plot

all_indices_N.quantiles_df_5 <- readRDS(file = "./outputs/Threat_maps/Overlaps/all_indices_N.quantiles_df_5.rds")

empty_rows <- data.frame(matrix(data = NA, nrow = 4, ncol = 4))
names(empty_rows) <- names(all_indices_N.quantiles_df_5)
# empty_rows[, 4] <- c("sp.richness", "sp.rarity", "MPD", "ring.richness")
empty_rows[, 4] <- c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity")
empty_rows[, 2] <- ""
empty_rows[, 1] <-  0
all_indices_N.quantiles_df_5 <- rbind(all_indices_N.quantiles_df_5, empty_rows)

# Without parenthesis in labels
all_indices_N.quantiles_df_5$labels <- str_remove(string = str_remove(string = all_indices_N.quantiles_df_5$labels, pattern = "\\("), pattern = "\\)")

str_remove(string = all_indices_N.quantiles_df_5$labels, pattern = "\\(.+\\)") # Remove all % lines
str_remove(string = all_indices_N.quantiles_df_5$labels, pattern = "[\\)|\\(]") # Does work, but should

# With only %
all_indices_N.quantiles_df_5$labels <- str_remove(string = all_indices_N.quantiles_df_5$labels, pattern = ".+[\n]")

# With group names on two lines
all_indices_N.quantiles_df_5$quantiles <- as.factor(rep(c("Risk areas\n(5%)", "Risk areas\n(25%)", "Refuges\n(25%)", "Refuges\n(5%)"), 5))

# Remove labels overlapped by the 25% and 5% lines
all_indices_N.quantiles_df_5$labels[c(2, 13)] <- ""

pdf(file = "./graphs/Overlaps/areas_overlap_quantiles_hotspot_5.pdf", height = 8, width = 12)

# levels(all_indices_N.quantiles_df_5$quantiles)

# Play with witdh and the extra empty level for N.quantiles to get a space inbetween group of bars

# internal_margins <- par()$mar
# par(mar = c(3.1, 3.5, 3.5, 2.1))
# external_margins <- par()$oma
# par(oma = c(0, 0, 3, 0))

g2 <- ggplot(all_indices_N.quantiles_df_5, aes(x = quantiles, y = N.quantiles, fill = indices)) +
  geom_col(position = position_dodge(width = 0.9), 
           # width = c(0.9,0.9,0.9,0.9,
           #           0.9,0.9,0.9,0.9,
           #           0.9,0.9,0.9,0.9,
           #           0.9,0.9,0.9,0.9, 
           #           0.1,0.1,0.1,0.1),
           show.legend = T) +
  labs(title = "", 
       y = "Communities", x = NULL) +
  # scale_x_discrete(limits = c("Risk areas (5%)", "Risk areas (25%)", "Refuges (25%)", "Refuges (5%)")) +
  scale_x_discrete(limits = c("Risk areas\n(25%)", "Refuges\n(25%)", "", "Risk areas\n(5%)", "Refuges\n(5%)")) +
  scale_fill_manual(name = "Indices",
                    # values = c("sp.richness" = "#FF4040", "sp.rarity" = "#80C2FF", "MPD" = "#FFD84D", "ring.richness" = "#7ACC29"),
                    # labels = c("sp.richness" = "Species richness", "sp.rarity" = "Species rarity", "MPD" = "MPD", "ring.richness" = "Mimicry richness")) +
                    values = c("sp.richness" = "#FF4040", "sp.rarity" = "#80C2FF", "ring.richness" = "#FFD84D", "ring.rarity" = "#7ACC29"),
                    labels = c("sp.richness" = "Species richness", "sp.rarity" = "Species rarity", "ring.richness" = "Mimicry richness", "ring.rarity" = "Mimicry rarity")) +
  # theme_classic() +
  scale_y_continuous(limits = c(0, 600), expand = c(0, 0), breaks = seq(0, 600, 200)) +
  # ylim = c(0, 3500) +
  # geom_rangeframe(# data = limits_df_5, aes(x = x, y = y),
  #                 size = 1.4, sides = "l") +
  theme(panel.background = element_rect(fill = NA),
        legend.position = c(0.82, 0.57),
        # legend.key = element_rect(fill = "white", colour = "white"),
        legend.title = element_text(size = 24, vjust = 3, face = "bold"),
        legend.text = element_text(size = 22, face = "plain"),
        # margin = margin(t = 0, unit = "pt")),
        legend.key.size = unit(2.5, 'lines'),
        # legend.spacing.y = unit(1,"cm"),
        # axis.line.x = element_line(color = NA),
        axis.line.y = element_line(color = "black", size = 1.3),
        axis.ticks = element_line(color = "black"),
        axis.ticks.x = element_line(size = NA),
        axis.ticks.y = element_line(size = 1.2),
        axis.ticks.length.y = unit(8, "pt"),
        axis.text = element_text(size = 22, face = "bold", color = "black"),
        axis.text.x = element_text(vjust = 0, size = 24, margin = margin(t = 8, b = 15)),
        axis.title = element_text(size = 28, face = "bold", color = "black"),
        # axis.title.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title.y = element_text(size = 28, margin = margin(l = 5, r = 10))) +
  geom_segment(aes(x = 3.5, y = 53.55, xend = 5.5, yend = 53.55), lwd = 1, color = "grey30", lty = "dashed") +
  geom_segment(aes(x = 0.5, y = 267.75, xend = 2.5, yend = 267.75), lwd = 1, color = "grey30", lty = "dashed") +
  
  # geom_segment(aes(x = 0.5, y = 1338.5, xend = 0.75, yend = 1338.5), lwd = 1, color = "grey30", lty = "dashed") +
  # geom_segment(aes(x = 1, y = 1338.5, xend = 1.55, yend = 1338.5), lwd = 1, color = "grey30", lty = "dashed") +
  # geom_segment(aes(x = 1.8, y = 1338.5, xend = 2.0, yend = 1338.5), lwd = 1, color = "grey30", lty = "dashed") +
  # geom_segment(aes(x = 2.48, y = 1338.5, xend = 2.5, yend = 1338.5), lwd = 1, color = "grey30", lty = "dashed") +
  
  geom_text(aes(label = labels), vjust=-0.8, size = 6, fontface = "bold", position = position_dodge(width = 0.9)) +
  
  annotate("text", x = 0.665, y = 290, size = 6, fontface = "bold", label = "23.6%") +
  annotate("text", x = 4.380, y = 94, size = 6, fontface = "bold", label = "6.7%") +
   
  annotate("text", x = 3, y = 585, size = 11, fontface = "bold", label = "Top 5% hotspots") +
  annotate("text", x = 0.7, y = 580, size = 11, fontface = "bold", label = "(b)") +
  
  coord_cartesian(clip = "off") # To allow to draw label outside of the plot area

print(g2)

par(mar = internal_margins, oma = external_margins)

dev.off()


### 13/ Bar plots for proportions of overlap with hotspots 25% and 5% aside ####

library(ggpubr)

pdf(file = "./graphs/Overlaps/areas_overlap_quantiles_hotspot_25_5_Vertical.pdf", height = 15, width = 11.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,4.5,2.1))
par(mfrow = c(2,1))

ggarrange(g1, g2, # List of plots
          ncol = 1, nrow = 2)

par(mar = internal_margins)
par(mfrow = c(1,1))

dev.off()

pdf(file = "./graphs/Overlaps/areas_overlap_quantiles_hotspot_25_5_Horizontal.pdf", height = 9, width = 23)

internal_margins <- par()$mar
par(mar = c(2.1,3.5,4.5,2.1))
par(mfrow = c(1,2))

ggarrange(g1, g2, # List of plots
          ncol = 2, nrow = 1)

par(mar = internal_margins)
par(mfrow = c(1,1))

dev.off()



### 14/ Randomisation test for overlaps ####

# 14.1/ Hotspots 25%

# Result of the test depend on sample size, so directly of community size.
N.hotspots <- sum(sp.richness_mask_rank_25[] > 0, na.rm = T)
N.hotspots <- 20000
N.hotspots <- 200

N.com <- N.hotspots*4

sampled_overlap <- rep(NA, 1000)
for (i in 1:999) {
  sample_i <- sample(x = c(rep(0, N.com-N.hotspots),
                           rep(1, N.hotspots)),
                     size = N.hotspots)
  
  overlap_i <- round(sum(sample_i)/N.hotspots, 3)*100
  
  sampled_overlap[i] <- overlap_i
}

hist(sampled_overlap)
