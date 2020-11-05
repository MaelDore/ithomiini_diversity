
###### Script 24: Richness maps: core-group vs. backbone group ######

# Only for Jaccard.80

###

# Inputs: 
# Probability map of species
# Phylogeny

# Outputs:
# Richness maps for core-group and backbone group

###


### 1/ Load stuff ####

# Clean environment
rm(list = ls())

library(raster)
library(ape)
library(picante)
library(geiger)
library(treespace)

library(rangeBuilder)

pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")
crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_15.rds"))


# Stack of species probabilities
All_sp_proba_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))

### 2/ Define groups ####

phylo.Ithomiini <- readRDS(file = "./input_data/Phylogenies/Final_phylogeny.rds")

pdf(file = "./input_data/phylo.Ithomiini.pdf", height = 50, width = 10)
plot(phylo.Ithomiini)
dev.off()

load(file = paste0("./input_data/list.sp.RData"))
list.sp <- list.sp$Sp_full

Backbone_genera <- c("Patricia", "Athesis", "Methona", "Elzunia", "Tithorea", "Aeria", "Scada", "Sais", "Mechanitis", "Forbestra", "Thyridia", "Olyras", "Paititia", "Eutresis", "Melinaea", "Athyrtis")

Backbone_sp <- str_subset(string = list.sp, pattern = paste(Backbone_genera, collapse="|"))
Coregroup_sp <- setdiff(list.sp, Backbone_sp)


### 3/ Extract layers and map richness

Backbone_proba_stack_Jaccard.80 <- subset(All_sp_proba_stack_Jaccard.80, subset = which(names(All_sp_proba_stack_Jaccard.80) %in% Backbone_sp))
Coregroup_proba_stack_Jaccard.80 <- subset(All_sp_proba_stack_Jaccard.80, subset = setdiff(1:388, which(names(All_sp_proba_stack_Jaccard.80) %in% Backbone_sp)))

saveRDS(Backbone_proba_stack_Jaccard.80, file = "./outputs/Indices_maps/By_groups/Backbone_proba_stack_Jaccard.80.rds", version = "2")
saveRDS(Coregroup_proba_stack_Jaccard.80, file = "./outputs/Indices_maps/By_groups/Coregroup_proba_stack_Jaccard.80.rds", version = "2")

Backbone.sp.richness_Jaccard.80 <- readAll(calc(Backbone_proba_stack_Jaccard.80, fun = sum))
Coregroup.sp.richness_Jaccard.80 <- readAll(calc(Coregroup_proba_stack_Jaccard.80, fun = sum))

saveRDS(Backbone.sp.richness_Jaccard.80, file = "./outputs/Indices_maps/By_groups/Backbone.sp.richness_Jaccard.80.rds", version = "2")
saveRDS(Coregroup.sp.richness_Jaccard.80, file = "./outputs/Indices_maps/By_groups/Coregroup.sp.richness_Jaccard.80.rds", version = "2")

# Proportion of backbone group species
prop_Backbone.sp.richness_Jaccard.80 <- Backbone.sp.richness_Jaccard.80/(Backbone.sp.richness_Jaccard.80 + Coregroup.sp.richness_Jaccard.80) * 100
saveRDS(prop_Backbone.sp.richness_Jaccard.80, file = "./outputs/Indices_maps/By_groups/prop_Backbone.sp.richness_Jaccard.80.rds", version = "2")


### 4/ Plot richness maps ####

# 4.1/ Backbone richness map ####

pdf(file = paste0("./maps/Indices_maps/By_groups/Backbone.sp.richness_Jaccard.80.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(Backbone.sp.richness_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Species richness for Backbone group"), 
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
addRasterLegend(Backbone.sp.richness_Jaccard.80, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.1, label = "Species")

par(mar = internal_margins)

dev.off()

# 4.2/ Core-group richness map ####

pdf(file = paste0("./maps/Indices_maps/By_groups/Coregroup.sp.richness_Jaccard.80.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(Coregroup.sp.richness_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Species richness for Core-group"), 
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
addRasterLegend(Coregroup.sp.richness_Jaccard.80, locs = seq(0, 100, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.1, label = "Species")

par(mar = internal_margins)

dev.off()

# 4.3/ Proportion of Backbone group species ####

prop_Backbone.sp.richness_Jaccard.80[prop_Backbone.sp.richness_Jaccard.80 > 60] <- 60

pdf(file = paste0("./maps/Indices_maps/By_groups/prop_Backbone.sp.richness_Jaccard.80.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(prop_Backbone.sp.richness_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Proportion of Backbone group species"), 
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
addRasterLegend(prop_Backbone.sp.richness_Jaccard.80, locs = seq(0, 60, 15), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -109, y = 4, font = 2, cex = 1.1, label = "Proportion (%)")

par(mar = internal_margins)

dev.off()

hist(prop_Backbone.sp.richness_Jaccard.80)

# 4.4/ Multiple plot ####

pdf(file = paste0("./maps/Indices_maps/By_groups/Backbone_vs_Coregroup.sp.richness_Jaccard.80.pdf"), height = 5.3, width = 13)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
par(mfrow = c(1,2))

# A/ Backbone group
image(Backbone.sp.richness_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Species richness for Backbone group"), 
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
addRasterLegend(Backbone.sp.richness_Jaccard.80, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.1, label = "Species")
legend(legend = "a", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))

# B/ Core-group
image(Coregroup.sp.richness_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Species richness for Core-group"), 
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
addRasterLegend(Coregroup.sp.richness_Jaccard.80, locs = seq(0, 100, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.1, label = "Species")
legend(legend = "b", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))

par(mfrow = c(1,1))
par(mar = internal_margins)

dev.off()

