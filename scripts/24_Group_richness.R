
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
bg_mask <- readRDS(file = "./input_data/Map_stuff/bg_mask.rds")

# Stack of species probabilities
All_sp_proba_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))

### 2/ Define groups ####

phylo.Ithomiini <- readRDS(file = "./input_data/Phylogenies/Final_phylogeny.rds")

pdf(file = "./input_data/phylo.Ithomiini.pdf", height = 50, width = 10)
plot(phylo.Ithomiini)
dev.off()

load(file = paste0("./input_data/list.sp.RData"))

Backbone_genera <- c("Patricia", "Athesis", "Methona", "Elzunia", "Tithorea", "Aeria", "Scada", "Sais", "Mechanitis", "Forbestra", "Thyridia", "Olyras", "Paititia", "Eutresis", "Melinaea", "Athyrtis")

Backbone_sp <- str_subset(string = list.sp$Sp_full, pattern = paste(Backbone_genera, collapse="|"))
Coregroup_sp <- setdiff(list.sp$Sp_full, Backbone_sp)


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


Backbone.sp.richness_Jaccard.80 <- readRDS(file = "./outputs/Indices_maps/By_groups/Backbone.sp.richness_Jaccard.80.rds")
Coregroup.sp.richness_Jaccard.80 <- readRDS(file = "./outputs/Indices_maps/By_groups/Coregroup.sp.richness_Jaccard.80.rds")

prop_Backbone.sp.richness_Jaccard.80 <- readRDS(file = "./outputs/Indices_maps/By_groups/prop_Backbone.sp.richness_Jaccard.80.rds")


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
legend(legend = "(a)", x = "bottomright", bty = "n",
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
legend(legend = "(b)", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))

par(mfrow = c(1,1))
par(mar = internal_margins)

dev.off()

# Try to find ways to merge those base plots with the Venn diagram (grob/gList object) for facet d in the multiple plot ####

# pdf(file = paste0("./maps/Indices_maps/By_groups/test.pdf"), height = 5.3, width = 13)
# 
# internal_margins <- par()$mar
# par(mar = c(3.1,3.5,3.5,2.1))
# par(mfrow = c(1,2))
# 
# cowplot::plot_grid(venn_plot_2, NULL,
#                    ncol = 2, nrow = 1,
#                    labels = c(NULL, "(b)"),
#                    label_x = 0.95, label_y = 0.05,
#                    label_size = 18)
# 
# # A/ Backbone group
# 
# plot_backbone_sp_richness <- function()
# {
#    image(Backbone.sp.richness_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Species richness for Backbone group"), 
#          cex.axis = 1.4, cex.main = 1.6, axis.args=list(cex.axis=1.4),
#          ylab = "", xlab = "",
#          legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
#          legend  = F )
#    plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
#    # plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
#    # plot(country_borders, lwd = 0.2, border = "grey20", color = "", add = T)
#    plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
#    scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
#    prettymapr::addnortharrow(scale = 0.55, padin = c(0.1, 0.1), text.col = "#00000000")
#    addRasterLegend(Backbone.sp.richness_Jaccard.80, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
#    addRasterLegend(Backbone.sp.richness_Jaccard.80, locs = seq(0, 20, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
#    graphics::text(x = -110, y = 4, font = 2, cex = 1.1, label = "Species")
#    legend(legend = "(a)", x = "bottomright", bty = "n",
#           text.font = 2, cex = 1.8, inset=c(0, 0))
# }
# 
# library("grid")
# library("ggplotify")
# 
# test <- as.grob(plot_backbone_sp_richness)
# grid.newpage()
# grid.draw(test)
# 
# ggtest <- as.ggplot(plot_backbone_sp_richness)
# 
# print(ggtest)
# 
# 
# 
# # # grid.newpage()
# # grid.draw(venn_plot_2_backbone)
# # grid.text(label = expression(bold("Backbone")), x = unit(0.5, "npc"), y = unit(0.85, "npc"), gp = gpar(fontsize = 20))
# # grid.newpage()
# 
# par(mfrow = c(1,1))
# par(mar = internal_margins)
# 
# dev.off()
# 
# ?pdf
# 
# pdf(file = paste0("./maps/Indices_maps/By_groups/test.pdf"), height = 5.3, width = 13)
# 
# cowplot::plot_grid(test, venn_plot_2,
#                    ncol = 2, nrow = 1,
#                    labels = c("", "(b)"),
#                    label_x = 0.95, label_y = 0.05,
#                    label_size = 18)
# 
# dev.off()
# 
# 
# library("multipanelfigure")
# 
# ?multi_panel_figure
# ?capture_base_plot
# ?fill_panel
# 
# figure_1 <- multi_panel_figure(rows = 1, columns = 2, height = 5.3, width = 10, unit = "inches", panel_label_type = "none")
# 
# figure_1 <- fill_panel(figure = figure_1, panel = venn_plot_2, row = 1, column = 1)
# figure_1 <- fill_panel(figure = figure_1, panel = capture_base_plot(plot_backbone_sp_richness()), row = 1, column = 2)
# print(figure_1)
# 
# save_multi_panel_figure(figure = figure_1, filename = "./maps/Indices_maps/By_groups/test.pdf")
# 
# pdf(file = paste0("./maps/Indices_maps/By_groups/test.pdf"), height = 5.3, width = 13)
# print(figure_1)
# dev.off()
# 
# pdf(file = paste0("./maps/Indices_maps/By_groups/test.pdf"), height = 5.3, width = 13)
# par(mfrow = c(1,2))
# 
# # plot.new()
# grid.draw(venn_plot_2)
# plot_backbone_sp_richness()
# 
# par(mfrow = c(1,1))
# dev.off()
# 
# plot(x = NULL)
# 

### 5/ Functions for final multiple plot with graticules and Mollweide projections ####

# 5.1/ Mollweide projection ####

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
   
   # Save new object in outputs folder
   # save(eval(parse(text = paste0(x_name, "_Mollweide"))), file = paste0("./outputs/Indices_maps/mollweide_projections/",x_name, "_Mollweide.RData"))
   saveRDS(new_map, file = paste0("./outputs/Indices_maps/mollweide_projections/",x_name, "_Mollweide.rds"))
}

Mollweide_shp_projection <-  function(x) # Shp to project
{
   x_name <- deparse(substitute(x)) # Get the name of the initial shp as a character string
   
   new_shp <- spTransform(x, CRSobj = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
   
   # Generate new object with "_Mollweide" suffix in the global environment
   eval(call("<<-", as.name(paste0(x_name, "_Mollweide")), new_shp))
}


# 5.2/ Plotting function ####

map_indices_Mollweide <- function(x,                                    # Raster to map
                                  color_palette = pal_bl_red_Mannion,   # Color palette
                                  main_title,                           # Main title
                                  main_title_cex = 1.4,                 # Main title size
                                  
                                  xlim = c(-4600, 4600),   # Limit of plot on x-axis (Longitude)
                                  ylim = c(-4450, 3400),    # Limit of plot on y-axis (Latitude)
                                  axis_cex = 1,             # Axes size
                                  
                                  xlab = "",                # X-axis label
                                  ylab = "",                # Y-axis label
                                  x_axis_breaks = c(-3930, -2170, -420, 1500, 3050),            # X-axis tick breaks
                                  y_axis_breaks = c(-3650, -2450, -1220, 0, 1230, 2445),        # Y-axis tick breaks
                                  x_axis_labels = c("120°E", "100°E", "80°E", "60°E", "40°E"),      # X-axis tick labels
                                  y_axis_labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"),  # Y-axis tick labels
                                  
                                  legend_title,             # Legend title
                                  legend_title_cex = 1.1,   # Legend title size
                                  legend_title_x = -3550,   # Legend title x position
                                  legend_title_y = 430,     # Legend title y position
                                  legend_cex = 1.2,         # Legend size
                                  legend_breaks,            # Legend tick positions
                                  legend_location = c(-4100, -3800, -3950, 0),  # Legend position
                                  
                                  scale_bar_position = c(-2600, -4000),  # Scale bar position
                                  
                                  arrow_scale = 0.45,           # North arrow size
                                  arrow_padin = c(0.15, 0.15),  # North arrow position adjustement
                                  
                                  facet_letter = "",                  # Small case letter for facet
                                  facet_letter_cex = 1.6,             # Size of small case letter for facet
                                  facet_letter_inset = c(0, -0.008))  # Position adjustment of small case letter for facet

{
   # Plot raster background without axis
   image(x, col = color_palette,
         xlim = xlim, ylim = ylim, axes = F,
         xlab = xlab, ylab = ylab)
   title(main = main_title, cex.main = main_title_cex, line = 1)
   
   # Generate axes with manual positioning of ticks
   axis(side = 1, at = x_axis_breaks, labels = x_axis_labels, cex.axis = axis_cex, lwd = 0.2, lwd.ticks = 1)
   axis(side = 2, at = y_axis_breaks, labels = y_axis_labels, cex.axis = axis_cex, lwd = 0.2, lwd.ticks = 1)
   
   # Add background, borders and graticules
   plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
   plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
   plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
   
   # Add scale bar in legend
   scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = scale_bar_position, label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
   prettymapr::addnortharrow(scale = arrow_scale, padin = arrow_padin, text.col = "#00000000")
   rangeBuilder::addRasterLegend(x, locs = legend_breaks, cex.axis = legend_cex, ramp = color_palette, ncolors = 200, border = T, location = legend_location)
   rangeBuilder::addRasterLegend(x, locs = legend_breaks, cex.axis = legend_cex, ramp = color_palette, ncolors = 200, border = T, location = legend_location)
   graphics::text(x = legend_title_x, y = legend_title_y, font = 2, cex = legend_title_cex, label = legend_title)
   
   # Add facet letter
   legend(legend = facet_letter, x = "bottomright", bty = "n",
          text.font = 2, cex = facet_letter_cex, inset = facet_letter_inset)
   
}

# map_indices_Mollweide(x = sp.richness_Mollweide,
#                       main_title = "Species richness",
#                       legend_title = "Species",
#                       legend_breaks = seq(0, 120, 20), 
#                       facet_letter = "a")

# 5.3/ Add continental null values to raster ####

add_continental_null_values <- function(x)
{
   continent_mask <- readRDS(file = paste0("./input_data/Map_stuff/continent_mask.rds"))
   
   y <- continent_mask  # Create final new raster from continental mask
   y[!is.na(x[])] <- x[!is.na(x[])]  # Add initial raster values
   
   return(y)
}

# 5.4/ Contrasting raster ####
contrasting_raster <- function(x, zmin, zmax)
{
   Community_mask <- readRDS(file = "./input_data/Map_stuff/Community_mask.rds")
   continent_mask <- readRDS(file = paste0("./input_data/Map_stuff/continent_mask.rds"))
   
   x[x[] <= zmin] <- zmin  # Fix low values
   x[x[] >= zmax] <- zmax  # Fix high values
   
   x <- mask(x, mask = Community_mask)  # Cut out values that are outside Ithomiini range
   
   y <- continent_mask + zmin  # Create final new raster from continental mask with baseline = zmin
   y[!is.na(x[])] <- x[!is.na(x[])]  # Add initial raster values
   
   return(y)
}


### 6/ Final multiple plot with graticules and Mollweide projections ####

# New color palette
pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")

# Load mask for continent borders, plot border, and grid
grid_Mollweide_out <- readRDS(file = "./input_data/Map_stuff/grid_Mollweide_out.rds")
large_bg_mask_Mollweide <- readRDS(file = "./input_data/Map_stuff/large_bg_mask_Mollweide.rds")
bbox_sp_Mollweide <- readRDS(file = "./input_data/Map_stuff/bbox_sp_Mollweide.rds")
Community_mask <- readRDS(file = "./input_data/Map_stuff/Community_mask.rds")
continent_mask <- readRDS(file = paste0("./input_data/Map_stuff/continent_mask.rds"))

# 6.1/ Add continental null values

prop_Backbone.sp.richness_Jaccard.80 <- add_continental_null_values(prop_Backbone.sp.richness_Jaccard.80)

plot(prop_Backbone.sp.richness_Jaccard.80)

# 6.2/ Contrast rasters

hist(Backbone.sp.richness_Jaccard.80)
hist(Coregroup.sp.richness_Jaccard.80)
hist(prop_Backbone.sp.richness_Jaccard.80)

prop_Backbone.sp.richness_Jaccard.80_contrasted <- contrasting_raster(x = prop_Backbone.sp.richness_Jaccard.80, zmin = 0, zmax = 60)

# 6.3/ Project into Mollweide projection

Mollweide_projection(Backbone.sp.richness_Jaccard.80)
Mollweide_projection(Coregroup.sp.richness_Jaccard.80)
Mollweide_projection(prop_Backbone.sp.richness_Jaccard.80)
Mollweide_projection(prop_Backbone.sp.richness_Jaccard.80_contrasted)

# Northern_Andes_shp <- as(Northern_Andes_shp, 'Spatial')
# Mollweide_shp_projection(Northern_Andes_shp)
# Central_Andes_shp <- as(Central_Andes_shp, 'Spatial')
# Mollweide_shp_projection(Central_Andes_shp)
# Western_Amazon_shp <- as(Western_Amazon_shp, 'Spatial')
# Mollweide_shp_projection(Western_Amazon_shp)

# 6.4/ Load and build the Venn diagram for facet d (Script_32_Subset_analyses)

load(file = paste0("./outputs/Subset_analyses/Andes_mask_lists.RData"))
Western_Amazon_mask_sp_list <- readRDS(file = paste0("./outputs/Subset_analyses/Western_Amazon_mask_sp_list.rds"))

venn_plot_2 <- draw.pairwise.venn(
   area1 = length(Western_Amazon_mask_sp_list),
   area2 = length(Andes_mask_sp_list),
   cross.area = length(intersect(Western_Amazon_mask_sp_list, Andes_mask_sp_list)),
   category = c("Upper\nAmazon", "Central &\nNorthern\nAndes"),
   fontfamily = "sans",
   col = c("blue", "tan4"),
   fill = c("dodgerblue", "bisque3"),
   alpha = 0.3,
   lty = 1,
   cex = 1.3,
   cat.cex = 1.3,
   cat.pos = c(-55, 60),
   cat.col = c("blue", "tan4"),
   cat.dist = c(0.08, 0.08),
   cat.fontfamily = "sans",
   ext.pos = 30,
   ext.dist = -0.07,
   ext.length = 0.60,
   ext.line.lwd = 2,
   ext.line.lty = 1,
   label.col = c("blue", "black", "tan4"),
   print.mode = c("raw", "percent"),
   sigdigs = 2,
   margin = 0.01
)
grid.newpage()
par(mfrow = c(2,2))
vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)
grid.draw(venn_plot_2)
grid.text(label = expression(bold("Species list overlaps across regions")), x = unit(0.5, "npc"), y = unit(1.2, "npc"), gp = gpar(fontsize = 15))
grid.text(label = expression(bold("(d)")), x = unit(0.95, "npc"), y = unit(0.05, "npc"), gp = gpar(fontsize = 17))


# 6.5/ Plot maps

# To add shp of Andes and Western Amazon
#
# Northern_Andes_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Northern_Andes_raw_shp.rds")
# Central_Andes_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Central_Andes_raw_shp.rds")
# Western_Lowlands_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Western_Lowlands_raw_shp.rds")
# Western_Amazon_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Western_Amazon_raw_shp.rds")
# 
# plot(Northern_Andes_shp)
# Northern_Andes_shp2_Mollweide <- remove.holes(gSimplify(Northern_Andes_shp_Mollweide, tol = 0.2, topologyPreserve = F))
# plot(Northern_Andes_shp2)
# 
# Central_Andes_shp2_Mollweide <- remove.holes(gSimplify(Central_Andes_shp_Mollweide, tol = 30.0, topologyPreserve = F))
# plot(Central_Andes_shp2_Mollweide)
# 
# Western_Amazon_shp2_Mollweide <- remove.holes(gSimplify(Western_Amazon_shp_Mollweide, tol = 60.0, topologyPreserve = F))
# plot(Western_Amazon_shp2_Mollweide)
# 
# plot(Western_Amazon_shp2_Mollweide, col = NA, border = "dodgerblue", lwd = 2, add = T)
# plot(Northern_Andes_shp2_Mollweide, col = NA, border = "tan4", lwd = 2, add = T)
# plot(Central_Andes_shp2_Mollweide, col = NA, border = "snow4", lwd = 2, add = T)
# 
# graphics::text(x = -570, y = 1500, font = 2, cex = 1.3, label = "CA", col = "#F8766D")
# graphics::text(x = -770, y = 10, font = 2, cex = 1.3, label = "EE", col = "#7CAE00")
# graphics::text(x = -630, y = -1340, font = 2, cex = 1.3, label = "EP", col = "#00BFC4")
# 
# segments(x0 = -610, y0 = 1350 , x1 = -700, y1 = 1100, col = "#F8766D", lwd = 4)
# segments(x0 = -620, y0 = -20 , x1 = -400, y1 = -100, col = "#7CAE00", lwd = 4)
# segments(x0 = -490, y0 = -1320 , x1 = -250, y1 = -1150, col = "#00BFC4", lwd = 4)


library(RGraphics)
library(gridExtra)
library(gridBase)

pdf(file = paste0("./maps/Indices_maps/By_groups/Backbone_vs_Coregroup_richness_&_prop.pdf"), height = 8, width = 8)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
par(mfrow = c(2,2))

map_indices_Mollweide(x = Backbone.sp.richness_Jaccard.80_Mollweide, color_palette = pal_bl_red_Mannion, 
                      main_title = "Species richness for Backbone group",
                      legend_title = "Species",
                      legend_breaks = seq(0, 20, 5), 
                      facet_letter = "(a)")

map_indices_Mollweide(x = Coregroup.sp.richness_Jaccard.80_Mollweide, color_palette = pal_bl_red_Mannion, 
                      main_title = "Species richness for Core-group",
                      legend_title = "Species",
                      legend_breaks = seq(0, 100, 20), 
                      facet_letter = "(b)")

map_indices_Mollweide(x = prop_Backbone.sp.richness_Jaccard.80_contrasted_Mollweide, color_palette = pal_bl_red_Mannion, 
                      main_title = "Proportion of Backbone group species",
                      legend_title_x = -2850,
                      legend_title_y = 600,
                      legend_title = "Proportion (%)",
                      legend_breaks = seq(0, 60, 15), 
                      facet_letter = "(c)")

# Facet (d) to complete with Venn Diagram from Script 32_Subset_analyses
plot.new()

vps <- baseViewports()
pushViewport(vps$inner, vps$figure, vps$plot)

grid.draw(venn_plot_2)
grid.text(label = expression(bold("Species list overlaps across regions")), x = unit(0.5, "npc"), y = unit(1.07, "npc"), gp = gpar(fontsize = 14))
grid.text(label = expression(bold("(d)")), x = unit(0.95, "npc"), y = unit(0.05, "npc"), gp = gpar(fontsize = 17))

popViewport(0)


par(mfrow = c(1,1))
par(mar = internal_margins)

dev.off()







