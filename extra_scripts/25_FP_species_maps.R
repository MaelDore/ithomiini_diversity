
###### Script 25: Species Fair-Proportion maps ######

# Map the absolute and proportional richness of communities in the most/least phylogenetically original species (based on species FP)
# Investigate the relationship between geographic rarity and phylogenetic originality

# Only for Jaccard.80

###

# Inputs: 
   # Species FP indices from script 14
   # Species probability stack from script 12

# Outputs:
   # Top and Low 25% FP species richness
      # Absolute and proportion
   # Scatter plot species FP ~ species rarity

###


### 1/ Load stuff ####

# Clean environment
rm(list = ls())

library(raster)

pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")
crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_15.rds"))
bg_mask <- readRDS(file = "./input_data/Map_stuff/bg_mask.rds")

# Load table of ES and FP indices of the 339 species in the phylogeny
load(file = "./outputs/Indices_Maps/Species_ES_FP_index_table.RData")
Phylo.sp_list <- as.character(Species_ES_FP_index_table$Sp.full)

# Stack of species probabilities
All_sp_proba_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))


### 2/ Extract layers  ####

Phylo_sp_proba_stack_Jaccard.80 <- All_sp_proba_stack_Jaccard.80[[Phylo.sp_list]]
saveRDS(Phylo_sp_proba_stack_Jaccard.80, file = "./outputs/Indices_maps/By_groups/Phylo_sp_proba_stack_Jaccard.80.rds", version = "2")

Phylo_sp.richness_Jaccard.80 <- readAll(calc(Phylo_sp_proba_stack_Jaccard.80, fun = sum))

# Reorder rows by species FP
Species_ES_FP_index_table_ordered <- Species_ES_FP_index_table %>% 
  arrange(FP)

# Nb of species to extract for Top and Low 25%
N.sp <- round(nrow(Species_ES_FP_index_table_ordered)/4, 0)

Top_25 <- as.character(Species_ES_FP_index_table_ordered$Sp.full[1:N.sp])
Low_25 <- as.character(Species_ES_FP_index_table_ordered$Sp.full[(nrow(Species_ES_FP_index_table_ordered)-N.sp+1):nrow(Species_ES_FP_index_table_ordered)])

Top_25_sp_proba_stack_Jaccard.80 <- All_sp_proba_stack_Jaccard.80[[Top_25]]
Low_25_sp_proba_stack_Jaccard.80 <- All_sp_proba_stack_Jaccard.80[[Low_25]]
saveRDS(Top_25_sp_proba_stack_Jaccard.80, file = "./outputs/Indices_maps/By_groups/Top_25_sp_proba_stack_Jaccard.80.rds", version = "2")
saveRDS(Low_25_sp_proba_stack_Jaccard.80, file = "./outputs/Indices_maps/By_groups/Low_25_sp_proba_stack_Jaccard.80.rds", version = "2")


### 3/ Compute richness maps ####

# Absolute maps
Top_25_FP_sp.richness_Jaccard.80 <- readAll(calc(Top_25_sp_proba_stack_Jaccard.80, fun = sum))
Low_25_FP_sp.richness_Jaccard.80 <- readAll(calc(Low_25_sp_proba_stack_Jaccard.80, fun = sum))

saveRDS(Top_25_FP_sp.richness_Jaccard.80, file = "./outputs/Indices_maps/By_groups/Backbone.sp.richness_Jaccard.80.rds", version = "2")
saveRDS(Low_25_FP_sp.richness_Jaccard.80, file = "./outputs/Indices_maps/By_groups/Coregroup.sp.richness_Jaccard.80.rds", version = "2")

# Relative maps
Top_25_FP_prop_sp.richness_Jaccard.80 <- Top_25_FP_sp.richness_Jaccard.80/Phylo_sp.richness_Jaccard.80*100
Low_25_FP_prop_sp.richness_Jaccard.80 <- Low_25_FP_sp.richness_Jaccard.80/Phylo_sp.richness_Jaccard.80*100

saveRDS(Top_25_FP_prop_sp.richness_Jaccard.80, file = "./outputs/Indices_maps/By_groups/Top_25_FP_prop_sp.richness_Jaccard.80.rds", version = "2")
saveRDS(Low_25_FP_prop_sp.richness_Jaccard.80, file = "./outputs/Indices_maps/By_groups/Low_25_FP_prop_sp.richness_Jaccard.80.rds", version = "2")

### 4/ Plot maps ####

# 4.1/ Top 25% FP map ####

# Absolute
pdf(file = paste0("./maps/Indices_maps/By_groups/Top_25_FP_sp.richness_Jaccard.80.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(Top_25_FP_sp.richness_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Species richness for most phylogenetically original species"), 
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
addRasterLegend(Top_25_FP_sp.richness_Jaccard.80, locs = seq(0, 35, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1.1, label = "Species")

par(mar = internal_margins)

dev.off()

# Proportion
pdf(file = paste0("./maps/Indices_maps/By_groups/Top_25_FP_prop_sp.richness_Jaccard.80.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

Top_25_FP_prop_sp.richness_Jaccard.80_contrasted <- Top_25_FP_prop_sp.richness_Jaccard.80
Top_25_FP_prop_sp.richness_Jaccard.80_contrasted[Top_25_FP_prop_sp.richness_Jaccard.80 > 50] <- 50

image(Top_25_FP_prop_sp.richness_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Proportion of most phylogenetically original species"), 
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
addRasterLegend(Top_25_FP_prop_sp.richness_Jaccard.80_contrasted, locs = seq(0, 50, 10), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.1, label = "Proportion (%)")

par(mar = internal_margins)

dev.off()

# 4.2/ Low 25% FP map ####

# Absolute
pdf(file = paste0("./maps/Indices_maps/By_groups/Low_25_FP_sp.richness_Jaccard.80.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

image(Low_25_FP_sp.richness_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Species richness for least phylogenetically original species"), 
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
addRasterLegend(Low_25_FP_sp.richness_Jaccard.80, locs = seq(0, 35, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1.1, label = "Species")

par(mar = internal_margins)

dev.off()

# Proportion
pdf(file = paste0("./maps/Indices_maps/By_groups/Low_25_FP_prop_sp.richness_Jaccard.80.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

Low_25_FP_prop_sp.richness_Jaccard.80_contrasted <- Low_25_FP_prop_sp.richness_Jaccard.80
Low_25_FP_prop_sp.richness_Jaccard.80_contrasted[Low_25_FP_prop_sp.richness_Jaccard.80 > 80] <- 80

image(Low_25_FP_prop_sp.richness_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Proportion of least phylogenetically original species"), 
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
addRasterLegend(Low_25_FP_prop_sp.richness_Jaccard.80_contrasted, locs = seq(0, 80, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.1, label = "Proportion (%)")

par(mar = internal_margins)

dev.off()

# 4.3/ Multiple plot

pdf(file = paste0("./maps/Indices_maps/By_groups/Species_FP_all_maps_Jaccard.80.pdf"), height = 11, width = 13)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
par(mfrow = c(2,2))

# A/ Top 25% - Absolute
image(Top_25_FP_sp.richness_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Species richness for most phylogenetically original species"), 
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
addRasterLegend(Top_25_FP_sp.richness_Jaccard.80, locs = seq(0, 35, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1.1, label = "Species")
legend(legend = "a", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))

# B/ Top 25% - Proportion
image(Top_25_FP_prop_sp.richness_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Proportion of most phylogenetically original species"), 
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
addRasterLegend(Top_25_FP_prop_sp.richness_Jaccard.80_contrasted, locs = seq(0, 50, 10), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.1, label = "Proportion (%)")
legend(legend = "b", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))

# C/ Low 25% - Absolute
image(Low_25_FP_sp.richness_Jaccard.80, col = pal_bl_red_Mannion, main = paste0("Species richness for least phylogenetically original species"), 
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
addRasterLegend(Low_25_FP_sp.richness_Jaccard.80, locs = seq(0, 35, 5), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -111, y = 4, font = 2, cex = 1.1, label = "Species")
legend(legend = "c", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))

# D/ Low 25% - Proportion
image(Low_25_FP_prop_sp.richness_Jaccard.80_contrasted, col = pal_bl_red_Mannion, main = paste0("Proportion of least phylogenetically original species"), 
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
addRasterLegend(Low_25_FP_prop_sp.richness_Jaccard.80_contrasted, locs = seq(0, 80, 20), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -110, y = 4, font = 2, cex = 1.1, label = "Proportion (%)")
legend(legend = "d", x = "bottomright", bty = "n",
       text.font = 2, cex = 1.8, inset=c(0, 0))

par(mfrow = c(1,1))
par(mar = internal_margins)

dev.off()


### 5/ Scatter plot species FP ~ species rarity ####

# Load species indices

sp.rarity_indices_Jaccard.80 <- readRDS(file = "./outputs/Indices_maps/sp.rarity_indices_Jaccard.80.rds")

Phylo_sp.rarity_indices_Jaccard.80 <- sp.rarity_indices_Jaccard.80[which(names(All_sp_proba_stack_Jaccard.80) %in% Phylo.sp_list)]

Species_ES_FP_index_table$Rarity <- Phylo_sp.rarity_indices_Jaccard.80

test <- cbind(names(All_sp_proba_stack_Jaccard.80), sp.range_size_Jaccard.80, sp.rarity_indices_Jaccard.80)

# Highlight species with the 4 highest FP (which happen to also be rather geographically rare)
Risk_species <- as.character(Species_ES_FP_index_table$Sp.full[tail(order(Species_ES_FP_index_table$FP), 4)])
Risk_species_df <- Species_ES_FP_index_table[Species_ES_FP_index_table$Sp.full %in% Risk_species, ]

# Plot

library(ggrepel)

pdf(file = paste0("./graphs/Species_FP_vs_Rarity.pdf"), height = 5.3, width = 6.5)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))

ggplot(Species_ES_FP_index_table, aes(x = FP, y = Rarity)) +
  geom_point() +
  geom_smooth() +
  labs(title = "Species Rarity ~ Phylogenetic originality (FP)") +
  coord_trans(y = "exp") +
  ylim(0, 1) +
  ggrepel::geom_label_repel(# data = Risk_species_df, # If you want to use only a subset of data to plot text
                            # aes(label = Sp.full),
                            
                            data = . %>% 
                              mutate(Sp.full = as.character(Sp.full)) %>% 
                              mutate(label = ifelse(Sp.full %in% Risk_species,
                                                    Sp.full, "")),
                            aes(label = label),
                            
                            alpha = 0.8,     # Transparency of the label and the text
                            fill = "white",  # Background label color
                            color = "black",           # Color is applied to the text and the label border
                            box.padding = 0.63,
                            segment.color = "black",   # Color of the line joining the point and the label. If NA, no segment is drawn.
                            segment.size = 1,    # Width of the line joining the point and the label
                            segment.alpha = 0.5) # Transparency of the line joining the point and the label

par(mar = internal_margins)

dev.off()


?geom_label_repel

# 6/ Map of most phylogenetically original species

pdf(file = paste0("./maps/By_species/Most_phylo_original_sp.pdf"), height = 10, width = 12)

internal_margins <- par()$mar
par(mar = c(3.1,3.5,3.5,2.1))
par(mfrow = c(2,2))

plot(All_sp_proba_stack_Jaccard.80[[Risk_species[1]]], main = paste(Risk_species[1], "\nFP =", round(Species_ES_FP_index_table$FP[Species_ES_FP_index_table$Sp.full == Risk_species[1]], 2), "My")) # Sais.rosalia
plot(All_sp_proba_stack_Jaccard.80[[Risk_species[2]]], main = paste(Risk_species[2], "\nFP =", round(Species_ES_FP_index_table$FP[Species_ES_FP_index_table$Sp.full == Risk_species[2]], 2), "My")) # Eutresis.hypereia
plot(All_sp_proba_stack_Jaccard.80[[Risk_species[3]]], main = paste(Risk_species[3], "\nFP =", round(Species_ES_FP_index_table$FP[Species_ES_FP_index_table$Sp.full == Risk_species[3]], 2), "My")) # Athyrtis.mechanitis
plot(All_sp_proba_stack_Jaccard.80[[Risk_species[4]]], main = paste(Risk_species[4], "\nFP =", round(Species_ES_FP_index_table$FP[Species_ES_FP_index_table$Sp.full == Risk_species[4]], 2), "My")) # Thyridia.psidii

par(mar = internal_margins)
par(mfrow = c(1, 1))

dev.off()