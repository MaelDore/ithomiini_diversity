
###### Script 27: Investigate Mimicry saturation effect ######

# Only for Jaccard.80

# Plot the relationship between species and mimicry richness by bioregions
# Test for siginificance of the relationship by generating randomly assembled communities within each bioregions


###

# Inputs: 
    # Indices maps from script 14a
    # WWF Ecoregions (Olson, 2001)

# Outputs:
   # Shp of bioregions with characteristic relationship (saturation effect)
   # Dataframe with indices values per bioregions
   # Scatter plot mim. richness ~ sp. richness with colored bioregions 
   # Stack of randomized communities and richnesses
   # Scatter plot mim. richness ~ sp. richness with colored bioregions + random values

###




# Clean environment
rm(list = ls())

### Load maps
sp.richness <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
ring.richness <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))

mean.ring.size <- sp.richness/ring.richness

plot(mean.ring.size)
plot(sp.richness)
plot(ring.richness)

# Save
save(mean.ring.size, file = "./outputs/Correlation_tests/mean.ring.size.RData")

# Export raster to use in QGIS
writeRaster(sp.richness, filename = "./outputs/Correlation_tests/sp.richness") # Don't put any extension because the function will write several files (.gri and .grd)
writeRaster(ring.richness, filename = "./outputs/Correlation_tests/ring.richness") # Don't put any extension because the function will write several files (.gri and .grd)
writeRaster(mean.ring.size, filename = "./outputs/Correlation_tests/mean.ring.size") # Don't put any extension because the function will write several files (.gri and .grd)

# Load map stuff
continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_15.rds"))
crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_15.rds"))
bg_mask <- readRDS(file = "./input_data/Map_stuff/bg_mask.rds")
bg_mask_pixel <- readRDS(file = "./input_data/Map_stuff/bg_mask_pixel.rds")

pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")
pal_bl_red_Mannion[1] <- "#EDEDED"


### 1/ Generate bioregion buffers ####

# Try to find which part of the scatter plot belong to which region

plot(sp.richness, zlim = c(25, 35))
plot(ring.richness, zlim = c(5, 7))

# Mata Atlantica
# 25/35 ; 5/7

# Central America = 60129 + 60167 + 60119 + 60130 + 60224
# 30/50 ; 8/10

# Andean slopes South = Peruvian Yungas (60153 ; 53)
# 75/125 ; 18/22

# Ecuador hotpot = Northern Andean Montane Forests (60121 ; 40)
# 50/125 ; 22/30


sp.richness_test <- sp.richness[sp.richness[] > 25 & sp.richness[] < 35]

sp.richness_test <- calc(x = sp.richness, fun = function(x){ x[(x > 25) & (x < 35)] <- NA; return(x)} )

# Get biomes

library(tidyverse)
library(sf)
library(rgeos)
library(alphahull)
library(spatialEco)
library(raster)

# Functions to convert alpha-hull into spatial objects
source("./functions/Alpha_functions.R", local = TRUE) 

Biomes_sf <- st_read("./input_data/WWF_official_teow/wwf_terr_ecos.shp")

xmin = -120 ; xmax = -30 ; ymin = -37 ; ymax = 28
e <- extent(c(xmin,xmax,ymin,ymax))

Neotropics <- Biomes_sf %>% 
  st_crop(e) %>% 
  filter(BIOME != 98)

Central_America <- Neotropics %>% 
  filter(ECO_ID %in% c(60129, 60167, 60119, 60130, 60224))
Central_America_shp <- summarize(Central_America)
Central_America_shp2 <- gSimplify(as(Central_America_shp, 'Spatial'), tol = 0.1, topologyPreserve = F)

plot(Central_America_shp)
plot(Central_America_shp2)

East_Ecuador_hotspot <- Neotropics %>% 
  filter(ECO_ID == 60121)
East_Ecuador_shp <- summarize(East_Ecuador_hotspot)
# East_Ecuador_shp2 <- gSimplify(as(East_Ecuador_shp, 'Spatial'), tol = 0.25, topologyPreserve = F)

East_Ecuador_shp2 <- ashape_to_SPLDF(ashape(unique(st_coordinates(East_Ecuador_shp)[,1:2]), alpha = 0.5)) # Try an alpha-hull first
East_Ecuador_shp2@proj4string <- as(East_Ecuador_shp, 'Spatial')@proj4string
East_Ecuador_shp2 <- gBuffer(spgeom = East_Ecuador_shp2, width = 0.01)
East_Ecuador_shp2 <- remove.holes(East_Ecuador_shp2)

plot(East_Ecuador_shp)
plot(East_Ecuador_shp2)

West_Ecuador_hotspot <- Neotropics %>% 
  filter(ECO_ID %in% c(60178, 60214))
West_Ecuador_shp <- summarize(West_Ecuador_hotspot)

Peruvian_hotspot <- Neotropics %>% 
  filter(ECO_ID == 60153)
Peruvian_shp <- summarize(Peruvian_hotspot)
# Peruvian_shp2 <- gSimplify(as(Peruvian_shp, 'Spatial'), tol = 0.3, topologyPreserve = F)

Peruvian_shp2 <- ashape_to_SPLDF(ashape(unique(st_coordinates(Peruvian_shp)[,1:2]), alpha = 0.5)) # Try an alpha-hull first
Peruvian_shp2@proj4string <- as(Peruvian_shp, 'Spatial')@proj4string
Peruvian_shp2 <- gBuffer(spgeom = Peruvian_shp2, width = 0.01)
Peruvian_shp2 <- remove.holes(Peruvian_shp2)

plot(Peruvian_shp)
plot(Peruvian_shp2)


Bolivian_hotspot <- Neotropics %>% 
  filter(ECO_ID == 60105)
Bolivian_shp <- summarize(Bolivian_hotspot)


# Merge by biome to extract all ecoregions for Mata Atlantica more easily
Neotropics_biomes <- Neotropics %>% 
  group_by(BIOME) %>% 
  summarize()  # Pour union les features d'après le grouping

plot(Neotropics_biomes[1,])

Forests <- Neotropics_biomes[1,] %>% 
  st_cast("POLYGON") %>% 
  mutate(area = st_area(.)) %>% 
  arrange(desc(area))

plot(Forests$geometry, col = c(1:30, "grey"))
plot(Forests[1,])
plot(Forests[2,])
plot(Forests[5,]) # Pernambuco section
plot(Forests[8,]) # Brazilia extent

plot(1:20, 1:20, col = c(1:20, "grey"))

Mata_Atlantica <- Forests[2,]
# Mata_Atlantica <- Forests[c(2,5,8),]
Mata_Atlantica_shp <- as(Mata_Atlantica, 'Spatial')
Mata_Atlantica_shp2 <- gSimplify(Mata_Atlantica_shp, tol = 0.3, topologyPreserve = F)

Mata_Atlantica_shp2 <- ashape_to_SPLDF(ashape(unique(st_coordinates(Mata_Atlantica)[,1:2]), alpha = 0.5)) # Try an alpha-hull first
Mata_Atlantica_shp2@proj4string <- as(Mata_Atlantica_shp, 'Spatial')@proj4string
Mata_Atlantica_shp2 <- gBuffer(spgeom = Mata_Atlantica_shp2, width = 0.01)
Mata_Atlantica_shp2 <- remove.holes(Mata_Atlantica_shp2)
Mata_Atlantica_shp2 <- remove.holes(Mata_Atlantica_shp2)

plot(Mata_Atlantica_shp)
plot(Mata_Atlantica_shp2)

crs(Mata_Atlantica_shp) <- crs(Mata_Atlantica_shp2) <- crs(Central_America_shp2)

# Option with the Terrabrasilis biome

Mata_Atlantica_sf <- st_read("./input_data/Map_stuff/Terrabrasilis_biome_limits/mata_atlantica_border/mata_atlantica_border.shp")


plot(Mata_Atlantica_sf)
Mata_Atlantica_shp3 <- as(Mata_Atlantica_sf, 'Spatial')
Mata_Atlantica_shp4 <- gSimplify(Mata_Atlantica_shp3, tol = 0.36, topologyPreserve = F)
plot(Mata_Atlantica_shp4)

# Merge both to add the Northern East portion, but keep the SouthWest one
Mata_Atlantica_shp5 <- gUnion(Mata_Atlantica_shp2, Mata_Atlantica_shp4)
plot(Mata_Atlantica_shp5, col = "limegreen", add = T)

### Save bioregions simplified shp ###
saveRDS(Central_America_shp, file = "./input_data/Map_stuff/Bioregions/Central_America_shp.rds", version = "2")
saveRDS(Central_America_shp2, file = "./input_data/Map_stuff/Bioregions/Central_America_shp2.rds", version = "2")
saveRDS(East_Ecuador_shp, file = "./input_data/Map_stuff/Bioregions/East_Ecuador_shp.rds", version = "2")
saveRDS(East_Ecuador_shp2, file = "./input_data/Map_stuff/Bioregions/East_Ecuador_shp2.rds", version = "2")
saveRDS(West_Ecuador_shp, file = "./input_data/Map_stuff/Bioregions/West_Ecuador_shp.rds", version = "2")
saveRDS(Peruvian_shp, file = "./input_data/Map_stuff/Bioregions/Peruvian_shp.rds", version = "2")
saveRDS(Peruvian_shp2, file = "./input_data/Map_stuff/Bioregions/Peruvian_shp2.rds", version = "2")
saveRDS(Bolivian_shp, file = "./input_data/Map_stuff/Bioregions/Bolivian_shp.rds", version = "2")
saveRDS(Mata_Atlantica_shp, file = "./input_data/Map_stuff/Bioregions/Mata_Atlantica_shp.rds", version = "2")
saveRDS(Mata_Atlantica_shp2, file = "./input_data/Map_stuff/Bioregions/Mata_Atlantica_shp2.rds", version = "2")
saveRDS(Mata_Atlantica_shp3, file = "./input_data/Map_stuff/Bioregions/Mata_Atlantica_shp3.rds", version = "2")
saveRDS(Mata_Atlantica_shp4, file = "./input_data/Map_stuff/Bioregions/Mata_Atlantica_shp4.rds", version = "2")
saveRDS(Mata_Atlantica_shp5, file = "./input_data/Map_stuff/Bioregions/Mata_Atlantica_shp5.rds", version = "2")

### Load bioregions simplified shp ###
Central_America_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Central_America_shp.rds")
Central_America_shp2 <- readRDS(file = "./input_data/Map_stuff/Bioregions/Central_America_shp2.rds")
East_Ecuador_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/East_Ecuador_shp.rds")
East_Ecuador_shp2 <- readRDS(file = "./input_data/Map_stuff/Bioregions/East_Ecuador_shp2.rds")
West_Ecuador_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/West_Ecuador_shp.rds")
Peruvian_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Peruvian_shp.rds")
Peruvian_shp2 <- readRDS(file = "./input_data/Map_stuff/Bioregions/Peruvian_shp2.rds")
Bolivian_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Bolivian_shp.rds")
Mata_Atlantica_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Mata_Atlantica_shp.rds")
Mata_Atlantica_shp5 <- readRDS(file = "./input_data/Map_stuff/Bioregions/Mata_Atlantica_shp5.rds")

### Plot the regions ###
pdf(file = paste0("./maps/Indices_maps/mean.ring.size.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))

mean.ring.size@data@min <- 0

image(mean.ring.size, col = pal_bl_red_Mannion)
plot(crop_mask_shp, lwd = 1, border = "grey80", add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)

plot(Central_America_shp2, col = NA, border = "red", lwd = 1.5, add = T)
plot(East_Ecuador_shp2, col = NA, border = "darkgreen", lwd = 1.5, add = T)
# plot(West_Ecuador_shp, col = NA, border = "cyan", add = T)
plot(Peruvian_shp2, col = NA, border = "blue", lwd = 1.5, add = T)
# plot(Bolivian_shp, col = NA, border = "darkgreen", add = T)
plot(Mata_Atlantica_shp5, col = NA, border = "purple", lwd = 1.5, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-100, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
rangeBuilder::addRasterLegend(mean.ring.size, locs = seq(0, 6, 2), cex.axis = 1.2, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 6.5, font = 2, cex = 1.2, label = "Mean\nring size")

legend(x = "bottomleft", legend = c("Central America", "East Ecuador", "East Peru", "Mata Atlantica"), 
       border = c("red", "darkgreen", "blue", "purple") , fill = NA,
       inset = c(0.15,0.15),
       bty = "n") 

par(mar = internal_margins)
dev.off()

### 2/ Generate df for plot ####

# Load stack of indices maps
all_indices_stack <- readRDS(file = "./outputs/Correlation_tests/all_indices_stack.rds")

# Get indices of pixels within regions
Central_America_mask <- mask(all_indices_stack, Central_America_shp)
Central_America_index <- which(!is.na(Central_America_mask@data@values[,1]))

East_Ecuador_mask <- mask(all_indices_stack, East_Ecuador_shp)
East_Ecuador_index <- which(!is.na(East_Ecuador_mask@data@values[,1]))

# West_Ecuador_mask <- mask(all_indices_stack, West_Ecuador_hotspot)
# West_Ecuador_index <- which(!is.na(West_Ecuador_mask@data@values[,1]))

Peruvian_mask <- mask(all_indices_stack, Peruvian_shp)
Peruvian_index <- which(!is.na(Peruvian_mask@data@values[,1]))

# Bolivian_mask <- mask(all_indices_stack, Bolivian_hotspot)
# Bolivian_index <- which(!is.na(Bolivian_mask@data@values[,1]))

Mata_Atlantica_mask <- mask(all_indices_stack, Mata_Atlantica_shp5)
Mata_Atlantica_index <- which(!is.na(Mata_Atlantica_mask@data@values[,1]))

# Save indices
save(Central_America_index, East_Ecuador_index, Peruvian_index, Mata_Atlantica_index, file = "./outputs/Correlation_tests/Bioregions_indices.RData")

# Load df with all values
# indices_df <- readRDS(file = "./outputs/Correlation_tests/Indices_df.rds")
indices_df <- tibble(sp.richness = sp.richness[], sp.mean.rarity = sp.mean.rarity[], MPD = MPD[], ring.richness = ring.richness[], mimicry.mean.rarity = mimicry.mean.rarity[])

# Add column with region as factor
indices_df$Region <- "Other"
indices_df$Region[Central_America_index] <- "Central_America"
indices_df$Region[East_Ecuador_index] <- "East_Ecuador"
# indices_df$Region[West_Ecuador_index] <- "West_Ecuador"
indices_df$Region[Peruvian_index] <- "East_Peru"
# indices_df$Region[Bolivian_index] <- "Bolivian"
indices_df$Region[Mata_Atlantica_index] <- "Mata_Atlantica"

indices_df$Region <- as.factor(indices_df$Region)

indices_df <- indices_df %>% 
  mutate(Region = fct_relevel(Region, 
                          "Other", 
                          # "Bolivian", 
                          "Central_America", 
                          "East_Ecuador", 
                          # "West_Ecuador", 
                          "East_Peru", 
                          "Mata_Atlantica"))

# Remove NA data
indices_df <- indices_df %>% filter(!is.na(indices_df$sp.richness))

# Save indices file
saveRDS(indices_df, file = "./outputs/Correlation_tests/Indices_df_bioregions.rds", version = "2")

# Remove other regions
# indices_df <- indices_df %>% drop_na()



### 3/ Plot scatterplot with different colors for each bioregions ####

indices_df <- readRDS(file = "./outputs/Correlation_tests/Indices_df_bioregions.rds")
limits_df <- tibble(x = c(0, 150), y = c(0, 30))
legend_df <- data.frame(fake_x = as.numeric(c(NA, NA, NA, NA)),
                        Region = c("CA", "EE", "EP", "BAF"))

library(tidyverse)
library(ggthemes)
library(scales)

pdf(file = "./graphs/Correlation_tests/Bioregions_Rich_Mim.pdf", height = 6, width = 8)
g1 <- ggplot(data = indices_df) +
  geom_point(data = indices_df[indices_df$Region == "Other",], aes(x = sp.richness, y = ring.richness, fill = "black", alpha = 0.1), show.legend = F) +
  geom_point(data = indices_df[indices_df$Region != "Other",], aes(x = sp.richness, y = ring.richness, col = Region), show.legend = F) +
  geom_smooth(aes(x = sp.richness, y = ring.richness), color = "#FF000080", size = 1.5, se = F) +
  geom_bar(data = legend_df, aes(x = fake_x, fill = Region), show.legend = T) +
  # annotate("text", x = max(indices_df$sp.richness), y = min(indices_df$ring.richness), hjust = 1, vjust = 0, label = "C", size = 7, fontface = 2) +
  
  # annotate("text", x = min(indices_df$sp.richness), y = max(indices_df$ring.richness), hjust = 0.010, vjust = 8.0, label = bquote(rho ~ "= 0.93"), size = 5, fontface = 2) +
  
  annotate("text", x = min(indices_df$sp.richness), y = max(indices_df$ring.richness), hjust = 0.05, vjust = 10.1, label = bquote(rho ~ "= 0.93"), size = 5, fontface = 2) +
  annotate("text", x = min(indices_df$sp.richness), y = max(indices_df$ring.richness), hjust = -0.360, vjust = 11.7, label = "GAM predict", size = 4.1, fontface = 2) +
  
  geom_segment(x = -1.1, y = 20.2, xend = 6, yend = 20.2, color = "#FF000080", alpha = 0.7, size = 1.5) +
  
  ylab(label = "Mimicry richness\n[Mimicry rings]") +
  xlab(label = "Species richness\n[Species]") +
  
  guides(alpha = F, color = F) +
  scale_x_continuous(breaks = seq(0, 150, 25)) +
  # scale_color_discrete(name = "Bioregions", labels = c("CA", "EE", "EP", "BAF"), values = c("Central_America" = "#F8766D", "East_Ecuador" = "#7CAE00", "East_Peru" = "#00BFC4", "Mata_Atlantica" = "#C77CFF")) +
  ggthemes::geom_rangeframe(data = limits_df, aes(x = x, y = y), size = 1.4) +
  # scale_fill_discrete(name = "Bioregions", labels = c("CA", "EE", "EP", "BAF"), limits = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
  scale_fill_discrete(name = "Bioregions", labels = c("Central America", "East Ecuador", "East Peru", "Atlantic Forest"), limits = c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")) +
  
  theme(panel.background = element_rect(fill = "white"),
        legend.position = c(0.139, 0.83),
        legend.background = element_rect(fill = NA, colour = NA),
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.title = element_text(size = 12, vjust = 3, face = "bold"),
        legend.text = element_text(size = 10, face = "bold"),
                      # margin = margin(t = 0, unit = "pt")),
        # legend.key.size = unit(1.5, 'lines'),
        # legend.spacing.y = unit(1,"cm"),
        axis.ticks = element_line(size = 1.2),
        axis.ticks.length = unit(5, "pt"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title.y = element_text(margin = margin(l = 5, r = 10)))
  

print(g1)
dev.off()

?scale_fill_discrete
c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
?element_text


### With inset_map

library(ggspatial)
library(cowplot)
library(sf)

crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_15.rds"))

load(file = "./input_data/Map_stuff/continent_mask_smooth.RData") # Load bg shp
load(file = "./input_data/Map_stuff/bbox_sp.RData")
bg_mask <- readRDS(file = "./input_data/Map_stuff/bg_mask.rds")

scale_bar <- st_crop(x = st_as_sf(bbox_sp), xmin = -115, xmax = -96, ymin = -33, ymax = -32.2)  

inset_map <- ggplot() + 
  theme_void() +
  geom_sf(data = st_as_sf(bbox_sp), fill = NA) +
  geom_sf(data = st_as_sf(continent_mask_smooth), fill = "#EDEDED") +
  geom_sf(data = st_as_sf(Central_America_shp2), fill = "#F8766D", col = "#F8766D") +  # "#F8766D"  "red"
  geom_sf(data = st_as_sf(East_Ecuador_shp2), fill = "#7CAE00", col = "#7CAE00") +  # "#7CAE00"  "darkgreen"
  geom_sf(data = st_as_sf(Peruvian_shp2), fill = "#00BFC4", col = "#00BFC4") +  # "#00BFC4"  "blue"
  geom_sf(data = st_as_sf(Mata_Atlantica_shp5), fill = "#C77CFF", col = "#C77CFF") +  # "#C77CFF"  "purple"
  geom_sf(data = st_as_sf(bg_mask), fill = "aliceblue", col = "grey20") +
  geom_sf(data = scale_bar, fill = "black", col = "black") # +
  # annotate("text", x = -105, y = -30, label = "2000 km", fontface = "bold", size = 4.5) +
  # ggspatial::annotation_north_arrow(height = unit(1.3, "cm"), width = unit(1.3, "cm"), pad_x = unit(1.4, "cm"), pad_y = unit(1.2, "cm"), location = "topright")

print(inset_map)

# Graph with inset map

final_inset_graph <- ggdraw() +
  draw_plot(g1) +
  draw_plot(inset_map, x = 0.54, y = 0.16, width = 0.43, height = 0.43) +
  annotate("text", x = 0.623, y = 0.235, label = "2000 km", fontface = "bold", size = 3) +
  ggspatial::annotation_north_arrow(height = unit(0.7, "cm"), width = unit(0.7, "cm"), pad_x = unit(1.2, "cm"), pad_y = unit(7.5, "cm"), location = "br", style = north_arrow_orienteering(text_col = "aliceblue"))

# pdf(file = "./graphs/Correlation_tests/Bioregions_Rich_Mim_with_map.pdf", height = 6, width = 8)
# print(final_inset_graph)
# dev.off()

pdf(file = "./graphs/Correlation_tests/Bioregions_Rich_Mim_with_map_and_GAM.pdf", height = 6, width = 8)
print(final_inset_graph)
dev.off()
   

# Version plot() à intégrer dans l'autre plot via Illustrator
plot(bbox_sp)
plot(continent_mask_smooth, col = "#EDEDED", add = T)
plot(Central_America_shp2, col = "red", border = "red", lwd = 1.5, add = T)
plot(East_Ecuador_shp2, col = "darkgreen", border = "darkgreen", lwd = 1.5, add = T)
plot(Peruvian_shp2, col = "blue", border = "blue", lwd = 1.5, add = T)
plot(Mata_Atlantica_shp5, col = "purple", border = "purple", lwd = 1.5, add = T)
plot(bg_mask, lwd = 0.8, border = "grey20", col = "aliceblue", add = T)
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-115, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.5, padin = c(0.4, 0.5))


### 4/ Randomization "test" to see if we retrieve the same pattern when regional communities assemble at random ####

# Load Summary table for OMU/unit
load(file = paste0("./input_data/list.models.RData"))

### 4.1/ Loop by OMU to stack all OMU proba maps
continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_15.rds"))
OMU.stack_Jaccard.80 <- stack(continent_mask) # Initiate stack
for (i in 1:nrow(list.models)) 
{
  # i <- 1
  
  unit <-  as.character(list.models$Tag.model[i]) # Load the unit name 
  
  if (list.models$initial_model_type[i] != "rasterized")
  {
    OMU_map <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_Jaccard_median_cropped_80.rds"))
  } 
  else 
  {
    OMU_map <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/rasterized_map_",unit,".rds"))
  }
  
  OMU.stack_Jaccard.80 <- addLayer(OMU.stack_Jaccard.80, OMU_map)
  
  if (i %% 10 == 0) {print(i)}
  
}

OMU.stack_Jaccard.80 <- dropLayer(OMU.stack_Jaccard.80, i = 1)
names(OMU.stack_Jaccard.80) <- list.models$Tag.model

saveRDS(OMU.stack_Jaccard.80, file = "./outputs/By_OMU/OMU.stack_Jaccard.80.rds", version = "2")

OMU.stack_Jaccard.80 <- readRDS(file = "./outputs/By_OMU/OMU.stack_Jaccard.80.rds")

### Choose the bioregion ###
bioreg <- "Central_America" ; cropping_shp <- Central_America_shp
bioreg <- "East_Ecuador" ; cropping_shp <- East_Ecuador_shp
bioreg <- "Peruvian" ; cropping_shp <- Peruvian_shp
bioreg <- "Mata_Atlantica" ; cropping_shp <- Mata_Atlantica_shp5


### Crop for specific bioregion ###
bioreg_OMU.stack_Jaccard.80 <- crop(OMU.stack_Jaccard.80, cropping_shp) # Crop to reduce raster size to minimum
bioreg_OMU.stack_Jaccard.80 <- mask(bioreg_OMU.stack_Jaccard.80, cropping_shp) # Mask communities outside of bioregion
presence_test <- cellStats(x = bioreg_OMU.stack_Jaccard.80, stat = sum, na.rm = T) # Check presence of OMU in the bioregion
bioreg_OMU.stack_Jaccard.80 <- subset(bioreg_OMU.stack_Jaccard.80, subset = which(presence_test != 0)) # Keep only layers for OMU present in the bioregion
nlayers(bioreg_OMU.stack_Jaccard.80) # Number of OMU
bioreg_Present_OMU <- names(bioreg_OMU.stack_Jaccard.80) # Extract OMU names

plot(bioreg_OMU.stack_Jaccard.80)

# Randomize OMU distri patterns (keep values but change them in space) ###
set.seed(42)

bioreg_OMU.stack_Jaccard.80_random <- bioreg_OMU.stack_Jaccard.80
for (i in 1:nlayers(bioreg_OMU.stack_Jaccard.80)) # Loop per OMU
{
  # i <- 1
  
  # plot(bioreg_OMU.stack_Jaccard.80[[i]])
  
  proba_data <- bioreg_OMU.stack_Jaccard.80[[i]][]
  noNA_commu <- !is.na(proba_data) # Get position of communities without NA = the ones you want to shuffle
  noNA_commu_data <- proba_data[noNA_commu] # Extract only cells with data
  noNA_commu_data_random <- sample(noNA_commu_data) # Shuffle data
  
  proba_data_random <- proba_data
  proba_data_random[noNA_commu] <- noNA_commu_data_random # Replace in original data vector
  bioreg_OMU.stack_Jaccard.80_random[[i]][] <- proba_data_random # Replace in raster layer
  
  # plot(bioreg_OMU.stack_Jaccard.80_random[[i]])
  
  print(i)
  
}

# Compare OMU stacks
plot(bioreg_OMU.stack_Jaccard.80) # Original
plot(bioreg_OMU.stack_Jaccard.80_random) # Randomized

# Function to compute probability of presence with multiple OMUs
aggreg_prob = function(x, na.rm) { 
  y <- 1-prod(1-x) # Probability of presence of species = probability of presence of at least one OMU = opposite of probability of absence of all OMU
  return(y) # Output
}

### 4.2/ Aggregate at species level ####

load(file = paste0("./input_data/list.sp.RData"))

bioreg_list.sp <- list.models$Sp_full[list.models$Tag.model %in% bioreg_Present_OMU] # Get species for each OMU
bioreg_list.sp_unique <- unique(bioreg_list.sp) # Get unique species list
bioreg_sp.stack_Jaccard.80_random <- bioreg_OMU.stack_Jaccard.80[[1]] # Initiate stack

for (i in 1:length(bioreg_list.sp_unique))# Loop for each unique species
{
  # i <- 1
  
  sp <- bioreg_list.sp_unique[i]
  
  OMU_list <- bioreg_Present_OMU[bioreg_list.sp ==  sp] # Get OMU for this species
  
  # Extract only layers of OMU for this species
  sp_bioreg.OMU.stack_Jaccard.80 <- subset(bioreg_OMU.stack_Jaccard.80_random, which(names(bioreg_OMU.stack_Jaccard.80_random) %in% OMU_list))
  
  # Compute probability of presence of species (at least one OMU)
  sp.cont_Jaccard.80 <- calc(sp_bioreg.OMU.stack_Jaccard.80, fun = aggreg_prob)
  
  # Add to final stack
  bioreg_sp.stack_Jaccard.80_random <- stack(bioreg_sp.stack_Jaccard.80_random, sp.cont_Jaccard.80)
  
  print(i)
}
bioreg_sp.stack_Jaccard.80_random <- dropLayer(bioreg_sp.stack_Jaccard.80_random, i = 1)
names(bioreg_sp.stack_Jaccard.80_random) <- bioreg_list.sp_unique
nlayers(bioreg_sp.stack_Jaccard.80_random) # Nb of species in the bioregion

bioreg_sp.richness_Jaccard.80_random <- calc(bioreg_sp.stack_Jaccard.80_random, fun = sum)

plot(bioreg_sp.richness_Jaccard.80_random)

sp.richness <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))

# Baseline without randomisation
test <- crop(sp.richness, cropping_shp)
test <- mask(test, cropping_shp)
plot(test)


### 4.3/ Aggregate at mimicry ring level ####

load(file = paste0("./input_data/list.models.RData"))
mimicry.list <- as.character(unique(list.models$Mimicry.model)) # 44 Mimicry rings

bioreg_list.ring <- as.character(list.models$Mimicry.model[list.models$Tag.model %in% bioreg_Present_OMU]) # Get mimicry ring for each OMU
bioreg_list.ring_unique <- unique(bioreg_list.ring) # Get unique mimicry ring list
bioreg_ring.stack_Jaccard.80_random <- bioreg_OMU.stack_Jaccard.80[[1]] # Initiate stack

for (i in 1:length(bioreg_list.ring_unique))# Loop for each unique mimicry ring
{
  # i <- 1
  
  ring <- bioreg_list.ring_unique[i]
  
  OMU_list <- bioreg_Present_OMU[bioreg_list.ring ==  ring] # Get OMU for this species
  
  # Extract only layers of OMU for this ring
  ring_bioreg.OMU.stack_Jaccard.80 <- subset(bioreg_OMU.stack_Jaccard.80_random, which(names(bioreg_OMU.stack_Jaccard.80_random) %in% OMU_list))
  
  # Compute probability of presence of mimicry ring (at least one OMU)
  ring.cont_Jaccard.80 <- calc(ring_bioreg.OMU.stack_Jaccard.80, fun = aggreg_prob)
  
  # Add to final stack
  bioreg_ring.stack_Jaccard.80_random <- stack(bioreg_ring.stack_Jaccard.80_random, ring.cont_Jaccard.80)
  
  print(i)
}
bioreg_ring.stack_Jaccard.80_random <- dropLayer(bioreg_ring.stack_Jaccard.80_random, i = 1)
names(bioreg_ring.stack_Jaccard.80_random) <- bioreg_list.ring_unique
nlayers(bioreg_ring.stack_Jaccard.80_random) # Nb of mimicry rings

bioreg_ring.richness_Jaccard.80_random <- calc(bioreg_ring.stack_Jaccard.80_random, fun = sum)

plot(bioreg_ring.richness_Jaccard.80_random)

ring.richness <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))

# Baseline without randomisation
test <- crop(ring.richness, cropping_shp)
test <- mask(test, cropping_shp)
plot(test)

# Save random index maps
saveRDS(bioreg_ring.richness_Jaccard.80_random, file = paste0("./outputs/Correlation_tests/",bioreg,"_ring.richness_Jaccard.80_random.rds"))
saveRDS(bioreg_sp.richness_Jaccard.80_random, file = paste0("./outputs/Correlation_tests/",bioreg,"_sp.richness_Jaccard.80_random.rds"))

### 5/ Plots with random values ####

### Load the random maps for species and mimicry richness

# Clean environment
rm(list = ls())

CA_sp.richness_Jaccard.80_random <- readRDS(file = paste0("./outputs/Correlation_tests/Central_America_sp.richness_Jaccard.80_random.rds"))
CA_ring.richness_Jaccard.80_random <- readRDS(file = paste0("./outputs/Correlation_tests/Central_America_ring.richness_Jaccard.80_random.rds"))
CA_random <- stack(CA_sp.richness_Jaccard.80_random, CA_ring.richness_Jaccard.80_random)

Ecuador_sp.richness_Jaccard.80_random <- readRDS(file = paste0("./outputs/Correlation_tests/East_Ecuador_sp.richness_Jaccard.80_random.rds"))
Ecuador_ring.richness_Jaccard.80_random <- readRDS(file = paste0("./outputs/Correlation_tests/East_Ecuador_ring.richness_Jaccard.80_random.rds"))
Ecuador_random <- stack(Ecuador_sp.richness_Jaccard.80_random, Ecuador_ring.richness_Jaccard.80_random)

Peruvian_sp.richness_Jaccard.80_random <- readRDS(file = paste0("./outputs/Correlation_tests/Peruvian_sp.richness_Jaccard.80_random.rds"))
Peruvian_ring.richness_Jaccard.80_random <- readRDS(file = paste0("./outputs/Correlation_tests/Peruvian_ring.richness_Jaccard.80_random.rds"))
Peruvian_random <- stack(Peruvian_sp.richness_Jaccard.80_random, Peruvian_ring.richness_Jaccard.80_random)

MA_sp.richness_Jaccard.80_random <- readRDS(file = paste0("./outputs/Correlation_tests/Mata_Atlantica_sp.richness_Jaccard.80_random.rds"))
MA_ring.richness_Jaccard.80_random <- readRDS(file = paste0("./outputs/Correlation_tests/Mata_Atlantica_ring.richness_Jaccard.80_random.rds"))
MA_random <- stack(MA_sp.richness_Jaccard.80_random, MA_ring.richness_Jaccard.80_random)


All_random_maps <- list(CA_sp.richness_Jaccard.80_random, CA_ring.richness_Jaccard.80_random, Ecuador_sp.richness_Jaccard.80_random, Ecuador_ring.richness_Jaccard.80_random,
                         Peruvian_sp.richness_Jaccard.80_random, Peruvian_ring.richness_Jaccard.80_random, MA_sp.richness_Jaccard.80_random, MA_ring.richness_Jaccard.80_random)

All_random_stacks <- list(CA_random, Ecuador_random, Peruvian_random, MA_random)
bioreg_list <- c("Central_America", "East_Ecuador", "East_Peru", "Mata_Atlantica")

# Extract data in tidy df
indices_df_random <- NULL
for(i in 1:length(All_random_stacks)) 
{
  # i <- 1
  indices_df_random <- tibble(sp.richness = All_random_stacks[[i]][][, 1], ring.richness = All_random_stacks[[i]][][, 2], Region = bioreg_list[i]) %>% 
    drop_na() %>% 
    bind_rows(indices_df_random, .)
}


### 5.1/ Plot scatterplot with different colors for each bioregions, only random ####
pdf(file = "./graphs/Correlation_tests/Bioregions_Rich_Mim_random.pdf", height = 6, width = 8)
g1 <- ggplot(data = indices_df_random) +
  geom_point(aes(x = sp.richness, y = ring.richness, col = Region)) +
  ylab(label = "Mimicry richness\n[Mimicry rings]") +
  xlab(label = "Species richness\n[Species]") +
  # scale_x_continuous(breaks = seq(0, 125, 25)) +
  # scale_colour_manual(values = c("Central_America" = "red", "East_Ecuador" = "royalblue", "East_Peru" = "purple", "Mata_Atlantica" = "darkgreen")) +
  theme_gray() +
  theme(panel.grid.minor.x = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.ticks = element_line(size = 1.2),
        axis.ticks.length = unit(5, "pt"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title.y = element_text(margin = margin(l = 5, r = 10)))

print(g1)
dev.off()



### 5.2/ Plot scatterplot with different colors for each bioregions, random + original ####

indices_df <- readRDS(file = "./outputs/Correlation_tests/Indices_df_bioregions.rds")
indices_df_random_orig <- bind_rows(indices_df, indices_df_random) %>% 
  mutate(type = c(rep("original", nrow(indices_df)), rep("random", nrow(indices_df_random))))

pdf(file = "./graphs/Correlation_tests/Bioregions_Rich_Mim_random_vs_orig.pdf", height = 6, width = 8)
g2 <- ggplot(data = indices_df_random_orig) +
  geom_point(data = indices_df_random_orig[indices_df_random_orig$Region == "Other",], aes(x = sp.richness, y = ring.richness, fill = "black", alpha = 0.1, col = NA), show.legend = F) +
  geom_point(data = indices_df_random_orig[indices_df_random_orig$Region != "Other" & indices_df_random_orig$type == "original",], aes(x = sp.richness, y = ring.richness, col = Region)) +
  geom_point(data = indices_df_random_orig[indices_df_random_orig$type == "random",], aes(x = sp.richness, y = ring.richness, fill = Region), col = "black", pch = 21, show.legend = F) +
  # annotate("text", x = max(indices_df$sp.richness), y = min(indices_df$ring.richness), hjust = 1, vjust = 0, label = "C", size = 7, fontface = 2) +
  annotate("text", x = min(indices_df$sp.richness), y = max(indices_df$ring.richness), hjust = 0, vjust = 0.8, label = bquote(rho ~ "= 0.93"), size = 5, fontface = 2) +
  ylab(label = "Mimicry richness\n[Mimicry rings]") +
  xlab(label = "Species richness\n[Species]") +
  scale_x_continuous(breaks = seq(0, 125, 25)) +
  # scale_colour_manual(values = c("Central_America" = "red", "East_Ecuador" = "royalblue", "East_Peru" = "purple", "Mata_Atlantica" = "darkgreen")) +
  scale_fill_manual(values = c("black" = NA, "Central_America" = "#F8766D", "East_Ecuador" = "#7CAE00", "East_Peru" = "#00BFC4", "Mata_Atlantica" = "#C77CFF")) +
  theme_gray() +
  theme(panel.grid.minor.x = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.ticks = element_line(size = 1.2),
        axis.ticks.length = unit(5, "pt"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title.y = element_text(margin = margin(l = 5, r = 10))) # +
  # geom_point(data = Bioreg_mean_rich, aes(x = Bioreg_mean_sp.rich, y = Bioreg_mean_ring.rich), size = 2, col = "black", fill = "red", pch = 21)

print(g2)
dev.off()




gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

gg_color_hue(n = 4)


Bioreg_mean_sp.rich <- by(indices_df$sp.richness, indices_df$Region, mean)
Bioreg_mean_ring.rich <- by(indices_df$ring.richness, indices_df$Region, mean)
Bioreg_mean_rich <- as.data.frame(cbind(Bioreg_mean_sp.rich, Bioreg_mean_ring.rich))[-1,]

