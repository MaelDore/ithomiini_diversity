
###### Script 28: Test correlation between higher mimicry richness than expected and terrain ruggedness ######

# Test correlation between higher mimicry richness than expected and terrain ruggedness

###

# Inputs: 
    # DEM (SRTM Dataset v.4.1; Farr et al., 2007) 
    # Species and mimicry richness maps from script 14a

# Outputs:
    # Map of terrain ruggedness index (TRI)
    # Boxplot of TRI per bioregions
    # Boxplot of GAM residuals per bioregions

###


# Effacer l'environnement
rm(list = ls())
gc()


library(raster)

# Change temp folder for raster files
rasterOptions(tmpdir = paste0("./temp"))
# To clean regularly temp folder
unlink(list.files(path = rasterOptions()$tmpdir, all.files = T, recursive = T, include.dirs = T, full.names = T), force = T, recursive = T) # Clean raster store in temp

# Default color palette in ggplot
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

# Load indices map
sp.richness <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
ring.richness <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))


#### 1/ Generate DEM with proper resolution ####

# Load elevation raster
DEM <- raster(x = paste0("./input_data/SRTM/SRTM_5min.tif"))

# SRTM 90m DEM Version 4. Load tiles
DEM_1 <- raster(x = paste0("./input_data/SRTM/cut_n00w060.tif"))
DEM_2 <- raster(x = paste0("./input_data/SRTM/cut_n00w090.tif"))
DEM_3 <- raster(x = paste0("./input_data/SRTM/cut_n00w120.tif"))
DEM_4 <- raster(x = paste0("./input_data/SRTM/cut_s30w060.tif"))
DEM_5 <- raster(x = paste0("./input_data/SRTM/cut_s30w090.tif"))
DEM_6 <- raster(x = paste0("./input_data/SRTM/cut_s30w120.tif"))
DEM_7 <- raster(x = paste0("./input_data/SRTM/cut_s60w060.tif"))
DEM_8 <- raster(x = paste0("./input_data/SRTM/cut_s60w090.tif"))

#  Assemblage of tiles
# DEM_SRTM_v4_90m <- merge(DEM_1, DEM_2)
DEM_SRTM_v4_90m <- DEM_1
for (i in 2:8) {
  DEM_to_add <- eval(parse(text = paste0("DEM_",i)))
  DEM_SRTM_v4_90m <- merge(DEM_SRTM_v4_90m, DEM_to_add)
}

DEM_SRTM_v4_90m@crs@projargs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Save full 90m resolution map
plot(DEM_SRTM_v4_90m)
saveRDS(DEM_SRTM_v4_90m, file = "./input_data/SRTM/DEM_SRTM_v4_90m.rds", version = "2")

head(DEM_SRTM_v4_90m[])


# Version with DEM at 1' resolution (1.8km x 1.8km)
# for (i in 1:8) {
for (i in 2:8) {
  DEM_to_aggregate <- eval(parse(text = paste0("DEM_",i)))
  
  eval(call("<-", as.name(paste0("DEM_",i,"_1m")), raster::aggregate(DEM_to_aggregate, fact = 20, fun = mean, expand = T, na.rm = T)))
  
  saveRDS(readAll(eval(parse(text = paste0("DEM_",i,"_1m")))), file = paste0("./input_data/SRTM/DEM_",i,"_1m.rds"), version = "2")
  gc() # Clean temporary files from virtual memory
  
}

plot(DEM_8_1m)

#  Assemblage of tiles for 1m resolution
# DEM_SRTM_v4_1m <- merge(DEM_1_1m, DEM_2_1m)
DEM_SRTM_v4_1m <- DEM_1_1m
for (i in 2:8) {
  DEM_to_add <- eval(parse(text = paste0("DEM_",i,"_1m")))
  DEM_SRTM_v4_1m <- merge(DEM_SRTM_v4_1m, DEM_to_add)
}

DEM_SRTM_v4_1m@crs@projargs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Save full 1m resolution map
plot(DEM_SRTM_v4_1m)
saveRDS(DEM_SRTM_v4_1m, file = "./input_data/SRTM/DEM_SRTM_v4_1m.rds", version = "2")


# #### 2/ Masking pixels outside of terrestrial lands
# 
# # Mask for 15min arc and final extent
# continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_15.rds"))
# 
# # Get a mask for terrestrial lands, if needed to discard water bodies.
# 
# library("rnaturalearth")
# 
# world_map <- ne_countries(scale = "medium",
#                           continent = NULL,     # To extract only a specific continent
#                           country = NULL,       # To extract a specific country
#                           returnclass = "sf")   # To get an sp object (with geom in list, an attributes in separated table) or an sf object (with geom nested in a sfc) as output
# world_map <- world_map$geometry # Extract only the geometry
# world_map_sp <- as(world_map, Class = "Spatial") # Convert worldmap to sp object
# 
# saveRDS(world_map_sp, file = "./input_data/Map_stuff/world_map_sp.rds", version = "2")
# 
# # Generate mask with extent and resolution of the DEM
# world_map_mask_90m <- rasterize(x = world_map_sp, y = DEM_SRTM_v4_90m, # Provide the grid to fill with CRS, bbox and resolution
#           field = 1, # How to fill non empty cells. With the value of a varibel in the df of the sp_obj, or directly with a fixed value ?
#           background = NA, # Value to use to fill empty cells
#           fun = mean) # What to do when multiple sp features overlap with a cell ?
# 
# saveRDS(world_map_mask_90m, file = "./input_data/Map_stuff/world_map_mask_90m.rds", version = "2")
# 
# 
# DEM_masked <- mask(DEM_SRTM_v4_90m, world_map_mask)
# 
# 
# x11()
# plot(world_map_mask)
# plot(DEM)
# plot(DEM_masked)


#### 3/ Compute Terrain Ruggedness Index (TRI) ####

library(spatialEco)

?tri

# # Version for 15' resolution
# DEM_TRI <- tri(r = DEM, s = 3, # size of the window
#                exact = TRUE, file.name = NULL)
# 
# plot(DEM_TRI)
# saveRDS(DEM_TRI, file = "./outputs/Indices_maps/DEM_TRI.rds", version = "2")
# 
# DEM_TRI_15m <- resample(DEM_TRI, Community_mask, method = "bilinear")
# saveRDS(DEM_TRI_15m, file = "./outputs/Indices_maps/DEM_TRI_15m.rds", version = "2")
# plot(DEM_TRI_15m)

# Version for 1' resolution, window size = 1 community
DEM_TRI_1m <- tri(r = DEM_SRTM_v4_1m, s = 15, # size of the window
               exact = TRUE, file.name = NULL)

plot(DEM_TRI_1m)
saveRDS(DEM_TRI_1m, file = "./outputs/Indices_maps/DEM_TRI_1m.rds", version = "2")

DEM_TRI_15m <- resample(DEM_TRI_1m, sp.richness, method = "bilinear")
saveRDS(DEM_TRI_15m, file = "./outputs/Indices_maps/DEM_TRI_15m.rds", version = "2")
plot(DEM_TRI_15m)


#### 4/ Mask communities outside of Ithomiini range ####

DEM_TRI <- readRDS(file = "./outputs/Indices_maps/DEM_TRI.rds")
DEM_TRI_15m <- readRDS(file = "./outputs/Indices_maps/DEM_TRI_15m.rds")

# Create a mask of communities with presence of Ithomiini
Community_mask <- sp.richness
Community_mask[!(Community_mask > 0)] <- NA
Community_mask[Community_mask > 0] <- 1

sum(Community_mask[], na.rm = T) # 21 415 communities
plot(Community_mask)
save(Community_mask, file = "input_data/Map_stuff/Community_mask.RData")

# Mask all rasters to constrain to Ithomiini communities
DEM_masked <- mask(DEM, world_map_mask)
DEM_TRI_masked <- mask(DEM_TRI_15m, Community_mask)
sp.richness_masked <- mask(sp.richness, Community_mask)
ring.richness_masked <- mask(ring.richness, Community_mask)

#### 5/ GAM model of mimicry richness ~ species richness => Extract residuals ####

library(mgcv)

sp.richness_values <- sp.richness_masked[which(!is.na(sp.richness_masked[]))]
ring.richness_values <- ring.richness_masked[which(!is.na(sp.richness_masked[]))]


mod_gam <- gam(ring.richness_values ~ s(sp.richness_values))
summary(mod_gam)

res_GAM <- residuals.gam(mod_gam, type = "pearson")

GAM_df <- data.frame(cbind(ring.richness_values, sp.richness_values, res_GAM))
saveRDS(GAM_df, file = "./outputs/Indices_maps/GAM_df.rds")

library(ggplot2)

GAM_df <- readRDS(file = "./outputs/Indices_maps/GAM_df.rds")

GAM_plot <-  ggplot(GAM_df, aes(x = sp.richness_values, y = ring.richness_values, color = res_GAM)) +
  geom_point() +
  scale_color_gradient2(midpoint = 0,  low = "darkblue", mid = "white",
                            high = "red")



print(GAM_plot)


#### 6/ Map residuals aside TRI ####

GAM_residuals_map <- sp.richness_masked
GAM_residuals_map[which(!is.na(sp.richness_masked[]))] <- res_GAM

saveRDS(GAM_residuals_map, file = "./outputs/Indices_maps/GAM_residuals_map.rds")

GAM_residuals_map <- readRDS(file = "./outputs/Indices_maps/GAM_residuals_map.rds")

par(mfrow=(c(1,2)))
plot(GAM_residuals_map)
plot(DEM_TRI_masked)

par(mfrow=(c(1,1)))
plot(x = DEM_TRI_masked, y = GAM_residuals_map, xlab = "TRI", ylab = "Residuals")

#### 7/ Extract TRI among bioregions ####

# ### Extract via shapes
# # Load shp
# Central_America_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Central_America_shp.rds")
# East_Ecuador_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/East_Ecuador_shp.rds")
# Peruvian_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Peruvian_shp.rds")
# Mata_Atlantica_shp5 <- readRDS(file = "./input_data/Map_stuff/Bioregions/Mata_Atlantica_shp5.rds")
# 
# # Remove communities above 3000m
# mask_3000m <- DEM_15m <- readRDS(file = paste0("./input_data/SRTM/Elevation_15.rds"))
# mask_3000m[DEM_15m[] > 3000] <- NA
# mask_3000m[DEM_15m[] <= 3000] <- 1
# 
# # With communities above 3000m
# # TRI_CA <- unlist(raster::extract(DEM_TRI_masked, Central_America_shp))
# # TRI_EE <- unlist(raster::extract(DEM_TRI_masked, East_Ecuador_shp))
# # TRI_EP <- unlist(raster::extract(DEM_TRI_masked, Peruvian_shp))
# # TRI_MA <- unlist(raster::extract(DEM_TRI_masked, Mata_Atlantica_shp5))
# 
# # Without communities above 3000m
# DEM_TRI_masked_3000m <- DEM_TRI_masked*mask_3000m
# plot(DEM_TRI_masked_3000m)
# 
# TRI_CA <- unlist(raster::extract(DEM_TRI_masked_3000m, Central_America_shp))
# TRI_EE <- unlist(raster::extract(DEM_TRI_masked_3000m, East_Ecuador_shp))
# TRI_EP <- unlist(raster::extract(DEM_TRI_masked_3000m, Peruvian_shp))
# TRI_MA <- unlist(raster::extract(DEM_TRI_masked_3000m, Mata_Atlantica_shp5))


### Extract via indices
load(file = "./outputs/Correlation_tests/Bioregions_indices.RData")

TRI_CA <- DEM_TRI_masked[Central_America_index]
TRI_EE <- DEM_TRI_masked[East_Ecuador_index]
TRI_EP <- DEM_TRI_masked[Peruvian_index]
TRI_MA <- DEM_TRI_masked[Mata_Atlantica_index]
TRI_others <- DEM_TRI_masked[][-c(Central_America_index, East_Ecuador_index, Peruvian_index, Mata_Atlantica_index)]
TRI_others <- TRI_others[!is.na(TRI_others)]


# Generate df
TRI_bioregions <- c(TRI_CA, TRI_EE, TRI_EP, TRI_MA)
Name_bioregions <- c(rep("Central_America", times = length(TRI_CA)), rep("East_Ecuador", times = length(TRI_EE)), rep("East_Peru", times = length(TRI_EP)), rep("Mata_Atlantica", times = length(TRI_MA)))
TRI_df <- data.frame(cbind(TRI_bioregions, Name_bioregions))
TRI_df$TRI_bioregions <- as.numeric(as.character(TRI_df$TRI_bioregions))

saveRDS(TRI_df, file = "./outputs/Indices_maps/TRI_df.rds", version = "2")

TRI_df <- readRDS(file = "./outputs/Indices_maps/TRI_df.rds")

boxplot(data = TRI_df, TRI_bioregions ~ Name_bioregions, col = gg_color_hue(4))

library(tidyverse)
library(ggsignif)
library(ggpubr)

pdf(file = "./graphs/TRI/Boxplot_TRI_Bioregions.pdf", height = 8, width = 10)


gg_boxplot <- ggplot(data = TRI_df, aes(x = Name_bioregions, y = TRI_bioregions)) +
  geom_boxplot(fill = gg_color_hue(4)) + 
  ggtitle("") +
  xlab("") + ylab("Terrain Ruggedness Index") +
  scale_x_discrete(labels = c("Central America", "East Ecuador", "East Peru", "BAF")) +
  scale_y_continuous(limits = c(0, 20000), expand = c(0, 0), breaks = seq(0, 20000, 5000)) +
  # stat_compare_means(mapping = aes(x = Name_bioregions, y = TRI_bioregions), comparisons = list(c("East_Ecuador", "East_Peru"), c("Central_America", "East_Ecuador")), hide.ns = FALSE,
  #                    label = "p.signif",  label.x = NULL, label.y = NULL) # +
  geom_signif(comparisons = list(c("East_Ecuador", "East_Peru") #,
                                 # c("Central_America", "East_Ecuador"),
                                 # c("Central America", "East_Peru"),
                                 # c("Central_America", "BAF"),
                                 # c("East_Ecuador", "BAF"),
                                 # c("East_Peru", "BAF")
                                 ),
              map_signif_level = T,
              textsize = 10,
              size = 1, margin_top = 0.15) +
  theme(panel.background = element_rect(fill = NA),
        legend.position = c(0.85, 0.55),
        # legend.key = element_rect(fill = "white", colour = "white"),
        legend.title = element_text(size = 15, vjust = 3, face = "bold"),
        legend.text = element_text(size = 14, face = "plain"),
        # margin = margin(t = 0, unit = "pt")),
        legend.key.size = unit(2, 'lines'),
        # legend.spacing.y = unit(1,"cm"),
        # axis.line.x = element_line(color = NA),
        axis.line.y = element_line(color = "black", size = 1.3),
        axis.ticks = element_line(color = "black"),
        axis.ticks.x = element_line(size = NA),
        axis.ticks.y = element_line(size = 1.2),
        axis.ticks.length.y = unit(8, "pt"),
        axis.text = element_text(size = 14, face = "bold", color = "black"),
        axis.text.x = element_text(vjust = -4, size = 14, margin = margin(t = -5, b = 20)),
        axis.title = element_text(size = 18, face = "bold", color = "black"),
        # axis.title.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title.y = element_text(margin = margin(l = 5, r = 10)))
  
print(gg_boxplot)

dev.off()


#### 8/ Test differences of TRI among bioregions ####

hist(TRI_CA)
shapiro.test(TRI_CA)

hist(TRI_EE)
shapiro.test(TRI_EE)

model_anova <- glm(TRI_bioregions ~ Name_bioregions, family = poisson(link="log"))
summary(model_anova)
par(mfrow = c(2,2))
plot(model_anova)

pairwise.t.test(x = TRI_bioregions, g = Name_bioregions, p.adjust.method = "bonferroni", pool.sd = F) # Tests post-hoc par paires sur un lm # See p.adjust() for information about methods of correction

wilcox.test(x = TRI_EE, y = TRI_EP, alternative = "greater", paired = F)
# W = 14757, p-value = 1

#### 9/ Plot TRI and bioregions ####

Central_America_shp2 <- readRDS(file = "./input_data/Map_stuff/Bioregions/Central_America_shp2.rds")
East_Ecuador_shp2 <- readRDS(file = "./input_data/Map_stuff/Bioregions/East_Ecuador_shp2.rds")
Peruvian_shp2 <- readRDS(file = "./input_data/Map_stuff/Bioregions/Peruvian_shp2.rds")
Mata_Atlantica_shp5 <- readRDS(file = "./input_data/Map_stuff/Bioregions/Mata_Atlantica_shp5.rds")



plot(DEM_TRI_masked)
plot(crop(DEM_TRI_masked, extent(c(-90, -60, -20, 10))))
plot(crop(DEM_TRI_masked_3000m, extent(c(-90, -60, -20, 10))))
plot(Central_America_shp2, border = gg_color_hue(4)[1], col = NA, lwd = 3, add = T)
plot(East_Ecuador_shp2, border = gg_color_hue(4)[2], col = NA, lwd = 3, add = T)
plot(Peruvian_shp2, border = gg_color_hue(4)[3], col = NA, lwd = 3, add = T)
plot(Mata_Atlantica_shp5, border = gg_color_hue(4)[4], col = NA, lwd = 3, add = T)


?extent

extent(c(-90, -70, -20, 10))


#### 10/ Test Spatial correlation between residuals and TRI ####

test_map <- DEM_TRI

test_map[] <- 0
test_map[Central_America_index] <- 1

plot(test_map)


#### 11/ Extract GAM residuals per bioregions and compare them #####

### Extract via shapes
# # Load shp
# Central_America_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Central_America_shp.rds")
# East_Ecuador_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/East_Ecuador_shp.rds")
# Peruvian_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Peruvian_shp.rds")
# Mata_Atlantica_shp5 <- readRDS(file = "./input_data/Map_stuff/Bioregions/Mata_Atlantica_shp5.rds")
# 
# GAM_resid_CA <- unlist(raster::extract(GAM_residuals_map, Central_America_shp))
# GAM_resid_EE <- unlist(raster::extract(GAM_residuals_map, East_Ecuador_shp))
# GAM_resid_EP <- unlist(raster::extract(GAM_residuals_map, Peruvian_shp))
# GAM_resid_MA <- unlist(raster::extract(GAM_residuals_map, Mata_Atlantica_shp5))

### Extract via indices
load(file = "./outputs/Correlation_tests/Bioregions_indices.RData")

GAM_resid_CA <- GAM_residuals_map[Central_America_index]
GAM_resid_EE <- GAM_residuals_map[East_Ecuador_index]
GAM_resid_EP <- GAM_residuals_map[Peruvian_index]
GAM_resid_MA <- GAM_residuals_map[Mata_Atlantica_index]
GAM_resid_others <- GAM_residuals_map[][-c(Central_America_index, East_Ecuador_index, Peruvian_index, Mata_Atlantica_index)]
GAM_resid_others <- GAM_resid_others[!is.na(GAM_resid_others)]

# Generate df
GAM_resid_bioregions <- c(GAM_resid_CA, GAM_resid_EE, GAM_resid_EP, GAM_resid_MA, GAM_resid_others)
Name_bioregions <- c(rep("Central_America", times = length(GAM_resid_CA)), rep("East_Ecuador", times = length(GAM_resid_EE)), rep("East_Peru", times = length(GAM_resid_EP)), rep("Mata_Atlantica", times = length(GAM_resid_MA)), rep("Others", times = length(GAM_resid_others)))
GAM_resid_df <- data.frame(cbind(GAM_resid_bioregions, Name_bioregions))
GAM_resid_df$GAM_resid_bioregions <- as.numeric(as.character(GAM_resid_df$GAM_resid_bioregions))
GAM_resid_df$Name_bioregions <- factor(GAM_resid_df$Name_bioregions, levels = c("Others", "Mata_Atlantica", "Central_America", "East_Ecuador", "East_Peru"))
saveRDS(GAM_resid_df, file = "./outputs/Indices_maps/GAM_resid_df.rds", version = "2")

GAM_resid_df <- readRDS(file = "./outputs/Indices_maps/GAM_resid_df.rds")



#### 12/ Boxplot of GAM residuals per bioregions #####

boxplot(data = GAM_resid_df, GAM_resid_bioregions ~ Name_bioregions, col = c("grey", gg_color_hue(4)))

library(tidyverse)
library(ggsignif)
library(ggpubr)

?geom_signif

# Select the color to keep the same than in other graph with bioregions
gg_color_hue(4)
color_list <- c("grey", "#C77CFF", "#F8766D", "#7CAE00", "#00BFC4")

pdf(file = "./graphs/TRI/Boxplot_GAM_resid_Bioregions.pdf", height = 8, width = 10)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,0,0,0), mar = c(4.5,4,3.5,1.3), xpd = T)

# Set margin parameters: Bottom, left, top, right



gg_boxplot <- ggplot(data = GAM_resid_df, aes(x = Name_bioregions, y = GAM_resid_bioregions)) +
  geom_boxplot(fill = color_list) + 
  ggtitle("") +
  xlab("") + ylab("GAM Residuals") +
  scale_x_discrete(labels = c("Others", "BAF", "Central America", "East Ecuador", "East Peru")) +
  # scale_y_continuous(limits = c(-8, 12), expand = c(0, 0), breaks = seq(-8, 12, 4)) +
  scale_y_continuous(limits = c(-6, 12), expand = c(0, 0), breaks = seq(-6, 12, 2)) +
  geom_hline(yintercept = 0, linetype = "dashed", 
             color = "black", size = 1) +
  # stat_compare_means(mapping = aes(x = Name_bioregions, y = GAM_resid_bioregions), comparisons = list(c("East_Ecuador", "East_Peru"), c("Central_America", "East_Ecuador")), hide.ns = FALSE,
  #                    label = "p.signif",  label.x = NULL, label.y = NULL) # +
  geom_signif(comparisons = list(c("East_Ecuador", "East_Peru"),
                                 c("Others", "Mata_Atlantica"),
                                 c("Others", "Central_America") ), #,
                                 # c("Central America", "East_Peru") #,
                                 # c("Central_America", "BAF") #,
                                 # c("East_Ecuador", "BAF") #,
                                 # c("East_Peru", "BAF")
  map_signif_level = T,
  textsize = 10,
  size = 1, margin_top = c(0.07, 0.07, 0.15)) +
  theme(panel.background = element_rect(fill = NA),
        legend.position = c(0.85, 0.55),
        # legend.key = element_rect(fill = "white", colour = "white"),
        legend.title = element_text(size = 15, vjust = 3, face = "bold"),
        legend.text = element_text(size = 14, face = "plain"),
        # margin = margin(t = 0, unit = "pt")),
        legend.key.size = unit(2, 'lines'),
        # legend.spacing.y = unit(1,"cm"),
        # axis.line.x = element_line(color = NA),
        axis.line.y = element_line(color = "black", size = 1.3),
        axis.ticks = element_line(color = "black"),
        axis.ticks.x = element_line(size = NA),
        axis.ticks.y = element_line(size = 1.2),
        axis.ticks.length.y = unit(8, "pt"),
        axis.text = element_text(size = 14, face = "bold", color = "black"),
        axis.text.x = element_text(vjust = -4, size = 16, margin = margin(t = -2, b = 20)),
        axis.title = element_text(size = 18, face = "bold", color = "black"),
        # axis.title.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title.y = element_text(margin = margin(l = 5, r = 10)))

print(gg_boxplot)

par(oma = original_ext_margins, mar = original_int_margins, xpd = F)  # Reload old_settings

dev.off()


#### 13/ Test differences of GAM residuals among bioregions ####

hist(GAM_resid_CA)
shapiro.test(GAM_resid_CA)

hist(GAM_resid_EE)
shapiro.test(GAM_resid_EE)

wilcox.test(x = GAM_resid_others, y = GAM_resid_CA, alternative = "greater", paired = F)
# W = 2132984, p-value << 0.001

wilcox.test(x = GAM_resid_others, y = GAM_resid_MA, alternative = "greater", paired = F)
# W = 27268335, p-value << 0.001

wilcox.test(x = GAM_resid_EE, y = GAM_resid_EP, alternative = "greater", paired = F)
# W = 32564, p-value << 0.001


#### Plot Residuals ~ TRI with bioregions colors for points


par(mfrow=(c(1,1)))
plot(x = DEM_TRI_masked, y = GAM_residuals_map, xlab = "TRI", ylab = "Residuals")

TRI_all <- c(TRI_df$TRI_bioregions, TRI_others)

plot(x = rev(TRI_all), y = rev(GAM_resid_df$GAM_resid_bioregions), xlab = "TRI", ylab = "Residuals", pch = 16, col = c("grey", gg_color_hue(4))[rev(GAM_resid_df$Name_bioregions)])

