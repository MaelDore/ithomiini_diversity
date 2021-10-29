
###### Script 23: Map bioregions ######

# Map bioregions used in the analyses (Figure S1.2)

###

# Inputs: 
   # WWF Ecoregions (Olson et al., 2001)

# Outputs:
    # Map of bioregions to improve in ArcGIS (Figure S1.2)

###


### 0/ Load stuff ####

# Clean environment
rm(list = ls())


# world_map_sp <- readRDS(file = "./input_data/Map_stuff/world_map_sp.rds")
# 
# plot(world_map_sp)
# 
# xmin = -120 ; xmax = -30 ; ymin = -37 ; ymax = 28
# e <- extent(c(xmin,xmax,ymin,ymax))
# 
# America_map_sp <- crop(world_map_sp, e)

library(sf)

load(file = "./input_data/Map_stuff/country_borders.RData") # Load country borders
# rivers <- st_read("./input_data/Map_stuff/Major_rivers/MajorRivers.shp") # Load rivers
# load(file = "./input_data/Map_stuff/bg_mask.RData") # Load bg shp

xmin = -120 ; xmax = -30 ; ymin = -37 ; ymax = 28
e <- extent(c(xmin,xmax,ymin,ymax))

DEM_SRTM_v4_1m <- readRDS(file = "./input_data/SRTM/DEM_SRTM_v4_1m.rds")
DEM_SRTM_v4_1m_cropped <- crop(DEM_SRTM_v4_1m, e)

pdf(file = paste0("./maps/Bioregions/Regions.pdf"), height = 12, width = 15)

plot(DEM_SRTM_v4_1m_cropped, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))
plot(country_borders, border = "#00000030", add = T)

dev.off()

##### 1/ Delimitate regions from WWF Ecoregions ####

library(sf)
library(tidyverse)
library(rgeos)
library(spatialEco)

### All ecoregions ####

Biomes_sf <- st_read("./input_data/WWF_official_teow/wwf_terr_ecos.shp")

xmin = -120 ; xmax = -30 ; ymin = -37 ; ymax = 28
e <- extent(c(xmin,xmax,ymin,ymax))

Neotropics <- Biomes_sf %>% 
  st_crop(e) %>%
  filter(BIOME != 98)

# All names used in the main text to locate stuff

### Guyana Shield ####
# 60173, 60125, 60707, 60147, 61401, 60182, 60124

Guyana_Shield <- Neotropics %>% 
  filter(ECO_ID %in% c(60173, 60125, 60707, 60147, 60182, 60124)) %>%
  summarize() 

Guyana_Shield_shp <- remove.holes(gSimplify(as(Guyana_Shield, 'Spatial'), tol = 0.6, topologyPreserve = F))

plot(DEM_SRTM_v4_1m_cropped, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))
plot(Guyana_Shield_shp, col = "purple", add = T)

saveRDS(Guyana_Shield_shp, file = "./input_data/Map_stuff/Bioregions/Guyana_Shield_shp.rds", version = "2")

### Caatinga ####
# 61304

Caatinga <- Neotropics %>% 
  filter(ECO_ID %in% c(61304)) %>%
  summarize() 

Caatinga_shp <- remove.holes(gSimplify(as(Caatinga, 'Spatial'), tol = 0.5, topologyPreserve = F))

plot(DEM_SRTM_v4_1m_cropped, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))
plot(Caatinga_shp, col = "red", add = T)

saveRDS(Caatinga_shp, file = "./input_data/Map_stuff/Bioregions/Caatinga_shp.rds", version = "2")


### Cerrado ####
# 60704

Cerrado <- Neotropics %>% 
  filter(ECO_ID %in% c(60704)) %>%
  summarize() 

Cerrado_shp <- remove.holes(gSimplify(as(Cerrado, 'Spatial'), tol = 1.2, topologyPreserve = F))

plot(DEM_SRTM_v4_1m_cropped, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))
plot(Cerrado_shp, col = "orange", add = T)

saveRDS(Cerrado_shp, file = "./input_data/Map_stuff/Bioregions/Cerrado_shp.rds", version = "2")


### Brazilian Atlantic Forest ####
# 60150, 60101, 60160, 60150, 60104, 60103, 60151, 60152
# Voir autre shp file tout fait

BAF <- Neotropics %>% 
  filter(ECO_ID %in% c(60150, 60101, 60160, 60150, 60104, 60103, 60151, 60152)) %>%
  summarize() 

BAF_shp <- remove.holes(gSimplify(as(BAF, 'Spatial'), tol = 0.5, topologyPreserve = F))

plot(DEM_SRTM_v4_1m_cropped, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))
plot(BAF, col = "darkgreen", add = T)
plot(BAF_shp, col = "limegreen", add = T)

# Only ecoregions in the tropical forests Biome
Mata_Atlantica_shp2 <- readRDS(file = "./input_data/Map_stuff/Bioregions/Mata_Atlantica_shp2.rds")
plot(Mata_Atlantica_shp2, col = "limegreen", add = T)

# From Terrabrasilis
Mata_Atlantica_shp4 <- readRDS(file = "./input_data/Map_stuff/Bioregions/Mata_Atlantica_shp4.rds")
plot(Mata_Atlantica_shp4, col = "limegreen", add = T)

# Merge both to add the Northern East portion, but keep the SouthWest one
Mata_Atlantica_shp5 <- gUnion(Mata_Atlantica_shp2, Mata_Atlantica_shp4)
plot(Mata_Atlantica_shp5, col = "limegreen", add = T)

saveRDS(Mata_Atlantica_shp5, file = "./input_data/Map_stuff/Bioregions/Mata_Atlantica_shp5.rds", version = "2")

### Central America ####
# 60129, 60167, 60119, 60130, 60224, 60209, 60303, 60154, 60181, 60235, 60230, 60211, 60309, ... Use Ithomiini range to get the last upper part

box_CA <- SpatialPolygons(list(Polygons(list(Polygon(coords = matrix(data = c(-116, 28, -79, 28, -79, 8, -116, 8, -116, 28),
                                                                                ncol = 2, byrow = T), hole = F))
                                                   , ID = 1))
                                     , proj4string = crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

plot(box_CA, col = "orange", add = T)

box_Caribbean_N <- SpatialPolygons(list(Polygons(list(Polygon(coords = matrix(data = c(-86, 28, -57, 28, -57, 17, -86, 17, -86, 28),
                                                                     ncol = 2, byrow = T), hole = F))
                                        , ID = 1))
                          , proj4string = crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(box_Caribbean_N, col = "orange", add = T)

Neotropics_shp <- as(Neotropics, 'Spatial')
Neotropics_shp@proj4string <- box_CA@proj4string

full_CA <- gDifference(gIntersection(Neotropics_shp, box_CA), box_Caribbean_N)

full_CA_shp <- remove.holes(gSimplify(full_CA, tol = 0.2, topologyPreserve = F))

plot(DEM_SRTM_v4_1m_cropped, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))
plot(full_CA_shp, col = "pink", add = T)

saveRDS(full_CA_shp, file = "./input_data/Map_stuff/Bioregions/full_CA_shp.rds", version = "2")


### Northern Andes ####
# 60175, 60118, 60136, 60221, 60109, 60207, 61006, 60145, 60225, 60121, 60159, 61007

Northern_Andes <- Neotropics %>% 
  filter(ECO_ID %in% c(60175, 60118, 60136, 60221, 60109, 60207, 61006, 60145, 60225, 60121, 60159, 61007)) %>%
  summarize() 

saveRDS(Northern_Andes, file = "./input_data/Map_stuff/Bioregions/Northern_Andes_raw_shp.rds", version = "2")

Northern_Andes_shp <- remove.holes(gSimplify(as(Northern_Andes, 'Spatial'), tol = 0.3, topologyPreserve = F))

plot(DEM_SRTM_v4_1m_cropped, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))
plot(Northern_Andes_shp, col = "orange", add = T)

saveRDS(Northern_Andes_shp, file = "./input_data/Map_stuff/Bioregions/Northern_Andes_shp.rds", version = "2")


### Central Andes ####
# 60153, 60223, 61003, 61002, 61001, 60105, 60206, 60165

Central_Andes <- Neotropics %>% 
  filter(ECO_ID %in% c(60153, 60223, 61003, 61002, 61001, 60105, 60206, 60165)) %>%
  summarize() 

saveRDS(Central_Andes, file = "./input_data/Map_stuff/Bioregions/Central_Andes_raw_shp.rds", version = "2")

Central_Andes_shp <- remove.holes(gSimplify(as(Central_Andes, 'Spatial'), tol = 0.5, topologyPreserve = F))

plot(DEM_SRTM_v4_1m_cropped, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))
plot(Central_Andes_shp, col = "orange", add = T)

saveRDS(Central_Andes_shp, file = "./input_data/Map_stuff/Bioregions/Central_Andes_shp.rds", version = "2")

### Western Lowlands ####
# 60232, 60214, 60178, 60115, 60137, 60229, 61308

Western_Lowlands <- Neotropics %>% 
  filter(ECO_ID %in% c(60232, 60214, 60178, 60115, 60137, 60229, 61308)) %>%
  summarize() 

saveRDS(Western_Lowlands, file = "./input_data/Map_stuff/Bioregions/Western_Lowlands_raw_shp.rds", version = "2")

Western_Lowlands_shp <- remove.holes(gSimplify(as(Western_Lowlands, 'Spatial'), tol = 0.2, topologyPreserve = F))

plot(DEM_SRTM_v4_1m_cropped, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))
plot(Western_Lowlands_shp, col = "limegreen", add = T)

saveRDS(Western_Lowlands_shp, file = "./input_data/Map_stuff/Bioregions/Western_Lowlands_shp.rds", version = "2")


### Coastal desert ####
# 61315

Coastal_desert <- Neotropics %>% 
  filter(ECO_ID %in% c(61315)) %>%
  summarize() 

Coastal_desert_shp <- remove.holes(gSimplify(as(Coastal_desert, 'Spatial'), tol = 0.2, topologyPreserve = F))

plot(DEM_SRTM_v4_1m_cropped, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))
plot(Coastal_desert_shp, col = "yellow", add = T)

saveRDS(Coastal_desert_shp, file = "./input_data/Map_stuff/Bioregions/Coastal_desert_shp.rds", version = "2")

### Llanos ####
# 60709

Llanos <- Neotropics %>% 
  filter(ECO_ID %in% c(60709)) %>%
  summarize() 

Llanos_shp <- remove.holes(gSimplify(as(Llanos, 'Spatial'), tol = 0.2, topologyPreserve = F))

plot(DEM_SRTM_v4_1m_cropped, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))
plot(Llanos_shp, col = "yellow", add = T)

saveRDS(Llanos_shp, file = "./input_data/Map_stuff/Bioregions/Llanos_shp.rds", version = "2")


### Chacos ####
# 60210, 60708

Chacos <- Neotropics %>% 
  filter(ECO_ID %in% c(60210, 60708)) %>%
  summarize() 

Chacos_shp <- remove.holes(gSimplify(as(Chacos, 'Spatial'), tol = 0.3, topologyPreserve = F))

plot(DEM_SRTM_v4_1m_cropped, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))
plot(Chacos_shp, col = "blue", add = T)

saveRDS(Chacos_shp, file = "./input_data/Map_stuff/Bioregions/Chacos_shp.rds", version = "2")

### Caribbean Islands ####
# 60123, 60120, 60902, 60304, 60213, 61402, 60131, 60218, 60127, 60305, 60215, 60155, 60226 + crop rectangle around islands

box_Florida <- SpatialPolygons(list(Polygons(list(Polygon(coords = matrix(data = c(-83, 28, -79, 28, -79, 25, -83, 25, -83, 28),
                                                                              ncol = 2, byrow = T), hole = F))
                                                 , ID = 1))
                                   , proj4string = crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(box_Florida, col = "orange", add = T)

box_Caribbean_S <- SpatialPolygons(list(Polygons(list(Polygon(coords = matrix(data = c(-64, 19, -57, 19, -57, 11, -64, 11, -64, 19),
                                                                              ncol = 2, byrow = T), hole = F))
                                                 , ID = 1))
                                   , proj4string = crs("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
plot(box_Caribbean_S, col = "orange", add = T)

box_Caribbean <- gDifference(gUnion(box_Caribbean_N, box_Caribbean_S), box_Florida)

plot(DEM_SRTM_v4_1m_cropped, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))
plot(box_Caribbean, col = "orange", add = T)

Caribbean_Islands <- gIntersection(Neotropics_shp, box_Caribbean)

Caribbean_Islands_shp <- remove.holes(gSimplify(Caribbean_Islands, tol = 0.05, topologyPreserve = F))

plot(DEM_SRTM_v4_1m_cropped, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))
plot(Caribbean_Islands_shp, col = "pink", add = T)

saveRDS(Caribbean_Islands_shp, file = "./input_data/Map_stuff/Bioregions/Caribbean_Islands_shp.rds", version = "2")

### Western/Upper Amazon ####
# 60142 (Napo), 60174 (Ucayali), 60128 (Iquitos), 60166 (SW), 60107 (Caqueta), 60163 (Solimões-Japura) 

Western_Amazon <- Neotropics %>% 
  filter(ECO_ID %in% c(60142, 60174, 60128, 60166, 60107, 60163)) %>%
  summarize() 

saveRDS(Western_Amazon, file = "./input_data/Map_stuff/Bioregions/Western_Amazon_raw_shp.rds", version = "2")

Western_Amazon_shp <- remove.holes(gSimplify(as(Western_Amazon, 'Spatial'), tol = 1.0, topologyPreserve = F))

plot(DEM_SRTM_v4_1m_cropped, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))
plot(Western_Amazon_shp, col = "blue", add = T)

saveRDS(Western_Amazon_shp, file = "./input_data/Map_stuff/Bioregions/Western_Amazon_shp.rds", version = "2")

### Lower Amazon ####
# 60132 (Solimões-Japura-Negro), 60143 (Negro Blanco), 60156 (Purus), 60133 (Jurua-Purus), 60157 (Purus-Madeira), 60135 (Madeira-Tapajos), 60140 (Mato Grosso),
# 60168 (Tapajos-Xingu), 60173 (Uatuma-Trombetas), 60180 (Xingu-Tocantins-Araguaia), 60138 (Marajo), 60170 (Tocantins/Pindare), 60139 (Maranhão), 60144

Lower_Amazon <- Neotropics %>% 
  filter(ECO_ID %in% c(60132, 60143, 60156, 60133, 60157, 60135, 60140, 60168, 60173, 60180, 60138, 60170, 60139, 60144)) %>%
  summarize() 

Lower_Amazon_shp <- remove.holes(gSimplify(as(Lower_Amazon, 'Spatial'), tol = 0.5, topologyPreserve = F))

plot(DEM_SRTM_v4_1m_cropped, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))
plot(Lower_Amazon_shp, col = "blue", add = T)

saveRDS(Lower_Amazon_shp, file = "./input_data/Map_stuff/Bioregions/Lower_Amazon_shp.rds", version = "2")


### Pantanal wetlands + Chiquitano dry forests
# 60907 (Pantanal flooded savannas), 60212 (Chiquitano dry forests)

Pantanal <- Neotropics %>% 
  filter(ECO_ID %in% c(60907, 60212)) %>%
  summarize() 

Pantanal_shp <- remove.holes(gSimplify(as(Pantanal, 'Spatial'), tol = 0.5, topologyPreserve = F))

plot(DEM_SRTM_v4_1m_cropped, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))
plot(Pantanal_shp, col = "blue", add = T)

saveRDS(Pantanal_shp, file = "./input_data/Map_stuff/Bioregions/Pantanal_shp.rds", version = "2")


### Pampas grasslands
# 60710 (Uruguayan savannas), 60801 (Espinal), 60803 (Humid Pampas)
# 60908 (Parana flodded savannas) # Can be added to join the two portions of the real pampas

Pampas <- Neotropics %>% 
  filter(ECO_ID %in% c(60710, 60801, 60803)) %>%
  summarize() # %>% 
  # st_cast("POLYGON") %>% # Divide Multipolygon into several polygons
  # slice(5) # Extract only the 5th one that fall into Ithomiini range

Pampas_shp <- remove.holes(gSimplify(as(Pampas, 'Spatial'), tol = 0.3, topologyPreserve = F))

plot(DEM_SRTM_v4_1m_cropped, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))
plot(Pampas_shp, col = "blue", add = T)

saveRDS(Pampas_shp, file = "./input_data/Map_stuff/Bioregions/Pampas_shp.rds", version = "2")

##### 2/ Plot map with all regions to delimitate on Illustrator ####

# Load all shp files

Guyana_Shield_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Guyana_Shield_shp.rds")
Caatinga_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Caatinga_shp.rds")
Cerrado_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Cerrado_shp.rds")
Mata_Atlantica_shp5 <- readRDS(file = "./input_data/Map_stuff/Bioregions/Mata_Atlantica_shp5.rds")
full_CA_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/full_CA_shp.rds")
Northern_Andes_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Northern_Andes_shp.rds")
Central_Andes_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Central_Andes_shp.rds")
Western_Lowlands_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Western_Lowlands_shp.rds")
Coastal_desert_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Coastal_desert_shp.rds")
Chacos_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Chacos_shp.rds")
Llanos_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Llanos_shp.rds")
Caribbean_Islands_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Caribbean_Islands_shp.rds")
Western_Amazon_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Western_Amazon_shp.rds")
Lower_Amazon_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Lower_Amazon_shp.rds")
Pantanal_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Pantanal_shp.rds")
Pampas_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Pampas_shp.rds")

# Save in a single file
save(Guyana_Shield_shp, Caatinga_shp, Cerrado_shp, Mata_Atlantica_shp5, full_CA_shp, Northern_Andes_shp, Central_Andes_shp, Western_Lowlands_shp, Coastal_desert_shp, Chacos_shp, Llanos_shp, Caribbean_Islands_shp, Western_Amazon_shp, Lower_Amazon_shp, Pantanal_shp, Pampas_shp, file = "./input_data/Map_stuff/Bioregions/All_bioregions_in_figure.RData")

# Load country borders
load(file = "./input_data/Map_stuff/country_borders.RData") 

# Plot map

pdf(file = paste0("./maps/Bioregions/Regions.pdf"), height = 12, width = 15)

plot(DEM_SRTM_v4_1m_cropped, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))
plot(country_borders, border = "#00000030", add = T)
plot(Guyana_Shield_shp, col = "bisque", add = T)
plot(Caatinga_shp, col = "chocolate1", add = T)
plot(Cerrado_shp, col = "darkgoldenrod1", add = T)
plot(Mata_Atlantica_shp5, col = "darkgreen", add = T)
plot(full_CA_shp, col = "firebrick1", add = T)
plot(Northern_Andes_shp, col = "azure2", add = T)
plot(Central_Andes_shp, col = "darkgrey", add = T)
plot(Western_Lowlands_shp, col = "limegreen", add = T)
plot(Coastal_desert_shp, col = "yellow", add = T)
plot(Llanos_shp, col = "darkviolet", add = T)
plot(Chacos_shp, col = "deeppink", add = T)
plot(Caribbean_Islands_shp, col = "pink", add = T)
plot(Western_Amazon_shp, col = "cyan", add = T)
plot(Lower_Amazon_shp, col = "blue", add = T)
plot(Pantanal_shp, col = "grey100", add = T)
plot(Pampas_shp, col = "tan4", add = T)

dev.off()


grid_Mollweide_out <- readRDS(file = "./input_data/Map_stuff/grid_Mollweide_out.rds")
large_bg_mask_Mollweide <- readRDS(file = "./input_data/Map_stuff/large_bg_mask_Mollweide.rds")
bbox_sp_Mollweide <- readRDS(file = "./input_data/Map_stuff/bbox_sp_Mollweide.rds")

DEM_SRTM_v4_1m_cropped_Mollweide <- projectRaster(from = DEM_SRTM_v4_1m_cropped, method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                                  
                                                  crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                                  alignOnly = F)

plot(DEM_SRTM_v4_1m_cropped_Mollweide, col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL))

pdf(file = paste0("./maps/Regions_background_Mollweide.pdf"), height = 12, width = 15)

image(DEM_SRTM_v4_1m_cropped_Mollweide, 
      col = gray.colors(100, start = 0.9, end = 0, gamma = 2.2, alpha = NULL),
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "",
      legend = F)

axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), padj = 0.5, labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 2, lwd = 0.2, lwd.ticks = 2)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 2, lwd = 0.2, lwd.ticks = 2)

plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)

scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-3900, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 2)
prettymapr::addnortharrow(scale = 0.9, padin = c(0.45, 0.45), text.col = "#00000000")

dev.off()


