
##### Explore diversity values for specific coordinates #####

library(raster)

sp_richness <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
ring_richness <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))

load(file = "./input_data/Map_stuff/country_borders.RData")
country_borders <- as(country_borders, "Spatial")



# Jatun Sacha = 
target_coords <- c(-77.61664182271451, -1.0869885282012113)
target_coords <- c(-77.5, -1.087)
target_coords <- c(-77.4, -1.087)

raster::cellFromXY(object = sp_richness, xy = target_coords)

raster::extract(x = sp_richness, y = as.matrix(t(target_coords)))

pdf(file = "./Jatun_Sacha.pdf", width = 10, height = 8)

plot(sp_richness)
plot(country_borders, add = T)
points(x = target_coords[1], y  = target_coords[2], col = "red", pch = 16, cex = 0.5)

dev.off()

?extract

table(Ithomiini_final$DB_orig)
