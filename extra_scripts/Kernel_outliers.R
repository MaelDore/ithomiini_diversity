

library(raster)

### Load occurrence dataset
load(file = "./input_data/Ithomiini_final.RData")

# Load mask for continent borders
res <- "15"
continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_", res, ".rds"))

load(file = "./input_data/Map_stuff/bg_mask.RData") # Load bg shp
load(file = "./input_data/Map_stuff/country_borders.RData") # Load country borders
rivers <- st_read("./input_data/Map_stuff/Major_rivers/MajorRivers.shp") # Load rivers

# New color palette
pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")

# Outer borders
crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".rds"))


# Save occurence map
pdf(file = paste0("./maps/All_occurrences.pdf"), height = 5.3, width = 6.5)
internal_margins <- par()$mar
par(mar = c(3.1, 3.5,3.5,2.1))
image(continent_mask, col = "antiquewhite", main = paste0("All occurrences"), 
      cex.axis = 1.4, cex.main = 1.4, axis.args=list(cex.axis=1.4),
      # asp = 1, colNA = "aliceblue",
      # xlim = tot.sp.richness_Jaccard.80@extent[1:2], ylim = tot.sp.richness_Jaccard.80@extent[3:4],
      ylab = "", xlab = "",
      legend.args=list(text="          Species", cex=1.2, line = 1, font = 2), 
      legend  = F)
points(Ithomiini_final$Longitude, Ithomiini_final$Latitude, pch = 16, cex = 0.5)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
plot(rivers, lwd = 0.05, col = "#75b2f6", alpha = 1, add = T)
plot(country_borders, lwd = 0.2, border = "red", color = "", add = T)
plot(bg_mask, lwd = 1, border = "grey20", col = "aliceblue", add = T)
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-95, -33), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1)
prettymapr::addnortharrow(scale = 0.75, padin = c(0.2, 0.2))
legend(legend = c("Borders", "Rivers"), lty = 1, col = c("red", "#75b2f6"), x = "bottomleft", bg = "white")
par(mar = internal_margins)
dev.off()


### Kernel stuff to detect outliers from Gomes et al., 2018

# Select unit to apply the kernel

# occ.unit <- Ithomiini_final # All occurrences
Tag.unit <- "Aeria.elara.EURIMEDIA.EURIMEDIA"
Tag.unit <- "Mechanitis.menapis.MAMERCUS.MAMERCUS"

occ.unit <- Ithomiini_final %>%
  dplyr::filter(Tag == Tag.unit)

plot(continent_mask)
points(occ.unit$Longitude, occ.unit$Latitude, pch = 16, cex = 0.5)

library(spatstat)

?density.ppp 
?ppp
?owin

# Need to transform the numeric mask in a logical matrix
window_mask <- continent_mask
window_mask[continent_mask[] == 0] <- T
window_mask[is.na(continent_mask[])] <- F

# Since windows must have rows in increasing order from top to bottom, need to reverse the rows
window_mask <- apply(matrix(data = as.logical(as.matrix(window_mask[])), ncol = continent_mask@ncols, byrow = T), 2, rev)

# Generate the window
window_obs <- owin(xrange = continent_mask@extent[1:2], yrange = continent_mask@extent[3:4], mask = as.matrix(window_mask[]))
str(window_obs) # Check xstep/ystep to be 0.25°

# Generate a point pattern dataset from occurrence specifiying the window of potential occurrences (= the continental mask)
ppp_obs <- ppp(x = occ.unit$Longitude, y = occ.unit$Latitude, window = window_obs) 

# Compute the kernel smoothed intensity function
density_kernel_obs <- density(x = ppp_obs)
# Weights can be added such as robustness of sites attribution !

# Plot the intensity function
image(density_kernel_obs$xcol,density_kernel_obs$yrow,t(density_kernel_obs$v))

hist(density_kernel_obs$v)

# The result of density.ppp is not a probability density, otherwise values at occurrences points would be 1, and of no use to detect outliers.
# It is an estimate of the intensity function of the point process that generated the point pattern data. 
# Intensity is the expected number of random points per unit area. The units of intensity are “points per unit area”. 
# Intensity is usually a function of spatial location, and it is this function which is estimated by density.ppp. 
# The integral of the intensity function over a spatial region gives the expected number of points falling in this region.

# A good threshold could be 1, or 0.5, since in means following the distribution of occurrences, 
# it is expected to have less than 1 (or 0.5) data points in that pixel. 
# But for a low sample prevalence, (i.e., few pixels with an occurrence) the intensity function will be way lower than 0.5

# Option of Gomes et al., 2018 = use a quantile
# It is not a quantile of occurrences but a quantile of pixels on the map. Another option would be to use a quantile of the occurrence points
# but then we will have the issue of always discarding a portion of occurrences even if they are fine...

quantiles <- quantile(x = density_kernel_obs$v, probs = seq(0, 1, 0.1),  na.rm = T) ; quantiles

# Threshold for outliers
threshold <- 1
threshold <- quantiles[8] # Remove 30% of the map surface

# Store output in a rasterLayer
density_kernal_map <- continent_mask
density_kernal_map[] <- as.vector(t(apply(density_kernel_obs$v, 2, rev)))

image(density_kernal_map)
points(occ.unit$Longitude, occ.unit$Latitude, pch = 16, cex = 0.5)

# Binarization of the kernel smoothed intensity function
bin_density_kernal_map <- density_kernal_map
bin_density_kernal_map[bin_density_kernal_map < threshold] <- 0
bin_density_kernal_map[bin_density_kernal_map >= threshold] <- 1

image(bin_density_kernal_map)
points(occ.unit$Longitude, occ.unit$Latitude, pch = 16, cex = 0.5)
