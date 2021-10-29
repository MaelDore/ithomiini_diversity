
##### Break_down the Venter's 2009 HF to each of its dimensions #####

# Clean environment
rm(list = ls())

# Load continental_mask as template for final resolution and extent (initial one is too big to handle)
continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_15.rds"))

##### 1/ Load each dimensions of the Human Footprint #####

### 1.1/ Built areas ####

Built_2009 <- raster("./input_data/HII/HumanFootprint_break_down/Maps/Built2009.tif")

Built_2009 <- readAll(Built_2009)
table(sample(Built_2009@data@values, 100000))

plot(Built_2009)
crs(Built_2009) # units in m
res(Built_2009) # 1km x 1km

# Built areas = 10
Built_2009@data@attributes

# Test of resampling of data to final extent, CRS and resolution

# Built_2009_WGS84 <- projectRaster(from = Built_2009, to = continent_mask, method = "ngb")

res(continent_mask)
Built_2009_WGS84 <- projectRaster(from = Built_2009, res = 0.25/3, crs = continent_mask@crs, method = "ngb")

plot(Built_2009_WGS84)

Built_2009_Neotropics_5min <- raster::resample(x = Built_2009_WGS84, y = continent_mask, method = "ngb")

plot(Built_2009_Neotropics_5min)


### 1.2/ Crop lands ####

Croplands_2005 <- raster("./input_data/HII/HumanFootprint_break_down/Maps/croplands2005.tif")

plot(Croplands_2005)
crs(Croplands_2005) # units in m
res(Croplands_2005) # 1km x 1km

# Crop lands = 7
Croplands_2005@data@attributes


### 1.3/ Pasture lands ####

Pasture_2009 <- raster("./input_data/HII/HumanFootprint_break_down/Maps/Pasture2009.tif")

plot(Pasture_2009)
crs(Pasture_2009) # units in m
res(Pasture_2009) # 1km x 1km

# Pasture lands = 0 to 4 based on % of cover per pixel
Pasture_2009@data@attributes


### 1.4/ Light pollution ####

Lights_2009 <- raster("./input_data/HII/HumanFootprint_break_down/Maps/Lights2009.tif")

plot(Lights_2009)
crs(Lights_2009) # units in m
res(Lights_2009) # 1km x 1km

# Light pollution = equal quintile (from 1994 reference) on a 0 - 10 scale
Lights_2009@data@attributes


### 1.5/ Navigable waterways ####

NavWater_2009 <- raster("./input_data/HII/HumanFootprint_break_down/Maps/NavWater2009.tif")

plot(NavWater_2009)
crs(NavWater_2009) # units in m
res(NavWater_2009) # 1km x 1km

# Closeness to navigable waterways = 4, then decaying exponentially out to 15km
table(sample(x = readAll(NavWater_2009)@data@values, size = 100000))

hist(sample(x = readAll(NavWater_2009)@data@values, size = 100000))

### 1.6/ Population density ####

PopDensity_2010 <- raster("./input_data/HII/HumanFootprint_break_down/Maps/Popdensity2010.tif")

plot(PopDensity_2010)
crs(PopDensity_2010) # units in m
res(PopDensity_2010) # 1km x 1km

# Pressure score = 10 for more than 1000/people/km²
# Pressure score = 3.333 * log(population density + 1) for less than 1000/people/km²
# Then pressure score is round on a scale of 10
PopDensity_2010@data@attributes


### 1.7/ Railways ####

Railways <- raster("./input_data/HII/HumanFootprint_break_down/Maps/Railways.tif")

plot(Railways)
crs(Railways) # units in m
res(Railways) # 1km x 1km

# Pressure score = 8 within a range of 0.5km
Railways@data@attributes

### 1.8/ Roads ####

Roads <- raster("./input_data/HII/HumanFootprint_break_down/Maps/Roads.tif")

plot(Roads)
crs(Roads) # units in m
res(Roads) # 1km x 1km

# Closeness to navigable waterways = 8 within 0.5km, then decaying exponentially out from a score of 4 at 0.5km to 0 at 15km
table(sample(x = readAll(Roads)@data@values, size = 100000))
hist(sample(x = readAll(Roads)@data@values, size = 100000))


##### 2/ Retrieve 2009's Human Footprint = sum of all the 8 dimensions ####

### 2.1/ Load the original from paper's data ####

HF_2009_original <- raster("./input_data/HII/HumanFootprint_break_down/Maps/HFP2009.tif")

plot(HF_2009_original)
crs(HF_2009_original) # units in m
res(HF_2009_original) # 1km x 1km

# Max = 50
table(sample(x = readAll(HF_2009_original)@data@values, size = 100000))
hist(sample(x = readAll(HF_2009_original)@data@values, size = 100000))

# Resample to Neotropics at 15min resolution
continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_15.rds"))
HF_2009_original_WGS84 <- projectRaster(from = HF_2009_original, res = 0.25, crs = continent_mask@crs, method = "bilinear")
HF_2009_original_Neotropics_15min <- raster::resample(x = HF_2009_original_WGS84, y = continent_mask, method = "bilinear")
plot(HF_2009_original_Neotropics_15min)

saveRDS(HF_2009_original_Neotropics_15min, file = "./input_data/HII/HumanFootprint_break_down/Maps/HF_2009_original_Neotropics_15min.rds")
HF_2009_original_Neotropics_15min <- readRDS(file = "./input_data/HII/HumanFootprint_break_down/Maps/HF_2009_original_Neotropics_15min.rds")


### 2.2/ Rebuilt the 2009's HF from dimensions ####

HF_2009_stack <- stack(Built_2009, Croplands_2005, Lights_2009, NavWater_2009, Pasture_2009, PopDensity_2010, Railways, Roads)
HF_2009_stack <- readAll(HF_2009_stack) 

saveRDS(HF_2009_stack, file = "./input_data/HII/HumanFootprint_break_down/Maps/HF_2009_stack.rds")
HF_2009_stack <- readRDS(file = "./input_data/HII/HumanFootprint_break_down/Maps/HF_2009_stack.rds")

plot(HF_2009_stack)

rm(Built_2009, Croplands_2005, Lights_2009, NavWater_2009, Pasture_2009, PopDensity_2010, Railways, Roads)

# Aggregate to smaller resolution because too big to handle...

continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_15.rds"))
# HF_2009_stack_WGS84 <- projectRaster(from = HF_2009_stack, to = continent_mask, method = "ngb")
HF_2009_stack_WGS84 <- projectRaster(from = HF_2009_stack, res = 0.25, crs = continent_mask@crs, method = "ngb")

HF_2009_stack_Neotropics_15min <- raster::resample(x = HF_2009_stack_WGS84, y = continent_mask, method = "ngb")
plot(HF_2009_stack_Neotropics_15min)

saveRDS(HF_2009_stack_Neotropics_15min, file = "./input_data/HII/HumanFootprint_break_down/Maps/HF_2009_stack_Neotropics_15min.rds")
HF_2009_stack_Neotropics_15min <- readRDS(file = "./input_data/HII/HumanFootprint_break_down/Maps/HF_2009_stack_Neotropics_15min.rds")


## Compute sum of all pressure scores

# Too long to do on the full initial dataset
HF_2009_rebuilt <- calc(x = HF_2009_stack, fun = sum)

# Compute only for the 15min resolution dataset in the Neotropics
HF_2009_rebuilt_Neotropics_15min <- calc(x = HF_2009_stack_Neotropics_15min, fun = sum)
saveRDS(HF_2009_rebuilt_Neotropics_15min, file = "./input_data/HII/HumanFootprint_break_down/Maps/HF_2009_rebuilt_Neotropics_15min.rds")
HF_2009_rebuilt_Neotropics_15min <- readRDS(file = "./input_data/HII/HumanFootprint_break_down/Maps/HF_2009_rebuilt_Neotropics_15min.rds")

plot(HF_2009_rebuilt_Neotropics_15min)
plot(HF_2009_original_Neotropics_15min)

# Correlation between ranks of value from original and rebuilt HF
cor.test(x = HF_2009_original_Neotropics_15min[], y = HF_2009_rebuilt_Neotropics_15min[], method = "spearman", alternative = "two.sided", na.rm = T)


##### 3/ Tests for correlations ####

### 3.1/ Within the dimensions ####

# Add the final HF too, to be able to evaluate which variable relates the most to the final cumulative index

HF_2009_stack_Neotropics_15min_full <- stack(HF_2009_original_Neotropics_15min, HF_2009_stack_Neotropics_15min)

# Reorganize order according to correlation to HF 2009
HF_2009_stack_Neotropics_15min_full <- raster::subset(x = HF_2009_stack_Neotropics_15min_full, subset = c(1, 7, 9, 3, 4, 6, 2, 8, 5))

names(HF_2009_stack_Neotropics_15min_full) <- c("Human Footprint", "Pop. density", "Roads", "Croplands", "Lights", "Pasture", "Built", "Railways", "Nav. Waterways")

plot(HF_2009_stack_Neotropics_15min_full)

Nb_dims <- nlayers(HF_2009_stack_Neotropics_15min_full)

rho_matrix <- rho_uncorrected_pvalues_matrix <- matrix(nrow = Nb_dims, ncol = Nb_dims)
for (i in 1:Nb_dims)
{
  # for (j in (i+1):Nb_dims) # Most efficient version with no useless double tests
  for (j in 1:Nb_dims) # Version to fill the all matrix
  {
    test <- cor.test(x = HF_2009_stack_Neotropics_15min_full[[i]][], y = HF_2009_stack_Neotropics_15min_full[[j]][], method = "spearman", alternative = "two.sided", na.rm = T)
    rho_matrix[i,j] <- round(test$estimate, 3)
    rho_uncorrected_pvalues_matrix[i,j] <- round(test$p.value, 3)
  }
}

rownames(rho_matrix) <- colnames(rho_matrix) <- names(HF_2009_stack_Neotropics_15min_full)
rownames(rho_uncorrected_pvalues_matrix) <- colnames(rho_uncorrected_pvalues_matrix) <- names(HF_2009_stack_Neotropics_15min_full)

rho_matrix
rho_uncorrected_pvalues_matrix


corrplot::corrplot(corr = rho_matrix, method = "circle",  
         type = "upper", # order = "hclust",
         addCoef.col = "black", # Ajout du coefficient de corrélation
         # labels = c("Human Footprint", "Pop. density", "Roads", "Croplands", "Lights", "Pasture", "Built", "Railways", "Nav. Waterways"),
         tl.col = "black", tl.srt = 90, font = 2, # Rotation des etiquettes de textes
         # Combiner avec le niveau de significativité
         # p.mat = rho_uncorrected_pvalues_matrix, sig.level = 0.05, insig = "blank", 
         number.cex = 0.8, number.font = 2, outline = F,
         # Cacher les coefficients de corrélation sur la diagonale
         diag = FALSE)




### 3.2/ Test correlation with and without Light Pollution ####

HF_2009_no_light_rebuilt_Neotropics_15min <- calc(x = raster::subset(x = HF_2009_stack_Neotropics_15min, subset = which(!(names(HF_2009_stack_Neotropics_15min) == "Lights2009"))), fun = sum)

plot(HF_2009_no_light_rebuilt_Neotropics_15min*2)
plot(HF_2009_rebuilt_Neotropics_15min*2)
plot(HF_2009_original_Neotropics_15min*2)

# Correlation between ranks of value from rebuilt HF with and without Light pollution
cor.test(x = HF_2009_rebuilt_Neotropics_15min[], y = HF_2009_no_light_rebuilt_Neotropics_15min[], method = "spearman", alternative = "two.sided", na.rm = T)

# Correlation between ranks of value from original HF and rebuilt HF without Light pollution
cor.test(x = HF_2009_original_Neotropics_15min[], y = HF_2009_no_light_rebuilt_Neotropics_15min[], method = "spearman", alternative = "two.sided", na.rm = T)


### 3.2/ Test correlation with and without Waterways ####

HF_2009_no_waterways_rebuilt_Neotropics_15min <- calc(x = raster::subset(x = HF_2009_stack_Neotropics_15min, subset = which(!(names(HF_2009_stack_Neotropics_15min) == "NavWater2009"))), fun = sum)

plot(HF_2009_no_waterways_rebuilt_Neotropics_15min*2)
plot(HF_2009_rebuilt_Neotropics_15min*2)

# Correlation between ranks of value from rebuilt HF with and without Light pollution
cor.test(x = HF_2009_rebuilt_Neotropics_15min[], y = HF_2009_no_waterways_rebuilt_Neotropics_15min[], method = "spearman", alternative = "two.sided", na.rm = T)

# Correlation between ranks of value from original HF and rebuilt HF without Light pollution
cor.test(x = HF_2009_original_Neotropics_15min[], y = HF_2009_no_waterways_rebuilt_Neotropics_15min[], method = "spearman", alternative = "two.sided", na.rm = T)


### 3.3/ Test correlation with and without Light pollution & Waterways ####

HF_2009_no_light_no_waterways_rebuilt_Neotropics_15min <- calc(x = raster::subset(x = HF_2009_stack_Neotropics_15min, subset = which(!(names(HF_2009_stack_Neotropics_15min) == c("Lights2009", "NavWater2009")))), fun = sum)

plot(HF_2009_no_light_no_waterways_rebuilt_Neotropics_15min*2)
plot(HF_2009_rebuilt_Neotropics_15min*2)

# Correlation between ranks of value from rebuilt HF with and without Light pollution
cor.test(x = HF_2009_rebuilt_Neotropics_15min[], y = HF_2009_no_light_no_waterways_rebuilt_Neotropics_15min[], method = "spearman", alternative = "two.sided", na.rm = T)

# Correlation between ranks of value from original HF and rebuilt HF without Light pollution
cor.test(x = HF_2009_original_Neotropics_15min[], y = HF_2009_no_light_no_waterways_rebuilt_Neotropics_15min[], method = "spearman", alternative = "two.sided", na.rm = T)



### 4/ Plot HF with and without light pollution and waterways ####

### 4.1/ Rescale to 0-100 ####

# Maximum for HF with light pollution = 50. Scaled on 100.
HF_full <- HF_2009_rebuilt_Neotropics_15min/50*100

# Maximum for HF without light pollution = 40. Scaled on 100
HF_no_light <- HF_2009_no_light_rebuilt_Neotropics_15min/40*100

# Maximum for HF without waterways = 46. Scaled on 100
HF_no_waterways <- HF_2009_no_waterways_rebuilt_Neotropics_15min/46*100

# Maximum for HF without light pollution and waterways = 36. Scaled on 100
HF_no_light_no_waterways <- HF_2009_no_light_no_waterways_rebuilt_Neotropics_15min/36*100

save(HF_full, HF_no_light, HF_no_waterways, HF_no_light_no_waterways, file = "input_data/HII/HumanFootprint_break_down/Maps/HF_comparison_layers_15min.RData")
load(file = "input_data/HII/HumanFootprint_break_down/Maps/HF_comparison_layers_15min.RData")


### 4.2/ Correlogram ###

HF_comparison <- stack(HF_full, HF_no_light, HF_no_waterways, HF_no_light_no_waterways)
names(HF_comparison) <- c("HF full", "HF no Lights", "HF no Waterways", "HF no Lights no Waterways")
Nb_dims <- nlayers(HF_comparison)

rho_matrix <- rho_uncorrected_pvalues_matrix <- matrix(nrow = Nb_dims, ncol = Nb_dims)
for (i in 1:Nb_dims)
{
  # for (j in (i+1):Nb_dims) # Most efficient version with no useless double tests
  for (j in 1:Nb_dims) # Version to fill the all matrix
  {
    test <- cor.test(x = HF_comparison[[i]][], y = HF_comparison[[j]][], method = "spearman", alternative = "two.sided", na.rm = T)
    rho_matrix[i,j] <- round(test$estimate, 3)
    rho_uncorrected_pvalues_matrix[i,j] <- round(test$p.value, 3)
  }
}

rownames(rho_matrix) <- colnames(rho_matrix) <- names(HF_comparison)
rownames(rho_uncorrected_pvalues_matrix) <- colnames(rho_uncorrected_pvalues_matrix) <- names(HF_comparison)

rho_matrix
rho_uncorrected_pvalues_matrix

corrplot::corrplot(corr = rho_matrix, method = "circle",  
                   type = "upper", # order = "hclust",
                   addCoef.col = "black", # Ajout du coefficient de corrélation
                   # labels = c("Human Footprint", "Pop. density", "Roads", "Croplands", "Lights", "Pasture", "Built", "Railways", "Nav. Waterways"),
                   tl.col = "black", tl.srt = 90, font = 2, # Rotation des etiquettes de textes
                   # Combiner avec le niveau de significativité
                   # p.mat = rho_uncorrected_pvalues_matrix, sig.level = 0.05, insig = "blank", 
                   number.cex = 1, number.font = 2, number.digits = 3, outline = F,
                   # Cacher les coefficients de corrélation sur la diagonale
                   diag = FALSE)

# ### 4.3/ Build quantiles shapes ####
# 
# ### Add shapes of top 5% and 25% quantiles !
# 
# # Make a function
# 
# raster <- HF_full
# quantiles_break <- c(0.75, 0.95)
# 
# build_top_ranking_shp <- function(raster, quantiles_break)
# {
#   quantiles_values <- quantile(x = raster, probs = quantiles_break, na.rm = TRUE)
#   
#   for (i in 1:length(quantiles_values))
#   {
#     # i <- 2
#     
#     quantile <- quantiles_values[i]
#     mask <- raster >= quantile
#     
#     shp <- raster::rasterToPolygons(x = mask, dissolve = T) # To merge contiguous polygons of the same field value
#     
#     plot(raster)
#     plot(mask)
#     plot(shp, add = T)
#     
#     
#   }
#   
#   return()
# }

### 4.4/ Load the stuff for nice plots ####

# res <- "5"
res <- "15"

library(tmaptools)
pal_grn_red_NA <- pal_grn_red <- rev(get_brewer_pal("RdYlGn", n = 200, plot = F))
pal_grn_red_NA[1] <- "#EDEDED"

# Load mask for continent borders
continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_", res, ".rds"))
crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".rds"))

# Load shp for plot border, and grid
grid_Mollweide_out <- readRDS(file = "./input_data/Map_stuff/grid_Mollweide_out.rds")
large_bg_mask_Mollweide <- readRDS(file = "./input_data/Map_stuff/large_bg_mask_Mollweide.rds")
bbox_sp_Mollweide <- readRDS(file = "./input_data/Map_stuff/bbox_sp_Mollweide.rds")

load(file = "./input_data/Map_stuff/country_borders.RData")
country_borders <- as(country_borders, "Spatial")


### 4.5/ Contrast maps ####

# Contrasting raster function
contrasting_raster_15min <- function(x, zmin, zmax)
{
  # ithomiini_range_15min <- readRDS(file = paste0("./input_data/Map_stuff/ithomiini_range_15min.rds"))
  continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_15.rds"))
  
  x[x[] <= zmin] <- zmin  # Fix low values
  x[x[] >= zmax] <- zmax  # Fix high values
  
  # x <- mask(x, mask = ithomiini_range_15min)  # Cut out values that are outside Ithomiini range
  
  y <- continent_mask + zmin  # Create final new raster from continental mask with baseline = zmin
  y[!is.na(x[])] <- x[!is.na(x[])]  # Add initial raster values
  
  return(y)
}

hist(HF_full)

HF_full_contrasted <- contrasting_raster_15min(x = HF_full, zmin = 0, zmax = 30)

plot(HF_full, col = pal_grn_red)
plot(HF_full_contrasted, col = pal_grn_red)

HF_no_light_contrasted <- contrasting_raster_15min(x = HF_no_light, zmin = 0, zmax = 30)
HF_no_waterways_contrasted <- contrasting_raster_15min(x = HF_no_waterways, zmin = 0, zmax = 30)
HF_no_light_no_waterways_contrasted <- contrasting_raster_15min(x = HF_no_light_no_waterways, zmin = 0, zmax = 30)


### 4.6/ Project into Mollweide ####

# 4.6.1/ Projection functions

# Project raster into Mollweide
Mollweide_projection <- function(x, name) # Raster to project
{
  x_name <- name # Get the name of the initial raster as a character string
  
  new_map <- projectRaster(from = x, 
                           method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                           crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                           alignOnly = F)
  
  # Generate new object with "_Mollweide" suffix in the global environment
  eval(call("<<-", as.name(paste0(x_name, "_Mollweide")), new_map))
}

# Project sp shape into Mollweide
Mollweide_shp_projection <-  function(x) # Shp to project
{
  x_name <- deparse(substitute(x)) # Get the name of the initial shp as a character string
  
  new_shp <- spTransform(x, CRSobj = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
  
  # Generate new object with "_Mollweide" suffix in the global environment
  eval(call("<<-", as.name(paste0(x_name, "_Mollweide")), new_shp))
}

# 4.6.2/ Project all maps

# Put all maps in a list

list_all_maps <- list(HF_full_contrasted, HF_no_light_contrasted, HF_no_waterways_contrasted, HF_no_light_no_waterways_contrasted)
names(list_all_maps) <- c("HF_full_contrasted", "HF_no_light_contrasted", "HF_no_waterways_contrasted", "HF_no_light_no_waterways_contrasted")

# Loop to generate projected maps
for (i in 1:length(list_all_maps))
{
  Mollweide_projection(x = list_all_maps[[i]], name = names(list_all_maps)[i])
}

# Project country_borders shp
Mollweide_shp_projection(country_borders)

#### 4.7/ Plot HF maps ####

# Function to map indices
{
  map_indices_Mollweide <- function(x,                                    # Raster to map
                                    color_palette = pal_bl_red_Mannion,   # Color palette
                                    main_title,                           # Main title
                                    main_title_cex = 1.5,                 # Main title size
                                    
                                    xlim = c(-4600, 4600),   # Limit of plot on x-axis (Longitude)
                                    ylim = c(-4450, 3400),    # Limit of plot on y-axis (Latitude)
                                    axis_cex = 1.4,             # Axes size
                                    
                                    xlab = "",                # X-axis label
                                    ylab = "",                # Y-axis label
                                    x_axis_breaks = c(-3930, -2170, -420, 1500, 3050),            # X-axis tick breaks
                                    y_axis_breaks = c(-3650, -2450, -1220, 0, 1230, 2445),        # Y-axis tick breaks
                                    x_axis_labels = c("120°E", "100°E", "80°E", "60°E", "40°E"),      # X-axis tick labels
                                    y_axis_labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"),  # Y-axis tick labels
                                    
                                    legend_title,             # Legend title
                                    legend_title_cex = 1.4,   # Legend title size
                                    legend_title_x = -3550,   # Legend title x position
                                    legend_title_y = 430,     # Legend title y position
                                    legend_cex = 1.4,         # Legend size
                                    legend_breaks,            # Legend tick positions
                                    legend_location = c(-4100, -3800, -3950, 0),  # Legend position
                                    
                                    scale_bar_position = c(-2600, -4000),  # Scale bar position
                                    
                                    arrow_scale = 0.55,           # North arrow size
                                    arrow_padin = c(0.15, 0.15),  # North arrow position adjustement
                                    
                                    facet_letter = "",                  # Small case letter for facet
                                    facet_letter_col = "black",         # Color of case letter for facet
                                    facet_letter_cex = 2.2,             # Size of small case letter for facet
                                    facet_letter_inset = c(0, -0.008))  # Position adjustment of small case letter for facet
  
  {
    # Plot raster background without axis
    image(x, col = color_palette,
          xlim = xlim, ylim = ylim, axes = F,
          xlab = xlab, ylab = ylab)
    title(main = main_title, cex.main = main_title_cex, line = 1)
    
    # Generate axes with manual positioning of ticks
    axis(side = 1, at = x_axis_breaks, labels = x_axis_labels, cex.axis = axis_cex, lwd = 0.2, lwd.ticks = 1, gap.axis = 0, padj = 0.5)
    axis(side = 2, at = y_axis_breaks, labels = y_axis_labels, cex.axis = axis_cex, lwd = 0.2, lwd.ticks = 1, gap.axis = 0)
    
    # Add background, borders and graticules
    plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
    plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
    plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
    plot(country_borders_Mollweide, lwd = 1, border = "#00000030", col = NA, add = T)
    
    # Add scale bar in legend
    scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = scale_bar_position, label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1.2)
    prettymapr::addnortharrow(scale = arrow_scale, padin = arrow_padin, text.col = "#00000000")
    rangeBuilder::addRasterLegend(x, locs = legend_breaks, cex.axis = legend_cex, ramp = color_palette, ncolors = 200, border = T, location = legend_location)
    rangeBuilder::addRasterLegend(x, locs = legend_breaks, cex.axis = legend_cex, ramp = color_palette, ncolors = 200, border = T, location = legend_location)
    graphics::text(x = legend_title_x, y = legend_title_y, font = 2, cex = legend_title_cex, label = legend_title)
    
    # Add facet letter
    legend(legend = facet_letter, x = "bottomright", bty = "n", text.col = facet_letter_col,
           text.font = 2, cex = facet_letter_cex, inset = facet_letter_inset)
    
  }
}

### Make a 4 facets plot

pdf(file = paste0("./supplementaries/HF_comparison_maps.pdf"), height = 10, width = 10)

internal_margins <- par()$mar
par(mar = c(3.1, 3.1, 2.7, 1.6))
par(mfrow = c(2,2))

# Panel A = Human Footprint full
map_indices_Mollweide(x = HF_full_contrasted_Mollweide,
                      main_title = "Human Footprint (HF)",
                      color_palette = pal_grn_red,
                      axis_cex = 1.3,
                      legend_title = "Human\nFootprint",
                      legend_title_x = -3350,
                      legend_title_y = 880,
                      legend_breaks = seq(0, 30, 5))
text(x = -2600, y = 0, labels = "-100", cex = 1.4)  # Add the contrast on the scale
text(x = -2600, y = 0, labels = "-100", cex = 1.4)  # Add the contrast on the scale
text(x = 4000, y = -3900, labels = "(a)", font = 2, cex = 1.8)  # Add the contrast on the scale

# Panel B = Human Footprint without Light Pollution
map_indices_Mollweide(x = HF_no_light_contrasted_Mollweide,
                      main_title = "HF without Light Pollution (LP)",
                      color_palette = pal_grn_red,
                      axis_cex = 1.3,
                      legend_title = "Human\nFootprint",
                      legend_title_x = -3350,
                      legend_title_y = 880,
                      legend_breaks = seq(0, 30, 5))
text(x = -2600, y = 0, labels = "-100", cex = 1.4)  # Add the contrast on the scale
text(x = -2600, y = 0, labels = "-100", cex = 1.4)  # Add the contrast on the scale
text(x = 4000, y = -3900, labels = "(b)", font = 2, cex = 1.8)  # Add the contrast on the scale

# Add correlation with full HF
text(x = 3450, y = 1800, labels = bquote(rho ~ "= 0.997"), font = 2, cex = 1.4)
text(x = 3450, y = 1800, labels = bquote(rho ~ "= 0.997"), font = 2, cex = 1.4)


# Panel C = Human Footprint without Navigable Waterways
map_indices_Mollweide(x = HF_no_waterways_contrasted_Mollweide,
                      main_title = "HF without Navigable Waterways (NW)",
                      color_palette = pal_grn_red,
                      axis_cex = 1.3,
                      legend_title = "Human\nFootprint",
                      legend_title_x = -3350,
                      legend_title_y = 880,
                      legend_breaks = seq(0, 30, 5))
text(x = -2600, y = 0, labels = "-100", cex = 1.4)  # Add the contrast on the scale
text(x = -2600, y = 0, labels = "-100", cex = 1.4)  # Add the contrast on the scale
text(x = 4000, y = -3900, labels = "(c)", font = 2, cex = 1.8)  # Add the contrast on the scale

# Add correlation with full HF
text(x = 3450, y = 1800, labels = bquote(rho ~ "= 0.991"), font = 2, cex = 1.4)
text(x = 3450, y = 1800, labels = bquote(rho ~ "= 0.991"), font = 2, cex = 1.4)


# Panel D = Human Footprint without Light Pollution and Navigable Waterways
map_indices_Mollweide(x = HF_no_light_no_waterways_contrasted_Mollweide,
                      main_title = "HF without LP and NW",
                      color_palette = pal_grn_red,
                      axis_cex = 1.3,
                      legend_title = "Human\nFootprint",
                      legend_title_x = -3350,
                      legend_title_y = 880,
                      legend_breaks = seq(0, 30, 5))
text(x = -2600, y = 0, labels = "-100", cex = 1.4)  # Add the contrast on the scale
text(x = -2600, y = 0, labels = "-100", cex = 1.4)  # Add the contrast on the scale
text(x = 4000, y = -3900, labels = "(d)", font = 2, cex = 1.8)  # Add the contrast on the scale

# Add correlation with full HF
text(x = 3450, y = 1800, labels = bquote(rho ~ "= 0.988"), font = 2, cex = 1.4)
text(x = 3450, y = 1800, labels = bquote(rho ~ "= 0.988"), font = 2, cex = 1.4)


par(mar = internal_margins)
par(mfrow = c(1,1))
dev.off()