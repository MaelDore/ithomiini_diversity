
### Additional maps ###


# Clean environment
rm(list = ls())

### 1/ Load stuff ####

# Packages
library(raster)
library(prettymapr)
library(rangeBuilder)

# New color palette
pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")

# Load mask for continent borders, plot border, and grid
grid_Mollweide_out <- readRDS(file = "./input_data/Map_stuff/grid_Mollweide_out.rds")
large_bg_mask_Mollweide <- readRDS(file = "./input_data/Map_stuff/large_bg_mask_Mollweide.rds")
bbox_sp_Mollweide <- readRDS(file = "./input_data/Map_stuff/bbox_sp_Mollweide.rds")
Community_mask <- readRDS(file = "./input_data/Map_stuff/Community_mask.rds")
continent_mask <- readRDS(file = paste0("./input_data/Map_stuff/continent_mask.rds"))

### Load maps
sp.richness <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
ring.richness <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))
vulnerability <- readRDS(file = paste0("./outputs/Indices_maps/vulnerability_Jaccard.80.rds"))

# Load mimicry ring richness stack
ring_rich_stack <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_Jaccard.80.RData"))

# Load mimicry ring proba stack
ring_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_Jaccard.80.RData"))

### 2/ Generate useful functions ####

# 2.1/ Mollweide projection ####

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

Mollweide_projection(sp.richness)

# 2.2/ Plotting function ####

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

# 2.3/ Add continental null values to raster ####

add_continental_null_values <- function(x)
{
  continent_mask <- readRDS(file = paste0("./input_data/Map_stuff/continent_mask.rds"))
  
  y <- continent_mask  # Create final new raster from continental mask
  y[!is.na(x[])] <- x[!is.na(x[])]  # Add initial raster values
  
  return(y)
}

# 2.4/ Contrasting raster ####
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

### 3/ Compute mean ring size ####

# 3.1/ Compute index
mean.ring.size <- sp.richness/ring.richness
mean.ring.size <- add_continental_null_values(mean.ring.size)
save(mean.ring.size, file = "./outputs/Indices_maps/mean.ring.size.RData")

# 3.2/ Project into Mollweide and save results
Mollweide_projection(mean.ring.size)

# 3.3/ Plot 
map_indices_Mollweide(x = mean.ring.size_Mollweide,
                      main_title = "Mean ring size",
                      legend_title = "Ring size",
                      legend_title_y = 530,
                      legend_breaks = seq(0, 6, 1), 
                      facet_letter = "a")

### 4/ Weighted mean ring size ####

# (weights = proba of presence of the ring)

# 4.1/ Compute index

?weighted.mean
# test <- overlay(x = ring_rich_stack, y = ring_proba_stack, fun = weighted.mean, na.rm = T, unstack=F)

ring_rich_brick <- brick(ring_rich_stack)
ring_proba_brick <- brick(ring_proba_stack)

View(ring_rich_brick[])

weighted_mean_ring_size <- vector()
for (i in 1:nrow(ring_rich_brick[]))  # Loop between communities/pixels
{
  # i <- 55000
  
  rich_vector <- ring_rich_brick@data@values[i,]   # Extract ring richness
  proba_vector <- ring_proba_brick@data@values[i,] # Extract ring proba of presence
  
  # Compute mean weighted by proba. of presence
  weighted_mean_ring_size[i] <- weighted.mean(x = rich_vector, w = proba_vector)
}
  
weighted_mean_ring_size <- add_continental_null_values(weighted_mean_ring_size)
save(weighted_mean_ring_size, file = "./outputs/Indices_maps/weighted_mean_ring_size.RData")
saveRDS(weighted_mean_ring_size, file = "./outputs/Indices_maps/weighted_mean_ring_size.rds")

# 4.2/ Project into Mollweide and save results
Mollweide_projection(weighted_mean_ring_size)

# 4.3/ Plot 
map_indices_Mollweide(x = weighted_mean_ring_size_Mollweide,
                      main_title = "Weighted mean ring size",
                      legend_title = "Ring size",
                      legend_title_y = 530,
                      legend_breaks = seq(0, 6, 1), 
                      facet_letter = "a")

### 5/ Max ring size ####

# 5.1/ Compute index
max_ring_size <- calc(ring_rich_stack, fun = max)*1
save(max_ring_size, file = "./outputs/Indices_maps/max_ring_size.RData")

# 5.2/ Project into Mollweide and save results
Mollweide_projection(max_ring_size)

# 5.3/ Plot 
map_indices_Mollweide(x = max_ring_size_Mollweide,
                      main_title = "Max ring size",
                      legend_title = "Ring size",
                      legend_title_y = 530,
                      legend_breaks = seq(0, 20, 5), 
                      facet_letter = "a")

### 6/ Quantiles 25%/50%/75% of ring size ####

# 6.1/ Compute indices

# Need to use rounded to then select quantiles otherwise, very unlikely present rings will pull down the index in rich regions

ring_rich_brick_rounded <- round(ring_rich_brick, digits = 0)
table(ring_rich_brick_rounded@data@values)

# Not working for some memory issues
# quantile_25_ring_size <- calc(x = ring_rich_stack_rounded, fun = quantile)

quantile_25_ring_size <- quantile_75_ring_size <- median_ring_size <- vector()
for (i in 1:nrow(ring_rich_brick_rounded[]))  # Loop between communities/pixels
{
  # i <- 52400
  
  # Extract rounded ring richness
  rich_vector <- ring_rich_brick_rounded@data@values[i,]   
  # Remove null values
  rich_vector_positive <- rich_vector[!rich_vector == 0] 
  
  # Extract quantiles 25%, median, 75%
  quantiles_i <- quantile(x = rich_vector_positive, probs = c(0.25, 0.5, 0.75), na.rm = T)
  
  # Store quantiles into vectors
  quantile_25_ring_size[i] <- quantiles_i[1]
  median_ring_size[i] <- quantiles_i[2]
  quantile_75_ring_size[i] <- quantiles_i[3]

}
# Add continental null values and put into a raster layer
quantile_25_ring_size <- add_continental_null_values(quantile_25_ring_size)
median_ring_size <- add_continental_null_values(median_ring_size)
quantile_75_ring_size <- add_continental_null_values(quantile_75_ring_size)

# save
save(quantile_25_ring_size, file = "./outputs/Indices_maps/quantile_25_ring_size.RData")
save(median_ring_size, file = "./outputs/Indices_maps/median_ring_size.RData")
save(quantile_75_ring_size, file = "./outputs/Indices_maps/quantile_75_ring_size.RData")

# 6.2/ Project into Mollweide and save results
Mollweide_projection(quantile_25_ring_size)
Mollweide_projection(median_ring_size)
Mollweide_projection(quantile_75_ring_size)

plot(median_ring_size)

# 6.3/ Plot
map_indices_Mollweide(x = quantile_25_ring_size_Mollweide,
                      main_title = "Quantile 25% ring size",
                      legend_title = "Ring size",
                      legend_title_y = 530,
                      legend_breaks = seq(0, 4, 1), 
                      facet_letter = "a")

map_indices_Mollweide(x = quantile_75_ring_size_Mollweide,
                      main_title = "Quantile 75% ring size",
                      legend_title = "Ring size",
                      legend_title_y = 530,
                      legend_breaks = seq(0, 10, 2), 
                      facet_letter = "a")

map_indices_Mollweide(x = median_ring_size_Mollweide,
                      main_title = "Median ring size",
                      legend_title = "Ring size",
                      legend_title_y = 530,
                      legend_breaks = seq(0, 8, 2), 
                      facet_letter = "a")

### 7/ Mean ring phenotypic rarity and prevalence ####

# Normalized_function

normalized_0_1 <- function(x, ...)
{
  if(max(x) == min(x)){
    y <- x/length(x)
  } else {
    y <- (x-min(x))/(max(x)-min(x))
  }
  return(y)
}


# 7.1/ Compute ring phenotypic rarity based on counts of OMUs per ring

# Extract order of ring layers in stack
ring_ordered <- data.frame(names(ring_proba_stack)) ; names(ring_ordered) <- "ring"

load(file = paste0("./input_data/list.models.RData"))

N_OMUs_per_ring <- list.models %>% 
  group_by(Mimicry.model) %>% 
  summarize(counts = n()) %>% 
  ungroup() %>% 
  mutate(normed_commonness = normalized_0_1(counts)) %>% 
  mutate(phenotypic_rarity = 1 - normed_commonness) %>% 
  left_join(x = ring_ordered, y = ., by = c("ring" = "Mimicry.model"))

save(N_OMUs_per_ring, file = paste0("./outputs/N_OMUs_per_ring.RData"))

hist(N_OMUs_per_ring$counts)

# Ideally, ring phenotypic rarity should be evaluated in the phenospace weighted by occurrences in the phylogeny. 
# Also phylogeny is not accounted. A way to account for phylogeny would be to count independent appearance of the patern instead of all occurrences

# 7.2/ Compute indices

# Weighted mean of phenotypic rarity

# Weights = probability of presence of the ring
# Values for rarity based on size of the ring in the whole tribe => 1 - normalized counts of OMUs per ring

mean_ring_phenotypic_rarity <- calc(x = ring_proba_stack, 
                                    fun = function(x) {weighted.mean(x = N_OMUs_per_ring$phenotypic_rarity, w = x)})*1

# Weighted mean of ring prevalence (size in the whole tribe)
mean_ring_prevalence <- calc(x = ring_proba_stack, 
                                 fun = function(x) {weighted.mean(x = N_OMUs_per_ring$counts, w = x)})*1

# Add null continental values
mean_ring_phenotypic_rarity <- add_continental_null_values(mean_ring_phenotypic_rarity)
mean_ring_prevalence <- add_continental_null_values(mean_ring_prevalence)

# Save
save(mean_ring_phenotypic_rarity, file = "./outputs/Indices_maps/mean_ring_phenotypic_rarity.RData")
save(mean_ring_prevalence, file = "./outputs/Indices_maps/mean_ring_prevalence.RData")

hist(mean_ring_phenotypic_rarity)
hist(mean_ring_prevalence)

# 7.3/ Project into Mollweide and save results
Mollweide_projection(mean_ring_phenotypic_rarity)
Mollweide_projection(mean_ring_prevalence)

# 7.4/ Plot
map_indices_Mollweide(x = mean_ring_phenotypic_rarity_Mollweide,
                      main_title = "Mean ring phenotypic rarity",
                      legend_title = "Phenotypic\nrarity",
                      legend_title_y = 730,
                      legend_breaks = seq(0, 0.8, 0.2), 
                      facet_letter = "a")

map_indices_Mollweide(x = mean_ring_prevalence_Mollweide,
                      main_title = "Mean ring prevalence",
                      legend_title = "Ring size",
                      legend_title_y = 530,
                      legend_breaks = seq(0, 80, 20), 
                      facet_letter = "a")

# 7.5/ Contrasted plot

# Weighted mean phenotypic rarity
hist(mean_ring_phenotypic_rarity) # 0.4 to 0.8

mean_ring_phenotypic_rarity_contrasted <- contrasting_raster(x = mean_ring_phenotypic_rarity, zmin = 0.4, zmax = 0.8)
Mollweide_projection(mean_ring_phenotypic_rarity_contrasted)

map_indices_Mollweide(x = mean_ring_phenotypic_rarity_contrasted_Mollweide,
                      main_title = "Mean ring phenotypic rarity",
                      legend_title = "Phenotypic\nrarity",
                      legend_title_y = 730,
                      legend_breaks = seq(0.4, 0.8, 0.1), 
                      facet_letter = "a")

# Mean ring prevalence
hist(mean_ring_prevalence) # 25 to 55

mean_ring_prevalence_contrasted <- contrasting_raster(x = mean_ring_prevalence, zmin = 25, zmax = 55)
Mollweide_projection(mean_ring_prevalence_contrasted)

map_indices_Mollweide(x = mean_ring_prevalence_contrasted_Mollweide,
                      main_title = "Mean ring prevalence",
                      legend_title = "Ring size",
                      legend_title_y = 530,
                      legend_breaks = seq(25, 55, 5), 
                      facet_letter = "a")

### 8/ OMUs mean phenotypic rarity ####

# 8.1/ Load stuff

# Load stack of OMU proba
OMU_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/All_OMU_stack_Jaccard.80.rds"))

# Load Summary table for OMUs
load(file = paste0("./input_data/list.models.RData"))
# Load summary table for rings
load(file = paste0("./outputs/N_OMUs_per_ring.RData"))

list.models <- list.models %>% 
  left_join(x = ., y = N_OMUs_per_ring[, c("ring", "phenotypic_rarity")], by = c("Mimicry.model" = "ring"))

# 8.2/ Compute indices

# Weighted mean of phenotypic rarity

# Weights = probability of presence of the OMUs
# Values for rarity based on size of the ring in the whole tribe => 1 - normalized counts of OMUs per ring

identical(list.models$Tag.model, names(OMU_proba_stack))

mean_OMU_phenotypic_rarity <- calc(x = OMU_proba_stack, 
                                    fun = function(x) {weighted.mean(x = list.models$phenotypic_rarity, w = x)})*1

# Add null continental values
mean_OMU_phenotypic_rarity <- add_continental_null_values(mean_OMU_phenotypic_rarity)

# Save
save(mean_OMU_phenotypic_rarity, file = "./outputs/Indices_maps/mean_OMU_phenotypic_rarity.RData")

# 8.3/ Project into Mollweide and save results
Mollweide_projection(mean_OMU_phenotypic_rarity)

# 8.4/ Plot
map_indices_Mollweide(x = mean_OMU_phenotypic_rarity_Mollweide,
                      main_title = "Mean OMU phenotypic rarity",
                      legend_title = "Phenotypic\nrarity",
                      legend_title_y = 730,
                      legend_breaks = seq(0, 0.8, 0.2), 
                      facet_letter = "a")

# 8.5/ Contrasted plot

# Weighted mean phenotypic rarity
hist(mean_OMU_phenotypic_rarity) # 0.4 to 0.7

mean_OMU_phenotypic_rarity_contrasted <- contrasting_raster(x = mean_OMU_phenotypic_rarity, zmin = 0.4, zmax = 0.7)
Mollweide_projection(mean_OMU_phenotypic_rarity_contrasted)

map_indices_Mollweide(x = mean_OMU_phenotypic_rarity_contrasted_Mollweide,
                      main_title = "Mean OMU phenotypic rarity",
                      legend_title = "Phenotypic\nrarity",
                      legend_title_y = 730,
                      legend_breaks = seq(0.4, 0.7, 0.1), 
                      facet_letter = "a")


### 9/ Community vulnerability ####

# 9.1/ Compute unweighted version (assuming all ring with N > 0.5 are present)

# Function to compute community vulnerability as the sum of rings vulnerability approximated as the inverse of the local richness of each mimicry ring
# A ring with only one species as a vulnerability of 1/1 = 1. A ring of 4 species as a vulnerability of 1/4 = 0.25.
# Final vulnerability is standardized by total ring richness such as a community mean vulnerability

vulnerability = function(x, na.rm) {
  x <- round(x, digits = 0) # Need to round ring richness to avoid inflation of value due to numerous rings with richness values close to 0 but not null
  if (sum(x, na.rm = T) > 0) { # Computed only if local mimicry richness >= 1 once rounded
    x <- x[x>0] # Remove all 0 values to avoid error with 1/0
    V <- sum(1/x, na.rm = T) # Compute non-standardized vulnerability = sum of mimicry ring vulnerability.
    V <- V/length(x) # Standardization by local mimicry ring richness
  }else{
    V <- NA
  }
  return(V) # Output
}

vulnerability <- calc(x = ring_rich_stack, fun = vulnerability)*1

vulnerability <- add_continental_null_values(vulnerability)
save(vulnerability, file = "./outputs/Indices_maps/vulnerability.RData")

# 9.2/ Compute weighted version (weighted by ring probability of presence)

ring_proba_brick <- brick(ring_proba_stack)
ring_rich_brick <- brick(ring_rich_stack)
ring_rich_brick_rounded <- round(ring_rich_brick, digits = 0)
table(ring_rich_brick_rounded@data@values)

weighted_vulnerability <- vector()
for (i in 1:nrow(ring_rich_brick_rounded[]))  # Loop between communities/pixels
{
  # i <- 25
  
  rich_rounded_vector <- ring_rich_brick_rounded@data@values[i,]   # Extract ring richness
  if (sum(rich_rounded_vector, na.rm = T) > 0) { # Computed only if local mimicry richness >= 1 once rounded
    rich_cleaned_vector <- rich_rounded_vector[rich_rounded_vector > 0] # Remove all 0 values to avoid error with 1/0
    ring_vulnerability <- 1/rich_cleaned_vector # Compute non-standardized vulnerabilities = inverse of mimicry ring richness
    
    proba_vector <- ring_proba_brick@data@values[i,] # Extract ring proba of presence
    proba_vector <- proba_vector[rich_rounded_vector > 0] # Remove all ring with richness = 0
    
    # Compute mean weighted by proba. of presence
    weighted_vulnerability[i] <- weighted.mean(x = ring_vulnerability, w = proba_vector)
    
  }else{
    weighted_vulnerability[i] <- NA
  }

}

weighted_vulnerability <- add_continental_null_values(weighted_vulnerability)
save(weighted_vulnerability, file = "./outputs/Indices_maps/weighted_vulnerability.RData")

# 9.3/ Project into Mollweide and save results
Mollweide_projection(vulnerability)
Mollweide_projection(weighted_vulnerability)

# 9.4/ Plot
map_indices_Mollweide(x = vulnerability_Mollweide,
                      main_title = "Community vulnerability",
                      legend_title = "Vulnerability",
                      legend_title_y = 530,
                      legend_breaks = seq(0, 1, 0.2), 
                      facet_letter = "a")

map_indices_Mollweide(x = weighted_vulnerability_Mollweide,
                      main_title = "Community weighted vulnerability",
                      legend_title = "Vulnerability",
                      legend_title_y = 530,
                      legend_breaks = seq(0, 1, 0.2), 
                      facet_letter = "a")


# 9.5/ Contrasted plot

# Weighted mean phenotypic rarity
hist(vulnerability) # 0.2 to 1.0

vulnerability_contrasted <- contrasting_raster(x = vulnerability, zmin = 0.2, zmax = 1.0)
Mollweide_projection(vulnerability_contrasted)

map_indices_Mollweide(x = vulnerability_contrasted_Mollweide,
                      main_title = "Community vulnerability",
                      legend_title = "Vulnerability",
                      legend_title_y = 530,
                      legend_breaks = seq(0.2, 1, 0.2), 
                      facet_letter = "a")

hist(weighted_vulnerability) # 0.2 to 1.0

weighted_vulnerability_contrasted <- contrasting_raster(x = weighted_vulnerability, zmin = 0.2, zmax = 1.0)
Mollweide_projection(weighted_vulnerability_contrasted)

map_indices_Mollweide(x = weighted_vulnerability_contrasted_Mollweide,
                      main_title = "Community weighted vulnerability",
                      legend_title = "Vulnerability",
                      legend_title_y = 530,
                      legend_breaks = seq(0.2, 1, 0.2), 
                      facet_letter = "a")

# 10/ Plot everything all together ####

pdf(file = paste0("./maps/Indices_maps/Index_maps_Mollweide_bonus.pdf"), height = 15.5, width = 12)

internal_margins <- par()$mar
par(mar = c(3.1, 3.5, 2.7, 2.1))
par(mfrow = c(4, 3))

map_indices_Mollweide(x = mean.ring.size_Mollweide,
                      main_title = "Mean ring size",
                      legend_title = "Ring size",
                      legend_title_y = 530,
                      legend_breaks = seq(0, 6, 1), 
                      facet_letter = "a")

map_indices_Mollweide(x = weighted_mean_ring_size_Mollweide,
                      main_title = "Weighted mean ring size",
                      legend_title = "Ring size",
                      legend_title_y = 530,
                      legend_breaks = seq(0, 6, 1), 
                      facet_letter = "b")

map_indices_Mollweide(x = max_ring_size_Mollweide,
                      main_title = "Max ring size",
                      legend_title = "Ring size",
                      legend_title_y = 530,
                      legend_breaks = seq(0, 20, 5), 
                      facet_letter = "c")

map_indices_Mollweide(x = quantile_25_ring_size_Mollweide,
                      main_title = "Quantile 25% ring size",
                      legend_title = "Ring size",
                      legend_title_y = 530,
                      legend_breaks = seq(0, 4, 1), 
                      facet_letter = "d")

map_indices_Mollweide(x = quantile_75_ring_size_Mollweide,
                      main_title = "Quantile 75% ring size",
                      legend_title = "Ring size",
                      legend_title_y = 530,
                      legend_breaks = seq(0, 10, 2), 
                      facet_letter = "e")

map_indices_Mollweide(x = median_ring_size_Mollweide,
                      main_title = "Median ring size",
                      legend_title = "Ring size",
                      legend_title_y = 530,
                      legend_breaks = seq(0, 8, 2), 
                      facet_letter = "f")


map_indices_Mollweide(x = mean_ring_phenotypic_rarity_contrasted_Mollweide,
                      main_title = "Mean ring phenotypic rarity",
                      legend_title = "Phenotypic\nrarity",
                      legend_title_y = 730,
                      legend_breaks = seq(0.4, 0.8, 0.1), 
                      facet_letter = "g")

map_indices_Mollweide(x = mean_ring_prevalence_contrasted_Mollweide,
                      main_title = "Mean ring prevalence",
                      legend_title = "Ring size",
                      legend_title_y = 530,
                      legend_breaks = seq(25, 55, 5), 
                      facet_letter = "h")

map_indices_Mollweide(x = mean_OMU_phenotypic_rarity_contrasted_Mollweide,
                      main_title = "Mean OMU phenotypic rarity",
                      legend_title = "Phenotypic\nrarity",
                      legend_title_y = 730,
                      legend_breaks = seq(0.4, 0.7, 0.1), 
                      facet_letter = "i")

map_indices_Mollweide(x = vulnerability_contrasted_Mollweide,
                      main_title = "Community vulnerability",
                      legend_title = "Vulnerability",
                      legend_title_y = 530,
                      legend_breaks = seq(0.2, 1, 0.2), 
                      facet_letter = "j")

map_indices_Mollweide(x = weighted_vulnerability_contrasted_Mollweide,
                      main_title = "Community weighted vulnerability",
                      legend_title = "Vulnerability",
                      legend_title_y = 530,
                      legend_breaks = seq(0.2, 1, 0.2), 
                      facet_letter = "k")

par(mar = internal_margins)
par(mfrow = c(1, 1))

dev.off()




