
##### Address concerns regarding Sampling effort #####


# Clean environment
rm(list = ls())

##### 1/ Map sampling sites in environmental space #####

# Load dataset
load(file = "./input_data/Ithomiini_final.RData")

# Load environmental stack
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_15.rds"))
env_data_df <- envData@data@values[complete.cases(envData@data@values), ]
dim(env_data_df)  # 24 986 environmental data = communities in the study area

occurrences_env <- raster::extract(x = envData, y = Ithomiini_final[, c("Longitude", "Latitude")])
occurrences_env <- occurrences_env[complete.cases(occurrences_env), ]
dim(occurrences_env)  # 27 795 occurrences (with environmental data available)

sampling_sites <- unique(Ithomiini_final[, c("Longitude_raster", "Latitude_raster")]) # 1834 grid-cell sampling sites

sites_env <- raster::extract(x = envData, y = sampling_sites)
sites_env <- sites_env[complete.cases(sites_env), ]
dim(sites_env)  # 1718 gird-cell with occurrences (with environmental data available)

full_data <- rbind(env_data_df, sites_env)
dim(full_data)

# Apply PCA on environmental variables to reduce into 2D
library("ade4")
library("factoextra")
library("Rcpp")

PCA_res <- FactoMineR::PCA(X = full_data,            # Background data of all environmental conditions available + occurrence data
                           ind.sup = (dim(env_data_df)[1] + 1):(dim(env_data_df)[1] + dim(occurrences_env)[1]),  # Indices of sample sites/occurrence data
                           graph = F)

eigenvalues_summary <- get_eigenvalue(PCA_res) ; eigenvalues_summary # To get EV ordered, % of explained variance and cumulative % of explained var
PC1_var <- round(eigenvalues_summary[1,2], 1)
PC2_var <- round(eigenvalues_summary[2,2], 1)
PCA_summary_ind <- get_pca_ind(PCA_res) # To extract results for individuals under a clearer format. Include cos2 and contrib
PCA_summary_var <- get_pca_var(PCA_res) # To extract results for variables under a clearer format. Include cos2 and contrib

# Extract coordinates
env_coords <- as.data.frame(PCA_res$ind$coord[, 1:2])
occ_coords <- as.data.frame(PCA_res$ind.sup$coord[, 1:2])
names(occ_coords) <- names(env_coords) <- c("PC1", "PC2")

# plot(env_coords, pch = 15, col = "#00000030")
plot(env_coords, pch = 15, col = "grey")
points(occ_coords, pch = 15, col = "red")


### ggplot style

pal_bl_red_Mannion <- readRDS("D:/Mael/R_projects/ithomiini_diversity/maps/pal_bl_red_Mannion.rds")

?geom_bin2d
?scale_fill_continuous

tmaptools::palette_explorer()

range(env_coords$PC1)
range(env_coords$PC2)

## In one plot with alpha transparency for the sample sites 

ggplot(data = env_coords, mapping = aes(x = PC1, y = PC2)) +
  geom_bin2d(bins = 50, binwidth = c(0.2, 0.16)) +
  
  xlim(-7, 3) +
  ylim(-5, 3) +
  xlab(paste0("PC1 (",PC1_var,"%)")) + 
  ylab(paste0("PC2 (",PC2_var,"%)")) +
  
  scale_fill_gradientn(colours = grey.colors(n = 200, start = 0.8, end = 0.2), na.value = NA) +
  labs(fill = "Env.       ") + # Need space because order of legend is determined by number of character in the title...
  
  ggnewscale::new_scale_fill() +
  geom_bin2d(data = occ_coords, aes(x = PC1, y = PC2),
             bins = 50, binwidth = c(0.2, 0.16), alpha = 0.7) +
  
  scale_fill_gradientn(colours = rev(tmaptools::get_brewer_pal("Spectral", n = 100, plot = F)), na.value = pal_bl_red_Mannion[1]) +
  # scale_fill_gradientn(colours = pal_bl_red_Mannion[150:200], na.value = pal_bl_red_Mannion[1]) +
  labs(fill = "Sampling") +
  
  scale_y_continuous(minor_breaks = seq(-7.05 , 3.95, 0.16)) +
  scale_x_continuous(minor_breaks = seq(-7 , 4, 0.2)) +
  labs(title = "Sampling effort in environmental space") +
  
  ggplot2::theme(panel.background = ggplot2::element_rect(fill = NA),
                 plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5, vjust = 3),
                 # legend.position = c(0.85, 0.55),
                 legend.title = ggplot2::element_text(size = 12, vjust = 3, face = "bold"),
                 legend.text = ggplot2::element_text(size = 10, face = "plain"),
                 panel.grid.major = element_line(colour = "grey50", size = 0.5),
                 panel.grid.minor = element_line(colour = "grey90", size = 0.5),
                 plot.margin = ggplot2::margin(t = 18, unit = "pt"),
                 # axis.line.x = ggplot2::element_line(color = NA),
                 axis.line = ggplot2::element_line(color = "black", size = 1.3),
                 axis.ticks = ggplot2::element_line(color = "black", size = 1.2),
                 axis.ticks.length = ggplot2::unit(8, "pt"),
                 axis.text = ggplot2::element_text(size = 14, face = "bold", color = "black"),
                 axis.text.x = ggplot2::element_text(vjust = -4, size = 14, margin = ggplot2::margin(t = -5, b = 20)),
                 axis.title = ggplot2::element_text(size = 16, face = "bold", color = "black"),
                 # axis.title.x = ggplot2::element_text(margin = margin(t = 10, b = 5)),
                 axis.title.y = ggplot2::element_text(margin = ggplot2::margin(l = 5, r = 10)))

## In two adjacent plots 

env_plot <- ggplot(data = env_coords, mapping = aes(x = PC1, y = PC2)) +
  geom_bin2d(bins = 50, binwidth = c(0.2, 0.16)) +
  
  xlim(-7, 3) +
  ylim(-5, 3) +
  xlab(paste0("PC1 (",PC1_var,"%)")) + 
  ylab(paste0("PC2 (",PC2_var,"%)")) +
  
  scale_fill_gradientn(colours = grey.colors(n = 200, start = 0.8, end = 0.2), na.value = NA) +
  labs(fill = "Env.") +
  
  scale_y_continuous(minor_breaks = seq(-7.05 , 3.95, 0.16)) +
  scale_x_continuous(minor_breaks = seq(-7 , 4, 0.2)) +
  labs(title = "Environmental space") +
  
  annotate("text", x = 2.5, y = -4.3, hjust = 0, vjust = 1,
           label = paste0("(a)"),
           size = 7, fontface = 2) +
  
  guides(fill = FALSE) +
  
  ggplot2::theme(panel.background = ggplot2::element_rect(fill = NA),
                 plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5, vjust = 3),
                 # legend.position = c(0.85, 0.55),
                 legend.title = ggplot2::element_text(size = 12, vjust = 3, face = "bold"),
                 legend.text = ggplot2::element_text(size = 12, face = "plain"),
                 panel.grid.major = element_line(colour = "grey50", size = 0.3),
                 panel.grid.minor = element_line(colour = "grey90", size = 0.3),
                 plot.margin = ggplot2::margin(t = 18, b = 10, unit = "pt"),
                 # axis.line.x = ggplot2::element_line(color = NA),
                 axis.line = ggplot2::element_line(color = "black", size = 1.3),
                 axis.ticks = ggplot2::element_line(color = "black", size = 1.2),
                 axis.ticks.length = ggplot2::unit(8, "pt"),
                 axis.text = ggplot2::element_text(size = 14, face = "bold", color = "black"),
                 axis.text.x = ggplot2::element_text(vjust = -4, size = 14, margin = ggplot2::margin(t = -5, b = 20)),
                 axis.title = ggplot2::element_text(size = 16, face = "bold", color = "black"),
                 # axis.title.x = ggplot2::element_text(margin = margin(t = 10, b = 5)),
                 axis.title.y = ggplot2::element_text(margin = ggplot2::margin(l = 5, r = 10))) 

print(env_plot)

occ_plot <- ggplot(data = env_coords, mapping = aes(x = PC1, y = PC2)) +
  geom_bin2d(bins = 50, binwidth = c(0.2, 0.16), alpha = 1) + # Full transparent. Only here to get the scale
  
  xlim(-7, 3) +
  ylim(-5, 3) +
  xlab(paste0("PC1 (",PC1_var,"%)")) + 
  # ylab(paste0("PC2 (",PC2_var,"%)")) +
  ylab("") +
  
  scale_fill_gradientn(colours = grey.colors(n = 200, start = 0.8, end = 0.2), na.value = NA) +
  labs(fill = "Env.           ") + # Need space because order of legend is determined by number of character in the title...
  
  ggnewscale::new_scale_fill() +
  geom_bin2d(data = occ_coords, aes(x = PC1, y = PC2),
             bins = 50, binwidth = c(0.2, 0.16), alpha = 1) + 
  
  # scale_fill_gradientn(colours = rev(tmaptools::get_brewer_pal("Spectral", n = 100, plot = F)), na.value = pal_bl_red_Mannion[1]) +
  # scale_fill_gradientn(colours = pal_bl_red_Mannion[150:200], na.value = pal_bl_red_Mannion[1]) +
  scale_fill_gradientn(colours = tmaptools::get_brewer_pal("YlOrRd", n = 100, plot = F), na.value = pal_bl_red_Mannion[1]) +
  labs(fill = "\n\nSampling") +
  
  scale_y_continuous(minor_breaks = seq(-7.05 , 3.95, 0.16)) +
  scale_x_continuous(minor_breaks = seq(-7 , 4, 0.2)) +
  labs(title = "Sampling effort in environmental space") +
  
  annotate("text", x = 2.5, y = -4.3, hjust = 0, vjust = 1,
           label = paste0("(b)"),
           size = 7, fontface = 2) +

  
  # guides(Env = guide_legend(order = 1),
  #        Sampling = guide_legend(order = 2)) +
  
  ggplot2::theme(panel.background = ggplot2::element_rect(fill = NA),
                 plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5, vjust = 3),
                 # legend.position = c(0.85, 0.55),
                 legend.title = ggplot2::element_text(size = 12, vjust = 3, face = "bold"),
                 legend.text = ggplot2::element_text(size = 12, face = "plain"),
                 panel.grid.major = element_line(colour = "grey50", size = 0.3),
                 panel.grid.minor = element_line(colour = "grey90", size = 0.3),
                 plot.margin = ggplot2::margin(t = 18, b = 10, unit = "pt"),
                 # axis.line.x = ggplot2::element_line(color = NA),
                 axis.line = ggplot2::element_line(color = "black", size = 1.3),
                 axis.ticks = ggplot2::element_line(color = "black", size = 1.2),
                 axis.ticks.length = ggplot2::unit(8, "pt"),
                 axis.text = ggplot2::element_text(size = 14, face = "bold", color = "black"),
                 axis.text.x = ggplot2::element_text(vjust = -4, size = 14, margin = ggplot2::margin(t = -5, b = 20)),
                 axis.title = ggplot2::element_text(size = 16, face = "bold", color = "black"),
                 # axis.title.x = ggplot2::element_text(margin = margin(t = 10, b = 5)),
                 axis.title.y = ggplot2::element_text(margin = ggplot2::margin(l = 5, r = 10))) 

print(occ_plot)



### Correlation test ###

# Extract data for environment
env_plot_data <- ggplot(data = env_coords, mapping = aes(x = PC1, y = PC2)) +
  geom_bin2d(bins = 50, binwidth = c(0.2, 0.16), alpha = 1) # Full transparent. Only here to get the scale
env_plot_data <- ggplot_build(env_plot_data)$data[[1]]

env_plot_data_counts <- env_plot_data %>% 
  select(x, y, count) %>% 
  mutate("PC1" = round(x, 2)) %>% 
  mutate("PC2" = round(y, 2)) %>% 
  select(PC1, PC2, count)

# Extract data for sampling
occ_plot_data <- ggplot(data = occ_coords, mapping = aes(x = PC1, y = PC2)) +
  geom_bin2d(bins = 50, binwidth = c(0.2, 0.16), alpha = 1) # Full transparent. Only here to get the scale
occ_plot_data <- ggplot_build(occ_plot_data)$data[[1]]

occ_plot_data_counts <- occ_plot_data %>% 
  select(x, y, count) %>% 
  mutate("PC1" = round(x, 2)) %>% 
  mutate("PC2" = round(y, 2)) %>% 
  select(PC1, PC2, count)

# Join data
counts_df <- full_join(x = env_plot_data_counts, y = occ_plot_data_counts, by = c("PC1", "PC2")) ####
counts_df <- counts_df %>% 
  rename("count_env" = "count.x") %>% 
  rename("count_occ" = "count.y") 
counts_df$count_occ[is.na(counts_df$count_occ)] <- 0

# Extract only occurrence within Ithomiini range
counts_df_range <- counts_df[counts_df$count_occ > 0, ]


## Compute t-test without correction for spatial autocorrelation

compute_clifford_correlation_test <- function(x, y, coords, method = c("spearman", "pearson", "both"))
{
  # Sample size # Must be equal for x and y
  N <- length(x)

  # Use Struges's formula to define the number of classes for correlograms
  Ht = round(1 + 3.322 * log(N), 0)

  if (method %in% c("pearson", "both"))
  {
    # For Pearson's r : uncorrected test
    pearson_test <- cor.test(x, y, alternative = "two.sided", method = "pearson")
    pearson_r <- pearson_test$estimate # Extract the Pearson's R
    pearson_t_uncorrected <- round(pearson_r * sqrt(N - 2) / sqrt(1 - pearson_r^2), digits = 3) #  Uncorrected t-stat for Pearson's r
    pearson_p_uncorrected <- round(2 * (1 - pt(q = abs(pearson_t_uncorrected), df = N - 2)), digits = 3) # Uncorrected p-value for two-sided test
    
    # Transform points into a polynomial surface trend
    pearson_x_surf <- spatial::surf.ls(0, coords[,1], coords[,2], x)
    pearson_y_surf <- spatial::surf.ls(0, coords[,1], coords[,2], y)
    
    # Run spatial autocorrelation correlogram
    pearson_x_correlog <- spatial::correlogram(krig = pearson_x_surf, nint = Ht, plotit = F)
    pearson_y_correlog <- spatial::correlogram(krig = pearson_y_surf, nint = Ht, plotit = F)
    
    # Choose the max distance class to look at spatial autocorrelation depending on the plot of the correlograms
    # After a certain distance, autocorrelation should be null, thus no need to include the next classes in the computation
    
    # Get automatically the first class with negative autocorrelation for both variables
    pearson_H <- max(which.max(pearson_x_correlog$y < 0), which.max(pearson_y_correlog$y < 0))
    
    # Record distance of absence of spatial autocorrelation
    pearson_threshold_dist <- round(max(pearson_x_correlog$x[pearson_H], pearson_y_correlog$x[pearson_H]), digits = 4)

    ## Estimate the variance of the sampling distribution (sampling error) of Pearson's r from the first H classes where autocorrelation is still present 
    pearson_nr1r2 <- pearson_x_correlog$cnt * pearson_x_correlog$y * pearson_y_correlog$y
    pearson_sr_hat <- sum(pearson_nr1r2[1:pearson_H]) * N^(-2)
    
    ## Compute corrected sample size
    pearson_ne_hat <- round(1 + 1 / pearson_sr_hat, digits = 3)
    
    ## Compute the corrected t-test with corrected sample size for the df
    pearson_t_corrected <- round(pearson_r * sqrt(pearson_ne_hat - 2) / sqrt(1 - pearson_r^2), digits = 3) #  Corrected t-stat
    pearson_p_corrected <- round(2 * (1 - pt(q = abs(pearson_t_corrected), df = pearson_ne_hat - 2)), digits = 3) # Corrected p-value for two-sided test
 
  }
  
  
  if (method %in% c("spearman", "both"))
  {
    # For Spearman's rho: uncorrected test
    spearman_test <- cor.test(x, y, alternative = "two.sided", method = "spearman")
    spearman_rho <- spearman_test$estimate # Extract the Spearman's rho
    spearman_t_uncorrected <- round(spearman_rho * sqrt(N - 2) / sqrt(1 - spearman_rho^2), digits = 3) #  Uncorrected t-stat for Spearman's rho
    spearman_p_uncorrected <- round(2 * (1 - pt(q = abs(spearman_t_uncorrected), df = N - 2)), digits = 3) # Uncorrected p-value for two-sided test
    
    # Compute ranks for Spearman's test
    x_rank <- rank(x)
    y_rank <- rank(y)
    
    # Transform ranks into a polynomial surface trend
    spearman_x_surf <- spatial::surf.ls(0, coords[,1], coords[,2], x_rank)
    spearman_y_surf <- spatial::surf.ls(0, coords[,1], coords[,2], y_rank)
    
    # Run spatial autocorrelation correlogram
    spearman_x_correlog <- spatial::correlogram(krig = spearman_x_surf, nint = Ht, plotit = F)
    spearman_y_correlog <- spatial::correlogram(krig = spearman_y_surf, nint = Ht, plotit = F)
    
    # Choose the max distance class to look at spatial autocorrelation depending on the plot of the correlograms
    # After a certain distance, autocorrelation should be null, thus no need to include the next classes in the computation
    
    # Get automatically the first class with negative autocorrelation for both variables
    spearman_H <- max(which.max(spearman_x_correlog$y < 0), which.max(spearman_y_correlog$y < 0))
    
    # Record distance of absence of spatial autocorrelation
    spearman_threshold_dist <- round(max(spearman_x_correlog$x[spearman_H], spearman_y_correlog$x[spearman_H]), digits = 4)
    
    ## Estimate the variance of the sampling distribution (sampling error) of Pearson's r from the first H classes where autocorrelation is still present 
    spearman_nr1r2 <- spearman_x_correlog$cnt * spearman_x_correlog$y * spearman_y_correlog$y
    spearman_sr_hat <- sum(spearman_nr1r2[1:spearman_H]) * N^(-2)
    
    ## Compute corrected sample size
    spearman_ne_hat <- round(1 + 1 / spearman_sr_hat, digits = 3)
    
    ## Compute the corrected t-test with corrected sample size for the df
    spearman_t_corrected <- round(spearman_rho * sqrt(spearman_ne_hat - 2) / sqrt(1 - spearman_rho^2), digits = 3) #  Corrected t-stat
    spearman_p_corrected <- round(2 * (1 - pt(q = abs(spearman_t_corrected), df = spearman_ne_hat - 2)), digits = 3) # Corrected p-value for two-sided test

  }
  
  # Built output
  if (method %in% c("pearson", "both"))
  {
    pearson_output <- data.frame("Coeff" = round(pearson_r, 3), 
                                 "T_stat_uncorrected" = pearson_t_uncorrected, 
                                 "P_value_uncorrected" = pearson_p_uncorrected, 
                                 "Threshold_autocorrelation_dist" = pearson_threshold_dist, 
                                 "T_stat_corrected" = pearson_t_corrected,
                                 "DF_corrected" = round(pearson_ne_hat, 1),
                                 "P_value_corrected" = pearson_p_corrected)
    row.names(pearson_output) <- "Pearson"
  }
  
  if (method == "pearson") { return(pearson_output) }

  if (method %in% c("spearman", "both"))
  {
    spearman_output <- data.frame("Coeff" = round(spearman_rho, 3), 
                                 "T_stat_uncorrected" = spearman_t_uncorrected, 
                                 "P_value_uncorrected" = spearman_p_uncorrected, 
                                 "Threshold_autocorrelation_dist" = spearman_threshold_dist, 
                                 "T_stat_corrected" = spearman_t_corrected,
                                 "DF_corrected" = round(spearman_ne_hat, 1),
                                 "P_value_corrected" = spearman_p_corrected)
    row.names(spearman_output) <- "Spearman"
  }
  
  if (method == "spearman") { return(spearman_output) }

  if (method == "both")
  {
    both_output <- rbind(pearson_output, spearman_output)
    return(both_output)
  }
  
}

spearman_corrected_test <- compute_clifford_correlation_test(x = counts_df_range$count_env,
                                                             y = counts_df_range$count_occ,
                                                             coords = counts_df_range[, c("PC1", "PC2")], 
                                                             method = "spearman")

### Add correlation test to the plot

print(occ_plot)

occ_plot <- occ_plot +
  annotate("rect", xmin = -6.2, xmax = -4.20, ymin = 1.75, ymax = 3.2,
           col = "grey20", fill = "white", alpha = 1) +
  annotate("text", x = -6.0, y = 3.0, hjust = 0, vjust = 1,
           label = bquote(rho ~ "= " ~ .(spearman_corrected_test$Coeff)),
           size = 4, fontface = 2) +
  annotate("text", x = -6.0, y = 2.6, hjust = 0, vjust = 1,
           label = paste0("df = ", spearman_corrected_test$DF_corrected),
           size = 4, fontface = 2) +
  annotate("text", x = -6.0, y = 2.2, hjust = 0, vjust = 1,
           label = paste0("p < 0.001 "),
           size = 4, fontface = 2)

### Generate final multiplot

library(gridExtra)
library(grid)

pdf(file = "./supplementaries/Sampling_map_env.pdf", width = 12, height = 6)

grid.arrange(
  grobs = list(env_plot, occ_plot), # List of ggplots
  widths = c(4.2, 0, 5),  # Width of columns
  heights = c(4),
  nrow = 1,
  ncol = 3,
  layout_matrix = rbind(c(1, NA, 2))) # Position of ggplots in the layout

dev.off()


# Correlation circle but with style! 
# May be added in a Viewport, but probably not so useful
fviz_pca_var(X = PCA_res, axes = c(1, 2), 
             geom.var = c("arrow", "point", "text"))




##### 2/ Thiessen polygon networks for sampling density #####

# See Jenkins et al., 2015 - Patterns of vertebrate diversity and protection in Brazil

### 2.1/ Build the Thiessen polygons ####

library("sf")
library("units")

### Load country border file
load(file = "./input_data/Map_stuff/country_borders.RData")

plot(country_borders)

# Transform into a multipolygon of land borders
country_borders_merged <- st_union(country_borders)[2]
plot(country_borders_merged)

# Cast into several polygon instead of a single multipolygon
country_borders_merged_poly <- st_cast(country_borders_merged, "POLYGON")

# Compute area
land_poly_sf <- st_sf(data.frame("Area" = set_units(st_area(country_borders_merged_poly), km^2)),
                      geometry = country_borders_merged_poly)

plot(land_poly_sf, col = land_poly_sf$Area)

# Remove holes
land_poly_sf <- nngeo::st_remove_holes(land_poly_sf)

# Save
saveRDS(object = land_poly_sf, file = "./input_data/Map_stuff/land_poly_sf.rds")

plot(land_poly_sf, col = land_poly_sf$Area)

# Filter by area for fun
land_poly_sf_filtered <- dplyr::filter(land_poly_sf, Area > set_units(1000, km^2))
plot(land_poly_sf_filtered, col = land_poly_sf$Area)


### Load occurrence data
occ_df <- Ithomiini_final[, c("Longitude", "Latitude")]
occ_sfg <- st_multipoint(x = as.matrix(occ_df))
occ_sfc <- st_sfc(occ_sfg, crs = 4326)

plot(occ_sfc)

# Data need to be projected because Voronoï thesselation works only on projected data

occ_Mollweide <- st_transform(x = occ_sfc, crs = CRS("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
country_borders_merged_Mollweide <- st_transform(x = country_borders_merged, crs = CRS("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))
land_poly_sf_Mollweide <- st_transform(x = land_poly_sf, crs = CRS("+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs"))

Voronoi_cells_Mollweide <- st_voronoi(x = occ_Mollweide, 
                                      envelope = NA) # Envelop does not work. Need a single POLYGON. Better to apply intersection afterwards
plot(Voronoi_cells_Mollweide, col = "white", border = "black")

Voronoi_cells_Mollweide_poly <- st_cast(Voronoi_cells_Mollweide)
plot(Voronoi_cells_Mollweide_poly)

Voronoi_cells <- st_transform(x = Voronoi_cells_Mollweide_poly, crs = 4326, aoi = c(-115,-37,-35,28))
plot(Voronoi_cells)

plot(occ_Mollweide, pch = 16, col = "red", cex = 0.5, add = T)

### Intersect to keep only terrestrial lands and nullify values beyond Ithomiini range

sp_richness <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
Ithomiini_range <- sp_richness > 0
Ithomiini_range_Mollweide <- projectRaster(from = Ithomiini_range, 
                                           method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                           crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                           alignOnly = F)
Ithomiini_range_Mollweide <- Ithomiini_range_Mollweide > 0
plot(Ithomiini_range_Mollweide)

Ithomiini_range_poly_sp_Mollweide <- raster::rasterToPolygons(x = Ithomiini_range_Mollweide, 
                                                              fun = function(x) { (x == 1) },
                                                              dissolve = T) # To merge contiguous polygons of the same field value

plot(Ithomiini_range_poly_sp_Mollweide)

library("sf")
library("units")

Ithomiini_range_poly_sf_Mollweide <- st_as_sf(Ithomiini_range_poly_sp_Mollweide)

Voronoi_cells_Mollweide_land_clip <- st_intersection(x = land_poly_sf_Mollweide, y = Voronoi_cells_Mollweide_poly)
Voronoi_cells_Mollweide_range_clip <- st_intersection(x = Ithomiini_range_poly_sf_Mollweide, y = Voronoi_cells_Mollweide_poly)

Voronoi_cells_Mollweide_range_clip$Area <- set_units(st_area(Voronoi_cells_Mollweide_range_clip), km^2)
plot(Voronoi_cells_Mollweide_range_clip)

Voronoi_cells_Mollweide_range_clip_sp <- as(Voronoi_cells_Mollweide_range_clip, "Spatial")
plot(Voronoi_cells_Mollweide_range_clip_sp, col = Voronoi_cells_Mollweide_range_clip_sp$Area)

saveRDS(Voronoi_cells_Mollweide_range_clip_sp, file = "./maps/Sampling_effort/Voronoi_cells_Mollweide_range_clip_sp.rds")

### Rasterize to have a background for legend
library("raster")

continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_15.rds"))
continent_mask_Mollweide <- projectRaster(from = continent_mask, 
                                          method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                                          crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                                          alignOnly = F)

Voronoi_cells_Mollweide_range_clip_raster <- rasterize(x = Voronoi_cells_Mollweide_range_clip_sp, y = continent_mask_Mollweide,
                                                       field = "Area", background = NA, fun = max)

temp <- continent_mask_Mollweide
temp[!is.na(Voronoi_cells_Mollweide_range_clip_raster@data@values)] <- Voronoi_cells_Mollweide_range_clip_raster@data@values[!is.na(Voronoi_cells_Mollweide_range_clip_raster@data@values)]
Voronoi_cells_Mollweide_range_clip_raster <- temp

plot(Voronoi_cells_Mollweide_range_clip_raster, col = c(pal_bl_red_Mannion[1],rev(pal_bl_red_Mannion[2:200])))
plot(Voronoi_cells_Mollweide_range_clip_raster, col = tmaptools::get_brewer_pal("Spectral", n = 112, plot = F)[8:107])

saveRDS(Voronoi_cells_Mollweide_range_clip_raster, file = "./maps/Sampling_effort/Voronoi_cells_Mollweide_range_clip_raster.rds")
Voronoi_cells_Mollweide_range_clip_raster <- readRDS(file = "./maps/Sampling_effort/Voronoi_cells_Mollweide_range_clip_raster.rds")


col_breaks <- BAMMtools::getJenksBreaks(var = Voronoi_cells_Mollweide_range_clip_sp$Area, k = 101)

col_factors_raster <- cut(x = Voronoi_cells_Mollweide_range_clip_raster@data@values, breaks = col_breaks)

temp <- Voronoi_cells_Mollweide_range_clip_raster
temp@data@values <- as.numeric(col_factors_raster)

plot(temp, col = rev(pal_bl_red_Mannion[-1]))
plot(temp, col = rev(tmaptools::get_brewer_pal("YlOrRd", n = 100, plot = F)))
plot(temp, col = rev(grey.colors(n = 200, start = 0.8, end = 0.2)))
plot(temp, col = tmaptools::get_brewer_pal("Spectral", n = 100, plot = F)[8:95])

col_factors_sp <- cut(x = Voronoi_cells_Mollweide_range_clip_sp$Area, breaks = col_breaks)


table(as.numeric(col_factors_sp))

plot(Voronoi_cells_Mollweide_range_clip_sp, border = NA,
     col = tmaptools::get_brewer_pal("Spectral", n = 112, plot = F)[8:107][as.numeric(col_factors_sp)])


### 2.2/ Map the sampling effort ####

# Packages
library(raster)
library(prettymapr)
library(rangeBuilder)
library(sf)

# New color palette
pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")

# Load mask for continent borders, plot border, and grid
grid_Mollweide_out <- readRDS(file = "./input_data/Map_stuff/grid_Mollweide_out.rds")
large_bg_mask_Mollweide <- readRDS(file = "./input_data/Map_stuff/large_bg_mask_Mollweide.rds")
bbox_sp_Mollweide <- readRDS(file = "./input_data/Map_stuff/bbox_sp_Mollweide.rds")

load(file = "./input_data/Map_stuff/country_borders.RData")
country_borders <- as(country_borders, "Spatial")

rivers <- st_read("./input_data/Map_stuff/Major_rivers/MajorRivers.shp") # Load rivers
rivers <- as(rivers, "Spatial")
xmin = -120 ; xmax = -30 ; ymin = -37 ; ymax = 28
e <- extent(c(xmin,xmax,ymin,ymax))
rivers <- crop(rivers, e)

# Load the shp, color scheme and raster of sampling effort
Voronoi_cells_Mollweide_range_clip_sp <- readRDS(file = "./maps/Sampling_effort/Voronoi_cells_Mollweide_range_clip_sp.rds")
col_breaks <- BAMMtools::getJenksBreaks(var = Voronoi_cells_Mollweide_range_clip_sp$Area, k = 101)
col_factors_sp <- cut(x = Voronoi_cells_Mollweide_range_clip_sp$Area, breaks = col_breaks)
Voronoi_cells_Mollweide_range_clip_raster <- readRDS(file = "./maps/Sampling_effort/Voronoi_cells_Mollweide_range_clip_raster.rds")


### Function to project shp object
Mollweide_shp_projection <-  function(x) # Shp to project
{
  x_name <- deparse(substitute(x)) # Get the name of the initial shp as a character string
  
  new_shp <- spTransform(x, CRSobj = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
  
  # Generate new object with "_Mollweide" suffix in the global environment
  eval(call("<<-", as.name(paste0(x_name, "_Mollweide")), new_shp))
}

# Apply projection
Mollweide_shp_projection(country_borders)
Mollweide_shp_projection(rivers)


### Plot sampling effort

pdf(file = paste0("./supplementaries/Sampling_effort.pdf"), height = 6, width = 6)

internal_margins <- par()$mar
par(mar = c(3.1, 3.1, 1.5, 1.6))

# Plot raster background without axis
image(x = continent_mask_Mollweide, col = pal_bl_red_Mannion[1],
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "")

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.4, lwd = 0.2, lwd.ticks = 1, gap.axis = 0, padj = 0.5)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.4, lwd = 0.2, lwd.ticks = 1, gap.axis = 0)

# Add background, borders and graticules
plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)

# Add sampling effort
plot(Voronoi_cells_Mollweide_range_clip_sp, border = NA, add = T,
     col = tmaptools::get_brewer_pal("Spectral", n = 112, plot = F)[8:107][as.numeric(col_factors_sp)])

# Add borders and rivers
# plot(country_borders_Mollweide, lwd = 1, border = "#FF000030", col = NA, add = T)
# plot(rivers_Mollweide, lwd = 1, col = "#6495ED60", add = T)

# Add occurrence points
# plot(Occ_shp_Mollweide, pch = 16, cex = 0.5, col = "#00000090", add = T)

# Add scale bar in legend
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-2500, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1.1)
prettymapr::addnortharrow(scale = 0.45, padin = c(0.15, 0.15), text.col = "#00000000")
# legend(legend = "Country borders", x = -4450, y = -2350, cex = 1.1, lty = 1, lwd = 2, col = "red", bty = "n")
# legend(legend = "Country borders", x = -4450, y = -2350, cex = 1.1, lty = 1, lwd = 2, col = "red", bty = "n")
# legend(legend = "Rivers", x = -4450, y = -2850, cex = 1.1, lty = 1, lwd = 2, col = "#6495ED", bty = "n")
# legend(legend = "Rivers", x = -4450, y = -2850, cex = 1.1, lty = 1, lwd = 2, col = "#6495ED", bty = "n")

# Voronoi_cells_Mollweide_range_clip_raster_inverse <- Voronoi_cells_Mollweide_range_clip_raster
# Voronoi_cells_Mollweide_range_clip_raster_inverse@data@values <- 1/Voronoi_cells_Mollweide_range_clip_raster@data@values
# 
# hist(Voronoi_cells_Mollweide_range_clip_raster@data@values)

# Add legend
rangeBuilder::addRasterLegend(Voronoi_cells_Mollweide_range_clip_raster, nTicks = 3,
                              labelDist = 500, # To remove the value from the plot because they do not fit the actual scale
                              cex.axis = 1.4, ramp = rev(tmaptools::get_brewer_pal("Spectral", n = 112, plot = F)[8:107]),
                              ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))

# Quartile 1%
quartile_1 <- round(((1/col_breaks)*(10^6))[101],1)
graphics::text(x = -3500, y = -3900, font = 1, cex = 1.2, label = quartile_1, adj = 0)

# Quartile 25%
quartile_25 <- round(((1/col_breaks)*(10^6))[75],1)
graphics::text(x = -3500, y = -2975, font = 1, cex = 1.2, label = quartile_25, adj = 0)

# Quartile 50 %
quartile_50 <- round(((1/col_breaks)*(10^6))[50],1)
graphics::text(x = -3500, y = -1975, font = 1, cex = 1.2, label = quartile_50, adj = 0)

# Quartile 75%
quartile_75 <- round(((1/col_breaks)*(10^6))[25],1)
graphics::text(x = -3500, y = -975, font = 1, cex = 1.2, label = quartile_75, adj = 0)

# Quartile 99%
quartile_99 <- round(((1/col_breaks)*(10^6))[2],0)
graphics::text(x = -3500, y = 0, font = 1, cex = 1.2, label = quartile_99, adj = 0)

graphics::text(x = -3000, y = 975, font = 2, cex = 1.2, label = "Sampling density")
graphics::text(x = -3150, y = 550, font = 2, cex = 1.2, label = expression(paste(bold("[10")^bold("-6"), bold(" km")^bold("-2"),bold("]"))))

# # Plot cities on the map
# 
# target_coords <- data.frame(t(c(-49.1, -5.34)))  # Maraba
# target_coords <- data.frame(t(c(-60.16, -3.04))) # Manaus
# target_coords <- data.frame(t(c(-54.74, -2.47))) # Santarem
# 
# target_coords_sp <- SpatialPoints(coords = target_coords, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# Mollweide_shp_projection(target_coords_sp)
# 
# plot(target_coords_sp_Mollweide, pch = 16, cex = 1, col = "black", add = T)

par(mar = internal_margins)

dev.off()


hist(Voronoi_cells_Mollweide_range_clip_raster@data@values)


##### 3/ Map occurrences and Sampling effort aside ####

### Load extra stuff for occurrence map

# Load contiental mask
continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_15.rds"))

# Generate spatial point object of occurrences 
Occ_shp <- SpatialPoints(coords = Ithomiini_final[, c("Longitude", "Latitude")],
                         proj4string = continent_mask@crs)
# Filter out occurrences falling in the sea
Occ_shp <- rgeos::gIntersection(spgeom1 = country_borders, spgeom2 = Occ_shp)

# Function to project raster
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
}

# Project stuff
Mollweide_projection(continent_mask)
Mollweide_shp_projection(Occ_shp)


### Plot final figure

pdf(file = paste0("./supplementaries/Figure_S1.1_Occurrrence_&_Sampling_effort.pdf"), height = 6, width = 12)

internal_margins <- par()$mar
par(mar = c(3.1, 3.1, 3.5, 1.6))
par(mfrow = c(1,2))

### Panel A = Occurrence map ###

# Plot raster background without axis
image(x = continent_mask_Mollweide, col = "antiquewhite",
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "",
      main = "Occurrence map", cex.main = 1.8)

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.4, lwd = 0.2, lwd.ticks = 1, gap.axis = 0, padj = 0.5)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.4, lwd = 0.2, lwd.ticks = 1, gap.axis = 0)

# Add background, borders and graticules
plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)
plot(country_borders_Mollweide, lwd = 1, border = "#FF000030", col = NA, add = T)
plot(rivers_Mollweide, lwd = 1, col = "#6495ED60", add = T)

# Add occurrence points
plot(Occ_shp_Mollweide, pch = 16, cex = 0.5, col = "#00000090", add = T)

# Add scale bar in legend
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-4000, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1.1)
prettymapr::addnortharrow(scale = 0.45, padin = c(0.15, 0.15), text.col = "#00000000")
legend(legend = "Country borders", x = -4450, y = -2350, cex = 1.1, lty = 1, lwd = 2, col = "red", bty = "n")
legend(legend = "Country borders", x = -4450, y = -2350, cex = 1.1, lty = 1, lwd = 2, col = "red", bty = "n")
legend(legend = "Rivers", x = -4450, y = -2850, cex = 1.1, lty = 1, lwd = 2, col = "#6495ED", bty = "n")
legend(legend = "Rivers", x = -4450, y = -2850, cex = 1.1, lty = 1, lwd = 2, col = "#6495ED", bty = "n")

# Add panel label
text(x = 4000, y = -3900, labels = "(a)", font = 2, cex = 1.8)


### Panel B = Sampling effort ###

# Plot raster background without axis
image(x = continent_mask_Mollweide, col = pal_bl_red_Mannion[1],
      xlim = c(-4600, 4600), ylim = c(-4450, 3400), axes = F,
      xlab = "", ylab = "",
      main = "Sampling effort", cex.main = 1.8)

# Generate axes with manual positioning of ticks
axis(side = 1, at = c(-3930, -2170, -420, 1500, 3050), labels = c("120°E", "100°E", "80°E", "60°E", "40°E"), cex.axis = 1.4, lwd = 0.2, lwd.ticks = 1, gap.axis = 0, padj = 0.5)
axis(side = 2, at = c(-3650, -2450, -1220, 0, 1230, 2445), labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"), cex.axis = 1.4, lwd = 0.2, lwd.ticks = 1, gap.axis = 0)

# Add background, borders and graticules
plot(large_bg_mask_Mollweide, lwd = 1, border = "grey20", col = "aliceblue", add = T)
plot(grid_Mollweide_out, lty = 92, col = "grey80", add = T)
plot(bbox_sp_Mollweide, lwd = 2, border = "black", col = NA, add = T)

# Add sampling effort
plot(Voronoi_cells_Mollweide_range_clip_sp, border = NA, add = T,
     col = tmaptools::get_brewer_pal("Spectral", n = 112, plot = F)[8:107][as.numeric(col_factors_sp)])

# Add borders and rivers
# plot(country_borders_Mollweide, lwd = 1, border = "#FF000030", col = NA, add = T)
# plot(rivers_Mollweide, lwd = 1, col = "#6495ED60", add = T)

# Add occurrence points
# plot(Occ_shp_Mollweide, pch = 16, cex = 0.5, col = "#00000090", add = T)

# Add scale bar in legend
scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = c(-2500, -4000), label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1.1)
prettymapr::addnortharrow(scale = 0.45, padin = c(0.15, 0.15), text.col = "#00000000")
# legend(legend = "Country borders", x = -4450, y = -2350, cex = 1.1, lty = 1, lwd = 2, col = "red", bty = "n")
# legend(legend = "Country borders", x = -4450, y = -2350, cex = 1.1, lty = 1, lwd = 2, col = "red", bty = "n")
# legend(legend = "Rivers", x = -4450, y = -2850, cex = 1.1, lty = 1, lwd = 2, col = "#6495ED", bty = "n")
# legend(legend = "Rivers", x = -4450, y = -2850, cex = 1.1, lty = 1, lwd = 2, col = "#6495ED", bty = "n")

# Voronoi_cells_Mollweide_range_clip_raster_inverse <- Voronoi_cells_Mollweide_range_clip_raster
# Voronoi_cells_Mollweide_range_clip_raster_inverse@data@values <- 1/Voronoi_cells_Mollweide_range_clip_raster@data@values
# 
# hist(Voronoi_cells_Mollweide_range_clip_raster@data@values)

# Add legend
rangeBuilder::addRasterLegend(Voronoi_cells_Mollweide_range_clip_raster, nTicks = 3,
                              labelDist = 500, # To remove the value from the plot because they do not fit the actual scale
                              cex.axis = 1.4, ramp = rev(tmaptools::get_brewer_pal("Spectral", n = 112, plot = F)[8:107]),
                              ncolors = 200, border = T, location = c(-4100, -3800, -3950, 0))

# Quartile 1%
quartile_1 <- round(((1/col_breaks)*(10^6))[101],1)
graphics::text(x = -3500, y = -3900, font = 1, cex = 1.2, label = quartile_1, adj = 0)

# Quartile 25%
quartile_25 <- round(((1/col_breaks)*(10^6))[75],1)
graphics::text(x = -3500, y = -2975, font = 1, cex = 1.2, label = quartile_25, adj = 0)

# Quartile 50 %
quartile_50 <- round(((1/col_breaks)*(10^6))[50],1)
graphics::text(x = -3500, y = -1975, font = 1, cex = 1.2, label = quartile_50, adj = 0)

# Quartile 75%
quartile_75 <- round(((1/col_breaks)*(10^6))[25],1)
graphics::text(x = -3500, y = -975, font = 1, cex = 1.2, label = quartile_75, adj = 0)

# Quartile 99%
quartile_99 <- round(((1/col_breaks)*(10^6))[2],0)
graphics::text(x = -3500, y = 0, font = 1, cex = 1.2, label = quartile_99, adj = 0)

graphics::text(x = -3000, y = 975, font = 2, cex = 1.2, label = "Sampling density")
graphics::text(x = -3150, y = 550, font = 2, cex = 1.2, label = expression(paste(bold("[10")^bold("-6"), bold(" km")^bold("-2"),bold("]"))))

# # Plot cities on the map
# 
# target_coords <- data.frame(t(c(-49.1, -5.34)))  # Maraba
# target_coords <- data.frame(t(c(-60.16, -3.04))) # Manaus
# target_coords <- data.frame(t(c(-54.74, -2.47))) # Santarem
# 
# target_coords_sp <- SpatialPoints(coords = target_coords, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
# Mollweide_shp_projection(target_coords_sp)
# 
# plot(target_coords_sp_Mollweide, pch = 16, cex = 1, col = "black", add = T)

# Add panel label
text(x = 4000, y = -3900, labels = "(b)", font = 2, cex = 1.8)


par(mar = internal_margins)
par(mfrow = c(1,1))

dev.off()



##### 4/ Sampling completness #####

library(vegan)

?vegan::specpool

### 4.1/ Built community matrix of species/OMU incidence per sites ####

generate_incidence_matrix <- function(df, units, sites)
{
  units <- df[, names(df) == units, drop = T]
  sites <- df[, names(df) == sites, drop = T]
  
  units_list <- as.character(unique(units))
  sites_list <- as.character(unique(sites))
  
  # Retrieve dimensions of community matrix
  nb_row <- length(units_list)
  nb_col <- length(sites_list)
  
  # Initialize incidence matrix
  incidence_matrix <- matrix(data = 0, nrow = nb_row, ncol = nb_col)
  rownames(incidence_matrix) <- units_list
  colnames(incidence_matrix) <- sites_list
  
  # Fill the matrix
  for (k in 1:nrow(df))
  {
    i <- which(units[k] == units_list)
    j <- which(sites[k] == sites_list)
    
    incidence_matrix[i, j] <- 1
  }
  
  return(incidence_matrix)
}

sp_incidence_matrix <- generate_incidence_matrix(df = Ithomiini_final, units = "Sp_full", sites = "sampling_ID")
OMU_incidence_matrix <- generate_incidence_matrix(df = Ithomiini_final, units = "Tag", sites = "sampling_ID")

##### 4.2/ Global test ####
completeness_global <- vegan::specpool(x = t(sp_incidence_matrix)) ; completeness_global
completeness_global$Species/completeness_global$boot*100

 

##### 4.3/ Per grid cells ####

# Chose resolution in degrees
# resolution <- 1  # 1 degrees = 111km
# resolution <- 5  # 5 degrees = 556km
resolution <- 10 # 10 degrees = 1113km

compute_and_map_sampling_completness_per_grid <- function (incidence_matrix, resolution, index)
{
  # Load continental mask
  continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_15.rds"))
  # Project to new resolution
  sample_grid <- projectRaster(from = continent_mask, res = resolution,
                               method = "ngb", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                               crs = continent_mask@crs, # If you have the CRS arguments
                               alignOnly = F)
  # plot(sample_grid)
  
  # Extract coordinates of uniques sites
  sites_coords <- unique(Ithomiini_final[, c("sampling_ID", "Longitude_raster", "Latitude_raster")])[, c("Longitude_raster", "Latitude_raster")]
  
  # Extract grid ID per sampling sites
  grid_group <- raster::cellFromXY(sample_grid, xy = sites_coords)
  
  # Compute and export completness df
  completness_df <<- vegan::specpool(x = t(incidence_matrix), pool = grid_group) ; completness_df # Estimations de richesse spécifique théorique à partir d'un échantillon en utilisant les fréquences = Chao1 diversity index
  print(completness_df)
  
  # Compute complexness index
  completness_index <- round(completness_df$Species/completness_df[, names(completness_df) == index]*100, 1)
  
  completness_raster <- sample_grid
  completness_raster@data@values[unique(grid_group)[order(unique(grid_group))]] <- completness_index
  
  return(completness_raster)
}

completness_raster_grid <- compute_and_map_sampling_completness_per_grid(incidence_matrix = sp_incidence_matrix,
                                                                    resolution = 10, index = "boot")
plot(completness_raster_grid)

mean(completness_df$Species/completness_df$boot*100) 
sd(completness_df$Species/completness_df$boot*100)
range(completness_df$Species/completness_df$boot*100)

map_within_range <- function(map, final_res)
{
  # Load continental mask
  continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_",final_res,".rds"))
  
  # Load ithomiini range
  ithomiini_range <- readRDS(file = paste0("./input_data/Map_stuff/ithomiini_range_",final_res,"min.rds"))
  
  # Upgrade resolution
  final_res_map <- raster::resample(x = map, y = continent_mask, method = "ngb")
  
  # Mask borders
  final_res_map <- raster::mask(x = final_res_map, mask = ithomiini_range)
  
  # Add continental borders
  temp <- continent_mask
  temp[!is.na(final_res_map[])] <- final_res_map[!is.na(final_res_map[])]
  final_res_map <- temp
  
  return(final_res_map)
}

completness_raster_grid_within_range <- map_within_range(map = completness_raster_grid, final_res = "15")

library(tmaptools)
pal_red_grn_NA <- pal_red_grn <- get_brewer_pal("RdYlGn", n = 200, plot = F)
pal_red_grn_NA[1] <- "#EDEDED"

plot(completness_raster_grid_within_range, col = pal_red_grn_NA)

saveRDS(completness_raster_grid_within_range, file = "./outputs/Completeness_analyses/completness_raster_grid10_bootstrap.rds")


##### 4.4/ Per bioregions ####


### Build raster of bioregions

continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_15.rds"))
load(file = "./input_data/Map_stuff/Bioregions/All_bioregions_in_figure.RData")

bioregions_shp_list <- list(Guyana_Shield_shp, Caatinga_shp, Cerrado_shp, Mata_Atlantica_shp5, full_CA_shp, Northern_Andes_shp, Central_Andes_shp, Western_Lowlands_shp, Coastal_desert_shp, Chacos_shp, Llanos_shp, Caribbean_Islands_shp, Western_Amazon_shp, Lower_Amazon_shp, Pantanal_shp, Pampas_shp)
bioregions_names <- c("Guyana Shield", "Caatinga", "Cerrado", "Atlantic Forest", "Central America", "Northern Andes", "Central Andes", "Western Lowlands", "Coastal Desert", "Chacos", "Llanos", "Caribbean Islands", "Upper Amazon", "Lower Amazon", "Pantanal", "Pampas")

bioregion_stack <- stack()
for (i in 1:length(bioregions_shp_list))
{
  # i <- 1
  
  bioregion_shp <- bioregions_shp_list[[i]]
  
  bioregion_layer <- rasterize(x = bioregion_shp, y = continent_mask, field = i)
  
  bioregion_stack <- stack(bioregion_stack, bioregion_layer)
}

plot(bioregion_layer)
plot(bioregion_stack)

# Superimpose bioregions in a single layer
# bioregion_raster <- calc(bioregion_stack, fun = function(x) {as.numeric(names(sort(table(x), decreasing = TRUE))[1])})
bioregion_raster <- calc(bioregion_stack, fun = pracma::Mode)
plot(bioregion_raster)

?raster::approxNA

test <- approxNA(x = stack(bioregion_raster, bioregion_raster))
plot(test)

# Fill the gaps
width = 15 # Width of the focal windows to apply to look for values to fill the gap
bioregion_raster <- focal(x = bioregion_raster, w = matrix(1,width,width), fun = pracma::Mode, 
                          pad = F, NAonly = T)

# Limit to ithomiini_range
ithomiini_range <- readRDS(file = paste0("./input_data/Map_stuff/ithomiini_range_15min.rds"))
bioregion_raster <- mask(x = bioregion_raster, mask = ithomiini_range)

plot(bioregion_raster, col = tmaptools::get_brewer_pal("Spectral", n = 16, plot = F))

saveRDS(bioregion_raster, file = "./input_data/Map_stuff/Bioregions/bioregion_raster_15min.rds")
bioregion_raster <- readRDS(file = "./input_data/Map_stuff/Bioregions/bioregion_raster_15min.rds")

### Extract border shp of bioregions
plot(bioregion_shp)

bioregions_shp <- raster::rasterToPolygons(x = bioregion_raster, dissolve = T) # To merge contiguous polygons of the same field value

saveRDS(bioregions_shp, file = "./input_data/Map_stuff/Bioregions/bioregions_shp_15min.rds")
bioregions_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/bioregions_shp_15min.rds")

### Compute completeness per bioregions


compute_and_map_sampling_completness_per_bioregions <- function (incidence_matrix, index)
{
  # Load continental mask
  continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_15.rds"))

  # Extract coordinates of uniques sites
  sites_coords <- unique(Ithomiini_final[, c("sampling_ID", "Longitude_raster", "Latitude_raster")])[, c("Longitude_raster", "Latitude_raster")]
  
  # Extract bioregion ID per sampling sites
  grid_group <- raster::extract(bioregion_raster, y = sites_coords, buffer = 50000, fun = pracma::Mode)
  sum(is.na(grid_group))
  table(grid_group)
  
  # Compute and export completeness df
  completness_df <- vegan::specpool(x = t(incidence_matrix), pool = grid_group) ; completness_df # Estimations de richesse spécifique théorique à partir d'un échantillon en utilisant les fréquences = Chao1 diversity index
  completness_df <- data.frame("Bioregion" = bioregions_names, completness_df)
  completness_df <<- completness_df
  print(completness_df)
  
  # Compute completeness index
  completness_index <- round(completness_df$Species/completness_df[, names(completness_df) == index]*100, 1)
  
  # Replace bioregions ID by completness values
  completness_raster <- subs(x = bioregion_raster, y = data.frame("ID" = 1:16, "completness_index" = completness_index), by = "ID", which = "completness_index")
  
  # Add continental borders
  temp <- continent_mask
  temp[!is.na(completness_raster[])] <- completness_raster[!is.na(completness_raster[])]
  completness_raster <- temp
  
  return(completness_raster)
}

completness_raster_bioregions <- compute_and_map_sampling_completness_per_bioregions(incidence_matrix = sp_incidence_matrix,
                                                                                     index = "boot")
completness_df
mean(completness_df$Species/completness_df$boot*100) 
sd(completness_df$Species/completness_df$boot*100)
range(completness_df$Species/completness_df$boot*100)

plot(completness_raster_bioregions, col = pal_red_grn_NA)
plot(bioregions_shp, border = "grey20", lwd = 0.5, add = T)

saveRDS(completness_raster_bioregions, file = "./outputs/Completeness_analyses/completness_raster_bioregions_bootstrap.rds")

##### 4.5/ Plot completeness #####

### Load stuff

completness_raster_grid_within_range <- readRDS(file = "./outputs/Completeness_analyses/completness_raster_grid10_bootstrap.rds")

completness_raster_bioregions <- readRDS(file = "./outputs/Completeness_analyses/completness_raster_bioregions_bootstrap.rds")
bioregions_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/bioregions_shp_15min.rds")

# res <- "5"
res <- "15"

# Load mask for continent borders
continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_", res, ".rds"))
crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".rds"))

# Load shp for plot border, and grid
grid_Mollweide_out <- readRDS(file = "./input_data/Map_stuff/grid_Mollweide_out.rds")
large_bg_mask_Mollweide <- readRDS(file = "./input_data/Map_stuff/large_bg_mask_Mollweide.rds")
bbox_sp_Mollweide <- readRDS(file = "./input_data/Map_stuff/bbox_sp_Mollweide.rds")

load(file = "./input_data/Map_stuff/country_borders.RData")
country_borders <- as(country_borders, "Spatial")

library(tmaptools)
pal_grn_red_NA <- pal_grn_red <- rev(get_brewer_pal("RdYlGn", n = 200, plot = F))
pal_grn_red_NA[1] <- "#EDEDED"

### Project stuff

# Project raster into Mollweide
Mollweide_projection <- function(x) # Raster to project
{
  x_name <- deparse(substitute(x))  # Get the name of the initial raster as a character string
  
  new_map <- projectRaster(from = x, 
                           method = "ngb", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                           crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                           alignOnly = F)
  
  # Generate new object with "_Mollweide" suffix in the global environment
  eval(call("<<-", as.name(paste0(x_name, "_Mollweide")), new_map))
}

Mollweide_projection(x = completness_raster_grid_within_range)
Mollweide_projection(x = completness_raster_bioregions)

# Project sp shape into Mollweide
Mollweide_shp_projection <-  function(x) # Shp to project
{
  x_name <- deparse(substitute(x)) # Get the name of the initial shp as a character string
  
  new_shp <- spTransform(x, CRSobj = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
  
  # Generate new object with "_Mollweide" suffix in the global environment
  eval(call("<<-", as.name(paste0(x_name, "_Mollweide")), new_shp))
}

Mollweide_shp_projection(bioregions_shp)
Mollweide_shp_projection(country_borders)


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

### Make a 2 facets plot

pdf(file = paste0("./supplementaries/Sampling_completeness.pdf"), height = 5, width = 10)

internal_margins <- par()$mar
par(mar = c(3.1, 3.1, 2.7, 1.6))
par(mfrow = c(1,2))

# Panel A = Sampling completeness per 10° grid
map_indices_Mollweide(x = completness_raster_grid_within_range_Mollweide,
                      main_title = "Sampling completeness in 10° grid",
                      main_title_cex = 1.4, 
                      color_palette = pal_red_grn_NA,
                      axis_cex = 1.3,
                      legend_title = "Sampling\n              Completeness [%]",
                      legend_title_cex = 1.2,
                      legend_title_x = -3450,
                      legend_title_y = 880,
                      legend_breaks = seq(0, 100, 20))
text(x = 4000, y = -3900, labels = "(a)", font = 2, cex = 1.8)  # Add the contrast on the scale

# Panel B = Human Footprint without Light Pollution
map_indices_Mollweide(x = completness_raster_bioregions_Mollweide,
                      main_title = "Sampling completeness per Bioregions",
                      main_title_cex = 1.4, 
                      color_palette = pal_red_grn_NA,
                      axis_cex = 1.3,
                      legend_title = "Sampling\n              Completeness [%]",
                      legend_title_cex = 1.2,
                      legend_title_x = -3450,
                      legend_title_y = 880,
                      legend_breaks = seq(0, 100, 20))
text(x = 4000, y = -3900, labels = "(b)", font = 2, cex = 1.8)  # Add the contrast on the scale

par(mar = internal_margins)
par(mfrow = c(1,1))
dev.off()


##### 4.6/ Extract table #####


