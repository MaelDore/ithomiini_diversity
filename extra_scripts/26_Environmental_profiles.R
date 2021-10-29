
###### Script 26: Environmental variable profile of indices ######

# Plot the evolution of the altitudinal gradient of different biodiversity indices

# Only for Jaccard.80

### Use only 5 indices

# 1/ Species richness
# 2/ Mean Geographic Species Rarity = Range-size Weighted Species Richness
# 3/ Mean pairwise Phylogenetic Distance
# 4/ Mimicry richness
# 5/ Mean Geographic Mimicry Rarity = Range-size Weighted Mimicry Richness

###

###

# Inputs: 
   # Indices maps (4 main indices) from script 14a
   # Elevation raster

# Outputs:
   # Scatter plots of communities diversity indices ~ elevation

###


### 1/ Load stuff ####

# Effacer l'environnement
rm(list = ls())

library(raster)

# Load indices map
sp.richness <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
sp.rarity <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_Leroy_Jaccard.80.rds")
MPD <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80.rds")
ring.richness <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))
ring.rarity <- readRDS(file = paste0("./outputs/Indices_maps/ring.mean.rarity_Leroy_Jaccard.80.rds"))


# Create a mask of communities with presence of Ithomiini
Community_mask <- sp.richness
Community_mask[!(Community_mask > 0)] <- NA

plot(Community_mask)


# Load environmental layers
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_15.rds"))

# Mask all rasters to constrain to Ithomiini communities
envData_masked <- mask(envData, Community_mask)
sp.richness_masked <- mask(sp.richness, Community_mask)
sp.rarity_masked <- mask(sp.rarity, Community_mask)
MPD_masked <- mask(MPD, Community_mask)
ring.richness_masked <- mask(ring.richness, Community_mask)
ring.rarity_masked <- mask(ring.rarity, Community_mask)

# Quick scatter plots
pairs(stack(envData_masked, sp.richness_masked, sp.rarity_masked, MPD_masked, ring.richness_masked, ring.rarity_masked))

# Better scatter plots
envData_df <- as_tibble(envData@data@values)
names(envData_df)[1:4] <- c("Tmean", 'Tvar', "Htot", "Hvar")
envData_df$Tmean <- envData_df$Tmean / 10  # Tmean in °C
envData_df$Tvar <- envData_df$Tvar / 1000  # Tvar in °C
envData_df$Htot <- envData_df$Htot / 1000  # Htot in % water mass/air mass
envData_df$Hvar <- envData_df$Hvar / 1000  # Hvar in % of variation (CV * 100)

### 2/ Build df with data for plotting ####

indices_df <- tibble(sp.richness = sp.richness[], sp.rarity = sp.rarity[], MPD = MPD[], mimicry.richness = ring.richness[], mimicry.rarity = ring.rarity[])
env_profile_df <- cbind(envData_df, indices_df)

# Remove NA data
env_profile_df <- env_profile_df %>% 
  filter(!is.na(sp.richness)) %>% 
  filter(!is.nan(MPD))

# Save file with all environmental and indices data
saveRDS(env_profile_df, file = "./outputs/Env_profile/env_profile_df.rds", version = "2")

library(tidyverse)
library(ggthemes)
library(scales)
library(ggpubr)

# Build a dataframe with limits for axis ticks
limits_df <- tibble(x = c(0, 5000), y = c(0, 150))
min_vars <- c(5, rep(0, ncol(env_profile_df)-1))
max_vars <- c(35, 6, 2.5, 40, 5000, 80, 150, 0.6, 50, 30, 0.8)
breaks_vars <- c(5, 1, 0.5, 10, 1000, 20, 30, 0.1, 10, 5, 0.2)
limits_df <- as.data.frame(rbind(min_vars, max_vars, breaks_vars))
names(limits_df) <- names(env_profile_df)

# Record labels
var_labels <- c("Tmean\n[°C]", "Tvar\n[°C]", "Htot\n[% water mass/air mass]", "Hvar\n[% CV]", "Elevation\n[meters]", "Forest cover\n[%]", 
                "Species richness\n[Species]", "Species rarity\n[Rarity index]", "MPD\n[My]", "Mimicry richness\n[Mimicry rings]", "Mimicry rarity\n[Rarity index]")

### 3/ Plotting functions  ####

# Plotting function with all labels
scatter_plot_env_profile_all_labs <- function (i, j) 
{
  # Extract names of variables
  x <- sym(names(env_profile_df)[i])
  y <- sym(names(env_profile_df)[j])
  
  # Extract limits for axis ticks
  min_x <- limits_df[1, i]
  max_x <- limits_df[2, i]
  seq_x <- limits_df[3, i]
  
  min_y <- limits_df[1, j]
  max_y <- limits_df[2, j]
  seq_y <- limits_df[3, j]
  
  # Plot
  g1 <- ggplot(data = env_profile_df) +
    geom_point(mapping = aes(x = !!x, y = !!y), fill = "black", alpha = 0.1, show.legend = F) +
    geom_smooth(mapping = aes(x = !!x, y = !!y), 
                method = "gam", formula = y ~ s(x, bs = "cs"),
                color = "#FF0000FF", size = 1.5, se = F) +
    
    xlab(label = var_labels[i]) +
    ylab(label = var_labels[j]) +
    
    guides(alpha = F, color = T) +
    scale_x_continuous(breaks = seq(min_x, max_x, seq_x)) +
    scale_y_continuous(breaks = seq(min_y, max_y, seq_y)) +
    geom_rangeframe(data = limits_df, mapping = aes(x = !!x, y = !!y), size = 1.4) +
    
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
  
  return(g1)
}

# Plotting function with labels only for first column and last row
scatter_plot_env_profile_custom_labs <- function (i, j, last_index, first_env) 
{
  # Extract names of variables
  x <- sym(names(env_profile_df)[i])
  y <- sym(names(env_profile_df)[j])
  
  # Extract limits for axis ticks
  min_x <- limits_df[1, i]
  max_x <- limits_df[2, i]
  seq_x <- limits_df[3, i]
  
  min_y <- limits_df[1, j]
  max_y <- limits_df[2, j]
  seq_y <- limits_df[3, j]
  
  # Plot
  g1 <- ggplot(data = env_profile_df) +
    geom_point(mapping = aes(x = !!x, y = !!y), fill = "black", alpha = 0.1, show.legend = F) +
    geom_smooth(mapping = aes(x = !!x, y = !!y), 
                method = "gam", formula = y ~ s(x, bs = "cs"),
                color = "#FF0000FF", size = 1.5, se = F) +
    
    xlab(label = NULL) +
    ylab(label = NULL) +
    
    guides(alpha = F, color = T) +
    scale_x_continuous(breaks = seq(min_x, max_x, seq_x)) +
    scale_y_continuous(breaks = seq(min_y, max_y, seq_y)) +
    geom_rangeframe(data = limits_df, mapping = aes(x = !!x, y = !!y), size = 1.4) +
    
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
  
  # Add x-label only if last row
  if(j == last_index) 
    { g1 <- g1 + xlab(label = var_labels[i]) }
  
  # Add y-label only if first column
  if(i == first_env)
    { g1 <- g1 + ylab(label = var_labels[j]) }
  
  return(g1)
}

# Function to generate all graphs within lists
generate_plot_list <- function (Biodiv_indices, Env_indices) 
{
  plot_list_all_labels <- plot_list_custom_labels <- list()
  k <- 1
  for (j in Biodiv_indices) # By row = by Indices first
  {
    for (i in Env_indices) # By column = by Env. variables second
    {
      # eval(call("<<-", as.name(paste0("g", i, "_",j)), scatter_plot_env_profile(i, j))) # Store them in the environment
      
      # scatter_plot_env_profile(i, j) # Only plot them
      
      # Store them in a lists
      plot_list_all_labels[[k]] <<- scatter_plot_env_profile_all_labs(i, j) 
      plot_list_custom_labels[[k]] <<- scatter_plot_env_profile_custom_labs(i, j,
                                                                            last_index = last(Biodiv_indices),
                                                                            first_env = first(Env_indices))
      
      k <- k + 1
    }
  }
}

### 4/ Plot ####

# Choose environmental variables
Env_list <- c("Tmean", "Tvar", "Htot", "Hvar", "Elevation", "Forests")
Env_indices <- match(Env_list, names(env_profile_df))

# Choose biodiversity indices
Biodiv_list <- c("sp.richness", "sp.rarity", "MPD", "mimicry.richness", "mimicry.rarity")
Biodiv_indices <- match(Biodiv_list, names(env_profile_df))

# Generate all graphs within lists
generate_plot_list(Biodiv_indices, Env_indices)

# Multiple plots
?ggarrange

n_indices <- length(Biodiv_list)
n_vars <- length(Env_list)

# Plot with all labels and save all graphs on a single plot
pdf(file = paste0("./graphs/Env_profiles/all_variables_all_indices_all_labels.pdf"),
    height = 4 * n_indices, width = 5 * n_vars)
ggpubr::ggarrange(plotlist = plot_list_all_labels, nrow = n_indices, ncol = n_vars)
dev.off()

# Plot with custom labels and save all graphs on a single plot
pdf(file = paste0("./graphs/Env_profiles/all_variables_all_indices_custom_labels.pdf"),
    height = 4 * n_indices, width = 5 * n_vars)
ggpubr::ggarrange(plotlist = plot_list_custom_labels, nrow = n_indices, ncol = n_vars)
dev.off()




### Only for Elevation

# Load elevation raster
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_15.rds"))
DEM <- envData[["Elevation"]]

# Mask all rasters to constrain to Ithomiini communities
DEM_masked <- mask(DEM, Community_mask)

# Plot scatter plots
pairs(stack(DEM_masked, sp.richness_masked))
pairs(stack(DEM_masked, sp.rarity_masked))
pairs(stack(DEM_masked, MPD_masked))
pairs(stack(DEM_masked, ring.richness_masked))


### Boxplot elevation per ring

library(raster)

elevation <- readRDS(file = "./input_data/Env_data/Select_env_15.rds")[[5]]
ring_richness_brick <- brick(readRDS(file = "./outputs/Indices_stacks/All_ring_rich_stack_Jaccard.80.RData"))

elevation_data <- elevation@data@values
ring_richness_data <- ring_richness_brick@data@values

### Only weighted mean, no variance

weighted_mean_elevation_per_ring <- c()
for(i in 1:ncol(ring_richness_data))
{
  weighted_mean_elevation_per_ring[i] <- weighted.mean(x = elevation_data, w = ring_richness_data[, i], na.rm = T)
}
names(weighted_mean_elevation_per_ring) <- names(ring_richness_brick)

weighted_mean_elevation_per_ring_df <- data.frame("Elevation" = weighted_mean_elevation_per_ring, "Mimicry" = as.factor(names(ring_richness_brick)))
weighted_mean_elevation_per_ring_df <- weighted_mean_elevation_per_ring_df[order(weighted_mean_elevation_per_ring),]
weighted_mean_elevation_per_ring_df$Mimicry <- factor(x = weighted_mean_elevation_per_ring_df$Mimicry, levels = weighted_mean_elevation_per_ring_df$Mimicry)

boxplot(formula = Elevation ~ Mimicry, data = weighted_mean_elevation_per_ring_df)

View(weighted_mean_elevation_per_ring_df)

saveRDS(object = weighted_mean_elevation_per_ring_df, file = "./outputs/Env_profile/weighted_mean_elevation_per_ring_df.rds")
weighted_mean_elevation_per_ring_df <- readRDS(file = "./outputs/Env_profile/weighted_mean_elevation_per_ring_df.rds")



### Boxplot with variance but binarization and no weighting

# Extract elevation data only for richness > threshold

threshold <-  0.5
threshold <- 0.8

ring_richness_data_bin <- ring_richness_data > threshold

elevation_data_per_ring <- mimicry_data_per_ring <- elevation_data_per_ring_new <- mimicry_data_per_ring_new <- c()
for (i in 1:ncol(ring_richness_data))
{
  # i <- 1
  
  ring <- names(ring_richness_brick)[i]
  elevation_data_per_ring_new <- elevation_data[(ring_richness_data_bin[,i]) & !is.na(ring_richness_data_bin[,i])]
  mimicry_data_per_ring_new <- rep(x = ring, times = length(elevation_data_per_ring_new))
    
  elevation_data_per_ring <- c(elevation_data_per_ring, elevation_data_per_ring_new)  
  mimicry_data_per_ring <- c(mimicry_data_per_ring, mimicry_data_per_ring_new) 
}
elevation_per_ring_df <- data.frame("Elevation" = elevation_data_per_ring, "Mimicry" = as.factor(mimicry_data_per_ring))

median_elevation_per_ring <- as.numeric(by(data = elevation_data_per_ring, INDICES = as.factor(mimicry_data_per_ring), FUN = median))
names(median_elevation_per_ring) <- levels(as.factor(mimicry_data_per_ring))

median_elevation_per_ring_ordered <- median_elevation_per_ring[order(median_elevation_per_ring)]
View(median_elevation_per_ring_ordered)

elevation_per_ring_df$Mimicry <- factor(x = elevation_per_ring_df$Mimicry, levels = names(median_elevation_per_ring)[order(median_elevation_per_ring)])

boxplot(formula = Elevation ~ Mimicry, data = elevation_per_ring_df)

saveRDS(object = median_elevation_per_ring_ordered, file = "./outputs/Env_profile/median_elevation_per_ring_ordered.rds")
saveRDS(object = elevation_per_ring_df, file = "./outputs/Env_profile/elevation_per_ring_df.rds")

