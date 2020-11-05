
###### Script 19: Threat curves ######

# Compute evolution of mean threat levels among proportions of most diverse/rare communities

# Only for Jaccard.80

### Use only 4 indices

# 1/ Species richness
# 2/ Mean Continuous Species Rarity = Range-size Weighted Species Richness
# 3/ Mimicry richness
# 4/ Mean Continuous Mimicry rarity = Range-size Weighted Mimicry Richness
# Bonus/ Mean pairwise Phylogenetic Distance (removed from the final analyses)

###


# Inputs: 
   # Indices maps
   # Threat map: HII (HF, Venter et al., 2016)

# Outputs:
   # Evolution of mean, median, quantile 75%, 90% and 95% of threat levels (HII) among proportions of most rich/rare communities
   # Boxplot of HII levels for top 5%, 10% and 20% of most rich/rare communities
   # Evolution of mean, median, quantile 75%, 90% and 95% of threat levels (HII) among groups of most rich/rare communities based on normalized index thresholds (0.95, 0.90, 0.75)
   # Boxplot of HII levels among groups of most rich/rare communities based on normalized index thresholds (0.95, 0.90, 0.75)
   # Evolution of mean normalized indices among proportions of most threatened communities

###



### 1/ Load stuff ####


# Effacer l'environnement
rm(list = ls())

library(raster)
library(tidyverse)
library(ggplot2)
library(gridExtra) # To plot gglots side-by-side

# Load indices map
sp.richness <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
sp.rarity <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_Jaccard.80.rds")
MPD <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80.rds")
ring.richness <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))
ring.rarity <- readRDS(file = paste0("./outputs/Indices_maps/mimicry.mean.rarity_Jaccard.80.rds"))

HII <- readRDS(file = "./input_data/HII/HII.rds")

# Create a mask of communities with presence of Ithomiini
Community_mask <- sp.richness
Community_mask[!(Community_mask > 0)] <- NA

# Extract number of community with Ithomiini
N.com <- sum(!is.na(Community_mask[]))

# Different nb of community for MPD because cannot be computed from community predicted with less than 1 species
N.com.MPD <- length(MPD[MPD > 0])


### 2/ Most diverse, by rank = % of communities ####

### 2.1/ Extract values ####

# Extract indices of communities with Ithomiini ordered by decreasing value of HII
sp.richness_indices <- order(sp.richness[], decreasing = T)[1:N.com]
sp.rarity_indices <- order(sp.rarity[], decreasing = T)[1:N.com]
MPD_indices <- order(MPD[], decreasing = T)[1:N.com.MPD]
ring.richness_indices <- order(ring.richness[], decreasing = T)[1:N.com]
ring.rarity_indices <- order(ring.rarity[], decreasing = T)[1:N.com]

# Initiate vectors and lists to store outputs
sp.richness_mean_values <- sp.rarity_mean_values <- MPD_mean_values <- ring.richness_mean_values <- ring.rarity_mean_values <- integer()
sp.richness_median_values <- sp.rarity_median_values <- MPD_median_values <- ring.richness_median_values <- ring.rarity_median_values <- integer()
sp.richness_75_values <- sp.rarity_75_values <- MPD_75_values <- ring.richness_75_values <- ring.rarity_75_values <-integer()
sp.richness_90_values <- sp.rarity_90_values <- MPD_90_values <- ring.richness_90_values <- ring.rarity_90_values <- integer()
sp.richness_95_values <- sp.rarity_95_values <- MPD_95_values <- ring.richness_95_values <- ring.rarity_95_values <- integer()
sp.richness_boxplot_values <- sp.rarity_boxplot_values <- MPD_boxplot_values <- ring.richness_boxplot_values <- ring.rarity_boxplot_values <- list()
sp.richness_map_indices <- sp.rarity_map_indices <- MPD_map_indices <- ring.richness_map_indices <- ring.rarity_map_indices <- list()
j <- 1

# Loop
for (i in 1:1000)
{
  # i <- 1
  
  ## Species richness
  
  # Compute nb of communities to extract
  N.extract <- round(length(sp.richness_indices) * i/1000)
  # Get indices of communities to extract
  sp.richness_indices_extract <- sp.richness_indices[1:N.extract]
  # Extract values and compute mean
  sp.richness_values <- HII[sp.richness_indices_extract]
  sp.richness_mean_values[i] <- round(mean(sp.richness_values, na.rm = T), 3)
  sp.richness_median_values[i] <- round(median(sp.richness_values, na.rm = T), 3)
  sp.richness_75_values[i] <- round(quantile(x = sp.richness_values, probs = 0.75, na.rm = T), 3)
  sp.richness_90_values[i] <- round(quantile(x = sp.richness_values, probs = 0.90, na.rm = T), 3)
  sp.richness_95_values[i] <- round(quantile(x = sp.richness_values, probs = 0.95, na.rm = T), 3)
  
  ## Species rarity
  
  # Compute nb of communities to extract
  N.extract <- round(length(sp.rarity_indices) * i/1000)
  # Get indices of communities to extract
  sp.rarity_indices_extract <- sp.rarity_indices[1:N.extract]
  # Extract values and compute mean
  sp.rarity_values <- HII[sp.rarity_indices_extract]
  sp.rarity_mean_values[i] <- round(mean(sp.rarity_values, na.rm = T), 3)
  sp.rarity_median_values[i] <- round(median(sp.rarity_values, na.rm = T), 3)  
  sp.rarity_75_values[i] <- round(quantile(x = sp.rarity_values, probs = 0.75, na.rm = T), 3)
  sp.rarity_90_values[i] <- round(quantile(x = sp.rarity_values, probs = 0.90, na.rm = T), 3)
  sp.rarity_95_values[i] <- round(quantile(x = sp.rarity_values, probs = 0.95, na.rm = T), 3)
  
  ## MPD
  
  # Compute nb of communities to extract
  N.extract <- round(length(MPD_indices) * i/1000)
  # Get indices of communities to extract
  MPD_indices_extract <- MPD_indices[1:N.extract]
  # Extract values and compute mean
  MPD_values <- HII[MPD_indices_extract]
  MPD_mean_values[i] <- round(mean(MPD_values, na.rm = T), 3)
  MPD_median_values[i] <- round(median(MPD_values, na.rm = T), 3)
  MPD_75_values[i] <- round(quantile(x = MPD_values, probs = 0.75, na.rm = T), 3)
  MPD_90_values[i] <- round(quantile(x = MPD_values, probs = 0.90, na.rm = T), 3)
  MPD_95_values[i] <- round(quantile(x = MPD_values, probs = 0.95, na.rm = T), 3)
  
  ## Ring richness
  
  # Compute nb of communities to extract
  N.extract <- round(length(ring.richness_indices) * i/1000)
  # Get indices of communities to extract
  ring.richness_indices_extract <- ring.richness_indices[1:N.extract]
  # Extract values and compute mean
  ring.richness_values <- HII[ring.richness_indices_extract]
  ring.richness_mean_values[i] <- round(mean(ring.richness_values, na.rm = T), 3)
  ring.richness_median_values[i] <- round(median(ring.richness_values, na.rm = T), 3)
  ring.richness_75_values[i] <- round(quantile(x = ring.richness_values, probs = 0.75, na.rm = T), 3)
  ring.richness_90_values[i] <- round(quantile(x = ring.richness_values, probs = 0.90, na.rm = T), 3)
  ring.richness_95_values[i] <- round(quantile(x = ring.richness_values, probs = 0.95, na.rm = T), 3)

  ## Mimicry rarity
  
  # Compute nb of communities to extract
  N.extract <- round(length(ring.rarity_indices) * i/1000)
  # Get indices of communities to extract
  ring.rarity_indices_extract <- ring.rarity_indices[1:N.extract]
  # Extract values and compute mean
  ring.rarity_values <- HII[ring.rarity_indices_extract]
  ring.rarity_mean_values[i] <- round(mean(ring.rarity_values, na.rm = T), 3)
  ring.rarity_median_values[i] <- round(median(ring.rarity_values, na.rm = T), 3)  
  ring.rarity_75_values[i] <- round(quantile(x = ring.rarity_values, probs = 0.75, na.rm = T), 3)
  ring.rarity_90_values[i] <- round(quantile(x = ring.rarity_values, probs = 0.90, na.rm = T), 3)
  ring.rarity_95_values[i] <- round(quantile(x = ring.rarity_values, probs = 0.95, na.rm = T), 3) 
  
  # Store stuff for specific thresholds (5%, 10%, and 25%)
  
  if(i %in% c(50, 100, 250)) {
  
    # Store all values for specific threshold 
    sp.richness_boxplot_values[[j]] <- sp.richness_values
    sp.rarity_boxplot_values[[j]] <- sp.rarity_values
    MPD_boxplot_values[[j]] <- MPD_values
    ring.richness_boxplot_values[[j]] <- ring.richness_values
    ring.rarity_boxplot_values[[j]] <- ring.rarity_values
    
    # Store indices for specific threshold to be able to mask maps
    sp.richness_map_indices[[j]] <- sp.richness_indices_extract
    sp.rarity_map_indices[[j]] <- sp.rarity_indices_extract
    MPD_map_indices[[j]] <- MPD_indices_extract
    ring.richness_map_indices[[j]] <- ring.richness_indices_extract
    ring.rarity_map_indices[[j]] <- ring.rarity_indices_extract
    
    j <- j + 1
    
  }
  
  if(i %% 100 == 0) {print(i)}  
}

# Save data for boxplots
most_diverse_rank_global_boxplot_values <- save(sp.richness_boxplot_values, sp.rarity_boxplot_values, MPD_boxplot_values, ring.richness_boxplot_values, ring.rarity_boxplot_values, file = "./outputs/Threat_curves/Most_diverse/most_diverse_rank_global_boxplot_values.rds", version = "2")

# Save data for map masks
most_diverse_rank_global_map_indices <- save(sp.richness_map_indices, sp.rarity_map_indices, MPD_map_indices, ring.richness_map_indices, ring.rarity_map_indices, file = "./outputs/Threat_curves/Most_diverse/most_diverse_rank_global_map_indices.rds", version = "2")


### 2.2/ Curves plots for HII ~ ranks ####

### 2.2.1/ mean HII ####

# Store results for mean in a tidy df
most_diverse_mean_df <- tibble(sp.richness = sp.richness_mean_values, 
                               sp.rarity = sp.rarity_mean_values, 
                               # MPD = MPD_mean_values, 
                               ring.richness = ring.richness_mean_values,
                               ring.rarity = ring.rarity_mean_values)

most_diverse_mean_df <- most_diverse_mean_df %>% 
  mutate(prop_com = (1:1000)/10) %>%    # Add a column for % of communities
  mutate(stat = "mean") %>%             # Add a column for stat computed
  filter(prop_com > 2) %>%              # Remove value for the top 2% because curves are idiosyncratic (nb of communities too low)
  pivot_longer(# cols = c("sp.richness", "sp.rarity", "MPD", "ring.richness"), # List of the name of the "columns" that are actually factor levels
               cols = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"), # List of the name of the "columns" that are actually factor levels
               names_to = "index",   # name of the new column that is the real variable for these factor level
               values_to = "HII", # name of the new columns that will get the values in the cells of the previous dataset
               values_drop_na = T) %>% # Drop rows for missing values ?
  mutate(index = factor(index, 
                        # levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness")
                        levels = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))) # Change factor level order

# Save results
saveRDS(most_diverse_mean_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_mean_df.rds", version = "2")

# Load results
most_diverse_mean_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_mean_df.rds")

# Plot resulting curves for mean HII
pdf(file = paste0("./graphs/Threat_curves/Most_diverse/Most_diverse_rank_global_mean.pdf"), height = 5.3, width = 6.5)
g <- ggplot(most_diverse_mean_df, aes(x = prop_com, y = HII)) +
  geom_line(aes(color = index), lwd = 2) +
  geom_hline(yintercept = mean(most_diverse_mean_df[most_diverse_mean_df$prop_com == 100, "HII", drop = T]), # Get the final mean HII for 100% of communities
             lwd = 1.2) +
  geom_vline(xintercept = c(5, 10, 25), # Draw the threshold for the boxplots
             col = "grey20", lty = 3, lwd = 0.3) +
  expand_limits(y = 0, x = 0) + # To tell to always keep those values inside the plot. # Set the origin for the graph
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, .05))) + # To tell not to add padding around the limits of the axis
  labs(y = "mean HII", x = "Proportion of communities",
       title = "Full map ; Mean HII") +
  theme_bw() +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness"))
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity"))
print(g)
dev.off()

### 2.1.2/ median HII ####

# Store results for mean in a tidy df
most_diverse_median_df <- tibble(sp.richness = sp.richness_median_values,
                                 sp.rarity = sp.rarity_median_values,
                                 # MPD = MPD_median_values,
                                 ring.richness = ring.richness_median_values,
                                 ring.rarity = ring.rarity_median_values)

most_diverse_median_df <- most_diverse_median_df %>% 
  mutate(prop_com = (1:1000)/10) %>%    # Add a column for % of communities
  mutate(stat = "median") %>%             # Add a column for stat computed
  filter(prop_com > 2) %>%              # Remove value for the top 2% because curves are idiosyncratic (nb of communities too low)
  pivot_longer(# cols = c("sp.richness", "sp.rarity", "MPD", "ring.richness"), # List of the name of the "columns" that are actually factor levels
               cols = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"), # List of the name of the "columns" that are actually factor levels
               names_to = "index",   # name of the new column that is the real variable for these factor level
               values_to = "HII", # name of the new columns that will get the values in the cells of the previous dataset
               values_drop_na = T) %>% # Drop rows for missing values ?
  mutate(index = factor(index,
                        # levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness")
                        levels = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))) # Change factor level order

# Save results
saveRDS(most_diverse_median_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_median_df.rds", version = "2")

# Load results
most_diverse_median_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_median_df.rds")

# Plot resulting curves for median HII
pdf(file = paste0("./graphs/Threat_curves/Most_diverse/Most_diverse_rank_global_median.pdf"), height = 5.3, width = 6.5)
g <- ggplot(most_diverse_median_df, aes(x = prop_com, y = HII)) +
  geom_line(aes(color = index), lwd = 2) +
  geom_hline(yintercept = median(most_diverse_median_df[most_diverse_median_df$prop_com == 100, "HII", drop = T]), # Get the final median HII for 100% of communities
             lwd = 1.2) +
  geom_vline(xintercept = c(5, 10, 25), # Draw the threshold for the boxplots
             col = "grey20", lty = 3, lwd = 0.3) +
  expand_limits(y = 0, x = 0) + # To tell to always keep those values inside the plot. # Set the origin for the graph
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, .05))) + # To tell not to add padding around the limits of the axis
  labs(y = "median HII", x = "Proportion of communities",
       title = "Full map ; Median HII") +
  theme_bw() +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness"))
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity"))
print(g)
dev.off()

### 2.1.3/ HII 75% ####

# Store results for mean in a tidy df
most_diverse_75_df <- tibble(sp.richness = sp.richness_75_values,
                             sp.rarity = sp.rarity_75_values,
                             # MPD = MPD_75_values,
                             ring.richness = ring.richness_75_values,
                             ring.rarity = ring.rarity_75_values)

most_diverse_75_df <- most_diverse_75_df %>% 
  mutate(prop_com = (1:1000)/10) %>%    # Add a column for % of communities
  mutate(stat = "75") %>%             # Add a column for stat computed
  filter(prop_com > 2) %>%              # Remove value for the top 2% because curves are idiosyncratic (nb of communities too low)
  pivot_longer(# cols = c("sp.richness", "sp.rarity", "MPD", "ring.richness"), # List of the name of the "columns" that are actually factor levels
               cols = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"), # List of the name of the "columns" that are actually factor levels
               names_to = "index",   # name of the new column that is the real variable for these factor level
               values_to = "HII", # name of the new columns that will get the values in the cells of the previous dataset
               values_drop_na = T) %>% # Drop rows for missing values ?
  mutate(index = factor(index,
                        # levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness")
                        levels = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))) # Change factor level order

# Save results
saveRDS(most_diverse_75_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_75_df.rds", version = "2")

# Load results
most_diverse_75_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_75_df.rds")

# Plot resulting curves for HII 75%
pdf(file = paste0("./graphs/Threat_curves/Most_diverse/Most_diverse_rank_global_75.pdf"), height = 5.3, width = 6.5)
g <- ggplot(most_diverse_75_df, aes(x = prop_com, y = HII)) +
  geom_line(aes(color = index), lwd = 2) +
  geom_hline(yintercept = quantile(x = most_diverse_75_df[most_diverse_75_df$prop_com == 100, "HII", drop = T], probs = 0.75), # Get the final quantile 75% HII for 100% of communities
             lwd = 1.2) +
  geom_vline(xintercept = c(5, 10, 25), # Draw the threshold for the boxplots
             col = "grey20", lty = 3, lwd = 0.3) +
  expand_limits(y = 0, x = 0) + # To tell to always keep those values inside the plot. # Set the origin for the graph
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, .05))) + # To tell not to add padding around the limits of the axis
  labs(y = "HII (quantile 0.75)", x = "Proportion of communities",
       title = "Full map ; HII 75%") +
  theme_bw() +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness"))
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity"))
print(g)
dev.off()

### 2.1.4/ HII 90% ####

# Store results for mean in a tidy df
most_diverse_90_df <- tibble(sp.richness = sp.richness_90_values,
                             sp.rarity = sp.rarity_90_values,
                             # MPD = MPD_90_values,
                             ring.richness = ring.richness_90_values,
                             ring.rarity = ring.rarity_90_values)

most_diverse_90_df <- most_diverse_90_df %>% 
  mutate(prop_com = (1:1000)/10) %>%    # Add a column for % of communities
  mutate(stat = "90") %>%             # Add a column for stat computed
  filter(prop_com > 2) %>%              # Remove value for the top 2% because curves are idiosyncratic (nb of communities too low)
  pivot_longer(# cols = c("sp.richness", "sp.rarity", "MPD", "ring.richness"), # List of the name of the "columns" that are actually factor levels
               cols = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"), # List of the name of the "columns" that are actually factor levels
               names_to = "index",   # name of the new column that is the real variable for these factor level
               values_to = "HII", # name of the new columns that will get the values in the cells of the previous dataset
               values_drop_na = T) %>% # Drop rows for missing values ?
  mutate(index = factor(index,
                        # levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness")
                        levels = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))) # Change factor level order

# Save results
saveRDS(most_diverse_90_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_90_df.rds", version = "2")

# Load results
most_diverse_90_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_90_df.rds")

# Plot resulting curves for HII 90%
pdf(file = paste0("./graphs/Threat_curves/Most_diverse/Most_diverse_rank_global_90.pdf"), height = 5.3, width = 6.5)
g <- ggplot(most_diverse_90_df, aes(x = prop_com, y = HII)) +
  geom_line(aes(color = index), lwd = 2) +
  geom_hline(yintercept = quantile(x = most_diverse_90_df[most_diverse_90_df$prop_com == 100, "HII", drop = T], probs = 0.90), # Get the final quantile 90% HII for 100% of communities
             lwd = 1.2) +
  geom_vline(xintercept = c(5, 10, 25), # Draw the threshold for the boxplots
             col = "grey20", lty = 3, lwd = 0.3) +
  expand_limits(y = 0, x = 0) + # To tell to always keep those values inside the plot. # Set the origin for the graph
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, .05))) + # To tell not to add padding around the limits of the axis
  labs(y = "HII (quantile 0.90)", x = "Proportion of communities",
       title = "Full map ; HII 90%") +
  theme_bw() +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness"))
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity"))
print(g)
dev.off()


### 2.1.5/ HII 95% ####

# Store results for mean in a tidy df
most_diverse_95_df <- tibble(sp.richness = sp.richness_95_values,
                             sp.rarity = sp.rarity_95_values,
                             # MPD = MPD_95_values,
                             ring.richness = ring.richness_95_values,
                             ring.rarity = ring.rarity_95_values)

most_diverse_95_df <- most_diverse_95_df %>% 
  mutate(prop_com = (1:1000)/10) %>%    # Add a column for % of communities
  mutate(stat = "95") %>%             # Add a column for stat computed
  filter(prop_com > 2) %>%              # Remove value for the top 2% because curves are idiosyncratic (nb of communities too low)
  pivot_longer(# cols = c("sp.richness", "sp.rarity", "MPD", "ring.richness"), # List of the name of the "columns" that are actually factor levels
               cols = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"), # List of the name of the "columns" that are actually factor levels
               names_to = "index",   # name of the new column that is the real variable for these factor level
               values_to = "HII", # name of the new columns that will get the values in the cells of the previous dataset
               values_drop_na = T) %>% # Drop rows for missing values ?
  mutate(index = factor(index,
                        # levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness")
                        levels = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))) # Change factor level order

# Save results
saveRDS(most_diverse_95_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_95_df.rds", version = "2")

# Load results
most_diverse_95_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_95_df.rds")

# Plot resulting curves for HII 95%
pdf(file = paste0("./graphs/Threat_curves/Most_diverse/Most_diverse_rank_global_95.pdf"), height = 5.3, width = 6.5)
g <- ggplot(most_diverse_95_df, aes(x = prop_com, y = HII)) +
  geom_line(aes(color = index), lwd = 2) +
  geom_hline(yintercept = quantile(x = most_diverse_95_df[most_diverse_95_df$prop_com == 100, "HII", drop = T], probs = 0.95), # Get the final quantile 95% HII for 100% of communities
             lwd = 1.2) +
  geom_vline(xintercept = c(5, 10, 25), # Draw the threshold for the boxplots
             col = "grey20", lty = 3, lwd = 0.3) +
  expand_limits(y = 0, x = 0) + # To tell to always keep those values inside the plot. # Set the origin for the graph
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, .05))) + # To tell not to add padding around the limits of the axis
  labs(y = "HII (quantile 0.95)", x = "Proportion of communities",
       title = "Full map ; HII 95%") +
  theme_bw() +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness"))
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity"))
print(g)
dev.off()

### 2.1.6/ all HII quantiles together (50%, 75%, 90%, 95%) ####

# Store results for mean in a tidy df
most_diverse_all_quantiles_df <- rbind(most_diverse_median_df, most_diverse_75_df, most_diverse_90_df, most_diverse_95_df)

# Save results
saveRDS(most_diverse_all_quantiles_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_all_quantiles_df.rds", version = "2")

# Load results
most_diverse_all_quantiles_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_all_quantiles_df.rds")

# Extract only median and quantile 95%
most_diverse_median_95_df <- most_diverse_all_quantiles_df %>% 
  filter(stat %in% c("median","95"))

saveRDS(most_diverse_median_95_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_median_95_df.rds", version = "2")

# Plot resulting curves for HII quantiles 50% and 95%
pdf(file = paste0("./graphs/Threat_curves/Most_diverse/Most_diverse_rank_global_median_95.pdf"), height = 5.3, width = 6.5)
g <- ggplot(most_diverse_median_95_df, aes(x = prop_com, y = HII)) +
  geom_line(aes(color = index, lty = stat), lwd = 2) +
  geom_hline(yintercept = quantile(x = most_diverse_median_95_df[(most_diverse_median_95_df$prop_com == 100), "HII", drop = T], probs = 0.95), # Get the final quantile 95% HII for 100% of communities
             lwd = 1.2) +
  geom_hline(yintercept = median(x = most_diverse_median_95_df[(most_diverse_median_95_df$prop_com == 100) & (most_diverse_median_95_df$stat == "median"), "HII", drop = T]), # Get the final median HII for 100% of communities
             lwd = 1.2) +
  geom_vline(xintercept = c(5, 10, 25), # Draw the threshold for the boxplots
             col = "grey20", lty = 3, lwd = 0.3) +
  expand_limits(y = 0, x = 0) + # To tell to always keep those values inside the plot. # Set the origin for the graph
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, .05))) + # To tell not to add padding around the limits of the axis
  labs(y = "HII (median & quantile 0.95)", x = "Proportion of communities",
       title = "Full map ; HII median and 95%") +
  theme_bw() +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness"))
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity"))
print(g)
dev.off()

### 2.3/ Boxplots for specific % of communities ####

# Prepare data

# Load data for boxplots
most_diverse_rank_global_boxplot_values <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_rank_global_boxplot_values.rds")

most_diverse_rank_global_boxplot_lists <- c(sp.richness_boxplot_values, 
                                            sp.rarity_boxplot_values, 
                                            # MPD_boxplot_values, 
                                            ring.richness_boxplot_values,
                                            ring.rarity_boxplot_values)

# Counts nb of observations recorded for each threshold and each index
n_com <- unlist(map(.x = most_diverse_rank_global_boxplot_lists, .f = length))

# Build the final tidy df
most_diverse_rank_global_boxplot_df <- tibble(HII = unlist(most_diverse_rank_global_boxplot_lists)) %>% 
  mutate(threshold = c(rep("5", n_com[1]),
                      rep("10", n_com[2]),
                      rep("25", n_com[3]),
                      rep("5", n_com[4]),
                      rep("10", n_com[5]),
                      rep("25", n_com[6]),
                      rep("5", n_com[7]),
                      rep("10", n_com[8]),
                      rep("25", n_com[9]),
                      rep("5", n_com[10]),
                      rep("10", n_com[11]),
                      rep("25", n_com[12]))) %>% 
  mutate(index = c(rep("sp.richness", sum(n_com[1:3])),
                   rep("sp.rarity", sum(n_com[4:6])),
                   # rep("MPD", sum(n_com[7:9])),
                   # rep("ring.richness", sum(n_com[10:12])),
                   rep("ring.richness", sum(n_com[7:9])),
                   rep("ring.rarity", sum(n_com[10:12])))) %>% 
  mutate(index = factor(index,
                        # levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness")
                        levels = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))) # Change factor level order

saveRDS(most_diverse_rank_global_boxplot_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_rank_global_boxplot_values.rds", version = "2")

### 2.3.1/ Boxplot for 5% of most diverse communities ####

# Load data
most_diverse_rank_global_boxplot_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_rank_global_boxplot_values.rds")

# Filter data
most_diverse_rank_global_boxplot_5_df <- most_diverse_rank_global_boxplot_df %>% 
  filter(threshold == "5")

# Compute quantiles 95 to add to the plot
quant_95 <- as.numeric(with(most_diverse_rank_global_boxplot_5_df, by(HII, INDICES = index, FUN = quantile, probs = 0.95)))
quant_95 <- tibble(quant_value = quant_95, 
                   # index = c("sp.richness", "sp.rarity", "MPD", "ring.richness")
                   index = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity"))
most_diverse_rank_global_boxplot_5_df <- most_diverse_rank_global_boxplot_5_df %>% 
  left_join(y = quant_95, by = "index") %>% 
  mutate(index = factor(index,
                        # levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness")
                        levels = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))) # Change factor level order

saveRDS(most_diverse_rank_global_boxplot_5_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_rank_global_boxplot_5_df.rds", version = "2")

# Plot the boxplot
pdf(file = paste0("./graphs/Threat_curves/Most_diverse/Most_diverse_rank_global_boxplot_5.pdf"), height = 5.3, width = 6.5)
g <- ggplot(most_diverse_rank_global_boxplot_5_df, aes(x = index, y = HII), show.legend = F) +
  geom_boxplot(aes(fill = index)) +
  geom_point(aes(x = index, y = quant_value, col = index), size = 5) +
  labs(y = "HII", x = "Index",
  title = "Full map ; 5% most rich/diverse/rare") +
  theme_bw() +
  # scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                    name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  
  guides(color = "none")
print(g)
dev.off()

### 2.3.2/ Boxplot for 10% of most diverse communities ####

# Load data
most_diverse_rank_global_boxplot_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_rank_global_boxplot_values.rds")

# Filter data
most_diverse_rank_global_boxplot_10_df <- most_diverse_rank_global_boxplot_df %>% 
  filter(threshold == "10")

# Compute quantiles 95 to add to the plot
quant_95 <- as.numeric(with(most_diverse_rank_global_boxplot_10_df, by(HII, INDICES = index, FUN = quantile, probs = 0.95)))
quant_95 <- tibble(quant_value = quant_95, 
                   # index = c("sp.richness", "sp.rarity", "MPD", "ring.richness")
                   index = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity"))
most_diverse_rank_global_boxplot_10_df <- most_diverse_rank_global_boxplot_10_df %>% 
  left_join(y = quant_95, by = "index") %>% 
  mutate(index = factor(index,
                        # levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness")
                        levels = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))) # Change factor level order

saveRDS(most_diverse_rank_global_boxplot_10_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_rank_global_boxplot_10_df.rds", version = "2")


# Plot the boxplot
pdf(file = paste0("./graphs/Threat_curves/Most_diverse/Most_diverse_rank_global_boxplot_10.pdf"), height = 5.3, width = 6.5)
g <- ggplot(most_diverse_rank_global_boxplot_10_df, aes(x = index, y = HII), show.legend = F) +
  geom_boxplot(aes(fill = index)) +
  geom_point(aes(x = index, y = quant_value, col = index), size = 5) +
  labs(y = "HII", x = "Index",
       title = "Full map ; 10% most rich/diverse/rare") +
  theme_bw() +
  # scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                    name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  guides(color = "none")
print(g)
dev.off()

### 2.3.3/ Boxplot for 25% of most diverse communities ####

# Load data
most_diverse_rank_global_boxplot_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_rank_global_boxplot_values.rds")

# Filter data
most_diverse_rank_global_boxplot_25_df <- most_diverse_rank_global_boxplot_df %>% 
  filter(threshold == "25")

# Compute quantiles 95 to add to the plot
quant_95 <- as.numeric(with(most_diverse_rank_global_boxplot_25_df, by(HII, INDICES = index, FUN = quantile, probs = 0.95)))
quant_95 <- tibble(quant_value = quant_95, 
                   # index = c("sp.richness", "sp.rarity", "MPD", "ring.richness")
                   index = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity"))
most_diverse_rank_global_boxplot_25_df <- most_diverse_rank_global_boxplot_25_df %>% 
  left_join(y = quant_95, by = "index") %>% 
  mutate(index = factor(index,
                        # levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness")
                        levels = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))) # Change factor level order

saveRDS(most_diverse_rank_global_boxplot_25_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_rank_global_boxplot_25_df.rds", version = "2")

# Plot the boxplot
pdf(file = paste0("./graphs/Threat_curves/Most_diverse/Most_diverse_rank_global_boxplot_25.pdf"), height = 5.3, width = 6.5)
g <- ggplot(most_diverse_rank_global_boxplot_25_df, aes(x = index, y = HII), show.legend = F) +
  geom_boxplot(aes(fill = index)) +
  geom_point(aes(x = index, y = quant_value, col = index), size = 5) +
  labs(y = "HII", x = "Index",
       title = "Full map ; 25% most rich/diverse/rare") +
  theme_bw() +
  # scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                    name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  guides(color = "none")
print(g)
dev.off()


### 2.4/ Multiple plot: curves and boxplots ####

most_diverse_median_95_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_median_95_df.rds")
most_diverse_rank_global_boxplot_5_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_rank_global_boxplot_5_df.rds")
most_diverse_rank_global_boxplot_10_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_rank_global_boxplot_10_df.rds")
most_diverse_rank_global_boxplot_25_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_rank_global_boxplot_25_df.rds")

# Plot the boxplot
pdf(file = paste0("./graphs/Threat_curves/Most_diverse/Most_diverse_rank_global_curves_boxplot.pdf"), height = 11, width = 13)

g1 <- ggplot(most_diverse_median_95_df, aes(x = prop_com, y = HII)) +
  geom_line(aes(color = index, lty = stat), lwd = 2) +
  geom_hline(yintercept = quantile(x = most_diverse_median_95_df[(most_diverse_median_95_df$prop_com == 100), "HII", drop = T], probs = 0.95), # Get the final quantile 95% HII for 100% of communities
             lwd = 1.2) +
  geom_hline(yintercept = median(x = most_diverse_median_95_df[(most_diverse_median_95_df$prop_com == 100) & (most_diverse_median_95_df$stat == "median"), "HII", drop = T]), # Get the final median HII for 100% of communities
             lwd = 1.2) +
  geom_vline(xintercept = c(5, 10, 25), # Draw the threshold for the boxplots
             col = "grey20", lty = 3, lwd = 0.3) +
  expand_limits(y = 0, x = 0) + # To tell to always keep those values inside the plot. # Set the origin for the graph
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, .05))) + # To tell not to add padding around the limits of the axis
  labs(y = "HII (median & quantile 0.95)", x = "Proportion of communities",
       title = "Full map ; HII median and 95%") +
  theme_bw() +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness"))
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity"))
  
g2 <- ggplot(most_diverse_rank_global_boxplot_5_df, aes(x = index, y = HII), show.legend = F) +
  geom_boxplot(aes(fill = index)) +
  geom_point(aes(x = index, y = quant_value, col = index), size = 5) +
  labs(y = "HII", x = "Index",
       title = "Full map ; 5% most rich/diverse/rare") +
  theme_bw() +
  # scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                    name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  guides(color = "none")

g3 <- ggplot(most_diverse_rank_global_boxplot_10_df, aes(x = index, y = HII), show.legend = F) +
  geom_boxplot(aes(fill = index)) +
  geom_point(aes(x = index, y = quant_value, col = index), size = 5) +
  labs(y = "HII", x = "Index",
       title = "Full map ; 10% most rich/diverse/rare") +
  theme_bw() +
  # scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                    name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  guides(color = "none")

g4 <- ggplot(most_diverse_rank_global_boxplot_25_df, aes(x = index, y = HII), show.legend = F) +
  geom_boxplot(aes(fill = index)) +
  geom_point(aes(x = index, y = quant_value, col = index), size = 5) +
  labs(y = "HII", x = "Index",
       title = "Full map ; 25% most rich/diverse/rare") +
  theme_bw() +
  # scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                    name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  guides(color = "none")

gridExtra::grid.arrange(g1, g2, g3, g4, ncol=2)


dev.off()




### 3/ Most diverse, by values. = threshold of normalized indices ####

### 3.1/ Extract values ####

# Load normalized indices map
sp.richness_norm <- readRDS(file = "./outputs/Threat_maps/sp.richness_norm.rds")
sp.rarity_norm <- readRDS(file = "./outputs/Threat_maps/sp.rarity_norm.rds")
MPD_norm <- readRDS(file = "./outputs/Threat_maps/MPD_norm.rds")
ring.richness_norm <- readRDS(file = "./outputs/Threat_maps/ring.richness_norm.rds")
ring.rarity_norm <- readRDS(file = "./outputs/Threat_maps/ring.rarity_norm.rds")

# Initiate vectors and lists to store outputs
sp.richness_mean_values <- sp.rarity_mean_values <- MPD_mean_values <- ring.richness_mean_values <- ring.rarity_mean_values <- integer()
sp.richness_median_values <- sp.rarity_median_values <- MPD_median_values <- ring.richness_median_values <- ring.rarity_median_values <- integer()
sp.richness_75_values <- sp.rarity_75_values <- MPD_75_values <- ring.richness_75_values <- ring.rarity_75_values <- integer()
sp.richness_90_values <- sp.rarity_90_values <- MPD_90_values <- ring.richness_90_values <- ring.rarity_90_values <- integer()
sp.richness_95_values <- sp.rarity_95_values <- MPD_95_values <- ring.richness_95_values <- ring.rarity_95_values <- integer()
sp.richness_boxplot_values <- sp.rarity_boxplot_values <- MPD_boxplot_values <- ring.richness_boxplot_values <- ring.rarity_boxplot_values <- list()
sp.richness_map_indices <- sp.rarity_map_indices <- MPD_map_indices <- ring.richness_map_indices <- ring.rarity_map_indices <- list()
j <- 1

# Loop
for (i in 1:1000)
{
  # i <- 1
  
  # Limit values to include a community in the pool
  threshold <- 1-(i/1000)
  
  ## Species richness
  
  # Get indices of communities to extract
  sp.richness_indices_extract <- which(sp.richness_norm[] > threshold)
  # Extract values and compute mean
  sp.richness_values <- HII[sp.richness_indices_extract]
  sp.richness_mean_values[i] <- round(mean(sp.richness_values, na.rm = T), 3)
  sp.richness_median_values[i] <- round(median(sp.richness_values, na.rm = T), 3)
  sp.richness_75_values[i] <- round(quantile(x = sp.richness_values, probs = 0.75, na.rm = T), 3)
  sp.richness_90_values[i] <- round(quantile(x = sp.richness_values, probs = 0.90, na.rm = T), 3)
  sp.richness_95_values[i] <- round(quantile(x = sp.richness_values, probs = 0.95, na.rm = T), 3)
  
  ## Species rarity
  
  # Get indices of communities to extract
  sp.rarity_indices_extract <- which(sp.rarity_norm[] > threshold)
  # Extract values and compute mean
  sp.rarity_values <- HII[sp.rarity_indices_extract]
  sp.rarity_mean_values[i] <- round(mean(sp.rarity_values, na.rm = T), 3)
  sp.rarity_median_values[i] <- round(median(sp.rarity_values, na.rm = T), 3)  
  sp.rarity_75_values[i] <- round(quantile(x = sp.rarity_values, probs = 0.75, na.rm = T), 3)
  sp.rarity_90_values[i] <- round(quantile(x = sp.rarity_values, probs = 0.90, na.rm = T), 3)
  sp.rarity_95_values[i] <- round(quantile(x = sp.rarity_values, probs = 0.95, na.rm = T), 3)
  
  ## MPD
  
  # Get indices of communities to extract
  MPD_indices_extract <- which(MPD_norm[] > threshold)
  # Extract values and compute mean
  MPD_values <- HII[MPD_indices_extract]
  MPD_mean_values[i] <- round(mean(MPD_values, na.rm = T), 3)
  MPD_median_values[i] <- round(median(MPD_values, na.rm = T), 3)
  MPD_75_values[i] <- round(quantile(x = MPD_values, probs = 0.75, na.rm = T), 3)
  MPD_90_values[i] <- round(quantile(x = MPD_values, probs = 0.90, na.rm = T), 3)
  MPD_95_values[i] <- round(quantile(x = MPD_values, probs = 0.95, na.rm = T), 3)
  
  ## Ring richness
  
  # Get indices of communities to extract
  ring.richness_indices_extract <- which(ring.richness_norm[] > threshold)
  # Extract values and compute mean
  ring.richness_values <- HII[ring.richness_indices_extract]
  ring.richness_mean_values[i] <- round(mean(ring.richness_values, na.rm = T), 3)
  ring.richness_median_values[i] <- round(median(ring.richness_values, na.rm = T), 3)
  ring.richness_75_values[i] <- round(quantile(x = ring.richness_values, probs = 0.75, na.rm = T), 3)
  ring.richness_90_values[i] <- round(quantile(x = ring.richness_values, probs = 0.90, na.rm = T), 3)
  ring.richness_95_values[i] <- round(quantile(x = ring.richness_values, probs = 0.95, na.rm = T), 3)
  
  ## Miicry rarity
  
  # Get indices of communities to extract
  ring.rarity_indices_extract <- which(ring.rarity_norm[] > threshold)
  # Extract values and compute mean
  ring.rarity_values <- HII[ring.rarity_indices_extract]
  ring.rarity_mean_values[i] <- round(mean(ring.rarity_values, na.rm = T), 3)
  ring.rarity_median_values[i] <- round(median(ring.rarity_values, na.rm = T), 3)  
  ring.rarity_75_values[i] <- round(quantile(x = ring.rarity_values, probs = 0.75, na.rm = T), 3)
  ring.rarity_90_values[i] <- round(quantile(x = ring.rarity_values, probs = 0.90, na.rm = T), 3)
  ring.rarity_95_values[i] <- round(quantile(x = ring.rarity_values, probs = 0.95, na.rm = T), 3)
  
  
  # Store stuff for specific thresholds (5%, 10%, and 25%)
  
  if(i %in% c(50, 100, 250)) {
    
    # Store all values for specific threshold 
    sp.richness_boxplot_values[[j]] <- sp.richness_values
    sp.rarity_boxplot_values[[j]] <- sp.rarity_values
    MPD_boxplot_values[[j]] <- MPD_values
    ring.richness_boxplot_values[[j]] <- ring.richness_values
    ring.rarity_boxplot_values[[j]] <- ring.rarity_values
    
    # Store indices for specific threshold to be able to mask maps
    sp.richness_map_indices[[j]] <- sp.richness_indices_extract
    sp.rarity_map_indices[[j]] <- sp.rarity_indices_extract
    MPD_map_indices[[j]] <- MPD_indices_extract
    ring.richness_map_indices[[j]] <- ring.richness_indices_extract
    ring.rarity_map_indices[[j]] <- ring.rarity_indices_extract
    
    j <- j + 1
    
  }
  
  if(i %% 100 == 0) {print(i)}  
}

# Save data for boxplots
most_diverse_norm_global_boxplot_values <- save(sp.richness_boxplot_values, sp.rarity_boxplot_values, MPD_boxplot_values, ring.richness_boxplot_values, ring.rarity_boxplot_values, file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_global_boxplot_values.rds", version = "2")

# Save data for map masks
most_diverse_norm_global_map_indices <- save(sp.richness_map_indices, sp.rarity_map_indices, MPD_map_indices, ring.richness_map_indices, ring.rarity_map_indices, file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_global_map_indices.rds", version = "2")


### 3.2/ Curves plots for HII ~ norm ####

### 3.2.1/ mean HII ####

# Store results for mean in a tidy df
most_diverse_norm_mean_df <- tibble(sp.richness = sp.richness_mean_values, 
                                    sp.rarity = sp.rarity_mean_values, 
                                    # MPD = MPD_mean_values, 
                                    ring.richness = ring.richness_mean_values,
                                    ring.rarity = ring.rarity_mean_values)

most_diverse_norm_mean_df <- most_diverse_norm_mean_df %>% 
  mutate(threshold = 1-(1:1000)/1000) %>%    # Add a column for threhold norm values
  mutate(stat = "mean") %>%             # Add a column for stat computed
  filter(threshold <= 1) %>%              # Remove value for the top thesholds because curves are idiosyncratic (nb of communities too low)
  pivot_longer(# cols = c("sp.richness", "sp.rarity", "MPD", "ring.richness"), # List of the name of the "columns" that are actually factor levels
               cols = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"), # List of the name of the "columns" that are actually factor levels
               names_to = "index",   # name of the new column that is the real variable for these factor level
               values_to = "HII", # name of the new columns that will get the values in the cells of the previous dataset
               values_drop_na = T) %>% # Drop rows for missing values ?
  mutate(index = factor(index, 
                        # levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness"))) # Change factor level order
                        levels = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))) # Change factor level order
# Save results
saveRDS(most_diverse_norm_mean_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_mean_df.rds", version = "2")

# Load results
most_diverse_norm_mean_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_mean_df.rds")

# Plot resulting curves for mean HII
pdf(file = paste0("./graphs/Threat_curves/Most_diverse/Most_diverse_norm_global_mean.pdf"), height = 5.3, width = 6.5)
g <- ggplot(most_diverse_norm_mean_df, aes(x = threshold, y = HII)) +
  geom_line(aes(color = index), lwd = 2) +
  geom_hline(yintercept = mean(most_diverse_norm_mean_df[most_diverse_norm_mean_df$threshold == 0, "HII", drop = T]), # Get the final mean HII for 100% of communities
             lwd = 1.2) +
  geom_vline(xintercept = c(0.95, 0.90, 0.75), # Draw the threshold for the boxplots
             col = "grey20", lty = 3, lwd = 0.3) +
  expand_limits(y = 0, x = 0) + # To tell to always keep those values inside the plot. # Set the origin for the graph
  scale_x_reverse(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, .05))) + # To tell not to add padding around the limits of the axis
  # scale_x_reverse() + # Reverse the x-axis to have an increase in the number of communities included at each iteration
  labs(y = "mean HII", x = "Normalized indices",
       title = "Full map ; Mean HII") +
  theme_bw() +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness"))
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity"))
print(g)
dev.off()

### 3.2.2/ median HII ####

# Store results for median in a tidy df
most_diverse_norm_median_df <- tibble(sp.richness = sp.richness_median_values,
                                      sp.rarity = sp.rarity_median_values,
                                      # MPD = MPD_median_values,
                                      ring.richness = ring.richness_median_values,
                                      ring.rarity = ring.rarity_median_values)

most_diverse_norm_median_df <- most_diverse_norm_median_df %>% 
  mutate(threshold = 1-(1:1000)/1000) %>%    # Add a column for threhold norm values
  mutate(stat = "median") %>%             # Add a column for stat computed
  filter(threshold <= 1) %>%              # Remove value for the top thesholds because curves are idiosyncratic (nb of communities too low)
  pivot_longer(# cols = c("sp.richness", "sp.rarity", "MPD", "ring.richness"), # List of the name of the "columns" that are actually factor levels
    cols = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"), # List of the name of the "columns" that are actually factor levels
    names_to = "index",   # name of the new column that is the real variable for these factor level
    values_to = "HII", # name of the new columns that will get the values in the cells of the previous dataset
    values_drop_na = T) %>% # Drop rows for missing values ?
  mutate(index = factor(index, 
                        # levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness"))) # Change factor level order
                        levels = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))) # Change factor level order

# Save results
saveRDS(most_diverse_norm_median_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_median_df.rds", version = "2")

# Load results
most_diverse_norm_median_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_median_df.rds")

# Plot resulting curves for median HII
pdf(file = paste0("./graphs/Threat_curves/Most_diverse/Most_diverse_norm_global_median.pdf"), height = 5.3, width = 6.5)
g <- ggplot(most_diverse_norm_median_df, aes(x = threshold, y = HII)) +
  geom_line(aes(color = index), lwd = 2) +
  geom_hline(yintercept = median(most_diverse_norm_median_df[most_diverse_norm_median_df$threshold == 0, "HII", drop = T]), # Get the final median HII for 100% of communities
             lwd = 1.2) +
  geom_vline(xintercept = c(0.95, 0.90, 0.75), # Draw the threshold for the boxplots
             col = "grey20", lty = 3, lwd = 0.3) +
  expand_limits(y = 0, x = 0) + # To tell to always keep those values inside the plot. # Set the origin for the graph
  scale_x_reverse(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, .05))) + # To tell not to add padding around the limits of the axis
  # scale_x_reverse() + # Reverse the x-axis to have an increase in the number of communities included at each iteration
  labs(y = "median HII", x = "Normalized indices",
       title = "Full map ; Median HII") +
  theme_bw() +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness"))
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity"))
print(g)
dev.off()


### 3.2.3/ HII 95% ####

# Store results for mean in a tidy df
most_diverse_norm_95_df <- tibble(sp.richness = sp.richness_95_values,
                                  sp.rarity = sp.rarity_95_values,
                                  # MPD = MPD_95_values,
                                  ring.richness = ring.richness_95_values,
                                  ring.rarity = ring.rarity_95_values)

most_diverse_norm_95_df <- most_diverse_norm_95_df %>% 
  mutate(threshold = 1-(1:1000)/1000) %>%    # Add a column for threhold norm values
  mutate(stat = "95") %>%             # Add a column for stat computed
  filter(threshold <= 1) %>%              # Remove value for the top thesholds because curves are idiosyncratic (nb of communities too low)
  pivot_longer(# cols = c("sp.richness", "sp.rarity", "MPD", "ring.richness"), # List of the name of the "columns" that are actually factor levels
    cols = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"), # List of the name of the "columns" that are actually factor levels
    names_to = "index",   # name of the new column that is the real variable for these factor level
    values_to = "HII", # name of the new columns that will get the values in the cells of the previous dataset
    values_drop_na = T) %>% # Drop rows for missing values ?
  mutate(index = factor(index, 
                        # levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness"))) # Change factor level order
                        levels = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))) # Change factor level order

# Save results
saveRDS(most_diverse_norm_95_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_95_df.rds", version = "2")

# Load results
most_diverse_norm_95_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_95_df.rds")

# Plot resulting curves for HII 95%
pdf(file = paste0("./graphs/Threat_curves/Most_diverse/Most_diverse_norm_global_95.pdf"), height = 5.3, width = 6.5)
g <- ggplot(most_diverse_norm_95_df, aes(x = threshold, y = HII)) +
  geom_line(aes(color = index), lwd = 2) +
  geom_hline(yintercept = quantile(most_diverse_norm_95_df[most_diverse_norm_95_df$threshold == 0, "HII", drop = T], probs = 0.95), # Get the final quantile 0.95 of HII for 100% of communities
             lwd = 1.2) +
  geom_vline(xintercept = c(0.95, 0.90, 0.75), # Draw the threshold for the boxplots
             col = "grey20", lty = 3, lwd = 0.3) +
  expand_limits(y = 0, x = 0) + # To tell to always keep those values inside the plot. # Set the origin for the graph
  scale_x_reverse(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, .05))) + # To tell not to add padding around the limits of the axis
  # scale_x_reverse() + # Reverse the x-axis to have an increase in the number of communities included at each iteration
  labs(y = "HII (quantile 0.95)", x = "Normalized indices",
       title = "Full map ; HII 95%") +
  theme_bw() +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness"))
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity"))
print(g)
dev.off()

### 3.2.4/ all HII quantiles together (50%, 95%) ####

# Store results for mean in a tidy df
most_diverse_norm_all_quantiles_df <- rbind(most_diverse_norm_median_df, most_diverse_norm_95_df)

# Save results
saveRDS(most_diverse_norm_all_quantiles_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_all_quantiles_df.rds", version = "2")

# Load results
most_diverse_norm_all_quantiles_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_all_quantiles_df.rds")

# Extract only median and quantile 95%
most_diverse_norm_median_95_df <- most_diverse_norm_all_quantiles_df %>% 
  filter(stat %in% c("median","95"))

saveRDS(most_diverse_norm_median_95_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_median_95_df.rds", version = "2")

# Plot resulting curves for HII quantiles 50% and 95%
pdf(file = paste0("./graphs/Threat_curves/Most_diverse/Most_diverse_norm_global_median_95.pdf"), height = 5.3, width = 6.5)
g <- ggplot(most_diverse_norm_median_95_df, aes(x = threshold, y = HII)) +
  geom_line(aes(color = index, lty = stat), lwd = 2) +
  geom_hline(yintercept = quantile(x = most_diverse_norm_median_95_df[(most_diverse_norm_median_95_df$threshold == 0), "HII", drop = T], probs = 0.95), # Get the final quantile 95% HII for 100% of communities
             lwd = 1.2) +
  geom_hline(yintercept = median(x = most_diverse_norm_median_95_df[(most_diverse_norm_median_95_df$threshold == 0) & (most_diverse_norm_median_95_df$stat == "median"), "HII", drop = T]), # Get the final median HII for 100% of communities
             lwd = 1.2) +
  geom_vline(xintercept = c(0.95, 0.90, 0.75), # Draw the threshold for the boxplots
             col = "grey20", lty = 3, lwd = 0.3) +
  expand_limits(y = 0, x = 0) + # To tell to always keep those values inside the plot. # Set the origin for the graph
  scale_x_reverse(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, .05))) + # To tell not to add padding around the limits of the axis
  labs(y = "HII (median & quantile 0.95)", x = "Normalized indices",
       title = "Full map ; HII median and 95%") +
  theme_bw() +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness"))
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity"))
print(g)
dev.off()

### 3.3/ Boxplots for specific normalized index values ####

# Prepare data

# Load data for boxplots
most_diverse_norm_global_boxplot_values <- load(file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_global_boxplot_values.rds")

most_diverse_norm_global_boxplot_lists <- c(sp.richness_boxplot_values, 
                                            sp.rarity_boxplot_values, 
                                            # MPD_boxplot_values, 
                                            ring.richness_boxplot_values,
                                            ring.rarity_boxplot_values)

# Counts nb of observations recorded for each threshold and each index
n_com <- unlist(map(.x = most_diverse_norm_global_boxplot_lists, .f = length))

# Build the final tidy df
most_diverse_norm_global_boxplot_df <- tibble(HII = unlist(most_diverse_norm_global_boxplot_lists)) %>% 
  mutate(threshold = c(rep("0.95", n_com[1]),
                       rep("0.90", n_com[2]),
                       rep("0.75", n_com[3]),
                       rep("0.95", n_com[4]),
                       rep("0.90", n_com[5]),
                       rep("0.75", n_com[6]),
                       rep("0.95", n_com[7]),
                       rep("0.90", n_com[8]),
                       rep("0.75", n_com[9]),
                       rep("0.95", n_com[10]),
                       rep("0.90", n_com[11]),
                       rep("0.75", n_com[12]))) %>% 
  mutate(index = c(rep("sp.richness", sum(n_com[1:3])),
                   rep("sp.rarity", sum(n_com[4:6])),
                   # rep("MPD", sum(n_com[7:9])),
                   # rep("ring.richness", sum(n_com[10:12])),
                   rep("ring.richness", sum(n_com[7:9])),
                   rep("ring.rarity", sum(n_com[10:12])))) %>% 
  mutate(index = factor(index, 
                        # levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness"))) # Change factor level order
                        levels = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))) # Change factor level order

saveRDS(most_diverse_norm_global_boxplot_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_global_boxplot_values.rds", version = "2")


### 3.3.1/ Boxplot for normalized index threshold = 0.95 ####

# Load data
most_diverse_norm_global_boxplot_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_global_boxplot_values.rds")

# Filter data
most_diverse_norm_global_boxplot_0.95_df <- most_diverse_norm_global_boxplot_df %>% 
  filter(threshold == "0.95")

# Compute quantiles 95 to add to the plot
quant_95 <- as.numeric(with(most_diverse_norm_global_boxplot_0.95_df, by(HII, INDICES = index, FUN = quantile, probs = 0.95)))
quant_95 <- tibble(quant_value = quant_95,
                   # index = c("sp.richness", "sp.rarity", "MPD", "ring.richness"))
                   index = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))
most_diverse_norm_global_boxplot_0.95_df <- most_diverse_norm_global_boxplot_0.95_df %>% 
  left_join(y = quant_95, by = "index") %>% 
  mutate(index = factor(index, 
                        # levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness"))) # Change factor level order
                        levels = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))) # Change factor level order

saveRDS(most_diverse_norm_global_boxplot_0.95_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_global_boxplot_0.95_df.rds", version = "2")

# Plot the boxplot
pdf(file = paste0("./graphs/Threat_curves/Most_diverse/Most_diverse_norm_global_boxplot_0.95.pdf"), height = 5.3, width = 6.5)
g <- ggplot(most_diverse_norm_global_boxplot_0.95_df, aes(x = index, y = HII), show.legend = F) +
  geom_boxplot(aes(fill = index)) +
  geom_point(aes(x = index, y = quant_value, col = index), size = 5) +
  geom_hline(yintercept = median(most_diverse_norm_median_df[most_diverse_norm_median_df$threshold == 0, "HII", drop = T]), # Get the final median HII for 100% of communities
             lty = 2, lwd = 0.8) +
  labs(y = "HII", x = "Index",
       title = "Full map ; Threshold normalized indices = 0.95") +
  theme_bw() +
  # scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                    name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  guides(color = "none")
print(g)
dev.off()

### 3.3.2/ Boxplot for normalized index threshold = 0.90 ####

# Load data
most_diverse_norm_global_boxplot_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_global_boxplot_values.rds")

# Filter data
most_diverse_norm_global_boxplot_0.90_df <- most_diverse_norm_global_boxplot_df %>% 
  filter(threshold == "0.90")

# Compute quantiles 95 to add to the plot
quant_95 <- as.numeric(with(most_diverse_norm_global_boxplot_0.90_df, by(HII, INDICES = index, FUN = quantile, probs = 0.95)))
quant_95 <- tibble(quant_value = quant_95,
                   # index = c("sp.richness", "sp.rarity", "MPD", "ring.richness"))
                   index = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))
most_diverse_norm_global_boxplot_0.90_df <- most_diverse_norm_global_boxplot_0.90_df %>% 
  left_join(y = quant_95, by = "index") %>% 
  mutate(index = factor(index, 
                        # levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness"))) # Change factor level order
                        levels = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))) # Change factor level order

saveRDS(most_diverse_norm_global_boxplot_0.90_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_global_boxplot_0.90_df.rds", version = "2")

# Plot the boxplot
pdf(file = paste0("./graphs/Threat_curves/Most_diverse/Most_diverse_norm_global_boxplot_0.90.pdf"), height = 5.3, width = 6.5)
g <- ggplot(most_diverse_norm_global_boxplot_0.90_df, aes(x = index, y = HII), show.legend = F) +
  geom_boxplot(aes(fill = index)) +
  geom_point(aes(x = index, y = quant_value, col = index), size = 5) +
  geom_hline(yintercept = median(most_diverse_norm_median_df[most_diverse_norm_median_df$threshold == 0, "HII", drop = T]), # Get the final median HII for 100% of communities
             lty = 2, lwd = 0.8) +
  labs(y = "HII", x = "Index",
       title = "Full map ; Threshold normalized indices = 0.90") +
  theme_bw() +
  # scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                    name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  guides(color = "none")
print(g)
dev.off()

### 3.3.3/ Boxplot for normalized index threshold = 0.75 ####

# Load data
most_diverse_norm_global_boxplot_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_global_boxplot_values.rds")

# Filter data
most_diverse_norm_global_boxplot_0.75_df <- most_diverse_norm_global_boxplot_df %>% 
  filter(threshold == "0.75")

# Compute quantiles 95 to add to the plot
quant_95 <- as.numeric(with(most_diverse_norm_global_boxplot_0.75_df, by(HII, INDICES = index, FUN = quantile, probs = 0.95)))
quant_95 <- tibble(quant_value = quant_95,
                   # index = c("sp.richness", "sp.rarity", "MPD", "ring.richness"))
                   index = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))
most_diverse_norm_global_boxplot_0.75_df <- most_diverse_norm_global_boxplot_0.75_df %>% 
  left_join(y = quant_95, by = "index") %>% 
  mutate(index = factor(index, 
                        # levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness"))) # Change factor level order
                        levels = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))) # Change factor level order

saveRDS(most_diverse_norm_global_boxplot_0.75_df, file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_global_boxplot_0.75_df.rds", version = "2")

# Plot the boxplot
pdf(file = paste0("./graphs/Threat_curves/Most_diverse/Most_diverse_norm_global_boxplot_0.75.pdf"), height = 5.3, width = 6.5)
g <- ggplot(most_diverse_norm_global_boxplot_0.75_df, aes(x = index, y = HII), show.legend = F) +
  geom_boxplot(aes(fill = index)) +
  geom_point(aes(x = index, y = quant_value, col = index), size = 5) +
  geom_hline(yintercept = median(most_diverse_norm_median_df[most_diverse_norm_median_df$threshold == 0, "HII", drop = T]), # Get the final median HII for 100% of communities
             lty = 2, lwd = 0.8) +
  labs(y = "HII", x = "Index",
       title = "Full map ; Threshold normalized indices = 0.75") +
  theme_bw() +
  # scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                    name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  guides(color = "none")
print(g)
dev.off()

### 3.4/ Multiple plot: curves and boxplots ####

most_diverse_norm_median_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_median_df.rds")
most_diverse_norm_global_boxplot_0.95_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_global_boxplot_0.95_df.rds")
most_diverse_norm_global_boxplot_0.90_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_global_boxplot_0.90_df.rds")
most_diverse_norm_global_boxplot_0.75_df <- readRDS(file = "./outputs/Threat_curves/Most_diverse/most_diverse_norm_global_boxplot_0.75_df.rds")

# Plot the boxplot
pdf(file = paste0("./graphs/Threat_curves/Most_diverse/Most_diverse_norm_global_curves_boxplot.pdf"), height = 11, width = 13)

g1 <- ggplot(most_diverse_norm_median_df, aes(x = threshold, y = HII)) +
  geom_line(aes(color = index), lwd = 2) +
  geom_hline(yintercept = median(most_diverse_norm_median_df[most_diverse_norm_median_df$threshold == 0, "HII", drop = T]), # Get the final median HII for 100% of communities
             lwd = 1.2) +
  geom_vline(xintercept = c(0.95, 0.90, 0.75), # Draw the threshold for the boxplots
             col = "grey20", lty = 3, lwd = 0.3) +
  expand_limits(y = 0, x = 0) + # To tell to always keep those values inside the plot. # Set the origin for the graph
  scale_x_reverse(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, .05))) + # To tell not to add padding around the limits of the axis
  # scale_x_reverse() + # Reverse the x-axis to have an increase in the number of communities included at each iteration
  labs(y = "median HII", x = "Normalized indices",
       title = "Full map ; Median HII") +
  theme_bw() +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness"))
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity"))
  
g2 <- ggplot(most_diverse_norm_global_boxplot_0.95_df, aes(x = index, y = HII), show.legend = F) +
  geom_boxplot(aes(fill = index)) +
  geom_point(aes(x = index, y = quant_value, col = index), size = 5) +
  geom_hline(yintercept = median(most_diverse_norm_median_df[most_diverse_norm_median_df$threshold == 0, "HII", drop = T]), # Get the final median HII for 100% of communities
             lty = 2, lwd = 0.8) +
  labs(y = "HII", x = "Index",
       title = "Full map ; Threshold normalized indices = 0.95") +
  theme_bw() +
  # scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                    name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  guides(color = "none")

g3 <- ggplot(most_diverse_norm_global_boxplot_0.90_df, aes(x = index, y = HII), show.legend = F) +
  geom_boxplot(aes(fill = index)) +
  geom_point(aes(x = index, y = quant_value, col = index), size = 5) +
  geom_hline(yintercept = median(most_diverse_norm_median_df[most_diverse_norm_median_df$threshold == 0, "HII", drop = T]), # Get the final median HII for 100% of communities
             lty = 2, lwd = 0.8) +
  labs(y = "HII", x = "Index",
       title = "Full map ; Threshold normalized indices = 0.90") +
  theme_bw() +
  # scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                    name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  guides(color = "none")

g4 <- ggplot(most_diverse_norm_global_boxplot_0.75_df, aes(x = index, y = HII), show.legend = F) +
  geom_boxplot(aes(fill = index)) +
  geom_point(aes(x = index, y = quant_value, col = index), size = 5) +
  geom_hline(yintercept = median(most_diverse_norm_median_df[most_diverse_norm_median_df$threshold == 0, "HII", drop = T]), # Get the final median HII for 100% of communities
             lty = 2, lwd = 0.8) +
  labs(y = "HII", x = "Index",
       title = "Full map ; Threshold normalized indices = 0.75") +
  theme_bw() +
  # scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness")) +
  scale_fill_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                    name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity")) +
  guides(color = "none")

gridExtra::grid.arrange(g1, g2, g3, g4, ncol=2)


dev.off()




### 4/ Version for most threatened communities ####


# Load normalized indices map
sp.richness_norm <- readRDS(file = "./outputs/Threat_maps/sp.richness_norm.rds")
sp.rarity_norm <- readRDS(file = "./outputs/Threat_maps/sp.rarity_norm.rds")
# MPD_norm <- readRDS(file = "./outputs/Threat_maps/MPD_norm.rds")
ring.richness_norm <- readRDS(file = "./outputs/Threat_maps/ring.richness_norm.rds")
ring.rarity_norm <- readRDS(file = "./outputs/Threat_maps/ring.rarity_norm.rds")

# Extract HII values for communities with Ithomiini
HII_masked <- mask(HII, Community_mask)

# Extract indices of communities with Ithomiini ordred by decreasing value of HII
HII_indices <- order(HII_masked[], decreasing = T)
HII_indices <- HII_indices[1:N.com]

sp.richness_norm_mean_values <- sp.rarity_norm_mean_values <- MPD_norm_mean_values <- ring.richness_norm_mean_values <- ring.rarity_norm_mean_values <- integer()
for (i in 1:1000)
{
  # i <- 1
  
  # Compute nb of communities to extract
  N.extract <- round(length(HII_indices) * i/1000)
  # Get indices of communities to extract
  Index_extract <- HII_indices[1:N.extract]
  
  # Extract values and compute mean
  sp.richness_norm_values <- sp.richness_norm[Index_extract]
  sp.richness_norm_mean_values[i] <- round(mean(sp.richness_norm_values, na.rm = T), 3)
  
  sp.rarity_norm_values <- sp.rarity_norm[Index_extract]
  sp.rarity_norm_mean_values[i] <- round(mean(sp.rarity_norm_values, na.rm = T), 3)
  
  MPD_norm_values <- MPD_norm[Index_extract]
  MPD_norm_mean_values[i] <- round(mean(MPD_norm_values, na.rm = T), 3)
  
  ring.richness_norm_values <- ring.richness_norm[Index_extract]
  ring.richness_norm_mean_values[i] <- round(mean(ring.richness_norm_values, na.rm = T), 3)
  
  ring.rarity_norm_values <- ring.rarity_norm[Index_extract]
  ring.rarity_norm_mean_values[i] <- round(mean(ring.rarity_norm_values, na.rm = T), 3)

  if(i %% 100 == 0) {print(i)}  
}

# Store results for mean in a tidy df
most_threatened_mean_df <- tibble(sp.richness = sp.richness_norm_mean_values,
                                  sp.rarity = sp.rarity_norm_mean_values,
                                  # MPD = MPD_norm_mean_values,
                                  ring.richness = ring.richness_norm_mean_values,
                                  ring.rarity = ring.rarity_norm_mean_values)

most_threatened_mean_df <- most_threatened_mean_df %>% 
  mutate(prop_com = (1:1000)/10) %>% 
  mutate(stat = "mean") %>%             # Add a column for stat computed
  filter(prop_com > 2) %>%              # Remove value for the top 2% because curves are idiosyncratic (nb of communities too low)
  pivot_longer(# cols = c("sp.richness", "sp.rarity", "MPD", "ring.richness"), # List of the name of the "columns" that are actually factor levels
               cols = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"), # List of the name of the "columns" that are actually factor levels
               names_to = "index",   # name of the new column that is the real variable for these factor level
               values_to = "values", # name of the new columns that will get the values in the cells of the previous dataset
               values_drop_na = T) %>% # Drop rows for missing values ?
  mutate(index = factor(index, 
                        # levels = c("sp.richness", "sp.rarity", "MPD", "ring.richness"))) # Change factor level order
                        levels = c("sp.richness", "sp.rarity", "ring.richness", "ring.rarity"))) # Change factor level order

# Save results
saveRDS(most_threatened_mean_df, file = "./outputs/Threat_curves/Most_threatened/most_threatened_mean_df.rds", version = "2")

# Load results
most_threatened_mean_df <- readRDS(file = "./outputs/Threat_curves/Most_threatened/most_threatened_mean_df.rds")

# Plot resulting curves for mean HII
pdf(file = paste0("./graphs/Threat_curves/Most_threatened/Most_threatened_rank_global_mean.pdf"), height = 5.3, width = 6.5)
g <- ggplot(most_threatened_mean_df, aes(x = prop_com, y = values)) +
  geom_line(aes(color = index), lwd = 2) +
  geom_hline(yintercept = last(sp.richness_norm_mean_values), # Get the final mean HII for 100% of communities
             lwd = 0.8) +
  geom_hline(yintercept = last(sp.rarity_norm_mean_values), # Get the final mean HII for 100% of communities
             lwd = 0.8) +
  # geom_hline(yintercept = last(MPD_norm_mean_values), # Get the final mean HII for 100% of communities
  #            lwd = 0.8) +
  geom_hline(yintercept = last(ring.richness_norm_mean_values), # Get the final mean HII for 100% of communities
             lwd = 0.8) +
  geom_hline(yintercept = last(ring.rarity_norm_mean_values), # Get the final mean HII for 100% of communities
             lwd = 0.8) +
  geom_vline(xintercept = c(5, 10, 25), # Draw the threshold for the boxplots
             col = "grey20", lty = 3, lwd = 0.3) +
  expand_limits(y = 0, x = 0) + # To tell to always keep those values inside the plot. # Set the origin for the graph
  scale_x_continuous(expand = expansion(mult = c(0, 0))) + scale_y_continuous(expand = expansion(mult = c(0, .05))) + # To tell not to add padding around the limits of the axis
  labs(y = "Mean normalized index", x = "Proportion of communities",
       title = "Full map ; Mean normalized indices") +
  theme_bw() +
  # scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "MPD" = "darkorange", "ring.richness" = "limegreen"),
  #                    name = "Indices", labels = c("Species richness", "Species rarity", "MPD", "Mimicry richness"))
  scale_color_manual(values = c("sp.richness" = "red", "sp.rarity" = "dodgerblue", "ring.richness" = "darkorange", "ring.rarity" = "limegreen"),
                     name = "Indices", labels = c("Species richness", "Species rarity", "Mimicry richness", "Mimicry rarity"))
print(g)
dev.off()



# Do also climatic anomalies and deforestation ???



### Aggregate Ecoregion to make bioregions

