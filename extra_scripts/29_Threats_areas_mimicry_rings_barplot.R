
###### Script 29: Barplots of mimicry rings ranked by threatened areas  ######

# Only for Jaccard.80

### Two types of bins for HII: quantiles and absolute values

###

# Inputs: 
   # Ring probability maps from script 13
   # Threat map: HII (HF, Venter et al., 2016)

# Outputs:
   # Histrogram of proportion of areas under threat per mimicry rings
       # Per quantiles 5 & 25% 
       # Per absolute values (1, 5 vs. 11, 23)
   # Barplot with rings ordered by perc. of threatened area
       # Per quantiles 5 & 25% 
       # Per absolute values (1, 5 vs. 11, 23)

###

### 1/ Load stuff ####

# Clean environment
rm(list = ls())

library(raster)
library(rangeBuilder)

# Load masks for HII quantiles
load(file = "./outputs/Threat_maps/HII_cropped/HII_masks_quantiles_files.rds") 

# Load masks for HII bins based on absolute values
load(file = "./outputs/Threat_maps/HII_cropped/HII_masks_bins_files.rds")

# Load mimicry ring probability maps
All_ring_proba_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_Jaccard.80.RData"))

plot(All_ring_proba_stack_Jaccard.80)

### 2/ Compute % areas under threats for each ring ####

threatened_rank_25 <- threatened_rank_5 <- threatened_high <- threatened_very_high <- NA
threatened_rank_25_perc <- threatened_rank_5_perc <- threatened_high_perc <- threatened_very_high_perc <- NA
for (i in 1:nlayers(All_ring_proba_stack_Jaccard.80))
{
  # Get total area
  ring.map <- All_ring_proba_stack_Jaccard.80[[i]]
  all_area <- sum(ring.map[], na.rm = T)
  
  # Extract for quantiles 25%
  threatened_rank_25[i] <- sum((HII_mask_rank_25*ring.map)[], na.rm = T)
  threatened_rank_25_perc[i] <- threatened_rank_25[i]/all_area*100
  
  # Extract for quantiles 5%
  threatened_rank_5[i] <- sum((HII_mask_rank_5*ring.map)[], na.rm = T)
  threatened_rank_5_perc[i] <- threatened_rank_5[i]/all_area*100
  
  # Extract for HII > 11
  threatened_high[i] <- sum((HII_mask_high*ring.map)[], na.rm = T)
  threatened_high_perc[i] <- threatened_high[i]/all_area*100
  
  # Extract for HII > 23
  threatened_very_high[i] <- sum((HII_mask_very_high*ring.map)[], na.rm = T)
  threatened_very_high_perc[i] <- threatened_very_high[i]/all_area*100
  
  print(i)
}

# Make regional categories for rings: Andean, Amazonian, Global, MA, CA
ring_regions_types <- c("Global", "Andean", "Global", "Andean", "Central America", "Andean", "Global", "Amazonian", "Amazonian", "Global", "Amazonian", "Amazonian", "Global", "Andean", "Global", "Amazonian", "Andean", "Andean", "Amazonian", "Andean", "Andean", "Andean", "Global", "Andean", "Andean", "Andean", "Mata Atlantica", "Andean", "Andean", "Andean", "Andean", "Andean", "Mata Atlantica", "Amazonian", "Amazonian", "Andean", "Andean", "Central America", "Andean", "Andean", "Andean", "Andean", "Andean", "Amazonian")

# Make a final table
ring_threats_perc_df <- data.frame(cbind(ring_name, ring_regions_types, threatened_rank_25_perc, threatened_rank_5_perc, threatened_high_perc, threatened_very_high_perc))
ring_threats_perc_df[, 3:6] <- apply(X = ring_threats_perc_df[, 3:6], MARGIN = 2, FUN = as.character)
ring_threats_perc_df[, 3:6] <- apply(X = ring_threats_perc_df[, 3:6], MARGIN = 2, FUN = as.numeric)
ring_threats_perc_df$ring_name <- factor(x = names(All_ring_proba_stack_Jaccard.80), levels = names(All_ring_proba_stack_Jaccard.80))

saveRDS(ring_threats_perc_df, file = "./outputs/Threat_maps/Overlaps/ring_threats_perc_df.rds", version = "2")

### 3/ Try histogram, but not really interesting... ####

hist(threatened_rank_25_perc)
hist(threatened_high_perc)
hist(threatened_very_high_perc)
hist(threatened_rank_5_perc)

### 4/ Try barplot with rings ordered by perc. of threatened area ####

library(tidyverse)

# 4.1/ For threat as Top 25% ####

pdf(file = "./graphs/Overlaps/ring_threats_barplot_rank_25.pdf", height = 6, width = 10)

ring_threats_perc_df <- ring_threats_perc_df %>%
  arrange(threatened_rank_25_perc)

ggplot(ring_threats_perc_df, aes(x = ring_name, y = threatened_rank_25_perc, fill = ring_regions_types)) +
  geom_col(show.legend = T) +
  labs(title = "Threat: Top 25%", y = "Area threatened (%)", x = NULL) +
  scale_fill_manual(name = "Regions",
                    values = c("Andean" = "firebrick1", "Amazonian" = "chartreuse3", "Global" = "darkgoldenrod1", "Mata Atlantica" = "steelblue1", "Central America" = "plum2")) +
  scale_x_discrete(limits = ring_threats_perc_df$ring_name) +
  theme_classic() +
  theme(title = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 5, 0)),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10))

dev.off()

# 4.2/ For threat as Top 5% ####

pdf(file = "./graphs/Overlaps/ring_threats_barplot_rank_5.pdf", height = 6, width = 10)

ring_threats_perc_df <- ring_threats_perc_df %>%
  arrange(threatened_rank_5_perc)

ggplot(ring_threats_perc_df, aes(x = ring_name, y = threatened_rank_5_perc, fill = ring_regions_types)) +
  geom_col(show.legend = T) +
  labs(title = "Threat: Top 5%", y = "Area threatened (%)", x = NULL) +
  scale_fill_manual(name = "Regions",
                    values = c("Andean" = "firebrick1", "Amazonian" = "chartreuse3", "Global" = "darkgoldenrod1", "Mata Atlantica" = "steelblue1", "Central America" = "plum2")) +
  scale_x_discrete(limits = ring_threats_perc_df$ring_name) +
  theme_classic() +
  theme(title = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 5, 0)),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10))

dev.off()

# 4.3/ For threat as HII > 11 ####

pdf(file = "./graphs/Overlaps/ring_threats_barplot_bin_high.pdf", height = 6, width = 10)

ring_threats_perc_df <- ring_threats_perc_df %>%
  arrange(threatened_high_perc)

ggplot(ring_threats_perc_df, aes(x = ring_name, y = threatened_high_perc, fill = ring_regions_types)) +
  geom_col(show.legend = T) +
  labs(title = "Threat: HII > 11", y = "Area threatened (%)", x = NULL) +
  scale_fill_manual(name = "Regions",
                    values = c("Andean" = "firebrick1", "Amazonian" = "chartreuse3", "Global" = "darkgoldenrod1", "Mata Atlantica" = "steelblue1", "Central America" = "plum2")) +
  scale_x_discrete(limits = ring_threats_perc_df$ring_name) +
  theme_classic() +
  theme(title = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 5, 0)),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10))

dev.off()

# 4.4/ For threat as HII > 23 ####

pdf(file = "./graphs/Overlaps/ring_threats_barplot_bin_very_high.pdf", height = 6, width = 10)

ring_threats_perc_df <- ring_threats_perc_df %>%
  arrange(threatened_very_high_perc)

ggplot(ring_threats_perc_df, aes(x = ring_name, y = threatened_very_high_perc, fill = ring_regions_types)) +
  geom_col(show.legend = T) +
  labs(title = "Threat: HII > 23", y = "Area threatened (%)", x = NULL) +
  scale_fill_manual(name = "Regions",
                    values = c("Andean" = "firebrick1", "Amazonian" = "chartreuse3", "Global" = "darkgoldenrod1", "Mata Atlantica" = "steelblue1", "Central America" = "plum2")) +
  scale_x_discrete(limits = ring_threats_perc_df$ring_name) +
  theme_classic() +
  theme(title = element_text(size = 14, hjust = 0.5, margin = margin(0, 0, 5, 0)),
        axis.text = element_text(size = 12, face = "bold", color = "black"),
        axis.title = element_text(size = 15),
        axis.line = element_line(color = "black"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size = 10))

dev.off()
