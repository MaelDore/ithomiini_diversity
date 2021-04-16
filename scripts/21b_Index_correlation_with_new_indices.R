
### Prepare stuff ####

# Effacer l'environnement
rm(list = ls())

library(raster)
library(rangeBuilder)

### Load maps

# Main figure
sp.richness <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
# sp.mean.rarity <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_Jaccard.80.rds")       # Linear weighting
sp.mean.rarity <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_Leroy_Jaccard.80.rds")   # Leroy's weighting
Faith.PD <- readRDS(file = "./outputs/Indices_Maps/PD.raster_Jaccard.80.rds")
ring.richness <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))
# mimicry.mean.rarity <- readRDS(file = "./outputs/Indices_maps/mimicry.mean.rarity_Jaccard.80.rds")     # Linear weighting
mimicry.mean.rarity <- readRDS(file = "./outputs/Indices_maps/ring.mean.rarity_Leroy_Jaccard.80.rds")    # Leroy's weighting
weighted_mean_ring_size <- readRDS(file = "./outputs/Indices_maps/weighted_mean_ring_size.rds")


# Extra indices for SM figure
sp.diversity <- readRDS(file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jost_Jaccard.80.rds"))
ring.diversity <- readRDS(file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jost_Jaccard.80.rds"))
ED <- readRDS(file = "./outputs/Indices_Maps/FP.raster_Jaccard.80.rds")
load(file = "./outputs/Indices_maps/max_ring_size.RData")

# Removed index
# MPD <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80.rds")

# Generate a mask only for communities within Ithominii range
Ithomiini_range_mask <- sp.richness > 0
Ithomiini_range_mask[Ithomiini_range_mask[] == 0] <-  NA
plot(Ithomiini_range_mask)

# Generate stack for analysis 
all_indices_stack <- stack(sp.richness, sp.mean.rarity, Faith.PD, ring.richness, mimicry.mean.rarity, weighted_mean_ring_size, sp.diversity, ring.diversity, ED, max_ring_size)
all_indices_stack <- mask(x = all_indices_stack, mask = Ithomiini_range_mask) # Remove data outside Ithomiini range
all_indices_stack <- virtualspecies::synchroniseNA(all_indices_stack) # Synchronize NA between layers

names(all_indices_stack) <- indices_list <- c("Sp. richness", "Sp. rarity", "Faith's PD", "Mim. richness", "Mim. rarity", "Mean ring size", "Sp. diversity", "Mim. diversity", "ED", "Max ring size")

save(all_indices_stack, file = "./outputs/Correlation_tests/all_indices_stack.RData", version = 2)
saveRDS(all_indices_stack, file = "./outputs/Correlation_tests/all_indices_stack.rds", version = 2)

pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")

plot(all_indices_stack, col = pal_bl_red_Mannion)

# Extract only selected indices for main text
main_indices_stack <- subset(all_indices_stack, 1:6)
main_indices_list <- indices_list[1:6]

### Choose data extraction specification (for corrected tests)

all_coords <- coordinates(all_indices_stack[[1]])

# Version with all communities
sample.all_index <- which(all_indices_stack[[1]][] > 0) 

# Version with subsampling
# set.seed <- 45512
# N.sample = 1000
# sample.all_index <- sample(x = which(Ithomiini_range_mask[] == 1), size = N.sample, replace = F)

# Get coordinates
sample.all_coords <- all_coords[sample.all_index,]

# Check position of sampled communities
plot(Ithomiini_range_mask)
points(sample.all_coords, pch = 16, cex = 0.5)

# Sample indices
sample_indices_df <- data.frame(raster::extract(x = all_indices_stack, y = sample.all_coords))

# Save stuff
save(sample.all_index, file =  paste0("./outputs/Correlation_tests/sample.all_index.RData"), version = "2")
saveRDS(sample.all_index, file =  paste0("./outputs/Correlation_tests/sample.all_index.rds"), version = "2")
save(sample.all_coords, file =  paste0("./outputs/Correlation_tests/sample.all_coords.RData"), version = "2")
saveRDS(sample.all_coords, file =  paste0("./outputs/Correlation_tests/sample.all_coords.rds"), version = "2")
save(sample_indices_df, file =  paste0("./outputs/Correlation_tests/sample_indices_df.RData"), version = "2")
saveRDS(sample_indices_df, file =  paste0("./outputs/Correlation_tests/sample_indices_df.rds"), version = "2")

# Load indices and coordinates
sample.all_index <- readRDS(file =  paste0("./outputs/Correlation_tests/sample.all_index.rds"))
sample.all_coords <- readRDS(file =  paste0("./outputs/Correlation_tests/sample.all_coords.rds"))
sample_indices_df <- readRDS(file =  paste0("./outputs/Correlation_tests/sample_indices_df.rds"))

# Load stack
all_indices_stack <- readRDS(file = "./outputs/Correlation_tests/all_indices_stack.rds")
main_indices_stack <- subset(all_indices_stack, 1:6)

### 1/ Approach with uncorrected Spearman's rank test ####

# With cor.test, easier, no need to synchronize the test
N.indices <- nlayers(all_indices_stack)

### 1.1/ Build a matrix of pairwise rho coefficients and associated p-values ####
rho_matrix <- rho_uncorrected_pvalues_matrix <- matrix(nrow = N.indices, ncol = N.indices)
for (i in 1:N.indices)
{
  # for (j in (i+1):N.indices) # Most efficient version with no useless double tests
  for (j in 1:N.indices) # Version to fill the all matrix
  {
    test <- cor.test(x = all_indices_stack[[i]][], y = all_indices_stack[[j]][], method = "spearman", alternative = "two.sided", na.rm = T)
    rho_matrix[i,j] <- round(test$estimate, 3)
    rho_uncorrected_pvalues_matrix[i,j] <- round(test$p.value, 3)
  }
}

rownames(rho_matrix) <- colnames(rho_matrix) <- indices_list
rownames(rho_uncorrected_pvalues_matrix) <- colnames(rho_uncorrected_pvalues_matrix) <- indices_list

# Save rho matrix, and uncorrected p-values matrix
saveRDS(rho_matrix, file =  paste0("./outputs/Correlation_tests/rho_matrix.rds"), version = "2")
saveRDS(rho_uncorrected_pvalues_matrix, file =  paste0("./outputs/Correlation_tests/rho_uncorrected_pvalues_matrix.rds"), version = "2")


### 1.2/ Prepare things for heatmap plot ####

rho_matrix <- readRDS(file =  paste0("./outputs/Correlation_tests/rho_matrix.rds"))
rho_matrix_main <- rho_matrix[1:6, 1:6]

# Compute manually distance from the rho indices = 1 - |absolute rho|
rho_dist <- 1-abs(rho_matrix)
rho_dist_main <- 1-abs(rho_matrix_main)

library(dendextend)

# Compute the dendrogram outside of heatmap() and set branch width manually
dd <- set(as.dendrogram(hclust(as.dist(rho_dist), method = "complete")), "branches_lwd", 2)
dd_main <- set(as.dendrogram(hclust(as.dist(rho_dist_main), method = "complete")), "branches_lwd", 2)

# Reorder leaves of the dendrogram manually
dd_reorder <- reorder(x = dd, wts = rowMeans(x = rho_dist))
dd_reorder_main <- reorder(x = dd_main, wts = rowMeans(x = rho_dist_main))

plot(dd_reorder)
str(dd_reorder)
dd_reorder <-  reorder(dd_reorder, wts = c(0.114, 0.4, 0.118, 0.186, 0.368, 0.179, 0.113, 0.15, 0.121, 0.169))
plot(dd_reorder)

plot(dd_reorder_main)
dd_reorder_main <-  reorder(dd_reorder_main, wts = c(6:3, 1, 2), agglo.FUN = mean)
plot(dd_reorder_main)

# Set color palette
heatmap_pal <- rev(tmaptools::get_brewer_pal("RdYlBu", n = 100))[10:95]
# heatmap_pal_main <- rev(tmaptools::get_brewer_pal("RdYlBu", n = 100))[5:95]

heatmap_pal

test <- heatmap_pal[round(min(abs(rho_matrix))*length(heatmap_pal)):length(heatmap_pal)]
test <- heatmap_pal[round(min(abs(rho_matrix_main))*length(heatmap_pal)):length(heatmap_pal)]

### 1.3/ Plot heatmap ####
pdf(file = "./graphs/Correlation_tests/Heatmap_without_legend.pdf", height = 6, width = 6)
original_ext_margins <- par()$oma
par(oma = c(4,0,0,4), xpd = NA)

heatmap(abs(rho_matrix_main), symm = T, Rowv = rev(dd_reorder_main), Colv = rev(dd_reorder_main), revC = T, 
        col = heatmap_pal,
        add.expr = segments(x0 = 0, y0 = 0.6, x1 = 0, y1 = 6.4, lwd = 2, col = "red", lty = 2) +
          segments(x0 = 0.6, y0 = 6.95, x1 = 6.25, y1 = 6.95, lwd = 2, col = "red", lty = 2))

par(oma = original_ext_margins, xpd = FALSE)
dev.off()

pdf(file = "./graphs/Correlation_tests/Heatmap_without_legend_full.pdf", height = 6, width = 6)
original_ext_margins <- par()$oma
par(oma = c(2,0,0,3), xpd = NA)

heatmap(abs(rho_matrix), symm = T, Rowv = rev(dd_reorder), Colv = rev(dd_reorder), revC = T, 
        col = heatmap_pal,
        add.expr = segments(x0 = -0.4, y0 = 0.6, x1 = -0.4, y1 = 10.4, lwd = 2, col = "red", lty = 2) +
                   segments(x0 = 0.6, y0 = 11.5, x1 = 10.3, y1 = 11.5, lwd = 2, col = "red", lty = 2))

par(oma = original_ext_margins, xpd = FALSE)
dev.off()

?heatmap

# Get the scale in another pdf
min_scale <- round(min(abs(rho_matrix)), digits = 1)
max_scale <- round(max(abs(rho_matrix)), digits = 1)
scale_raster <- sp.richness
scale_raster[scale_raster < min_scale] <- min_scale
scale_raster[scale_raster > max_scale] <- max_scale

pdf(file = "./graphs/Correlation_tests/Legend_for_heatmap.pdf", height = 8, width = 8)
plot(scale_raster, col = heatmap_pal)
addRasterLegend(scale_raster, locs = seq(min_scale, 1, 0.1), cex.axis = 1.2, ramp = heatmap_pal, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 5.5, font = 2, cex = 1.4, label = bquote("|"~rho~"|"))
dev.off()


# # With another package for more options to make the plot pretty, including a legend
# 
# library(pheatmap)
# 
# ?pheatmap
# 
# pheatmap(mat = abs(rho_matrix))


### 1.4/ Plot pairwise plots of index correlation ####
pairs(all_indices_stack)

pdf(file = "./graphs/Correlation_tests/Pairs_graph.pdf", height = 12, width = 12)
pairs(all_indices_stack)
dev.off()

# Plot pair graphs without predicts
plot(sp.mean.rarity[] ~ sp.richness[],
     ylab = "Species mean rarity",
     xlab = "Species richness",
     pch = 16)

# Build a tidy df for ggplot2
library(tidyverse)

n <- length(sp.richness[])
indices_df <- tibble(sp.richness = sp.richness[], sp.mean.rarity = sp.mean.rarity[], MPD = MPD[], ring.richness = ring.richness[], mimicry.mean.rarity = mimicry.mean.rarity[])

save(indices_df, file = "./outputs/Correlation_tests/Indices_df.RData", version = 2)
saveRDS(indices_df, file = "./outputs/Correlation_tests/Indices_df.rds", version = 2)

indices_df <- readRDS(file = "./outputs/Correlation_tests/Indices_df.rds")

indices_df <- indices_df %>% drop_na()

g1 <- ggplot(data = indices_df, aes(x = sp.richness, y = sp.mean.rarity)) +
  geom_point() +
  geom_smooth(color = "red", size = 1.5, se = F) +
  annotate("text", x = max(indices_df$sp.richness), y = min(indices_df$sp.mean.rarity), hjust = 1, vjust = 0, label = "(a)", size = 7, fontface = 2) +
  annotate("text", x = max(indices_df$sp.richness), y = max(indices_df$sp.mean.rarity) + 0.035, hjust = 0.95, vjust = 1, label = bquote(rho ~ "= 0.47"), size = 5, fontface = 2) +
  # annotate("rect", xmin = 135, xmax = 145, ymin = -0.02, ymax = 0.05, fill = "white", col = NA) +
  # annotate("text", label = "A", x = 140, y = 0.02, size = 7, colour = "black", fontface = 2) +
  ylab(label = "Species geographic rarity \n[Rarity index]") +
  xlab(label = "Species richness\n[Species]") +
  scale_x_continuous(breaks = seq(0, 125, 25)) +  # Ticks from 0-10, every .25 +
  scale_y_continuous(breaks = seq(0, 0.6, 0.2), limits = c(0, 0.6)) + 
  theme_gray() +
  theme(panel.grid.minor.x=element_blank(),
        panel.border = element_rect(fill = NA),
        axis.ticks = element_line(size = 1.2),
        axis.ticks.length = unit(5, "pt"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title.y = element_text(margin = margin(l = 5, r = 10)))

print(g1)


g2 <- ggplot(data = indices_df, aes(x = sp.richness, y = MPD)) +
  geom_point() +
  geom_smooth(color = "red", size = 1.5, se = F) +
  annotate("text", x = max(indices_df$sp.richness), y = min(indices_df$MPD), hjust = 1, vjust = 0, label = "(b)", size = 7, fontface = 2) +
  annotate("text", x = max(indices_df$sp.richness), y = max(indices_df$MPD), hjust = 0.9, vjust = 0.8, label = bquote(rho ~ "= 0.05"), size = 5, fontface = 2) +
  # annotate("rect", xmin = 135, xmax = 145, ymin = 36, ymax = 37, fill = "white", col = NA) +
  # annotate("text", label = "A", x = 140, y = 36.5, size = 7, colour = "black", fontface = 2) +
  ylab(label = "Mean pairwise Phylogenetic Distance\n[Evolutionary time in My]") +
  xlab(label = "Species richness\n[Species]") +
  scale_x_continuous(breaks = seq(0, 125, 25)) +  # Ticks from 0-10, every .25 +
  theme_gray() +
  theme(panel.grid.minor.x=element_blank(),
        panel.border = element_rect(fill = NA),
        axis.ticks = element_line(size = 1.2),
        axis.ticks.length = unit(5, "pt"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title.y = element_text(margin = margin(l = 5, r = 10)))

print(g2)

g3 <- ggplot(data = indices_df, aes(x = sp.richness, y = ring.richness)) +
  geom_point() +
  geom_smooth(color = "red", size = 1.5, se = F) +
  annotate("text", x = max(indices_df$sp.richness), y = min(indices_df$ring.richness), hjust = 1, vjust = 0, label = "(c)", size = 7, fontface = 2) +
  annotate("text", x = min(indices_df$sp.richness), y = max(indices_df$ring.richness), hjust = 0, vjust = 1, label = bquote(rho ~ "= 0.93"), size = 5, fontface = 2) +
  ylab(label = "Mimicry richness\n[Mimicry rings]") +
  xlab(label = "Species richness\n[Species]") +
  scale_x_continuous(breaks = seq(0, 125, 25)) +
  theme_gray() +
  theme(panel.grid.minor.x = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.ticks = element_line(size = 1.2),
        axis.ticks.length = unit(5, "pt"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title.y = element_text(margin = margin(l = 5, r = 10)))

print(g3)

g4 <- ggplot(data = indices_df, aes(x = ring.richness, y = mimicry.mean.rarity)) +
  geom_point() +
  geom_smooth(color = "red", size = 1.5, se = F) +
  annotate("text", x = max(indices_df$ring.richness), y = min(indices_df$mimicry.mean.rarity), hjust = 1, vjust = 0, label = "(d)", size = 7, fontface = 2) +
  annotate("text", x = max(indices_df$ring.richness), y = max(indices_df$mimicry.mean.rarity) + 0.055, hjust = 0.95, vjust = 0, label = bquote(rho ~ "= 0.65"), size = 5, fontface = 2) +
  ylab(label = "Mimicry geographic rarity\n[Rarity index]") +
  xlab(label = "Mimicry richness\n[Mimicry rings]") +
  scale_x_continuous(breaks = seq(0, 30, 5)) +  
  scale_y_continuous(breaks = seq(0, 0.75, 0.25), limits = c(0, 0.75)) +  
  theme_gray() +
  theme(panel.grid.minor.x = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.ticks = element_line(size = 1.2),
        axis.ticks.length = unit(5, "pt"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title.y = element_text(margin = margin(l = 5, r = 10)))

print(g4)

g5 <- ggplot(data = indices_df, aes(x = sp.richness, y = mimicry.mean.rarity)) +
  geom_point() +
  geom_smooth(color = "red", size = 1.5, se = F) +
  annotate("text", x = max(indices_df$sp.richness), y = min(indices_df$mimicry.mean.rarity), hjust = 1, vjust = 0, label = "(b)", size = 7, fontface = 2) +
  annotate("text", x = max(indices_df$sp.richness), y = max(indices_df$mimicry.mean.rarity) + 0.055, hjust = 0.95, vjust = -0.1, label = bquote(rho ~ "= 0.61"), size = 5, fontface = 2) +
  ylab(label = "Mimicry geographic rarity\n[Rarity index]") +
  xlab(label = "Species richness\n[Species]") +
  scale_x_continuous(breaks = seq(0, 125, 25)) + 
  scale_y_continuous(breaks = seq(0, 0.75, 0.25), limits = c(0, 0.75)) +  
  theme_gray() +
  theme(panel.grid.minor.x = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.ticks = element_line(size = 1.2),
        axis.ticks.length = unit(5, "pt"),
        axis.text = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 14, face = "bold"),
        axis.title.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title.y = element_text(margin = margin(l = 5, r = 10)))

print(g5)

library(gridExtra)

# Plot the selected pairwise plots. With MPD ~ Species richness

pdf(file = "./graphs/Correlation_tests/Selected_pairs_graph_smoothed_MPD.pdf", height = 10, width = 10)
# pdf(file = "./graphs/Correlation_tests/Selected_pairs_graph.pdf", height = 10, width = 10)
grid.arrange(g1, g2, g3, g4,
  widths = c(10, 0.4, 10),  # Width of columns
  heights = c(10, 0.2, 10),
  nrow = 3,
  ncol = 3,
  layout_matrix = rbind(c(1, NA, 2),  # Position of ggplots in the layout
                        c(NA, NA, NA),
                        c(3, NA, 4)))
dev.off()

# Plot the selected pairwise plots. With Mimicry rarity ~ Species richness

pdf(file = "./graphs/Correlation_tests/Selected_pairs_graph_smoothed.pdf", height = 10, width = 10)
# pdf(file = "./graphs/Correlation_tests/Selected_pairs_graph.pdf", height = 10, width = 10)
grid.arrange(g1, g5, g3, g4,
             widths = c(10, 0.4, 10),  # Width of columns
             heights = c(10, 0.2, 10),
             nrow = 3,
             ncol = 3,
             layout_matrix = rbind(c(1, NA, 2),  # Position of ggplots in the layout
                                   c(NA, NA, NA),
                                   c(3, NA, 4)))
dev.off()

###### 2/ Corrected tests accounting for spatial autocorrelation ######

### 2.1/ Pearson's rho corrected test version ####

# Clifford correction: see Chapter 11.2 in Spatial Data Analysis in Ecology and Agriculture Using R

### 2.1.1/ modified.ttest function in SpatialPack

# Version pour tester un Pearson's r avec des corrections sur l'effectif size (et les ddl) et la variance de r_hat (utilisée pour estimer l'effective size).

library(SpatialPack)

# Fonction toute faite pour un Pearson's r
# Problem = pas adaptable pour un Spearman car le code source est compilé sous C !

?modified.ttest
modified.ttest(x = sample.sp.richness, y = sample.ring.richness, coords = sample.all_coords, nclass = NULL)

# Non-modified version
cor.test(sample.sp.richness, sample.ring.richness)

### 2.1.1/ Manual implementation of Clifford's correction for Pearson's r and Spearman's rho

# Only difference for Spearman = compute the autocorrelation coefficients on ranks rather than values to estimate the effective sample size,
# and use the Spearman's rho stat instead of the Pearson's r (also follow a t-stats with dl = N'- 2, but computed differently)

## Simulate data
# library(spdep)
# 
# # Set a level of spatial autocorrelation
# rho <- 0.6
# 
# # Generate a list of neighbors from a matrix
# nlist <- cell2nb(14, 14)
# 
# # Compute the matrix used to generate simultaneous autoregressive random variable from a list of neighbor and a given spatial autocorrelation value (rho) then generate them
# IrWinv <- invIrM(nlist, rho)
# set.seed(123)
# Y1 <- IrWinv %*% rnorm(14^2)
# Y2 <- IrWinv %*% rnorm(14^2)
# 
# # Remove the edge cells to avoid edge effects
# Y1samp <- matrix(Y1, nrow = 14, byrow = TRUE)[3:12,3:12]
# Y2samp <- matrix(Y2, nrow = 14, byrow = TRUE)[3:12,3:12]

## Or use indices data
Y1samp <- sample_indices_df$Sp..richness
Y1samp <- sample_indices_df$Sp..diversity
Y1samp <- sample_indices_df$Sp..rarity
Y1samp <- sample_indices_df$Faith.s.PD
Y1samp <- sample_indices_df$Sum.of.ED
Y1samp <- sample_indices_df$MPD
Y1samp <- sample_indices_df$Mim..richness
Y1samp <- sample_indices_df$Mim..diversity
Y1samp <- sample_indices_df$Mim..rarity

Y1samp <- sample_indices_df$Sp..richness
Y2samp <- sample_indices_df$Mim..richness

## Compute t-test without correction for spatial autocorrelation

# Pearson's r
cortest <- cor.test(Y1samp, Y2samp, alternative = "two.sided", method = "pearson")
cortest
r <- cortest$estimate # Extract the Pearson's R
print(t.uncorr <- cortest$statistic , digits = 3) #  Uncorrected t-stat for Pearson's r
print(t.uncorr <- r * sqrt(N - 2) / sqrt(1 - r^2), digits = 3) #  Uncorrected t-stat for Pearson's r # Manual computation

# Spearman's rho
cortest_sp <- cor.test(Y1samp, Y2samp, alternative = "two.sided", method = "spearman")
cortest_sp
r_sp <- cortest_sp$estimate # Extract the Spearman's rho
print(t.uncorr_sp <- r_sp * sqrt(N - 2) / sqrt(1 - r_sp^2), digits = 3) #  Uncorrected t-stat for Spearman's rho


## Clifford's correction (1989)
library(spatial)

# Get coordinates of points
# coords.xy <- expand.grid(1:10,10:1)
coords.xy <- sample.all_coords

# Extract values in a vector (useless for indices which ar already in vector. Only for simulated data stored in matrices)
Y1.vec <- as.vector(t(Y1samp))
Y2.vec <- as.vector(t(Y2samp))
N <- length(Y1.vec) # Store nb of observations = uncorrected sample size

# Compute ranks for Spearman's test
Y1.rank <- rank(Y1.vec)
Y2.rank <- rank(Y2.vec)

# Transform points into a surface for Pearson's
Y1.surf <- surf.ls(0,coords.xy[,1], coords.xy[,2], Y1.vec)
Y2.surf <- surf.ls(0,coords.xy[,1], coords.xy[,2], Y2.vec)

# Transform points into a surface for Spearman's
Y1.surf_sp <- surf.ls(0,coords.xy[,1], coords.xy[,2], Y1.rank)
Y2.surf_sp <- surf.ls(0,coords.xy[,1], coords.xy[,2], Y2.rank)

# Compute correlograms for Ht distance classes based on a trend surface
?correlogram

Ht <- 70 # Choose the number of classes arbitrary (make classes of around 1° width)

# Ht <- as.integer(1.5 + 3.3 * log10(N)) # Use Struges's formula to define the number of classes

# # Use distance (class witdh) to define nb of classes
# dist_bandwidth <- 1
# max_dist <- max(dist(coords.xy))
# Ht <- round(max_dist/dist_bandwidth , 0)

# For Pearson's r, based on values
r1 <- correlogram(Y1.surf, Ht)
r2 <- correlogram(Y2.surf, Ht)

# For Spearman's rho, based on ranks
r1_sp <- correlogram(Y1.surf_sp, Ht)
r2_sp <- correlogram(Y2.surf_sp, Ht)

H <- 18 # Choose the max distance class to look at spatial autocorrelation depending on the plot fo the correlograms
# After a certain distance, autocorrelation should be null, thus no need to include the next classes in the computation

# Get automatically the first class with negative autocorrelation, for this class of smaller, for both indices 
# H <- max(which.max(r1$y < 0), which.max(r2$y < 0))
# H <- max(which.max(r1_sp$y < 0), which.max(r2_sp$y < 0))

# Record distance of absence of spatial autocorrelation (in decimal degrees)
print(threshold_dist <- r1$x[H], digits = 4)
print(threshold_dist <- r1_sp$x[H], digits = 4)


## Estimate the variance of the sampling distribution (sampling error) of Pearson's r from the first H classes where autocorrelation is still present 

# Pearson's r
nr1r2 <- r1$cnt * r1$y * r2$y
print(sr.hat <- sum(nr1r2[1:H]) * N^(-2), digits = 3)

# Spearman's rho
nr1r2_sp <- r1_sp$cnt * r1_sp$y * r2_sp$y
print(sr.hat_sp <- sum(nr1r2_sp[1:H]) * N^(-2), digits = 3)

## Compute corrected sample size
print(ne.hat <- 1 + 1 / sr.hat, digits = 3) # Pearson's r
print(ne.hat_sp <- 1 + 1 / sr.hat_sp, digits = 3) # Spearman's rho

## Compute the corrected t-test with corrected sample size for the df

# Pearson's r
print(t.corr <- r * sqrt(ne.hat - 2) / sqrt(1 - r^2), digits = 3) #  Corrected t-stat
print(p.corr <- 2 * (1 - pt(q = abs(t.corr), df = ne.hat - 2)), digits = 3) # Corrected p-value

# Spearman's rho
print(t.corr_sp <- r_sp * sqrt(ne.hat_sp - 2) / sqrt(1 - r_sp^2), digits = 3) #  Corrected t-stat
print(p.corr_sp <- 2 * (1 - pt(q = abs(t.corr_sp), df = ne.hat_sp - 2)), digits = 3) # Corrected p-value


### 2.1.2/ Automatic script for Clifford's corrected t-test for Spearman's rho ####

library(spatial)

N.indices <- ncol(sample_indices_df)
coords.xy <- sample.all_coords # Get coordinates of points

## Compute autocorrelogram of each variables and store autocorrelation coefficients used to estimate the variance of the sampling distribution of Spearman's rho, and effective sample size afterward
autocorr_coefs <- list(NULL)
for (i in 1:N.indices)
{
  # Get the indices values
  Y <- sample_indices_df[ ,i]
  
  # Compute ranks
  Y.rank <- rank(Y)
  
  # Transform points into a surface
  Y.surf_sp <- surf.ls(0, coords.xy[,1], coords.xy[,2], Y.rank)
  
  # Compute correlograms for Ht distance classes based on a trend surface
  Ht <- 70 # Choose the number of classes arbitrary (make classes of around 1° width)
  r_sp <- correlogram(Y.surf_sp, Ht)
  
  # Store autocorrelation coefficients
  autocorr_coefs[[i]] <- r_sp$y
  
  # Get counts of pairs by distance classes
  pairs_cnt <- r_sp$cnt
  
  print(i)
}

## Build a matrix of pairwise rho coefficients and associated corrected p-values
rho_matrix <- rho_corrected_pvalues_matrix <- rho_corrected_tstat_matrix <- rho_corrected_df_matrix <- matrix(nrow = N.indices, ncol = N.indices)
for (i in 1:N.indices)
{
  # for (j in (i+1):v) # Most efficient version with no useless double tests
  for (j in 1:N.indices) # Version to fill the all matrix
  {
    # Get the indices values
    Y1 <- sample_indices_df[ ,i]
    Y2 <- sample_indices_df[ ,j]
    N <- length(Y1) # Store nb of observations = uncorrected sample size
    
    # Get the rho coefficient form the uncorrected t-test function
    test <- cor.test(x = Y1, y = Y2, method = "spearman", alternative = "two.sided", na.rm = T)
    rho_matrix[i,j] <- r_sp <- round(test$estimate, 3)
    
    # Choose the max distance class to look at spatial autocorrelation depending on the plot fo the correlograms
    # After a certain distance, autocorrelation should be null, thus no need to include the next classes in the computation
    H <- 18 
    
    # Estimate the variance of the sampling distribution (sampling error) of Spearman's rho from the first H classes where autocorrelation is still present 
    nr1r2_sp <- pairs_cnt * autocorr_coefs[[i]] * autocorr_coefs[[j]]
    sr.hat_sp <- sum(nr1r2_sp[1:H]) / N^2
    
    # Compute corrected sample size
    ne.hat_sp <- 1 + 1 / sr.hat_sp
    
    # Compute the corrected t-test with corrected sample size for the df
    t.corr_sp <- r_sp * sqrt(ne.hat_sp - 2) / sqrt(1 - r_sp^2) #  Corrected t-stat
    p.corr_sp <- 2 * (1 - pt(q = abs(t.corr_sp), df = ne.hat_sp - 2)) # Corrected p-value
    
    # Save the corrected df
    rho_corrected_df_matrix[i,j] <- round((ne.hat_sp - 2), 1)
    # Save the t-stat
    rho_corrected_tstat_matrix[i,j] <- round(t.corr_sp, 2)
    # Save the p-value
    rho_corrected_pvalues_matrix[i,j] <- round(p.corr_sp, 3)
    
  }
  
  print(i)
}

rownames(rho_matrix) <- colnames(rho_matrix) <- indices_list
rownames(rho_corrected_pvalues_matrix) <- colnames(rho_corrected_pvalues_matrix) <- indices_list
rownames(rho_corrected_tstat_matrix) <- colnames(rho_corrected_tstat_matrix) <- indices_list
rownames(rho_corrected_df_matrix) <- colnames(rho_corrected_df_matrix) <- indices_list

rho_matrix
rho_corrected_pvalues_matrix
rho_corrected_tstat_matrix
rho_corrected_df_matrix

# Get the order of indices used in the full dendrogram
hclust_reorder <- as.hclust(dd_reorder)
str(hclust_reorder)
rev(hclust_reorder$labels[hclust_reorder$order])
indices_order <- rev(hclust_reorder$order)

# Reorder all matrix to follow the order of indices in the full dendrogram
rho_matrix <- rho_matrix[indices_order, indices_order]
rho_corrected_pvalues_matrix <- rho_corrected_pvalues_matrix[indices_order, indices_order]
rho_corrected_tstat_matrix <- rho_corrected_tstat_matrix[indices_order, indices_order]
rho_corrected_df_matrix <- rho_corrected_df_matrix[indices_order, indices_order]

# Save final test outputs
saveRDS(rho_matrix, file =  paste0("./outputs/Correlation_tests/rho_matrix.rds"), version = "2")
saveRDS(rho_corrected_pvalues_matrix, file =  paste0("./outputs/Correlation_tests/rho_corrected_pvalues_matrix.rds"), version = "2")
saveRDS(rho_corrected_tstat_matrix, file =  paste0("./outputs/Correlation_tests/rho_corrected_tstat_matrix.rds"), version = "2")
saveRDS(rho_corrected_df_matrix, file =  paste0("./outputs/Correlation_tests/rho_corrected_df_matrix.rds"), version = "2")

### Generate table with rho in lower triangle and p-value in upper triangle
rho_pvalue_matrix <- matrix(nrow = N.indices, ncol = N.indices)

lower.tri.mask <- lower.tri(rho_pvalue_matrix) # Generate a mask to extract and fill only lower triangle
upper.tri.mask <- upper.tri(rho_pvalue_matrix) # Generate a mask to extract and fill only upper triangle

rho_pvalue_matrix[lower.tri.mask] <- rho_matrix[lower.tri.mask] # Extract and fill lower triangle with rho indices
rho_pvalue_matrix[upper.tri.mask] <- rho_corrected_pvalues_matrix[upper.tri.mask]  # Extract and fill lower triangle with corrected p-values

rownames(rho_pvalue_matrix) <- colnames(rho_pvalue_matrix) <- indices_list[indices_order] # Provide names for rows and columns

rho_pvalue_matrix

# Save/export result
saveRDS(rho_pvalue_matrix, file =  paste0("./outputs/Correlation_tests/rho_pvalue_matrix.rds"), version = "2")
write.csv2(x = rho_pvalue_matrix, file = "./outputs/Correlation_tests/rho_pvalue_matrix.csv")

### Generate table with t-stat in lower triangle and corrected df in upper triangle
tstat_df_matrix <- matrix(nrow = N.indices, ncol = N.indices)

lower.tri.mask <- lower.tri(tstat_df_matrix) # Generate a mask to extract and fill only lower triangle
upper.tri.mask <- upper.tri(tstat_df_matrix) # Generate a mask to extract and fill only upper triangle

tstat_df_matrix[lower.tri.mask] <- rho_corrected_tstat_matrix[lower.tri.mask] # Extract and fill lower triangle with rho indices
tstat_df_matrix[upper.tri.mask] <- rho_corrected_df_matrix[upper.tri.mask]  # Extract and fill lower triangle with corrected p-values

rownames(tstat_df_matrix) <- colnames(tstat_df_matrix) <- indices_list[indices_order] # Provide names for rows and columns

tstat_df_matrix

# Save/export result
saveRDS(tstat_df_matrix, file =  paste0("./outputs/Correlation_tests/tstat_df_matrix.rds"), version = "2")
write.csv2(x = tstat_df_matrix, file = "./outputs/Correlation_tests/tstat_df_matrix.csv")


### 3/ Tjostheim spatial non-parametric correlation coefficient based on ranks and distances in space between observations of close ranks. ####
# In practice, it need to specify a distance function to compute distances from coordinates, and it includes suitable modifications of the coordinates (???)

library(SpatialPack)

?cor.spatial
SpatialPack::cor.spatial(x = sample_indices_df$Sp..richness, y = sample_indices_df$Mim..richness, coords = sample.all_coords)

### 3.1/ Build a matrix of pairwise rho coefficients and associated p-values ####
N.indices <- ncol(sample_indices_df)

Tjostheim_matrix <- matrix(nrow = N.indices, ncol = N.indices)
for (i in 1:N.indices)
{
  # for (j in (i+1):v) # Most efficient version with no useless double tests
  for (j in 1:N.indices) # Version to fill the all matrix
  {
    # Get the indices values
    Y1 <- sample_indices_df[ ,i]
    Y2 <- sample_indices_df[ ,j]
    
    # Compute the indices and store them
    Tjostheim_matrix[i,j] <- SpatialPack::cor.spatial(x = Y1, y = Y2, coords = sample.all_coords)
  }
}

rownames(Tjostheim_matrix) <- colnames(Tjostheim_matrix) <- indices_list

### 3.2/ Prepare things for heatmap plot ####

# Compute manually distance from the rho indices = 1 - |absolute rho|
Tjostheim_dist <- 1-abs(Tjostheim_matrix)

library(dendextend)

# Compute the dendrogram outside of heatmap() and set branch width manually
dd <- set(as.dendrogram(hclust(as.dist(Tjostheim_dist), method = "complete")), "branches_lwd", 2)

# Reorder leaves of the dendrogram manually
dd_reorder <- reorder(x = dd, wts = rowMeans(x = Tjostheim_dist))

plot(dd_reorder)

# Set color palette
heatmap_pal <- rev(tmaptools::get_brewer_pal("RdYlBu", n = 100))[10:95]


### 3.3/ Plot heatmap ####
pdf(file = "./graphs/Correlation_tests/Heatmap_Tjostheim_without_legend.pdf", height = 6, width = 6)
original_ext_margins <- par()$oma
par(oma = c(2,0,0,4))

heatmap(abs(Tjostheim_matrix), symm = T, Rowv = rev(dd_reorder), Colv = rev(dd_reorder), revC = T, col = heatmap_pal)

par(oma = original_ext_margins)
dev.off()

?heatmap

# Get the scale in another pdf
scale_raster <- sp.richness > 0
# scale_raster[1] <- -1 # For centered scale

pdf(file = "./graphs/Correlation_tests/Legend_for_heatmap.pdf", height = 8, width = 8)
plot(scale_raster, col = heatmap_pal)
addRasterLegend(scale_raster, locs = seq(0, 1, 0.2), cex.axis = 1.2, ramp = heatmap_pal, ncolors = 200, border = T, location = c(-115, -112, -33, 0))
graphics::text(x = -112, y = 5.5, font = 2, cex = 1.4, label = bquote("|"~rho~"|"))
dev.off()


# Set color palette to display lower values for -1 and max value for 1
min <- min(Tjostheim_matrix)
N_col_full <- length(heatmap_pal)
N_col <-  round((1 - min)/2 * N_col_full, 0)
# Use only portion of the color palette representing the portion of values range
heatmap_pal_centered_scale <- heatmap_pal[(N_col_full-N_col):N_col_full]

heatmap(Tjostheim_matrix, symm = T, Rowv = rev(dd_reorder), Colv = rev(dd_reorder), revC = T, col = heatmap_pal_centered_scale)
# Limit : values in the heatmap are not the one used to build the dendrogram...

# Get the centered scale in another pdf
scale_raster <- sp.richness > 0
scale_raster[1] <- -1 # For centered scale

pdf(file = "./graphs/Correlation_tests/Legend_for_heatmap_centered.pdf", height = 8, width = 8)
plot(scale_raster, col = heatmap_pal)
addRasterLegend(scale_raster, locs = seq(-1, 1, 0.5), cex.axis = 1.2, digits = 1, ramp = heatmap_pal, ncolors = 200, border = T, location = c(-110, -107, -33, 0))
graphics::text(x = -107, y = 6.5, font = 2, cex = 1.3, label = "Spatial\ncorrelation")
dev.off()


# # With another package for more options to make the plot pretty, including a legend
# 
# library(pheatmap)
# 
# ?pheatmap
# 
# pheatmap(mat = abs(rho_matrix))



### Extraction of coordinates inside Ithomiini range from raster

sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/sp.richness_Jaccard.80.rds"))

Ithomiini_range_mask <- sp.richness_Jaccard.80 > 0


all_coords <- coordinates(Ithomiini_range_mask)

### Choose subsampling index

N.sample = 1000

set.seed <- 45512

sample.all_index <- which(Ithomiini_range_mask[] == 1) # Version with all communities
sample.all_index <- sample(x = which(Ithomiini_range_mask[] == 1), size = N.sample, replace = F)
sample.all_coords <- all_coords[sample.all_index,]

plot(Ithomiini_range_mask)
points(sample.all_coords, pch = 16, cex = 0.5, add = T)

save(sample.all_index, file =  paste0("./outputs/Correlation_tests/sample.all_index.RData"), version = "2")
saveRDS(sample.all_index, file =  paste0("./outputs/Correlation_tests/sample.all_index.rds"), version = "2")
save(sample.all_coords, file =  paste0("./outputs/Correlation_tests/sample.all_coords.RData"), version = "2")
saveRDS(sample.all_coords, file =  paste0("./outputs/Correlation_tests/sample.all_coords.rds"), version = "2")


### Project coords on plan to use in nlme corStructure functions

sample.all_coords_Spobj <- SpatialPoints(sample.all_coords, proj4string = Ithomiini_range_mask@crs)
sample.all_coords_Spobj <- spTransform(sample.all_coords_Spobj, CRS("+proj=utm +zone=20 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")) 
proj_sample.all_coords <- sample.all_coords_Spobj@coords

save(proj_sample.all_coords, file =  paste0("./outputs/Correlation_tests/proj_sample.all_coords.RData"), version = "2")
saveRDS(proj_sample.all_coords, file =  paste0("./outputs/Correlation_tests/proj_sample.all_coords.rds"), version = "2")


### Compute geographic distance (if find a way to use it in correlation strucutre or autocorrelogram)
?distm
library(geosphere)

# For autocorrelation
dist_all <- distm(sample.all_coords, fun = distGeo) # Génère la matrice symétrique de distance par calcul trigo.
dist_all <- as.dist(dist_all/1000) # Transforme la matrice en matrice triangulaire inférieure pour éviter les doublons (n(n-1)/2 éléments). # Conversion en km
save(dist_all, file = paste0("./outputs/Correlation_tests/dist_all.RData"), version = "2")
saveRDS(dist_all, file = paste0("./outputs/Correlation_tests/dist_all.rds"), version = "2")


##### Approche GLS avec structure de covariance modélisée (Devictor et al., 2014 => corExp)
# Besoin de subsampling sinon la matrice de distance est trop fat...

# Load maps
sp.richness <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
sp.diversity.compatible_Jost_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/sp.diversity.compatible_Jost_Jaccard.80.rds"))
sp.mean.rarity_Jaccard.80 <- readRDS(file = "./outputs/Indices_maps/sp.mean.rarity_Jaccard.80.rds")
PD.raster_Jaccard.80 <- readRDS(file = "./outputs/Indices_Maps/PD.raster_Jaccard.80.rds")
ED.raster_Jaccard.80 <- readRDS(file = "./outputs/Indices_Maps/ED.raster_Jaccard.80.rds")
MPD.raster_Jaccard.80_contrasted <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80_contrasted.rds")
ring.richness <- readRDS(file = paste0("./outputs/Indices_maps/ring.richness_Jaccard.80.rds"))
ring.diversity.compatible_Jost_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/ring.diversity.compatible_Jost_Jaccard.80.rds"))
mimicry.mean.rarity_Jaccard.80 <- readRDS(file = "./outputs/Indices_maps/mimicry.mean.rarity_Jaccard.80.rds")


# Extract the samples
sample.sp.richness <- sp.richness[sample.all_index]
sample.ring.richness <- ring.richness[sample.all_index]

library(nlme)

?gls
?corClasses
?corSpatial

?corExp # Relation exponentielle (courbe log)
?corLin # Relation linéaire (droite)
?corGaus # Relation sinusoidale (telle la cumulative df d'une loi Normale)
?corSpher # Relation exponentielle avec plateau
?corRatio # ???


### Version quadratique

m0 <- gls(sample.ring.richness ~ sample.sp.richness + I(sample.sp.richness^2))
summary(m0)
anova(m0, type = "marginal")
resid.model <- residuals(m0, type = "response")
?Variogram
vario0 <- Variogram(resid.model, distance = dist_all, resType = "pearson", maxDist = 3000)
?plot.Variogram
plot(vario0, smooth = TRUE, ylim = c(0, 1.2))

load(file =  paste0(internal.wd, "/Maps_exploitation/proj_sample_coords.RData"))

data.gls <- as.data.frame(cbind(sample.EDD, sample.sp.richness))

for (i in 1) {
  m1 <- gls(data = data.gls, sample.EDD ~ sample.sp.richness + I(sample.sp.richness^2), correlation = corExp(form = ~ proj_sample.all_coords[,1] + proj_sample.all_coords[,2], nugget = TRUE))
  print("m1")
  m2 <- gls(data = data.gls, sample.EDD ~ sample.sp.richness + I(sample.sp.richness^2), correlation = corGaus(form = ~ proj_sample.all_coords[,1] + proj_sample.all_coords[,2], nugget = TRUE))
  print("m2")
  m3 <- gls(data = data.gls, sample.EDD ~ sample.sp.richness + I(sample.sp.richness^2), correlation = corSpher(form = ~ proj_sample.all_coords[,1] + proj_sample.all_coords[,2], nugget = TRUE))
  print("m3")
  m4 <- gls(data = data.gls, sample.EDD ~ sample.sp.richness + I(sample.sp.richness^2), correlation = corLin(form = ~ proj_sample.all_coords[,1] + proj_sample.all_coords[,2], nugget = TRUE))
  print("m4")
  m5 <- gls(data = data.gls, sample.EDD ~ sample.sp.richness + I(sample.sp.richness^2), correlation = corRatio(form = ~ proj_sample.all_coords[,1] + proj_sample.all_coords[,2], nugget = TRUE))
  print("m5")
}
save(m0, m1, m2, m3, m4, m5, file =  paste0(internal.wd, "/Maps_exploitation/Models/Models_EDD_vs_sp.richness_quad.RData"))
load(file =  paste0(internal.wd, "/Maps_exploitation/Models/Models_EDD_vs_sp.richness_quad.RData"))

# Find the best model
AIC(m0, m1, m2, m3, m4, m5)
which.min(AIC(m0, m1, m2, m3, m4, m5)[,2])

Cor_EDD_Sp.richness_AIC_df <- data.frame(c("No structure", "corExp", "corGaus", "corSpher", "corLin", "corRatio")) ; names(Cor_EDD_Sp.richness_AIC_df) <- "Correlation Structure"
Cor_EDD_Sp.richness_AIC_df$AIC <- AIC(m0, m1, m2, m3, m4, m5)
save(Cor_EDD_Sp.richness_AIC_df, file =  paste0(internal.wd, "/Maps_exploitation/Models/Cor_EDD_Sp.richness_AIC_df.RData"))
load(file =  paste0(internal.wd, "/Maps_exploitation/Models/Cor_EDD_Sp.richness_AIC_df.RData"))

?Variogram
# Check the fit of the model to the semi-variogram
vario_best <- Variogram(m1, form = ~ proj_sample.all_coords[,1] + proj_sample.all_coords[,2], resType = "pearson", maxDist = 4000)
plot(vario_best, smooth = FALSE, ylim = c(0, 1.2))

# Plot Variogram on normalized residuals that should show no trend if things went well
vario_norm <- Variogram(m1, form = ~ proj_sample.all_coords[,1] + proj_sample.all_coords[,2], resType = "normalized", maxDist = 3000)
plot(vario_norm, smooth = FALSE, ylim = c(0, 1.2))

### Plot the correlogram of residuals before and after taking spatial autocorrelation into consideration

library(ncf)

resid.model <- resid(m0, type = "response") # Sans prise en compte de la strucutre spatiale
resid.model <- resid(m1, type = "response") # Avec prise en compte de la strucutre spatiale

hist(resid.model)

# bins <- round(x= 1 + 3.3*log(x = max(dist), base = 2), digits = 0) # Nombre de bins selon la règle de Struges
# increment <- round(max(dist)/bins,0) ; increment
increment <- 200

?correlog

correlogram <- ncf::correlog(x = proj_sample.all_coords[,1], y = proj_sample.all_coords[,2], z = resid.model, latlon=F, increment = increment, resamp=1000) # Génère le corrélogramme en calculant des distance de grand cercle
par(mfcol=c(2,2)) ; plot(correlogram_m1)  ; abline(h=0)
plot(correlogram$mean.of.class,correlogram$correlation, type="b", col = c("black","red")[(correlogram$p>0.05)+1], ylim = c(-1,1), pch=16, cex=1, lwd=1.5, xlab="Distance (km)", ylab="Moran's I", cex.lab=1.5, cex.axis=1.2) ; abline(h=0)
plot(correlogram$mean.of.class[1:15],correlogram$correlation[1:15], type="b", col = c("black","red")[(correlogram$p>0.05)+1], ylim = c(-1,1), pch=16, cex=1, lwd=1.5, xlab="Distance (km)", ylab="Moran's I", cex.lab=1.5, cex.axis=1.2) ; abline(h=0)
plot(correlogram$mean.of.class[1:10],correlogram$correlation[1:10], type="b", col = c("black","red")[(correlogram$p>0.05)+1], ylim = c(-1,1), pch=16, cex=1, lwd=1.5, xlab="Distance (km)", ylab="Moran's I", cex.lab=1.5, cex.axis=1.2) ; abline(h=0)

# correlogram_m0 <- correlogram
save(correlogram_m0, file = paste0(internal.wd,"/Maps_exploitation/Correlog/Correlog_SP.richness_EDD_Classic_LM.RData"))
# correlogram_m1 <- correlogram
save(correlogram_m1, file = paste0(internal.wd,"/Maps_exploitation/Correlog/Correlog_SP.richness_EDD_GLS_Spatial.RData"))

# Plot Correlog Classic LM
jpeg(filename = paste0(internal.wd,"/Maps_exploitation/Correlog/Correlog_SP.richness_EDD_Classic_LM.jpeg"), quality =100)
par(mfcol=c(1,1))
plot(main = "Correlogram of residuals \n Species Richness ~ Faith's PD \n (Classic LM)", x = correlogram_m0$mean.of.class[1:20], y = correlogram_m0$correlation[1:20], type="b", col = c("black","red")[(correlogram_m0$p>0.05)+1], ylim = c(-1,1), pch=16, cex=1, lwd=1.5, xlab="Distance (km)", ylab="Moran's I", cex.lab=1.5, cex.axis=1.2) ; abline(h=0)
dev.off()

# Plot Correlog GLS Spatial
jpeg(filename = paste0(internal.wd,"/Maps_exploitation/Correlog/Correlog_SP.richness_EDD_GLS_Spatial.jpeg"), quality =100)
par(mfcol=c(1,1))
plot(main = "Correlogram of residuals \n Species Richness ~ Faith's PD \ (GLS Spatial)", x = correlogram_m1$mean.of.class[1:20], y = correlogram_m1$correlation[1:20], type="b", col = c("black","red")[(correlogram_m1$p>0.05)+1], ylim = c(-1,1), pch=16, cex=1, lwd=1.5, xlab="Distance (km)", ylab="Moran's I", cex.lab=1.5, cex.axis=1.2) ; abline(h=0)
dev.off()

### Compute R and F-Stats

summary(m1)
p <- summary(m1)$tTable[2,4]
# F-value = SCM/ddlM / SCR/ddlR
F_value <- ((sum((predict(m1)-mean(sample.EDD))^2))/2)/((sum((sample.EDD-predict(m1))^2))/997)

# Compute regular (?) pseudo R² based on the traditional definition in OLS context => R² = 1 - SSR/SST = SSM/SST
pseudo.R2 <- 1 - (sum((sample.EDD-predict(m1))^2)/sum((sample.EDD-mean(sample.EDD))^2))
# Compute other R² applicable only in for a linear model
R2 <- cor(sample.EDD,predict(m1))^2

### Plot the model
jpeg(filename = paste0(internal.wd, "/Maps_exploitation/Correl_Plots/Correl_Plot_EDD_sp.richness_quad.RData"), quality = 100)
plot(sample.EDD ~ sample.sp.richness, main = "Correlation Plot \n Faith's PD ~ Species Richness \n (Quadratic model)", ylab = "Faith's Phylogenetic Distance", xlab = "Species Richness")
# abline(coef(m1), lwd = 2, col = "red") # For linear model
order.index <- order(sample.sp.richness)
lines(y = predict(m1)[order.index], x = sample.sp.richness[order.index], col = "red", lwd = 2)
legend(legend = c(paste0("Pseudo R² = ", round(pseudo.R2,2)), paste0("F = ", round(F_value, 0)), paste0("ddl = 2, 997") , paste0("p = 0.001")), x = "bottomright")
dev.off()

### Test to see if the predictions need the correlation structure need to be taken in account. Answer : NO !

load(file = paste0(internal.wd, "/Maps_exploitation/Ithomiini_range_raster.RData"))
?predict.gls

# Generate newdata df
newdata <-  data.frame(sp.richness[]) ; names(newdata) <- "sample.sp.richness"
newdata$sample.EDD <- PD.raster[]
newdata <- newdata[!is.na(Ithomiini.range.raster[]),] # Remove NA

# Extract projected coords of all non NA pixels
full_coords <- coordinates(Ithomiini.range.raster)[!is.na(Ithomiini.range.raster[]),]
full_coords_Spobj <- SpatialPoints(full_coords, proj4string = Ithomiini.range.raster@crs)
full_coords_Spobj <- spTransform(full_coords_Spobj, CRS("+proj=utm +zone=20 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")) # Exemple de projection : "+proj=longlat" pour le WGS84
proj_full_coords <- full_coords_Spobj@coords

?predict.gls
pred_net.EDD <- predict(object = m1, newdata = newdata, correlation = corExp(form = ~ full_coords[,1] + full_coords[,2]))
pred_net.EDD_2 <- predict(object = m1, newdata = newdata)
sum(pred_net.EDD_2 - pred_net.EDD)

pred.raster <- Ithomiini.range.raster
pred.raster[!is.na(pred.raster[])] <- pred_net.EDD

# Aucune différence ! La structure de correlation n'est prise en compte que pour l'estimation des paramètres, mais pas pour les prédictions par la suite...

### Direct predictions of net values from model coefs

pred.raster <- coef(m1)[1] + coef(m1)[2]*sp.richness + coef(m1)[3]*sp.richness^2
plot(pred.raster, col = pal_bl_red)
plot(PD.raster, col = pal_bl_red)
resid.raster <- pred.raster - PD.raster
resid.raster[is.na(Ithomiini.range.raster[])] <- NA
plot(resid.raster, col = pal_bl_red)

# Plot the contrasted version
load(file = paste0(internal.wd, "/crop_mask_shp.RData"))
resid.raster[resid.raster<(-50)] <- -50
resid.raster[resid.raster>50] <- 50
jpeg(filename = paste0(internal.wd, "/Maps_exploitation/Net_Faith_PD.jpeg"), quality = 100)
plot(resid.raster, col = pal_bl_red,  main = "Net Faith's Phylogenetic Diversity")
plot(crop_mask_shp, add = T)
dev.off()


### Models GLS linéaires (pas d'effets quadratiques testés)

# Faire les 45 combinaisons pour récupérer le r et plotter une Heatmap

load(file =  paste0(internal.wd, "/Maps_exploitation/proj_sample_coords.RData"))
load(file =  paste0(internal.wd, "/Maps_exploitation/sample.IN.OUT_index.RData"))

# Load maps of all 11 indices
load(file = paste0(internal.wd,"/Com.sp.richness.RData"))
load(file = paste0(internal.wd,"/Com.sp.diversity.compatible_Jost.RData"))
load(file = paste0(internal.wd,"/Com.sp.mean.rarity.RData"))

load(file = paste0(internal.wd,"/Com.MPD.RData"))
load(file = paste0(internal.wd,"/Com.PD.RData"))
load(file = paste0(internal.wd,"/Com.ED.RData"))

load(file = paste0(internal.wd,"/Com.ring.richness.RData"))
load(file = paste0(internal.wd,"/Com.ring.diversity.compatible_Jost.RData"))
load(file = paste0(internal.wd,"/Com.mimicry.mean.rarity.RData"))

load(file = paste0(internal.wd,"/Com.vulnerability.RData"))

load(file = paste0(internal.wd, "/Maps_exploitation/HII.raster.RData"))

# Generate df for object name and label name

Cor_indices_df <- as.data.frame(c("sp.richness", "sp.diversity.compatible_Jost", "mean.sp.rarity", "MPD.raster", "PD.raster", "ED.raster", "tot.ring.richness", "tot.ring.diversity.compatible_Jost", "mimicry.mean.rarity", "vulnerability", "HII.raster"))
names(Cor_indices_df) <- "Object"
Cor_indices_df$Main_label <- c("Species Richness", "Species Diversity", "Mean Species Rarity", "Mean Phylogenetic Pairwise Distance", "Faith's Phylogenetc Diversity", "Sum of Fair-Proportion", "Mimicry Richness", "Mimicry Diversity", "Mean Mimicry Rarity", "Vulnerability", "Human Influence Index")
Cor_indices_df$Y_label <- c("Species Richness", "Jost's Effective Species Richness", "Mean Species Rarity", "Mean Phylogenetic Pairwise Distance", "Faith's Phylogenetc Diversity", "Sum of Fair-Proportion", "Nb of Mimicry rings", "Jost's Effective Mimicry Richness", "Mean Mimicry Rarity", "Vulnerability", "HII (%)")

# Create a loop that test, compute F/R²/p and plot for all possible combinations

combinaisons <- t(combn(x = 1:11, m = 2, simplify = T))
R_coef_mat <- matrix(nrow = 11, ncol = 11, data = NA)


library(nlme)

for (i in 1:55) {
  # for (i in 1:nrow(combinaisons)) {
  comb <- combinaisons[i,]
  j <- comb[1] # Index A
  k <- comb[2] # Index B
  
  
  # Extract sample data from raster of index A
  sample.A <- eval(parse(text = paste0(Cor_indices_df$Object[j], "[sample.all_index]")))
  # Extract sample data from raster of index B
  sample.B <- eval(parse(text = paste0(Cor_indices_df$Object[k], "[sample.all_index]")))
  
  data.gls <- as.data.frame(cbind(sample.A, sample.B, proj_sample.all_coords))
  data.gls <- na.omit(data.gls)
  
  # Modeling
  # gls.model <- gls(data = data.gls, sample.A ~ sample.B, correlation = corExp(form = ~ x + y, nugget = TRUE), na.action = na.omit)
  # save(gls.model, file =  paste0(internal.wd, "/Maps_exploitation/Models/Models_", Cor_indices_df$Object[j], "_vs_" , Cor_indices_df$Object[k], ".RData"))
  load(file =  paste0(internal.wd, "/Maps_exploitation/Models/Models_", Cor_indices_df$Object[j], "_vs_" , Cor_indices_df$Object[k], ".RData"))
  
  # Modeling H0
  # gls.model.H0 <- gls(data = data.gls, sample.A ~ 1, correlation = corExp(form = ~ x + y, nugget = TRUE), na.action = na.omit)
  # save(gls.model.H0, file =  paste0(internal.wd, "/Maps_exploitation/Models/Models_H0_", Cor_indices_df$Object[j], "_vs_" , Cor_indices_df$Object[k], ".RData"))
  # load(file =  paste0(internal.wd, "/Maps_exploitation/Models/Models_H0_", Cor_indices_df$Object[j], "_vs_" , Cor_indices_df$Object[k], ".RData"))
  
  ### Compute p-value, r, R² and F-Stats
  
  # F-value = SCM/ddlM / SCR/ddlR
  F_value <- ((sum((predict(gls.model)-mean(data.gls$sample.A))^2))/2)/((sum((data.gls$sample.A-predict(gls.model))^2))/997)
  p <- df(x = F_value, df1 = 2, df2 = 997)
  
  # T-value de la pente et p-value associée
  T_value <- summary(gls.model)$tTable[2,3]
  p_T <- summary(gls.model)$tTable[2,4]
  
  dt(x = 0.6023832, df = 2)
  
  # Test du modèle par likelihood
  # Problème de la paramétrisation via RELM ? Voir cours STAR...
  # Delta.lik <- gls.model$logLik - gls.model.H0$logLik
  # p_Khi <- dchisq(x = Delta.lik, df = 1)
  
  # Compute other R², applicable only for a linear model
  R2 <- cor(data.gls$sample.A, predict(gls.model))^2
  # Get the sign of the relationship
  if (coef(gls.model)[2]<0) {
    R2 <- (-1)*R2
  }
  R_coef_mat[j,k] <- R2
  R_coef_mat[k,j] <- R2
  save(R_coef_mat, file = paste0(internal.wd, "/Maps_exploitation/R_coef_mat.RData"))
  
  ### Plot the model A~B
  jpeg(filename = paste0(internal.wd, "/Maps_exploitation/Correl_Plots/Correl_Plot_", Cor_indices_df$Object[k],"_", Cor_indices_df$Object[j],".jpeg"), quality = 100)
  original_ext_margins <- par()$oma
  par(oma = c(0,4,0,0))
  plot(data.gls$sample.A ~ data.gls$sample.B,
       main = paste0("Correlation Plot \n", Cor_indices_df$Main_label[j], " ~ ", Cor_indices_df$Main_label[k]),
       ylab = Cor_indices_df$Y_label[j], xlab = Cor_indices_df$Y_label[k],
       cex.axis = 1.7, cex.main = 1.5, cex.lab = 1.7, axis.args=list(cex.axis=1.7))
  abline(coef(gls.model), lwd = 2, col = "red") # For linear model
  legend(legend = c(paste0("R² = ", abs(round(R2,2))), paste0("F = ", round(F_value, 0)), paste0("p = ", round(p,3)), paste0("t = ", round(T_value, 2)), paste0("p = ", round(p_T, 3))), x = "bottomright")
  par(oma = original_ext_margins)
  dev.off()
  
  ### Plot the model B~A
  jpeg(filename = paste0(internal.wd, "/Maps_exploitation/Correl_Plots/Correl_Plot_", Cor_indices_df$Object[j],"_", Cor_indices_df$Object[k],".jpeg"), quality = 100)
  original_ext_margins <- par()$oma
  par(oma = c(0,4,0,0))
  plot(data.gls$sample.B ~ data.gls$sample.A, 
       main = paste0("Correlation Plot \n", Cor_indices_df$Main_label[k], " ~ ", Cor_indices_df$Main_label[j]), 
       ylab = Cor_indices_df$Y_label[k], xlab = Cor_indices_df$Y_label[j],
       cex.axis = 1.7, cex.main = 1.5, cex.lab = 1.7, axis.args=list(cex.axis=1.7))
  abline(a = (-coef(gls.model)[1])/coef(gls.model)[2], b = 1/coef(gls.model)[2], lwd = 2, col = "red") # For linear model
  legend(legend = c(paste0("R² = ", abs(round(R2,2))), paste0("F = ", round(F_value, 0)), paste0("p = ", round(p,3)), paste0("t = ", round(T_value, 2)), paste0("p = ", round(p_T, 3))), x = "bottomright")
  par(oma = original_ext_margins)
  dev.off()
  
  cat(paste0(Sys.time(), " - Model n°", i, "\n"))
}

# R_coef_mat[9,10] <- R_coef_mat[10,9] <- 0.05

rownames(R_coef_mat) <- c("Species Richness", "Species Diversity", "Mean Species Rarity", "MPD", "Faith's PD", "Sum of Fair-Proportion", "Mimicry Richness", "Mimicry Diversity", "Mean Mimicry Rarity", "Vulnerability", "Human Influence Index")
colnames(R_coef_mat) <- c("Species Richness", "Species Diversity", "Mean Species Rarity", "MPD", "Faith's PD", "Sum of Fair-Proportion", "Mimicry Richness", "Mimicry Diversity", "Mean Mimicry Rarity", "Vulnerability", "Human Influence Index")

diag(R_coef_mat) <- 1

save(R_coef_mat, file = paste0(internal.wd, "/Maps_exploitation/R_coef_mat.RData"))
load(file = paste0(internal.wd, "/Maps_exploitation/R_coef_mat.RData"))

# R² heat map
library(corrgram)
jpeg(filename = paste0(internal.wd, "/Maps_exploitation/Correl_Plots/Heat_Map_Rsq.RData"), quality = 100)
corrgram(x = R_coef_mat, type = "cor", order = F, lower.panel = panel.shade,  upper.panel = panel.cor, diag.panel = NULL, col.regions = colorRampPalette(c("navy", "royalblue", "white", "orange", "red"))) # pour plotter les correlations entre variables de manière stylée
dev.off()

# Pearson's r heatmap
neg.index <- R_coef_mat<0
Pearson_R_coef_mat <- sqrt(abs(R_coef_mat))
Pearson_R_coef_mat[neg.index] <- -Pearson_R_coef_mat[neg.index]
Pearson_R_coef_mat

# Pearson's r heat map
library(corrgram)
jpeg(filename = paste0(internal.wd, "/Maps_exploitation/Correl_Plots/Heat_Map_R_Pearson.RData"), quality = 100)
corrgram(x = Pearson_R_coef_mat, type = "cor", order = F, lower.panel = panel.shade,  upper.panel = panel.cor, diag.panel = NULL, col.regions = colorRampPalette(c("navy", "royalblue", "white", "orange", "red"))) # pour plotter les correlations entre variables de manière stylée
dev.off()

# Dendrogram des Indices

# To compute on real R²
# R_coef_mat <- abs(R_coef_mat)

jpeg(filename = paste0(internal.wd, "/Maps_exploitation/Correl_Plots/Heat_Map_R_Sq_Dendro.jpeg"), quality = 100)
original_ext_margins <- par()$oma
par(oma = c(7,0,0,7))
heatmap(x = R_coef_mat, hclustfun = hclust, symm = T, col = pal_bl_red)
par(oma = original_ext_margins)
dev.off()

?hclust # Method par défaut = "complete linkage"

# Scale for the heatmap
Ithomiini.range.raster@data@values[1] <- 0
tiff(filename = paste0(internal.wd, "/Maps_exploitation/Correl_Plots/Heat_Map_Scale.tiff"))
plot(Ithomiini.range.raster, col = pal_bl_red,
     axis.args=list(cex.axis=1.7))
dev.off()


##### GAM approach (Zupan et al., 2014)

# Utilise les coordonnées géographiques comme prédicteurs inclus avec un smooth
?gam

x_cor <- coordinates(PD.raster)[,1]
y_cor <- coordinates(PD.raster)[,2]

library(mgcv)
model <- mgcv::gam(PD.raster[] ~ sp.richness[] + s(x_cor, y_cor)) # s() pour smooth la variable explicative
summary(model) ; plot(model) # Permet de visualiser la forme de la relation sous forme de courbe lissée et de détecter d'éventuel effet quadratique à tester en modèle LM/GLM/MM/GMM classiques via anova(mod1,mod2)


##### Test correlation avec functions package raster => ne prend pas en compte la strucutre spatiale !!!
load(file = paste0(internal.wd,"/Com.sp.richness.RData"))

load(file = paste0(internal.wd,"/Com.PD.RData"))
cor.test(x = sp.richness[], y = PD.raster[], method = "spearman", exact = F)
# Ne prend pas en compte la strucutre spatiale...

?corLocal # Pour faire des tests au sein d'une moving window afin d'obtenir un raster qui montre la correlation locale
# ngb = size of moving windows (odd number)
local.cor <- corLocal(x = sp.richness, y = PD.raster, method = "spearman", ngb = 5, test = T)
plot(local.cor)
local.cor.p <- local.cor[[2]]
local.cor.p[local.cor.p>0.05] <- 1
local.cor.p[local.cor.p<=0.05] <- 0
plot(local.cor.p)