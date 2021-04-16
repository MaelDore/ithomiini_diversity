
####################### Test and plot correlation ring size ~ ring richness ############################

### Create summary table for mimicry rings

### Compute correlation

### Plot correlation


### 1/ Load stuff

# Effacer l'environnement
rm(list = ls())

# Load summary table for mimicry rings
load(file = paste0("./outputs/N_OMUs_per_ring.RData"))
list.rings <- N_OMUs_per_ring

# Load probability stack for mimicry rings
ring_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_Jaccard.80.RData"))

### 2/ Compute range size and geographic rarity

list.rings$range_size <- NA
list.rings$geographic_rarity <- NA
for (i in 1:nlayers(ring_proba_stack)) {
  # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
  # Compute the estimated number of pixels occupied by the ring as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
  list.rings$range_size[i] <- sum(ring_proba_stack[[i]]@data@values, na.rm = T)*774.51/1000
}
list.rings$geographic_rarity <- 1 - (list.rings$range_size/max(list.rings$range_size)) # Weights by max ring extent = rarity indices

# Compute range size standardized by number of species
list.rings$range_size_std <- list.rings$range_size / list.rings$counts

### 3/ Compute mean range size for OMU/species among each ring

load(file = paste0("./input_data/list.models.RData"))
OMU_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/All_OMU_stack_Jaccard.80.rds"))

identical(list.models$Tag.model, names(OMU_proba_stack))

for (i in 1:nlayers(OMU_proba_stack)) {
  # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
  # Compute the estimated number of pixels occupied by the ring as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
  list.models$OMU_range_size[i] <- sum(OMU_proba_stack[[i]]@data@values, na.rm = T)*774.51/1000
}
list.models$geographic_rarity <- 1 - (list.models$OMU_range_size/max(list.models$OMU_range_size)) # Weights by max ring extent = rarity indices

# save(list.models, file = paste0("./input_data/list.models.RData"))

list.rings <- list.models %>% 
  group_by(Mimicry.model) %>% 
  summarize(mean_range_size_per_sp = mean(OMU_range_size)) %>% 
  ungroup() %>% 
  left_join(x = list.rings, y = ., by = c("ring" = "Mimicry.model"))

# save(list.rings, file = paste0("./outputs/list.rings.RData"))

### 3/ Test for correlation between size and range

# Not normally distributed
hist(list.rings$range_size)
hist(list.rings$range_size_std)
hist(list.rings$mean_range_size_per_sp)
hist(list.rings$counts)

# Raw range
cor.test(x = list.rings$range_size, y = list.rings$counts, method = "spearman", alternative = "two.sided")
# rho = 0.758, p-value << 0.001

# Range std per richness
cor.test(x = list.rings$range_size_std, y = list.rings$counts, method = "spearman", alternative = "two.sided")
# rho = -0.071, p-value = 0.64

# Mean range per sp
cor.test(x = list.rings$mean_range_size_per_sp, y = list.rings$counts, method = "spearman", alternative = "two.sided")
# rho = 0.316, p-value = 0.036


### 4/ Plot correlation

### Correlation plot function

{
correlation_plot <- function(data_df,
                             x_var,
                             y_var,
                             
                             main = "",
                             x_label,
                             y_label,
                             
                             x_transfo = "identity",
                             y_transfo = "identity",
                             
                             facet_label = "",
                             facet_label_size = 7,
                             facet_label_hjust = 0.5,
                             facet_label_vjust = 0,

                             rho_label = "",
                             rho_size = 5,
                             rho_hjust = 0.94,
                             rho_vjust = -3.3,
                             
                             p_label = "",
                             p_size = 5,
                             p_hjust = 0.94,
                             p_vjust = -2.5
                             
                             )
{
  # x_var_full <- eval(parse(text = paste0(data_df,"$",x_var)))
  # y_var_full <- eval(parse(text = paste0(data_df,"$",x_var)))
  
  ggplot(data = data_df, aes(x = x_var, y = y_var)) +
    geom_point() +
    geom_smooth(color = "red", size = 1.5, se = F) +
    annotate("text", x = max(x_var), y = min(y_var), hjust = facet_label_hjust, vjust = facet_label_vjust, label = facet_label, size = facet_label_size, fontface = 2) +
    annotate("text", x = max(x_var), y = min(y_var), hjust = rho_hjust, vjust = rho_vjust, label = rho_label, size = rho_size, fontface = 2) +
    annotate("text", x = max(x_var), y = min(y_var), hjust = p_hjust, vjust = p_vjust, label = p_label, size = p_size, fontface = 2) +
    # annotate("rect", xmin = 135, xmax = 145, ymin = -0.02, ymax = 0.05, fill = "white", col = NA) +
    # annotate("text", label = "A", x = 140, y = 0.02, size = 7, colour = "black", fontface = 2) +
    ggtitle(label = main) +
    xlab(label = x_label) +
    ylab(label = y_label) +
    scale_x_continuous(trans = x_transfo) + 
    scale_y_continuous(trans = y_transfo) + 
    # scale_y_continuous(breaks = seq(0, 1, 0.25), limits = c(0, 1)) + 
    theme_gray() +
    theme(panel.grid.minor.x=element_blank(),
          panel.border = element_rect(fill = NA),
          axis.ticks = element_line(size = 1.2, color = "black"),
          axis.ticks.length = unit(5, "pt"),
          axis.text = element_text(size = 12, face = "bold", color = "black"),
          axis.title = element_text(size = 14, face = "bold"),
          axis.title.x = element_text(margin = margin(t = 10, b = 5)),
          axis.title.y = element_text(margin = margin(l = 5, r = 10)))
}
}

### 4.1/ Raw range ~ size

corr_ring_size_range <- correlation_plot(
  data_df = list.rings,
  x_var = list.rings$counts, 
  y_var = list.rings$range_size,
  main = "Correlation range ~ size: ring range",
  x_label = "Ring size  [species]", 
  y_label = "Ring range (log scale)\n[10^3 km²]",
  y_transfo = "log10", 
  facet_label = "a",
  rho_label = bquote(rho ~ "= 0.758"),
  p_label = "p-value < 0.001")
                
# save(corr_ring_size_range, file = paste0("./outputs/Correlation_tests/corr_ring_size_range.RData"))
# load(file = paste0("./outputs/Correlation_tests/corr_ring_size_range.RData")

pdf(file = "./graphs/Correlation_tests/corr_ring_size_range.pdf", height = 6, width = 8)
print(corr_ring_size_range)
dev.off()

### 4.2/ Std range ~ size

corr_ring_size_std_range <- correlation_plot(
  data_df = list.rings,
  x_var = list.rings$counts, 
  y_var = list.rings$range_size_std,
  main = "Correlation range ~ size: ring standardized range",
  x_label = "Ring size  [species]", 
  y_label = "Ring range std per species\n(log scale) [10^3 km²]",
  y_transfo = "log10", 
  facet_label = "b",
  rho_label = bquote(rho ~ "= -0.071"),
  p_label = "p-value = 0.64")

# save(corr_ring_size_std_range, file = paste0("./outputs/Correlation_tests/corr_ring_size_std_range.RData"))
# load(file = paste0("./outputs/Correlation_tests/corr_ring_size_std_range.RData")

pdf(file = "./graphs/Correlation_tests/corr_ring_size_std_range.pdf", height = 6, width = 8)
print(corr_ring_size_std_range)
dev.off()

### 4.3/ Mean sp range ~ size

corr_ring_size_mean_sp_range <- correlation_plot(
  data_df = list.rings,
  x_var = list.rings$counts, 
  y_var = list.rings$mean_range_size_per_sp,
  main = "Correlation range ~ size: mean sp range per ring",
  x_label = "Ring size  [species]", 
  y_label = "Mean sp range (log scale)\n[10^3 km²]",
  y_transfo = "log10", 
  facet_label = "c",
  rho_label = bquote(rho ~ "= 0.316"),
  p_label = "p-value = 0.036")

# save(corr_ring_size_mean_sp_range, file = paste0("./outputs/Correlation_tests/corr_ring_size_mean_sp_range.RData"))
# load(file = paste0("./outputs/Correlation_tests/corr_ring_size_mean_sp_range.RData")

pdf(file = "./graphs/Correlation_tests/corr_ring_size_mean_sp_range.pdf", height = 6, width = 8)
print(corr_ring_size_mean_sp_range)
dev.off()

library(gridExtra)

pdf(file = "./graphs/Correlation_tests/corr_ring_size_all.pdf", height = 12, width = 14)

grid.arrange(corr_ring_size_range, corr_ring_size_std_range, corr_ring_size_mean_sp_range,
             widths = c(10, 0.4, 10),  # Width of columns
             heights = c(10, 0.2, 10),
             nrow = 3,
             ncol = 3,
             layout_matrix = rbind(c(1, NA, 2),  # Position of ggplots in the layout
                                   c(NA, NA, NA),
                                   c(3, NA, 4)))

dev.off()
