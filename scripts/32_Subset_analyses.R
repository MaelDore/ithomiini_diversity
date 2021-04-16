###### Script 32: Is Upper Amazon species pool a subset of the Andean one?  ######

# Compute overlap between Upper Amazon species pool and Andean species pool

###

# Inputs: 
# Species probability stack
# OMU probability stack
# Bioregions borders

# Outputs:
# Venn Diagram of species pools

###


# Clean environment
rm(list = ls())
gc()

### 1/ Load stuff ####

library(raster)
library(sf)

# Load stack of species probability of presence
sp_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))
OMU_proba_stack <- readRDS(file = paste0("./outputs/Indices_stacks/All_OMU_stack_Jaccard.80.rds"))

# Load regions shp files
Northern_Andes_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Northern_Andes_raw_shp.rds")
Central_Andes_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Central_Andes_raw_shp.rds")
Western_Lowlands_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Western_Lowlands_raw_shp.rds")
Western_Amazon_shp <- readRDS(file = "./input_data/Map_stuff/Bioregions/Western_Amazon_raw_shp.rds")

plot(Western_Amazon_shp)
plot(Northern_Andes_shp)
plot(Central_Andes_shp)
plot(Western_Lowlands_shp)

### 2/ Create mask for each bioregions ####

Western_Amazon_mask <- rasterize(x = Western_Amazon_shp, y = sp_proba_stack[[1]], # Provide the grid to fill with CRS, bbox and resolution
                                 field = 1, background = NA, # Value to use to fill empty cells
                                 fun = max)
Northern_Andes_mask <- rasterize(x = Northern_Andes_shp, y = sp_proba_stack[[1]], # Provide the grid to fill with CRS, bbox and resolution
                                 field = 1, background = NA, # Value to use to fill empty cells
                                 fun = max)
Central_Andes_mask <- rasterize(x = Central_Andes_shp, y = sp_proba_stack[[1]], # Provide the grid to fill with CRS, bbox and resolution
                                field = 1, background = NA, # Value to use to fill empty cells
                                fun = max)
Western_Lowlands_mask <- rasterize(x = Western_Lowlands_shp, y = sp_proba_stack[[1]], # Provide the grid to fill with CRS, bbox and resolution
                                field = 1, background = NA, # Value to use to fill empty cells
                                fun = max)

save(Western_Amazon_mask, Northern_Andes_mask, Central_Andes_mask, Western_Lowlands_mask, file = paste0("./outputs/Subset_analyses/all_regions_masks.RData"))

plot(Western_Amazon_mask)
plot(Northern_Andes_mask)
plot(Central_Andes_mask)
plot(Western_Lowlands_mask)

### 3/ Extract species/OMU list for each bioregions ####

### 3.1/ Function to extract species P/A or list from raster brick and generate associated objects

extract_species_list <- function(stack, mask, type)
{
  # type <- c("sp", "OMU")
  
  mask_name <- deparse(substitute(mask)) # Get the name of the mask region as a character string
  
  brick <- readAll(brick(stack)) # Extract data of probability of presences in a matrix format
  
  com_mat <- brick@data@values[!is.na(mask[]), ] # Extract only community within the bioregion/mask
  com_bin <- apply(X = com_mat, MARGIN = 2, FUN = function(x) {x >= 0.5}) # Apply filter to consider only pixels with p >= 0.5
  unit_bin <- apply(X = com_bin, MARGIN = 2, FUN = any) # Mark as presence only unit with at least one pixel with p >= 0.5
  unit_bin[is.na(unit_bin)] <- F # Replace NA by FALSE
  unit_list <- names(stack)[unit_bin] # Extract unit names
  
  # Generate new object with "type_bin" suffix in the global environment
  eval(call("<<-", as.name(paste0(mask_name,"_", type, "_bin")), unit_bin))
  # Generate new object with "type_bin" suffix in the global environment
  eval(call("<<-", as.name(paste0(mask_name,"_", type, "_list")), unit_list))
  
  # Save unit bin and list
  saveRDS(unit_bin, file = paste0("./outputs/Subset_analyses/",mask_name, "_", type, "_bin.rds"))
  saveRDS(unit_list, file = paste0("./outputs/Subset_analyses/",mask_name, "_", type, "_list.rds"))
  
}

### 3.2/ For Species ####

# A species belong to a region if it as at least one cell with (p > 0.5)

extract_species_list(stack = sp_proba_stack, mask = Western_Amazon_mask, type = "sp")
extract_species_list(stack = sp_proba_stack, mask = Northern_Andes_mask, type = "sp")
extract_species_list(stack = sp_proba_stack, mask = Central_Andes_mask, type = "sp")
extract_species_list(stack = sp_proba_stack, mask = Western_Lowlands_mask, type = "sp")

### 3.3/ For OMU ####

# An OMU belongs to a region if it as at least one cell with (p > 0.5) 

extract_species_list(stack = OMU_proba_stack, mask = Western_Amazon_mask, type = "OMU")
extract_species_list(stack = OMU_proba_stack, mask = Northern_Andes_mask, type = "OMU")
extract_species_list(stack = OMU_proba_stack, mask = Central_Andes_mask, type = "OMU")
extract_species_list(stack = OMU_proba_stack, mask = Western_Lowlands_mask, type = "OMU")

### 4/ Get core-group and backbone differentiation ####

### 4.1/ Load phylogeny ####

library(ape)
library(geiger)

phylo.Ithomiini <- readRDS(file = "./input_data/Phylogenies/Final_phylogeny.rds")

pdf(file = paste0("./graphs/quick_phylo.pdf"), height = 50, width = 12)
plot(phylo.Ithomiini)
nodelabels()
dev.off()

# Core-group = Nod 390

### 4.2/ Extract names of species in the phylogeny belonging to core-group and backbone ####

core_group_phylo_sp <- tips(phylo.Ithomiini, node = 390) # Retrieve the set of species descending from nod 390
backbone_group_phylo_sp <- setdiff(phylo.Ithomiini$tip.label, core_group_phylo_sp)

core_group_genera <- unique(str_split(string = core_group_phylo_sp, pattern = "\\.", simplify = T)[, 1])
backbone_group_genera <- unique(str_split(string = backbone_group_phylo_sp, pattern = "\\.", simplify = T)[, 1])

### 4.3/ Update list.models and list.sp ####

load(file = "./input_data/list.models.RData")

list.models$Phylo_group[list.models$Genus %in% core_group_genera] <- "core-group"
list.models$Phylo_group[list.models$Genus %in% backbone_group_genera] <- "backbone"

table(list.models$Phylo_group)
# 670 Core-group 
# 113 Backbone

# save(list.models, file = "./input_data/list.models.RData")

load(file = "./input_data/list.sp.RData")

list.sp$Phylo_group[list.sp$Genus %in% core_group_genera] <- "core-group"
list.sp$Phylo_group[list.sp$Genus %in% backbone_group_genera] <- "backbone"

table(list.sp$Phylo_group)
# 330 Core-group ???
# 58 Backbone ???

# save(list.sp, file = "./input_data/list.sp.RData")

### 4.4/ Extract full list of OMU and sp per group

core_group_sp <- list.sp$Sp_full[list.sp$Phylo_group == "core-group"]
backbone_group_sp <- list.sp$Sp_full[list.sp$Phylo_group == "backbone"]

### 5/ Extract list of species per phylo-group per bioregions ####

Western_Amazon_mask_sp_list_coregroup <- intersect(Western_Amazon_mask_sp_list, core_group_sp)
Western_Amazon_mask_sp_list_backbone <- intersect(Western_Amazon_mask_sp_list, backbone_group_sp)

Northern_Andes_mask_sp_list_coregroup <- intersect(Northern_Andes_mask_sp_list, core_group_sp)
Northern_Andes_mask_sp_list_backbone <- intersect(Northern_Andes_mask_sp_list, backbone_group_sp)

Central_Andes_mask_sp_list_coregroup <- intersect(Central_Andes_mask_sp_list, core_group_sp)
Central_Andes_mask_sp_list_backbone <- intersect(Central_Andes_mask_sp_list, backbone_group_sp)

### 6/ Plot Venn diagrams ####

library(VennDiagram)
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)

?venn.diagram

# Make the plot

{
venn_diag_4 <- venn.diagram(
  x = list(Western_Amazon_mask_sp_list, Northern_Andes_mask_sp_list, Central_Andes_mask_sp_list, Western_Lowlands_mask_sp_list),
  category.names = c("Western Amazon" , "Northern Andes" , "Central Andes", "Western Lowlands"),
  filename = "./graphs/Venn_diag/All_4_regions.tiff",
  imagetype = "tiff" ,
  output = T,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1, # ector giving the width of each circle's circumference
  col = c("#440154ff", '#21908dff', '#fde725ff', "red"),
  fill = c("#440154ff", '#21908dff', '#fde725ff', "red"),
  alpha = 0.3,
  cex = 0.5, # Vector giving the size for each area label
  fontfamily = "sans",
  cat.cex = 0.3, # Vector giving the size for each category name
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135, -135), # 	Vector giving the position (in degrees) of each category name along the circle, with 0 at 12 o'clock
  cat.dist = c(0.055, 0.055, 0.085, 0.085), # Vector giving the distance (in npc units) of each category name from the edge of the circle (can be negative)
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#fde725ff', "red"), # Vector giving the colour for each category name
  rotation.degree = 0,
  margin = 0.2
)

grid.draw(venn.plot)


venn_diag_3 <- venn.diagram(
  x = list(Western_Amazon_mask_sp_list, Northern_Andes_mask_sp_list, Central_Andes_mask_sp_list),
  category.names = c("Western Amazon" , "Northern Andes" , "Central Andes"),
  filename = "./graphs/Venn_diag/Andes_2_Upper_Amazon_regions.tiff",
  imagetype = "tiff" ,
  output = T,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1, # ector giving the width of each circle's circumference
  col = c("#440154ff", '#21908dff', '#fde725ff'),
  fill = c("#440154ff", '#21908dff', '#fde725ff'),
  alpha = 0.3,
  cex = 0.5, # Vector giving the size for each area label
  fontfamily = "sans",
  cat.cex = 0.3, # Vector giving the size for each category name
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135), # 	Vector giving the position (in degrees) of each category name along the circle, with 0 at 12 o'clock
  cat.dist = c(0.055, 0.055, 0.085), # Vector giving the distance (in npc units) of each category name from the edge of the circle (can be negative)
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff', '#fde725ff'), # Vector giving the colour for each category name
  rotation.degree = 0,
  margin = 0.2
)

# Merge the 2 Andean list of species
Andes_mask_sp_list <- union(Northern_Andes_mask_sp_list, Central_Andes_mask_sp_list)

venn_diag_2 <- venn.diagram(
  x = list(Western_Amazon_mask_sp_list, Andes_mask_sp_list),
  category.names = c("Western Amazon" , "Central & Northern Andes"),
  filename = "./graphs/Venn_diag/Andes_Upper_Amazon_regions.tiff",
  imagetype = "tiff" ,
  output = T,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1, # Vector giving the width of each circle's circumference
  col = c("#440154ff", '#21908dff'),
  fill = c("#440154ff", '#21908dff'),
  alpha = 0.3,
  cex = 0.5, # Vector giving the size for each area label
  fontfamily = "sans",
  cat.cex = 0.3, # Vector giving the size for each category name
  cat.default.pos = "outer",
  cat.pos = c(-27, 27), # 	Vector giving the position (in degrees) of each category name along the circle, with 0 at 12 o'clock
  cat.dist = c(0.055, 0.055), # Vector giving the distance (in npc units) of each category name from the edge of the circle (can be negative)
  cat.fontfamily = "sans",
  cat.col = c("#440154ff", '#21908dff'), # Vector giving the colour for each category name
  rotation.degree = 0,
  margin = 0.2
)
}

# Manual drawing

### 6.1/ Complete lists ####

# Merge the 2 Andean list of species
Andes_mask_sp_list <- union(Northern_Andes_mask_sp_list, Central_Andes_mask_sp_list)
Andes_mask_sp_list_coregroup <- union(Northern_Andes_mask_sp_list_coregroup, Central_Andes_mask_sp_list_coregroup)
Andes_mask_sp_list_backbone <- union(Northern_Andes_mask_sp_list_backbone, Central_Andes_mask_sp_list_backbone)

save(Andes_mask_sp_list, Andes_mask_sp_list_coregroup, Andes_mask_sp_list_backbone, file = paste0("./outputs/Subset_analyses/Andes_mask_lists.RData"))

pdf(file = paste0("./graphs/Venn_diag/2_regions.pdf"), height = 12, width = 12)

grid.newpage()
venn_plot_2 <- draw.pairwise.venn(
  area1 = length(Western_Amazon_mask_sp_list),
  area2 = length(Andes_mask_sp_list),
  cross.area = length(intersect(Western_Amazon_mask_sp_list, Andes_mask_sp_list)),
  category = c("Uppern\nAmazon", "Central &\nNorthern\nAndes"),
  fontfamily = "sans",
  col = c("blue", "tan4"),
  fill = c("dodgerblue", "bisque3"),
  alpha = 0.3,
  lty = 1,
  cex = 2.6,
  cat.cex = 2.8,
  cat.pos = c(-55, 60),
  cat.col = c("blue", "tan4"),
  cat.dist = c(0.08, 0.08),
  cat.fontfamily = "sans",
  ext.pos = 30,
  ext.dist = -0.07,
  ext.length = 0.60,
  ext.line.lwd = 2,
  ext.line.lty = 1,
  label.col = c("blue", "black", "tan4"),
  print.mode = c("raw", "percent"),
  sigdigs = 2,
  margin = 0.07
)
grid.draw(venn_plot_2)
grid.text(label = expression(bold("Species list overlaps across regions")), x = unit(0.5, "npc"), y = unit(0.90, "npc"), gp = gpar(fontsize = 35))
grid.text(label = expression(bold("(d)")), x = unit(0.9, "npc"), y = unit(0.18, "npc"), gp = gpar(fontsize = 35))
grid.newpage()

dev.off()


### 6.2/ Coregroup ####

pdf(file = paste0("./graphs/Venn_diag/2_regions_coregroup.pdf"), height = 12, width = 12)

grid.newpage()
venn_plot_2_coregroup <- draw.pairwise.venn(
  area1 = length(Western_Amazon_mask_sp_list_coregroup),
  area2 = length(Andes_mask_sp_list_coregroup),
  cross.area = length(intersect(Western_Amazon_mask_sp_list_coregroup, Andes_mask_sp_list_coregroup)),
  category = c("Uppern\nAmazon", "Central\n& Northern\nAndes"),
  fontfamily = "sans",
  col = c("blue", "tan4"),
  fill = c("dodgerblue", "bisque3"),
  alpha = 0.3,
  lty = 1,
  cex = 1.8,
  cat.cex = 2.2,
  cat.pos = c(-55, 55),
  cat.col = c("blue", "tan4"),
  cat.dist = 0.10,
  cat.fontfamily = "sans",
  ext.pos = 30,
  ext.dist = -0.07,
  ext.length = 0.60,
  ext.line.lwd = 2,
  ext.line.lty = 1,
  label.col = c("blue", "black", "tan4"),
  print.mode = c("raw", "percent"),
  sigdigs = 2,
  margin = 0.1
)
grid.draw(venn_plot_2_coregroup)
grid.text(label = expression(bold("Coregroup")), x = unit(0.5, "npc"), y = unit(0.85, "npc"), gp = gpar(fontsize = 20))
grid.newpage()

dev.off()


### 6.3/ Backbone ####

pdf(file = paste0("./graphs/Venn_diag/2_regions_backbone.pdf"), height = 12, width = 12)

grid.newpage()
venn_plot_2_backbone <- draw.pairwise.venn(
  area1 = length(Western_Amazon_mask_sp_list_backbone),
  area2 = length(Andes_mask_sp_list_backbone),
  cross.area = length(intersect(Western_Amazon_mask_sp_list_backbone, Andes_mask_sp_list_backbone)),
  category = c("Uppern\nAmazon", "Central &\nNorthern Andes"),
  fontfamily = "sans",
  col = c("blue", "tan4"),
  fill = c("dodgerblue", "bisque3"),
  alpha = 0.3,
  lty = 1,
  cex = 1.8,
  cat.cex = 2.2,
  cat.pos = c(-55, 55),
  cat.col = c("blue", "tan4"),
  cat.dist = 0.10,
  cat.fontfamily = "sans",
  ext.pos = 30,
  ext.dist = -0.07,
  ext.length = 0.60,
  ext.line.lwd = 2,
  ext.line.lty = 1,
  label.col = c("blue", "black", "tan4"),
  print.mode = c("raw", "percent"),
  sigdigs = 2,
  margin = 0.1
)
grid.draw(venn_plot_2_backbone)
grid.text(label = expression(bold("Backbone")), x = unit(0.5, "npc"), y = unit(0.85, "npc"), gp = gpar(fontsize = 20))
grid.newpage()

dev.off()

### 6.4/ Plot all 3 diagram together ####

# With cowplot::plot_grid()

?plot_grid

pdf(file = paste0("./graphs/Venn_diag/all_diags.pdf"), height = 12, width = 12)

cowplot::plot_grid(venn_plot_2, NULL, venn_plot_2_coregroup, venn_plot_2_backbone,
                   ncol = 2, nrow = 2,
                   labels = c("(a)", "(b)", "(c)", "(d)"),
                   label_x = 0.95, label_y = 0.05,
                   label_size = 18)
dev.off()


# With grid::grid.arrange()

# library("grid")
# 
# ?grid.arrange
# 
# grid.arrange(venn_plot_2)
# 
# grid.arrange(grobs = list(venn_plot_2, venn_plot_2, venn_plot_2_coregroup, venn_plot_2_backbone), nrow = 2, ncol = 2, widths = c(3, 2))
# 
# facet_plot <- arrangeGrob(grobs = list(venn_plot_2, venn_plot_2_coregroup, venn_plot_2_backbone, venn_plot_2), nrow = 2, ncol = 2, widths = c(3, 2))
# grid.arrange(facet_plot)
# grid.draw(facet_plot)


### 6.5/ Plot full list in facet (d) position

pdf(file = paste0("./graphs/Venn_diag/full_diags_facet_d.pdf"), height = 12, width = 12)

cowplot::plot_grid(NULL, NULL, NULL, venn_plot_2,
                   ncol = 2, nrow = 2)
grid.text(label = expression(bold("Species list overlaps across regions")), x = unit(0.75, "npc"), y = unit(0.45, "npc"), gp = gpar(fontsize = 15))
grid.text(label = expression(bold("(d)")), x = unit(0.95, "npc"), y = unit(0.05, "npc"), gp = gpar(fontsize = 18))

dev.off()