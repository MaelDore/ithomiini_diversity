
##### Script 24: Create final distribution maps with phylogeny as inset #####

# Create final distribution maps with phylogeny as inset, for the online archive (10.5281/zenodo.4673446)

# 4 types of maps:
   # 783 OMU distribution maps
   # 388 species distribution maps
   # 44 mimicry ring distributions maps
   # 44 mimicry ring richness maps

###

# Input: 
  
  # Ithominii phylogeny (Chazot et al., 2019)
  # OMU/species/mimicry ring maps from scripts 11/12/13

# Outputs:
  
  # Pretty phylogeny with subtribes (Figure S3.4)
  # Summary table for subtribes
  # Custom phylogeny to use as inset in maps
  # Distribution map with embedded phylogeny for each OMU/species/mimicry ring (10.5281/zenodo.4673446)
  # Richness map with embedded phylogeny for each mimicry ring (10.5281/zenodo.4673446) 

###

### 1/ Load stuff ####

# Clean environnement
rm(list = ls())

# Load Packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("treeio")
BiocManager::install("ggtree")

library(tidyverse)
library(dplyr)

library(treeio)    # To read tree from any phylogenetic format
library(ggtree)    # To plot tree using the tidyverse grammar
library(tidytree)  # To manipulate tlb_tree tibble that store info on tree topology as list of edges

library(ape)
library(picante)
library(geiger)
library(treespace)

# phylo.Ithomiini <- readRDS(file = "./input_data/Phylogenies/Final_phylogeny.rds")
phylo.Ithomiini <- readRDS(file = "./input_data/Phylogenies/Final_phylogeny_with_updated_names.rds")
load(file = paste0("./input_data/list.sp.RData"))

# ### Clean sp name
# list.sp$Abs_phlyo <- F
# list.sp$Abs_phlyo[str_which(string = list.sp$Species, pattern = "AbsPhylo")] <- T
# 
# list.sp$Species <- str_replace(string = list.sp$Species, pattern = "AbsPhylo", replacement = "")
# list.sp$Sp_full <- str_replace(string = list.sp$Sp_full, pattern = "AbsPhylo", replacement = "")
# 
# save(list.sp, file = paste0("./input_data/list.sp.RData"))

# ### Add Subtribes to list.sp
# Ithomiini_taxonomy_with_subtribes <- xlsx::read.xlsx("./input_data/Ithomiini_taxonomy_with_subtribes.xlsx", sheetIndex = 1)
# Ithomiini_taxonomy_with_subtribes$sp_full <- paste(Ithomiini_taxonomy_with_subtribes$Genus, Ithomiini_taxonomy_with_subtribes$Species, sep = ".")
# 
# list.sp$Subtribe[match(Ithomiini_taxonomy_with_subtribes$sp_full, list.sp$sp_full_new)] <- as.character(Ithomiini_taxonomy_with_subtribes$Subtribe)
# 
# save(list.sp, file = paste0("./input_data/list.sp.RData"))


# # Prepare df of species (tips) data
# list_sp_for_tree <- list.sp %>% 
#   slice(match(phylo.Ithomiini$tip.label, list.sp$Sp_full)) %>%  # Remove from list species that are not in the phylogeny. And reorder df to match phylogeny tips
#   rename("label" = "Sp_full")  # Rename Sp_full as 'label' to include them into the datatree obj
# 
# # Transform into a treedata object (treeio package)
# phylo.Ithomiini_treedata <- phylo.Ithomiini %>% 
#   as.treedata() %>% 
#   full_join(x = ., y = list_sp_for_tree, by = "label") # Associate the list_sp table data to each tip
# 
# save(phylo.Ithomiini_treedata, file = "./phylogenies/phylo.Ithomiini_treedata.RData")

load(file = "./phylogenies/phylo.Ithomiini_treedata.RData")

### 2/ Plot first phylogeny with label and nod and branch indices ####
pdf(file = "./phylogenies/Full_phylo_rectangular_all_labels.pdf", height = 50, width = 20)
Ithomiini_phylo_plot <- ggtree(phylo.Ithomiini_treedata,
                               layout = "rectangular") +
  geom_tiplab(align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
               linetype = 0, 
               size = 2, 
               offset = 0.2,    # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
               hjust = 0) +    # To center, left-align, right-aligned text
  geom_label(aes(label = node),
             fill = 'steelblue', size = 1) +
  geom_label(aes(x = branch, label = node),
             fill = 'lightgreen', size = 1) 

print(Ithomiini_phylo_plot)
dev.off()


### 3/ Extract genera and subtribe info ####

## 3.1/ Find the MRCA nod for each genera ####

Genus <- as.character(unique(phylo.Ithomiini_treedata@extraInfo$Genus))
Genus <- Genus[!is.na(Genus)] # Remove NA
list_genera <- as.data.frame(Genus)

for (i in seq_along(list_genera$Genus)) 
{
  list_genera$Genus[i]
  
  # Get list of taxa from string match within the tips
  # list_sp_for_each_genera <- phylo.Ithomiini_treedata@phylo$tip.label[str_which(string = phylo.Ithomiini_treedata@phylo$tip.label, pattern = list_genera[i])]
  
  # Get list of taxa from string match within the data associated with the tips
  list_sp_for_each_genera <- phylo.Ithomiini_treedata@phylo$tip.label[phylo.Ithomiini_treedata@extraInfo$Genus %in% list_genera$Genus[i]]
  
  if(length(list_sp_for_each_genera) > 1) # Check if there is more than one species for this Genus
  {
    # Store nod for MRCA of Genera
    list_genera$MRCA_nod[i] <- phylobase::MRCA(phylo.Ithomiini_treedata@phylo, tip = list_sp_for_each_genera)
  } else { # Otherwise, get directly the tip index
    list_genera$MRCA_nod[i] <- which(phylo.Ithomiini_treedata@phylo$tip.label == list_sp_for_each_genera)
  }

}

## 3.2/ Find the MRCA nod for each of the 10 sub-tribe ####
Ithomiini_taxonomy_with_subtribes <- xlsx::read.xlsx("./input_data/Databases/Ithomiini_taxonomy_with_subtribes.xlsx", sheetIndex = 1)
Ithomiini_taxonomy_with_subtribes$sp_full <- paste(Ithomiini_taxonomy_with_subtribes$Genus, Ithomiini_taxonomy_with_subtribes$Species, sep = ".")

list.sp$Subtribe[match(Ithomiini_taxonomy_with_subtribes$sp_full, list.sp$sp_full_new)] <- as.character(Ithomiini_taxonomy_with_subtribes$Subtribe)
list.models <- left_join(x = list.models, y = Ithomiini_taxonomy_with_subtribes[, c("sp_full", "Subtribe")], by = c("sp_full_new" = "sp_full"))

# save(list.models, file = "./input_data/list.models.RData")
# save(list.sp, file = "./input_data/list.sp.RData")

Subtribe <- as.character(unique(Ithomiini_taxonomy_with_subtribes$Subtribe))
list.subtribes <-  as.data.frame(Subtribe)

for (i in seq_along(list.subtribes$Subtribe)) 
{
  list.subtribes$Subtribe[i]
  
  # Get list of taxa from string match within the tips
  # list_sp_for_each_genera <- phylo.Ithomiini_treedata@phylo$tip.label[str_which(string = phylo.Ithomiini_treedata@phylo$tip.label, pattern = list.subtribes[i])]
  
  # Get list of taxa from string match within the data associated with the tips
  list_sp_for_each_subtribe <- phylo.Ithomiini_treedata@phylo$tip.label[phylo.Ithomiini_treedata@extraInfo$Subtribe %in% list.subtribes$Subtribe[i]]
  
  if(length(list_sp_for_each_subtribe) > 1) # Check if there is more than one species for this Subtribe
  {
    # Store nod for MRCA of Genera
    list.subtribes$MRCA_nod[i] <- phylobase::MRCA(phylo.Ithomiini_treedata@phylo, tip = list_sp_for_each_subtribe)
  } else { # Otherwise, get directly the tip index
    list.subtribes$MRCA_nod[i] <- which(phylo.Ithomiini_treedata@phylo$tip.label == list_sp_for_each_subtribe)
  }
  
}

# Reorder following crown age in phylogeny
match(c("Melinaeina", "Mechanitina", "Methonina", "Tithoreina", "Athesitina", "Ithomiina", "Napeogenina", "Oleriina", "Godyridina", "Dircennia"), list.subtribes$Subtribe)
list.subtribes <- list.subtribes[match(c("Melinaeina", "Mechanitina", "Methonina", "Tithoreina", "Athesitina", "Ithomiina", "Napeogenina", "Oleriina", "Godyridina", "Dircennia"), list.subtribes$Subtribe), ]

# Add 3-letters symbol
list.subtribes$subtribe_symbol <- c("MEL", "MEC", "MET", "TIT", "ATH", "ITH", "NAP", "OLE", "GOD", "DIR")

save(list.subtribes, file = paste0("./input_data/list.subtribes.RData"))


### 4/ Plot phylogeny with labels ####

## 4.1/ Rectangular with all nod/branch labels ####

pdf(file = "./phylogenies/Full_phylo_rectangular_all_labels.pdf", height = 50, width = 20)
Ithomiini_phylo_plot <- ggtree(phylo.Ithomiini_treedata,
                               layout = "rectangular") +
  geom_tiplab(align = T,     # To align all tip label rather than move them aside each tips (useful for non-ultrametric tree where tips can be all over the x-axis)
              linetype = 0, 
              size = 2, 
              offset = 0.3,    # To move to the exterior. Useful if you wish to put mulitple layer of labels with different offsets.
              hjust = 0) +    # To center, left-align, right-aligned text
  geom_label(aes(label = node),
             fill = 'steelblue', size = 1) +
  geom_label(aes(x = branch, label = node),
             fill = 'lightgreen', size = 1)

for (i in seq_along(list.subtribes$Subtribe)) {
  Ithomiini_phylo_plot <- Ithomiini_phylo_plot +
    # Add a label to the median tip. Add a bar alongside all the clade tips
    geom_cladelabel(node = list.subtribes$MRCA_nod[i],
                    label = list.subtribes$Subtribe[i], 
                    align = T,                     # For non-ultrametric tree, to align veritcally all label, whever the tip finish
                    geom = 'text',     #  Tip of label, with or without rectangle
                    fill = NA,            # Color of the label background
                    color = 'black',                 # Color of the text and the bar
                    angle = 270,
                    hjust = 'center',
                    offset = 2.2,
                    offset.text = 0.2,
                    barsize = 2)             # Width of the bar
}

# Ithomiini_phylo_plot <- Ithomiini_phylo_plot +
#   layout_circular()

print(Ithomiini_phylo_plot)
dev.off()

Ithomiini_phylo_plot <- Ithomiini_phylo_plot +
  layout_circular()
print(Ithomiini_phylo_plot)

## 4.2/ Circular layout with colored subtribe labels around the tree ####

pdf(file = "./phylogenies/Full_phylo_circular_all_subtribes.pdf", height = 20, width = 20)
Ithomiini_phylo_plot <- ggtree(phylo.Ithomiini_treedata,
                               layout = "circular", size = 1.7)

angles <- c(277, 295, -52, 314, -39, -21,  22,  80, 145, 226)
offset.texts <- c(4.0, 1.5, 7.0, 4.0, 1.5, 1.5, 2.0, 2.0, 2.0, 2.0)

for (i in seq_along(list.subtribes$Subtribe)) {
  Ithomiini_phylo_plot <- Ithomiini_phylo_plot +
    # Add a label to the median tip. Add a bar alongside all the clade tips
      geom_cladelabel(node = list.subtribes$MRCA_nod[i],
                      label = list.subtribes$Subtribe[i], 
                      align = T,                       # For non-ultrametric tree, to align veritcally all label, whever the tip finish
                      geom = 'text',                   #  Tip of label, with or without rectangle
                      fill = NA,                       # Color of the label background
                      color =  tmaptools::get_brewer_pal("Spectral", n = 10, plot = F)[i],                 # Color of the text and the bar
                      fontsize = 10,
                      angle = angles[i],
                      # angle = 0,
                      hjust = 'center',
                      offset = 0.5,
                      offset.text = offset.texts[i],
                      barsize = 2) 

}

print(Ithomiini_phylo_plot)

dev.off()


## 4.3/ Circular layout with black subtribe labels with horizontal names ####

pdf(file = "./phylogenies/Full_phylo_circular_all_subtribes_final.pdf", height = 16, width = 20)
Ithomiini_phylo_plot <- ggtree(phylo.Ithomiini_treedata,
                               layout = "circular", size = 1.5)

list.subtribes$label <- paste0(list.subtribes$Subtribe, "\n(", list.subtribes$subtribe_symbol, ")")
# list.subtribes$custom_label <- as.character(list.subtribes$Subtribe)
# list.subtribes$custom_label[c(3:6, 10)] <- paste0(list.subtribes$Subtribe, "  (", list.subtribes$subtribe_symbol, ")")[c(3:6, 10)] 
offset.texts <- c(8.0, 6.0, 6.5, 6.5, 7.0, 6.0, 5.8, 6.5, 6.5, 7.0)
hjust <- c(0.5, 0.2, 0.2, 0.2, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)

for (i in seq_along(list.subtribes$Subtribe)) 
{
  Ithomiini_phylo_plot <- Ithomiini_phylo_plot +
    # Add a label to the median tip. Add a bar alongside all the clade tips
    geom_cladelabel(node = list.subtribes$MRCA_nod[i],
                    # label = list.subtribes$Subtribe[i], 
                    label = paste0('bold(',list.subtribes$Subtribe[i],')'), parse = T,
                    # label = paste0('bold(',list.subtribes$custom_label[i],')'), parse = T,
                    # label = list.subtribes$label[i],
                    # label = paste0('bold(',list.subtribes$label[i],')'), parse = T,
                    # label = bquote(atop(bold(.(list.subtribes$Subtribe[i])), ~ "(" ~ .(list.subtribes$subtribe_symbol[i]) ~ ")")), parse = T,
                    # label = paste0('bold(',list.subtribes$Subtribe[i],',', list.subtribes$subtribe_symbol[i],')'), parse = T,
                    align = T,                       # For non-ultrametric tree, to align vertically all label, whever the tip finish
                    geom = 'text',                   #  Tip of label, with or without rectangle
                    fill = NA,                       # Color of the label background
                    color =  c(tmaptools::get_brewer_pal("Spectral", n = 10, plot = F)[i], "black"),                 # Color of the text and the bar
                    fontsize = 12,
                    # fontface = 2,
                    angle = 0,
                    hjust = hjust[i],
                    offset = 2,
                    offset.text = offset.texts[i],
                    barsize = 25)
}
age <- max(Ithomiini_phylo_plot$data$x) - 5
while (age > 0)
{ 
  Ithomiini_phylo_plot <- Ithomiini_phylo_plot +
     geom_vline(xintercept = age, linetype = 92, size = 1, color = 'grey80')
  age  <- age - 5
}
# geom_vline(xintercept = 10, linetype = 29, size = 0.2, color = 'grey80') +
# geom_vline(xintercept = 20, linetype = 29, size = 0.2, color = 'grey80') +
# geom_vline(xintercept = root_age, linetype = 29, size = 0.2, color = 'grey80')

# Adapt margins
Ithomiini_phylo_plot <- Ithomiini_phylo_plot +
  theme(plot.margin = unit(c(-80, -20, -100, -20), "mm"),
        plot.background = element_rect(fill = "transparent",colour = NA),
        panel.border = element_rect(color = "transparent", fill = NA),
        panel.background = element_rect(color = NA, fill = "transparent"))

print(Ithomiini_phylo_plot)

dev.off()

?geom_cladelabel

### Tests to highlight one species ####

# Choose the species of focus
focus_sp <- list.sp$sp_full_new[2]

# Extract nod ID from taxa name
get_taxa_ID <- function(focus_taxa, phylo_treedata) 
{
  index <- which(phylo_treedata@phylo$tip.label == focus_taxa)
  return(index)
}

# Get list of ancestral branches from root to focal nod
taxa_ID <- get_taxa_ID(focus_taxa = focus_sp, phylo_treedata = phylo.Ithomiini_treedata)

?ancestor
ancestor(.data = phylo.Ithomiini_treedata, .node = taxa_ID)

branches_list <- ancestor(.data = phylo.Ithomiini_treedata, .node = 1)
branches_list <- c(branches_list, taxa_ID)
branches_color <- rep(x = NA, times = nrow(phylo.Ithomiini_treedata@extraInfo))
branches_color[branches_list] <- "red"

# Circular layout
pdf(file = "./phylogenies/Full_phylo_circular_one_sp.pdf", height = 20, width = 20)
Ithomiini_phylo_plot <- ggtree(phylo.Ithomiini_treedata,
                               layout = "circular", size = 1.7,
                               # color = branches_color
                               )

Ithomiini_phylo_plot <- Ithomiini_phylo_plot +
  geom_tree(color = branches_color, size = 10)

for (i in seq_along(list.subtribes$Subtribe)) {
  Ithomiini_phylo_plot <- Ithomiini_phylo_plot +
    # Add a label to the median tip. Add a bar alongside all the clade tips
    geom_cladelabel(node = list.subtribes$MRCA_nod[i],
                    label = NA, 
                    align = T,                       # For non-ultrametric tree, to align vertically all label, wherever the tip finish
                    geom = 'text',                   # Tip of label, with or without rectangle
                    fill = NA,                       # Color of the label background
                    color =  tmaptools::get_brewer_pal("Spectral", n = 10, plot = F)[i],                 # Color of the text and the bar
                    fontsize = 10,
                    angle = 0,
                    hjust = 'center',
                    offset = 0.5,
                    offset.text = 0,
                    barsize = 15) 
}

print(Ithomiini_phylo_plot)

dev.off()


### 5/ Generate phylo_inset for all OMUs ####

# When the species is not in the phylogeny. Find any species from the same Genus as replacement
taxonomy_update <- xlsx::read.xlsx(file = "./input_data/Taxonomy_update_2021.xlsx", sheetIndex = 1, stringsAsFactors = F)
load(file = "./input_data/list.models.RData")

# # Update species name in the list.model file
# 
# # Clean list model from old species names
# list.models$Updated_name <- F
# list.models$Genus_new <- as.character(list.models$Genus)
# list.models$Species_new <- as.character(list.models$Species)
# 
# for(i in 1:nrow(taxonomy_update))
# {
#   match_indices <- which((list.models$Genus %in% taxonomy_update$Genus_old[i]) & (list.models$Species %in% taxonomy_update$Species_old[i]))
#   list.models$Updated_name[match_indices] <- T
#   list.models$Genus_new[match_indices] <- taxonomy_update$Genus_new[i]
#   list.models$Species_new[match_indices] <- taxonomy_update$Species_new[i]
#   
# }
# 
# list.models$Species_new <- str_replace(string = list.models$Species_new, pattern = "AbsPhylo", replacement = "")
# list.models$Species_new[(list.models$Species_new == "nsp4")] <- "bifurcata"
# list.models$Updated_name[(list.models$Species_new == "bifurcata")] <- T
# list.models$sp_full_new <- paste(as.character(list.models$Genus_new), as.character(list.models$Species_new), sep = ".")
# list.models$Genus_new <- as.factor(list.models$Genus_new)
# list.models$Species_new <- as.factor(list.models$Species_new)
# 
# # save(list.models, file = paste0("./input_data/list.models.RData"))

# Clean environnement
rm(list = ls())


# Save tree_plot as a basis
load(file = "./phylogenies/phylo.Ithomiini_treedata.RData")
load(file = "./input_data/list.subtribes.RData")

# Generate base plot
Ithomiini_phylo_plot <- ggtree(phylo.Ithomiini_treedata,
                               layout = "circular", size = 1.7)

for (i in seq_along(list.subtribes$Subtribe)) {
  Ithomiini_phylo_plot <- Ithomiini_phylo_plot +
    # Add a label to the median tip. Add a bar alongside all the clade tips
    geom_cladelabel(node = list.subtribes$MRCA_nod[i],
                    label = NA, 
                    align = T,                       # For non-ultrametric tree, to align vertically all label, wherever the tip finish
                    geom = 'text',                   # Tip of label, with or without rectangle
                    fill = NA,                       # Color of the label background
                    color =  tmaptools::get_brewer_pal("Spectral", n = 10, plot = F)[i],                 # Color of the text and the bar
                    fontsize = 10,
                    angle = 0,
                    hjust = 'center',
                    offset = 0.5,
                    offset.text = 0,
                    barsize = 15) 
}

# print(Ithomiini_phylo_plot)

save(Ithomiini_phylo_plot, file = "./phylogenies/Ithomiini_phylo_plot.RData")

# Function to plot each OMU inset

generate_phylo_inset <- function(OMU)
{
  # # # Needed but we don't want to load them every time we use the function
  # load(file = "./phylogenies/phylo.Ithomiini_treedata.RData")
  # load(file = "./input_data/list.models.RData")
  
  # Choose the species of focus
  focus_sp <- list.models$sp_full_new[list.models$Tag.model == OMU]
  focus_genus <- as.character(list.models$Genus_new[list.models$Tag.model == OMU])
  
  # Check if the species is in the phylogeny
  if ((focus_sp %in% phylo.Ithomiini_treedata@phylo$tip.label))
  {
    # Extract nod number
    taxa_ID <- get_taxa_ID(focus_taxa = focus_sp, phylo_treedata = phylo.Ithomiini_treedata)
  } else  {
    # Replace by a species within the same genus
    taxa_ID <- which.max(focus_genus == phylo.Ithomiini_treedata@extraInfo$Genus)
  }
  
  # Get list of ancestral branches from root to focal nod
  # ?ancestor
  
  branches_list <- ancestor(.data = phylo.Ithomiini_treedata, .node = taxa_ID)
  branches_list <- c(branches_list, taxa_ID)
  branches_color <- rep(x = NA, times = nrow(phylo.Ithomiini_treedata@extraInfo))
  branches_color[branches_list] <- "red"
  
  load(file = "./phylogenies/Ithomiini_phylo_plot.RData")

  Ithomiini_phylo_plot_taxa <- Ithomiini_phylo_plot +
    geom_tree(color = branches_color, size = 10)
  
  return(Ithomiini_phylo_plot_taxa)
}

# Test with one OMU
OMU <- list.models$Tag.model[2]

Ithomiini_phylo_plot <- generate_phylo_inset(OMU = OMU)
print(Ithomiini_phylo_plot)

# Generate inset ggplot for all OMU
for (i in seq_along(list.models$Tag.model)) 
{
  unit <- list.models$Tag.model[i]
  
  Ithomiini_phylo_plot <- generate_phylo_inset(OMU = unit)
  
  # Plot
  pdf(file = paste0("./phylogenies/By_OMU/plots/Phylo_inset_",unit,".pdf"), height = 20, width = 20)
  print(Ithomiini_phylo_plot)
  dev.off()
  
  # Save
  saveRDS(Ithomiini_phylo_plot, file = paste0("./phylogenies/By_OMU/outputs/Phylo_inset_",unit,".rds"))
  
  if (i %% 10 == 0)
  {
    cat(paste0("\n", Sys.time()," ------ Plotting ", unit, " = Unit N°",i,"/",length(list.models$Tag.model)," ------\n")) 
  }
  
}

### 6/ Plot pretty distribution map with occurrences and inset for the phylogeny ####

## 6.0/ Load stuff ####

load(file = "./input_data/list.subtribes.RData")

pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")

grid_Mollweide_out <- readRDS(file = "./input_data/Map_stuff/grid_Mollweide_out.rds")
large_bg_mask_Mollweide <- readRDS(file = "./input_data/Map_stuff/large_bg_mask_Mollweide.rds")
bbox_sp_Mollweide <- readRDS(file = "./input_data/Map_stuff/bbox_sp_Mollweide.rds")

load(file = "./input_data/Map_stuff/country_borders.RData")
country_borders <- as(country_borders, "Spatial")

Mollweide_shp_projection <-  function(x) # Shp to project
{
  x_name <- deparse(substitute(x)) # Get the name of the initial shp as a character string
  
  new_shp <- spTransform(x, CRSobj = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
  
  # Generate new object with "_Mollweide" suffix in the global environment
  eval(call("<<-", as.name(paste0(x_name, "_Mollweide")), new_shp))
}

Mollweide_shp_projection(country_borders)

Mollweide_projection <- function(x) # Raster to project
{
  new_map <- projectRaster(from = x, 
                           method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                           crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                           alignOnly = F)
  return(new_map)
}

# Function to plot a unit map
{
  map_unit_Mollweide <- function(x,                                       # Raster to map
                                 unit_pts,                                # Shp with data points of ocurrences
                                    color_palette = pal_bl_red_Mannion,   # Color palette
                                    main_title,                           # Main title
                                    main_title_cex = 1.5,                 # Main title size
                                    
                                    xlim = c(-4600, 4600),    # Limit of plot on x-axis (Longitude)
                                    ylim = c(-4450, 3400),    # Limit of plot on y-axis (Latitude)
                                    axis_cex = 1.4,           # Axes size
                                    
                                    xlab = "",                # X-axis label
                                    ylab = "",                # Y-axis label
                                    x_axis_breaks = c(-3930, -2170, -420, 1500, 3050),            # X-axis tick breaks
                                    y_axis_breaks = c(-3650, -2450, -1220, 0, 1230, 2445),        # Y-axis tick breaks
                                    x_axis_labels = c("120°E", "100°E", "80°E", "60°E", "40°E"),      # X-axis tick labels
                                    y_axis_labels = c("30°S", "20°S", "10°S", "0°", "10°N", "20°N"),  # Y-axis tick labels
                                    
                                    occ_cex = 0.5,         # Size of occurrences
                                    occ_col = "#00000080",  # Color of occurrences
                                 
                                    legend_title,             # Legend title
                                    legend_title_cex = 1.4,   # Legend title size
                                    legend_title_x = -3550,   # Legend title x position
                                    legend_title_y = 430,     # Legend title y position
                                    legend_cex = 1.4,         # Legend size
                                    legend_breaks,            # Legend tick positions
                                    legend_location = c(-4100, -3800, -3950, 0),  # Legend position
                                    
                                    scale_bar_position = c(-2600, -4000),  # Scale bar position
                                    
                                    arrow_scale = 0.45,           # North arrow size
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
    
    # Add occurrence data
    plot(unit_pts, add = TRUE, col = occ_col, pch = 16, cex = occ_cex)
    
    # Add occurrence legend
    legend(legend = NA, pch = 16, cex = 1, x = -2700, y = -2600, bty = "n")
    text(labels = "Occurrences", cex = 1.05, x = -950, y = -2970, bty = "n", font = 2)

    
    # Add scale bar in legend
    scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = scale_bar_position, label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1.2)
    prettymapr::addnortharrow(pos = "bottomright", scale = arrow_scale, padin = arrow_padin, text.col = "#00000000")
    rangeBuilder::addRasterLegend(x, locs = legend_breaks, cex.axis = legend_cex, ramp = color_palette, ncolors = 200, border = T, location = legend_location)
    rangeBuilder::addRasterLegend(x, locs = legend_breaks, cex.axis = legend_cex, ramp = color_palette, ncolors = 200, border = T, location = legend_location)
    graphics::text(x = legend_title_x, y = legend_title_y, font = 2, cex = legend_title_cex, label = legend_title)
    
    # Add facet letter
    legend(legend = facet_letter, x = "bottomright", bty = "n", text.col = facet_letter_col,
           text.font = 2, cex = facet_letter_cex, inset = facet_letter_inset)
    
  }
}

# Function to plot a sp map

{
    map_sp_Mollweide <- function(x,                                    # Raster to map
                                 unit_pts,                             # Shp with data points of occurrences
                                 color_palette = pal_bl_red_Mannion,   # Color palette
                                 main_title,                           # Main title
                                 main_title_cex = 1.5,                 # Main title size
                                 
                                 xlim = c(-4600, 4600),    # Limit of plot on x-axis (Longitude)
                                 ylim = c(-4450, 3400),    # Limit of plot on y-axis (Latitude)
                                 axis_cex = 1.4,           # Axes size
                                 
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
                                 
                                 arrow_scale = 0.45,           # North arrow size
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
    
    # Add occurrence data
    occ_palette <- tmaptools::get_brewer_pal("Set1", n = 8, plot = F)
    plot(unit_pts, add = TRUE, col = "black", pch = 16, cex = 0.7)
    plot(unit_pts, add = TRUE, col = occ_palette[as.factor(unit_pts$Mimicry)], pch = 16, cex = 0.3)

    
    # Extract mimicry ring names
    mimicry_list <- levels(as.factor(unit_pts$Mimicry))
    N.mimicry <- length(mimicry_list)
    
    # Add occurrence legend
    text(labels = "Patterns", cex = 1.0, x = -1700, y = -2900 + 320*N.mimicry, bty = "n", font = 2)
    for (i in seq_along(mimicry_list)) 
    {
      legend(legend = mimicry_list[i], text.font = 2, pch = 16, cex = 0.8, pt.cex = 1.2, col = "black",
             x = -2700, y = -2750 + 320*N.mimicry - 320*i, bty = "n")
      legend(legend = mimicry_list[i], text.font = 2, pch = 16, cex = 0.8, col = occ_palette[i],
             x = -2700, y = -2750 + 320*N.mimicry - 320*i, bty = "n") 
    }

    # Add scale bar in legend
    scalebar(d = 2000, type = "line", lwd = 4, divs = 4, xy = scale_bar_position, label = c("", "2000 km", ""), adj = c(0.5, -0.8), font = 2, cex = 1.2)
    prettymapr::addnortharrow(pos = "bottomright", scale = arrow_scale, padin = arrow_padin, text.col = "#00000000")
    rangeBuilder::addRasterLegend(x, locs = legend_breaks, cex.axis = legend_cex, ramp = color_palette, ncolors = 200, border = T, location = legend_location)
    rangeBuilder::addRasterLegend(x, locs = legend_breaks, cex.axis = legend_cex, ramp = color_palette, ncolors = 200, border = T, location = legend_location)
    graphics::text(x = legend_title_x, y = legend_title_y, font = 2, cex = legend_title_cex, label = legend_title)
    
    # Add facet letter
    legend(legend = facet_letter, x = "bottomright", bty = "n", text.col = facet_letter_col,
           text.font = 2, cex = facet_letter_cex, inset = facet_letter_inset)
    
  }
}


## 6.1/ Per OMU ####

### Functions

# Function to retrieve type of model and parse name
extract_OMU_spec <- function(unit)
{
  load(file = "./input_data/list.rings.RData")
  
  ID <- which(list.models$Tag.model == unit)
  species_name <- paste(list.models$Genus_new[ID], list.models$Species_new[ID])
  
  # Update pattern name
  pattern <- as.character(list.models$Mimicry.model[ID])
  pattern <- list.rings$Mimicry_new[which(list.rings$Mimicry == pattern)]
  
  # Retrieve model type
  model_type <- as.character(list.models$initial_model_type[ID])
  
  model_type_list <- c("rasterized", "restricted", "complete")
  model_type_output <- c("No SDM (N < 6)", "Restricted SDM (5 < N < 31)", "Complete SDM (N > 30)")
  
  model_type_output <- model_type_output[which(model_type_list == model_type)]
  
  return(c(species_name, pattern, model_type_output))
}

# Function to generate the phylo inset adapted to the dimension of the final map. Version for pdf output (tiny)
generate_phylo_inset_tiny <- function(OMU)
{
  # # # Needed but we don't want to load them every time we use the function
  # load(file = "./phylogenies/phylo.Ithomiini_treedata.RData")
  # load(file = "./input_data/list.models.RData")
  # load(file = "./input_data/list.subtribes.RData")
  
  # Choose the species of focus
  focus_sp <- list.models$sp_full_new[list.models$Tag.model == OMU]
  focus_genus <- as.character(list.models$Genus_new[list.models$Tag.model == OMU])
  focus_subtribe <- as.character(list.models$Subtribe[list.models$Tag.model == OMU])
  
  # Check if the species is in the phylogeny
  if ((focus_sp %in% phylo.Ithomiini_treedata@phylo$tip.label))
  {
    # Extract nod number
    taxa_ID <- get_taxa_ID(focus_taxa = focus_sp, phylo_treedata = phylo.Ithomiini_treedata)
  } else  {
    # Replace by a species within the same genus
    taxa_ID <- which.max(focus_genus == phylo.Ithomiini_treedata@extraInfo$Genus)
  }
  
  # Get list of ancestral branches from root to focal nod
  # ?ancestor
  
  branches_list <- ancestor(.data = phylo.Ithomiini_treedata, .node = taxa_ID)
  branches_list <- c(branches_list, taxa_ID)
  branches_color <- rep(x = NA, times = nrow(phylo.Ithomiini_treedata@extraInfo))
  branches_color[branches_list] <- "red"
  
  # Define subtribe symbol and position
  
  # subtribe_list <- c("Melinaeina", "Mechanitina", "Methonina", "Tithoreina", "Athesitina", "Ithomiina", "Napeogenina", "Oleriina", "Godyridina", "Dircennia")
  # subtribe_symbol <- c("MEL", "MEC", "MET", "TIT", "ATH", "ITH", "NAP", "OLE", "GOD", "DIR")
  # subtribe_position_x <- c(1, 1, 1, 50, 1, 1, 1, 1, 1, 1)
  # subtribe_position_y <- c(1, 1, 1, 30, 1, 1, 1, 1, 1, 1)
  
  # Get info on subtribe
  subtribe_index <- which(list.subtribes$Subtribe %in% focus_subtribe)
  
  # load(file = "./phylogenies/Ithomiini_phylo_plot.RData")
  load(file = "./phylogenies/Ithomiini_phylo_plot_tiny.RData")
  
  Ithomiini_phylo_plot_taxa <- Ithomiini_phylo_plot_tiny +
    geom_tree(color = branches_color, size = 1) +
    geom_cladelabel(node = list.subtribes$MRCA_nod[subtribe_index],
                    label = paste0('bold(',list.subtribes$subtribe_symbol[subtribe_index],')'), parse = T,
                    # label = list.subtribes$subtribe_symbol[subtribe_index],
                    align = T,                       # For non-ultrametric tree, to align vertically all label, wherever the tip finish
                    geom = 'text',                   # Tip of label, with or without rectangle
                    fill = NA,                       # Color of the label background
                    color =  c(NA, tmaptools::get_brewer_pal("Spectral", n = 10, plot = F)[subtribe_index]),                 # Color of the text and the bar
                    # fontface = 2,
                    fontsize = 3.2,
                    angle = 0,
                    hjust = 'center',
                    offset = 22,
                    offset.text = 0,
                    barsize = 0) +
    # xlim(NA, 40) # +
    theme(plot.margin = unit(c(-3, -0, -3, -3), "mm"),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.border = element_rect(color = "transparent", fill = NA),
          panel.background = element_rect(color = NA, fill = "transparent"))
  
  # theme(clip = list(panel="off"))
  
  # annotate(geom = "text",
  #          cex = 4,
  #          fontface = 2,
  #          label = subtribe_symbol[subtribe_index],
  #          color = tmaptools::get_brewer_pal("Spectral", n = 10, plot = F)[subtribe_index],
  #          x = subtribe_position_x[subtribe_index],
  #          y = subtribe_position_y[subtribe_index])
  
  # print(Ithomiini_phylo_plot_tiny)
  # print(Ithomiini_phylo_plot_taxa)
  
  return(Ithomiini_phylo_plot_taxa)
}


###

# Load summary df
load(file = "./input_data/list.models.RData")

# Load stack of SDM outputs
All_OMU_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_OMU_stack_Jaccard.80.rds"))

# Generate base plot for pdf output (tiny)

{
  # load(file = "./phylogenies/phylo.Ithomiini_treedata.RData")
  # load(file = "./input_data/list.subtribes.RData")
  # 
  # Ithomiini_phylo_plot_tiny <- ggtree(phylo.Ithomiini_treedata,
  #                                     layout = "circular", size = 0.1)
  # 
  # for (i in seq_along(list.subtribes$Subtribe)) {
  #   Ithomiini_phylo_plot_tiny <- Ithomiini_phylo_plot_tiny +
  #     # Add a label to the median tip. Add a bar alongside all the clade tips
  #     geom_cladelabel(node = list.subtribes$MRCA_nod[i],
  #                     label = NA, 
  #                     align = T,                       # For non-ultrametric tree, to align vertically all label, wherever the tip finish
  #                     geom = 'text',                   # Tip of label, with or without rectangle
  #                     fill = NA,                       # Color of the label background
  #                     color =  tmaptools::get_brewer_pal("Spectral", n = 10, plot = F)[i],                 # Color of the text and the bar
  #                     fontsize = 10,
  #                     angle = 0,
  #                     hjust = 'center',
  #                     offset = 4,
  #                     offset.text = 0,
  #                     barsize = 3) 
  # }
  # 
  # print(Ithomiini_phylo_plot_tiny)
  # 
  # save(Ithomiini_phylo_plot_tiny, file = "./phylogenies/Ithomiini_phylo_plot_tiny.RData")
}

# Load base phylo plot for PDF output (tiny version)
load(file = "./phylogenies/Ithomiini_phylo_plot_tiny.RData")

### Plot map for each OMU one by one

for (i in seq_along(list.models$Tag.model))
{
  unit <- as.character(list.models$Tag.model[i])
  
  # Load unit occurrence points
  load(paste0("./input_data/Species_data/15/Spatial_Points_Objects/occurrences_", unit,".RData")) # Load Spatial object with occurrences
  Mollweide_shp_projection(unit.points)
  
  # Retrieve type of model and parse name
  unit_spec <- extract_OMU_spec(unit)
  unit_title <- paste(as.character(list.models$sp_full_new[i]), unit_spec[2], sep = ".")
  
  # Plot the map
  
  pdf(file = paste0("./maps/By_OMU/", unit, "/Final_map_with_phylo_", unit, ".pdf"), height = 6, width = 6)
  
  map <- Mollweide_projection(All_OMU_stack_Jaccard.80[[unit]])
  map@data@max <- 1
  
  # Map the OMU habitat suitability
  map_unit_Mollweide(x = map,
                     unit_pts = unit.points_Mollweide,
                     color_palette = c("#EDEDED", tmaptools::get_brewer_pal("RdYlGn", n = 199, plot = F)),
                     main_title = paste0(unit_spec[1], " (", unit_spec[2], ")\n", unit_spec[3]),
                     scale_bar_position = c(-2400, -4000),
                     legend_breaks = seq(0, 1, 0.2),
                     legend_title = "Habitat\n   suitability",
                     legend_title_x = -3400,
                     legend_title_y = 850)
  
  # Generate the phylo inset
  Ithomiini_phylo_plot_taxa <- generate_phylo_inset_tiny(OMU = unit)
  # print(Ithomiini_phylo_plot_taxa)
  
  # Ithomiini_phylo_plot_taxa <- ggplotify::as.grob(Ithomiini_phylo_plot_taxa)
  
  library(grid)
  # Create an appropriate viewport to add the ggplot as inset.
  vp.BottomRight <- viewport(height = unit(.22, "npc"), width = unit(0.22, "npc"), 
                             just = c("left","top"), 
                             y = 0.868, x = 0.705,
                             clip = "on")
  
  # plot the ggplot using the print command
  print(Ithomiini_phylo_plot_taxa, vp = vp.BottomRight)
  # print(Ithomiini_phylo_plot, vp = vp.BottomRight)
  
  dev.off()
  
  # Copy in Supplementary folder
  file.copy(from = paste0("./maps/By_OMU/", unit, "/Final_map_with_phylo_", unit, ".pdf"), 
            to = paste0("./supplementaries/Maps/OMU_maps/Range_map_", unit_title, ".pdf"),
            overwrite = T)
  
  if (i %% 10 == 0)
  {
    cat(paste0("\n", Sys.time()," ------ Plotting ", unit_title, " = Unit N°",i,"/",length(list.models$Tag.model)," ------\n")) 
  }
  
}

### Plot all in one pdf with each page as an OMU map

list.models <- list.models[order(as.character(list.models$Tag.model)),]

pdf(file = paste0("./maps/By_OMU/Final_map_with_phylo_all_OMU.pdf"), height = 6, width = 6)
for (i in seq_along(list.models$Tag.model))
# for (i in 1:3)
{
  unit <- as.character(list.models$Tag.model[i])
  
  # Load unit occurrence points
  load(paste0("./input_data/Species_data/15/Spatial_Points_Objects/occurrences_", unit,".RData")) # Load Spatial object with occurrences
  Mollweide_shp_projection(unit.points)
  
  # Retrieve type of model and parse name
  unit_spec <- extract_OMU_spec(unit)
  unit_title <- paste(as.character(list.models$sp_full_new[i]), unit_spec[2], sep = ".")
  
  # Plot the map
  map <- Mollweide_projection(All_OMU_stack_Jaccard.80[[unit]])
  map@data@max <- 1
  
  map_unit_Mollweide(x = map,
                     unit_pts = unit.points_Mollweide,
                     color_palette = c("#EDEDED", tmaptools::get_brewer_pal("RdYlGn", n = 199, plot = F)),
                     main_title = paste0(unit_spec[1], " (", unit_spec[2], ")\n", unit_spec[3]),
                     scale_bar_position = c(-2400, -4000),
                     legend_breaks = seq(0, 1, 0.2),
                     legend_title = "Habitat\n   suitability",
                     legend_title_x = -3400,
                     legend_title_y = 850)
  
  # Generate the phylo inset
  Ithomiini_phylo_plot_taxa <- generate_phylo_inset_tiny(OMU = unit)
  # print(Ithomiini_phylo_plot_taxa)
  
  # Ithomiini_phylo_plot_taxa <- ggplotify::as.grob(Ithomiini_phylo_plot_taxa)
  
  library(grid)
  # Create an appropriate viewport to add the ggplot as inset.
  vp.BottomRight <- viewport(height = unit(.22, "npc"), width = unit(0.22, "npc"), 
                             just = c("left","top"), 
                             y = 0.868, x = 0.705,
                             clip = "on")
  
  # plot the ggplot using the print command
  print(Ithomiini_phylo_plot_taxa, vp = vp.BottomRight)
  # print(Ithomiini_phylo_plot, vp = vp.BottomRight)
  
  if (i %% 10 == 0)
  {
    cat(paste0("\n", Sys.time()," ------ Plotting ", unit_title, " = Unit N°",i,"/",length(list.models$Tag.model)," ------\n")) 
  }
}
dev.off()

# Copy in Supplementary folder
file.copy(from = paste0("./maps/By_OMU/Final_map_with_phylo_all_OMU.pdf"), 
          to = paste0("./supplementaries/Maps/OMU_maps/0_Range_maps_all_OMU.pdf"),
          overwrite = T)


## 6.2/ Per species ####

### To improve if time:
# For species with no SDM, the legend of Habitat suitability is fucked up
# Put species name in italic

load(file = "./input_data/list.sp.RData")
load(file = "./input_data/list.models.RData")

### Functions

# Function to retrieve type of model and parse name

sp <- list.sp$Sp_full[9]

extract_sp_spec <- function(sp)
{
  ID <- which(list.sp$Sp_full == sp)
  species_name <- paste(list.sp$Genus_new[ID], list.sp$Species_new[ID])
  mimicry_nb <- list.sp$N.mimicry[ID]
  
  mimicry_list <- as.character(list.models$Mimicry.model[which(list.models$Sp_full == sp)])
  unit_list <- as.character(list.models$Tag.model[which(list.models$Sp_full == sp)])
  
  return(list(species_name, mimicry_nb, mimicry_list, unit_list))
}

# Function to generate occurrence spatial dataframe for sp
extract_sp_occ <- function(sp)
{
  # load(file = "./input_data/list.sp.RData")
  # load(file = "./input_data/list.models.RData")
  
  unit_list <- as.character(list.models$Tag.model[which(list.models$Sp_full == sp)])
  mimicry_list <- as.character(list.models$Mimicry.model[which(list.models$Sp_full == sp)])
  
  full_sp_occ <- NA
  for(i in seq_along(unit_list))
  {
    unit <- unit_list[i]
    
    # Load unit occurrence points
    load(paste0("./input_data/Species_data/15/Spatial_Points_Objects/occurrences_", unit,".RData")) # Load Spatial object with occurrences
    
    # Add mimicry ring identity to the df
    unit.points@data$Mimicry <- mimicry_list[i]
    
    # Merge all occurrence data in one
    if (i == 1) 
    {
      full_sp_occ <- unit.points
    } else {
      full_sp_occ <- rbind(full_sp_occ, unit.points)
    }
  }
  
  return(full_sp_occ)
}


# Function to plot each OMU inset. Version for pdf output (tiny)

i <- 9
sp <- as.character(list.sp$Sp_full[i]) 

# Function to generate the phylo inset adapted to the dimension of the final map
generate_phylo_inset_tiny_for_sp <- function(sp)
{
  # # # Needed but we don't want to load them every time we use the function
  # load(file = "./phylogenies/phylo.Ithomiini_treedata.RData")
  # load(file = "./input_data/list.sp.RData")
  # load(file = "./input_data/list.subtribes.RData")
  
  # Choose the species of focus
  focus_sp <- sp
  focus_genus <- as.character(list.sp$Genus_new[list.sp$Sp_full == sp])
  focus_subtribe <- as.character(list.sp$Subtribe[list.sp$Sp_full == sp])
  
  # Check if the species is in the phylogeny
  if ((focus_sp %in% phylo.Ithomiini_treedata@phylo$tip.label))
  {
    # Extract nod number
    taxa_ID <- get_taxa_ID(focus_taxa = focus_sp, phylo_treedata = phylo.Ithomiini_treedata)
  } else  {
    # Replace by a species within the same genus
    taxa_ID <- which.max(focus_genus == phylo.Ithomiini_treedata@extraInfo$Genus)
  }
  
  # Get list of ancestral branches from root to focal nod
  # ?ancestor
  
  branches_list <- ancestor(.data = phylo.Ithomiini_treedata, .node = taxa_ID)
  branches_list <- c(branches_list, taxa_ID)
  branches_color <- rep(x = NA, times = nrow(phylo.Ithomiini_treedata@extraInfo))
  branches_color[branches_list] <- "red"
  
  # Get info on subtribe
  subtribe_index <- which(list.subtribes$Subtribe %in% focus_subtribe)
  
  # load(file = "./phylogenies/Ithomiini_phylo_plot.RData")
  load(file = "./phylogenies/Ithomiini_phylo_plot_tiny.RData")
  
  Ithomiini_phylo_plot_taxa <- Ithomiini_phylo_plot_tiny +
    geom_tree(color = branches_color, size = 1) +
    geom_cladelabel(node = list.subtribes$MRCA_nod[subtribe_index],
                    label = paste0('bold(',list.subtribes$subtribe_symbol[subtribe_index],')'), parse = T,
                    # label = list.subtribes$subtribe_symbol[subtribe_index],
                    align = T,                       # For non-ultrametric tree, to align vertically all label, wherever the tip finish
                    geom = 'text',                   # Tip of label, with or without rectangle
                    fill = NA,                       # Color of the label background
                    color =  c(NA, tmaptools::get_brewer_pal("Spectral", n = 10, plot = F)[subtribe_index]),                 # Color of the text and the bar
                    # fontface = 2,
                    fontsize = 3.2,
                    angle = 0,
                    hjust = 'center',
                    offset = 22,
                    offset.text = 0,
                    barsize = 0) +
    # xlim(NA, 40) # +
    theme(plot.margin = unit(c(-3, -0, -3, -3), "mm"),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.border = element_rect(color = "transparent", fill = NA),
          panel.background = element_rect(color = NA, fill = "transparent"))
  
  # print(Ithomiini_phylo_plot_tiny)
  # print(Ithomiini_phylo_plot_taxa)
  
  return(Ithomiini_phylo_plot_taxa)
}

Ithomiini_phylo_plot_taxa <- generate_phylo_inset_tiny_for_sp(sp = sp)

###

# Load summary df
load(file = "./input_data/list.models.RData")
load(file = "./input_data/list.sp.RData")

# Load stack of SDM outputs
All_sp_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80.rds"))

# Generate base plot for pdf output (tiny)

{
  # load(file = "./phylogenies/phylo.Ithomiini_treedata.RData")
  # load(file = "./input_data/list.subtribes.RData")
  # 
  # Ithomiini_phylo_plot_tiny <- ggtree(phylo.Ithomiini_treedata,
  #                                     layout = "circular", size = 0.1)
  # 
  # for (i in seq_along(list.subtribes$Subtribe)) {
  #   Ithomiini_phylo_plot_tiny <- Ithomiini_phylo_plot_tiny +
  #     # Add a label to the median tip. Add a bar alongside all the clade tips
  #     geom_cladelabel(node = list.subtribes$MRCA_nod[i],
  #                     label = NA, 
  #                     align = T,                       # For non-ultrametric tree, to align vertically all label, wherever the tip finish
  #                     geom = 'text',                   # Tip of label, with or without rectangle
  #                     fill = NA,                       # Color of the label background
  #                     color =  tmaptools::get_brewer_pal("Spectral", n = 10, plot = F)[i],                 # Color of the text and the bar
  #                     fontsize = 10,
  #                     angle = 0,
  #                     hjust = 'center',
  #                     offset = 4,
  #                     offset.text = 0,
  #                     barsize = 3) 
  # }
  # 
  # print(Ithomiini_phylo_plot_tiny)
  # 
  # save(Ithomiini_phylo_plot_tiny, file = "./phylogenies/Ithomiini_phylo_plot_tiny.RData")
}

# Load base phylo plot for PDF output (tiny version)
load(file = "./phylogenies/Ithomiini_phylo_plot_tiny.RData")

### Plot map for each sp one by one ####

# i <- 1

for (i in seq_along(list.sp$Sp_full))
{
  sp <- as.character(list.sp$Sp_full[i])
  sp_title <- as.character(list.sp$sp_full_new[i])
  
  # Retrieve name and list of mimicry rings and OMU
  unit_spec <- extract_sp_spec(sp)
  
  sp_occ <- extract_sp_occ(sp)
  Mollweide_shp_projection(sp_occ)
  
  # Plot the map
  
  pdf(file = paste0("./maps/By_species/", sp, "/Final_map_with_phylo_", sp_title, ".pdf"), height = 6, width = 6)
  
  map <- Mollweide_projection(All_sp_stack_Jaccard.80[[sp]])
  map@data@max <- 1
  
  # Main title choice
  if(unit_spec[2] == 1)  # In case of monomorphic species. Pattern must be singular
  {
    main_title <- paste0(unit_spec[1], "\n(", unit_spec[2], " pattern)")
  } else {              # In case of polymorphic species. Patterns must be plural
    main_title <- paste0(unit_spec[1], "\n(", unit_spec[2], " patterns)")
  }
  
  
  # Map the OMU habitat suitability
  map_sp_Mollweide(x = map,
                     unit_pts = sp_occ_Mollweide,
                     color_palette = c("#EDEDED", tmaptools::get_brewer_pal("RdYlGn", n = 199, plot = F)),
                     
                     main_title = main_title,
                     
                     #### Try to put title for species in italic
                     # substitute(paste(italic('p value'), " = 0.01")))
                     # paste0('bold(',list.subtribes$subtribe_symbol[subtribe_index],')'), parse = T,
                     # expression(italic(.(unit_spec[1])))
                     # main_title = parse(text = paste0('italic(',unit_spec[1], ')')),
                     # main_title = bquote(atop(italic(.(unit_spec[1])), ~ "(" ~ .(unit_spec[2]) ~ "patterns)")),
                     
                     scale_bar_position = c(-2400, -4000),
                     legend_breaks = seq(0, 1, 0.2),
                     legend_title = "Habitat\n   suitability",
                     legend_title_x = -3400,
                     legend_title_y = 850)
  

  
  
  # Generate the phylo inset
  Ithomiini_phylo_plot_taxa <- generate_phylo_inset_tiny_for_sp(sp = sp)
  # print(Ithomiini_phylo_plot_taxa)
  
  # Ithomiini_phylo_plot_taxa <- ggplotify::as.grob(Ithomiini_phylo_plot_taxa)
  
  library(grid)
  # Create an appropriate viewport to add the ggplot as inset.
  vp.BottomRight <- viewport(height = unit(.22, "npc"), width = unit(0.22, "npc"), 
                             just = c("left","top"), 
                             y = 0.868, x = 0.705,
                             clip = "on")
  
  # plot the ggplot using the print command
  print(Ithomiini_phylo_plot_taxa, vp = vp.BottomRight)
  # print(Ithomiini_phylo_plot, vp = vp.BottomRight)
  
  dev.off()
  
  # Copy in Supplementary folder
  file.copy(from = paste0("./maps/By_species/", sp, "/Final_map_with_phylo_", sp_title, ".pdf"), 
            to = paste0("./supplementaries/Maps/Species_maps/Range_maps_", sp_title, ".pdf"),
            overwrite = T)
  
  if (i %% 10 == 0)
  {
    cat(paste0("\n", Sys.time()," ------ Plotting ", sp_title, " = Species N°",i,"/",length(list.sp$Sp_full)," ------\n")) 
  }
  
}

### Plot all in one pdf with each page as a sp map ####

# Reorder in alphabetic order
list.sp <- list.sp[order(as.character(list.sp$Sp_full)),]

pdf(file = paste0("./maps/By_species/Final_map_with_phylo_all_species.pdf"), height = 6, width = 6)
for (i in seq_along(list.sp$Sp_full))
{
  sp <- as.character(list.sp$Sp_full[i])
  sp_title <- as.character(list.sp$sp_full_new[i])
  
  # Retrieve name and list of mimicry rings and OMU
  unit_spec <- extract_sp_spec(sp)
  
  sp_occ <- extract_sp_occ(sp)
  Mollweide_shp_projection(sp_occ)
  
  # Plot the map
  map <- Mollweide_projection(All_sp_stack_Jaccard.80[[sp]])
  map@data@max <- 1
  
  # Main title choice
  if(unit_spec[2] == 1)  # In case of monomorphic species. Pattern must be singular
  {
    main_title <- paste0(unit_spec[1], "\n(", unit_spec[2], " pattern)")
  } else {              # In case of polymorphic species. Patterns must be plural
    main_title <- paste0(unit_spec[1], "\n(", unit_spec[2], " patterns)")
  }
  
  # Map the OMU habitat suitability
  map_sp_Mollweide(x = map,
                   unit_pts = sp_occ_Mollweide,
                   color_palette = c("#EDEDED", tmaptools::get_brewer_pal("RdYlGn", n = 199, plot = F)),
                   
                   main_title = main_title,
                   
                   #### Try to put title for species in italic
                   # substitute(paste(italic('p value'), " = 0.01")))
                   # paste0('bold(',list.subtribes$subtribe_symbol[subtribe_index],')'), parse = T,
                   # expression(italic(.(unit_spec[1])))
                   # main_title = parse(text = paste0('italic(',unit_spec[1], ')')),
                   # main_title = bquote(atop(italic(.(unit_spec[1])), ~ "(" ~ .(unit_spec[2]) ~ "patterns)")),
                   
                   scale_bar_position = c(-2400, -4000),
                   legend_breaks = seq(0, 1, 0.2),
                   legend_title = "Habitat\n   suitability",
                   legend_title_x = -3400,
                   legend_title_y = 850)
  
 # Generate the phylo inset
  Ithomiini_phylo_plot_taxa <- generate_phylo_inset_tiny_for_sp(sp = sp)
  # print(Ithomiini_phylo_plot_taxa)
  
  # Ithomiini_phylo_plot_taxa <- ggplotify::as.grob(Ithomiini_phylo_plot_taxa)
  
  library(grid)
  # Create an appropriate viewport to add the ggplot as inset.
  vp.BottomRight <- viewport(height = unit(.22, "npc"), width = unit(0.22, "npc"), 
                             just = c("left","top"), 
                             y = 0.868, x = 0.705,
                             clip = "on")
  
  # plot the ggplot using the print command
  print(Ithomiini_phylo_plot_taxa, vp = vp.BottomRight)
  # print(Ithomiini_phylo_plot, vp = vp.BottomRight)
  
  if (i %% 10 == 0)
  {
    cat(paste0("\n", Sys.time()," ------ Plotting ", sp_title, " = Species N°",i,"/",length(list.sp$Sp_full)," ------\n")) 
  }
  
}
dev.off()

# Copy in Supplementary folder
file.copy(from = paste0("./maps/By_species/Final_map_with_phylo_all_species.pdf"), 
          to = paste0("./supplementaries/Maps/Species_maps/0_Range_maps_all_species.pdf"),
          overwrite = T)

### 6.3/ Per mimicry ring ####

# Per ring, with a heatmap around the phylogeny showing where the taxa are found (or just the branches for those taxa ?)

### With a butterfly as illustration ?!

# Should I add all occurrences ? 

load(file = "./input_data/list.sp.RData")
load(file = "./input_data/list.models.RData")

# Create a summary df for mimicry rings

list.rings <- data.frame(Mimicry = unique(list.models$Mimicry.model))
list.rings$Mimicry <- list.rings$Mimicry[order(list.rings$Mimicry)]
list.rings$ID_ring <- seq_along(list.rings$Mimicry)                            # Mimicry ring ID
list.rings$N.species <- table(list.models$Mimicry.model)[list.rings$Mimicry]   # Number of species per rings
list.rings <- list.models %>% 
  group_by(Mimicry.model) %>% 
  summarise(Total.obs = sum(N.obs_15m_used)) %>% 
  left_join(x = list.rings, y = ., by = c("Mimicry" = "Mimicry.model"))
list.rings$N.obs <- list.rings_temp$Total.obs                                  # Number of observations

# Clean mimicry name
list.rings$Mimicry_new <- as.character(list.rings$Mimicry)
list.rings$Mimicry_new[which(list.rings$Mimicry_new == "BANJANAM")] <- "BANJANA-M"
list.rings$Mimicry_new[which(list.rings$Mimicry_new == "THABENAF")] <- "THABENA-F"
list.rings$Mimicry_new[which(list.rings$Mimicry_new == "TICIDAM")] <- "TICIDA-M"

save(list.rings, file = "./input_data/list.rings.RData")

load(file = "./input_data/list.rings.RData")

### Functions

# Function to retrieve type of model and parse name

i <- 11
ring <- list.rings$Mimicry[i]

# Function to generate occurrence spatial dataframe for sp
extract_ring_occ <- function(ring)
{
  # load(file = "./input_data/list.rings.RData")
  # load(file = "./input_data/list.models.RData")
  
  unit_list <- as.character(list.models$Tag.model[which(list.models$Mimicry.model == ring)])
  species_list <- as.character(list.models$Sp_full[which(list.models$Mimicry.model == ring)])
  
  full_ring_occ <- NA
  for(i in seq_along(unit_list))
  {
    unit <- unit_list[i]
    
    # Load unit occurrence points
    load(paste0("./input_data/Species_data/15/Spatial_Points_Objects/occurrences_", unit,".RData")) # Load Spatial object with occurrences
    
    # Add species identity to the df
    unit.points@data$Species <- species_list[i]
    
    # Merge all occurrence data in one
    if (i == 1) 
    {
      full_ring_occ <- unit.points
    } else {
      full_ring_occ <- rbind(full_ring_occ, unit.points)
    }
  }
  
  return(full_ring_occ)
}

phylo.Ithomiini$tip.label

# Load base phylo plot for PDF output (tiny version)
load(file = "./phylogenies/Ithomiini_phylo_plot_tiny.RData")

identical(phylo.Ithomiini$tip.label, Ithomiini_phylo_plot_tiny$data$label[Ithomiini_phylo_plot_tiny$data$isTip])

# Generate a data.frame of mimicry ring belonging
heatmapData <- data.frame(label = Ithomiini_phylo_plot_tiny$data$label[Ithomiini_phylo_plot_tiny$data$isTip])
row.names(heatmapData) <- heatmapData$label
for (i in seq_along(list.rings$Mimicry))
{
  ring <- as.character(list.rings$Mimicry[i])
  species_list <- as.character(list.models$Sp_full[which(list.models$Mimicry.model == ring)])
  ring_membership <- heatmapData$label %in% species_list
  heatmapData <- cbind(heatmapData, ring_membership)
}
heatmapData <- heatmapData[, -1]
names(heatmapData) <- as.character(list.rings$Mimicry)

save(heatmapData, file = "./input_data/ring_membership_table.RData")

load(file = "./input_data/ring_membership_table.RData")

i <- 11
ring <- list.rings$Mimicry[i]

# Function to generate the phylo inset adapted to the dimension of the final map
# With a heatmap for mimicry membership
generate_phylo_inset_tiny_for_ring_heatmap <- function(ring)
{
  # # # Needed but we don't want to load them every time we use the function
  # load(file = "./phylogenies/phylo.Ithomiini_treedata.RData")
  # load(file = "./input_data/list.rings.RData")
  # load(file = "./input_data/list.subtribes.RData")
  load(file = "./input_data/ring_membership_table.RData")
  
  species_list <- as.character(list.models$Sp_full[which(list.models$Mimicry.model == ring)])
  
  ## Only useful is one ring has no species in the phylogeny
  # # Check if the species is in the phylogeny
  # if ((focus_sp %in% phylo.Ithomiini_treedata@phylo$tip.label))
  # {
  #   # Extract nod number
  #   taxa_ID <- get_taxa_ID(focus_taxa = focus_sp, phylo_treedata = phylo.Ithomiini_treedata)
  # } else  {
  #   # Replace by a species within the same genus
  #   taxa_ID <- which.max(focus_genus == phylo.Ithomiini_treedata@extraInfo$Genus)
  # }
  
  # # Get list of ancestral branches from root to focal nod
  # # ?ancestor
  # 
  # branches_list <- ancestor(.data = phylo.Ithomiini_treedata, .node = taxa_ID)
  # branches_list <- c(branches_list, taxa_ID)
  # branches_color <- rep(x = NA, times = nrow(phylo.Ithomiini_treedata@extraInfo))
  # branches_color[branches_list] <- "red"
  
  # # Get info on subtribe
  # subtribe_index <- which(list.subtribes$Subtribe %in% focus_subtribe)
  
  # Extract the df for heatmap
  heatmapData_ring <- heatmapData[, which(names(heatmapData) == ring), drop = FALSE]
  
  # load(file = "./phylogenies/Ithomiini_phylo_plot.RData")
  load(file = "./phylogenies/Ithomiini_phylo_plot_tiny.RData")
  
  Ithomiini_phylo_plot_taxa <- gheatmap(p = Ithomiini_phylo_plot_tiny,                   # A tree plot (not a data_tree_obj)
           data = heatmapData_ring,    # Must be a dataframe with row = tips. Columns = data. Cell values as character. Tip.label as row.names !!!
           width = 0.1,                  # Total width of heatmap, compare to width of tree
           offset = 5,                # Offset to allow to add multiple heatmaps with different offset aside. (Easier to add multiple variable in a single heatmap)
           color = "gray10",            # Color of heatmap cell border
           colnames = FALSE,
           hjust = 0.5,                # Horizontal adjustment for column names
           font.size = 3,              # Font size for columns names  
           legend_title = NA) +        # Legend title for the entire heatmap        
    scale_fill_manual(values = c(NA, "black")) + # To set manually the colors
    theme(legend.position = "none") +
  
  # Ithomiini_phylo_plot_taxa <- Ithomiini_phylo_plot_taxa +
  #   geom_tree(color = branches_color, size = 1) +
    
    # xlim(NA, 40) # +
    theme(plot.margin = unit(c(0, -0, 0, 0), "mm"),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.border = element_rect(color = "transparent", fill = NA),
          panel.background = element_rect(color = NA, fill = "transparent"))
  
  # print(Ithomiini_phylo_plot_tiny)
  # print(Ithomiini_phylo_plot_taxa)
  
  return(Ithomiini_phylo_plot_taxa)
}


# Function to generate the phylo inset adapted to the dimension of the final map
# With geom_tippoint to show ring membership
generate_phylo_inset_tiny_for_ring_tippoint <- function(ring)
{
  # # # Needed but we don't want to load them every time we use the function
  # load(file = "./phylogenies/phylo.Ithomiini_treedata.RData")
  # load(file = "./input_data/list.rings.RData")
  # load(file = "./input_data/list.subtribes.RData")
  load(file = "./input_data/ring_membership_table.RData")
  
  # species_list <- as.character(list.models$Sp_full[which(list.models$Mimicry.model == ring)])
  
  ## Only useful is one ring has no species in the phylogeny
  # # Check if the species is in the phylogeny
  # if ((focus_sp %in% phylo.Ithomiini_treedata@phylo$tip.label))
  # {
  #   # Extract nod number
  #   taxa_ID <- get_taxa_ID(focus_taxa = focus_sp, phylo_treedata = phylo.Ithomiini_treedata)
  # } else  {
  #   # Replace by a species within the same genus
  #   taxa_ID <- which.max(focus_genus == phylo.Ithomiini_treedata@extraInfo$Genus)
  # }
  
  # # Get list of ancestral branches from root to focal nod
  # # ?ancestor
  # 
  # branches_list <- ancestor(.data = phylo.Ithomiini_treedata, .node = taxa_ID)
  # branches_list <- c(branches_list, taxa_ID)
  # branches_color <- rep(x = NA, times = nrow(phylo.Ithomiini_treedata@extraInfo))
  # branches_color[branches_list] <- "red"
  
  # # Get info on subtribe
  # subtribe_index <- which(list.subtribes$Subtribe %in% focus_subtribe)
  
  # Extract the df for heatmap
  species_list_phylo <- row.names(heatmapData)[heatmapData[, which(names(heatmapData) == ring)]]
  
  # load(file = "./phylogenies/Ithomiini_phylo_plot.RData")
  load(file = "./phylogenies/Ithomiini_phylo_plot_tiny.RData")
  
  Ithomiini_phylo_plot_taxa <- Ithomiini_phylo_plot_tiny +
    geom_tippoint(aes(subset = label %in% species_list_phylo),
                  shape = 21, fill = "red", col = "black", size = 1.5) +

    # scale_fill_manual(values = c(NA, "black")) + # To set manually the colors
    # theme(legend.position = "none") +
    
    # Ithomiini_phylo_plot_taxa <- Ithomiini_phylo_plot_taxa +
    #   geom_tree(color = branches_color, size = 1) +
    
    # xlim(NA, 40) # +
    theme(plot.margin = unit(c(3, 3, 3, 3), "mm"),
          plot.background = element_rect(fill = "transparent",colour = NA),
          panel.border = element_rect(color = "transparent", fill = NA),
          panel.background = element_rect(color = NA, fill = "transparent"))
  
  # print(Ithomiini_phylo_plot_tiny)
  # print(Ithomiini_phylo_plot_taxa)
  
  return(Ithomiini_phylo_plot_taxa)
}

Ithomiini_phylo_plot_taxa <- generate_phylo_inset_tiny_for_ring_heatmap(ring = ring)
Ithomiini_phylo_plot_taxa <- generate_phylo_inset_tiny_for_ring_tippoint(ring = ring)

###

# Load summary df
load(file = "./input_data/list.models.RData")
load(file = "./input_data/list.sp.RData")

# Load stack of SDM outputs
All_ring_proba_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_Jaccard.80.RData"))
All_ring_rich_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_Jaccard.80.RData"))

# Load base phylo plot for PDF output (tiny version)
load(file = "./phylogenies/Ithomiini_phylo_plot_tiny.RData")

### Plot range map and richness map for each ring one by one ####

# i <- 3

for (i in seq_along(list.rings$Mimicry))
{
  ring <- as.character(list.rings$Mimicry[i])
  ring_title <- as.character(list.rings$Mimicry_new[i])
  
  ring_occ <- extract_ring_occ(ring)
  Mollweide_shp_projection(ring_occ)
  
  # Plot the range map
  
  pdf(file = paste0("./maps/By_ring/Range_maps/Final_map_with_phylo_", ring, ".pdf"), height = 6, width = 6)
  
  map <- Mollweide_projection(All_ring_proba_stack_Jaccard.80[[ring]])
  map@data@max <- 1
  
  # Map the ring habitat suitability/range
  map_unit_Mollweide(x = map,
                   unit_pts = ring_occ_Mollweide,
                   color_palette = c("#EDEDED", tmaptools::get_brewer_pal("RdYlGn", n = 199, plot = F)),
                   
                   main_title = paste0("Mimicry ring ", ring_title, "\n(", list.rings$N.species[i], " species)"),
                   
                   #### Try to put title for species in italic
                   # substitute(paste(italic('p value'), " = 0.01")))
                   # paste0('bold(',list.subtribes$subtribe_symbol[subtribe_index],')'), parse = T,
                   # expression(italic(.(unit_spec[1])))
                   # main_title = parse(text = paste0('italic(',unit_spec[1], ')')),
                   # main_title = bquote(atop(italic(.(unit_spec[1])), ~ "(" ~ .(unit_spec[2]) ~ "patterns)")),
                   
                   occ_cex = 0.5, occ_col = "black",
                   
                   scale_bar_position = c(-2400, -4000),
                   legend_breaks = seq(0, 1, 0.2),
                   legend_title = "Habitat\n   suitability",
                   legend_title_x = -3400,
                   legend_title_y = 850)
  
  # Generate the phylo inset
  # Ithomiini_phylo_plot_taxa <- generate_phylo_inset_tiny_for_ring_heatmap(ring = ring)
  Ithomiini_phylo_plot_taxa <- generate_phylo_inset_tiny_for_ring_tippoint(ring = ring)
  # print(Ithomiini_phylo_plot_taxa)
  
  # Ithomiini_phylo_plot_taxa <- ggplotify::as.grob(Ithomiini_phylo_plot_taxa)
  
  library(grid)
  # Create an appropriate viewport to add the ggplot as inset.
  vp.BottomRight <- viewport(height = unit(.22, "npc"), width = unit(0.22, "npc"), 
                             just = c("left","top"), 
                             y = 0.868, x = 0.705,
                             clip = "on")
  
  # plot the ggplot using the print command
  print(Ithomiini_phylo_plot_taxa, vp = vp.BottomRight)
  # print(Ithomiini_phylo_plot, vp = vp.BottomRight)
  
  dev.off()
  
  # Copy in Supplementary folder
  file.copy(from = paste0("./maps/By_ring/Range_maps/Final_map_with_phylo_", ring, ".pdf"), 
            to = paste0("./supplementaries/Maps/Mimicry_ring_maps/Range_maps/Range_map_", ring, ".pdf"),
            overwrite = T)
  
  # Plot the richness map
  
  pdf(file = paste0("./maps/By_ring/Richness_maps/Final_map_with_phylo_", ring, ".pdf"), height = 6, width = 6)
  
  map <- Mollweide_projection(All_ring_rich_stack_Jaccard.80[[ring]])

  # Define legend breaks depending on max value
  max <- map@data@max
  if (max <= 1) {map@data@max <- 1; legend_breaks <- seq(0, 1, 0.2)}
  if ((max > 1) & (max <= 2)) {map@data@max <- 2; legend_breaks <- seq(0, 2, 0.5)}
  if ((max > 2) & (max <= 5)) {map@data@max <- ceiling(map@data@max); legend_breaks <- seq(0, ceiling(map@data@max), 1)}
  if ((max > 5) & (max <= 10)) {map@data@max <- ceiling(map@data@max); legend_breaks <- seq(0, ceiling(map@data@max), 2)}
  if ((max > 10) & (max <= 15)) {map@data@max <- ceiling(map@data@max); legend_breaks <- seq(0, ceiling(map@data@max), 3)}
  if ((max > 15)) {legend_breaks <- seq(0, max, 5)}
  
  # Map the ring richness map
  map_unit_Mollweide(x = map,
                     unit_pts = ring_occ_Mollweide,
                     color_palette = pal_bl_red_Mannion,
                     
                     main_title = paste0("Mimicry ring ", ring_title, "\n(", list.rings$N.species[i], " species)"),
                     
                     #### Try to put title for species in italic
                     # substitute(paste(italic('p value'), " = 0.01")))
                     # paste0('bold(',list.subtribes$subtribe_symbol[subtribe_index],')'), parse = T,
                     # expression(italic(.(unit_spec[1])))
                     # main_title = parse(text = paste0('italic(',unit_spec[1], ')')),
                     # main_title = bquote(atop(italic(.(unit_spec[1])), ~ "(" ~ .(unit_spec[2]) ~ "patterns)")),
                     
                     occ_cex = 0.35,
                     
                     scale_bar_position = c(-2400, -4000),
                     legend_breaks = legend_breaks,
                     legend_title = "Number of\nspecies ",
                     legend_title_x = -3200,
                     legend_title_y = 850)
  
  # Generate the phylo inset
  # Ithomiini_phylo_plot_taxa <- generate_phylo_inset_tiny_for_ring_heatmap(ring = ring)
  Ithomiini_phylo_plot_taxa <- generate_phylo_inset_tiny_for_ring_tippoint(ring = ring)
  # print(Ithomiini_phylo_plot_taxa)
  
  # Ithomiini_phylo_plot_taxa <- ggplotify::as.grob(Ithomiini_phylo_plot_taxa)
  
  library(grid)
  # Create an appropriate viewport to add the ggplot as inset.
  vp.BottomRight <- viewport(height = unit(.22, "npc"), width = unit(0.22, "npc"), 
                             just = c("left","top"), 
                             y = 0.868, x = 0.705,
                             clip = "on")
  
  # plot the ggplot using the print command
  print(Ithomiini_phylo_plot_taxa, vp = vp.BottomRight)
  # print(Ithomiini_phylo_plot, vp = vp.BottomRight)
  
  dev.off()
  
  # Copy in Supplementary folder
  file.copy(from = paste0("./maps/By_ring/Richness_maps/Final_map_with_phylo_", ring, ".pdf"), 
            to = paste0("./supplementaries/Maps/Mimicry_ring_maps/Richness_maps/Richness_map_", ring, ".pdf"),
            overwrite = T)
  
  # Check run
  if (i %% 5 == 0)
  {
    cat(paste0("\n", Sys.time()," ------ Plotting ", ring_title, " = Mimicry ring N°",i,"/",length(list.rings$Mimicry)," ------\n")) 
  }
  
}


### Plot all in one pdf with each page as a mimicry ring range map ####

# Reorder in alphabetic order
list.rings <- list.rings[order(as.character(list.rings$Mimicry_new)),]

pdf(file = paste0("./maps/By_ring/Range_maps/0_Final_range_map_with_phylo_all_rings.pdf"), height = 6, width = 6)

for (i in seq_along(list.rings$Mimicry))
{
  ring <- as.character(list.rings$Mimicry[i])
  ring_title <- as.character(list.rings$Mimicry_new[i])
  
  ring_occ <- extract_ring_occ(ring)
  Mollweide_shp_projection(ring_occ)
  
  # Plot the map
  
  map <- Mollweide_projection(All_ring_proba_stack_Jaccard.80[[ring]])
  map@data@max <- 1
  
  map_unit_Mollweide(x = map,
                     unit_pts = ring_occ_Mollweide,
                     color_palette = c("#EDEDED", tmaptools::get_brewer_pal("RdYlGn", n = 199, plot = F)),
                     
                     main_title = paste0("Mimicry ring ", ring_title, "\n(", list.rings$N.species[i], " species)"),
                     
                     #### Try to put title for species in italic
                     # substitute(paste(italic('p value'), " = 0.01")))
                     # paste0('bold(',list.subtribes$subtribe_symbol[subtribe_index],')'), parse = T,
                     # expression(italic(.(unit_spec[1])))
                     # main_title = parse(text = paste0('italic(',unit_spec[1], ')')),
                     # main_title = bquote(atop(italic(.(unit_spec[1])), ~ "(" ~ .(unit_spec[2]) ~ "patterns)")),
                     
                     occ_cex = 0.5, occ_col = "black",
                     
                     scale_bar_position = c(-2400, -4000),
                     legend_breaks = seq(0, 1, 0.2),
                     legend_title = "Habitat\n   suitability",
                     legend_title_x = -3400,
                     legend_title_y = 850)
  
  # Generate the phylo inset
  # Ithomiini_phylo_plot_taxa <- generate_phylo_inset_tiny_for_ring_heatmap(ring = ring)
  Ithomiini_phylo_plot_taxa <- generate_phylo_inset_tiny_for_ring_tippoint(ring = ring)
  # print(Ithomiini_phylo_plot_taxa)
  
  # Ithomiini_phylo_plot_taxa <- ggplotify::as.grob(Ithomiini_phylo_plot_taxa)
  
  library(grid)
  # Create an appropriate viewport to add the ggplot as inset.
  vp.BottomRight <- viewport(height = unit(.22, "npc"), width = unit(0.22, "npc"), 
                             just = c("left","top"), 
                             y = 0.868, x = 0.705,
                             clip = "on")
  
  # plot the ggplot using the print command
  print(Ithomiini_phylo_plot_taxa, vp = vp.BottomRight)
  # print(Ithomiini_phylo_plot, vp = vp.BottomRight)
  
  if (i %% 5 == 0)
  {
    cat(paste0("\n", Sys.time()," ------ Plotting ", ring_title, " = Mimicry ring N°",i,"/",length(list.rings$Mimicry)," ------\n")) 
  }
}
dev.off()

# Copy in Supplementary folder
file.copy(from = paste0("./maps/By_ring/Range_maps/0_Final_range_map_with_phylo_all_rings.pdf"), 
          to = paste0("./supplementaries/Maps/Mimicry_ring_maps/Range_maps/0_Range_maps_all_rings.pdf"),
          overwrite = T)

### Plot all in one pdf with each page as a mimicry richness map ####

# Reorder in alphabetic order
list.rings <- list.rings[order(as.character(list.rings$Mimicry_new)),]

pdf(file = paste0("./maps/By_ring/Richness_maps/0_Final_richness_map_with_phylo_all_rings.pdf"), height = 6, width = 6)

for (i in seq_along(list.rings$Mimicry))
{
  ring <- as.character(list.rings$Mimicry[i])
  ring_title <- as.character(list.rings$Mimicry_new[i])
  
  ring_occ <- extract_ring_occ(ring)
  Mollweide_shp_projection(ring_occ)
  
  # Plot the richness map
  
  map <- Mollweide_projection(All_ring_rich_stack_Jaccard.80[[ring]])
  
  # Define legend breaks depending on max value
  max <- map@data@max
  if (max <= 1) {map@data@max <- 1; legend_breaks <- seq(0, 1, 0.2)}
  if ((max > 1) & (max <= 2)) {map@data@max <- 2; legend_breaks <- seq(0, 2, 0.5)}
  if ((max > 2) & (max <= 5)) {map@data@max <- ceiling(map@data@max); legend_breaks <- seq(0, ceiling(map@data@max), 1)}
  if ((max > 5) & (max <= 10)) {map@data@max <- ceiling(map@data@max); legend_breaks <- seq(0, ceiling(map@data@max), 2)}
  if ((max > 10) & (max <= 15)) {map@data@max <- ceiling(map@data@max); legend_breaks <- seq(0, ceiling(map@data@max), 3)}
  if ((max > 15)) {legend_breaks <- seq(0, max, 5)}
  
  # Map the ring richness map
  map_unit_Mollweide(x = map,
                     unit_pts = ring_occ_Mollweide,
                     color_palette = pal_bl_red_Mannion,
                     
                     main_title = paste0("Mimicry ring ", ring_title, "\n(", list.rings$N.species[i], " species)"),
                     
                     #### Try to put title for species in italic
                     # substitute(paste(italic('p value'), " = 0.01")))
                     # paste0('bold(',list.subtribes$subtribe_symbol[subtribe_index],')'), parse = T,
                     # expression(italic(.(unit_spec[1])))
                     # main_title = parse(text = paste0('italic(',unit_spec[1], ')')),
                     # main_title = bquote(atop(italic(.(unit_spec[1])), ~ "(" ~ .(unit_spec[2]) ~ "patterns)")),
                     
                     occ_cex = 0.35,
                     
                     scale_bar_position = c(-2400, -4000),
                     legend_breaks = legend_breaks,
                     legend_title = "Number of\nspecies ",
                     legend_title_x = -3200,
                     legend_title_y = 850)
  
  # Generate the phylo inset
  # Ithomiini_phylo_plot_taxa <- generate_phylo_inset_tiny_for_ring_heatmap(ring = ring)
  Ithomiini_phylo_plot_taxa <- generate_phylo_inset_tiny_for_ring_tippoint(ring = ring)
  # print(Ithomiini_phylo_plot_taxa)
  
  # Ithomiini_phylo_plot_taxa <- ggplotify::as.grob(Ithomiini_phylo_plot_taxa)
  
  library(grid)
  # Create an appropriate viewport to add the ggplot as inset.
  vp.BottomRight <- viewport(height = unit(.22, "npc"), width = unit(0.22, "npc"), 
                             just = c("left","top"), 
                             y = 0.868, x = 0.705,
                             clip = "on")
  
  # plot the ggplot using the print command
  print(Ithomiini_phylo_plot_taxa, vp = vp.BottomRight)
  # print(Ithomiini_phylo_plot, vp = vp.BottomRight)
  
  if (i %% 5 == 0)
  {
    cat(paste0("\n", Sys.time()," ------ Plotting ", ring_title, " = Mimicry ring N°",i,"/",length(list.rings$Mimicry)," ------\n")) 
  }
}
dev.off()

# Copy in Supplementary folder
file.copy(from = paste0("./maps/By_ring/Richness_maps/0_Final_richness_map_with_phylo_all_rings.pdf"), 
          to = paste0("./supplementaries/Maps/Mimicry_ring_maps/Richness_maps/0_Richness_maps_all_rings.pdf"),
          overwrite = T)

