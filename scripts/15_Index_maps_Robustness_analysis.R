
##### Script 14b: Robustness anlayses of index maps 

### Compute again index maps but with methodological changes #####

# Only species from phylogeny
# Only OMU with niche model
# Jaccard.95 = with buffer set to quantile 95%
# TSS.80 = with TSS instead of Jaccard for submodel selection

### Version with 6 Indices + 3 Zoom on Andes

# A/ Species richness
# B/ Mean Species geographic rarity
# C/ Faith's Phylogenetic Diversity

# D/ Mimicry richness
# E/ Mean mimicry geographic rarity
# F/ Weighted mean ring size

###


# Inputs:
   # Stack of species 'probability' of presence from script 12
   # Stack of mimicry ring 'probability' of presence from script 13
   # Stack of mimicry ring richness of presence from script 13
   # Phylogeny

# Outputs:
   # Alternative stacks for OMU/species/ring with other methodological choices
   # Alternative index maps with other methodological choices

###



# Clean environment
rm(list = ls())

### 0/ Load stuff ####

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
Andes_grid_Mollweide_out <- readRDS(file = "./input_data/Map_stuff/Andes_grid_Mollweide_out.rds")

load(file = "./input_data/Map_stuff/country_borders.RData")
country_borders <- as(country_borders, "Spatial")
plot(country_borders)

# Load Summary table for OMU/unit and for Species
load(file = paste0("./input_data/list.sp.RData"))
load(file = paste0("./input_data/list.models.RData"))

### Load directly the phylogeny with the good labels
phylo.Ithomiini <- readRDS(file = "./input_data/Phylogenies/Final_phylogeny.rds")

# Extract only the 339 species included in the phylogeny from list.sp and the sp.stack
list.sp$In_phylo <- list.sp$Sp_full %in% phylo.Ithomiini$tip.label
list.sp_phylo <- list.sp[list.sp$In_phylo,]

save(list.sp_phylo, file = paste0("./input_data/list.sp_phylo.RData"))
saveRDS(list.sp_phylo, file = paste0("./input_data/list.sp_phylo.rds"))

# Extract the 719 OMU included in the phylogeny from list.models
list.models$In_phylo <- list.models$Sp_full %in% phylo.Ithomiini$tip.label
list.models_phylo <- list.models[list.models$In_phylo,]

save(list.models_phylo, file = paste0("./input_data/list.models_phylo.RData"))
saveRDS(list.models_phylo, file = paste0("./input_data/list.models_phylo.rds"))

# Load the sp.list with only species from the phylogeny
load(file = paste0("./input_data/list.sp_phylo.RData"))
load(file = paste0("./input_data/list.models_phylo.RData"))


### 1/ Generation of species proba stack  ####

### Load directly the complete stack of species proba
All_sp_proba_stack_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.80_phylo.rds"))

### 1.1/ For TSS.80 ####
All_sp_proba_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80.rds"))

### 1.2/ For Jaccard.95 ####
All_sp_proba_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.95.rds"))

### 1.3/ For phylo only ####

# Keep only species in the phylogeny
All_sp_proba_stack_Jaccard.80_phylo <- All_sp_proba_stack_Jaccard.80[[phylo.Ithomiini$tip.label]]

# Save new stacks
save(All_sp_proba_stack_Jaccard.80_phylo, file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_Jaccard.80_phylo.RData"), version = "2")
saveRDS(All_sp_proba_stack_Jaccard.80_phylo, file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_Jaccard.80_phylo.rds"), version = "2")

### 1.4/ For no_binary. Need to aggregate at species level again, only with modeled OMUs ####

# Load Summary table for OMU/unit and for Species
load(file = paste0("./input_data/list.sp.RData"))
load(file = paste0("./input_data/list.models.RData"))

# Make a list of modeled and non-modeled units (occ.unit = with only occurrences points)
modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]

length(unique(modeled_OMU$Sp_full)) # 325 species with at least one modeled OMU

# Function to compute probability of presence with multiple OMUs
aggreg_prob = function(x, na.rm) { 
  y <- 1-prod(1-x) # Probability of presence of species = probability of presence of at least one OMU = opposite of probability of absence of all OMU
  return(y) # Output
}

# Load continent mask as template
load(file = paste0("./input_data/Env_data/continent_mask_15.RData"))

### Loop by species
All_sp_proba_stack_no_binary <- stack(continent_mask)
k <- 1
for (i in 1:nrow(list.sp)) 
{
  # i <- 23
  
  # Load a random layer to initiate the stacks (to remove later)
  sp.stack_no_binary <- stack(continent_mask) 
  
  # Load Species name
  sp <- as.character(list.sp$Sp_full[i]) 
  
  OMU_list <- modeled_OMU[modeled_OMU$Sp_full == sp, ]
  if (nrow(OMU_list) > 0) {      # If at least one OMU was modeled
    
    for (j in 1:nrow(OMU_list))  # For each OMU
    { 
      unit <-  as.character(OMU_list$Tag.model[j]) # Load the unit name                  
      
      # Jaccard.80: Load the continuous maps and stack them all
      Ensemble_Jaccard_median_cropped_80 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_Jaccard_median_cropped_80.rds"))
      sp.stack_no_binary <- addLayer(sp.stack_no_binary, Ensemble_Jaccard_median_cropped_80)
    
    }
    
    # Drop the useless first layer used to initiate the stack
    sp.stack_no_binary <- dropLayer(sp.stack_no_binary, i = 1)
    
    ### Compute probability map per species ####
    
    sp_proba_no_binary <- calc(sp.stack_no_binary, fun = aggreg_prob)
    
    # plot(sp_proba_no_binary, main = paste0(sp,"\nJaccard.80"))
    
    # Make sure to save the data, not just the link to the temp file!
    sp_proba_no_binary <- readAll(sp_proba_no_binary)
    
    save(sp_proba_no_binary, file = paste0("./outputs/By_species/",sp,"/sp_proba_no_binary_",sp,".RData"), version = "2")
    saveRDS(sp_proba_no_binary, file = paste0("./outputs/By_species/",sp,"/sp_proba_no_binary_",sp,".rds"), version = "2")
    
    # Stack final species probability maps
    All_sp_proba_stack_no_binary <- addLayer(All_sp_proba_stack_no_binary, sp_proba_no_binary)
    names(All_sp_proba_stack_no_binary)[k + 1] <- sp
    
    k <- k + 1 
  }
  
  print(i)
  
}

# Drop the useless first layer used to initiate the stack
All_sp_proba_stack_no_binary <- dropLayer(All_sp_proba_stack_no_binary, i = 1)

nlayers(All_sp_proba_stack_no_binary) # 325 species with at least one modeled OMUs
# plot(All_sp_proba_stack_no_binary)

# Make sure to save the data, not just the link to the temp file!
All_sp_proba_stack_no_binary <- readAll(All_sp_proba_stack_no_binary)

save(All_sp_proba_stack_no_binary, file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_no_binary.RData"), version = "2")
saveRDS(All_sp_proba_stack_no_binary, file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_no_binary.rds"), version = "2")


##### 2/ Generation of Mimicry ring stacks  ####

### 2.1/ For TSS.80 ####
All_ring_proba_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_TSS.80.RData"))

### 2.2/ For Jaccard.95 ####
All_ring_proba_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_Jaccard.95.RData"))

### 2.3/ For OMU of species found in the phylogeny (only_phylo)

# Need to rebuild mimicry ring maps only with OMU found in the phylogeny

# Load Summary table for OMU/unit and for Species
load(file = paste0("./input_data/list.sp_phylo.RData"))
load(file = paste0("./input_data/list.models_phylo.RData"))

# Mimicry list
mimicry.list <- as.character(unique(list.models_phylo$Mimicry.model)) # 44 Mimicry rings

# Make a list of modeled and non-modeled units (occ.unit = with only occurrences points)
modeled_OMU <- list.models_phylo[!is.na(list.models_phylo$Model_ID), ]
rasterized_OMU <- list.models_phylo[is.na(list.models_phylo$Model_ID), ]

# Function to compute probability of presence with multiple OMUs
aggreg_prob = function(x, na.rm) { 
  y <- 1-prod(1-x) # Probability of presence of species = probability of presence of at least one OMU = opposite of probability of absence of all OMU
  return(y) # Output
}

# Load continent mask as template
load(file = paste0("./input_data/Env_data/continent_mask_15.RData"))

All_ring_proba_stack_phylo <- All_ring_rich_stack_phylo <- stack(continent_mask) # 1e temp layer used to initiate stack, to remove afterwards
### Loop by mimicry ring
for (i in 1:length(mimicry.list)) 
{
  # i <- 10
  
  # Load a random layer to initiate the stacks (to remove later)
  ring_stack_phylo <- stack(continent_mask) 
  
  # Load Mimicry ring name
  ring <- as.character(mimicry.list[i])
  
  ### 2.3.1/ Make stacks of OMU maps ####
  
  OMU_list <- modeled_OMU[modeled_OMU$Mimicry.model == ring, ]
  if (nrow(OMU_list) > 0) {      # If at least one OMU was modeled
    
    for (j in 1:nrow(OMU_list))  # For each OMU
    { 
      unit <-  as.character(OMU_list$Tag.model[j]) # Load the unit name                  
      
      # Jaccard.80: Load the continuous maps and stack them all
      Ensemble_Jaccard_median_cropped_80 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_Jaccard_median_cropped_80.rds")) 
      ring_stack_phylo <- addLayer(ring_stack_phylo, Ensemble_Jaccard_median_cropped_80)   
      
    }
  }
  
  Binaries <- rasterized_OMU[rasterized_OMU$Mimicry.model == ring, ]
  if (nrow(Binaries) > 0) {  # If at least one OMU was rasterized
    
    for (j in 1:nrow(Binaries)) {  # For each  OMU
      
      unit <-  as.character(Binaries$Tag.model[j]) # Load the unit name 
      
      rasterized_map <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/rasterized_map_",unit,".rds"))
      
      # Jaccard.80: Load the binary maps and stack them all
      ring_stack_phylo <- addLayer(ring_stack_phylo, rasterized_map)
      
    }
  }
  
  # Drop the useless first layer used to initiate the stack
  ring_stack_phylo <- dropLayer(ring_stack_phylo, i = 1)

  ### 2.3.2/ Compute probability map per ring ####
  
  ring_proba_phylo <- calc(ring_stack_phylo, fun = aggreg_prob)
  # plot(ring.cont_Jaccard.80, main = paste0(ring,"\nJaccard.80"))
  
  # Make sure to save the data, not just the link to the temp file!
  ring_proba_phylo <- readAll(ring_proba_phylo)

  save(ring_proba_phylo, file = paste0("./outputs/Mimicry_rings_proba/Phylo_only/Ring_proba_phylo_",ring,".RData"), version = "2")
  saveRDS(ring_proba_phylo, file = paste0("./outputs/Mimicry_rings_proba/Phylo_only/Ring_proba_phylo_",ring,".rds"), version = "2")

  ### 2.3.3/ Compute richness map per ring (nb of OMU) ####
  
  ring_rich_phylo <- calc(ring_stack_phylo, fun = sum)
  # plot(ring_rich_phylo, main = paste0(ring,"\nJaccard.80"))
  
  # Make sure to save the data, not just the link to the temp file!
  ring_rich_phylo <- readAll(ring_rich_phylo)
  
  save(ring_rich_phylo, file = paste0("./outputs/Mimicry_ring_richness/Phylo_only/Ring_rich_phylo_",ring,".RData"), version = "2")
  saveRDS(ring_rich_phylo, file = paste0("./outputs/Mimicry_ring_richness/Phylo_only/Ring_rich_phylo_",ring,".rds"), version = "2")

  ### 2.3.4/ Build mimicry stacks for proba and richness ####
  
  All_ring_proba_stack_phylo <- addLayer(All_ring_proba_stack_phylo, ring_proba_phylo)
  All_ring_rich_stack_phylo <- addLayer(All_ring_rich_stack_phylo, ring_rich_phylo)
  
  print(i)
  
}

# Drop the useless first layer used to initiate the stack
All_ring_proba_stack_phylo <- dropLayer(All_ring_proba_stack_phylo, i = 1)
All_ring_rich_stack_phylo <- dropLayer(All_ring_rich_stack_phylo, i = 1)

# Name layers with ring names
names(All_ring_proba_stack_phylo) <- names(All_ring_rich_stack_phylo) <- mimicry.list

# nlayers(All_ring_proba_stack_phylo) # 44 rings in the final stack

# # Plot results
# plot(All_ring_proba_stack_phylo)
# plot(All_ring_rich_stack_phylo)

# Save stacks
save(All_ring_proba_stack_phylo, file = paste0("./outputs/Indices_stacks/phylo_only/All_ring_proba_stack_phylo.RData"), version = "2")
saveRDS(All_ring_proba_stack_phylo, file = paste0("./outputs/Indices_stacks/phylo_only/All_ring_proba_stack_phylo.rds"), version = "2")

save(All_ring_rich_stack_phylo, file = paste0("./outputs/Indices_stacks/phylo_only/All_ring_rich_stack_phylo.RData"), version = "2")
saveRDS(All_ring_rich_stack_phylo, file = paste0("./outputs/Indices_stacks/phylo_only/All_ring_rich_stack_phylo.rds"), version = "2")


### 2.4/ For OMU of species found in the phylogeny (only_phylo)

# Need to rebuild mimicry ring maps only with OMU found in the phylogeny

# Load Summary table for OMU/unit and for Species
load(file = paste0("./input_data/list.sp.RData"))
load(file = paste0("./input_data/list.models.RData"))

# Mimicry list
mimicry.list <- as.character(unique(list.models$Mimicry.model)) # 44 Mimicry rings

# Make a list of modeled and non-modeled units (occ.unit = with only occurrences points)
modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]

# Function to compute probability of presence with multiple OMUs
aggreg_prob = function(x, na.rm) { 
  y <- 1-prod(1-x) # Probability of presence of species = probability of presence of at least one OMU = opposite of probability of absence of all OMU
  return(y) # Output
}

# Load continent mask as template
load(file = paste0("./input_data/Env_data/continent_mask_15.RData"))

All_ring_proba_stack_no_binary <- All_ring_rich_stack_no_binary <- stack(continent_mask) # 1e temp layer used to initiate stack, to remove afterwards
k <- 1
### Loop by mimicry ring
for (i in 1:length(mimicry.list)) 
{
  # i <- 10
  
  # Load a random layer to initiate the stacks (to remove later)
  ring_stack_no_binary <- stack(continent_mask) 
  
  # Load Mimicry ring name
  ring <- as.character(mimicry.list[i])
  
  ### 2.4.1/ Make stacks of OMU maps ####
  
  OMU_list <- modeled_OMU[modeled_OMU$Mimicry.model == ring, ]
  if (nrow(OMU_list) > 0) {      # If at least one OMU was modeled
    
    for (j in 1:nrow(OMU_list))  # For each OMU
    { 
      unit <-  as.character(OMU_list$Tag.model[j]) # Load the unit name                  
      
      # Jaccard.80: Load the continuous maps and stack them all
      Ensemble_Jaccard_median_cropped_80 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_Jaccard_median_cropped_80.rds")) 
      ring_stack_no_binary <- addLayer(ring_stack_no_binary, Ensemble_Jaccard_median_cropped_80)   
      
    }
    
    # Drop the useless first layer used to initiate the stack
    ring_stack_no_binary <- dropLayer(ring_stack_no_binary, i = 1)
    
    ### 2.4.2/ Compute probability map per ring ####
    
    ring_proba_no_binary <- calc(ring_stack_no_binary, fun = aggreg_prob)
    # plot(ring_proba_no_binary, main = paste0(ring,"\nNo binaries"))
    
    # Make sure to save the data, not just the link to the temp file!
    ring_proba_no_binary <- readAll(ring_proba_no_binary)
    
    save(ring_proba_no_binary, file = paste0("./outputs/Mimicry_rings_proba/no_binary/Ring_proba_no_binary_",ring,".RData"), version = "2")
    saveRDS(ring_proba_no_binary, file = paste0("./outputs/Mimicry_rings_proba/no_binary/Ring_proba_no_binary_",ring,".rds"), version = "2")
    
    ### 2.4.3/ Compute richness map per ring (nb of OMU) ####
    
    ring_rich_no_binary <- calc(ring_stack_no_binary, fun = sum)
    # plot(ring_rich_no_binary, main = paste0(ring,"\nJaccard.80"))
    
    # Make sure to save the data, not just the link to the temp file!
    ring_rich_no_binary <- readAll(ring_rich_no_binary)
    
    save(ring_rich_no_binary, file = paste0("./outputs/Mimicry_ring_richness/no_binary/Ring_rich_no_binary_",ring,".RData"), version = "2")
    saveRDS(ring_rich_no_binary, file = paste0("./outputs/Mimicry_ring_richness/no_binary/Ring_rich_no_binary_",ring,".rds"), version = "2")
    
    ### 2.4.4/ Build mimicry stacks for proba and richness ####
    
    All_ring_proba_stack_no_binary <- addLayer(All_ring_proba_stack_no_binary, ring_proba_no_binary)
    All_ring_rich_stack_no_binary <- addLayer(All_ring_rich_stack_no_binary, ring_rich_no_binary)
    
    names(All_ring_proba_stack_no_binary)[k + 1] <- names(All_ring_rich_stack_no_binary)[k + 1] <- ring
    
    k <- k + 1
  }
  
  print(i)
  
}

# Drop the useless first layer used to initiate the stack
All_ring_proba_stack_no_binary <- dropLayer(All_ring_proba_stack_no_binary, i = 1)
All_ring_rich_stack_no_binary <- dropLayer(All_ring_rich_stack_no_binary, i = 1)

# nlayers(All_ring_proba_stack_no_binary) # 42 rings in the final stack that have at least one modeled OMU

# # Plot results
# plot(All_ring_proba_stack_no_binary)
# plot(All_ring_rich_stack_no_binary)

# Save stacks
save(All_ring_proba_stack_no_binary, file = paste0("./outputs/Indices_stacks/no_binary/All_ring_proba_stack_no_binary.RData"), version = "2")
saveRDS(All_ring_proba_stack_no_binary, file = paste0("./outputs/Indices_stacks/no_binary/All_ring_proba_stack_no_binary.rds"), version = "2")

save(All_ring_rich_stack_no_binary, file = paste0("./outputs/Indices_stacks/no_binary/All_ring_rich_stack_no_binary.RData"), version = "2")
saveRDS(All_ring_rich_stack_no_binary, file = paste0("./outputs/Indices_stacks/no_binary/All_ring_rich_stack_no_binary.rds"), version = "2")




##### 3/ Compute new index maps ####

### 3.0/ Choose subset and generate indices maps

# Choose your alternative set for robustness analyses
set <- "only_phylo"
set <- "no_binary"
set <- "TSS.80"
set <- "Jaccard.95"

# Load new stacks
if (set == "only_phylo")
{
  All_sp_proba_stack_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_Jaccard.80_phylo.rds"))
  All_ring_proba_stack_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_ring_proba_stack_phylo.rds"))
  All_ring_rich_stack_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_ring_rich_stack_phylo.rds"))
}

if (set == "no_binary")
{
  All_sp_proba_stack_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_Jaccard.80_no_binary.rds"))
  All_ring_proba_stack_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_ring_proba_stack_no_binary.rds"))
  All_ring_rich_stack_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_ring_rich_stack_no_binary.rds"))
}

if (set == "TSS.80")
{
  All_sp_proba_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80.rds"))
  All_ring_proba_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_TSS.80.RData"))
  All_ring_rich_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_TSS.80.RData"))
}

if (set == "Jaccard.95")
{
  All_sp_proba_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.95.rds"))
  All_ring_proba_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_Jaccard.95.RData"))
  All_ring_rich_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_Jaccard.95.RData"))
}

### 3.1/ Species richness ######

# Function
compute_sp.richness <- function(proba_stack)
{
  sp.richness <- readAll(calc(proba_stack, fun = sum))
  return(sp.richness)
}

# Computation
sp.richness_phylo <- compute_sp.richness(All_sp_proba_stack_phylo)
sp.richness_no_binary <- compute_sp.richness(All_sp_proba_stack_no_binary)
sp.richness_TSS.80 <- compute_sp.richness(All_sp_proba_stack_TSS.80)
sp.richness_Jaccard.95 <- compute_sp.richness(All_sp_proba_stack_Jaccard.95)

# Save
save(sp.richness_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.richness_phylo.RData"), version = "2")
saveRDS(sp.richness_phylo, file = paste0("./outputs/Indices_maps/phylo_only/sp.richness_phylo.rds"), version = "2")

save(sp.richness_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.richness_no_binary.RData"), version = "2")
saveRDS(sp.richness_no_binary, file = paste0("./outputs/Indices_maps/no_binary/sp.richness_no_binary.rds"), version = "2")

save(sp.richness_TSS.80, file = paste0("./outputs/Indices_maps/TSS.80/sp.richness_TSS.80.RData"), version = "2")
saveRDS(sp.richness_TSS.80, file = paste0("./outputs/Indices_maps/TSS.80/sp.richness_TSS.80.rds"), version = "2")

save(sp.richness_Jaccard.95, file = paste0("./outputs/Indices_maps/Jaccard.95/sp.richness_Jaccard.95.RData"), version = "2")
saveRDS(sp.richness_Jaccard.95, file = paste0("./outputs/Indices_maps/Jaccard.95/sp.richness_Jaccard.95.rds"), version = "2")


### 3.2/ Mean Species geographic rarity ####

# Contrasting raster function
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

# Compute species mean geography function
compute_sp.mean.rarity <- function(proba_stack, sp.richness)
{
  # Compute weights based on geographic ranges
  sp.range_size <- NA
  for (i in 1:nlayers(proba_stack)) {
    # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
    # Compute the estimated number of pixels occupied by the species as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
    sp.range_size[i] <- sum(proba_stack[[i]]@data@values, na.rm = T)*774.51/1000
  }
  
  library(Rarity)
  # ?rWeights
  
  # Build a matrix of species assemblage with continuous outputs from SDM (remove empty communities)
  sp_proba_brick <- readAll(brick(proba_stack))
  sp_assemblage <- t(tidyr::drop_na(as.data.frame(sp_proba_brick@data@values)))
  sp_assemblage <- sp_assemblage[, colSums(sp_assemblage) >= 1]
  
  # Rarity weights based on an inverse exponential function with inflexion point calibrated as the rarity threshold as such that communities host 25% of rare species in average
  Leroy_weights_df <- rWeights(occData = sp.range_size, wMethods = "W", rCutoff = "Leroy",
                               normalised = T, rounding = 5, assemblages = sp_assemblage)
  sp.rarity_Leroy_weights <- Leroy_weights_df$W
  
  # hist(sp.rarity_Leroy_weights)
  # plot(sp.rarity_Leroy_weights ~ sp.range_size)
  
  # Apply weights
  sp_rarity_weighted_stack <- proba_stack # Generate the stack to fill
  for (i in 1:nlayers(sp_rarity_weighted_stack)){
    # Multiply species probability of presence by species rarity indices
    sp_rarity_weighted_stack@layers[[i]]@data@values <- proba_stack@layers[[i]]@data@values * sp.rarity_Leroy_weights[i]
  }
  
  # Sum = Richness weighted by rarity based on Range-size
  sp.rarity <- readAll(calc(sp_rarity_weighted_stack, fun = sum)) 
  # Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each species)
  sp.mean.rarity <- sp.rarity/sp.richness
  
  # Add 0 values for empty continental pixels (transformed into min values as 0.05)
  continent_mask <- readRDS(file = paste0("./input_data/Map_stuff/continent_mask.rds"))
  sp.mean.rarity <- sp.mean.rarity + 0.005
  temp <- continent_mask
  temp[!is.na(sp.mean.rarity[])] <- sp.mean.rarity[!is.na(sp.mean.rarity[])]
  sp.mean.rarity <- temp
  
  # Need to contrast values above 0.5
  # hist(sp.mean.rarity_phylo@data@values)
  
  sp.mean.rarity_contrasted <- contrasting_raster(x = sp.mean.rarity, zmin = 0, zmax = 0.5)
  
  # Quick plot
  plot(sp.mean.rarity, col = pal_bl_red_Mannion, main = paste0("Species mean geographic rarity\n", set))
  plot(sp.mean.rarity_contrasted, col = pal_bl_red_Mannion, main = paste0("Species mean geographic rarity\n", set, " (contrasted)"))
  
  return(list(sp.range_size, sp.rarity_Leroy_weights, sp.mean.rarity, sp.mean.rarity_contrasted))
}

# Computation
sp.mean.rarity_list_phylo <- compute_sp.mean.rarity(proba_stack = All_sp_proba_stack_phylo,
                                                    sp.richness = sp.richness_phylo)
sp.mean.rarity_list_no_binary <- compute_sp.mean.rarity(proba_stack = All_sp_proba_stack_no_binary,
                                                        sp.richness = sp.richness_no_binary)
sp.mean.rarity_list_TSS.80 <- compute_sp.mean.rarity(proba_stack = All_sp_proba_stack_TSS.80,
                                                     sp.richness = sp.richness_TSS.80)
sp.mean.rarity_list_Jaccard.95 <- compute_sp.mean.rarity(proba_stack = All_sp_proba_stack_Jaccard.95,
                                                     sp.richness = sp.richness_Jaccard.95)

# Retrieve outputs
sp.range_size_phylo <- sp.mean.rarity_list_phylo[[1]]
sp.rarity_Leroy_weights_phylo <- sp.mean.rarity_list_phylo[[2]]
sp.mean.rarity_phylo <- sp.mean.rarity_list_phylo[[3]]
sp.mean.rarity_phylo_contrasted <- sp.mean.rarity_list_phylo[[4]]

sp.range_size_no_binary <- sp.mean.rarity_list_no_binary[[1]]
sp.rarity_Leroy_weights_no_binary <- sp.mean.rarity_list_no_binary[[2]]
sp.mean.rarity_no_binary <- sp.mean.rarity_list_no_binary[[3]]
sp.mean.rarity_no_binary_contrasted <- sp.mean.rarity_list_no_binary[[4]]

sp.range_size_TSS.80 <- sp.mean.rarity_list_TSS.80[[1]]
sp.rarity_Leroy_weights_TSS.80 <- sp.mean.rarity_list_TSS.80[[2]]
sp.mean.rarity_TSS.80 <- sp.mean.rarity_list_TSS.80[[3]]
sp.mean.rarity_TSS.80_contrasted <- sp.mean.rarity_list_TSS.80[[4]]

sp.range_size_Jaccard.95 <- sp.mean.rarity_list_Jaccard.95[[1]]
sp.rarity_Leroy_weights_Jaccard.95 <- sp.mean.rarity_list_Jaccard.95[[2]]
sp.mean.rarity_Jaccard.95 <- sp.mean.rarity_list_Jaccard.95[[3]]
sp.mean.rarity_Jaccard.95_contrasted <- sp.mean.rarity_list_Jaccard.95[[4]]

# Save species ranges 
save(sp.range_size_phylo, file = "./outputs/Indices_maps/phylo_only/sp.range_size_phylo.RData", version = "2")
saveRDS(sp.range_size_phylo, file = "./outputs/Indices_maps/phylo_only/sp.range_size_phylo.rds", version = "2")

save(sp.range_size_no_binary, file = "./outputs/Indices_maps/no_binary/sp.range_size_no_binary.RData", version = "2")
saveRDS(sp.range_size_no_binary, file = "./outputs/Indices_maps/no_binary/sp.range_size_no_binary.rds", version = "2")

save(sp.range_size_TSS.80, file = "./outputs/Indices_maps/TSS.80/sp.range_size_TSS.80.RData", version = "2")
saveRDS(sp.range_size_TSS.80, file = "./outputs/Indices_maps/TSS.80/sp.range_size_TSS.80.rds", version = "2")

save(sp.range_size_Jaccard.95, file = "./outputs/Indices_maps/Jaccard.95/sp.range_size_Jaccard.95.RData", version = "2")
saveRDS(sp.range_size_Jaccard.95, file = "./outputs/Indices_maps/Jaccard.95/sp.range_size_Jaccard.95.rds", version = "2")

# Save rarity weights
save(sp.rarity_Leroy_weights_phylo, file = "./outputs/Indices_maps/phylo_only/sp.rarity_Leroy_weights_phylo.RData", version = "2")
saveRDS(sp.rarity_Leroy_weights_phylo, file = "./outputs/Indices_maps/phylo_only/sp.rarity_Leroy_weights_phylo.rds", version = "2")

save(sp.rarity_Leroy_weights_no_binary, file = "./outputs/Indices_maps/no_binary/sp.rarity_Leroy_weights_no_binary.RData", version = "2")
saveRDS(sp.rarity_Leroy_weights_no_binary, file = "./outputs/Indices_maps/no_binary/sp.rarity_Leroy_weights_no_binary.rds", version = "2")

save(sp.rarity_Leroy_weights_TSS.80, file = "./outputs/Indices_maps/TSS.80/sp.rarity_Leroy_weights_TSS.80.RData", version = "2")
saveRDS(sp.rarity_Leroy_weights_TSS.80, file = "./outputs/Indices_maps/TSS.80/sp.rarity_Leroy_weights_TSS.80.rds", version = "2")

save(sp.rarity_Leroy_weights_Jaccard.95, file = "./outputs/Indices_maps/Jaccard.95/sp.rarity_Leroy_weights_Jaccard.95.RData", version = "2")
saveRDS(sp.rarity_Leroy_weights_Jaccard.95, file = "./outputs/Indices_maps/Jaccard.95/sp.rarity_Leroy_weights_Jaccard.95.rds", version = "2")

# Save mean rarity maps
save(sp.mean.rarity_phylo, file = "./outputs/Indices_maps/phylo_only/sp.mean.rarity_phylo.RData", version = "2")
saveRDS(sp.mean.rarity_phylo, file = "./outputs/Indices_maps/phylo_only/sp.mean.rarity_phylo.rds", version = "2")
save(sp.mean.rarity_phylo_contrasted, file = "./outputs/Indices_maps/phylo_only/sp.mean.rarity_phylo_contrasted.RData", version = "2")
saveRDS(sp.mean.rarity_phylo_contrasted, file = "./outputs/Indices_maps/phylo_only/sp.mean.rarity_phylo_contrasted.rds", version = "2")

save(sp.mean.rarity_no_binary, file = "./outputs/Indices_maps/no_binary/sp.mean.rarity_no_binary.RData", version = "2")
saveRDS(sp.mean.rarity_no_binary, file = "./outputs/Indices_maps/no_binary/sp.mean.rarity_no_binary.rds", version = "2")
save(sp.mean.rarity_no_binary_contrasted, file = "./outputs/Indices_maps/no_binary/sp.mean.rarity_no_binary_contrasted.RData", version = "2")
saveRDS(sp.mean.rarity_no_binary_contrasted, file = "./outputs/Indices_maps/no_binary/sp.mean.rarity_no_binary_contrasted.rds", version = "2")

save(sp.mean.rarity_TSS.80, file = "./outputs/Indices_maps/TSS.80/sp.mean.rarity_TSS.80.RData", version = "2")
saveRDS(sp.mean.rarity_TSS.80, file = "./outputs/Indices_maps/TSS.80/sp.mean.rarity_TSS.80.rds", version = "2")
save(sp.mean.rarity_TSS.80_contrasted, file = "./outputs/Indices_maps/TSS.80/sp.mean.rarity_TSS.80_contrasted.RData", version = "2")
saveRDS(sp.mean.rarity_TSS.80_contrasted, file = "./outputs/Indices_maps/TSS.80/sp.mean.rarity_TSS.80_contrasted.rds", version = "2")

save(sp.mean.rarity_Jaccard.95, file = "./outputs/Indices_maps/Jaccard.95/sp.mean.rarity_Jaccard.95.RData", version = "2")
saveRDS(sp.mean.rarity_Jaccard.95, file = "./outputs/Indices_maps/Jaccard.95/sp.mean.rarity_Jaccard.95.rds", version = "2")
save(sp.mean.rarity_Jaccard.95_contrasted, file = "./outputs/Indices_maps/Jaccard.95/sp.mean.rarity_Jaccard.95_contrasted.RData", version = "2")
saveRDS(sp.mean.rarity_Jaccard.95_contrasted, file = "./outputs/Indices_maps/Jaccard.95/sp.mean.rarity_Jaccard.95_contrasted.rds", version = "2")


### 3.3/ Phylogenetic diversity ####

library(ape)
library(picante)
library(geiger)
library(treespace)

### Create function to aggregate probabilities to higher hierachical level (aggregate pixel, or go up on a phylogenetic tree)
aggreg_prob = function(x, na.rm) 
{ 
  y <- 1-prod(1-x) 
  return(y) # Output
}


### 3.3.1/ Prepare stack with only species in the phylogeny

### Load directly the phylogeny with the good labels
phylo.Ithomiini <- readRDS(file = "./input_data/Phylogenies/Final_phylogeny.rds")

### Load directly the stacks of sp probas
All_sp_proba_stack_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_Jaccard.80_phylo.rds"))
All_sp_proba_stack_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_no_binary.rds"))
All_sp_proba_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80.rds"))
All_sp_proba_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.95.rds"))

### Need to extract only the species in the phylogeny

# Load Summary table for OMU/unit and for Species
load(file = paste0("./input_data/list.sp.RData"))

# Extract only the 339 species included in the phylogeny from list.sp and the sp.stack
list.sp_phylo <- list.sp[list.sp$Sp_full %in% phylo.Ithomiini$tip.label,]

# Extract only the 339 species included in the phylogeny from the stacks of sp probas
All_sp_proba_stack_phylo <- All_sp_proba_stack_phylo[[phylo.Ithomiini$tip.label]] # Useless since it already contains only species from phylogeny
All_sp_proba_stack_no_binary_phylo <- All_sp_proba_stack_no_binary[[phylo.Ithomiini$tip.label]]
All_sp_proba_stack_TSS.80_phylo <- All_sp_proba_stack_TSS.80[[phylo.Ithomiini$tip.label]]
All_sp_proba_stack_Jaccard.95_phylo <- All_sp_proba_stack_Jaccard.95[[phylo.Ithomiini$tip.label]]

nlayers(All_sp_proba_stack_phylo) # 339 species in the phylogeny
nlayers(All_sp_proba_stack_no_binary_phylo) # 300 species in the phylogeny with at least one modeled OMU
nlayers(All_sp_proba_stack_TSS.80_phylo) # 339 species in the phylogeny
nlayers(All_sp_proba_stack_Jaccard.95_phylo) # 339 species in the phylogeny

# Save stacks
save(All_sp_proba_stack_phylo, file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_phylo.RData"), version = "2")
saveRDS(All_sp_proba_stack_phylo, file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_phylo.rds"), version = "2")
save(All_sp_proba_stack_no_binary_phylo, file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_no_binary_phylo.RData"), version = "2")
saveRDS(All_sp_proba_stack_no_binary_phylo, file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_no_binary_phylo.rds"), version = "2")
save(All_sp_proba_stack_TSS.80_phylo, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80_phylo.RData"), version = "2")
saveRDS(All_sp_proba_stack_TSS.80_phylo, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80_phylo.rds"), version = "2")
save(All_sp_proba_stack_Jaccard.95_phylo, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.95_phylo.RData"), version = "2")
saveRDS(All_sp_proba_stack_Jaccard.95_phylo, file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.95_phylo.rds"), version = "2")

# 3.3.2/ Function to compute phylogenetic diversity

# For each community, compute probability of presence of the edge as the probability for at least one descending species to be in the community

compute_phylogenetic_diversity <- function(phylo, proba_stack, set)
{
  # Generate matrix to store info on edges/branches
  branches = matrix(NA, nrow(phylo$edge), ncol = 4) 
  names(branches) <- c("Starting nod", "Ending nod", "Length", "Proba_presence")
  branches[, 1:2] = phylo$edge # Retrieve starting and ending node 
  branches[, 3] = round(phylo$edge.length, 4) # Retrieve edge length
  
  sp.proba.mat <- readAll(brick(proba_stack))@data@values
  
  # Make a loop per community
  all.PD <- NA
  for (k in 1:nrow(sp.proba.mat)) {
    sp.proba.com <- sp.proba.mat[k,] # Extract the proba for the k community
    
    # Compute PD only if no NA is present in the community
    if (any(is.na(sp.proba.com))) {
      
      all.PD[k] <- NA # if NA present, PD = NA
      
    } else {
      
      # Compute probability of presence of each edge/branch (i.e., a species descending from this branch) in this community
      for (i in 1:nrow(branches)) {
        leaves.node = tips(phylo, branches[i, 2]) # Retrieve the set of species descending from branch i
        index <- which(colnames(sp.proba.mat) %in% leaves.node) # Find position of the species in the stack/matrix
        prob_edge <- aggreg_prob(sp.proba.com[index]) # Compute proba of presence of the edge
        branches[i, 4] <- prob_edge # Store info
      }
      PD <- round(sum(branches[,3]*branches[,4]),4) # Compute PD as the weighted sum of the edge length (weighted by the probabilty of presence of this edge in the community)
      all.PD[k]  <- PD
    }
    
    # Show k every 100 iterations and save a backup
    if (k %% 1000 == 0) {
      cat(paste0(Sys.time(), " - ", k," on ",nrow(sp.proba.mat),"\n"))
      save(all.PD, file = "./outputs/Indices_Maps/PD_backup.RData", version = "2")
    }
  }
  
  # Write PD values in a raster
  PD.raster <- proba_stack[[1]]
  PD.raster@data@values <-  all.PD
  
  # Repair issue with max value
  PD.raster@data@max <- max(PD.raster[], na.rm = T)
  
  # Quick plot
  plot(PD.raster, col = pal_bl_red_Mannion, main = paste0("Phylogenetic diversity\n", set))
  
  return(PD.raster)
}

# 3.3.3/ Compute indices #

# Load directly the sp proba stack with only species in the phylogeny
All_sp_proba_stack_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_sp_proba_stack_phylo.rds"))
All_sp_proba_stack_no_binary_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_sp_proba_stack_no_binary_phylo.rds"))
All_sp_proba_stack_TSS.80_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_TSS.80_phylo.rds"))
All_sp_proba_stack_Jaccard.95_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/All_sp_proba_stack_Jaccard.95_phylo.rds"))

# Load directly the phylogeny with the good labels
phylo.Ithomiini <- readRDS(file = "./input_data/Phylogenies/Final_phylogeny.rds")

PD_phylo <- compute_phylogenetic_diversity(phylo = phylo.Ithomiini,
                                           proba_stack = All_sp_proba_stack_phylo,
                                           set = "only_phylo")

PD_no_binary <- compute_phylogenetic_diversity(phylo = phylo.Ithomiini,
                                               proba_stack = All_sp_proba_stack_no_binary_phylo,
                                               set = "no_binary")

PD_TSS.80 <- compute_phylogenetic_diversity(phylo = phylo.Ithomiini,
                                            proba_stack = All_sp_proba_stack_TSS.80_phylo,
                                            set = "TSS.80")

PD_Jaccard.95 <- compute_phylogenetic_diversity(phylo = phylo.Ithomiini,
                                                proba_stack = All_sp_proba_stack_Jaccard.95_phylo,
                                                set = "Jaccard.95")

# Save
save(PD_phylo, file = "./outputs/Indices_Maps/phylo_only/PD_phylo.RData", version = "2")
saveRDS(PD_phylo, file = "./outputs/Indices_Maps/phylo_only/PD_phylo.rds", version = "2")

save(PD_no_binary, file = "./outputs/Indices_Maps/no_binary/PD_no_binary.RData", version = "2")
saveRDS(PD_no_binary, file = "./outputs/Indices_Maps/no_binary/PD_no_binary.rds", version = "2")

save(PD_TSS.80, file = "./outputs/Indices_Maps/TSS.80/PD_TSS.80.RData", version = "2")
saveRDS(PD_TSS.80, file = "./outputs/Indices_Maps/TSS.80/PD_TSS.80.rds", version = "2")

save(PD_Jaccard.95, file = "./outputs/Indices_Maps/Jaccard.95/PD_Jaccard.95.RData", version = "2")
saveRDS(PD_Jaccard.95, file = "./outputs/Indices_Maps/Jaccard.95/PD_Jaccard.95.rds", version = "2")


### 3.4/ Mimicry richness ####

# Function
compute_ring.richness <- function(proba_stack)
{
  ring.richness <- readAll(calc(proba_stack, fun = sum))
  return(ring.richness)
}

# Computation
ring.richness_phylo <- compute_ring.richness(All_ring_proba_stack_phylo)
ring.richness_no_binary <- compute_ring.richness(All_ring_proba_stack_no_binary)
ring.richness_TSS.80 <- compute_ring.richness(All_ring_proba_stack_TSS.80)
ring.richness_Jaccard.95 <- compute_ring.richness(All_ring_proba_stack_Jaccard.95)

# Quick plot
plot(ring.richness_phylo, col = pal_bl_red_Mannion, main = paste0("Mimicry richness\nphylo"))
plot(ring.richness_no_binary, col = pal_bl_red_Mannion, main = paste0("Mimicry richness\nno_binary"))
plot(ring.richness_TSS.80, col = pal_bl_red_Mannion, main = paste0("Mimicry richness\nTSS.80"))
plot(ring.richness_Jaccard.95, col = pal_bl_red_Mannion, main = paste0("Mimicry richness\nJaccard.95"))

# Save
save(ring.richness_phylo, file = paste0("./outputs/Indices_maps/phylo_only/ring.richness_phylo.RData"), version = "2")
saveRDS(ring.richness_phylo, file = paste0("./outputs/Indices_maps/phylo_only/ring.richness_phylo.rds"), version = "2")

save(ring.richness_no_binary, file = paste0("./outputs/Indices_maps/no_binary/ring.richness_no_binary.RData"), version = "2")
saveRDS(ring.richness_no_binary, file = paste0("./outputs/Indices_maps/no_binary/ring.richness_no_binary.rds"), version = "2")

save(ring.richness_TSS.80, file = paste0("./outputs/Indices_maps/TSS.80/ring.richness_TSS.80.RData"), version = "2")
saveRDS(ring.richness_TSS.80, file = paste0("./outputs/Indices_maps/TSS.80/ring.richness_TSS.80.rds"), version = "2")

save(ring.richness_Jaccard.95, file = paste0("./outputs/Indices_maps/Jaccard.95/ring.richness_Jaccard.95.RData"), version = "2")
saveRDS(ring.richness_Jaccard.95, file = paste0("./outputs/Indices_maps/Jaccard.95/ring.richness_Jaccard.95.rds"), version = "2")


### 3.5/ Mean mimicry ring geographic rarity ####

# Contrasting raster function
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

# Compute species mean geography function
compute_ring.mean.rarity <- function(proba_stack, ring.richness)
{
  # Compute weights based on geographic ranges
  ring.range_size <- NA
  for (i in 1:nlayers(proba_stack)) {
    # 1 pixel at 15min resolution ≈ 27.83 x 27.83 km² = 774.51 km²
    # Compute the estimated number of pixels occupied by the mimicr ring as the sum of probabilities and multiply by 774.51/1000 to have values in 10^3 km²
    ring.range_size[i] <- sum(proba_stack[[i]]@data@values, na.rm = T)*774.51/1000
  }
  
  library(Rarity)
  # ?rWeights
  
  # Build a matrix of mimicry ring assemblage with continuous outputs from SDM (remove empty communities)
  ring_proba_brick <- readAll(brick(proba_stack))
  ring_assemblage <- t(tidyr::drop_na(as.data.frame(ring_proba_brick@data@values)))
  ring_assemblage <- ring_assemblage[, colSums(ring_assemblage) >= 1]
  
  # Rarity weights based on an inverse exponential function with inflexion point calibrated as the rarity threshold as such that communities host 25% of rare mimicry rings in average
  Leroy_weights_df <- rWeights(occData = ring.range_size, wMethods = "W", rCutoff = "Leroy",
                               normalised = T, rounding = 5, assemblages = ring_assemblage)
  ring.rarity_Leroy_weights <- Leroy_weights_df$W
  
  # hist(ring.rarity_Leroy_weights)
  # plot(ring.rarity_Leroy_weights ~ ring.range_size)
  
  # Apply weights
  ring_rarity_weighted_stack <- proba_stack # Generate the stack to fill
  for (i in 1:nlayers(ring_rarity_weighted_stack)){
    # Multiply mimicry ring probability of presence by mimicry ring rarity indices
    ring_rarity_weighted_stack@layers[[i]]@data@values <- proba_stack@layers[[i]]@data@values * ring.rarity_Leroy_weights[i]
  }
  
  # Sum = Richness weighted by rarity based on Range-size
  ring.rarity <- readAll(calc(ring_rarity_weighted_stack, fun = sum)) 
  # Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each mimicry ring)
  ring.mean.rarity <- ring.rarity/ring.richness
  
  # Add 0 values for empty continental pixels (transformed into min values as 0.05)
  continent_mask <- readRDS(file = paste0("./input_data/Map_stuff/continent_mask.rds"))
  ring.mean.rarity <- ring.mean.rarity + 0.005
  temp <- continent_mask
  temp[!is.na(ring.mean.rarity[])] <- ring.mean.rarity[!is.na(ring.mean.rarity[])]
  ring.mean.rarity <- temp
  
  # Need to contrast values above 0.5
  # hist(ring.mean.rarity_phylo@data@values)
  
  ring.mean.rarity_contrasted <- contrasting_raster(x = ring.mean.rarity, zmin = 0, zmax = 0.5)
  
  # Quick plot
  plot(ring.mean.rarity, col = pal_bl_red_Mannion, main = paste0("Mimicry ring mean geographic rarity\n", set))
  plot(ring.mean.rarity_contrasted, col = pal_bl_red_Mannion, main = paste0("Mimicry ring mean geographic rarity\n", set, " (contrasted)"))
  
  return(list(ring.range_size, ring.rarity_Leroy_weights, ring.mean.rarity, ring.mean.rarity_contrasted))
}

# Computation
ring.mean.rarity_list_phylo <- compute_ring.mean.rarity(proba_stack = All_ring_proba_stack_phylo,
                                                        ring.richness = ring.richness_phylo)
ring.mean.rarity_list_no_binary <- compute_ring.mean.rarity(proba_stack = All_ring_proba_stack_no_binary,
                                                            ring.richness = ring.richness_no_binary)
ring.mean.rarity_list_TSS.80 <- compute_ring.mean.rarity(proba_stack = All_ring_proba_stack_TSS.80,
                                                         ring.richness = ring.richness_TSS.80)
ring.mean.rarity_list_Jaccard.95 <- compute_ring.mean.rarity(proba_stack = All_ring_proba_stack_Jaccard.95,
                                                             ring.richness = ring.richness_Jaccard.95)

# Retrieve outputs
ring.range_size_phylo <- ring.mean.rarity_list_phylo[[1]]
ring.rarity_Leroy_weights_phylo <- ring.mean.rarity_list_phylo[[2]]
ring.mean.rarity_phylo <- ring.mean.rarity_list_phylo[[3]]
ring.mean.rarity_phylo_contrasted <- ring.mean.rarity_list_phylo[[4]]

ring.range_size_no_binary <- ring.mean.rarity_list_no_binary[[1]]
ring.rarity_Leroy_weights_no_binary <- ring.mean.rarity_list_no_binary[[2]]
ring.mean.rarity_no_binary <- ring.mean.rarity_list_no_binary[[3]]
ring.mean.rarity_no_binary_contrasted <- ring.mean.rarity_list_no_binary[[4]]

ring.range_size_TSS.80 <- ring.mean.rarity_list_TSS.80[[1]]
ring.rarity_Leroy_weights_TSS.80 <- ring.mean.rarity_list_TSS.80[[2]]
ring.mean.rarity_TSS.80 <- ring.mean.rarity_list_TSS.80[[3]]
ring.mean.rarity_TSS.80_contrasted <- ring.mean.rarity_list_TSS.80[[4]]

ring.range_size_Jaccard.95 <- ring.mean.rarity_list_Jaccard.95[[1]]
ring.rarity_Leroy_weights_Jaccard.95 <- ring.mean.rarity_list_Jaccard.95[[2]]
ring.mean.rarity_Jaccard.95 <- ring.mean.rarity_list_Jaccard.95[[3]]
ring.mean.rarity_Jaccard.95_contrasted <- ring.mean.rarity_list_Jaccard.95[[4]]

# Save mimicry ring ranges 
save(ring.range_size_phylo, file = "./outputs/Indices_maps/phylo_only/ring.range_size_phylo.RData", version = "2")
saveRDS(ring.range_size_phylo, file = "./outputs/Indices_maps/phylo_only/ring.range_size_phylo.rds", version = "2")

save(ring.range_size_no_binary, file = "./outputs/Indices_maps/no_binary/ring.range_size_no_binary.RData", version = "2")
saveRDS(ring.range_size_no_binary, file = "./outputs/Indices_maps/no_binary/ring.range_size_no_binary.rds", version = "2")

save(ring.range_size_TSS.80, file = "./outputs/Indices_maps/TSS.80/ring.range_size_TSS.80.RData", version = "2")
saveRDS(ring.range_size_TSS.80, file = "./outputs/Indices_maps/TSS.80/ring.range_size_TSS.80.rds", version = "2")

save(ring.range_size_Jaccard.95, file = "./outputs/Indices_maps/Jaccard.95/ring.range_size_Jaccard.95.RData", version = "2")
saveRDS(ring.range_size_Jaccard.95, file = "./outputs/Indices_maps/Jaccard.95/ring.range_size_Jaccard.95.rds", version = "2")

# Save rarity weights
save(ring.rarity_Leroy_weights_phylo, file = "./outputs/Indices_maps/phylo_only/ring.rarity_Leroy_weights_phylo.RData", version = "2")
saveRDS(ring.rarity_Leroy_weights_phylo, file = "./outputs/Indices_maps/phylo_only/ring.rarity_Leroy_weights_phylo.rds", version = "2")

save(ring.rarity_Leroy_weights_no_binary, file = "./outputs/Indices_maps/no_binary/ring.rarity_Leroy_weights_no_binary.RData", version = "2")
saveRDS(ring.rarity_Leroy_weights_no_binary, file = "./outputs/Indices_maps/no_binary/ring.rarity_Leroy_weights_no_binary.rds", version = "2")

save(ring.rarity_Leroy_weights_TSS.80, file = "./outputs/Indices_maps/TSS.80/ring.rarity_Leroy_weights_TSS.80.RData", version = "2")
saveRDS(ring.rarity_Leroy_weights_TSS.80, file = "./outputs/Indices_maps/TSS.80/ring.rarity_Leroy_weights_TSS.80.rds", version = "2")

save(ring.rarity_Leroy_weights_Jaccard.95, file = "./outputs/Indices_maps/Jaccard.95/ring.rarity_Leroy_weights_Jaccard.95.RData", version = "2")
saveRDS(ring.rarity_Leroy_weights_Jaccard.95, file = "./outputs/Indices_maps/Jaccard.95/ring.rarity_Leroy_weights_Jaccard.95.rds", version = "2")

# Save mean rarity maps
save(ring.mean.rarity_phylo, file = "./outputs/Indices_maps/phylo_only/ring.mean.rarity_phylo.RData", version = "2")
saveRDS(ring.mean.rarity_phylo, file = "./outputs/Indices_maps/phylo_only/ring.mean.rarity_phylo.rds", version = "2")
save(ring.mean.rarity_phylo_contrasted, file = "./outputs/Indices_maps/phylo_only/ring.mean.rarity_phylo_contrasted.RData", version = "2")
saveRDS(ring.mean.rarity_phylo_contrasted, file = "./outputs/Indices_maps/phylo_only/ring.mean.rarity_phylo_contrasted.rds", version = "2")

save(ring.mean.rarity_no_binary, file = "./outputs/Indices_maps/no_binary/ring.mean.rarity_no_binary.RData", version = "2")
saveRDS(ring.mean.rarity_no_binary, file = "./outputs/Indices_maps/no_binary/ring.mean.rarity_no_binary.rds", version = "2")
save(ring.mean.rarity_no_binary_contrasted, file = "./outputs/Indices_maps/no_binary/ring.mean.rarity_no_binary_contrasted.RData", version = "2")
saveRDS(ring.mean.rarity_no_binary_contrasted, file = "./outputs/Indices_maps/no_binary/ring.mean.rarity_no_binary_contrasted.rds", version = "2")

save(ring.mean.rarity_TSS.80, file = "./outputs/Indices_maps/TSS.80/ring.mean.rarity_TSS.80.RData", version = "2")
saveRDS(ring.mean.rarity_TSS.80, file = "./outputs/Indices_maps/TSS.80/ring.mean.rarity_TSS.80.rds", version = "2")
save(ring.mean.rarity_TSS.80_contrasted, file = "./outputs/Indices_maps/TSS.80/ring.mean.rarity_TSS.80_contrasted.RData", version = "2")
saveRDS(ring.mean.rarity_TSS.80_contrasted, file = "./outputs/Indices_maps/TSS.80/ring.mean.rarity_TSS.80_contrasted.rds", version = "2")

save(ring.mean.rarity_Jaccard.95, file = "./outputs/Indices_maps/Jaccard.95/ring.mean.rarity_Jaccard.95.RData", version = "2")
saveRDS(ring.mean.rarity_Jaccard.95, file = "./outputs/Indices_maps/Jaccard.95/ring.mean.rarity_Jaccard.95.rds", version = "2")
save(ring.mean.rarity_Jaccard.95_contrasted, file = "./outputs/Indices_maps/Jaccard.95/ring.mean.rarity_Jaccard.95_contrasted.RData", version = "2")
saveRDS(ring.mean.rarity_Jaccard.95_contrasted, file = "./outputs/Indices_maps/Jaccard.95/ring.mean.rarity_Jaccard.95_contrasted.rds", version = "2")


### 3.6/ Weighted mean ring size ####

add_continental_null_values <- function(x)
{
  continent_mask <- readRDS(file = paste0("./input_data/Map_stuff/continent_mask.rds"))
  
  y <- continent_mask  # Create final new raster from continental mask
  y[!is.na(x[])] <- x[!is.na(x[])]  # Add initial raster values
  
  return(y)
}

# Function 

compute_weighted_mean_ring_size <- function(ring_proba_stack, ring_rich_stack, set)
{
  ring_proba_brick <- brick(ring_proba_stack)
  ring_rich_brick <- brick(ring_rich_stack)
  
  weighted_mean_ring_size <- vector()
  for (i in 1:nrow(ring_rich_brick[]))  # Loop between communities/pixels
  {
    # i <- 55000
    
    rich_vector <- ring_rich_brick@data@values[i,]   # Extract ring richness
    proba_vector <- ring_proba_brick@data@values[i,] # Extract ring proba of presence
    
    # Compute mean weighted by proba. of presence
    weighted_mean_ring_size[i] <- weighted.mean(x = rich_vector, w = proba_vector)
  }
  
  # Add NULL value to continental pixels
  weighted_mean_ring_size <- add_continental_null_values(weighted_mean_ring_size)
  
  # Quick plot
  plot(weighted_mean_ring_size, col = pal_bl_red_Mannion, main = paste0("Mean ring size\n", set))
  
  return(weighted_mean_ring_size)
}

# Load stuff
All_ring_proba_stack_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_ring_proba_stack_phylo.rds"))
All_ring_rich_stack_phylo <- readRDS(file = paste0("./outputs/Indices_stacks/phylo_only/All_ring_rich_stack_phylo.rds"))

All_ring_proba_stack_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_ring_proba_stack_no_binary.rds"))
All_ring_rich_stack_no_binary <- readRDS(file = paste0("./outputs/Indices_stacks/no_binary/All_ring_rich_stack_no_binary.rds"))

All_ring_proba_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_TSS.80.RData"))
All_ring_rich_stack_TSS.80 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_TSS.80.RData"))

All_ring_proba_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_proba_stack_Jaccard.95.RData"))
All_ring_rich_stack_Jaccard.95 <- readRDS(file = paste0("./outputs/Indices_stacks/All_ring_rich_stack_Jaccard.95.RData"))


# Computation
weighted_mean_ring_size_phylo <- compute_weighted_mean_ring_size(ring_proba_stack = All_ring_proba_stack_phylo,
                                                                 ring_rich_stack = All_ring_rich_stack_phylo,
                                                                 set = "only_phylo")
weighted_mean_ring_size_no_binary <- compute_weighted_mean_ring_size(ring_proba_stack = All_ring_proba_stack_no_binary,
                                                                     ring_rich_stack = All_ring_rich_stack_no_binary,
                                                                     set = "no_binary")
weighted_mean_ring_size_TSS.80 <- compute_weighted_mean_ring_size(ring_proba_stack = All_ring_proba_stack_TSS.80,
                                                                  ring_rich_stack = All_ring_rich_stack_TSS.80,
                                                                  set = "TSS.80")
weighted_mean_ring_size_Jaccard.95 <- compute_weighted_mean_ring_size(ring_proba_stack = All_ring_proba_stack_Jaccard.95,
                                                                      ring_rich_stack = All_ring_rich_stack_Jaccard.95,
                                                                      set = "Jaccard.95")
# Save
save(weighted_mean_ring_size_phylo, file = paste0("./outputs/Indices_maps/phylo_only/weighted_mean_ring_size_phylo.RData"), version = "2")
saveRDS(weighted_mean_ring_size_phylo, file = paste0("./outputs/Indices_maps/phylo_only/weighted_mean_ring_size_phylo.rds"), version = "2")

save(weighted_mean_ring_size_no_binary, file = paste0("./outputs/Indices_maps/no_binary/weighted_mean_ring_size_no_binary.RData"), version = "2")
saveRDS(weighted_mean_ring_size_no_binary, file = paste0("./outputs/Indices_maps/no_binary/weighted_mean_ring_size_no_binary.rds"), version = "2")

save(weighted_mean_ring_size_TSS.80, file = paste0("./outputs/Indices_maps/TSS.80/weighted_mean_ring_size_TSS.80.RData"), version = "2")
saveRDS(weighted_mean_ring_size_TSS.80, file = paste0("./outputs/Indices_maps/TSS.80/weighted_mean_ring_size_TSS.80.rds"), version = "2")

save(weighted_mean_ring_size_Jaccard.95, file = paste0("./outputs/Indices_maps/Jaccard.95/weighted_mean_ring_size_Jaccard.95.RData"), version = "2")
saveRDS(weighted_mean_ring_size_Jaccard.95, file = paste0("./outputs/Indices_maps/Jaccard.95/weighted_mean_ring_size_Jaccard.95.rds"), version = "2")


### 4/ Map indices based on phylogeny only ####

# A/ Species richness
# B/ Mean Species geographic rarity
# C/ Faith's Phylogenetic Diversity

# D/ Mimicry richness
# E/ Mean mimicry geographic rarity
# F/ Weighted mean ring size

### 4.0/ Load map stuff ####

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

load(file = "./input_data/Map_stuff/country_borders.RData")
country_borders <- as(country_borders, "Spatial")
# plot(country_borders)

# Choose your alternative set for robustness analyses
set <- "only_phylo"
set <- "no_binary"
set <- "TSS.80"
set <- "Jaccard.95"

### 4.1/ Load directly the index maps ####

if (set == "only_phylo")
{
  sp.richness <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/sp.richness_phylo.rds"))
  sp.mean.rarity_contrasted <- readRDS(file = "./outputs/Indices_maps/phylo_only/sp.mean.rarity_phylo_contrasted.rds")
  Faith_PD <- readRDS(file = "./outputs/Indices_Maps/PD.raster_Jaccard.80.rds")
  ring.richness <- readRDS(file = paste0("./outputs/Indices_maps/phylo_only/ring.richness_phylo.rds"))
  ring.mean.rarity <- readRDS(file = "./outputs/Indices_maps/phylo_only/ring.mean.rarity_phylo.rds")
  weighted_mean_ring_size <- readRDS(file = "./outputs/Indices_maps/phylo_only/weighted_mean_ring_size_phylo.rds")
}

if (set == "no_binary")
{
  sp.richness <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/sp.richness_no_binary.rds"))
  sp.mean.rarity_contrasted <- readRDS(file = "./outputs/Indices_maps/no_binary/sp.mean.rarity_no_binary_contrasted.rds")
  Faith_PD <- readRDS(file = "./outputs/Indices_Maps/no_binary/PD_no_binary.rds")
  ring.richness <- readRDS(file = paste0("./outputs/Indices_maps/no_binary/ring.richness_no_binary.rds"))
  ring.mean.rarity <- readRDS(file = "./outputs/Indices_maps/no_binary/ring.mean.rarity_no_binary.rds")
  weighted_mean_ring_size <- readRDS(file = "./outputs/Indices_maps/no_binary/weighted_mean_ring_size_no_binary.rds")
}

if (set == "TSS.80")
{
  sp.richness <- readRDS(file = paste0("./outputs/Indices_maps/TSS.80/sp.richness_TSS.80.rds"))
  sp.mean.rarity_contrasted <- readRDS(file = "./outputs/Indices_maps/TSS.80/sp.mean.rarity_TSS.80_contrasted.rds")
  Faith_PD <- readRDS(file = "./outputs/Indices_Maps/TSS.80/PD_TSS.80.rds")
  ring.richness <- readRDS(file = paste0("./outputs/Indices_maps/TSS.80/ring.richness_TSS.80.rds"))
  ring.mean.rarity <- readRDS(file = "./outputs/Indices_maps/TSS.80/ring.mean.rarity_TSS.80.rds")
  weighted_mean_ring_size <- readRDS(file = "./outputs/Indices_maps/TSS.80/weighted_mean_ring_size_TSS.80.rds")
}

if (set == "Jaccard.95")
{
  sp.richness <- readRDS(file = paste0("./outputs/Indices_maps/Jaccard.95/sp.richness_Jaccard.95.rds"))
  sp.mean.rarity_contrasted <- readRDS(file = "./outputs/Indices_maps/Jaccard.95/sp.mean.rarity_Jaccard.95_contrasted.rds")
  Faith_PD <- readRDS(file = "./outputs/Indices_Maps/Jaccard.95/PD_Jaccard.95.rds")
  ring.richness <- readRDS(file = paste0("./outputs/Indices_maps/Jaccard.95/ring.richness_Jaccard.95.rds"))
  ring.mean.rarity <- readRDS(file = "./outputs/Indices_maps/Jaccard.95/ring.mean.rarity_Jaccard.95.rds")
  weighted_mean_ring_size <- readRDS(file = "./outputs/Indices_maps/Jaccard.95/weighted_mean_ring_size_Jaccard.95.rds")
}

### 4.2/ Plotting and projection functions ####

Mollweide_shp_projection <-  function(x) # Shp to project
{
  x_name <- deparse(substitute(x)) # Get the name of the initial shp as a character string
  
  new_shp <- spTransform(x, CRSobj = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
  
  # Generate new object with "_Mollweide" suffix in the global environment
  eval(call("<<-", as.name(paste0(x_name, "_Mollweide")), new_shp))
}


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


Mollweide_projection(sp.richness)
Mollweide_projection(sp.mean.rarity_contrasted)
Mollweide_projection(Faith_PD)
Mollweide_projection(ring.richness)
Mollweide_projection(ring.mean.rarity)
Mollweide_projection(weighted_mean_ring_size)

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

### 4.3/ Plot and export ####
species_richness_legend_breaks <- seq(0, 120, 20)
if(set == "no_binary") { species_richness_legend_breaks <- seq(0, 80, 20) }

ring_rarity_legend_breaks <- seq(0, 0.6, 0.2)
if(set == "Jaccard.95") { ring_rarity_legend_breaks <- seq(0, 0.4, 0.1) }

pdf(file = paste0("./supplementaries/Index_maps_",set,".pdf"), height = 8, width = 12)

internal_margins <- par()$mar
# par(mar = c(3.1, 3.5, 3.5, 2.1))
par(mar = c(3.1, 3.1, 2.7, 1.6))
par(mfrow = c(2, 3))

# A/ Species richness ####

map_indices_Mollweide(x = sp.richness_Mollweide,
                      main_title = "Species richness",
                      legend_title = "Species",
                      legend_title_x = -3650,
                      legend_title_y = 430,
                      legend_breaks = species_richness_legend_breaks, 
                      facet_letter = "(a)")

# B/ Species mean geographic rarity ####

map_indices_Mollweide(x = sp.mean.rarity_contrasted_Mollweide,
                      main_title = "Species geographic rarity",
                      legend_title = "Rarity\nindex",
                      legend_title_x = -3650,
                      legend_title_y = 870,
                      legend_breaks = seq(0, 0.5, 0.1), 
                      facet_letter = "(b)")

# C/ Faith's PD ####

map_indices_Mollweide(x = Faith_PD_Mollweide,
                      main_title = "Phylogenetic diversity",
                      legend_title = "Evolutionary\nTime (My)",
                      legend_title_x = -3150,
                      legend_title_y = 670,
                      legend_breaks = seq(0, 800, 200), 
                      facet_letter = "(c)")

# D/ Mimicry richness ####

map_indices_Mollweide(x = ring.richness_Mollweide,
                      main_title = "Mimicry richness",
                      legend_title = "Mimicry\nrings",
                      legend_title_x = -3670,
                      legend_title_y = 670,
                      legend_breaks = seq(0, 25, 5), 
                      facet_letter = "(d)")

# E/ Mimicry ring mean geographic rarity ####

map_indices_Mollweide(x = ring.mean.rarity_Mollweide,
                      main_title = "Mimicry geographic rarity",
                      legend_title = "Rarity\nindex",
                      legend_title_x = -3650,
                      legend_title_y = 670,
                      legend_breaks = ring_rarity_legend_breaks, 
                      facet_letter = "(e)")

# F/ Mean ring size ####

map_indices_Mollweide(x = weighted_mean_ring_size_Mollweide,
                      main_title = "Mean mimicry ring size",
                      legend_title = "Ring\nsize",
                      legend_title_x = -3750,
                      legend_title_y = 670,
                      legend_breaks = seq(0, 6, 1), 
                      facet_letter = "(f)")

par(mar = internal_margins)
par(mfrow = c(1, 1))

dev.off()


