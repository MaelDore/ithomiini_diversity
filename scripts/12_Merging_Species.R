

##### Script 12: Merging species #####

# Create maps at species level by computing "probability" of presence of any of the OMU of the species

# Generate a continuous map for each 4 option between Jaccard/TSS and Buffer.80/95

##### 

# Inputs:
   # Stack of continuous EM and binary maps for each OMU, per species from script 11

# Outputs:
   # Species continuous map for each 4 option between Jaccard/TSS and Buffer.80/95
   # Summary of composition of species stack: with/without rasterized OMU, only rasterized OMU

#####


### Prepare stuff

# Effacer l'environnement
rm(list = ls())

library(raster)

# Choose the resolution
# res <-  "5"
res <- "15"

# Load environmental stack to use as reference for extent and CRS and generate a mask for continental borders
library(raster)
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))
continent_mask <- envData[[1]]
continent_mask <- calc(continent_mask, fun = function(x, ...) {x*0}, na.rm = F)

# Load Summary table for OMU/unit and for Species
load(file = paste0("./input_data/list.sp.RData"))
load(file = paste0("./input_data/list.models.RData"))

# Make a list of modeled and non-modeled units (occ.unit = with only occurrences points)
modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]
rasterized_OMU <- list.models[is.na(list.models$Model_ID), ]

# Function to compute probability of presence with multiple OMUs
aggreg_prob = function(x, na.rm) { 
  y <- 1-prod(1-x) # Probability of presence of species = probability of presence of at least one OMU = opposite of probability of absence of all OMU
  return(y) # Output
}

list.sp$one.Binaries <- F # To record which species include at least one binaries maps
list.sp$all.Binaries <- F # To record which species include only binaries maps

### Loop by species
for (i in 1:nrow(list.sp)) 
{
  # i <- 23
  
  # Load a random layer to initiate the stacks (to remove later)
  sp.stack_Jaccard.80 <- sp.stack_Jaccard.95 <- sp.stack_TSS.80 <- sp.stack_TSS.95 <- stack(continent_mask) 
  
  # Load Species name
  sp <- as.character(list.sp$Sp_full[i]) 
  
  # Create directory to store outputs per species
  if(!dir.exists(paste0("./outputs/By_species/",sp))) { # Test if the folder exists already or not
    dir.create(paste0("./outputs/By_species/",sp), recursive = T) # Create folder if absent
  }
  
  # Create directory to store maps per species
  if(!dir.exists(paste0("./maps/By_species/",sp))) { # Test if the folder exists already or not
    dir.create(paste0("./maps/By_species/",sp), recursive = T) # Create folder if absent
  }
  
  OMU_list <- modeled_OMU[modeled_OMU$Sp_full == sp, ]
  if (nrow(OMU_list) > 0) {      # If at least one OMU was modeled
    
    for (j in 1:nrow(OMU_list))  # For each OMU
    { 
      unit <-  as.character(OMU_list$Tag.model[j]) # Load the unit name                  
      
      # Jaccard.80: Load the continuous maps and stack them all
      Ensemble_Jaccard_median_cropped_80 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_Jaccard_median_cropped_80.rds"))
      sp.stack_Jaccard.80 <- addLayer(sp.stack_Jaccard.80, Ensemble_Jaccard_median_cropped_80)

      # Jaccard.95: Load the continuous maps and stack them all
      Ensemble_Jaccard_median_cropped_95 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_Jaccard_median_cropped_95.rds"))
      sp.stack_Jaccard.95 <- addLayer(sp.stack_Jaccard.95, Ensemble_Jaccard_median_cropped_95)
      
      # TSS.80: Load the continuous maps and stack them all
      Ensemble_TSS_median_cropped_80 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_TSS_median_cropped_80.rds")) 
      sp.stack_TSS.80 <- addLayer(sp.stack_TSS.80, Ensemble_TSS_median_cropped_80)   
      
      # TSS.95: Load the continuous maps and stack them all
      Ensemble_TSS_median_cropped_95 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_TSS_median_cropped_95.rds")) 
      sp.stack_TSS.95 <- addLayer(sp.stack_TSS.95, Ensemble_TSS_median_cropped_95)
    }
  }
  
  Binaries <- rasterized_OMU[rasterized_OMU$Sp_full == sp,]
  if (nrow(Binaries) > 0) {  # If at least one OMU was rasterized
    
    list.sp$one.Binaries[list.sp$Sp_full == sp] <- T
    
    for (j in 1:nrow(Binaries)) {  # For each  OMU
      
      unit <-  as.character(Binaries$Tag.model[j]) # Load the unit name 
      
      # Create directory to store outputs per OMU
      if(!dir.exists(paste0("./outputs/By_OMU/",unit))) { # Test if the folder exists already or not
        dir.create(paste0("./outputs/By_OMU/",unit), recursive = T) # Create folder if absent
      }
      
      # Create directory to store maps per OMU
      if(!dir.exists(paste0("./maps/By_OMU/",unit))) { # Test if the folder exists already or not
        dir.create(paste0("./maps/By_OMU/",unit), recursive = T) # Create folder if absent
      }
      
      # # Load the rasterized occurrence map
      # Raster_Occ <- readRDS(file = paste0("./input_data/Species_data/",res,"/Binary_Rasters_from_Occurrences/",unit,".rds"))
      # 
      # # Add the continental borders as 0 values for empty pixels (oceans are NA)
      # rasterized_map <- continent_mask
      # rasterized_map[!is.na(Raster_Occ[])] <- 1
      # 
      # save(rasterized_map, file = paste0("./outputs/By_OMU/",unit,"/rasterized_map_",unit,".RData"), version = "2")
      # saveRDS(rasterized_map, file = paste0("./outputs/By_OMU/",unit,"/rasterized_map_",unit,".rds"), version = "2")
      # 
      # pdf(file = paste0("./maps/By_OMU/",unit,"/rasterized_map_",unit,".pdf"), height = 6, width = 7)
      # plot(rasterized_map, main = unit)
      # dev.off()
      # # Copy in sp folder
      # file.copy(from = paste0("./maps/By_OMU/",unit,"/rasterized_map_",unit,".pdf"), to = paste0("./maps/By_species/",sp,"/rasterized_map_",unit,".pdf"), overwrite = T)
      
      # Load the rasterized map, with continental borders
      rasterized_map <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/rasterized_map_",unit,".rds"))
      
      # Jaccard.80: Load the binary maps and stack them all
      sp.stack_Jaccard.80 <- addLayer(sp.stack_Jaccard.80, rasterized_map)

      # Jaccard.95: Load the binary maps and stack them all
      sp.stack_Jaccard.95 <- addLayer(sp.stack_Jaccard.95, rasterized_map)
      
      # TSS.80: Load the binary maps and stack them all
      sp.stack_TSS.80 <- addLayer(sp.stack_TSS.80, rasterized_map)   
      
      # TSS.95: Load the binary maps and stack them all
      sp.stack_TSS.95 <- addLayer(sp.stack_TSS.95, rasterized_map)
      
    }
  }
  
  # Note case with only rasterized OMU
  if ((nrow(Binaries) > 0) & (nrow(OMU_list) == 0)) {
    list.sp$all.Binaries[list.sp$Sp_full == sp] <- T
  }
  
  # Drop the useless first layer used to initiate the stack
  sp.stack_Jaccard.80 <- dropLayer(sp.stack_Jaccard.80, i = 1)
  sp.stack_Jaccard.95 <- dropLayer(sp.stack_Jaccard.95, i = 1)
  sp.stack_TSS.80 <- dropLayer(sp.stack_TSS.80, i = 1)
  sp.stack_TSS.95 <- dropLayer(sp.stack_TSS.95, i = 1)
  
  # Compute probability map per Species
  
  sp.cont_Jaccard.80 <- calc(sp.stack_Jaccard.80, fun = aggreg_prob)
  sp.cont_Jaccard.95 <- calc(sp.stack_Jaccard.95, fun = aggreg_prob)
  sp.cont_TSS.80 <- calc(sp.stack_TSS.80, fun = aggreg_prob)
  sp.cont_TSS.95 <- calc(sp.stack_TSS.95, fun = aggreg_prob)
  
  plot(sp.cont_TSS.95)
  
  save(sp.cont_Jaccard.80, file = paste0("./outputs/By_species/",sp,"/cont_Jaccard.80_",sp,".RData"), version = "2")
  saveRDS(sp.cont_Jaccard.80, file = paste0("./outputs/By_species/",sp,"/cont_Jaccard.80_",sp,".rds"), version = "2")
  save(sp.cont_Jaccard.95, file = paste0("./outputs/By_species/",sp,"/cont_Jaccard.95_",sp,".RData"), version = "2")
  saveRDS(sp.cont_Jaccard.95, file = paste0("./outputs/By_species/",sp,"/cont_Jaccard.95_",sp,".rds"), version = "2")
  save(sp.cont_TSS.80, file = paste0("./outputs/By_species/",sp,"/cont_TSS.80_",sp,".RData"), version = "2")
  saveRDS(sp.cont_TSS.80, file = paste0("./outputs/By_species/",sp,"/cont_TSS.80_",sp,".rds"), version = "2")
  save(sp.cont_TSS.95, file = paste0("./outputs/By_species/",sp,"/cont_TSS.95_",sp,".RData"), version = "2")
  saveRDS(sp.cont_TSS.95, file = paste0("./outputs/By_species/",sp,"/cont_TSS.95_",sp,".rds"), version = "2")
  
  pdf(file = paste0("./maps/By_species/",sp,"/cont_Jaccard.80_",sp,".pdf"), height = 6, width = 7)
  plot(sp.cont_Jaccard.80, main = sp)
  dev.off()

  pdf(file = paste0("./maps/By_species/",sp,"/cont_Jaccard.95_",sp,".pdf"), height = 6, width = 7)
  plot(sp.cont_Jaccard.95, main = sp)
  dev.off()
  
  pdf(file = paste0("./maps/By_species/",sp,"/cont_TSS.80_",sp,".pdf"), height = 6, width = 7)
  plot(sp.cont_TSS.80, main = sp)
  dev.off()
  
  pdf(file = paste0("./maps/By_species/",sp,"/cont_TSS.95_",sp,".pdf"), height = 6, width = 7)
  plot(sp.cont_TSS.95, main = sp)
  dev.off()
  
  print(i)
  
}

# save(list.sp, file = paste0(internal.wd,"/list.sp.RData"))

sum(list.sp$one.Binaries) # 170 species among the 388 have at least one binary map used
sum(list.sp$all.Binaries) # 63 species among the 388 are just modeled under binary maps


##### Generate stack of all OMUs outputs for all options


# Load a random layer to initiate the stacks (to remove later)
All_OMU_stack_Jaccard.80 <- All_OMU_stack_Jaccard.95 <- All_OMU_stack_TSS.80 <- All_OMU_stack_TSS.95 <- stack(continent_mask)

### Loop by OMU
for (i in 1:nrow(list.models)) 
{
  # i <- 23
  
  # Load OMU name
  unit <- as.character(list.models$Tag.model[i]) 
  
  if (list.models$initial_model_type[i] != "rasterized") # For modeled OMUs
  { 
    # Jaccard.80: Load the continuous maps and stack them all
    Ensemble_Jaccard_median_cropped_80 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_Jaccard_median_cropped_80.rds"))
    All_OMU_stack_Jaccard.80 <- addLayer(All_OMU_stack_Jaccard.80, Ensemble_Jaccard_median_cropped_80)
    
    # Jaccard.95: Load the continuous maps and stack them all
    Ensemble_Jaccard_median_cropped_95 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_Jaccard_median_cropped_95.rds"))
    All_OMU_stack_Jaccard.95 <- addLayer(All_OMU_stack_Jaccard.95, Ensemble_Jaccard_median_cropped_95)
    
    # TSS.80: Load the continuous maps and stack them all
    Ensemble_TSS_median_cropped_80 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_TSS_median_cropped_80.rds")) 
    All_OMU_stack_TSS.80 <- addLayer(All_OMU_stack_TSS.80, Ensemble_TSS_median_cropped_80)   
    
    # TSS.95: Load the continuous maps and stack them all
    Ensemble_TSS_median_cropped_95 <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_TSS_median_cropped_95.rds")) 
    All_OMU_stack_TSS.95 <- addLayer(All_OMU_stack_TSS.95, Ensemble_TSS_median_cropped_95)
  
  } else  # For rasterized OMU
  { 
    # Load the unique rasterized binary map for this OMU
    rasterized_map <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/rasterized_map_",unit,".rds"))
    
    # Jaccard.80: Load the binary maps and stack them all
    All_OMU_stack_Jaccard.80 <- addLayer(All_OMU_stack_Jaccard.80, rasterized_map)
    
    # Jaccard.95: Load the binary maps and stack them all
    All_OMU_stack_Jaccard.95 <- addLayer(All_OMU_stack_Jaccard.95, rasterized_map)
    
    # TSS.80: Load the binary maps and stack them all
    All_OMU_stack_TSS.80 <- addLayer(All_OMU_stack_TSS.80, rasterized_map)   
    
    # TSS.95: Load the binary maps and stack them all
    All_OMU_stack_TSS.95 <- addLayer(All_OMU_stack_TSS.95, rasterized_map)
  }
  
  if (i %% 10 == 0)
  {print(i)}
  
}
  
# Drop the useless first layer used to initiate the stack
All_OMU_stack_Jaccard.80 <- dropLayer(All_OMU_stack_Jaccard.80, i = 1)
All_OMU_stack_Jaccard.95 <- dropLayer(All_OMU_stack_Jaccard.95, i = 1)
All_OMU_stack_TSS.80 <- dropLayer(All_OMU_stack_TSS.80, i = 1)
All_OMU_stack_TSS.95 <- dropLayer(All_OMU_stack_TSS.95, i = 1)

# Rename layers with OMU names
names(All_OMU_stack_Jaccard.80) <- names(All_OMU_stack_Jaccard.95) <- names(All_OMU_stack_TSS.80) <- names(All_OMU_stack_TSS.95) <- as.character(list.models$Tag.model)

plot(All_OMU_stack_Jaccard.80)

save(All_OMU_stack_Jaccard.80, file = paste0("./outputs/Indices_stacks/All_OMU_stack_Jaccard.80.RData"), version = "2")
saveRDS(All_OMU_stack_Jaccard.80, file = paste0("./outputs/Indices_stacks/All_OMU_stack_Jaccard.80.rds"), version = "2")
save(All_OMU_stack_Jaccard.95, file = paste0("./outputs/Indices_stacks/All_OMU_stack_Jaccard.95.RData"), version = "2")
saveRDS(All_OMU_stack_Jaccard.95, file = paste0("./outputs/Indices_stacks/All_OMU_stack_Jaccard.95.rds"), version = "2")
save(All_OMU_stack_TSS.80, file = paste0("./outputs/Indices_stacks/All_OMU_stack_TSS.80.RData"), version = "2")
saveRDS(All_OMU_stack_TSS.80, file = paste0("./outputs/Indices_stacks/All_OMU_stack_TSS.80.rds"), version = "2")
save(All_OMU_stack_TSS.95, file = paste0("./outputs/Indices_stacks/All_OMU_stack_TSS.95.RData"), version = "2")
saveRDS(All_OMU_stack_TSS.95, file = paste0("./outputs/Indices_stacks/All_OMU_stack_TSS.95.rds"), version = "2")  

  

