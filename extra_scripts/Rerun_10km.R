
##### Run all the SDM workflow for a resolution of 10km = 5min #####

# Clean environment
rm(list = ls())


####### Script 1: Import environmental stack #######

# Load libraries
library("raster")

# Chose resolution
res <- "5"
# res <- "15"

Select_env <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))

plot(Select_env)


##### Script 2: Filtering of occurrences #####

### Load table to track model outputs
load(file = "./input_data/list.models.RData")

# Save the version for 5min resolution
saveRDS(list.models, file = "./input_data/list.models_5min.rds")

### 1/ Add initial_model_type to list.models ####

list.models$initial_model_type <- NA
for (i in 1:nrow(list.models)) {
  
  N.obs <- list.models$N.obs_5m_used[i]
  
  if (N.obs < 6) {
    list.models$initial_model_type[i] <- "rasterized"
    
  } else {
    
    if (N.obs < 30) {
      list.models$initial_model_type[i] <- "restricted"
      
    } else {
      list.models$initial_model_type[i] <- "complete"
      
    }
  }
}

save(list.models, file = "./input_data/list.models_5min.RData")



#### Script 3: Generation of PsA datas ####

# Generate set of Pseudo-absences (PsA) for each OMU
# Nb of sets based on model type, based on number of occurrences available
# Pseudo-absences drawn randomly from other samples sites, with a weighting-distance scheme, and a minimum buffer around occurrences


# Inputs:

# Environmental stack from script 01
# Occurrence data for each OMU (Ithomiini_final)

# Outputs:

# Matrix of distances among sample sites
# Table of sampling site coordinates
# Spatial objects with occurrences and PsA coordinates
# PsA table for each OMU to use in biomod2

########## Only 1 PsA set for quick modeling in 5min resolution #############


##### 1/ Compute distance matrix between sampling sites ####

# Remove environment
rm(list = ls())

# Load packages
library("foreach")

# Need to define resolution of final models right now, to filter the proper number of occurrences using the final modeling grid of all outputs

# Chose resolution
res <- "5"
# res <- "15"

# Load associated environmental stack to use as grid to rasterize
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))

# Load occurrence dataset and extract occurrences coordinates
load(file = "./input_data/Databases/Ithomiini_final.RData")
occ.coords <- Ithomiini_final[, c("Longitude","Latitude")] 

# Rasterization of occurrences
raster.occ_full <- rasterize(x = occ.coords, # Occurrence coordinates
                             y = envData, # Env stack grid to apply
                             field = 1, # Value to store in the raster
                             fun = function(x, ...) mean(x)) # Apply mean to keep value = 1 in case of multiple occurrences in one pixel

cell_indices <- vector(mode = "integer", length = nrow(occ.coords))
for (i in 1:nrow(occ.coords)) {
  cell_indices[i] <- cellFromXY(raster.occ_full, occ.coords[i,])
}

# Use ID from 1 to last instead of cell indices
sampling_ID <- as.factor(as.numeric(as.factor(cell_indices)))
# 1834 sampling sites for 15min resolution
# 2843 sampling sites for 5min resolution

# Get coordinates of pixels associated with each occurrence
xy <- foreach (i = seq_along(cell_indices), .combine = rbind) %do% {
  xyFromCell(envData, cell_indices[i])
}

raster_infos <- as.tibble(cbind(cell_indices, sampling_ID, xy))
raster_infos <- raster_infos %>% 
  dplyr::rename(Longitude_raster_5 = x, Latitude_raster_5 = y, cell_indices_5 = cell_indices, sampling_ID_5 = sampling_ID)

# Store infos on occurrences on raster grid in df of occurrences
if (length(dplyr::intersect(names(Ithomiini_final), names(raster_infos))) > 0) {
  Ithomiini_final <- Ithomiini_final %>% 
    select(-names(raster_infos))
}
Ithomiini_final <- Ithomiini_final %>%
  bind_cols(raster_infos)

# save(Ithomiini_final, file = "./input_data/Databases/Ithomiini_final.Rdata")


### Generate table of sampling sites coordinates
sampling.sites_coords_with_index_and_cell_nb <- raster_infos %>% 
  unique() %>% 
  arrange(sampling_ID_5)

save(sampling.sites_coords_with_index_and_cell_nb, file = paste0("./input_data/Sampling_sites/sampling.sites_coords_with_index_and_cell_nb_",res,".RData"))

load(file = paste0("./input_data/Sampling_sites/sampling.sites_coords_with_index_and_cell_nb_",res,".RData"))

sampling.sites_coords <- sampling.sites_coords_with_index_and_cell_nb %>% 
  select(Longitude_raster_5, Latitude_raster_5)

save(sampling.sites_coords, file = paste0("./input_data/Sampling_sites/sampling.sites_coords_",res,".RData"))


### Compute distance matrix of all sampling sites
Sampling.sites.Dist = geosphere::distm(x = sampling.sites_coords)/1000 # Geometric distance on WGS84 ellipsoid, in km

save(Sampling.sites.Dist, file = paste0("./input_data/Sampling_sites/Sampling.sites.Dist_",res,".RData"))


##### 2/ Get a record of sampling sites that do not have associated env. data to discard them from the potential PsA pool ####

load(file = paste0("./input_data/Sampling_sites/Sampling.sites.Dist_",res,".RData"))

# Chose resolution
res <- "5"
# res <- "15"

# Load associated environmental stack to use as grid to rasterize
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))

sampling.sites_coords$env_data <- NA
for (i in 1:nrow(sampling.sites_coords)) {
  
  # Extract env values for each sampling sites
  test <- raster::extract(envData, sampling.sites_coords[i, 1:2])
  
  # If any env values is missing, note to not use this sampling site as potential PsA site
  sampling.sites_coords$env_data[i] <- all(!is.na(test))
}

sum(!sampling.sites_coords$env_data) # 25 sampling sites out of 2843 are not available to draw PsA because they lack env data

save(sampling.sites_coords, file = paste0("./input_data/Sampling_sites/sampling.sites_coords_",res,".RData"))



##### 3/ Generate PA tables: Parallelized version #####

# Changes =

# Add library loading in the loop
# Modify lines that assign value in object outside of the loop : list.models
#    N.obs et dérivés
# Generate a return(list()) for these objects that are not saved at each iteration, typically infos you want to store in a vector/list/table for each iteration
# Add script to generate and close the cluster
# foreach ... %dopar%


# Remove environment
rm(list = ls())

### 3.1/ Set up the cluster ####

library(foreach)
library(doParallel)

# Detect the number of threads
numCores <- detectCores()
# registerDoParallel(numCores)  # To set the nb if threads to use when using %dopar%

cl <- makeCluster(numCores) # Generate the virtual cluster from the number of threads
registerDoParallel(cl)      # To set the virtual cluster when using %dopar%

# getDoParWorkers() # Pour vérifier le nombre de coeurs enregistrés

### 3.2/ Load stuff useful for all units (outside the loop) ####

### Load table to track model outputs
load(file = "./input_data/list.models_5min.RData")
unit.list <- list.models$Tag.model
sp.list <- list.models$Sp_full

### Load occurrence dataset
load(file = "./input_data/Databases/Ithomiini_final.RData")

### Select Env data resolution 

# Chose resolution
res <- "5"
# res <- "15"

# Load environmental stack
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))

### Load sampling sites infos for PsA Generation
load(file = paste0("./input_data/Sampling_sites/Sampling.sites.Dist_",res,".RData"))
load(file = paste0("./input_data/Sampling_sites/sampling.sites_coords_",res,".RData"))

source("./functions/progcombine.R")

# Set seed to ensure repeatability of random PsA draws
set.seed(1)


### 3.3/ Loop to generate PsA table for each units/OMU ####

temp <- foreach (k = seq_along(unit.list), .combine = progcombine_rbind(nreps = length(unit.list)), .init = NULL) %dopar% {
  # temp <- foreach (k = 1:5) %dopar% {
  # k <- 1

# for (k in 1:length(unit.list)) 
# {
  
  unit <- as.character(unit.list[k])
  N.obs <- list.models$N.obs_5m_used[k]
  
  ##### 3.3.1/ PsA generation ####
  
  ### Only for OMU/unit that will be modeled, to save time
  
  if (N.obs < 6) { 
    
    # Warning message for OMU with too low sample size to try a SDM
    cat(paste0(Sys.time(), " - No Pseudo-Absences Generation for ", unit, " = Unit N°",k,"\n"))
    
  } else { # PsA Generation for OMU with at least 6 obs.
    
    cat("\n", as.character(Sys.time()), "----- PsA generation starts for", unit, "= Unit N°",k,"\n")
    
    # load libraries inside the loop for parallelized computation
    
    library(tidyverse)
    library(raster)
    
    # Load occurrences coordiantes and spatial object
    unit.occurrence <- Ithomiini_final[Ithomiini_final$Tag.model == unit, c("Longitude","Latitude")]
    load(file = paste0("./input_data/Species_data/",res,"/Spatial_Points_Objects/occurrences_", unit,".RData"))  
    
    # Get species name for output paths
    sp <- sp.list[k]
    
    # Retrieve sampling sites ID for this OMU
    unit.sampling_ID <- unique(Ithomiini_final[Ithomiini_final$Tag.model == unit, "sampling_ID_5"])
    
    # Extract the minimal distance of sampling sites, to any of the occurrence points of this OMU
    if (length(unit.points) == 1) { # Case with a simple occurrence point
      Closest.dist <- Sampling.sites.Dist[unit.sampling_ID,] 
    } else { # Case with multiple occurrence points
      Closest.dist <- apply(X = Sampling.sites.Dist[unit.sampling_ID,], MARGIN = 2, FUN = min)
    }
    
    # Select only sampling points outside the buffer zone to draw potential PsA
    # Also discard PsA that do not have env data associated with (because of small errors in GPS coordiantes for points lying close to coastlines)
    buffer <- 111.32 # Buffer = 1° = 111.32km
    index.potential.PsA <- which((Closest.dist > buffer) & sampling.sites_coords$env_data)
    
    # Compute weights for probability of drawing a site
    weights <-  Closest.dist[index.potential.PsA]
    # weights <-  1/(weights*weights) # Weights are inversely proportional to the square of distance such as 1/d²
    weights <-  1/(weights) # Weights are inversly proportional to distance such as 1/d
    # weights <- (weights - min(weights)) / (max(weights) - min(weights)) # Normalise to 0 - 1 range
    
    # Generate temporary PA.table to store indices of sampling sites
    
    if (N.obs < 30) {
      # n.PsA.set <- 10 # 10 PsA for OMU with limited N.obs
      n.PsA.set <- 1 # 1 PsA only for quick 10m resolution
    } else { 
      # n.PsA.set <- 3  # 3 PsA sets for OMU with enough N.obs to do CV
      n.PsA.set <- 1  # 1 PsA only for quick 10m resolution
    } 
    
    PA.indices.table <- data.frame(matrix(nrow = N.obs, ncol = n.PsA.set))
    for (i in 1:n.PsA.set) {
      PA.indices.table[,i] <- sample(x = index.potential.PsA, size = N.obs, prob = weights, replace = F)
    }
    PA.list <- as.vector(as.matrix(PA.indices.table)) # Concatenate in a single vector all the indices of sample sites to use as PsA
    
    PA_set <- as.factor(as.vector(matrix(data = 1:n.PsA.set, byrow = T, nrow = N.obs, ncol = n.PsA.set))) # Generate a vector of PsA set attribution
    ggplot.table <- tibble("indices" = PA.list, "PA_set" = PA_set, "dist" = Closest.dist[PA.list], "weights" = 1/Closest.dist[PA.list])
    
    ### 3.3.2/ Retrieve coordinates of PsA and save in sp object ####
    Occ.coord <- unit.points@coords
    PsA.coord <- sampling.sites_coords[PA.list, 1:2]
    names(Occ.coord) <- names(PsA.coord) <- c("x", "y")
    All.coord <- rbind(Occ.coord,PsA.coord)
    
    # Append the spatial object for this unit, adding the PsA
    
    unit.points.PsA <- SpatialPointsDataFrame(All.coord, # Coordinates of occurrences data and PsA
                                              data = data.frame(Presences_Pseudo.Absences = c(rep(x = 1, times = nrow(Occ.coord)), rep(x = NA, times = nrow(PsA.coord)))), # 1 = Presence point, NA = PsA
                                              proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")) # Define CRS
    
    save(unit.points.PsA, file = paste0("./input_data/Species_data/",res,"/Spatial_Points_Objects/occ_&_PsA_", unit,".RData"))
    saveRDS(unit.points.PsA, file = paste0("./input_data/Species_data/",res,"/Spatial_Points_Objects/occ_&_PsA_", unit,".rds"))
    
    ### 3.3.5/ Generate the associated PA.table for Biomod2 ####
    
    # Fill a df with F value. One row = one occurrence or PsA. One column = One PSA run
    PA.table <- data.frame(matrix(data = F, nrow = N.obs*(n.PsA.set+1), ncol = n.PsA.set))
    
    # Fill occurrences point rows with T
    PA.table[1:N.obs,] <- T # All presences need to be taken at all time
    
    # Fill PsA point rows with T, only for their PsA run
    for (i in 1:n.PsA.set) {
      PA.table[(N.obs+((i-1)*N.obs)+1):(N.obs+(i*N.obs)), i] <- T # Take only the 10 PsA associated with the random sampling
    }
    
    names(PA.table) <- paste0("PA", 1:n.PsA.set)
    
    save(PA.table, file = paste0("./input_data/Species_data/",res,"/PA.Tables/PA.table_", unit,".RData"))
    saveRDS(PA.table, file = paste0("./input_data/Species_data/",res,"/PA.Tables/PA.table_", unit,".rds"))
    
    cat(paste0(Sys.time(), " - Pseudo-Absences Generation for ", unit, " = Unit N°",k," - Done\n"))
    
  }
  
}

# Close clusters generated with registerDoParallel(numCores)
stopImplicitCluster()
# Close clusters generated with makeCluster 
stopCluster() 


##### Script 4: Generation of Cross-validation blocks for OMU with N ≥ 30 ######

######## No CV for quick 10m resolution #########

{

# # Remove environment
# rm(list = ls())
# 
# ### 1/ Load stuff useful for all units (outside the loop) ####
# 
# ### Load table to track model outputs
# load(file = "./input_data/list.models_5min.RData")
# unit.list <- list.models$Tag.model
# sp.list <- list.models$Sp_full
# 
# ### Select Env data resolution 
# 
# 
# # Chose resolution
# res <- "5"
# # res <- "15"
# 
# # Load environmental stack
# envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))
# 
# # Set seed to ensure repeatability of random PsA draws
# set.seed(1)
# 
# # Choose the number of CV folds to separate the occurrences and PsA in between
# nb.of.folds = 3
# 
# 
# ##### 2/ Loop to generate cross validation blocks  ######
# 
# for (k in 1:length(unit.list))
# {
#   
#   unit <- as.character(unit.list[k])
#   sp <- sp.list[k]
#   
#   N.obs <- list.models$N.obs_5m_used[k]
#   
#   # Do the CV only for units with at least 30 occurrences
#   if (N.obs < 30) {
#     
#     # Warning message for OMU with less than 30 obs
#     cat("\n", as.character(Sys.time()), "----- No need for CV blocks generation for", unit, "= Unit N°",k,"\n")
#     
#   } else { # CV design for OMU with at least 30 obs
#     
#     cat("\n", as.character(Sys.time()), "----- CV blocks generation starts for", unit, "= Unit N°",k,"\n")
#     
#     # load libraries (needed inside the loop when using parallelization)
#     
#     library(raster)
#     library(sf)
#     library(tmap)
#     library(tmaptools)
#     library(blockCV)
#     
#     # Load spatial points with PsA for this unit/OMU
#     load(file = paste0("./input_data/Species_data/",res,"/Spatial_Points_Objects/occ_&_PsA_", unit,".RData"))
#     
#     # Transform spatial points to an sf object
#     occ_PsA_data <- sf::st_as_sf(unit.points.PsA, coords = c("x", "y"), crs = raster::crs(envData))
#     
#     # Need to create as many folds as there are PsA sets
#     # Load PA.table
#     load(file = paste0("./input_data/Species_data/",res,"/PA.Tables/PA.table_", unit,".RData"))
#     
#     # 2.1/ Generate empty array to record calibration and validation lines (calib.lines in Biomod2) ####
#     
#     # NA = not use in this run. T = calibration/learning set. F = validation/test set.
#     # Start with array full of NA.
#     calib_array <- array(dim = c(nrow(PA.table),   # Number of occurrences and PsA
#                                  nb.of.folds,      # Number of CV folds
#                                  ncol(PA.table)),  # Number of PsA_set
#                          dimnames = list(NULL,                              # No name for point rows
#                                          paste0("_RUN", 1:nb.of.folds),     # Need to add "_" as in BioMod2 calib.lines
#                                          colnames(PA.table))                # Keep PA set names from PA.table
#     )
#     
#     for (pa in colnames(PA.table)) { # Generate CV table for each PA.run
#       
#       cur.occ_PsA_data <- occ_PsA_data[which(PA.table[, pa]), ]  # Get only lines for the current PsA run
#       
#       # Transform PsA into Absences (0) to get the spatialBlock function working and obtain balanced folds with as many Presence and PsA (seen as Absences)
#       cur.occ_PsA_data$Presences_Pseudo.Absences[which(is.na(cur.occ_PsA_data$Presences_Pseudo.Absences))] <- 0
#       
#       # 2.2/ Run the fold repartition process ####
#       cur.blocks <- spatialBlock(speciesData = cur.occ_PsA_data,           # Sf_object with data
#                                  species = "Presences_Pseudo.Absences",    # Column with Presence (1), Absence (0), and PsA (NA)
#                                  rasterLayer = envData[["Elevation"]],     
#                                  theRange = 100000,                        # Block range set to 100km
#                                  k = nb.of.folds,                          # nb of CV folds to generate
#                                  selection = "random",                     # Assign folds to each blocks randomly
#                                  iteration = 100,                          # Nb of iterations allowed to try to find the best sorting which balance Presence and PsA among folds
#                                  showBlocks = F,                           # Don't plot automatically
#                                  biomod2Format = TRUE)                     # Get an output in the Biomod2 format
#       
#       # Results : train_0 = Nb of PsA for calibration/learning
#       #           train_1 = Nb of presences for calibration/learning
#       #           test_0 = Nb of PsA for evaluation/learning
#       #           test_1 = Nb of presences for evaluation/learning
#       
#       # Fill calib.array with the output of spatialBlock
#       
#       calib_array[which(PA.table[, pa]), , pa] <- cur.blocks$biomodTable # Only the lines for this PsA set model
#       
#       # Save folds infos for this PsA set
#       saveRDS(cur.occ_PsA_data, file = paste0("./input_data/Species_data/",res,"/CV_blocks/CV_folds_data/CV_folds_data_", unit,"_",pa,".rds"))
#       
#     }
#     
#     # save CV infos for this OMU
#     save(calib_array, file = paste0("./input_data/Species_data/",res,"/CV_blocks/calib_arrays/calib_array_", unit,".RData"))
#     saveRDS(calib_array, file = paste0("./input_data/Species_data/",res,"/CV_blocks/calib_arrays/calib_array_", unit,".rds"))
#     
#     cat("\n", as.character(Sys.time()), "----- CV blocks generation for", unit, "= Unit N°",k,"- Done\n")
#     
#   }
# }

}

##### Script 5: Models Calibration and Projection ######

# Format data under BioMod2 standards
# Set parameters of algorithms
# Calibrate SDMs
# Project SDMs results under current climatic conditions


# Inputs:

# Environmental stack from script 01
# Summary table for OMU models (list.models)
# Spatial objects with occurrences and PsA from script 03
# Corrected environmental stacks from script 02
# PA table from script 03
# CV scheme from script 04

# Outputs:

# Calibrated models
# Maps of current environmental suitability

#####

####### Only for RF, a single PsA fold, and no CV ! #########

# Remove environment
rm(list = ls())


### Load stuff useful for all units (outside the loop) ####

### Select Env data resolution 
res <- "5"
# res <- "15"

# Load environmental stack
library(raster)
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))

# Load Full unit list
load(file = paste0("./input_data/list.models_5min.RData"))

# Generate new column in the summary table for units to keep tracks of models results

# list.models[, c("Model_ID", "Start_Time", "End_Time", "Models_success", "Models_select_0.5", "Mean_TSS", "Mean_TSS_select_0.5")] <- NA

# New_vars_for_models_summary <- data.frame(
#   Model_ID = rep(NA, nrow(list.models)),              # Internal Biomod model ID
#   Start_Time = rep(NA, nrow(list.models)),            # Starting Date & Time of last modeling
#   End_Time = rep(NA, nrow(list.models)),              # Ending Date & Time of last modeling
#   Models_success = rep(NA, nrow(list.models)),        # Nb of models not failing to converge
#   Models_select_0.5 = rep(NA, nrow(list.models)),     # Nb of models selected for TSS >= 0.5
#   Mean_TSS = rep(NA, nrow(list.models)),              # Mean TSS among all models
#   Mean_TSS_select_0.5 = rep(NA, nrow(list.models)))   # Mean TSS among selected models for TSS >= 0.5
#                                          
# list.models <- cbind(list.models, New_vars_for_models_summary)


### Select a reduced number of units to model ####

index_model_to_compute <- 1:nrow(list.models)
# index_model_to_compute <- which(list.models$N.obs_15m_used > 20)
# index_model_to_compute <- c(14)
# index_model_to_compute <- 1:10

############## Modif between scripts ######################
# index_model_to_compute <- c(1:250)
# index_model_to_compute <- c(251:500)
# index_model_to_compute <- c(501:783)
###########################################################


##### Loop to do models Calibration, Execution, Projection, and Ensemble  ######

# Parallelized version, libraries need to be load in the loop

library(biomod2)

# Save initial working directory
initial.wd <- getwd()

library(raster)
rasterOptions(tmpdir = paste0(initial.wd,"/temp")) # Modify folder to store temporary files for this session (Goal = avoid problem with paths too long)

# Set seed to ensure repetability of algorithm optimization
set.seed(1)

k <- 1 # Index to keep track of number of iterations when not parallelized, and run on limited nb of units
for (i in index_model_to_compute) { # Unparallelized version

  # i <- 4
  # i <- 6
  
  # Load unit/OMU name
  unit <- as.character(list.models$Tag.model[i]) 
  model_type <- list.models$initial_model_type[i]
  
  if (model_type == "rasterized") {
    
    cat(paste("\n", Sys.time(),"------ NO MODELING for", unit, "because of low N.obs \n")) # Warning message for OMU with no SDM
    
  } else {
    
    cat(paste0("\n", Sys.time()," ------ Starting process for ", unit, " = Unit N°",k,"/",length(index_model_to_compute)," ------\n")) # Affichage le temps au début de la procédure
    
    # Unparallelized version
    list.models$Start_Time[list.models$Tag.model == unit] <- as.character(Sys.time()) # Store date and time of starting of last modeling process
    
    # Load Spatial object with occurrences and PsA
    load(file = paste0("./input_data/Species_data/",res,"/Spatial_Points_Objects/occ_&_PsA_", unit,".RData"))
    
    # Load environmental stack with corrections for extrapoled data points
    load(file = paste0("./input_data/Species_data/",res,"/Env_Stacks/Env_stack_",unit,".RData"))
    unit.env
    
    # Load PA table for biomod2
    load(file = paste0("./input_data/Species_data/",res,"/PA.Tables/PA.table_", unit,".RData"))
    
    # Create directory to store model outputs
    if(!dir.exists(paste0("./models/",res,"/", unit))) { # Test if the folder exists already or not
      dir.create(paste0("./models/",res,"/", unit), recursive = T) # Create folder if absent
    }
    
    # Set wd as the models folder since Biomod save automatically outputs in a folder named as the modelled unit name
    setwd(paste0(initial.wd,"/models/",res,"/"))
    
    ################# 1/ Formating data for biomod2 ######################
    
    # Starting warning
    cat(paste0("\n", Sys.time(),"------ Formating input data for ", unit, " = Unit N°",k,"/",length(index_model_to_compute)," ------\n")) 
    
    formated.input.data <- BIOMOD_FormatingData(resp.var = unit.points.PsA,    # Spatial object with occurrences and PsA
                                                expl.var = unit.env,           # Env stack
                                                resp.name = unit,              # Unit/OMU name
                                                PA.nb.rep = ncol(PA.table),    # Nb of PsA sets
                                                PA.nb.absences = sum(unit.points.PsA@data, na.rm = T), # Nb of occurrences = Nb of PsA per set
                                                PA.strategy = "user.defined",  # To allow to provide a custom PA_table
                                                PA.table = PA.table)           # Provide the custom PA_table
    
    save(formated.input.data, file = paste0(initial.wd,"/models/",res,"/",unit,"/formated.input.data.RData"))
    saveRDS(formated.input.data, file = paste0(initial.wd,"/models/",res,"/",unit,"/formated.input.data.rds"))
    
    # Plot allow to check for the location of data in each PA set
    # x11()
    # plot(formated.input.data) 
    
    ################# 2/ Model parametrization and execution ######################
    
    # Starting warning
    cat(paste0("\n", Sys.time(),"------ Modeling for ", unit, " = Unit N°",k,"/",length(index_model_to_compute)," ------\n"))
    
    # Check what is going on with the do.full or DataSplit = 100 pour l'évaluation par la suite et la gestion du nom des modèles
    # Faire un if pour 2 setup différents !
    
    # DataSplitTable erase NbRunEva, DataSplit et do.full
    
    # Quick version: 1 PsA sets, No CV, only RF, for all model types !
    
    model.runs <- BIOMOD_Modeling(data = formated.input.data,
                                  models =  c('RF'), # Algorithm: RF = Random Forest, GBM = Gradient Boosting Model, ANN = Artifical Neural Network
                                  DataSplit = 100, # Proportion of data to use for calibration
                                  NbRunEval = 1, # Nb of evaluation run = 1 because we use 100% of data for calibration. No CV possible for low sample size (N < 30)
                                  VarImport = 0, # Nb of data permutation to evaluate variable importance
                                  models.eval.meth = c("TSS"), # Evaluation metric = True Skill Statistic
                                  do.full.models = TRUE, # Compute model with all data. RUN = "Full"
                                  SaveObj = T, # Save everything generated during the modeling process
                                  Prevalence = NULL, # Weights P/A. Poids relatifs d'une absence par rapport à une présence. Si 0.5 => poids répartis également selon le nb de P/A => somme des poids des P et A est identiques. Lorsque qu'utilisation de pseudo-absences, PR et PA équivalentes par défaut => Prévalence = NULL (= 0.5 par défaut). Si Prévalence > 0.5 => more weights for PR. If Prévalence < 0.5 => more weights for Absences.                                
                                  rescal.all.models = FALSE, # No rescaling of output needed with the last version
                                  modeling.id = , # Empty = generate a new ID each time
                                  models.options = NULL # Default options for RF
    )

    # Save models runs results
    save(model.runs, file = paste0(initial.wd,"/models/",res,"/",unit,"/model.runs.RData")) 
    saveRDS(model.runs, file = paste0(initial.wd,"/models/",res,"/",unit,"/model.runs.rds"))
    
    ### Explore model.runs
    # model.runs
    # str(model.runs)
    # get_variables_importance(model.runs) # Retrieve variable importance stored in a 4D array: importance[variable, Algo, RUN, PA]
    # get_evaluations(model.runs) # Retrieve evaluation metrics stored in a 5D array: eval[metric, infos, Algo, RUN, PA]
    # get_predictions(model.runs, as.data.frame = F) # Get predictions (output proba. for each evaluated data) in a 4D array: pred[obs, Algo, RUN, PA] (in a dataframe if as.data.frame = T)
    # get_formal_data(model.runs, subinfo = 'resp.var') # Vector of input data of Presence (1) /Absences (0)/PsA (NA)
    # get_formal_data(model.runs, subinfo = 'expl.var') # df of environnmental variables
    # get_built_models(model.runs) # Vector of convergent model names
    
    
    # # Store infos from modeling results: Unparallelized version
    list.models$Model_ID[list.models$Tag.model == unit] <- model.runs@modeling.id # Keep model ID
    list.models$Models_success[list.models$Tag.model == unit] <- length(get_built_models(model.runs)) # Save nb of converging models
    
    # Project everything, but do not do the Ensemble. We need to compute Jaccard and select proper models that pass our evaluation check prior to do the Ensemble manually (also allow to compute sd)
    
    
    if (!(length(get_built_models(model.runs)) > 0)) { # Stop process if no model has converged
      
      # Display warning if no model has converged
      cat(paste0("\n", Sys.time()," ------ No successful models for ", unit, " = Unit N°",k,"/",length(index_model_to_compute)," ------\n"))
      
    } else { # Projection of all models that have converged. # No ensembling because proper evaluation of models quality is done outside of BioMod2
      
      ### Save infos on number of models with TSS ≥ 0.5 just to compare with our manual evaluation outside of Biomod2, and have a first quick view of models quality.
      
      # Proper version when script run in local. Seems to produce error on the cluster (problm with reshape2 ?)
      # eval.runs <- subset(x=reshape2::melt(get_evaluations(model.runs, as.data.frame = F)), subset = (Var2 == "Testing.data"))[, c(3:6)]
      
      eval.runs <- get_evaluations(model.runs, as.data.frame = T)
      eval.runs$value <- eval.runs$Testing.data
      
      # Store info on Biomod2 Evaluation with TSS: Unparallelized version
      list.models$Mean_TSS_full[list.models$Tag.model == unit] <- round(mean(eval.runs$value, na.rm = T),3) # Mean value for TSS among all converged models
      list.models$Models_select_TSS_0.5[list.models$Tag.model == unit] <- sum(eval.runs$value >= 0.5, na.rm = T) # Nb of model with TSS ≥ 0.5
      list.models$Mean_TSS_select_TSS_0.5[list.models$Tag.model == unit] <- round(mean(eval.runs$value[eval.runs$value >= 0.5], na.rm = T),3) # Mean value for TSS among all models with TSS ≥ 0.5
      
      ########## 3/ Model projection ##############
      
      cat(paste0("\n", Sys.time()," ------ Projection of all models for ", unit, " = Unit N°",k,"/",length(index_model_to_compute)," ------\n"))
      
      # Projection of each "submodel" = RUN * PA * Algo
      projection.runs <- BIOMOD_Projection(modeling.output = model.runs, # Output of BIOMOD_Modeling
                                           new.env = stack(envData), # Environnmental stack to use for the projection
                                           proj.name = "Current", # Name of the projection based on name of the climatic model
                                           selected.models = 'all', # Use all submodels
                                           binary.meth = NULL, # No thresholding for binary maps. We will do it outside of Biomod2
                                           omit.na = T, # Do not predict pixel with incomplete env data
                                           on_0_1000 = T, # Predictions store on 0 to 1000 scale for LOT OF memory saving
                                           compress = TRUE, # Compress output files
                                           build.clamping.mask = TRUE, # To build a mask showing which predictions are made out of the environmental range out of the data used for calibration => i.e., show the area of extrapolation
                                           output.format = ".RData", # Chose projection output format to save on the HD. Save under .grd quick for vizualisation in Explorer, and QGIS compatibility
                                           do.stack = T, # Store all projected layers (one per model) into a single stack rather than separated
                                           keep.in.memory = F) # To keep only the link to the hard-drive copy of the layer in the output object rather than all the projections
      
      # Save projection runs
      save(projection.runs, file = paste0(initial.wd,"/models/",res,"/",unit,"/proj_Current/projection.runs.RData"))
      saveRDS(projection.runs, file = paste0(initial.wd,"/models/",res,"/",unit,"/proj_Current/projection.runs.rds"))
      
      # Reload .RData and Save final stack under .rds to be able to load it in any named object
      load(file = paste0(initial.wd,"/models/",res,"/",unit,"/proj_Current/proj_Current_",unit,".RData")) 
      # load(projection.runs@proj@link) # Does the same but retrieve the path to the .RData from the projection.runs object
      saveRDS(eval(parse(text = paste0("proj_Current_",unit))), file = paste0(initial.wd,"/models/",res,"/",unit,"/proj_Current/proj_Current_",unit,".rds"))
      rm(list = paste0("proj_Current_",unit)) # Remove it from the env.
      
      # Reload .RData and Save Clamp mask under .rds to be able to load it in any named object
      load(file = paste0(initial.wd,"/models/",res,"/",unit,"/proj_Current/proj_Current_ClampingMask.RData")) 
      saveRDS(eval(parse(text = paste0("proj_Current_",unit,"_ClampingMask"))), file = paste0(initial.wd,"/models/",res,"/",unit,"/proj_Current/proj_Current_",unit,"_ClampingMask.rds"))
      rm(list = paste0("proj_Current_",unit,"_ClampingMask")) # Remove it from the env.
      
      ##### Explore content
      # str(projection.runs)
      # plot(projection.runs,str.grep = "key") # To plot only projection from model with "key" in their name 
      # projection.runs@proj@link # Link to the final .RData 
      # load(projection.runs@proj@link) # Load the raster stack of projections
      #
      # proj_stack <- readRDS(file = paste0(initial.wd,"/models/",res,"/",unit,"/proj_Current/proj_Current_",unit,".rds"))
      # plot(proj_stack[[1:4]])
      # plot(proj_stack[[grep(pattern = "PA1", x = names(proj_stack))]])
      #####
      
      # Clean the RAM and temp folder to avoid issues
      rm(list = c("eval.runs", "formated.input.data", "model.runs", "unit.env", "unit.points.PsA"))
      unlink(list.files(path = rasterOptions()$tmpdir, all.files = T, recursive = T, include.dirs = T, full.names = T), force = T, recursive = T) # Clean raster store in temp
      
      cat(paste0("\n", Sys.time()," ------ Projections for ", unit, " = Unit N°",k,"/",length(index_model_to_compute)," - Finished ------\n"))
      
    } 
    
    # Store date and time of the end of modeling process: Unparallelized version
    list.models$End_Time[list.models$Tag.model == unit] <- as.character(Sys.time())
    
    # print(warnings()) # To control potential errors
    
  }    
  
  # Set wd back to initial working directory
  setwd(initial.wd)
  
  cat(paste0("\n", Sys.time()," ------ Process over for ", unit, " = Unit N°",k,"/",length(index_model_to_compute)," ------\n")) # Affichage le temps au début de la procédure
  
  # Unparallelized version
  save(list.models, file = "./input_data/list.models_5min.RData")
  
  k <- k + 1
  
}

setwd(initial.wd) # Pour reset le working directory dans le PC

# Unparallelized version
save(list.models, file = "./input_data/list.models_5min.RData")



##### Script 6: Response plot #####

# Not needed for quick modeling at 10km resolution


##### Script 7: Manual evaluation of sub.models #####

# Evaluate the quality of models using optimized Jaccard indices and TSS
# Estimate the threshold that maximize those evaluation metrics

#####

# Inputs:

# Summary table for OMU models (list.models)
# Model results from Script 05.1

# Outputs:

# Jaccard and TSS optimized values, with their associated threshold, for each model

#####  


# Effacer l'environnement
rm(list = ls())

# setwd("D:/Mael/R_projects/ithomiini_diversity/")

# Load libraries

library(biomod2)

# Select Env data resolution 
res <- "5"
# res <- "15"

### 1/ Create empty arrays to gather the evaluation data ####

# Load example OMUs/units to use as template to create the empty array that will gather the evaluation data for all units

# Load Summary table for OMU models
load(file = paste0("./input_data/list.models_5min.RData"))
# modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]
modeled_OMU <- list.models[!(list.models$initial_model_type == "rasterized"), ]

eval_stats <- c("Jaccard", "TSS")

# Create an array to store evaluation for all units
all.units_evaluation <- array(dim = c(nrow(modeled_OMU),  # Nb of units/OMUs
                                      length(eval_stats)),   # Nb of evaluation stats = 2 
                              dimnames = list(unit = modeled_OMU$Tag.model,   # Unit = OMU
                                              eval_stats = eval_stats))       # Names of evaluation stats = Jaccard, TSS      

# Create a similar array to store the cut-off value that maximize the evaluation stat for each sub_model. Usefull for binarization
all.units_cutoffs <- all.units_evaluation

# Save empty files
save(all.units_evaluation, file = "./controls/evaluations/5/all.units_evaluation.RData", version = 2)
saveRDS(all.units_evaluation, file = "./controls/evaluations/5/all.units_evaluation.rds", version = 2)
save(all.units_cutoffs, file = "./controls/evaluations/5/all.units_cutoffs.RData", version = 2)
saveRDS(all.units_cutoffs, file = "./controls/evaluations/5/all.units_cutoffs.rds", version = 2)


### 2/ Loop to evaluate each model with Jaccard and TSS ####

# Effacer l'environnement
rm(list = ls())

# setwd("D:/Mael/R_projects/ithomiini_current/")

# Load libraries
library(biomod2)

# Select Env data resolution 
res <- "5"
# res <- "15"

# Load Summary table for OMU models
load(file = paste0("./input_data/list.models_5min.RData"))
# modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]
modeled_OMU <- list.models[!(list.models$initial_model_type == "rasterized"), ]

# modeled_OMU  <- list.models[1:2,]

# Load summary files to fill
all.units_evaluation <- readRDS(file = "./controls/evaluations/5/all.units_evaluation.rds")
all.units_cutoffs <- readRDS(file = "./controls/evaluations/5/all.units_cutoffs.rds")

# Set up directory
initial.wd <- getwd()
setwd(paste0(initial.wd,"/models/",res,"/"))

for (i in 1:nrow(modeled_OMU)) 
{
  # i <- 1
  
  unit <- as.character(modeled_OMU$Tag.model[i])
  model_type <- modeled_OMU$initial_model_type[i]
  
  # cat(paste0("\n", Sys.time()," ------ Evaluation starts for ", unit, " = Unit N°",i,"/",nrow(modeled_OMU)," ------\n"))
  
  # Retrieve calibration lines to know which data point to use as test set for each sub_model
  model.runs <- readRDS(file = paste0(unit,"/model.runs.rds"))
  
  # Retrieve foramted data, especially the infos on P/PsA
  formated.input.data <- readRDS(file = paste0(unit,"/formated.input.data.rds"))
  # load(model.runs@formated.input.data@link) ; formated.input.data <- data  # Retrieve formated.input.data for infos on PA runs
  input_data <- formated.input.data@data.species # Presence = 1 ; Absence = 0 ; Pseudo-absence = NA
  
  # Retrieve calib.lines. TRUE = Points used for calibration (70%) ; FALSE = Points used for evaluation (30%) ; Points not used for this submodel
  load(model.runs@calib.lines@link)          
  
  # Load models.prediction ; 4D array saving predicts => models.prediction[points, Algo, CV_run, PA_set]
  load(model.runs@models.prediction@link)
  # Predicts go from 0 to 1000
  
  # dimnames(models.prediction)
  
  # Retrieve input values of data use for evaluation
  obs.data <- input_data    
  obs.data[is.na(obs.data)] <- 0 # /!\ Transformation of pseudo-absences in absences needed to generate the confusion matrix
  
  # Retrieve predictions
  cur.preds <- models.prediction[, 1, 1, 1]
  
  # Compute Jaccard and TSS index for all threshold between 0 et 1 (= 1000)
  jaccard.test <- TSS.test <-  NULL # Object to store the optimum thresholds for each metric
  
  
  for(cutoff in seq(0, 1000, by = 1)) # We test all threshold from 0 to 1000 with a step of 1.
  {
    pred.pa <- cur.preds                               # Predictions of the model for the test set
    pred.pa[pred.pa < cutoff] <- 0                     # We assigne a 0 (Absence) for all predicts under the cutoff 
    pred.pa[pred.pa >= cutoff] <- 1                    # We assign a 1 (presence) for all predicts equal or above the cutoff 
    
    # We compare predicts and inputs = confusion matrix
    TP <- sum(obs.data == 1 & pred.pa == 1)   # True positives
    FN <- sum(obs.data == 1 & pred.pa == 0)   # False negatives
    FP <- sum(obs.data == 0 & pred.pa == 1)   # False positives
    TN <- sum(obs.data == 0 & pred.pa == 0)   # True negatives
    
    sensitivity <- round(TP / (TP + FN), 3) # Ability to detect presence among presences points
    specificity <- round(TN / (TN + FP), 3)  # Ability to detect absences among absences points
    
    jaccard <- round(TP / (TP + FP + FN), 3)        # Compute Jaccard index
    TSS <- sensitivity + specificity - 1  # Compute TSS
    
    # Store all infos regarding the cut.off test for Jaccard in a df
    jaccard.test <- rbind.data.frame(jaccard.test,                  
                                     data.frame(cutoff = cutoff,
                                                TP = TP,
                                                FN = FN,
                                                FP = FP,
                                                jaccard = jaccard))
    
    # Store all infos regarding the cut.off test for TSS in a df
    TSS.test <- rbind.data.frame(TSS.test,                  
                                 data.frame(cutoff = cutoff,
                                            sensitivity = sensitivity,
                                            specificity = specificity,
                                            TSS = TSS))
  }
  
  # Find the best cut-off = threshold that maximize the stat
  # We take a mean because several values can lead to equal stat
  all.units_cutoffs[unit, "Jaccard"] <- mean(jaccard.test$cutoff[which(jaccard.test$jaccard == max(jaccard.test$jaccard))])
  all.units_cutoffs[unit, "TSS"] <- mean(TSS.test$cutoff[which(TSS.test$TSS == max(TSS.test$TSS))])
  
  # Save the eval stat value. Will be used to discard bad models prior ensemble, and eventually for binarization
  all.units_evaluation[unit, "Jaccard"] <- jaccard.test$jaccard[
    which(jaccard.test$cutoff == round(all.units_cutoffs[unit, "Jaccard"]))] # We round the mean to get an existing cut-off value
  
  all.units_evaluation[unit, "TSS"] <- TSS.test$TSS[
  which(TSS.test$cutoff == round(all.units_cutoffs[unit, "TSS"]))] # We round the mean to get an existing cut-off value

  if (i %% 10 == 0) 
  { 
    cat(paste0("\n", Sys.time()," ------ Evaluation done for ", unit, " = Unit N°",i,"/",nrow(modeled_OMU)," ------\n"))
  }
  
}

# setwd("D:/Mael/R_projects/ithomiini_current/")
setwd(initial.wd) # Back to initial wd

# Save cut-off data
save(all.units_cutoffs, file = "./controls/evaluations/5/all.units_cutoffs.RData")
saveRDS(all.units_cutoffs, file = "./controls/evaluations/5/all.units_cutoffs.rds")

# Save evaluations
save(all.units_evaluation, file = "./controls/evaluations/5/all.units_evaluation.RData")
saveRDS(all.units_evaluation, file = "./controls/evaluations/5/all.units_evaluation.rds")


##### Script 8: Selection of valid models #####

# Only need to check if Jaccard/TSS are not so bad
# If too bad, must use tha raster because we only have one model per OMU...

hist(all.units_evaluation[, "Jaccard"])
table(all.units_evaluation[, "Jaccard"])

hist(all.units_evaluation[, "TSS"])
table(all.units_evaluation[, "TSS"])

##### Script 9: Variable Importance #####

# No use for quick 10m resolution

##### Script 10: Ensemble #####

# Not needed since we have only one model per OMU


##### Script 11: Clipping maps #####

# Use alpha-hull around occurrences points and a buffer proportional to the quality of presence sampling to clip the maps
# Clip model projection for each OMU

#####

### Inputs:
# Buffer for continental borders
# Buffer for Andean regions
# OMU occurrence coordinates (to draw the alpha-hull) from Script 01
# Model projection from script 05


### Outputs:
# Buffer value based on 80%/95% quantile of occurrences distance to closest point
# Alpha-hull for each OMU
# Plots of final clipped EM maps

#####


# Effacer l'environnement
rm(list = ls())

# Load libraries
library(raster)
library(rgeos)
library(rgdal)
library(alphahull)
library(geosphere)
library(spatialEco)

# Load Summary table for OMU models
load(file = paste0("./input_data/list.models_5min.RData"))
# modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]
modeled_OMU <- list.models[!(list.models$initial_model_type == "rasterized"), ]


# Chose the resolution
res <-  "5"
# res <-  "15"

# Load environmental stack to use as reference for extent and CRS
library(raster)
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))

# Load mask for continent borders
# land_poly_sf <- readRDS(file = "./input_data/Map_stuff/land_poly_sf.rds")
# continent_mask <- rasterize(x = land_poly_sf, y = envData,
#                             field = 1, background = NA, fun = max)
# saveRDS(continent_mask, file = paste0("./input_data/Env_data/continent_mask_", res, ".rds"))

continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_", res, ".rds"))
plot(continent_mask)

# Load buffers for Andes
Andes_East <- readOGR(dsn = paste0("./input_data/biogeo_masks"), layer = "neo_e")
Andes_West <- readOGR(dsn = paste0("./input_data/biogeo_masks"), layer = "neo_w")
# # Rename CRS with proper order for arguments to avoid warnings
# Andes_East@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")
# Andes_West@proj4string <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")

# Rasterize the buffers
Andes_East.rasterized <- raster::rasterize(x = Andes_East, y = envData, field = 1)
Andes_West.rasterized <- raster::rasterize(x = Andes_West, y = envData, field = 1)
# plot(Andes_East.rasterized)
# plot(Andes_West.rasterized, add =T)

# Projection in Mollweide to allow intersection function to work properly
Andes_East <-  spTransform(Andes_East , CRSobj = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
Andes_West <-  spTransform(Andes_West , CRSobj = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")

# plot(Andes_East)
# plot(Andes_West, add = T)

# Create the variable to keep track of location of OMU regarding Andes slopes
list.models$Andes_side <-  NA 

# Functions to convert alpha-hull into spatial objects
source("./functions/Alpha_functions.R", local = TRUE) 


############## Modif between scripts ######################
index_model_to_compute <- c(1:585)
index_model_to_compute <- c(403:585)
# index_model_to_compute <- c(1:200)
# index_model_to_compute <- c(201:400)
# index_model_to_compute <- c(401:585)
###########################################################

### Loop for all OMU

# i <- 161
# Units en erreur à utiliser avec un alpha-shape plutôt que alpha-hull
# for (i in c(161, 173, 178, 207, 221, 237, 249, 270, 280, 293, 318, 341, 388, 444, 480, 494, 511, 528, 548)) {

### Loop for all modeled OMU/unit
# for (i in 1:nrow(modeled_OMU)) 
for (i in index_model_to_compute) 
{
  # i <- 1
  
  unit <- as.character(modeled_OMU$Tag.model[i])
  sp <- as.character(modeled_OMU$Sp_full[i])
  
  cat(paste0("\n", Sys.time()," ------ Starts for ", unit, " = Unit N°",i,"/",nrow(modeled_OMU)," ------\n"))
  
  ### Load unit occurrence points
  load(paste0("./input_data/Species_data/",res,"/Spatial_Points_Objects/occurrences_", unit,".RData")) # Load Spatial object with occurrences
  
  # Load projection map
  OMU_full_proj <- readRDS(file = paste0("./models/",res,"/",unit,"/proj_Current/proj_Current_",unit,".rds"))
  
  # plot(OMU_full_proj)
  
  ### 1/ Compute the buffer distance as the 80% or 95% quantile the distance to the closest occurrence points among occurrence points for each OMU ####
  
  # To be seen as an estimation of the uncertainty surrounding the completeness of our survey.
  # Well surveyed OMU should have short distances between occurrences, and no isolated points. Thus the distribution map can be clipped tighly around points
  # Unsufficiently surveyed OMU will display isolated points, thus the buffer must be wide to encompass the uncertainty around the distribution limits of this OMU
  # Use 95% quantile instead of maximum distance to avoid single isolated point to drive the size of the buffer. Allow 5% of outliers
  
  Coord_unit <- unit.points@coords
  
  Occ.dist = geosphere::distm(x = Coord_unit)/1000 # Geometric distance on WGS84 ellipsoid, in km
  
  # Remove diagonal (distance to thmeselves = 0 km)
  diag(Occ.dist) <- NA
  # Compute distance to closest other occurrence point
  Occ.dist_min <- apply(X = Occ.dist, MARGIN = 2, FUN = min, na.rm = T)
  
  # Get the 80% and 95% quantile and max.distance
  Uncertainty.80 <- round(quantile(x = Occ.dist_min, probs = 0.8),1)
  Uncertainty.95 <- round(quantile(x = Occ.dist_min, probs = 0.95),1)
  Max.dist <- round(max(Occ.dist_min),1)
  
  # Save in list.models
  list.models$Uncertainty.80[list.models$Tag.model == unit] <- Uncertainty.80
  list.models$Uncertainty.95[list.models$Tag.model == unit] <- Uncertainty.95
  list.models$Max.dist[list.models$Tag.model == unit] <- Max.dist
  
  
  ### 2/ Generate alpha-hull buffer ####
  
  # Temporary projection of the spatial object of occurrences in Mollweide
  proj.points <-spTransform(unit.points, CRSobj = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
  
  tryCatch(
    {
      Hull <- ahull_to_SPLDF(ahull(proj.points@coords, alpha = 1000)) # Try an alpha-hull first
    },
    error = function(e) # If fails, do an alpha shape instead
    {
      cat(paste0("Alpha-hull failed for ", unit, ", use alpha-shape instead \n"))
      cat("ERROR :",conditionMessage(e), "\n")   # Display the error message but do not stop the function
      Hull <<- ashape_to_SPLDF(ashape(proj.points@coords, alpha = 1000))
    })
  
  plot(Hull)
  
  # Hull <- ahull(proj.points@coords, alpha = 1000000) # alpha diameter = 1000km
  # Hull <- ahull_to_SPLDF(x = Hull)
  
  # For cases that cannot handle alpha-hull, use alpha-shape
  # Hull <- ashape(proj.points@coords, alpha = 1000000)
  # Hull <- ashape_to_SPLDF(x = Hull)
  
  # Save crs of the newly created alpha-hull
  Hull@proj4string <- proj.points@proj4string
  
  # Extend the shape with a buffer based on the quality of presence sampling (Uncertainty.80)
  # As a minimum, use the same buffer distance used to draw PsA = 111.32 km
  Buffer_80 <- gBuffer(spgeom = Hull, width = max(111.32, Uncertainty.80))

  # Remove holes
  Buffer_80_unholed <- spatialEco::remove.holes(Buffer_80)

  # Retransposition in WGS84
  Buffer_80.WGS84 <- spTransform(Buffer_80_unholed, CRSobj = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

  # Rasterization
  Buffer_80.rasterized <- rasterize(x = Buffer_80.WGS84, y = continent_mask, field = 1)
  # plot(Buffer_80.rasterized, add = T)

  # Save final buffers
  save(Buffer_80.rasterized, file = paste0("./outputs/By_OMU/",unit,"/Buffer_80_5.RData"), version = "2")
  saveRDS(Buffer_80.rasterized, file = paste0("./outputs/By_OMU/",unit,"/Buffer_80_5.rds"), version = "2")

  ### 3/ Clip the maps ####
  
  ## 3.1/ Clip with the buffers ####
  
  # Load buffers
  Buffer_80.rasterized <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Buffer_80_5.rds"))

  # Apply buffers
  OMU_map_cropped_80 <- raster::mask(OMU_full_proj, Buffer_80.rasterized)

  # plot(OMU_map_cropped_80)
  # plot(unit.points, pch = 16, add = T)
  
  ## 3.2/ Check within which Andean slope lies the occurrences, and crop out the slope where no occurrences have been recorded ####
  
  # Checkif presence points fall into the Andes slopes (used projected object in Mercator)
  test_East <- gIntersection(spgeom1 = Andes_East, spgeom2 = proj.points, checkValidity = 2)
  test_West <- gIntersection(spgeom1 = Andes_West, spgeom2 = proj.points, checkValidity = 2)
  
  # Clip out the slope where no occurrences have been recorded
  if (length(test_East)*length(test_West) == 0) { # Case with occurrences only on one side
    
    if (length(test_West)==0) { # Case with occurrences on the East side
      
      OMU_map_cropped_80 <- mask(OMU_map_cropped_80, Andes_East.rasterized)
      
      list.models$Andes_side[which(list.models$Tag.model == unit)] <- "East"
      
    } else { # Case with occurrences on the West side
      
      OMU_map_cropped_80 <- mask(OMU_map_cropped_80, Andes_West.rasterized)
      
      list.models$Andes_side[which(list.models$Tag.model == unit)] <- "West"
      
    }
    
  } else { # Case with occurrences on both sides = no additionnal cropping
    
    list.models$Andes_side[which(list.models$Tag.model == unit)] <- "Both"
  }
  
  ### Add 0 value for terrestrial pixels to be able to distinguish from oceans
  
  OMU_map_cropped_80_borders <- continent_mask*0
  OMU_map_cropped_80_borders@data@values[!is.na(OMU_map_cropped_80@data@values)] <- OMU_map_cropped_80@data@values[!is.na(OMU_map_cropped_80@data@values)]
  # plot(OMU_map_cropped_80_borders)

  # Save individually
  saveRDS(OMU_map_cropped_80_borders, file = paste0("./outputs/By_OMU/",unit,"/OMU_map_5.rds"))
    
  # Save list.models
  save(list.models, file = "./input_data/list.models_5min.RData")
  
  cat(paste0("\n", Sys.time()," ------ Done for ", unit, " = Unit N°",i,"/",nrow(modeled_OMU)," ------\n"))
  
} 

# Save Summary table for OMU models
save(list.models, file = "./input_data/list.models_5min.RData")



##### Script 12: Merging species #####

# Create maps at species level by computing "probability" of presence of any of the OMU of the species


##### 

# Inputs:
   # Stack of continuous maps for each OMU, per species from script 11

# Outputs:
   # Species continuous map
   # Stack of all OMU continuous maps

#####


### 1/ Prepare stuff ####

# Clean environment
rm(list = ls())

library(raster)

# Choose the resolution
res <-  "5"
# res <- "15"

# Load environmental stack to use as reference for extent and CRS and generate a mask for continental borders
library(raster)
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))
continent_mask <- envData[[1]]
continent_mask <- calc(continent_mask, fun = function(x, ...) {x*0}, na.rm = F)

# Load Summary table for OMU/unit and for Species
load(file = paste0("./input_data/list.sp_5min.RData"))
load(file = paste0("./input_data/list.models_5min.RData"))

# Make a list of modeled and non-modeled units (occ.unit = with only occurrences points)
# modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]
modeled_OMU <- list.models[!(list.models$initial_model_type == "rasterized"), ]
rasterized_OMU <- list.models[(list.models$initial_model_type == "rasterized"), ]

# Function to compute probability of presence with multiple OMUs
aggreg_prob = function(x, na.rm) { 
  y <- 1-prod(1-x) # Probability of presence of species = probability of presence of at least one OMU = opposite of probability of absence of all OMU
  return(y) # Output
}

list.sp$one.Binaries_5 <- F # To record which species include at least one binaries maps
list.sp$all.Binaries_5 <- F # To record which species include only binaries maps

### 2/ Loop by species to build species maps ####
for (i in 1:nrow(list.sp)) 
{
  # i <- 23
  
  # Load a random layer to initiate the stacks (to remove later)
  sp_stack_5min <- stack(continent_mask) 
  
  # Load Species name
  sp <- as.character(list.sp$Sp_full[i]) 
  
  # Create directory to store outputs per species
  if(!dir.exists(paste0("./outputs/By_species/",sp,"/5"))) { # Test if the folder exists already or not
    dir.create(paste0("./outputs/By_species/",sp,"/5"), recursive = T) # Create folder if absent
  }
  
  # Create directory to store maps per species
  if(!dir.exists(paste0("./maps/By_species/",sp,"/5"))) { # Test if the folder exists already or not
    dir.create(paste0("./maps/By_species/",sp,"/5"), recursive = T) # Create folder if absent
  }
  
  OMU_list <- modeled_OMU[modeled_OMU$Sp_full == sp, ]
  if (nrow(OMU_list) > 0) {      # If at least one OMU was modeled
    
    for (j in 1:nrow(OMU_list))  # For each OMU
    { 
      unit <-  as.character(OMU_list$Tag.model[j]) # Load the unit name                  
      
      # Load the continuous map
      OMU_map <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/OMU_map_5.rds"))
      # Convert into proba
      OMU_map <- OMU_map/1000
      
      # Stack them if there are several
      sp_stack_5min <- addLayer(sp_stack_5min, OMU_map)
      
    }
  }
  
  Binaries <- rasterized_OMU[rasterized_OMU$Sp_full == sp,]
  if (nrow(Binaries) > 0) {  # If at least one OMU was rasterized
    
    list.sp$one.Binaries_5[list.sp$Sp_full == sp] <- T
    
    for (j in 1:nrow(Binaries)) {  # For each  OMU
      
      unit <-  as.character(Binaries$Tag.model[j]) # Load the unit name 
      
      # Create directory to store outputs per OMU
      if(!dir.exists(paste0("./outputs/By_OMU/",unit,"/5"))) { # Test if the folder exists already or not
        dir.create(paste0("./outputs/By_OMU/",unit,"/5"), recursive = T) # Create folder if absent
      }
      
      # Create directory to store maps per OMU
      if(!dir.exists(paste0("./maps/By_OMU/",unit,"/5"))) { # Test if the folder exists already or not
        dir.create(paste0("./maps/By_OMU/",unit,"/5"), recursive = T) # Create folder if absent
      }
      
      # Load the rasterized occurrence map
      Raster_Occ <- readRDS(file = paste0("./input_data/Species_data/",res,"/Binary_Rasters_from_Occurrences/",unit,".rds"))

      # Add the continental borders as 0 values for empty pixels (oceans are NA)
      rasterized_map <- continent_mask
      rasterized_map[!is.na(Raster_Occ[])] <- 1

      save(rasterized_map, file = paste0("./outputs/By_OMU/",unit,"/5/rasterized_map_",unit,".RData"))
      saveRDS(rasterized_map, file = paste0("./outputs/By_OMU/",unit,"/5/rasterized_map_",unit,".rds"))

      # Load the rasterized map, with continental borders
      rasterized_map <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/5/rasterized_map_",unit,".rds"))
      
      # Stack  them if there are several
      sp_stack_5min <- addLayer(sp_stack_5min, rasterized_map)
      
    }
  }
  
  # Note case with only rasterized OMU
  if ((nrow(Binaries) > 0) & (nrow(OMU_list) == 0)) {
    list.sp$all.Binaries_5[list.sp$Sp_full == sp] <- T
  }
  
  # Drop the useless first layer used to initiate the stack
  sp_stack_5min <- dropLayer(sp_stack_5min, i = 1)

  # Compute probability map per Species
  
  sp_prob_map <- calc(sp_stack_5min, fun = aggreg_prob)

  plot(sp_prob_map)
  
  # Make sure to save the data, not just the link to the temp file!
  sp_prob_map <- readAll(sp_prob_map)
  
  save(sp_prob_map, file = paste0("./outputs/By_species/",sp,"/5/sp_prob_map_",sp,".RData"), version = "2")
  saveRDS(sp_prob_map, file = paste0("./outputs/By_species/",sp,"/5/sp_prob_map_",sp,".rds"), version = "2")

  pdf(file = paste0("./maps/By_species/",sp,"/5/sp_prob_map_",sp,".pdf"), height = 6, width = 7)
  plot(sp_prob_map, main = sp)
  dev.off()
  
  save(list.sp, file = paste0("./input_data/list.sp_5min.RData"))
  
  cat(paste0("\n", Sys.time()," ------ Done for ", sp, " = Species N°",i,"/",nrow(list.sp)," ------\n"))
  
}

# save(list.sp, file = paste0(internal.wd,"/list.sp_5min.RData"))

sum(list.sp$one.Binaries_5) # 157 species among the 388 have at least one binary map used
sum(list.sp$all.Binaries_5) # 54 species among the 388 are just modeled under binary maps



##### 3/ Generate stack of all OMUs outputs for all options ####


# Load a random layer to initiate the stacks (to remove later)
All_OMU_stack_5min <- stack(continent_mask)

### Loop by OMU
for (i in 1:nrow(list.models)) 
{
  # i <- 4
  
  # Load OMU name
  unit <- as.character(list.models$Tag.model[i]) 
  
  if (list.models$initial_model_type[i] != "rasterized") # For modeled OMUs
  { 
    # Load the continuous map
    OMU_map <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/OMU_map_5.rds"))
    # Convert into proba
    OMU_map <- OMU_map/1000
    # Stack them all
    All_OMU_stack_5min <- addLayer(All_OMU_stack_5min, OMU_map)
    
  } else  # For rasterized OMU
  { 
    # Load the rasterized map, with continental borders, and stack them all
    rasterized_map <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/5/rasterized_map_",unit,".rds"))
    All_OMU_stack_5min <- addLayer(All_OMU_stack_5min, rasterized_map)
    
  }
  
  if (i %% 10 == 0)
  {print(i)}
  
}

# Drop the useless first layer used to initiate the stack
All_OMU_stack_5min <- dropLayer(All_OMU_stack_5min, i = 1)

# Rename layers with OMU names
names(All_OMU_stack_5min) <- as.character(list.models$Tag.model)

plot(All_OMU_stack_5min)

save(All_OMU_stack_5min, file = paste0("./outputs/Indices_stacks/5/All_OMU_stack_5min.RData"))
saveRDS(All_OMU_stack_5min, file = paste0("./outputs/Indices_stacks/5/All_OMU_stack_5min.rds"))


##### 4/ Generate stack of all species #####

All_sp_proba_stack_5min <- stack(continent_mask) # 1e temp layer used to initiate stack, to remove afterwards
# i <- 1

for (i in 1:nrow(list.sp)) # By sp
{ 
  # Load sp name
  sp <- as.character(list.sp$Sp_full[i]) 
  
  # Load Sp continuous probability map 
  sp_prob_map <- readRDS(file = paste0("./outputs/By_species/",sp,"/5/sp_prob_map_",sp,".rds"))

  # Build stack
  All_sp_proba_stack_5min <- addLayer(All_sp_proba_stack_5min, sp_prob_map)

  if (i %% 10 == 0) {print(i)}
}

# Drop the useless first layer used to initiate the stack
All_sp_proba_stack_5min <- dropLayer(All_sp_proba_stack_5min, i = 1)

# Name layers with sp names
names(All_sp_proba_stack_5min) <- list.sp$Sp_full

# nlayers(All_sp_proba_stack_5min) # 388 species in the final stack

# Ensure data is store in the environment
All_sp_proba_stack_5min <- readAll(All_sp_proba_stack_5min)

# Save stacks
save(All_sp_proba_stack_5min, file = paste0("./outputs/Indices_stacks/5/All_sp_proba_stack_5min.RData"))
saveRDS(All_sp_proba_stack_5min, file = paste0("./outputs/Indices_stacks/5/All_sp_proba_stack_5min.rds"))


##### Script 13: Merging Mimicry rings #####

# Create map at mimicry ring level by computing "probability" of presence of any of the OMU of each mimicry ring

##### 

# Inputs:
    # Stack of maps for each OMU, per mimicry ring from script 11

# Outputs:
   # Mimicry ring probability continous map => aggregated with aggreg_proba function
   # Mimicry ring richness map => aggregated with sum

#####


### Prepare stuff

# Clean environment
rm(list = ls())

library(raster)

# Choose the resolution
res <-  "5"
# res <- "15"

# Load environmental stack to use as reference for extent and CRS and generate a mask for continental borders
library(raster)
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))
continent_mask <- envData[[1]]
continent_mask <- calc(continent_mask, fun = function(x, ...) {x*0}, na.rm = F)
# save(continent_mask, file = paste0("./input_data/Env_data/continent_mask_", res, ".RData"))
# saveRDS(continent_mask, file = paste0("./input_data/Env_data/continent_mask_", res, ".rds"))


# Load Summary table for OMU/unit and for Species
load(file = paste0("./input_data/list.sp_5min.RData"))
load(file = paste0("./input_data/list.models_5min.RData"))

# Mimicry list
mimicry.list <- as.character(unique(list.models$Mimicry.model)) # 44 Mimicry rings

# Make a list of modeled and non-modeled units (occ.unit = with only occurrences points)
# modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]
modeled_OMU <- list.models[!(list.models$initial_model_type == "rasterized"), ]
rasterized_OMU <- list.models[(list.models$initial_model_type == "rasterized"), ]

# Function to compute probability of presence with multiple OMUs
aggreg_prob = function(x, na.rm) { 
  y <- 1-prod(1-x) # Probability of presence of species = probability of presence of at least one OMU = opposite of probability of absence of all OMU
  return(y) # Output
}


### Loop by mimicry ring
for (i in 1:length(mimicry.list)) 
{
  # i <- 10
  
  # Load a random layer to initiate the stacks (to remove later)
  ring_stack_5min <- stack(continent_mask) 
  
  # Load Mimicry ring name
  ring <- as.character(mimicry.list[i])
  
  ### 1/ Make stacks of OMU maps ####
  
  OMU_list <- modeled_OMU[modeled_OMU$Mimicry.model == ring, ]
  if (nrow(OMU_list) > 0) {      # If at least one OMU was modeled
    
    for (j in 1:nrow(OMU_list))  # For each OMU
    { 
      unit <-  as.character(OMU_list$Tag.model[j]) # Load the unit name                  
      
      # Load the continuous map
      OMU_map <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/OMU_map_5.rds"))
      # Convert to proba
      OMU_map <- OMU_map/1000
      
      # Stack them if there are several
      ring_stack_5min <- addLayer(ring_stack_5min, OMU_map)
      
    }
  }
  
  Binaries <- rasterized_OMU[rasterized_OMU$Mimicry.model == ring, ]
  if (nrow(Binaries) > 0) {  # If at least one OMU was rasterized
    
    for (j in 1:nrow(Binaries)) {  # For each  OMU
      
      unit <-  as.character(Binaries$Tag.model[j]) # Load the unit name 
      
      # Load the rasterized map
      rasterized_map <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/5/rasterized_map_",unit,".rds"))
      
      # Stack  them if there are several
      ring_stack_5min <- addLayer(ring_stack_5min, rasterized_map)
      
    }
  }
  
  # Drop the useless first layer used to initiate the stack
  ring_stack_5min <- dropLayer(ring_stack_5min, i = 1)

  ### 2/ Compute probability map per ring ####
  
  ring_prob_map <- calc(ring_stack_5min, fun = aggreg_prob)

  plot(ring_prob_map, main = paste0(ring,"\nResolution: 5min"))
  
  # Make sure to save the data, not just the link to the temp file!
  ring_prob_map <- readAll(ring_prob_map)

  save(ring_prob_map, file = paste0("./outputs/Mimicry_rings_proba/5/ring_prob_map_",ring,".RData"))
  saveRDS(ring_prob_map, file = paste0("./outputs/Mimicry_rings_proba/5/ring_prob_map_",ring,".rds"))

  ### 3/ Compute richness map per ring (nb of OMU) ####
  
  ring_richness_map <- calc(ring_stack_5min, fun = sum)

  plot(ring_richness_map, main = paste0(ring,"\nResolution: 5min"))
  
  # Make sure to save the data, not just the link to the temp file!
  ring_richness_map <- readAll(ring_richness_map)

  save(ring_richness_map, file = paste0("./outputs/Mimicry_ring_richness/5/ring_richness_map_",ring,".RData"))
  saveRDS(ring_richness_map, file = paste0("./outputs/Mimicry_ring_richness/5/ring_richness_map_",ring,".rds"))

  print(i)
  
}


##### 4/ Generate stack of all mimicry rings proba and richness #####

# Mimicry list
mimicry.list <- as.character(unique(list.models$Mimicry.model)) # 44 Mimicry rings

All_ring_proba_stack_5min <- All_ring_richness_stack_5min <- stack(continent_mask) # 1e temp layer used to initiate stack, to remove afterwards
# i <- 1

for (i in 1:length(mimicry.list)) # By ring
{ 
  # Load ring name
  ring <- as.character(mimicry.list[i]) 
  
  # Load ring continuous probability map 
  ring_prob_map <- readRDS(file = paste0("./outputs/Mimicry_rings_proba/5/ring_prob_map_",ring,".rds"))
  
  # Load ring richness map
  ring_richness_map <- readRDS(file = paste0("./outputs/Mimicry_ring_richness/5/ring_richness_map_",ring,".rds"))
  
  # Build stacks
  All_ring_proba_stack_5min <- addLayer(All_ring_proba_stack_5min, ring_prob_map)
  All_ring_richness_stack_5min <- addLayer(All_ring_richness_stack_5min, ring_richness_map)
  
  if (i %% 10 == 0) {print(i)}
}

# Drop the useless first layer used to initiate the stack
All_ring_proba_stack_5min <- dropLayer(All_ring_proba_stack_5min, i = 1)
All_ring_richness_stack_5min <- dropLayer(All_ring_richness_stack_5min, i = 1)

# Name layers with sp names
names(All_ring_proba_stack_5min) <- mimicry.list
names(All_ring_richness_stack_5min) <- mimicry.list

# Ensure data is store in the environment
All_ring_proba_stack_5min <- readAll(All_ring_proba_stack_5min)
All_ring_richness_stack_5min <- readAll(All_ring_richness_stack_5min)

# Save stacks
save(All_ring_proba_stack_5min, file = paste0("./outputs/Indices_stacks/5/All_ring_proba_stack_5min.RData"))
saveRDS(All_ring_proba_stack_5min, file = paste0("./outputs/Indices_stacks/5/All_ring_proba_stack_5min.rds"))
save(All_ring_richness_stack_5min, file = paste0("./outputs/Indices_stacks/5/All_ring_richness_stack_5min.RData"))
saveRDS(All_ring_richness_stack_5min, file = paste0("./outputs/Indices_stacks/5/All_ring_richness_stack_5min.rds"))




##### Script 14a: Compute indices maps #####

# Compute indices map from stacks

###

# Inputs:
   # Stack of species 'probability' of presence from script 12
   # Stack of mimicry ring 'probability' of presence from script 13
   # Stack of mimicry ring richness of presence from script 13
   # Phylogeny

# Outputs:
   # Biodiversity index maps.

###

### Prepare stuff

# Clean environment
rm(list = ls())

# Load stacks of species and rings
All_sp_proba_stack_5min <- readRDS(file = paste0("./outputs/Indices_stacks/5/All_sp_proba_stack_5min.rds"))
All_ring_proba_stack_5min <- readRDS(file = paste0("./outputs/Indices_stacks/5/All_ring_proba_stack_5min.rds"))
All_ring_richness_stack_5min <- readRDS(file = paste0("./outputs/Indices_stacks/5/All_ring_richness_stack_5min.rds"))

# Function to add continental null values to a raster
add_continental_null_values <- function(x, res)
{
  continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_", res, ".rds"))

  y <- continent_mask  # Create final new raster from continental mask
  y[!is.na(x[])] <- x[!is.na(x[])]  # Add initial raster values
  
  return(y)
}




##### 1/ Species richness #####

plot(All_sp_proba_stack_5min)

# Compute index
sp_richness_5min <- calc(All_sp_proba_stack_5min, fun = sum)

# Save
save(sp_richness_5min, file = paste0("./outputs/Indices_maps/5/sp_richness_5min.RData"))
saveRDS(sp_richness_5min, file = paste0("./outputs/Indices_maps/5/sp_richness_5min.rds"))

plot(sp_richness_5min, col = pal_bl_red_Mannion)

# Get ithomiini range mask to use function for applying contrast
ithomiini_range_5min <- sp_richness_5min > 0
ithomiini_range_5min@data@values[ithomiini_range_5min@data@values == 0] <- NA
plot(ithomiini_range_5min)

saveRDS(ithomiini_range_5min, file = paste0("./input_data/Map_stuff/ithomiini_range_5min.rds"))


##### 2/ Mimicry richness #####

# Compute index
ring_richness_5min <- readAll(calc(All_ring_proba_stack_5min, fun = sum))

plot(ring_richness_5min, col = pal_bl_red_Mannion)

# Save
save(ring_richness_5min, file = paste0("./outputs/Indices_maps/5/ring_richness_5min.RData"))
saveRDS(ring_richness_5min, file = paste0("./outputs/Indices_maps/5/ring_richness_5min.rds"))


##### 3/ Continuous Species Rarity = Range-size Weighted Species Richness #####

# Function to compute index
compute_geographic_rarity <- function (proba_stack)
{
  library(Rarity)
  
  # Compute expected richness
  richness <- calc(x = proba_stack, fun = sum)
  
  # Extract gridcell area in km²
  grid_area_raster <- raster::area(x = richness, na.rm = F, weights = F)
  area_cell <- round(mean(grid_area_raster@data@values, na.rm = T), 2)
  
  # Compute range size for for weights based on geographic ranges
  range_size <- NA
  for (i in 1:nlayers(proba_stack)) {
    # Compute the estimated number of pixels occupied by the unit as the sum of probabilities and multiply by cell area/1000 to have values in 10^3 km²
    range_size[i] <- sum(proba_stack[[i]]@data@values, na.rm = T)*area_cell/1000
  }
  
  # Generate incidence/abundance matrix of community assemblage (rows = units, col = communities)
  proba_brick <- readAll(brick(proba_stack))
  assemblage <- t(tidyr::drop_na(as.data.frame(proba_brick@data@values)))
  assemblage <- assemblage[, colSums(assemblage) >= 1]
  
  # Rarity weights based on an inverse exponential function with inflexion point calibrated as the rarity threshold as such that communities host 25% of rare species in average
  Leroy_weights_df <- rWeights(occData = range_size, wMethods = "W", rCutoff = "Leroy",
                               normalised = T, rounding = 5, assemblages = assemblage)
  Leroy_weights <- Leroy_weights_df$W
  
  # Apply weights
  weighted_proba_stack <- proba_stack # Generate the stack to fill
  for (i in 1:nlayers(weighted_proba_stack)){
    # Multiply species probability of presence by species rarity indices
    weighted_proba_stack@layers[[i]]@data@values <- proba_stack@layers[[i]]@data@values * Leroy_weights[i]
  }
  
  # Sum = Richness weighted by rarity based on Range-size
  geo_rarity_sum <- readAll(calc(weighted_proba_stack, fun = sum)) 
  # Divided by local richness = Mean rarity indices in the community (weighted by the probability of presence of each unit)
  geo_rarity_mean <- geo_rarity_sum/richness
  
  # Add continental values
  geo_rarity_mean <- add_continental_null_values(x = geo_rarity_mean, res = "5")
  
  return(geo_rarity_mean)
}

# Compute index
sp_geo_rarity_mean_5min <- readAll(compute_geographic_rarity(All_sp_proba_stack_5min))

plot(sp_geo_rarity_mean_5min, col = pal_bl_red_Mannion)

# Need to add a few values to lower values to avoid being colored as null outer-range values
sp_mean_geo_rarity_Mollweide[(sp_mean_geo_rarity_Mollweide[] > 0) & (sp_mean_geo_rarity_Mollweide[] <= 0.005)] <- 0.005


# Save
save(sp_geo_rarity_mean_5min, file = paste0("./outputs/Indices_maps/5/sp_geo_rarity_mean_5min.RData"))
saveRDS(sp_geo_rarity_mean_5min, file = paste0("./outputs/Indices_maps/5/sp_geo_rarity_mean_5min.rds"))


##### 4/ Continuous Mimicry rarity = Range size weighted mimicry richness #####

# Compute index
ring_geo_rarity_mean_5min <- readAll(compute_geographic_rarity(All_ring_proba_stack_5min))

plot(ring_geo_rarity_mean_5min, col = pal_bl_red_Mannion)

# Save
save(ring_geo_rarity_mean_5min, file = paste0("./outputs/Indices_maps/5/ring_geo_rarity_mean_5min.RData"))
saveRDS(ring_geo_rarity_mean_5min, file = paste0("./outputs/Indices_maps/5/ring_geo_rarity_mean_5min.rds"))


##### 5/ Weighted Mean ring Size ####

# Function to compute weighted mean ring size
compute_weighted_mean_ring_size <- function (ring_proba_stack, ring_rich_stack)
{
  # Convert to Raster brick
  ring_proba_brick <- readAll(brick(ring_proba_stack))
  ring_rich_brick <- readAll(brick(ring_rich_stack))
  
  weighted_mean_ring_size <- vector()
  for (i in 1:nrow(ring_proba_brick[]))  # Loop between communities/pixels
  {
    # i <- 550000
    # i <- 650000
    
    rich_vector <- ring_rich_brick@data@values[i,]   # Extract ring richness
    proba_vector <- ring_proba_brick@data@values[i,] # Extract ring proba of presence
    
    # Compute mean weighted by proba. of presence
    weighted_mean_ring_size[i] <- weighted.mean(x = rich_vector, w = proba_vector)
    
    if(i %% 100000 == 0)
    {
      cat(paste0("Iteration n°", i, "/", nrow(ring_proba_brick[]), "\n"))
    }
  }
  
  # Add value to a raster with continental borders
  weighted_mean_ring_size <- add_continental_null_values(x = weighted_mean_ring_size, res = "5")
  
  return(weighted_mean_ring_size)
}

# Compute index
weighted_mean_ring_size_5min <- compute_weighted_mean_ring_size(ring_proba_stack = All_ring_proba_stack_5min, ring_rich_stack = All_ring_richness_stack_5min)

plot(weighted_mean_ring_size_5min, col = pal_bl_red_Mannion)

# Save
save(weighted_mean_ring_size_5min, file = paste0("./outputs/Indices_maps/5/weighted_mean_ring_size_5min.RData"))
saveRDS(weighted_mean_ring_size_5min, file = paste0("./outputs/Indices_maps/5/weighted_mean_ring_size_5min.rds"))


##### 6/ Faith's Phylogenetic Diversity #####

library(ape)
# library(picante)
library(geiger)
# library(treespace)

### Create function to aggregate probabilities to higher hierachical level (aggregate pixel, or go up on a phylogenetic tree)
aggreg_prob = function(x, na.rm) { 
  y <- 1-prod(1-x) 
  return(y) # Output
}

### Function to compute Faith's phylogenetic diversity
compute_faith_phylogenetic_diversity <- function (sp_proba_stack, phylo)
{
  # Extract only the species included in the phylogeny from the stacks of sp probas
  sp_proba_stack_phylo <- sp_proba_stack[[phylo$tip.label]]
  
  # Generate matrix to store info on edges/branches
  branches = matrix(NA, nrow(phylo$edge), ncol = 4) 
  names(branches) <- c("Starting nod", "Ending nod", "Length", "Proba_presence")
  branches[, 1:2] = phylo$edge # Retrieve starting and ending node 
  branches[, 3] = round(phylo$edge.length, 4) # Retrieve edge length
  
  sp_proba_brick <- readAll(brick(sp_proba_stack_phylo))

  # Make a loop per community
  PD_data <- NA
  for (k in 1:nrow(sp_proba_brick[]))
  {
    # k <- 50000
    # k <- 593788
    
    # Extract the proba for the k community
    sp_proba_com <- sp_proba_brick@data@values[k,] 
    
    # Compute PD only if no NA is present in the community
    if (any(is.na(sp_proba_com)))
    {
      PD_data[k] <- NA # if NA present, PD = NA
      
    } else {
      
      # Compute probability of presence of each edge/branch (i.e., a species descending from this branch) in this community
      for (i in 1:nrow(branches)) 
      {
        leaves.node = geiger::tips(phylo, branches[i, 2]) # Retrieve the set of species descending from branch i
        index <- which(colnames(sp_proba_brick[]) %in% leaves.node) # Find position of the species in the stack/matrix
        prob_edge <- aggreg_prob(sp_proba_com[index]) # Compute proba of presence of the edge
        branches[i, 4] <- prob_edge # Store info
      }
      PD <- round(sum(branches[,3]*branches[,4]),4) # Compute PD as the weighted sum of the edge length (weighted by the probabilty of presence of this edge in the community)
      PD_data[k]  <- PD
    }
    
    # Show k every 1000 iterations and save a backup
    if (k %% 1000 == 0) {
      cat(paste0(Sys.time(), " - ", k," on ",nrow(sp_proba_brick[]),"\n"))
      save(PD_data, file = "./PD_data_backup.RData")
    }
  }
  
  # Write PD values in a raster
  PD_raster <- sp_proba_stack_phylo[[1]]*0
  PD_raster@data@values <- PD_data
  # plot(PD_raster)
  
  # Repair issue with max value
  PD_raster@data@max <- max(PD_raster[], na.rm = T)
  
  return(PD_raster)
  
}


# Load the phylogeny
phylo.Ithomiini <- readRDS(file = "./input_data/Phylogenies/Final_phylogeny.rds")

# Compute index
Faith_PD_5min <- compute_faith_phylogenetic_diversity(sp_proba_stack = All_sp_proba_stack_5min, phylo = phylo.Ithomiini)

plot(Faith_PD_5min, col = pal_bl_red_Mannion)

# Save
save(Faith_PD_5min, file = paste0("./outputs/Indices_maps/5/Faith_PD_5min.RData"))
saveRDS(Faith_PD_5min, file = paste0("./outputs/Indices_maps/5/Faith_PD_5min.rds"))



##### Script 14c: Plot indices maps #####

### Plot nice maps for 5min resolution

### Version with 6 Indices (Figure S5.14)

# A/ Species richness
# B/ Mean Species geographic rarity
# C/ Faith's Phylogenetic Diversity

# D/ Mimicry richness
# E/ Mean mimicry geographic rarity
# F/ Weighted mean ring size

#####

# Inputs
   # Index maps from Script 14a

# Outputs
   # Clean maps with Mollweide projection of the 6 main indices (Figure S5.14)

#####

# Clean environment
rm(list = ls())

# Load library
library(raster)
library(rangeBuilder)
library(gplots)

##### 1/ Load useful stuff for map #####

# Color palette for plot
pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")

# Choose resolution
res <- "5"
# res <- "15"

# Load mask for continent borders
continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_", res, ".rds"))
# crop_mask_shp <- rasterToPolygons(x = continent_mask, dissolve = T, digits = 4)
# save(crop_mask_shp, file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".RData"), version = "2")
# saveRDS(crop_mask_shp, file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".rds"), version = "2")
crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".rds"))

# Load shp for plot border, and grid
grid_Mollweide_out <- readRDS(file = "./input_data/Map_stuff/grid_Mollweide_out.rds")
large_bg_mask_Mollweide <- readRDS(file = "./input_data/Map_stuff/large_bg_mask_Mollweide.rds")
bbox_sp_Mollweide <- readRDS(file = "./input_data/Map_stuff/bbox_sp_Mollweide.rds")

load(file = "./input_data/Map_stuff/country_borders.RData")
country_borders <- as(country_borders, "Spatial")
# plot(country_borders)

### Load maps 
sp_richness <- readRDS(file = paste0("./outputs/Indices_maps/5/sp_richness_5min.rds"))
sp_mean_geo_rarity <- readRDS(file = "./outputs/Indices_maps/5/sp_geo_rarity_mean_5min.rds")   # Leroy's weighting
Faith_PD <- readRDS(file = paste0("./outputs/Indices_maps/5/Faith_PD_5min.rds"))
ring_richness <- readRDS(file = paste0("./outputs/Indices_maps/5/ring_richness_5min.rds"))
ring_mean_geo_rarity <- readRDS(file = "./outputs/Indices_maps/5/ring_geo_rarity_mean_5min.rds")   # Leroy's weighting
weighted_mean_ring_size <- readRDS(file = paste0("./outputs/Indices_maps/5/weighted_mean_ring_size_5min.rds"))


##### 2/ Contrast indices when needed #####

# 2.1/ Contrasting raster function
contrasting_raster_5min <- function(x, zmin, zmax)
{
  ithomiini_range_5min <- readRDS(file = paste0("./input_data/Map_stuff/ithomiini_range_5min.rds"))
  continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_5.rds"))
  
  x[x[] <= zmin] <- zmin  # Fix low values
  x[x[] >= zmax] <- zmax  # Fix high values
  
  x <- mask(x, mask = ithomiini_range_5min)  # Cut out values that are outside Ithomiini range
  
  y <- continent_mask + zmin  # Create final new raster from continental mask with baseline = zmin
  y[!is.na(x[])] <- x[!is.na(x[])]  # Add initial raster values
  
  return(y)
}

# 2.2/ Apply contrast to rarity indices

ithomiini_range_5min <- readRDS(file = paste0("./input_data/Map_stuff/ithomiini_range_5min.rds"))

hist(sp_mean_geo_rarity)

sp_mean_geo_rarity_contrasted <- contrasting_raster_5min(x = sp_mean_geo_rarity, zmin = 0, zmax = 0.5)
sp_mean_geo_rarity_contrasted[!is.na(ithomiini_range_5min[]) & (sp_mean_geo_rarity_contrasted[] < 0.005)] <- sp_mean_geo_rarity_contrasted[!is.na(ithomiini_range_5min[]) & (sp_mean_geo_rarity_contrasted[] < 0.005)] + 0.005

plot(sp_mean_geo_rarity, col = pal_bl_red_Mannion)
plot(sp_mean_geo_rarity_contrasted, col = pal_bl_red_Mannion)

hist(ring_mean_geo_rarity)

ring_mean_geo_rarity_contrasted <- contrasting_raster_5min(x = ring_mean_geo_rarity, zmin = 0, zmax = 0.75)
ring_mean_geo_rarity_contrasted[!is.na(ithomiini_range_5min[]) & (ring_mean_geo_rarity_contrasted[] < 0.005)] <- ring_mean_geo_rarity_contrasted[!is.na(ithomiini_range_5min[]) & (ring_mean_geo_rarity_contrasted[] < 0.005)] + 0.005

plot(ring_mean_geo_rarity, col = pal_bl_red_Mannion)
plot(ring_mean_geo_rarity_contrasted, col = pal_bl_red_Mannion)

##### 2/ Project into Mollweide #####

### 2.1/ Projection functions

# Project raster into Mollweide
Mollweide_projection <- function(x, name) # Raster to project
{
  new_map <- projectRaster(from = x, 
                           method = "bilinear", # Method for interpolation => "ngb" = nearest neighbor for qualitative (or discrete) variables . "bilinear" = for quantitative variables
                           crs = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs", # If you have the CRS arguments
                           alignOnly = F)
  
  # Generate new object with "_Mollweide" suffix in the global environment
  eval(call("<<-", as.name(paste0(name, "_Mollweide")), new_map))
}

# Project sp shape into Mollweide
Mollweide_shp_projection <-  function(x) # Shp to project
{
  x_name <- deparse(substitute(x)) # Get the name of the initial shp as a character string
  
  new_shp <- spTransform(x, CRSobj = "+proj=moll +lon_0=-75 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=km +no_defs")
  
  # Generate new object with "_Mollweide" suffix in the global environment
  eval(call("<<-", as.name(paste0(x_name, "_Mollweide")), new_shp))
}

# 2.2/ Project all maps

# Put all maps in a list

list_all_maps <- list(sp_richness, sp_mean_geo_rarity_contrasted, Faith_PD, ring_richness, ring_mean_geo_rarity_contrasted, weighted_mean_ring_size)
names(list_all_maps) <- c("sp_richness", "sp_mean_geo_rarity_contrasted", "Faith_PD", "ring_richness", "ring_mean_geo_rarity_contrasted", "weighted_mean_ring_size")

# Loop to generate projected maps
for (i in 1:length(list_all_maps))
{
  Mollweide_projection(x = list_all_maps[[i]], name = names(list_all_maps)[i])
}

# Project country_borders shp
Mollweide_shp_projection(country_borders)


##### 3/ Plotting functions #####



# Function to nicely plot indices in Mollweide
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


##### 4/ Final pdf for Figure S5.14 #####

pdf(file = paste0("./supplementaries/Figure_S5.14_Indices_maps_5min.pdf"), height = 8, width = 12)

internal_margins <- par()$mar
# par(mar = c(3.1, 3.5, 3.5, 2.1))
par(mar = c(3.1, 3.1, 2.7, 1.6))
par(mfrow = c(2, 3))

# A/ Species richness ####

map_indices_Mollweide(x = sp_richness_Mollweide,
                      main_title = "Species richness",
                      legend_title = "Species",
                      legend_title_x = -3650,
                      legend_title_y = 430,
                      legend_breaks = seq(0, 120, 20), 
                      facet_letter = "(a)")

# B/ Species mean geographic rarity ####

map_indices_Mollweide(x = sp_mean_geo_rarity_contrasted_Mollweide,
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

map_indices_Mollweide(x = ring_richness_Mollweide,
                      main_title = "Mimicry richness",
                      legend_title = "Mimicry\nrings",
                      legend_title_x = -3670,
                      legend_title_y = 670,
                      legend_breaks = seq(0, 25, 5), 
                      facet_letter = "(d)")

# E/ Mimicry ring mean geographic rarity ####

map_indices_Mollweide(x = ring_mean_geo_rarity_contrasted_Mollweide,
                      main_title = "Mimicry geographic rarity",
                      legend_title = "Rarity\nindex",
                      legend_title_x = -3650,
                      legend_title_y = 670,
                      legend_breaks = seq(0, 0.6, 0.2), 
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


