
##### Script 4: Generation of Cross-validation blocks for OMU with N ≥ 30 ######

# Generate a 3-fold CV scheme based on spatial-blocks to evaluate models 

# Only for OMU with N ≥ 30


# Inputs:
   
   # Summary table for OMU models (list.models)
   # Spatial object wit occurrences and PsA coordinates from Script 03
   # PA tables from script 03

# Outputs:

   # CV scheme for each OMU with N ≥ 30 as calib.lines in Biomod2 (calib_array)
   # PDF plots with repartition of occurrences and PsA between folds




# Remove environment
rm(list = ls())


### 1/ Set up the cluster ####

library(foreach)
library(doParallel)

# Detect the number of threads
numCores <- detectCores()
# registerDoParallel(numCores)  # To set the nb of threads to use when using %dopar%

cl <- makeCluster(numCores) # Generate the virtual cluster from the number of threads
registerDoParallel(cl)      # To set the virtual cluster when using %dopar%

# getDoParWorkers() # Pour vérifier le nombre de coeurs enregistrés

### 2/ Load stuff useful for all units (outside the loop) ####

### Load table to track model outputs
load(file = "./input_data/list.models.RData")
unit.list <- list.models$Tag.model
sp.list <- list.models$Sp_full

### Select Env data resolution 

# Chose resolution
# res <- "5"
res <- "15"

# Load environmental stack
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))

# Set seed to ensure repetability of random PsA draws
set.seed(158340)

# Choose the number of CV folds to separate the occurrences and PsA in between
nb.of.folds = 3

source("./functions/progcombine.R")

set.seed(1383)

##### 3/ Loop to generate cross validation blocks  ######

# Parallelized version

temp <- foreach (k = seq_along(unit.list), .combine = progcombine_rbind(nreps = length(unit.list)), .init = NULL) %dopar% {
  # temp <- foreach (k = 1:5) %dopar% {
  # k <- 1
  
  unit <- as.character(unit.list[k])
  sp <- sp.list[k]
  
  N.obs <- list.models$N.obs_15m_used[k]
  
  # Do the CV only for units with at least 30 occurrences
  if (N.obs < 30) {
    
    # Warning message for OMU with less than 30 obs
    cat("\n", as.character(Sys.time()), "----- No need for CV blocks generation for", unit, "= Unit N°",k,"\n")
    
  } else { # CV design for OMU with at least 30 obs
    
    cat("\n", as.character(Sys.time()), "----- CV blocks generation starts for", unit, "= Unit N°",k,"\n")
    
    # load libraries
    
    library(raster)
    library(sf)
    library(tmap)
    library(tmaptools)
    library(blockCV)
    
    # Load spatial points with PsA for this unit/OMU
    load(file = paste0("./input_data/Species_data/",res,"/Spatial_Points_Objects/occ_&_PsA_", unit,".RData"))
    
    # Transform spatial points to an sf object
    occ_PsA_data <- sf::st_as_sf(unit.points.PsA, coords = c("x", "y"), crs = raster::crs(envData))
    
    # Need to create as many folds as there are PsA sets
    # Load PA.table
    load(file = paste0("./input_data/Species_data/",res,"/PA.Tables/PA.table_", unit,".RData"))
    
    # 3.1/ Generate empty array to record calibration and validation lines (calib.lines in Biomod2) ####
    
    # NA = not use in this run. T = calibration/learning set. F = validation/test set.
    # Start with array full of NA.
    calib_array <- array(dim = c(nrow(PA.table),   # Number of occrrences and PsA
                                 nb.of.folds,      # Number of CV folds
                                 ncol(PA.table)),  # Number of PsA_set
                         dimnames = list(NULL,                              # No name for point rows
                                         paste0("_RUN", 1:nb.of.folds),     # Need to add "_" as in BioMod2 calib.lines
                                         colnames(PA.table))                # Keep PA set names from PA.table
    )
    
    for (pa in colnames(PA.table)) { # Generate CV table for each PA.run
      
      cur.occ_PsA_data <- occ_PsA_data[which(PA.table[, pa]), ]  # Get only lines for the current PsA run
      
      # Transform PsA into Absences (0) to get the spatialBlock function working and obtain balanced folds with as many Presence and PsA (seen as Absences)
      cur.occ_PsA_data$Presences_Pseudo.Absences[which(is.na(cur.occ_PsA_data$Presences_Pseudo.Absences))] <- 0
      
      # 3.2/ Run the fold repartition process ####
      cur.blocks <- spatialBlock(speciesData = cur.occ_PsA_data,           # Sf_object with data
                                 species = "Presences_Pseudo.Absences",    # Column with Presence (1), Absence (0), and PsA (NA)
                                 rasterLayer = envData[["Elevation"]],     
                                 theRange = 100000,                        # Block range set to 100km
                                 k = nb.of.folds,                          # nb of CV folds to generate
                                 selection = "random",                     # Assign folds to each blocks randomly
                                 iteration = 100,                          # Nb of iterations allowed to try to find the best sorting which balance Presence and PsA among folds
                                 showBlocks = F,                           # Don't plot automatically
                                 biomod2Format = TRUE)                     # Get an output in the Biomod2 format
      
      # Results : train_0 = Nb of PsA for calibration/learning
      #           train_1 = Nb of presences for calibration/learning
      #           test_0 = Nb of PsA for evaluation/learning
      #           test_1 = Nb of presences for evaluation/learning
      
      # Fill calib.array with the output of spatialBlock
      
      calib_array[which(PA.table[, pa]), , pa] <- cur.blocks$biomodTable # Only the lines for this PsA set model
      
      # 3.3/ Save a plot of this CV repartition ####
      
      shape_vector <- as.factor((-1*cur.occ_PsA_data$Presences_Pseudo.Absences) + 2)
      cur.occ_PsA_data <- cbind(cur.occ_PsA_data, foldID = as.factor(cur.blocks$foldID), Type = shape_vector)
      
      pdf(file = paste0("./input_data/Species_data/",res,"/CV_blocks/block_maps/block_map_", unit,"_", pa, ".pdf"), height = 8, width = 10)
      
      block_map <- tm_shape(envData[["Elevation"]]) + tm_raster(palette = terrain.colors(10), legend.show = F) + tm_layout(legend.position = c("left", "bottom")) +
        tm_shape(cur.occ_PsA_data) + tm_dots(shape = "Type", col = "foldID", palette = c("black", "dodgerblue", "red"), size = 0.4, border.col = "black", border.lwd = 2, legend.show = F) +
        tm_layout(title = paste0("Fold repartition for spatial CV of ", unit, " for ", pa)) +
        tm_add_legend(type = "fill", labels = c(" 0 to 1,000"," 1,000 to 2,000", " 2,000 to 3,000", " 3,000 to 4,000", " 4,000 to 5,000", " 5,000 to 6,000"), col = terrain.colors(6), border.col = NA, title = "Elevation") +
        tm_add_legend(type = "symbol", shape = 20, labels = c("1","2","3"), col = c("black", "dodgerblue", "red"), title = "CV fold", size = 0.8) +
        tm_add_legend(type = "symbol", shape = c(21, 22) , labels = c(" Presence"," Pseudo-absence"), col = "black", title = "Type", size = 0.5)
      
      print(block_map)
      
      dev.off()
      
      # Save folds infos for this PsA set
      saveRDS(cur.occ_PsA_data, file = paste0("./input_data/Species_data/",res,"/CV_blocks/CV_folds_data/CV_folds_data_", unit,"_",pa,".rds"))
      
    }
    
    # save CV infos for this OMU
    save(calib_array, file = paste0("./input_data/Species_data/",res,"/CV_blocks/calib_arrays/calib_array_", unit,".RData"))
    saveRDS(calib_array, file = paste0("./input_data/Species_data/",res,"/CV_blocks/calib_arrays/calib_array_", unit,".rds"))
    
    cat("\n", as.character(Sys.time()), "----- CV blocks generation for", unit, "= Unit N°",k,"- Done\n")
    
  }
}

# Close clusters generated with registerDoParallel(numCores)
# stopImplicitCluster()
# Close clusters generated with makeCluster 
stopCluster() 
