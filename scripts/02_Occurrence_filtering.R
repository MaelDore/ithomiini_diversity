

##### Script 2: Filtering of occurrences #####

# Apply raster grid to merge occurrences in the same grid cell (spatial filtering)

# Remove occurrences without environnemntal data associated to their location
    # But try to retrieve data in a small range (max 100km) to keep a maximum of occurrences
    # Rationale is taht many occurrences with no env. data are just coastal points falling just aside the pixels with valid env. data.

# Generate spatial objects with corrected coordinates to correspond to grid cell centroids

# Stack final environmental layer, corrected for occurrences in error when possible


###

# Inputs:

  # Occurrence points for each OMU (Ithomiini_final)
  # Global env. stack from script 01
  # Summary table to track model output (for each OMU) = list.models

# Outputs:

  # Final env. stack for each OMU
  # Rasterize map of occurrences for each OMU (Binary maps)
  # PDF map of occurrences, with errors and corrections
  # Info on Nb of occurrences, model type (complete, restricted, rasterized)
  # PDF plot comparing final nb of occurrences per OMU depending on resolution used

### 2 Versions of the script: Without parallelization, and with parallelization ###



##### 1/ Without parallelisation #####

# Remove environment
rm(list = ls())

# Load libraries
library(raster)

### Load table to track model outputs
load(file = "./input_data/list.models.RData")
unit.list <- list.models$Tag.model
sp.list <- list.models$Sp_full

# list.models$N.obs_5m <- list.models$N.obs_5m_used <- list.models$N.obs_15m <- list.models$N.obs_15m_used <- NA


### Load occurrence dataset
load(file = "./input_data/Databases/Ithomiini_final.RData")

### Select Env data resolution 

# Chose resolution
res <- "5"
res <- "15"

# Load environmental stack
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))


### Loop per OMU

for (k in seq_along(unit.list)) {
  # for (k in 350:369) {
  # k <- 25

  unit <- as.character(unit.list[k])
  cat(as.character(Sys.time()), "----- Start for", unit, "= Unit N°",k,"\n")
  unit.coords <- Ithomiini_final[(Ithomiini_final$Tag.model == unit), c("Longitude","Latitude")] 
  
  ### 1.1/ Rasterize occurrences ####
  unit.occ <- rasterize(x = unit.coords, # Occurrence coordinates
                        y = envData, # Env stack grid to apply
                        field = 1, # Value to store in the raster
                        fun = function(x, ...) mean(x)) # Apply mean to keep value = 1 in case of multiple occurrences in one pixel
  
  unit.env <- stack(unit.occ, envData) # Add the rasterized occurrence layer to the env stack
  
  # Count occurrences after rasterization and before correction for NA values
  if (res == "5") {
    list.models$N.obs_5m[k] <- sum(unit.occ[], na.rm = T)
  } else {
    list.models$N.obs_15m[k] <- sum(unit.occ[], na.rm = T)
  }
  
  
  # Check if occurrences match env data
  coorXY <- xyFromCell(envData, 1:ncell(envData)) # Coordinates of centroid of pixels
  unit.env.df <- getValues(unit.env) # Values of the Env stack, in a df, per layers
  
  ### 1.2/ Correction of occurrences without env data by looking around for values in a limited range ####
  
  Err.test <- 0 # Register failure, to start with
  
  # Check if any occurrence pixel as an NA value for envrionnmental layer (just need to check for the first one since all NA of env layers are synchronized)
  if(any(is.na(unit.env.df[, "bio1"]) & !is.na(unit.env.df[, "layer"]))) {
  
    Err.test <- 1 # Register presence of errors
    index.error <- which(is.na(unit.env.df[, "bio1"]) & !is.na(unit.env.df[, "layer"])) # Retreive indices for pixels with errors
    
    # Display how many occurrence points are outside the bioclimatic mask
    cat(length(index.error), "occurrences of", unit ,"are in pixels without environnmental values\n") 
    
    coorXY_errors <- coorXY[index.error,] # Retrieve coordinates of pixels in error
    res.error <- data.frame(index.error)  # Create a df to store result of research
    
    # Retrieve env infos in neighboring pixels around the occurrence point
    for (i in seq_along(index.error)) {
      index <- index.error[i]
      Error.occ <- unit.env.df[index,]  # Extract occurence line
      buffer <- 5000                    # 1st buffer = 5km
      while (anyNA(Error.occ) && (buffer<=50000)) { # Enlarge buffer while no data has beeen retrieve or 50km radius is reached
        if (length(index.error) == 1) { # Different format when only once error
          Error.occ <- unlist(raster::extract(x = envData, y = t(as.matrix(coorXY_errors)), buffer = buffer, fun = mean))
        }else{
          Error.occ <- unlist(raster::extract(x = envData, y = t(as.matrix(coorXY_errors[i,])), buffer = buffer, fun = mean))
        }
        buffer <- buffer + 5000 # New buffer = + 5km
      }
      
      if (!anyNA(Error.occ)) { # When research for data succeeded
        
        cat("Buffer",index,"for",unit,"=",(buffer-5000)/1000,"km\n") # Display buffer final radius
        unit.env.df[index,2:ncol(unit.env.df)] <- Error.occ # Write env infos in the recorded location of the occurrence
        res.error$Res[i] <- 1 # Register the success
        
      }else{ # When research for data failed
        
        cat("Occurrence",index,"for",unit,"has been removed because no env. data was available in a 50km radius\n")
        unit.env.df[index, 1] <- NA # Delete the occurrence
        res.error$Res[i] <- 0 # Register the failure
      }
      
    }
    
    for (j in 1:nlayers(unit.env)) {
      unit.env@layers[[j]]@data@values <- unit.env.df[,j] # Replace data in the EnvData stack
    }  
    
  } else { # Case when all occurrences have env. data
    
    cat(unit,": All occurrences present Env. data\n")
    
  }
  
  # Save the modified occurrence layer independently
  Binary_Raster_Occ <- unit.env[[1]] # Name the layer of usable Occurrences with a generic name
  save(Binary_Raster_Occ, file = paste0("./input_data/Species_data/",res,"/Binary_Rasters_from_Occurrences/",unit,".RData")) # Save in .RData
  saveRDS(Binary_Raster_Occ, file = paste0("./input_data/Species_data/",res,"/Binary_Rasters_from_Occurrences/",unit,".rds")) # Save in .rds
  
  # Compute the new number of occurrences and extract their coordinates
  coorXY <- coorXY[-which(is.na(unit.env.df[, "bio1"])), ] # Coordinates of pixels with env data
  unit.env.df <- unit.env.df[-which(is.na(unit.env.df[, "bio1"])), ] # Extract only lines with env data
  
  # Store info on nb of occurrences after correction
  if (res == "5") {
    list.models$N.obs_5m_used[list.models$Tag.model == unit] <- length(which(unit.env.df[, 1] == 1))
  } else {
    list.models$N.obs_15m_used[list.models$Tag.model == unit] <- length(which(unit.env.df[, 1] == 1))
  }
  
  cat(unit, "Number of pixels of presence:", length(which(unit.env.df[, 1] == 1)), "\n")  # Nb of valid presence points (with env. values associated)
  # cat(unit, "Number of pixels of absence:", length(which(unit.env.df[, 1] == 0)), "\n") # Nb of valid absence points (with env. values associated)
  
  coorXY_valid <- coorXY[which(!is.na(unit.env.df[, 1])), ] # Coordinates of pixels with env. values and occurrences data = the ones we want to keep
  if (length(coorXY_valid) == 2) { # Format correction, when only one occurrence
    coorXY_valid <- data.frame(t(as.matrix(coorXY_valid)))
  }
  
  ### 1.3/ Generate sp object (package sp) with only valid occurrences, and coordinates of pixel centroids ####
  
  unit.points <- SpatialPointsDataFrame(coorXY_valid, # Coordinates of pixels with env. values and occurrences data = the ones we want to keep
                                        data = data.frame(Occurrence = unit.env.df[which(!is.na(unit.env.df[, 1])), 1]), # Extract values of occurrences (0/1)
                                        proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")) # Define the CRS
  
  ### 1.4/ Plot occurrence correction results per units ####
  
  pdf(file = paste0("./input_data/Species_data/",res,"/Occurrence_correction_map_",unit,".pdf"), height = 8, width = 8) 
  original_ext_margins <- par()$oma ; original_int_margins <- par()$mar
  par(oma = c(0,2,0,3), mar = c(5.1,2.1,4.1,2.1), xpd=NA)
  par(mfrow = c(1,1)); plot(envData[[1]], las = 1, main = unit)
  points(unit.coords[, c(1:2)], cex = 1, pch = 16) # Plot original occurrence points
  plot(unit.points, add = T, cex = 0.6, pch = 16, col = "blue") # Add the one with env data
  if (Err.test) { # Case when there was points with no data
    if (length(index.error) == 1) { # Case with only one error to correct
      if (res.error$Res) { # Case when the point was corrected
        points(coorXY_errors, cex = 0.6, pch = 16, col = "orange") # Add retrieved data points
      }else{ #  Case when the point was not corrected
        points(coorXY_errors, cex = 0.6, pch = 16, col = "red") # Add failed data points
      } # Case with multiple errors to correct
    }else{
      points(coorXY_errors[res.error$Res,], cex = 0.6, pch = 16, col = "orange") # Add retrieved data points
      points(coorXY_errors[!res.error$Res,], cex = 0.6, pch = 16, col = "red") # Add failed data points
    }
  }
  legend(legend = c("OK","Extrapoled","Deleted"), pch = 16, pt.cex = 1.3, col = c("blue", "orange", "red"), x = "bottomleft", cex = 0.9, bty ="o")
  par(mar = original_ext_margins, oma = original_ext_margins, xpd = F)
  dev.off()
  
  
  
  ## Copy in the Species folder for Maps
  sp <- sp.list[k]
  # Check if output folder exists
  if (!dir.exists(paths = paste0("./maps/By_species/",sp))) {
    dir.create(path = paste0("./maps/By_species/",sp), recursive = T) # Create it if not existing
  } 
  
  # Copy Occurrence map in Species folder
  file.copy(from = paste0("./input_data/Species_data/",res,"/Occurrence_correction_map_",unit,".pdf"), to = paste0("./maps/By_species/",sp,"/Occurrence_correction_map_",unit,".pdf"), overwrite = T)
  
  # Save final Spatial Object
  save(unit.points, file = paste0("./input_data/Species_data/",res,"/Spatial_Points_Objects/occurrences_", unit,".RData")) # Sauvegarde du nouvel objet spatial des P/A
  
  # Save environmental stack with modifs for extrapolated data
  unit.env <- dropLayer(x = unit.env, i = 1) # Remove Occurrence layer
  save(unit.env, file = paste0("./input_data/Species_data/",res,"/Env_Stacks/Env_stack_",unit,".RData")) 
  
  cat(paste0(Sys.time(), " - Occurences Filtering for ", unit, " = Unit N°",k," - Done\n"))
  
}


# save(list.models, file = "./input_data/list.models.RData")




##### 2/ Parallelized version #####

# Changes =

# Add library loading in the loop
# Modify lines that assign value in object outside of the loop : list.models
#    N.obs et dérivés
# Generate a return(list()) for these objects that are not saved at each iteration, typically infos you want to store in a vector/list/table for each iteration
# Add script to geenrate and close the cluster
# foreach ... %dopar%


# Remove environment
rm(list = ls())

### 2.1/ Set up the cluster ####

library(foreach)
library(doParallel)

# Detect the number of threads
numCores <- detectCores()
# registerDoParallel(numCores)  # To set the nb if threads to use when using %dopar%

cl <- makeCluster(numCores) # Generate the virtual cluster from the number of threads
registerDoParallel(cl)      # To set the virtual cluster when using %dopar%

# getDoParWorkers() # Pour vérifier le nombre de coeurs enregistrés


### Load stuff useful for all units (outside the loop) ####

### Load table to track model outputs
load(file = "./input_data/list.models.RData")
unit.list <- list.models$Tag.model
sp.list <- list.models$Sp_full

# list.models$N.obs_5m <- list.models$N.obs_5m_used <- list.models$N.obs_15m <- list.models$N.obs_15m_used <- list.models$NA


### Load occurrence dataset
load(file = "./input_data/Ithomiini_final.RData")

### Select Env data resolution 

# Chose resolution
# res <- "5"
res <- "15"

# Load environmental stack
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))


### Rasterize occurrences, retrieve env data, generate spatial object


# for (k in seq_along(unit.list)) {
temp <- foreach (k = 2:3) %dopar% {
  # k <- 25
  
  # Load libraries
  library(raster)
  
  unit <- as.character(unit.list[k])
  cat(as.character(Sys.time()), "----- Start for", unit, "= Unit N°",k,"\n")
  unit.coords <- Ithomiini_final[(Ithomiini_final$Tag.model == unit), c("Longitude","Latitude")] 
  
  ### 2.2/ Rasterization of occurrences ####
  unit.occ <- rasterize(x = unit.coords, # Occurrence coordinates
                        y = envData, # Env stack grid to apply
                        field = 1, # Value to store in the raster
                        fun = function(x, ...) mean(x)) # Apply mean to keep value = 1 in case of multiple occurrences in one pixel
  
  unit.env <- stack(unit.occ, envData) # Add the rasterized occurrence layer to the env stack
  
  # Count occurrences after rasterization and before correction for NA values
  if (res == "5") {
    N.obs_5m <- sum(unit.occ[], na.rm = T)
    N.obs_15m <- NA
  } else {
    N.obs_5m <- NA
    N.obs_15m <- sum(unit.occ[], na.rm = T)
  }
  
  # Check if occurrences match env data
  coorXY <- xyFromCell(envData, 1:ncell(envData)) # Coordinates of centroid of pixels
  unit.env.df <- getValues(unit.env) # Values of the Env stack, in a df, per layers
  
  ### 2.3/ Correction of occurrences without env data by looking around for values in a limited range ####
  
  Err.test <- 0 # Register failure, to start with
  
  # Check if any occurrence pixel as an NA value for envrionnmental layer (just need to check for the first one since all NA of env layers are synchronized)
  if(any(is.na(unit.env.df[, "bio1"]) & !is.na(unit.env.df[, "layer"]))) {
    
    Err.test <- 1 # Register presence of errors
    index.error <- which(is.na(unit.env.df[, "bio1"]) & !is.na(unit.env.df[, "layer"])) # Retreive indices for pixels with errors
    
    # Display how many occurrence points are outside the bioclimatic mask
    cat(length(index.error), "occurrences of", unit ,"are in pixels without environnmental values\n") 
    
    coorXY_errors <- coorXY[index.error,] # Retrieve coordinates of pixels in error
    res.error <- data.frame(index.error)  # Create a df to store result of research
    
    # Retrieve env infos in neighboring pixels around the occurrence point
    for (i in seq_along(index.error)) {
      index <- index.error[i]
      Error.occ <- unit.env.df[index,]  # Extract occurence line
      buffer <- 5000                    # 1st buffer = 5km
      while (anyNA(Error.occ) && (buffer<=50000)) { # Enlarge buffer while no data has beeen retrieve or 50km radius is reached
        if (length(index.error) == 1) { # Different format when only once error
          Error.occ <- unlist(raster::extract(x = envData, y = t(as.matrix(coorXY_errors)), buffer = buffer, fun = mean))
        }else{
          Error.occ <- unlist(raster::extract(x = envData, y = t(as.matrix(coorXY_errors[i,])), buffer = buffer, fun = mean))
        }
        buffer <- buffer + 5000 # New buffer = + 5km
      }
      
      if (!anyNA(Error.occ)) { # When research for data succeeded
        
        cat("Buffer",index,"for",unit,"=",(buffer-5000)/1000,"km\n") # Display buffer final radius
        unit.env.df[index,2:ncol(unit.env.df)] <- Error.occ # Write env infos in the recorded location of the occurrence
        res.error$Res[i] <- 1 # Register the success
        
      }else{ # When research for data failed
        
        cat("Occurrence",index,"for",unit,"has been removed because no env. data was available in a 50km radius\n")
        unit.env.df[index, 1] <- NA # Delete the occurrence
        res.error$Res[i] <- 0 # Register the failure
      }
      
    }
    
    for (j in 1:nlayers(unit.env)) {
      unit.env@layers[[j]]@data@values <- unit.env.df[,j] # Replace data in the EnvData stack
    }  
    
  } else { # Case when all occurrences have env. data
    
    cat(unit,": All occurrences present Env. data\n")
    
  }
  
  # Save the modified occurrence layer independently
  Binary_Raster_Occ <- unit.env[[1]] # Name the layer of usable Occurrences with a generic name
  save(Binary_Raster_Occ, file = paste0("./input_data/Species_data/",res,"/Binary_Rasters_from_Occurrences/",unit,".RData")) # Save in .RData
  saveRDS(Binary_Raster_Occ, file = paste0("./input_data/Species_data/",res,"/Binary_Rasters_from_Occurrences/",unit,".rds")) # Save in .rds
  
  # Compute the new number of occurrences and extract their coordinates
  coorXY <- coorXY[-which(is.na(unit.env.df[, "bio1"])), ] # Coordinates of pixels with env data
  unit.env.df <- unit.env.df[-which(is.na(unit.env.df[, "bio1"])), ] # Extract only lines with env data
  
  # Store info on nb of occurrences after correction
  if (res == "5") {
    N.obs_5m_used <- length(which(unit.env.df[, 1] == 1))
    N.obs_15m_used <-  NA
  } else {
    N.obs_5m_used <- NA
    N.obs_15m_used <- length(which(unit.env.df[, 1] == 1))
  }
  
  cat(unit, "Number of pixels of presence:", length(which(unit.env.df[, 1] == 1)), "\n")  # Nb of valid presence points (with env. values associated)
  # cat(unit, "Number of pixels of absence:", length(which(unit.env.df[, 1] == 0)), "\n") # Nb of valid absence points (with env. values associated)
  
  coorXY_valid <- coorXY[which(!is.na(unit.env.df[, 1])), ] # Coordinates of pixels with env. values and occurrences data = the ones we want to keep
  if (length(coorXY_valid) == 2) { # Format correction, when only one occurrence
    coorXY_valid <- data.frame(t(as.matrix(coorXY_valid)))
  }
  
  ### 2.4/ Generate sp object (package sp) with only valid occurrences, and coordinates of pixel centroids ####
  
  unit.points <- SpatialPointsDataFrame(coorXY_valid, # Coordinates of pixels with env. values and occurrences data = the ones we want to keep
                                        data = data.frame(Occurrence = unit.env.df[which(!is.na(unit.env.df[, 1])), 1]), # Extract values of occurrences (0/1)
                                        proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0")) # Define the CRS
  
  ### 2.5/ Plot occurrence correction results per units ####
  
  pdf(file = paste0("./input_data/Species_data/",res,"/Occurrence_correction_maps/Occurrence_correction_map_",unit,".pdf"), height = 8, width = 8) 
  original_ext_margins <- par()$oma ; original_int_margins <- par()$mar
  par(oma = c(0,2,0,3), mar = c(5.1,2.1,4.1,2.1), xpd=NA)
  par(mfrow = c(1,1)); plot(envData[[1]], las = 1, main = unit)
  points(unit.coords[, c(1:2)], cex = 1, pch = 16) # Plot original occurrence points
  plot(unit.points, add = T, cex = 0.6, pch = 16, col = "blue") # Add the one with env data
  if (Err.test) { # Case when there was points with no data
    if (length(index.error) == 1) { # Case with only one error to correct
      if (res.error$Res) { # Case when the point was corrected
        points(coorXY_errors, cex = 0.6, pch = 16, col = "orange") # Add retrieved data points
      }else{ #  Case when the point was not corrected
        points(coorXY_errors, cex = 0.6, pch = 16, col = "red") # Add failed data points
      } # Case with multiple errors to correct
    }else{
      points(coorXY_errors[res.error$Res,], cex = 0.6, pch = 16, col = "orange") # Add retrieved data points
      points(coorXY_errors[!res.error$Res,], cex = 0.6, pch = 16, col = "red") # Add failed data points
    }
  }
  legend(legend = c("OK","Extrapoled","Deleted"), pch = 16, pt.cex = 1.3, col = c("blue", "orange", "red"), x = "bottomleft", cex = 0.9, bty ="o")
  par(mar = original_ext_margins, oma = original_ext_margins, xpd = F)
  dev.off()
  
  
  
  ## Copy in the Species folder for Maps
  sp <- sp.list[k]
  # Check if output folder exists
  if (!dir.exists(paths = paste0("./maps/By_species/",sp))) {
    dir.create(path = paste0("./maps/By_species/",sp), recursive = T) # Create it if not existing
  } 
  
  # Copy Occurrence map in Species folder
  file.copy(from = paste0("./input_data/Species_data/",res,"/Occurrence_correction_maps/Occurrence_correction_map_",unit,".pdf"), to = paste0("./maps/By_species/",sp,"/Occurrence_correction_map_",unit,".pdf"), overwrite = T)
  
  # Save final Spatial Object
  save(unit.points, file = paste0("./input_data/Species_data/",res,"/Spatial_Points_Objects/occurrences_", unit,".RData")) # Sauvegarde du nouvel objet spatial des P/A
  
  # Save environnmental stack with modifs for extrapolated data
  unit.env <- dropLayer(x = unit.env, i = 1) # Remove Occurrence layer
  save(unit.env, file = paste0("./input_data/Species_data/",res,"/Env_Stacks/Env_stack_",unit,".RData")) 
  
  # Return objects linked to the index of the loop
  return(list(index = k, N.obs_5m = N.obs_5m, N.obs_5m_used = N.obs_5m_used, N.obs_15m = N.obs_15m, N.obs_15m_used = N.obs_15m_used))
  
  cat(paste0(Sys.time(), " - Occurences Filtering for ", unit, " = Unit N°",k," - Done\n"))
  
}

# Close clusters generated with registerDoParallel(numCores)
# stopImplicitCluster()
# Close clusters generated with makeCluster 
stopCluster() 


### 2.6/ Manage outputs and save them in the summary table: list.models ####

# For whatever reason, the reshaping of the output as a row.bind df does not work, so need to do it manually afterwards
temp <- dplyr::as_tibble(temp)

# Insert results in the summary table
for (i in seq_along(temp)) {
  list.models$N.obs_5m[temp$index[i]] <- temp$N.obs_5m[i]
  list.models$N.obs_15m[temp$index[i]] <- temp$N.obs_15m[i]
  list.models$N.obs_5m_used[temp$index[i]] <- temp$N.obs_5m_used[i]
  list.models$N.obs_15m_used[temp$index[i]] <- temp$N.obs_15m_used[i]
}

# save(list.models, file = "./input_data/list.models.RData")


# Try with a custom function to deal with as many variables as wanted... but it does not work for some reason... ####

actualize_summary_table <- function(df_parallel, df_summary, vars, index = "index") {
  
  for (j in seq_along(vars)) {
    var <- vars[j]
    
    for (i in seq_along(df_parallel)) {
      df_summary[df_parallel[i, names(df_parallel) == index], names(df_summary) == var] <- df_parallel[i, names(df_summary) == var]
    }
  }
  
  return(df_summary)
}

test <- actualize_summary_table(temp, list.models, vars = c("N.obs_5m", "N.obs_15m", "N.obs_5m_used", "N.obs_15m_used"))


### 3/ Plot histogram of occurrences availability, depending of the resolution chosen ####

### 3.1/ Play around to define threshold for modeling types ####

hist(list.models$N.obs_5m_used[list.models$N.obs_5m_used<20], breaks  = seq(0,20,2))
sum((list.models$N.obs_5m_used >= 10) & (list.models$N.obs_5m_used < 30))
sum((list.models$N.obs_5m_used >= 5) & (list.models$N.obs_5m_used < 30)) 
sum(list.models$N.obs_5m_used < 6)
sum(list.models$N.obs_5m_used < 10)
sum(list.models$N.obs_5m_used >= 30)

table(list.models$initial_model_type)

hist(list.models$N.obs_15m_used[list.models$N.obs_15m_used<20], breaks  = seq(0,20,2))
sum((list.models$N.obs_15m_used >= 10) & (list.models$N.obs_15m_used < 30))
restricted_15 <- sum((list.models$N.obs_15m_used >= 6) & (list.models$N.obs_15m_used < 30)) 
rasterized_15 <- sum(list.models$N.obs_15m_used < 6) ; rastesterized_15
sum(list.models$N.obs_15m_used < 10)
complete_15 <- sum(list.models$N.obs_15m_used >= 30)

perc_complete_15 <- round(complete_15/nrow(list.models)*100,1)
perc_restricted_15 <- round(restricted_15/nrow(list.models)*100,1)
perc_rasterized_15 <- round(rasterized_15/nrow(list.models)*100,1)

obs_data <- data.frame(N.obs = c(list.models$N.obs_5m_used, list.models$N.obs_15m_used), Type = c(rep("5m", times = nrow(list.models)), rep("15m", times = nrow(list.models))))

library(tidyverse)

# Histogram
obs_data %>% 
  filter(obs_data$N.obs <= 20) %>% 
  ggplot(aes(x = N.obs, fill = Type)) + geom_histogram(position = "dodge", binwidth = 1) +
  geom_vline(xintercept = c(5, 10))

# Density plot
obs_data %>% 
  filter(obs_data$N.obs <= 30) %>% 
  ggplot(aes(x = N.obs, fill = Type, color = Type)) + geom_density(alpha = 0.3) +
  geom_vline(xintercept = c(5, 10))


pdf(file = paste0("./input_data/Species_data/Comparison_occurrences_btw_res.pdf"), height = 8, width = 8)
# Bar plot with customed categories
obs_data_barplot <- obs_data
obs_data_barplot$N.obs <- obs_data_barplot$N.obs %>%
  cut(breaks = c(0, 6, 10, 30), right = F)
obs_data_barplot %>% 
  ggplot(aes(x = N.obs, fill = Type, color = Type)) + geom_bar(position = "dodge") +
  geom_text(stat = 'count', aes(label = stat(count)), color = "black", nudge_y = 10)
dev.off()

table(list.models$N.obs_15m_used)

scales::show_col(colors(), labels = FALSE)

### 3.2/ Plot bar chart for occurrences at resolution 15min ###

obs_data_15 <- data.frame(N.obs = list.models$N.obs_15m_used)
obs_data_15$model_type <- as.factor(cut(x = obs_data_15$N.obs, breaks = c(0, 5, 30, 500), right = F))

rasterized_sp <- length(unique(list.models$Sp_full[obs_data_15$model_type == "[0,5)"]))
restricted_sp <- length(unique(list.models$Sp_full[obs_data_15$model_type == "[5,30)"]))
complete_sp <- length(unique(list.models$Sp_full[obs_data_15$model_type == "[30,500)"]))

modeled_sp <- length(unique(list.models$Sp_full[obs_data_15$model_type != "[0,5)"]))

# Function to increase vertical spacing between legend keys

draw_key_polygon3 <- function(data, params, size) {
  lwd <- min(data$size, min(size) / 4)
  
  grid::rectGrob(
    width = grid::unit(0.6, "npc"),
    height = grid::unit(0.6, "npc"),
    gp = grid::gpar(
      col = data$colour,
      fill = alpha(data$fill, data$alpha),
      lty = data$linetype,
      lwd = lwd * .pt,
      linejoin = "mitre"
    ))
}

# Register new key drawing function, 
# The effect is global & persistent throughout the R session
GeomBar$draw_key = draw_key_polygon3

pdf(file = "./supplementaries/Distribution_occurrences.pdf", width = 8, height = 6)

# Histogram
obs_data_15 %>%
  # filter(obs_data$N.obs <= 50) %>%
  ggplot(aes(x = N.obs, fill = model_type)) + 
  geom_histogram(binwidth = 3, position = "identity") +
  # geom_vline(xintercept = c(5, 30)) +
  
  scale_fill_manual(name = "Model type",
                    breaks = c("[0,5)", "[5,30)", "[30,500)"),
                    labels = c(paste0('"Rasterized" models (N < 6)\n', rasterized_15 ,' OMUs (', perc_rasterized_15,'%)'),
                               paste0('"Restricted" models (6 \u2264 N < 30)\n', restricted_15 ,' OMUs (', perc_restricted_15,'%)'),
                               paste0('"Complete" models (N \u2265 30)\n', complete_15 ,' OMUs (', perc_complete_15,'%)')),
                    values = c("#FF0000BB", "#FFA500BB", "#32CD32BB")) +
  
  ylab("Counts") +
  xlab("Number of occurrences per OMU") +
  
  labs(fill = "Model type") +
  
  # title(main = "Distribution of occurrences per OMU") +
  
  # # Rectangle for "rasterized" OMU
  # annotate("rect", xmin = -1.5, xmax = 5, ymin = 0, ymax = 252,
  #          col = NA, fill = "red", alpha = 0.5) +
  # 
  # # Rectangle for "restricted" OMU
  # annotate("rect", xmin = 5, xmax = 30, ymin = 0, ymax = 200,
  #          col = NA, fill = "orange", alpha = 0.5) +
  # 
  # # Rectangle for "restricted" OMU
  # annotate("rect", xmin = 30, xmax = 400, ymin = 0, ymax = 50,
  #          col = NA, fill = "limegreen", alpha = 0.5) +
  
  ### Manual legend
  
  # # Rectangle for "rasterized" OMU
  # annotate("rect", xmin = 100, xmax = 140, ymin = 220, ymax = 240,
  #          col = "grey80", fill = "red", alpha = 1) +
  
  ggplot2::theme(panel.background = ggplot2::element_rect(fill = NA),
                 plot.title = ggplot2::element_text(size = 18, face = "bold", hjust = 0.5, vjust = 3),
                 legend.position = c(0.75, 0.55),
                 legend.title = ggplot2::element_text(size = 15, vjust = 2, hjust = 0.05, face = "bold"),
                 legend.text = ggplot2::element_text(size = 12, face = "bold", hjust = 0.5, margin = margin(l = -10, unit = "pt")),
                 legend.key = element_rect(size = 50, color = NA, fill = NA),
                 legend.key.size = unit(40, 'pt'),
                 # plot.margin = ggplot2::margin(t = 18, unit = "pt"),
                 # axis.line.x = ggplot2::element_line(color = NA),
                 axis.line = ggplot2::element_line(color = "black", size = 1.3),
                 axis.ticks = ggplot2::element_line(color = "black", size = 1.2),
                 axis.ticks.length = ggplot2::unit(8, "pt"),
                 axis.text = ggplot2::element_text(size = 16, face = "bold", color = "black"),
                 axis.text.x = ggplot2::element_text(vjust = -4, size = 14, margin = ggplot2::margin(t = -5, b = 20)),
                 axis.title = ggplot2::element_text(size = 16, face = "bold", color = "black"),
                 # axis.title.x = ggplot2::element_text(margin = margin(t = 10, b = 5)),
                 axis.title.y = ggplot2::element_text(margin = ggplot2::margin(l = 5, r = 10)))

dev.off()


### 4/ Add initial_model_type ####

# Decide of initial modeling type depending of the N.obs_used => Add a column for that "complete" (N >= 30), "restricted" (N > ? et N < 30), "rasterized" (N > ?)
# Threshold between "restricted" and "raserized" to decide afterwards regarding the evaluation metrics results. Still possible to discard low N.obs models even if they "succeeded"

# Do it manually afterwards, because the loop must also run with resolution = 5m

list.models$initial_model_type <- NA
for (i in 1:nrow(list.models)) {
  
  N.obs <- list.models$N.obs_15m_used[i]
  
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

# save(list.models, file = "./input_data/list.models.RData")


