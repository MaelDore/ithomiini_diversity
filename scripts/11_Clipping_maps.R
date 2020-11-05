
##### Script 11: Clipping maps #####

# Use alpha-hull around occurrences points and a buffer proportionnal to the quality of presence sampling to clip the maps
# Clip EM maps and individual submodel projections, for each OMU

#####

# Inputs:
   # Buffer for continental borders
   # Buffer for Andean regions
   # OMU occurrence coordinates (to draw the alpha-hull)
   # Stack of EM for each OMU
   # Stack of submodels for each OMU

# Outputs:
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
load(file = paste0("./input_data/list.models.RData"))
modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]

# Chose the resolution
res <-  "15"

# Load environmental stack to use as reference for extent and CRS
library(raster)
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))

# Load mask for continent borders
continent_mask <- readRDS(file = paste0("./input_data/Env_data/continent_mask_", res, ".rds"))

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

# Projection in Mercador (EPSG:3857) to allow intersection function to work properly
Andes_East <-  spTransform(Andes_East , CRSobj = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")
Andes_West <-  spTransform(Andes_West , CRSobj = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")

# plot(Andes_East)
# plot(Andes_West, add = T)

# Create the variable to keep track of location of OMU regarding Andes slopes
list.models$Andes_side <-  NA 

# Functions to convert alpha-hull into spatial objects
source("./functions/Alpha_functions.R", local = TRUE) 


############## Modif between scripts ######################
index_model_to_compute <- c(1:563)
# index_model_to_compute <- c(1:200)
# index_model_to_compute <- c(201:400)
# index_model_to_compute <- c(401:563)
###########################################################

### Loop for all OMU

# i <- 161
# Units en erreur à utiliser avec un alpha-shape plutôt que alpha-hull
# for (i in c(161, 173, 178, 207, 221, 237, 249, 270, 280, 293, 318, 341, 388, 444, 480, 494, 511, 528, 548)) {

### Loop for all modeled OMU/unit
# for (i in 1:nrow(modeled_OMU)) 
for (i in index_model_to_compute) 
{
  # i <- 69
  
  unit <- as.character(modeled_OMU$Tag.model[i])
  sp <- as.character(modeled_OMU$Sp_full[i])
  
  cat(paste0("\n", Sys.time()," ------ Starts for ", unit, " = Unit N°",i,"/",nrow(modeled_OMU)," ------\n"))
  
  ### Load unit occurrence points
  load(paste0("./input_data/Species_data/",res,"/Spatial_Points_Objects/occurrences_", unit,".RData")) # Load Spatial object with occurrences

  # Load Ensembles (continous, ensemble, binary, committee averaging) for Jaccard indices
  Ensemble_Jaccard_median <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_Jaccard_median.rds"))
  Ensemble_Jaccard_bin <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_Jaccard_bin.rds"))
  Ensemble_Jaccard_sd <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_Jaccard_sd.rds"))
  Ensemble_Jaccard_CA <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_Jaccard_CA.rds"))

  # all_EM_Jaccard <- stack(Ensemble_Jaccard_median, Ensemble_Jaccard_sd, Ensemble_Jaccard_bin, Ensemble_Jaccard_CA)
  # names(all_EM_Jaccard) <- c("Continuous", "Incertitude (sd)", "Binary", "Committee Averaging")
  # save(all_EM_Jaccard, file = paste0("./outputs/By_OMU/",unit,"/all_EM_Jaccard.RData"), version = "2")
  # saveRDS(all_EM_Jaccard, file = paste0("./outputs/By_OMU/",unit,"/all_EM_Jaccard.rds"), version = "2")

  # Load Ensembles (continous, ensemble, binary, committee averaging) for TSS indices
  Ensemble_TSS_median <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_TSS_median.rds"))
  Ensemble_TSS_bin <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_TSS_bin.rds"))
  Ensemble_TSS_sd <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_TSS_sd.rds"))
  Ensemble_TSS_CA <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Ensemble_TSS_CA.rds"))

  # all_EM_TSS <- stack(Ensemble_TSS_median, Ensemble_TSS_sd, Ensemble_TSS_bin, Ensemble_TSS_CA)
  # names(all_EM_TSS) <- c("Continuous", "Incertitude (sd)", "Binary", "Committee Averaging")
  # save(all_EM_TSS, file = paste0("./outputs/By_OMU/",unit,"/all_EM_TSS.RData"), version = "2")
  # saveRDS(all_EM_TSS, file = paste0("./outputs/By_OMU/",unit,"/all_EM_TSS.rds"), version = "2")

  # Load submodels stacks
  Select_Jaccard_stack_cont <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Select_Jaccard_stack_cont.rds"))
  Select_TSS_stack_cont <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Select_TSS_stack_cont.rds"))
  Select_Jaccard_stack_bin <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Select_Jaccard_stack_bin.rds"))
  Select_TSS_stack_bin <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Select_TSS_stack_bin.rds"))
  
  
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

  # Temporary projection of the spatial object of occurrences in Mercador (EPSG:3857)
  proj.points <-  spTransform(unit.points, CRSobj = "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")

  tryCatch(
   {
     Hull <- ahull_to_SPLDF(ahull(proj.points@coords, alpha = 1000000)) # Try an alpha-hull first
   },
   error = function(e) # If fails, do an alpha shape instead
     {
       cat(paste0("Alpha-hull failed for ", unit, ", use alpha-shape instead \n"))
       cat("ERROR :",conditionMessage(e), "\n")   # Display the error message but do not stop the function
       Hull <<- ashape_to_SPLDF(ashape(proj.points@coords, alpha = 1000000))
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
  Buffer_80 <- gBuffer(spgeom = Hull, width = max(111.32*1000, Uncertainty.80*1000))
  Buffer_95 <- gBuffer(spgeom = Hull, width = max(111.32*1000, Uncertainty.95*1000))

  # Remove holes
  Buffer_80_unholed <- spatialEco::remove.holes(Buffer_80)
  Buffer_95_unholed <- spatialEco::remove.holes(Buffer_95)

  # Retransposition in WGS84
  Buffer_80.WGS84 <- spTransform(Buffer_80_unholed, CRSobj = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))
  Buffer_95.WGS84 <- spTransform(Buffer_95_unholed, CRSobj = CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"))

  # plot(Ensemble_Jaccard_median)
  # points(unit.points, col = "red", cex = 1, pch = 16)
  # plot(Buffer_80.WGS84, border = "dodgerblue", lwd = 2, col = NA, add =T)
  # plot(Buffer_95.WGS84, lwd = 2, add =T)

  # Rasterization
  Buffer_80.rasterized <- rasterize(x = Buffer_80.WGS84, y = Ensemble_Jaccard_median, field = 1)
  Buffer_95.rasterized <- rasterize(x = Buffer_95.WGS84, y = Ensemble_Jaccard_median, field = 1)
  # plot(Buffer_80.rasterized, add = T)
  # plot(Buffer_95.rasterized, add = T)

  # Save final buffers
  save(Buffer_80.rasterized, file = paste0("./outputs/By_OMU/",unit,"/Buffer_80.RData"), version = "2")
  saveRDS(Buffer_80.rasterized, file = paste0("./outputs/By_OMU/",unit,"/Buffer_80.rds"), version = "2")
  save(Buffer_95.rasterized, file = paste0("./outputs/By_OMU/",unit,"/Buffer_95.RData"), version = "2")
  saveRDS(Buffer_95.rasterized, file = paste0("./outputs/By_OMU/",unit,"/Buffer_95.rds"), version = "2")

  
  ### 3/ Clip the EM and submodels maps ####
  
  ## 3.1/ Clip all EM and submodels with the buffers ####
  
  # Load buffers
  Buffer_80.rasterized <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Buffer_80.rds"))
  Buffer_95.rasterized <- readRDS(file = paste0("./outputs/By_OMU/",unit,"/Buffer_95.rds"))

  ## 3.1.1/ EM maps

  # Make a stack of EMs, both Eval metrics
  full_stack <- stack(Ensemble_Jaccard_median, Ensemble_Jaccard_sd, Ensemble_Jaccard_bin, Ensemble_Jaccard_CA, Ensemble_TSS_median, Ensemble_TSS_sd, Ensemble_TSS_bin, Ensemble_TSS_CA)

  # Apply buffers
  full_stack_cropped_80 <- mask(full_stack, Buffer_80.rasterized)
  full_stack_cropped_95 <- mask(full_stack, Buffer_95.rasterized)

  # Make a full stack for all Eval metric and buffers
  full_stack_cropped_by_buffers <- stack(full_stack_cropped_80, full_stack_cropped_95)

  # Create names for the layers
  Ensemble <- rep("Ensemble", 16)
  Eval <- c(rep("Jaccard", 4), rep("TSS", 4))
  Type <- c("median", "sd", "bin", "CA")
  Crop <- rep("cropped", 16)
  Buffer <- c(rep("80", 8), rep("95", 8))

  layer_names_df <- data.frame(Ensemble = Ensemble, Eval = Eval, Type = Type, Crop = Crop, Buffer = Buffer, stringsAsFactors = F)
  layer_names <- apply(X = layer_names_df, MARGIN = 1 , FUN = paste, collapse = "_")

  names(full_stack_cropped_by_buffers) <- layer_names

  ## 3.1.2/ For submodels maps
  
  # Apply buffers
  Select_Jaccard_stack_cont_cropped_80 <- mask(Select_Jaccard_stack_cont, Buffer_80.rasterized)
  Select_TSS_stack_cont_cropped_80 <- mask(Select_TSS_stack_cont, Buffer_80.rasterized)
  Select_Jaccard_stack_cont_cropped_95 <- mask(Select_Jaccard_stack_cont, Buffer_95.rasterized)
  Select_TSS_stack_cont_cropped_95 <- mask(Select_TSS_stack_cont, Buffer_95.rasterized)
  
  Select_Jaccard_stack_bin_cropped_80 <- mask(Select_Jaccard_stack_bin, Buffer_80.rasterized)
  Select_TSS_stack_bin_cropped_80 <- mask(Select_TSS_stack_bin, Buffer_80.rasterized)
  Select_Jaccard_stack_bin_cropped_95 <- mask(Select_Jaccard_stack_bin, Buffer_95.rasterized)
  Select_TSS_stack_bin_cropped_95 <- mask(Select_TSS_stack_bin, Buffer_95.rasterized)
  
  ## 3.2/ Check within which Andean slope lies the occurrences, and crop out the slope where no occurrences have been recorded ####
  
  # 3.2.1/ For EM maps

  # Checkif presence points fall into the Andes slopes (used projected object in Mercator)
  test_East <- gIntersection(spgeom1 = Andes_East, spgeom2 = proj.points, checkValidity = 2)
  test_West <- gIntersection(spgeom1 = Andes_West, spgeom2 = proj.points, checkValidity = 2)

  ### Plot
  # plot(Andes_East)
  # points(test_East, pch = 16, col ="limegreen")
  # points(test_West, pch = 16, col ="red")
  ###

  # Clip out the slope where no occurrences have been recorded
  if (length(test_East)*length(test_West) == 0) { # Case with occurrences only on one side

    if (length(test_West)==0) { # Case with occurrences on the East side

      full_stack_cropped_by_buffers_and_slopes <- mask(full_stack_cropped_by_buffers, Andes_East.rasterized)

      list.models$Andes_side[which(list.models$Tag.model == unit)] <- "East"

    } else { # Case with occurrences on the West side

      full_stack_cropped_by_buffers_and_slopes <- mask(full_stack_cropped_by_buffers, Andes_West.rasterized)

      list.models$Andes_side[which(list.models$Tag.model == unit)] <- "West"

    }

  } else { # Case with occurrences on both sides = no additionnal cropping

    full_stack_cropped_by_buffers_and_slopes <- full_stack_cropped_by_buffers

    list.models$Andes_side[which(list.models$Tag.model == unit)] <- "Both"
  }

  ### Add 0 value for terrestrial pixels to be able to distinguish from oceans

  full_stack_cropped_by_buffers_with_borders <- full_stack_cropped_by_buffers_and_slopes
  for (j in 1:nlayers(full_stack_cropped_by_buffers))
  {
    layer_name <- names(full_stack_cropped_by_buffers)[j]

    # Extract values from layer
    temp_data <- full_stack_cropped_by_buffers_and_slopes[[j]]@data@values
    temp_data[which(is.na(temp_data))] <- 0 # Add 0 instead of NA

    temp_layer <- Ensemble_Jaccard_median
    temp_layer@data@values <- temp_data
    names(temp_layer) <- layer_name

    # Mask with a layer of environnemental variable which have continental borders
    temp_layer <- mask(temp_layer, envData[[1]])

    # Store in new stack
    full_stack_cropped_by_buffers_with_borders[[j]] <- temp_layer

    # Save individually
    saveRDS(temp_layer, file = paste0("./outputs/By_OMU/",unit,"/",layer_name,".rds"), version = "2")

  }

  # Extract stacks per Eval metric and Buffer
  all_EM_Jaccard_buffer_80 <- full_stack_cropped_by_buffers_with_borders[[1:4]]
  all_EM_TSS_buffer_80 <- full_stack_cropped_by_buffers_with_borders[[5:8]]
  all_EM_Jaccard_buffer_95 <- full_stack_cropped_by_buffers_with_borders[[9:12]]
  all_EM_TSS_buffer_95 <- full_stack_cropped_by_buffers_with_borders[[13:16]]

  # Save stacks
  save(all_EM_Jaccard_buffer_80, file = paste0("./outputs/By_OMU/",unit,"/all_EM_Jaccard_buffer_80.RData"), version = "2")
  saveRDS(all_EM_Jaccard_buffer_80, file = paste0("./outputs/By_OMU/",unit,"/all_EM_Jaccard_buffer_80.rds"), version = "2")
  save(all_EM_Jaccard_buffer_95, file = paste0("./outputs/By_OMU/",unit,"/all_EM_Jaccard_buffer_95.RData"), version = "2")
  saveRDS(all_EM_Jaccard_buffer_95, file = paste0("./outputs/By_OMU/",unit,"/all_EM_Jaccard_buffer_95.rds"), version = "2")
  save(all_EM_TSS_buffer_80, file = paste0("./outputs/By_OMU/",unit,"/all_EM_TSS_buffer_80.RData"), version = "2")
  saveRDS(all_EM_TSS_buffer_80, file = paste0("./outputs/By_OMU/",unit,"/all_EM_TSS_buffer_80.rds"), version = "2")
  save(all_EM_TSS_buffer_95, file = paste0("./outputs/By_OMU/",unit,"/all_EM_TSS_buffer_95.RData"), version = "2")
  saveRDS(all_EM_TSS_buffer_95, file = paste0("./outputs/By_OMU/",unit,"/all_EM_TSS_buffer_95.rds"), version = "2")

  ## 3.2.2/ For Submodel maps
  
  # Get Andes side info
  Andes_side <- list.models$Andes_side[which(list.models$Tag.model == unit)]
  
  # Apply buffers
  if (Andes_side == "East") {
    
    Select_Jaccard_stack_cont_cropped_80 <- mask(Select_Jaccard_stack_cont_cropped_80, Andes_East.rasterized)
    Select_TSS_stack_cont_cropped_80 <- mask(Select_TSS_stack_cont_cropped_80, Andes_East.rasterized)
    Select_Jaccard_stack_cont_cropped_95 <- mask(Select_Jaccard_stack_cont_cropped_95, Andes_East.rasterized)
    Select_TSS_stack_cont_cropped_95 <- mask(Select_TSS_stack_cont_cropped_95, Andes_East.rasterized)
    
    Select_Jaccard_stack_bin_cropped_80 <- mask(Select_Jaccard_stack_bin_cropped_80, Andes_East.rasterized)
    Select_TSS_stack_bin_cropped_80 <- mask(Select_TSS_stack_bin_cropped_80, Andes_East.rasterized)
    Select_Jaccard_stack_bin_cropped_95 <- mask(Select_Jaccard_stack_bin_cropped_95, Andes_East.rasterized)
    Select_TSS_stack_bin_cropped_95 <- mask(Select_TSS_stack_bin_cropped_95, Andes_East.rasterized)
    
  } else {
    if (Andes_side == "West") {
      
      Select_Jaccard_stack_cont_cropped_80 <- mask(Select_Jaccard_stack_cont_cropped_80, Andes_West.rasterized)
      Select_TSS_stack_cont_cropped_80 <- mask(Select_TSS_stack_cont_cropped_80, Andes_West.rasterized)
      Select_Jaccard_stack_cont_cropped_95 <- mask(Select_Jaccard_stack_cont_cropped_95, Andes_West.rasterized)
      Select_TSS_stack_cont_cropped_95 <- mask(Select_TSS_stack_cont_cropped_95, Andes_West.rasterized)
      
      Select_Jaccard_stack_bin_cropped_80 <- mask(Select_Jaccard_stack_bin_cropped_80, Andes_West.rasterized)
      Select_TSS_stack_bin_cropped_80 <- mask(Select_TSS_stack_bin_cropped_80, Andes_West.rasterized)
      Select_Jaccard_stack_bin_cropped_95 <- mask(Select_Jaccard_stack_bin_cropped_95, Andes_West.rasterized)
      Select_TSS_stack_bin_cropped_95 <- mask(Select_TSS_stack_bin_cropped_95, Andes_West.rasterized)
    }
  }
  
  ### Add 0 value for terrestrial pixels to be able to distinguish from oceans
  Select_Jaccard_stack_cont_cropped_80 <- merge(x = continent_mask, y = Select_Jaccard_stack_cont_cropped_80)
  Select_Jaccard_stack_cont_cropped_95 <- merge(x = continent_mask, y = Select_Jaccard_stack_cont_cropped_95)
  Select_TSS_stack_cont_cropped_80 <- merge(x = continent_mask, y = Select_TSS_stack_cont_cropped_80)
  Select_TSS_stack_cont_cropped_95 <- merge(x = continent_mask, y = Select_TSS_stack_cont_cropped_95)
  
  Select_Jaccard_stack_bin_cropped_80 <- merge(x = continent_mask, y = Select_Jaccard_stack_bin_cropped_80)
  Select_Jaccard_stack_bin_cropped_95 <- merge(x = continent_mask, y = Select_Jaccard_stack_bin_cropped_95)
  Select_TSS_stack_bin_cropped_80 <- merge(x = continent_mask, y = Select_TSS_stack_bin_cropped_80)
  Select_TSS_stack_bin_cropped_95 <- merge(x = continent_mask, y = Select_TSS_stack_bin_cropped_95)
  
  
  # Save stacks as .RData/rds
  save(Select_Jaccard_stack_cont_cropped_80, file = paste0("./outputs/By_OMU/",unit,"/Select_Jaccard_stack_cont_cropped_80.RData"), version = "2")
  saveRDS(Select_Jaccard_stack_cont_cropped_80, file = paste0("./outputs/By_OMU/",unit,"/Select_Jaccard_stack_cont_cropped_80.rds"), version = "2")
  save(Select_TSS_stack_cont_cropped_80, file = paste0("./outputs/By_OMU/",unit,"/Select_TSS_stack_cont_cropped_80.RData"), version = "2")
  saveRDS(Select_TSS_stack_cont_cropped_80, file = paste0("./outputs/By_OMU/",unit,"/Select_TSS_stack_cont_cropped_80.rds"), version = "2")
  
  save(Select_Jaccard_stack_cont_cropped_95, file = paste0("./outputs/By_OMU/",unit,"/Select_Jaccard_stack_cont_cropped_95.RData"), version = "2")
  saveRDS(Select_Jaccard_stack_cont_cropped_95, file = paste0("./outputs/By_OMU/",unit,"/Select_Jaccard_stack_cont_cropped_95.rds"), version = "2")
  save(Select_TSS_stack_cont_cropped_95, file = paste0("./outputs/By_OMU/",unit,"/Select_TSS_stack_cont_cropped_95.RData"), version = "2")
  saveRDS(Select_TSS_stack_cont_cropped_95, file = paste0("./outputs/By_OMU/",unit,"/Select_TSS_stack_cont_cropped_95.rds"), version = "2")
  
  save(Select_Jaccard_stack_bin_cropped_80, file = paste0("./outputs/By_OMU/",unit,"/Select_Jaccard_stack_bin_cropped_80.RData"), version = "2")
  saveRDS(Select_Jaccard_stack_bin_cropped_80, file = paste0("./outputs/By_OMU/",unit,"/Select_Jaccard_stack_bin_cropped_80.rds"), version = "2")
  save(Select_TSS_stack_bin_cropped_80, file = paste0("./outputs/By_OMU/",unit,"/Select_TSS_stack_bin_cropped_80.RData"), version = "2")
  saveRDS(Select_TSS_stack_bin_cropped_80, file = paste0("./outputs/By_OMU/",unit,"/Select_TSS_stack_bin_cropped_80.rds"), version = "2")
  
  save(Select_Jaccard_stack_bin_cropped_95, file = paste0("./outputs/By_OMU/",unit,"/Select_Jaccard_stack_bin_cropped_95.RData"), version = "2")
  saveRDS(Select_Jaccard_stack_bin_cropped_95, file = paste0("./outputs/By_OMU/",unit,"/Select_Jaccard_stack_bin_cropped_95.rds"), version = "2")
  save(Select_TSS_stack_bin_cropped_95, file = paste0("./outputs/By_OMU/",unit,"/Select_TSS_stack_bin_cropped_95.RData"), version = "2")
  saveRDS(Select_TSS_stack_bin_cropped_95, file = paste0("./outputs/By_OMU/",unit,"/Select_TSS_stack_bin_cropped_95.rds"), version = "2")
  
  
  ### 4/ Plot the final cropped EM ####

  # Function adding the points to the plot on each layer
  add_points <- function() {
    plot(unit.points, add = TRUE, col = "red", pch = 16, cex = 0.3)
  }

  # 4.1/ Plot the 4 EM for 1 Eval metric and 1 Buffer, on a single pdf page
  all_EM_list <- list(all_EM_Jaccard_buffer_80, all_EM_Jaccard_buffer_95, all_EM_TSS_buffer_80, all_EM_TSS_buffer_95)
  all_EM_names_list <- c("all_EM_Jaccard_buffer_80", "all_EM_Jaccard_buffer_95", "all_EM_TSS_buffer_80", "all_EM_TSS_buffer_95")

  for (j in 1:length(all_EM_list))
  {
    EM_stack <- all_EM_list[[j]]
    names(EM_stack) <- c("Continuous", "Incertitude (sd)", "Binary", "Committee Averaging")
    EM_stack_name <- all_EM_names_list[j]
    pdf(file = paste0("./maps/By_OMU/",unit,"/",EM_stack_name,"_",unit,".pdf"), height = 6, width = 7)
    internal_margins <- par()$mar
    external_margins <- par()$oma
    par(mar = c(5.1,6.1,9.1,3.1))
    par(oma = c(1,2,4,3))
    plot(EM_stack, addfun = add_points)
    title(main = unit, outer = T)
    par(mar = internal_margins, oma = external_margins)
    dev.off()
    # Copy in species folder
    # sp <- as.character(modeled_OMU$Sp_full[i])
    file.copy(from = paste0("./maps/By_OMU/",unit,"/",EM_stack_name,"_",unit,".pdf"), to = paste0("./maps/By_species/",sp,"/",EM_stack_name,"_",unit,".pdf"), overwrite = T)
  }

  # 4.2/ Plot the full projection and the two cropped one, for continuous and binary maps

  # Jaccard
  Final_stack <- stack(Ensemble_Jaccard_median,   # Full continuous map
                       Ensemble_Jaccard_bin,      # Full binary map
                       all_EM_Jaccard_buffer_80[["Ensemble_Jaccard_median_cropped_80"]],  # Buffer 80 continuous map
                       all_EM_Jaccard_buffer_80[["Ensemble_Jaccard_bin_cropped_80"]],     # Buffer 80 binary map
                       all_EM_Jaccard_buffer_95[["Ensemble_Jaccard_median_cropped_95"]],  # Buffer 95 continuous map
                       all_EM_Jaccard_buffer_95[["Ensemble_Jaccard_bin_cropped_95"]])     # Buffer 95 binary map

  names(Final_stack) <- c("Continuous Map", "Binary Map", "Continuous Map\nBuffer 80%", "Binary Map\nBuffer 80%", "Continuous Map\nBuffer 95%", "Binary Map\nBuffer 95%")

  pdf(file = paste0("./maps/By_OMU/",unit,"/All_buffers_cont_bin_Jaccard_",unit,".pdf"), height = 9, width = 7)
  internal_margins <- par()$mar
  external_margins <- par()$oma
  par(mar = c(5.1,6.1,9.1,3.1))
  par(oma = c(1,2,4,3))
  plot(Final_stack, addfun = add_points, nc = 2)
  title(main = unit, outer = T)
  par(mar = internal_margins, oma = external_margins)
  dev.off()
  # Copy in species folder
  # sp <- as.character(modeled_OMU$Sp_full[i])
  file.copy(from = paste0("./maps/By_OMU/",unit,"/All_buffers_cont_bin_Jaccard_",unit,".pdf"), to = paste0("./maps/By_species/",sp,"/All_buffers_cont_bin_Jaccard_",unit,".pdf"), overwrite = T)

  # TSS
  Final_stack <- stack(Ensemble_TSS_median,   # Full continuous map
                       Ensemble_TSS_bin,      # Full binary map
                       all_EM_TSS_buffer_80[["Ensemble_TSS_median_cropped_80"]],  # Buffer 80 continuous map
                       all_EM_TSS_buffer_80[["Ensemble_TSS_bin_cropped_80"]],     # Buffer 80 binary map
                       all_EM_TSS_buffer_95[["Ensemble_TSS_median_cropped_95"]],  # Buffer 95 continuous map
                       all_EM_TSS_buffer_95[["Ensemble_TSS_bin_cropped_95"]])     # Buffer 95 binary map

  names(Final_stack) <- c("Continuous Map", "Binary Map", "Continuous Map\nBuffer 80%", "Binary Map\nBuffer 80%", "Continuous Map\nBuffer 95%", "Binary Map\nBuffer 95%")

  pdf(file = paste0("./maps/By_OMU/",unit,"/All_buffers_cont_bin_TSS_",unit,".pdf"), height = 9, width = 7)
  internal_margins <- par()$mar
  external_margins <- par()$oma
  par(mar = c(5.1,6.1,9.1,3.1))
  par(oma = c(1,2,4,3))
  plot(Final_stack, addfun = add_points, nc = 2)
  title(main = unit, outer = T)
  par(mar = internal_margins, oma = external_margins)
  dev.off()
  # Copy in species folder
  # sp <- as.character(modeled_OMU$Sp_full[i])
  file.copy(from = paste0("./maps/By_OMU/",unit,"/All_buffers_cont_bin_TSS_",unit,".pdf"), to = paste0("./maps/By_species/",sp,"/All_buffers_cont_bin_TSS_",unit,".pdf"), overwrite = T)

  cat(paste0("\n", Sys.time()," ------ Done for ", unit, " = Unit N°",i,"/",nrow(modeled_OMU)," ------\n"))
 
} 

# Save Summary table for OMU models
# save(list.models, file = paste0("./input_data/list.models.RData"), version = "2")
# saveRDS(list.models, file = paste0("./input_data/list.models.rds"), version = "2")


# plot(list.models$Uncertainty.80 ~ list.models$Uncertainty.95)
# 
# hist(list.models$Uncertainty.80)
# hist(list.models$Uncertainty.95)
