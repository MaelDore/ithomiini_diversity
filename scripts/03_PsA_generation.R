

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
  # PDF plot of probability distribution of selected PsA for each set, for each OMU
  # PDF map of PsA set for each OMU



##### 1/ Compute distance matrix between sampling sites ####

# Remove environment
rm(list = ls())

# Need to define resolution of final models right now, to filter the proper number of occurrences using the final modeling grid of all outputs

# Chose resolution
# res <- "5"
res <- "15"

# Load associated environmental stack to use as grid to rasterize
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))

# Load occurrence dataset and extract occurrences coordinates
load(file = "./input_data/Ithomiini_final.RData")
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
# 1834 sampling sites

# Get coordinates of pixels associated with each occurrence
xy <- foreach (i = seq_along(cell_indices), .combine = rbind) %do% {
  xyFromCell(envData, cell_indices[i])
}
  
raster_infos <- as.tibble(cbind(cell_indices, sampling_ID, xy))
raster_infos <- raster_infos %>% 
  dplyr::rename(Longitude_raster = x, Latitude_raster = y)

# Store infos on occurrences on raster grid in df of occurrences
if (length(dplyr::intersect(names(Ithomiini_final), names(raster_infos))) > 0) {
  Ithomiini_final <- Ithomiini_final %>% 
    select(-names(raster_infos))
}
Ithomiini_final <- Ithomiini_final %>%
  bind_cols(raster_infos)

# save(Ithomiini_final, file = "./input_data/Ithomiini_final.Rdata")


### Generate table of sampling sites coordinates
sampling.sites_coords_with_index_and_cell_nb <- raster_infos %>% 
  unique() %>% 
  arrange(sampling_ID)

save(sampling.sites_coords_with_index_and_cell_nb, file = paste0("./input_data/Sampling_sites/sampling.sites_coords_with_index_and_cell_nb_",res,".RData"))

sampling.sites_coords <- sampling.sites_coords_with_index_and_cell_nb %>% 
  select(Longitude_raster, Latitude_raster)

save(sampling.sites_coords, file = paste0("./input_data/Sampling_sites/sampling.sites_coords_",res,".RData"))


### Compute distance matrix of all sampling sites
Sampling.sites.Dist = geosphere::distm(x = sampling.sites_coords)/1000 # Geometric distance on WGS84 ellipsoid, in km

save(Sampling.sites.Dist, file = paste0("./input_data/Sampling_sites/Sampling.sites.Dist_",res,".RData"))


##### 2/ Get a record of sampling sites that do not have associated env. data to discard them from the potential PsA pool ####

load(file = "./input_data/Sampling_sites/sampling.sites_coords.RData")

# Chose resolution
# res <- "5"
res <- "15"

# Load associated environmental stack to use as grid to rasterize
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))

for (i in 1:nrow(sampling.sites_coords)) {
  
  # Extract env values for each sampling sites
  test <- raster::extract(envData, sampling.sites_coords[i, 1:2])
  
  # If any env values is missing, note to not use this sampling site as potential PsA site
  sampling.sites_coords$env_data[i] <- all(!is.na(test))
}

sum(!sampling.sites_coords$env_data) # 116 sampling sites out of 1834 are not available to draw PsA because they lack env data

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
load(file = "./input_data/list.models.RData")
unit.list <- list.models$Tag.model
sp.list <- list.models$Sp_full

### Load occurrence dataset
load(file = "./input_data/Databases/Ithomiini_final.RData")

### Select Env data resolution 

# Chose resolution
# res <- "5"
res <- "15"

# Load environmental stack
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))

### Load sampling sites infos for PsA Generation
load(file = paste0("./input_data/Sampling_sites/Sampling.sites.Dist_",res,".RData"))
load(file = paste0("./input_data/Sampling_sites/sampling.sites_coords_",res,".RData"))

source("./functions/progcombine.R")

# Set seed to ensure repeatability of random PsA draws
set.seed(158340)


### 3.3/ Loop to generate PsA table for each units/OMU ####

temp <- foreach (k = seq_along(unit.list), .combine = progcombine_rbind(nreps = length(unit.list)), .init = NULL) %dopar% {
# temp <- foreach (k = 1:5) %dopar% {
  # k <- 1
  
  unit <- as.character(unit.list[k])
  N.obs <- list.models$N.obs_15m_used[k]
  
  ##### 3.3.1/ PsA generation ####
  
  ### Only for OMU/unit that will be modeled, to save time
  
  if (N.obs < 6) { 
    
    # Warning message for OMU with too low sample size to try a SDM
    cat(paste0(Sys.time(), " - No Pseudo-Absences Generation for ", unit, " = Unit N°",k,"\n"))
    
  } else { # PsA Generation for OMU with at least 6 obs.
    
    cat("\n", as.character(Sys.time()), "----- PsA generation starts for", unit, "= Unit N°",k,"\n")
    
    # load libraries
    
    library(tidyverse)
    library(raster)
    
    # Load occurrences coordiantes and spatial object
    unit.occurrence <- Ithomiini_final[Ithomiini_final$Tag.model == unit, c("Longitude","Latitude")]
    load(file = paste0("./input_data/Species_data/",res,"/Spatial_Points_Objects/occurrences_", unit,".RData"))  
    
    # Get species name for output paths
    sp <- sp.list[k]
    
    # Retrieve sampling sites ID for this OMU
    unit.sampling_ID <- unique(Ithomiini_final[Ithomiini_final$Tag.model == unit, "sampling_ID"])
    
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
      n.PsA.set <- 10 # 10 PsA for OMU with limited N.obs
    } else { 
      n.PsA.set <- 3  # 3 PsA sets for OMU with enough N.obs to do CV
    } 
    
    PA.indices.table <- data.frame(matrix(nrow = N.obs, ncol = n.PsA.set))
    for (i in 1:n.PsA.set) {
      PA.indices.table[,i] <- sample(x = index.potential.PsA, size = N.obs, prob = weights, replace = F)
    }
    PA.list <- as.vector(as.matrix(PA.indices.table)) # Concatenate in a single vector all the indices of sample sites to use as PsA
    
    PA_set <- as.factor(as.vector(matrix(data = 1:n.PsA.set, byrow = T, nrow = N.obs, ncol = n.PsA.set))) # Generate a vector of PsA set attribution
    ggplot.table <- tibble("indices" = PA.list, "PA_set" = PA_set, "dist" = Closest.dist[PA.list], "weights" = 1/Closest.dist[PA.list])
    
    ### 3.3.2/ Plot probability distribution of selected PsA for each set ####
    
    pdf(file = paste0("./input_data/Species_data/",res,"/Pseudo_absences/Proba_distri_PsA_",unit,".pdf"), height = 8, width = 8)
    g <- ggplot.table %>% 
      ggplot(aes(x = dist, y = weights, colour = PA_set)) +
      geom_point(size = 3, alpha = 0.5, fill = NA, position = position_jitter(width = max(ggplot.table$dist)/50, height = max(ggplot.table$weights)/50, seed = 5)) +
      guides(colour = guide_legend(override.aes = list(alpha = 1))) +
      labs(title = paste0("Probability distribution of PsA for ", unit),
           x = "Distance in km", y = "Weights")
    print(g)
    dev.off()  
    
    # Copy in the sp folder
    file.copy(from = paste0("./input_data/Species_data/",res,"/Pseudo_absences/Proba_distri_PsA_",unit,".pdf"), to = paste0("./maps/By_species/",sp,"/Proba_distri_PsA_",unit,".pdf"), overwrite = T)
    
    
    ### 3.3.3/ Plot selected PsA of each set on a map ####
    
    pdf(file = paste0("./input_data/Species_data/",res,"/Pseudo_absences/PsA_map_",unit,".pdf"), height = 8, width = 8) 
    plot(envData[["Elevation"]], main = unit)
    for (i in 1:n.PsA.set) {
      points(sampling.sites_coords[PA.indices.table[,i], 1:2], col = c("red","darkgreen","blue", "gold","cyan","limegreen","purple","grey","orange","deeppink")[i], pch = 1, cex = 0.7)
      points(unit.points, col = "red", cex = 1.5, pch = 16)
      points(unit.points, cex = 1, pch = 16)
    }
    legend(legend = c("Presences",paste0(rep("PA",n.PsA.set),1:n.PsA.set)), pch = c(16, rep(1,n.PsA.set)), 
           col = c("black","red","darkgreen", "blue","gold","cyan","limegreen","purple","grey","orange","deeppink")[1:n.PsA.set],
           x = "bottomleft", pt.cex = 1, cex = 0.9, bty ="o")
    dev.off()
    
    # Copy in the sp folder for Maps
    file.copy(from =  paste0("./input_data/Species_data/",res,"/Pseudo_absences/PsA_map_",unit,".pdf"), to = paste0("./maps/By_species/",sp,"/PsA_map_",unit,".jpeg"), overwrite = T)
    
    
    ### 3.3.4/ Retrieve coordinates of PsA and save in sp object ####
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
# stopImplicitCluster()
# Close clusters generated with makeCluster 
stopCluster() 
