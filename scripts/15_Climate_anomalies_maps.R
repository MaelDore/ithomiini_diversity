

##### Script 15: Build climatic anomalies maps to assess threats on Ithomiini diversity #####

# Build maps of climate threats based on multivariate standardized climatic anomalies 

# Data from Bioclim 2.1

# 2 scenarios from the CMIP6 : SSP2 - RCP 4.5  and SSP5 - RCP 8.5

# 8 to 6 Global Climate Models (GCM) 

# 3 ways to compute anomalies
  # Standardized by baseline variance in space
  # Standardized by raw anomalies variance in space
  # With normalized anomalies between -1 and 1

# Multivariate anomalies are computed as euclidian distance in the multivariate space


#####

# Inputs:
  # Current climate conditions
  # Climatic projections under SSP2-RCP4.5 and SSP5-RCP8.5, for 6-8 GCMs

# Outputs:
  # Climatic anomalies EM maps. Median and incertitude (sd). Global, Andes and Mata Atlantica

#####


# Remove environment
rm(list = ls())

# Load libraries
library(raster)

### Select Env data resolution 
# res <- "5"
res <- "15"

# ### Create baseline (average on 1970-2000 = 1985)
# clim_baseline <- NA
# for (i in c(1, 4, 12, 15)) {
#   
#   if(i == 1) 
#   { 
#     clim_baseline <- stack(paste0("./input_data/Bioclim/wc2.1_5m_bio/wc2.1_5m_bio_", i,".tif")) # Create the stack
#     
#   } else { # Add layer iteratively
#       
#     clim_baseline <- addLayer(clim_baseline, raster(paste0("./input_data/Bioclim/wc2.1_5m_bio/wc2.1_5m_bio_", i,".tif")))
#   }
# }
# names(clim_baseline) <- c("Tmean", "Tvar", "Ptot", "Pvar") # Create layer name
# plot(clim_baseline)
# 
# # Aggregate and crop to final maps resolution
# 
# # Load environmental stack to use as template for resolution and extent
# envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))
# clim_baseline <- resample(x = clim_baseline, y = envData, method = "bilinear") # Ajust extent and resolution to match the one of our SDM maps
# plot(clim_baseline)
# 
# # Save baseline
# saveRDS(clim_baseline, file = "./input_data/Env_data/anomalies/clim_baseline.rds", version = "2")

### Load baseline (average on 1970-2000 = 1985)
clim_baseline <- readRDS(file = "./input_data/Env_data/anomalies/clim_baseline.rds")

### Load environmental stack to use as template for resolution, extent, and terrestrial boundaries
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))


library(tmaptools)

pal_grn_red <- rev(get_brewer_pal("RdYlGn", n = 200))

### For 2070 ####

# Set Time period
period <- "2061-2080"

# Per SSP-RCP
SSP_list <- c("ssp245", "ssp585")
# SSP <- SSP_list[1]

# Per GCM
# GCM_list <- c("BCC-CSM2-MR", "CNRM-CM6-1", "CNRM-ESM2-1", "CanESM5", "IPSL-CM6A-LR", "MIROC-ES2L", "MIROC6", "MRI-ESM2-0") ; GCM_set <- 8
GCM_list <- c("BCC-CSM2-MR", "CNRM-CM6-1", "CNRM-ESM2-1", "CanESM5", "IPSL-CM6A-LR", "MIROC6") ; GCM_set <- 6
# GCM <- GCM_list[1]

# ID for GCM_set

for (SSP in SSP_list) 
{
  
  cat(paste("\n------------- ", Sys.time(),"Start process for all GCMs for SSP =", SSP,"------------- \n"))
  
  for (GCM in GCM_list) 
  {
    cat(paste("\n------ ", Sys.time(),"Start process for GCM", GCM,"------ \n"))
    
    # Create directory to store maps by GCM
    if(!dir.exists(paste0("./maps/anomalies/",SSP,"/",GCM))) { # Test if the folder exists already or not
      dir.create(paste0("./maps/anomalies/",SSP,"/",GCM), recursive = T) # Create folder if absent
    }
    
    
    ### 1/ Extract stack of projected variables ####
    
    # Load projection stack for this SSP and this GCM
    clim_stack <- stack(x = paste0("./input_data/Bioclim/wc2.1_5m_",SSP,"/wc2.1_5m_bioc_",GCM,"_",SSP,"_",period,".tif"))
    # nlayers(clim_stack)
    # names(clim_stack)
    
    # Extract variables used in SDM
    clim_stack <- subset(x = clim_stack, subset = c(1, 4, 12, 15) , drop = F)
    names(clim_stack) <- c("Tmean", "Tvar", "Ptot", "Pvar")
    # plot(clim_stack)
    
    # Aggregate and crop to final maps resolution
    clim_stack <- resample(x = clim_stack, y = envData, method = "bilinear") # Ajust extent and resolution to match the one of our SDM maps
    # plot(clim_stack)
    
    # Remove pixels not computed for indices maps
    # clim_stack <- mask(x = clim_stack, mask = envData[[1]])
    
    # Save stack for this period, this SSP and this GCM
    saveRDS(clim_stack, file = paste0("./input_data/Env_data/anomalies/clim_stack_",period,"_",SSP,"_",GCM,".rds"), version = "2")
    
    ### 2/ Compute anomalies compare to baseline ####
    
    # Raw anomalies 
    anomalies_raw_stack <- clim_stack - clim_baseline
    
    # Plot Raw anomalies
    pdf(file = paste0("./maps/anomalies/",SSP,"/",GCM,"/raw_anomalies_",period,"_",SSP,"_",GCM,".pdf"), height = 8, width = 8)
    internal_margins <- par()$mar
    external_margins <- par()$oma
    par(mar = c(5.1,6.1,9.1,3.1))
    par(oma = c(1,2,4,3))
    plot(anomalies_raw_stack, col = pal_grn_red)
    title(main = paste0("Raw anomalies ",period,"\n",SSP," ",GCM), outer = T)
    par(mar = internal_margins, oma = external_margins)
    dev.off()
    
    # Standardized anomalies based on baseline variability (sd)
    baseline_sd_list <- cellStats(clim_baseline, sd)
    
    anomalies_base_std_stack <- anomalies_raw_stack/baseline_sd_list
    
    # Standardized anomalies based on anomalies sd (don't change the mean, only scale by dividing by sd)
    # Limits: will erase the differences between scenarios since higher anomalies tend to be more variable...
    
    # Do not use that because it will apply the root mean square instead of the sd !
    # anomalies_std_stack <- scale(anomalies_raw_stack, center = F, scale = T) 
    
    # Do the standardization "manually" instead
    sd_list <- cellStats(anomalies_raw_stack, sd)
    anomalies_std_stack <- anomalies_raw_stack/sd_list
    
    # plot(anomalies_std_stack, col = pal_grn_red)
    
    # Normalized anomalies (don't change the mean, only divide by the max absolute value to give equal weights for each variable)
    anomalies_norm_stack <- anomalies_raw_stack / cellStats(abs(anomalies_raw_stack), max)
    # plot(anomalies_norm_stack, col = pal_grn_red)
    
    # Save anomalies stacks
    save(anomalies_raw_stack, anomalies_base_std_stack, anomalies_std_stack, anomalies_norm_stack, file = paste0("./outputs/anomalies/anomalies_stacks_",period,"_",SSP,"_",GCM,".RData"), version = "2")
    # # Save sd values to allow eventual common standardization between the SSPs later
    # saveRDS(sd_list, file = paste0("./outputs/anomalies/sd_list_",period,"_",SSP,"_",GCM,".RData"), version = "2")
    
    ### 3/ Compute global anomalies as euclidian distance to origin (i.e., the baseline) in 4D ####
    
    # Since the baseline is the origin of the multivariate space of anomalies, the euclidian distance to the baseline is sqrt(sum(coordinates²))
    
    # For all baseline-standardized variables
    anomalies_base_std_all <- calc(x = anomalies_base_std_stack, 
                              fun = function(x) 
                              { y <- sqrt(sum(x*x)) ; return(y)} )
    names(anomalies_base_std_all) <- "Base std anomalies - All var"
    
    # For baseline-standardized Tmean and Ptot only
    anomalies_base_std_no_var <- calc(x = subset(anomalies_base_std_stack, c("Tmean", "Ptot")), 
                                 fun = function(x) 
                                 { y <- sqrt(sum(x*x)) ; return(y)} )
    names(anomalies_base_std_no_var) <- "Base std anomalies - Tmean & Ptot"
    
    # For all standardized variables
    anomalies_std_all <- calc(x = anomalies_std_stack, 
                              fun = function(x) 
                              { y <- sqrt(sum(x*x)) ; return(y)} )
    names(anomalies_std_all) <- "Std anomalies - All var"
    
    # For standardized Tmean and Ptot only
    anomalies_std_no_var <- calc(x = subset(anomalies_std_stack, c("Tmean", "Ptot")), 
                                 fun = function(x) 
                                 { y <- sqrt(sum(x*x)) ; return(y)} )
    names(anomalies_std_no_var) <- "Std anomalies - Tmean & Ptot"
    
    # For normalized 
    anomalies_norm_all <- calc(x = anomalies_norm_stack, 
                               fun = function(x) 
                               { y <- sqrt(sum(x*x)) ; return(y)} )
    names(anomalies_norm_all) <- "Normalized anomalies - All var"
    
    # For normalized Tmean and Ptot only
    anomalies_norm_no_var <- calc(x = subset(anomalies_norm_stack, c("Tmean", "Ptot")), 
                                  fun = function(x) 
                                  { y <- sqrt(sum(x*x)) ; return(y)} )
    names(anomalies_norm_no_var) <- "Normalized anomalies - Tmean & Ptot"
    
    # Save global anomalies
    save(anomalies_base_std_all, anomalies_base_std_no_var, anomalies_std_all, anomalies_std_no_var, anomalies_norm_all, anomalies_norm_no_var, file = paste0("./outputs/anomalies/global_anomalies_",period,"_",SSP,"_",GCM,".RData"), version = "2")
    saveRDS(anomalies_base_std_all, file = paste0("./outputs/anomalies/anomalies_base_std_all_",period,"_",SSP,"_",GCM,".rds"), version = "2")
    saveRDS(anomalies_base_std_no_var, file = paste0("./outputs/anomalies/anomalies_base_std_no_var_",period,"_",SSP,"_",GCM,".rds"), version = "2")
    saveRDS(anomalies_std_all, file = paste0("./outputs/anomalies/anomalies_std_all_",period,"_",SSP,"_",GCM,".rds"), version = "2")
    saveRDS(anomalies_std_no_var, file = paste0("./outputs/anomalies/anomalies_std_no_var_",period,"_",SSP,"_",GCM,".rds"), version = "2")
    saveRDS(anomalies_norm_all, file = paste0("./outputs/anomalies/anomalies_norm_all_",period,"_",SSP,"_",GCM,".rds"), version = "2")
    saveRDS(anomalies_norm_no_var, file = paste0("./outputs/anomalies/anomalies_norm_no_var_",period,"_",SSP,"_",GCM,".rds"), version = "2")
    
    # 4/ Plot global anomalies ####
    
    global_anomalies_stack <- stack(anomalies_base_std_all, anomalies_base_std_no_var, anomalies_std_all, anomalies_std_no_var, anomalies_norm_all, anomalies_norm_no_var)
    
    pdf(file = paste0("./maps/anomalies/",SSP,"/",GCM,"/global_anomalies_",period,"_",SSP,"_",GCM,".pdf"), height = 12, width = 8)
    internal_margins <- par()$mar
    external_margins <- par()$oma
    par(mar = c(5.1,6.1,9.1,3.1))
    par(oma = c(1,2,4,3))
    plot(global_anomalies_stack, col = pal_grn_red, nc = 2)
    title(main = paste0("Multivariate Anomalies ",period,"\n",SSP," ",GCM), outer = T)
    par(mar = internal_margins, oma = external_margins)
    dev.off()
  }
  
  # 5/ Stack all GCM, extract median and sd ####
  
  # 5.1/ Ensemble for Base std Anomalies ####
  SSP_anomalies_base_std_all_stack <-  NA
  for (i in 1:length(GCM_list)) 
  {
    GCM <- GCM_list[i]
    
    if(i == 1) 
    { 
      SSP_anomalies_base_std_all_stack <- stack(readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_all_",period,"_",SSP,"_",GCM,".rds"))) # Create the stack
      
    } else { # Add layer iteratively
      
      SSP_anomalies_base_std_all_stack <- addLayer(SSP_anomalies_base_std_all_stack, readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_all_",period,"_",SSP,"_",GCM,".rds")))
    }
  }
  
  # plot(SSP_anomalies_base_std_all_stack, col = pal_grn_red)
  
  # Compute median Ensemble
  SSP_anomalies_base_std_all_median <- calc(x = SSP_anomalies_base_std_all_stack, fun = median, na.rm = T)
  names(SSP_anomalies_base_std_all_median) <- "Base std anomalies - All var"
  
  plot(SSP_anomalies_base_std_all_median, col = pal_grn_red)
  
  # hist(SSP_anomalies_base_std_all_median)
  
  # Quick contrasted plot by changing zlim with min = 0.5
  # plot(SSP_anomalies_base_std_all_median, col = pal_grn_red, zlim = c(0.5, max(SSP_anomalies_base_std_all_median[], na.rm = T)))
  
  if (GCM_set == 8) {
    if (SSP == "ssp245") {min_cut <- 0.4; max_cut <- 1.1} else {min_cut <- 0.5 ; max_cut <- 1.1}
  } else  {
    if (SSP == "ssp245") {min_cut <- 0.4; max_cut <- 0.7} else {min_cut <- 0.6 ; max_cut <- 1.15}
  }
  
  
  # Contrast by mergin values < cut
  SSP_anomalies_base_std_all_median_contrasted <- calc(x = SSP_anomalies_base_std_all_median, fun = function(x){ x[x < min_cut] <- min_cut ; x[x > max_cut] <- max_cut ; return(x)} )
  names(SSP_anomalies_base_std_all_median_contrasted) <- "Contrasted base std anomalies - All var"

  # plot(SSP_anomalies_base_std_all_median_contrasted, col = pal_grn_red)

  # Compute Incertitude (sd)
  SSP_anomalies_base_std_all_sd <- calc(x = SSP_anomalies_base_std_all_stack, fun = sd, na.rm = T)
  names(SSP_anomalies_base_std_all_sd) <- "Base std anomalies - All var"

  # hist(SSP_anomalies_base_std_all_sd)
  
  # No need for contrast

  # Contrast by mergin values < 0 (No contrast)
  SSP_anomalies_base_std_all_sd_contrasted <- calc(x = SSP_anomalies_base_std_all_sd, fun = function(x){ x[x < 0] <- 0 ; return(x)} )
  names(SSP_anomalies_base_std_all_sd_contrasted) <- "Contrasted base std anomalies - All var"

  # plot(SSP_anomalies_base_std_all_sd_contrasted, col = pal_grn_red)
  
  # Save ensemble median, sd, (and their contrasted version) for this SSP
  saveRDS(SSP_anomalies_base_std_all_median, file = paste0("./outputs/anomalies/anomalies_base_std_all_median_",period,"_",SSP,".rds"), version = "2")
  saveRDS(SSP_anomalies_base_std_all_median_contrasted, file = paste0("./outputs/anomalies/anomalies_base_std_all_median_contrasted_",period,"_",SSP,".rds"), version = "2")
  saveRDS(SSP_anomalies_base_std_all_sd, file = paste0("./outputs/anomalies/anomalies_base_std_all_sd_",period,"_",SSP,".rds"), version = "2")
  saveRDS(SSP_anomalies_base_std_all_sd_contrasted, file = paste0("./outputs/anomalies/anomalies_base_std_all_sd_contrasted_",period,"_",SSP,".rds"), version = "2")
  
  
  # 5.2/ Ensemble for Base std Anomalies - Only Tmean & Ptot ####
  SSP_anomalies_base_std_no_var_stack <-  NA
  for (i in 1:length(GCM_list)) 
  {
    GCM <- GCM_list[i]
    
    if(i == 1) 
    { 
      SSP_anomalies_base_std_no_var_stack <- stack(readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_no_var_",period,"_",SSP,"_",GCM,".rds"))) # Create the stack
      
    } else { # Add layer iteratively
      
      SSP_anomalies_base_std_no_var_stack <- addLayer(SSP_anomalies_base_std_no_var_stack, readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_no_var_",period,"_",SSP,"_",GCM,".rds")))
    }
  }
  
  # plot(SSP_anomalies_base_std_no_var_stack, col = pal_grn_red)
  
  # Compute median Ensemble
  SSP_anomalies_base_std_no_var_median <- calc(x = SSP_anomalies_base_std_no_var_stack, fun = median, na.rm = T)
  names(SSP_anomalies_base_std_no_var_median) <- "Base std anomalies - Tmean & Ptot"
  
  # plot(SSP_anomalies_base_std_no_var_median, col = pal_grn_red)
  
  # hist(SSP_anomalies_base_std_no_var_median)
  
  if (GCM_set == 8) {
    if (SSP == "ssp245") {min_cut <- 0.3 ; max_cut <- 10} else {min_cut <- 0.5 ; max_cut <- 10}
  } else {
    if (SSP == "ssp245") {min_cut <- 0.3 ; max_cut <- 0.7} else {min_cut <- 0.5 ; max_cut <- 1.05}
  }
  
  
  # Quick contrasted plot by changing zlim with min = 0.5
  # plot(SSP_anomalies_base_std_no_var_median, col = pal_grn_red, zlim = c(0.5, max(SSP_anomalies_base_std_no_var_median[], na.rm = T)))
  
  # Contrast by mergin values < 0.5
  SSP_anomalies_base_std_no_var_median_contrasted <- calc(x = SSP_anomalies_base_std_no_var_median, fun = function(x){ x[x < min_cut] <- min_cut ; x[x > max_cut] <- max_cut ; return(x)} )
  names(SSP_anomalies_base_std_no_var_median_contrasted) <- "Contrasted base std anomalies - Tmean & Ptot"

  # plot(SSP_anomalies_base_std_no_var_median_contrasted, col = pal_grn_red)
  
  # Compute Incertitude (sd)
  SSP_anomalies_base_std_no_var_sd <- calc(x = SSP_anomalies_base_std_no_var_stack, fun = sd, na.rm = T)
  names(SSP_anomalies_base_std_no_var_sd) <- "Base std anomalies - Tmean & Ptot"
  
  # hist(SSP_anomalies_base_std_no_var_sd)
  
  # plot(SSP_anomalies_base_std_no_var_sd, col = pal_grn_red)
  
  ifelse(SSP == "ssp245", cut <- 0.3, cut <- 0.4)
  
  # Contrasted plot by changing zlim with max = 0.4
  # plot(SSP_anomalies_base_std_no_var_sd, col = pal_grn_red, zlim = c(min(SSP_anomalies_base_std_no_var_sd[], na.rm = T), 0.4))
  
  # Contrast by mergin values > 0.4
  SSP_anomalies_base_std_no_var_sd_contrasted <- calc(x = SSP_anomalies_base_std_no_var_sd, fun = function(x){ x[x > cut] <- cut ; return(x)} )
  names(SSP_anomalies_base_std_no_var_sd_contrasted) <- "Contrasted base std anomalies - Tmean & Ptot"
  
  # plot(SSP_anomalies_base_std_no_var_sd_contrasted, col = pal_grn_red)
  
  # Save ensemble median, sd, and their contrasted version for this SSP
  saveRDS(SSP_anomalies_base_std_no_var_median, file = paste0("./outputs/anomalies/anomalies_base_std_no_var_median_",period,"_",SSP,".rds"), version = "2")
  saveRDS(SSP_anomalies_base_std_no_var_median_contrasted, file = paste0("./outputs/anomalies/anomalies_base_std_no_var_median_contrasted_",period,"_",SSP,".rds"), version = "2")
  saveRDS(SSP_anomalies_base_std_no_var_sd, file = paste0("./outputs/anomalies/anomalies_base_std_no_var_sd_",period,"_",SSP,".rds"), version = "2")
  saveRDS(SSP_anomalies_base_std_no_var_sd_contrasted, file = paste0("./outputs/anomalies/anomalies_base_std_no_var_sd_contrasted_",period,"_",SSP,".rds"), version = "2")
  
  
  # 5.3/ Ensemble for Std Anomalies
  SSP_anomalies_std_all_stack <-  NA
  for (i in 1:length(GCM_list)) 
  {
    GCM <- GCM_list[i]
    
    if(i == 1) 
    { 
      SSP_anomalies_std_all_stack <- stack(readRDS(file = paste0("./outputs/anomalies/anomalies_std_all_",period,"_",SSP,"_",GCM,".rds"))) # Create the stack
      
    } else { # Add layer iteratively
      
      SSP_anomalies_std_all_stack <- addLayer(SSP_anomalies_std_all_stack, readRDS(file = paste0("./outputs/anomalies/anomalies_std_all_",period,"_",SSP,"_",GCM,".rds")))
    }
  }
  
  # plot(SSP_anomalies_std_all_stack, col = pal_grn_red)
  
  # Compute median Ensemble
  SSP_anomalies_std_all_median <- calc(x = SSP_anomalies_std_all_stack, fun = median, na.rm = T)
  names(SSP_anomalies_std_all_median) <- "Std anomalies - All var"
  
  # hist(SSP_anomalies_std_all_median)
  
  # Contrast by mergin values < 4
  SSP_anomalies_std_all_median_contrasted <- calc(x = SSP_anomalies_std_all_median, fun=function(x){ x[x < 4] <- 4 ; return(x)} )
  names(SSP_anomalies_std_all_median_contrasted) <- "Contrasted std anomalies - All var"
  
  # plot(SSP_anomalies_std_all_median_contrasted, col = pal_grn_red)
  
  # Compute Incertitude (sd)
  SSP_anomalies_std_all_sd <- calc(x = SSP_anomalies_std_all_stack, fun = sd, na.rm = T)
  names(SSP_anomalies_std_all_sd) <- "Std anomalies - All var"
  
  # hist(SSP_anomalies_std_all_sd)
  
  # Contrast by mergin values > 2.0
  SSP_anomalies_std_all_sd_contrasted <- calc(x = SSP_anomalies_std_all_sd, fun=function(x){ x[x > 2.0] <- 2.0 ; return(x)} )
  names(SSP_anomalies_std_all_sd_contrasted) <- "Contrasted std anomalies - All var"
  
  # plot(SSP_anomalies_std_all_sd_contrasted, col = pal_grn_red)
  
  # Save ensemble median, sd, and their contrasted version for this SSP
  saveRDS(SSP_anomalies_std_all_median, file = paste0("./outputs/anomalies/anomalies_std_all_median_",period,"_",SSP,".rds"), version = "2")
  saveRDS(SSP_anomalies_std_all_median_contrasted, file = paste0("./outputs/anomalies/anomalies_std_all_median_contrasted_",period,"_",SSP,".rds"), version = "2")
  saveRDS(SSP_anomalies_std_all_sd, file = paste0("./outputs/anomalies/anomalies_std_all_sd_",period,"_",SSP,".rds"), version = "2")
  saveRDS(SSP_anomalies_std_all_sd_contrasted, file = paste0("./outputs/anomalies/anomalies_std_all_sd_contrasted_",period,"_",SSP,".rds"), version = "2")
  
  
  # 5.4/ Ensemble for Std Anomalies - Only Tmean & Ptot
  SSP_anomalies_std_no_var_stack <-  NA
  for (i in 1:length(GCM_list)) 
  {
    GCM <- GCM_list[i]
    
    if(i == 1) 
    { 
      SSP_anomalies_std_no_var_stack <- stack(readRDS(file = paste0("./outputs/anomalies/anomalies_std_no_var_",period,"_",SSP,"_",GCM,".rds"))) # Create the stack
      
    } else { # Add layer iteratively
      
      SSP_anomalies_std_no_var_stack <- addLayer(SSP_anomalies_std_no_var_stack, readRDS(file = paste0("./outputs/anomalies/anomalies_std_no_var_",period,"_",SSP,"_",GCM,".rds")))
    }
  }
  
  # plot(SSP_anomalies_std_no_var_stack, col = pal_grn_red)
  
  # Compute median Ensemble
  SSP_anomalies_std_no_var_median <- calc(x = SSP_anomalies_std_no_var_stack, fun = median, na.rm = T)
  names(SSP_anomalies_std_no_var_median) <- "Std anomalies - Tmean & Ptot"
  
  # hist(SSP_anomalies_std_no_var_median)
  
  # Contrast by mergin values < 3.5
  SSP_anomalies_std_no_var_median_contrasted <- calc(x = SSP_anomalies_std_no_var_median, fun=function(x){ x[x < 3.5] <- 3.5 ; return(x)} )
  names(SSP_anomalies_std_no_var_median_contrasted) <- "Contrasted std anomalies - Tmean & Ptot"
  
  # plot(SSP_anomalies_std_no_var_median_contrasted, col = pal_grn_red)
  
  # Compute Incertitude (sd)
  SSP_anomalies_std_no_var_sd <- calc(x = SSP_anomalies_std_no_var_stack, fun = sd, na.rm = T)
  names(SSP_anomalies_std_no_var_sd) <- "Std anomalies - Tmean & Ptot"
  
  # hist(SSP_anomalies_std_no_var_sd)
  
  # Contrast by mergin values > 2.0
  SSP_anomalies_std_no_var_sd_contrasted <- calc(x = SSP_anomalies_std_no_var_sd, fun=function(x){ x[x > 2.0] <- 2.0 ; return(x)} )
  names(SSP_anomalies_std_no_var_sd_contrasted) <- "Contrasted std anomalies - Tmean & Ptot"
  
  # plot(SSP_anomalies_std_no_var_sd_contrasted, col = pal_grn_red)
  
  # Save ensemble median, sd, and their contrasted version for this SSP
  saveRDS(SSP_anomalies_std_no_var_median, file = paste0("./outputs/anomalies/anomalies_std_no_var_median_",period,"_",SSP,".rds"), version = "2")
  saveRDS(SSP_anomalies_std_no_var_median_contrasted, file = paste0("./outputs/anomalies/anomalies_std_no_var_median_contrasted_",period,"_",SSP,".rds"), version = "2")
  saveRDS(SSP_anomalies_std_no_var_sd, file = paste0("./outputs/anomalies/anomalies_std_no_var_sd_",period,"_",SSP,".rds"), version = "2")
  saveRDS(SSP_anomalies_std_no_var_sd_contrasted, file = paste0("./outputs/anomalies/anomalies_std_no_var_sd_contrasted_",period,"_",SSP,".rds"), version = "2")
  
  
  # 5.5/ Ensemble for Normalized Anomalies
  SSP_anomalies_norm_all_stack <-  NA
  for (i in 1:length(GCM_list)) 
  {
    GCM <- GCM_list[i]
    
    if(i == 1) 
    { 
      SSP_anomalies_norm_all_stack <- stack(readRDS(file = paste0("./outputs/anomalies/anomalies_norm_all_",period,"_",SSP,"_",GCM,".rds"))) # Create the stack
      
    } else { # Add layer iteratively
      
      SSP_anomalies_norm_all_stack <- addLayer(SSP_anomalies_norm_all_stack, readRDS(file = paste0("./outputs/anomalies/anomalies_norm_all_",period,"_",SSP,"_",GCM,".rds")))
    }
  }
  
  # plot(SSP_anomalies_norm_all_stack, col = pal_grn_red)
  
  # Compute median Ensemble
  SSP_anomalies_norm_all_median <- calc(x = SSP_anomalies_norm_all_stack, fun = median, na.rm = T)
  names(SSP_anomalies_norm_all_median) <- "Normalized anomalies - All var"
  
  # plot(SSP_anomalies_norm_all_median, col = pal_grn_red)
  
  # Compute Incertitude (sd)
  SSP_anomalies_norm_all_sd <- calc(x = SSP_anomalies_norm_all_stack, fun = sd, na.rm = T)
  names(SSP_anomalies_norm_all_sd) <- "Normalized anomalies - All var"
  
  # plot(SSP_anomalies_norm_all_sd, col = pal_grn_red)
  
  # Save ensemble median and sd for this SSP
  saveRDS(SSP_anomalies_norm_all_median, file = paste0("./outputs/anomalies/anomalies_norm_all_median_",period,"_",SSP,".rds"), version = "2")
  saveRDS(SSP_anomalies_norm_all_sd, file = paste0("./outputs/anomalies/anomalies_norm_all_sd_",period,"_",SSP,".rds"), version = "2")
  
  # 5.6/ Ensemble for Normalized Anomalies - Only Tmean & Ptot
  SSP_anomalies_norm_no_var_stack <-  NA
  for (i in 1:length(GCM_list)) 
  {
    GCM <- GCM_list[i]
    
    if(i == 1) 
    { 
      SSP_anomalies_norm_no_var_stack <- stack(readRDS(file = paste0("./outputs/anomalies/anomalies_norm_no_var_",period,"_",SSP,"_",GCM,".rds"))) # Create the stack
      
    } else { # Add layer iteratively
      
      SSP_anomalies_norm_no_var_stack <- addLayer(SSP_anomalies_norm_no_var_stack, readRDS(file = paste0("./outputs/anomalies/anomalies_norm_no_var_",period,"_",SSP,"_",GCM,".rds")))
    }
  }
  
  # plot(SSP_anomalies_norm_no_var_stack, col = pal_grn_red)
  
  # Compute median Ensemble
  SSP_anomalies_norm_no_var_median <- calc(x = SSP_anomalies_norm_no_var_stack, fun = median, na.rm = T)
  names(SSP_anomalies_norm_no_var_median) <- "Normalized anomalies - Tmean & Ptot"
  
  # plot(SSP_anomalies_norm_no_var_median, col = pal_grn_red)
  # plot(SSP_anomalies_norm_no_var_median, col = pal_grn_red, zlim = c(0.4,0.85))
  
  # Compute Incertitude (sd)
  SSP_anomalies_norm_no_var_sd <- calc(x = SSP_anomalies_norm_no_var_stack, fun = sd, na.rm = T)
  names(SSP_anomalies_norm_no_var_sd) <- "Normalized anomalies - Tmean & Ptot"
  
  # plot(SSP_anomalies_norm_no_var_sd, col = pal_grn_red)
  
  # Save ensemble median and sd for this SSP
  saveRDS(SSP_anomalies_norm_no_var_median, file = paste0("./outputs/anomalies/anomalies_norm_no_var_median_",period,"_",SSP,".rds"), version = "2")
  saveRDS(SSP_anomalies_norm_no_var_sd, file = paste0("./outputs/anomalies/anomalies_norm_no_var_sd_",period,"_",SSP,".rds"), version = "2")  
  
  
  # 6/ Plot and save maps for this SSP-RCP in a multiple pdf ####
  
  # Make stacks
  SSP_anomalies_median_stack <- stack(SSP_anomalies_base_std_all_median, SSP_anomalies_base_std_no_var_median, SSP_anomalies_base_std_all_median_contrasted, SSP_anomalies_base_std_no_var_median_contrasted, SSP_anomalies_std_all_median, SSP_anomalies_std_no_var_median, SSP_anomalies_std_all_median_contrasted, SSP_anomalies_std_no_var_median_contrasted, SSP_anomalies_norm_all_median, SSP_anomalies_norm_no_var_median)
  SSP_anomalies_sd_stack <- stack(SSP_anomalies_base_std_all_sd, SSP_anomalies_base_std_no_var_sd, SSP_anomalies_base_std_all_sd_contrasted, SSP_anomalies_base_std_no_var_sd_contrasted, SSP_anomalies_std_all_sd, SSP_anomalies_std_no_var_sd, SSP_anomalies_std_all_sd_contrasted, SSP_anomalies_std_no_var_sd_contrasted, SSP_anomalies_norm_all_sd, SSP_anomalies_norm_no_var_sd)
  
  # Save stacks
  saveRDS(SSP_anomalies_median_stack, file = paste0("./outputs/anomalies/SSP_anomalies_median_stack_",period,"_",SSP,".rds"), version = "2")
  saveRDS(SSP_anomalies_sd_stack, file = paste0("./outputs/anomalies/SSP_anomalies_sd_stack_",period,"_",SSP,".rds"), version = "2")  
  
  # Multiple pages plot
  pdf(file = paste0("./maps/anomalies/",SSP,"/global_anomalies_",period,"_",SSP,".pdf"), height = 12, width = 8)
  internal_margins <- par()$mar
  external_margins <- par()$oma
  par(mar = c(5.1,6.1,9.1,3.1))
  par(oma = c(1,2,6,3))
  
  # First page = Ensemble median
  plot(SSP_anomalies_median_stack, col = pal_grn_red, nc = 2)
  title(main = paste0("Multivariate Anomalies\nEnsemble (median)\n",period," ",SSP), outer = T)
  
  # Second page = Incertitudes (sd)
  plot(SSP_anomalies_sd_stack, col = pal_grn_red, nc = 2)
  title(main = paste0("Multivariate Anomalies\nIncertitude (sd)\n",period," ",SSP), outer = T)
  
  par(mar = internal_margins, oma = external_margins)
  dev.off()
  
}


### 7/ Retrieve maps from each SSP-RCP and do a multiple PDF with original scales ####

# # Color palet for plot
# cool = rainbow(50, s = 1, v = 1, start=rgb2hsv(col2rgb('yellow'))[1], end=rgb2hsv(col2rgb('blue'))[1])
# warm = rainbow(50, s = 1, v= 1, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
# pal_bl_red  = c(rev(cool), rev(warm))

# New color palette
pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")

# Load stack for SSP2-RCP4.5
SSP245_anomalies_median_stack <- readRDS(file = paste0("./outputs/anomalies/SSP_anomalies_median_stack_",period,"_ssp245.rds"))
SSP245_anomalies_sd_stack <- readRDS(file = paste0("./outputs/anomalies/SSP_anomalies_sd_stack_",period,"_ssp245.rds"))  

# Load stack for SSP5-RCP8.5
SSP585_anomalies_median_stack <- readRDS(file = paste0("./outputs/anomalies/SSP_anomalies_median_stack_",period,"_ssp585.rds"))
SSP585_anomalies_sd_stack <- readRDS(file = paste0("./outputs/anomalies/SSP_anomalies_sd_stack_",period,"_ssp585.rds"))  


pdf(file = paste0("./maps/anomalies/All_ssp_global_anomalies_",period,".pdf"), height = 12, width = 8)
internal_margins <- par()$mar
external_margins <- par()$oma
par(mar = c(5.1,6.1,9.1,3.1))
par(oma = c(1,2,6,3))

# First page = Ensemble median for SSP245
plot(SSP245_anomalies_median_stack, col = pal_grn_red, nc = 2)
title(main = paste0("Multivariate Anomalies\nEnsemble (median)\n",period," SSP2-RCP4.5"), outer = T)

# Second page = Ensemble median for SSP245
plot(SSP585_anomalies_median_stack, col = pal_grn_red, nc = 2)
title(main = paste0("Multivariate Anomalies\nEnsemble (median)\n",period," SSP5-RCP8.5"), outer = T)

# Third page = Incertitudes (sd) for SSP585
plot(SSP245_anomalies_sd_stack, col = pal_grn_red, nc = 2)
title(main = paste0("Multivariate Anomalies\nIncertitude (sd)\n",period," SSP2-RCP4.5"), outer = T)

# Forth page = Incertitudes (sd) for SSP585
plot(SSP585_anomalies_sd_stack, col = pal_grn_red, nc = 2)
title(main = paste0("Multivariate Anomalies\nIncertitude (sd)\n",period," SSP5-RCP8.5"), outer = T)

par(mar = internal_margins, oma = external_margins)
dev.off()

### Version for 8 GCMs

# 8/ Retrieve maps from each SSP-RCP and do a multiple PDF with a common scale to show differences among SSP-RCP ####

# Only usefull fore the baseline-standardized maps. All the others have fairly similar scales between SSP-RCP

# Get scale for Base std anomalies

SSP245_anomalies_base_std_all_median <- readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_all_median_",period,"_ssp245.rds"))
SSP585_anomalies_base_std_all_median <- readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_all_median_",period,"_ssp585.rds"))

zmin <- min(min(SSP245_anomalies_base_std_all_median[], na.rm = T), min(SSP585_anomalies_base_std_all_median[], na.rm = T))
zmax <- max(max(SSP245_anomalies_base_std_all_median[], na.rm = T), max(SSP585_anomalies_base_std_all_median[], na.rm = T))

hist(c(SSP245_anomalies_base_std_all_median[], SSP585_anomalies_base_std_all_median[]))

hist(SSP245_anomalies_base_std_all_median, col = "dodgerblue3", angle = 135, density = 10,  xlim = c(0.1, 1.2))
hist(SSP585_anomalies_base_std_all_median, col = "red", angle = 45, density = 10, add = T)

# Quick plots
plot(SSP245_anomalies_base_std_all_median, col = pal_grn_red)
plot(SSP245_anomalies_base_std_all_median, col = pal_grn_red, zlim = c(zmin, zmax))
plot(SSP245_anomalies_base_std_all_median, col = pal_grn_red, zlim = c(0.4, 1.1))
plot(SSP245_anomalies_base_std_all_median, col = pal_grn_red, zlim = c(0, 1.1))
plot(SSP245_anomalies_base_std_all_median, col = pal_grn_red, zlim = c(0, zmax))

plot(SSP585_anomalies_base_std_all_median, col = pal_grn_red)
plot(SSP585_anomalies_base_std_all_median, col = pal_grn_red, zlim = c(zmin, zmax))
plot(SSP585_anomalies_base_std_all_median, col = pal_grn_red, zlim = c(0.4, 1.1))
plot(SSP585_anomalies_base_std_all_median, col = pal_grn_red, zlim = c(0, 1.1))
plot(SSP585_anomalies_base_std_all_median, col = pal_grn_red, zlim = c(0, zmax))


# Contrasted maps
SSP245_anomalies_base_std_all_median_contrasted_0.4_1.1 <- calc(x = SSP245_anomalies_base_std_all_median, fun = function(x){ x[x < 0.4] <- 0.4 ; x[x > 1.1] <- 1.1 ; return(x)} )
SSP585_anomalies_base_std_all_median_contrasted_0.4_1.1 <- calc(x = SSP585_anomalies_base_std_all_median, fun = function(x){ x[x < 0.4] <- 0.4 ; x[x > 1.1] <- 1.1 ; return(x)} )
names(SSP245_anomalies_base_std_all_median_contrasted_0.4_1.1) <- names(SSP585_anomalies_base_std_all_median_contrasted_0.4_1.1) <- "Contrasted Scale"

SSP245_anomalies_base_std_all_median_contrasted_0_1.1 <- calc(x = SSP245_anomalies_base_std_all_median, fun = function(x){ x[x > 1.1] <- 1.1 ; return(x)} )
SSP585_anomalies_base_std_all_median_contrasted_0_1.1 <- calc(x = SSP585_anomalies_base_std_all_median, fun = function(x){ x[x > 1.1] <- 1.1 ; return(x)} )

# Save contrasted maps
saveRDS(SSP245_anomalies_base_std_all_median_contrasted_0.4_1.1, file = paste0("./outputs/anomalies/SSP245_anomalies_base_std_all_median_contrasted_0.4_1.1.rds"), version = "2")
saveRDS(SSP585_anomalies_base_std_all_median_contrasted_0.4_1.1, file = paste0("./outputs/anomalies/SSP585_anomalies_base_std_all_median_contrasted_0.4_1.1.rds"), version = "2")
saveRDS(SSP245_anomalies_base_std_all_median_contrasted_0_1.1, file = paste0("./outputs/anomalies/SSP245_anomalies_base_std_all_median_contrasted_0_1.1.rds"), version = "2")
saveRDS(SSP585_anomalies_base_std_all_median_contrasted_0_1.1, file = paste0("./outputs/anomalies/SSP585_anomalies_base_std_all_median_contrasted_0_1.1.rds"), version = "2")

# Load independantly contrasted maps and stack them
SSP245_anomalies_base_std_all_median_contrasted <- readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_all_median_contrasted_",period,"_ssp245.rds"))
SSP585_anomalies_base_std_all_median_contrasted <- readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_all_median_contrasted_",period,"_ssp585.rds"))
SSP_anomalies_base_std_all_median_indep_contrasted <- stack(SSP245_anomalies_base_std_all_median_contrasted, SSP585_anomalies_base_std_all_median_contrasted)
names(SSP_anomalies_base_std_all_median_indep_contrasted) <- c("SSP2-RCP4.5", "SSP5-RCP8.5")

# Stack stuff per scaling
SSP_anomalies_base_std_all_median_stack <- stack(SSP245_anomalies_base_std_all_median, SSP585_anomalies_base_std_all_median)
SSP_anomalies_base_std_all_median_contrasted_0.4_1.1_stack <- stack(SSP245_anomalies_base_std_all_median_contrasted_0.4_1.1, SSP585_anomalies_base_std_all_median_contrasted_0.4_1.1)
SSP_anomalies_base_std_all_median_contrasted_0_1.1_stack <- stack(SSP245_anomalies_base_std_all_median_contrasted_0_1.1, SSP585_anomalies_base_std_all_median_contrasted_0_1.1)
names(SSP_anomalies_base_std_all_median_stack) <- names(SSP_anomalies_base_std_all_median_contrasted_0.4_1.1_stack) <- names(SSP_anomalies_base_std_all_median_contrasted_0_1.1_stack) <- c("SSP2-RCP4.5", "SSP5-RCP8.5")

# Make a multiple pages PDF. One page = one scale.

pdf(file = paste0("./maps/anomalies/All_ssp_multiscales_base_std_anomalies_",period,".pdf"), height = 5, width = 10)
internal_margins <- par()$mar
external_margins <- par()$oma
par(mar = c(5.1,4.1,4.1,2.1))
par(oma = c(0.5,0.5,3,0.5))

# First page = Free Scale
plot(SSP_anomalies_base_std_all_median_stack, col = pal_grn_red, nc = 2)
title(main = paste0("\nMultivariate baseline-std anomalies\n",period,"\nFree scale"), outer = T)

# Second page = Independantly contrasted
plot(SSP_anomalies_base_std_all_median_indep_contrasted, col = pal_grn_red, nc = 2)
title(main = paste0("\nMultivariate baseline-std anomalies\n",period,"\nIndep. contrasted"), outer = T)

# Third page = All included (min-max)
plot(SSP_anomalies_base_std_all_median_stack, col = pal_grn_red, nc = 2, zlim = c(zmin, zmax))
title(main = paste0("\nMultivariate baseline-std anomalies\n",period,"\nAll included (min-max)"), outer = T)

# Fourth page = Contrasted
plot(SSP_anomalies_base_std_all_median_contrasted_0.4_1.1_stack, col = pal_grn_red, nc = 2, zlim = c(0.4, 1.1))
title(main = paste0("\nMultivariate baseline-std anomalies\n",period,"\nContrasted"), outer = T)

# Fifh page = Zero to contrasted max
plot(SSP_anomalies_base_std_all_median_contrasted_0_1.1_stack, col = pal_grn_red, nc = 2, zlim = c(0, 1.1))
title(main = paste0("\nMultivariate baseline-std anomalies\n",period,"\nZero to contrasted max"), outer = T)

# Sixth page = Zero to max
plot(SSP_anomalies_base_std_all_median_stack, col = pal_grn_red, nc = 2, zlim = c(0, zmax))
title(main = paste0("\nMultivariate baseline-std anomalies\n",period,"\nZero to max"), outer = T)


par(mar = internal_margins, oma = external_margins)
dev.off()


### Version for 6 GCMs

# 9/ Retrieve maps from each SSP-RCP and do a multiple PDF with a common scale to show differences among SSP-RCP ####

# Only usefull fore the baseline-standardized maps. All the others have fairly similar scales between SSP-RCP

# Get scale for Base std anomalies

SSP245_anomalies_base_std_all_median <- readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_all_median_",period,"_ssp245.rds"))
SSP585_anomalies_base_std_all_median <- readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_all_median_",period,"_ssp585.rds"))

zmin <- min(min(SSP245_anomalies_base_std_all_median[], na.rm = T), min(SSP585_anomalies_base_std_all_median[], na.rm = T))
zmax <- max(max(SSP245_anomalies_base_std_all_median[], na.rm = T), max(SSP585_anomalies_base_std_all_median[], na.rm = T))

hist(c(SSP245_anomalies_base_std_all_median[], SSP585_anomalies_base_std_all_median[]))

hist(SSP245_anomalies_base_std_all_median, col = "dodgerblue3", angle = 135, density = 10,  xlim = c(0.1, 1.2))
hist(SSP585_anomalies_base_std_all_median, col = "red", angle = 45, density = 10, add = T)

# Quick plots
plot(SSP245_anomalies_base_std_all_median, col = pal_grn_red)
plot(SSP245_anomalies_base_std_all_median, col = pal_grn_red, zlim = c(zmin, zmax))
plot(SSP245_anomalies_base_std_all_median, col = pal_grn_red, zlim = c(0.4, 1.15))
plot(SSP245_anomalies_base_std_all_median, col = pal_grn_red, zlim = c(0, 1.15))
plot(SSP245_anomalies_base_std_all_median, col = pal_grn_red, zlim = c(0, zmax))

plot(SSP585_anomalies_base_std_all_median, col = pal_grn_red)
plot(SSP585_anomalies_base_std_all_median, col = pal_grn_red, zlim = c(zmin, zmax))
plot(SSP585_anomalies_base_std_all_median, col = pal_grn_red, zlim = c(0.4, 1.15))
plot(SSP585_anomalies_base_std_all_median, col = pal_grn_red, zlim = c(0, 1.15))
plot(SSP585_anomalies_base_std_all_median, col = pal_grn_red, zlim = c(0, zmax))


# Contrasted maps
SSP245_anomalies_base_std_all_median_contrasted_0.4_1.15 <- calc(x = SSP245_anomalies_base_std_all_median, fun = function(x){ x[x < 0.4] <- 0.4 ; x[x > 1.15] <- 1.15 ; return(x)} )
SSP585_anomalies_base_std_all_median_contrasted_0.4_1.15 <- calc(x = SSP585_anomalies_base_std_all_median, fun = function(x){ x[x < 0.4] <- 0.4 ; x[x > 1.15] <- 1.15 ; return(x)} )
names(SSP245_anomalies_base_std_all_median_contrasted_0.4_1.15) <- names(SSP585_anomalies_base_std_all_median_contrasted_0.4_1.1) <- "Contrasted Scale"

SSP245_anomalies_base_std_all_median_contrasted_0_1.15 <- calc(x = SSP245_anomalies_base_std_all_median, fun = function(x){ x[x > 1.15] <- 1.15 ; return(x)} )
SSP585_anomalies_base_std_all_median_contrasted_0_1.15 <- calc(x = SSP585_anomalies_base_std_all_median, fun = function(x){ x[x > 1.15] <- 1.15 ; return(x)} )

# Save contrasted maps
saveRDS(SSP245_anomalies_base_std_all_median_contrasted_0.4_1.15, file = paste0("./outputs/anomalies/SSP245_anomalies_base_std_all_median_contrasted_0.4_1.15.rds"), version = "2")
saveRDS(SSP585_anomalies_base_std_all_median_contrasted_0.4_1.15, file = paste0("./outputs/anomalies/SSP585_anomalies_base_std_all_median_contrasted_0.4_1.15.rds"), version = "2")
saveRDS(SSP245_anomalies_base_std_all_median_contrasted_0_1.15, file = paste0("./outputs/anomalies/SSP245_anomalies_base_std_all_median_contrasted_0_1.15.rds"), version = "2")
saveRDS(SSP585_anomalies_base_std_all_median_contrasted_0_1.15, file = paste0("./outputs/anomalies/SSP585_anomalies_base_std_all_median_contrasted_0_1.15.rds"), version = "2")

# Load independantly contrasted maps and stack them
SSP245_anomalies_base_std_all_median_contrasted <- readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_all_median_contrasted_",period,"_ssp245.rds"))
SSP585_anomalies_base_std_all_median_contrasted <- readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_all_median_contrasted_",period,"_ssp585.rds"))
SSP_anomalies_base_std_all_median_indep_contrasted <- stack(SSP245_anomalies_base_std_all_median_contrasted, SSP585_anomalies_base_std_all_median_contrasted)
names(SSP_anomalies_base_std_all_median_indep_contrasted) <- c("SSP2-RCP4.5", "SSP5-RCP8.5")

# Stack stuff per scaling
SSP_anomalies_base_std_all_median_stack <- stack(SSP245_anomalies_base_std_all_median, SSP585_anomalies_base_std_all_median)
SSP_anomalies_base_std_all_median_contrasted_0.4_1.15_stack <- stack(SSP245_anomalies_base_std_all_median_contrasted_0.4_1.15, SSP585_anomalies_base_std_all_median_contrasted_0.4_1.15)
SSP_anomalies_base_std_all_median_contrasted_0_1.15_stack <- stack(SSP245_anomalies_base_std_all_median_contrasted_0_1.15, SSP585_anomalies_base_std_all_median_contrasted_0_1.15)
names(SSP_anomalies_base_std_all_median_stack) <- names(SSP_anomalies_base_std_all_median_contrasted_0.4_1.15_stack) <- names(SSP_anomalies_base_std_all_median_contrasted_0_1.15_stack) <- c("SSP2-RCP4.5", "SSP5-RCP8.5")

# Make a multiple pages PDF. One page = one scale.

pdf(file = paste0("./maps/anomalies/All_ssp_multiscales_base_std_anomalies_",period,".pdf"), height = 5, width = 10)
internal_margins <- par()$mar
external_margins <- par()$oma
par(mar = c(5.1,4.1,4.1,2.1))
par(oma = c(0.5,0.5,3,0.5))

# First page = Free Scale
plot(SSP_anomalies_base_std_all_median_stack, col = pal_grn_red, nc = 2)
title(main = paste0("\nMultivariate baseline-std anomalies\n",period,"\nFree scale"), outer = T)

# Second page = Independantly contrasted
plot(SSP_anomalies_base_std_all_median_indep_contrasted, col = pal_grn_red, nc = 2)
title(main = paste0("\nMultivariate baseline-std anomalies\n",period,"\nIndep. contrasted"), outer = T)

# Third page = All included (min-max)
plot(SSP_anomalies_base_std_all_median_stack, col = pal_grn_red, nc = 2, zlim = c(zmin, zmax))
title(main = paste0("\nMultivariate baseline-std anomalies\n",period,"\nAll included (min-max)"), outer = T)

# Fourth page = Contrasted
plot(SSP_anomalies_base_std_all_median_contrasted_0.4_1.15_stack, col = pal_grn_red, nc = 2, zlim = c(0.4, 1.15))
title(main = paste0("\nMultivariate baseline-std anomalies\n",period,"\nContrasted"), outer = T)

# Fifh page = Zero to contrasted max
plot(SSP_anomalies_base_std_all_median_contrasted_0_1.15_stack, col = pal_grn_red, nc = 2, zlim = c(0, 1.15))
title(main = paste0("\nMultivariate baseline-std anomalies\n",period,"\nZero to contrasted max"), outer = T)

# Sixth page = Zero to max
plot(SSP_anomalies_base_std_all_median_stack, col = pal_grn_red, nc = 2, zlim = c(0, zmax))
title(main = paste0("\nMultivariate baseline-std anomalies\n",period,"\nZero to max"), outer = T)


par(mar = internal_margins, oma = external_margins)
dev.off()


# 10/ Do a contrasted map of the Andes ####

# Color palet for diversity indices
# cool = rainbow(49, s = 1, v = 1, start=rgb2hsv(col2rgb('yellow'))[1], end=rgb2hsv(col2rgb('blue'))[1])
# warm = rainbow(50, s = 1, v= 1, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
# pal_bl_red  = c(rev(cool), rev(warm))
# pal_bl_red <- c(gplots::col2hex("grey93"), pal_bl_red)

# New color palet for diversity indices
pal_bl_red_Mannion <- readRDS(file = "./maps/pal_bl_red_Mannion.rds")

library(rangeBuilder)

# Load independantly contrasted maps and stack them
SSP245_anomalies_base_std_all_median_contrasted <- readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_all_median_contrasted_",period,"_ssp245.rds"))
SSP585_anomalies_base_std_all_median_contrasted <- readRDS(file = paste0("./outputs/anomalies/anomalies_base_std_all_median_contrasted_",period,"_ssp585.rds"))
SSP_anomalies_base_std_all_median_indep_contrasted <- stack(SSP245_anomalies_base_std_all_median_contrasted, SSP585_anomalies_base_std_all_median_contrasted)
names(SSP_anomalies_base_std_all_median_indep_contrasted) <- c("SSP2-RCP4.5", "SSP5-RCP8.5")

### Load environmental stack to use as template for resolution and extent
res <- "15"
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))

# Get elevation contours
Andes_ext <- extent(c(-90, -59, -15, 14))
Close_Andes_ext <- extent(c(-90, -67.2, -15, 14))
DEM <- envData[["Elevation"]]
Andes_DEM <- crop(DEM, Andes_ext)
Close_Andes_DEM <- crop(DEM, Close_Andes_ext)

# Function to use in stack plots. Issue = cannot change palette between plots
contour_fun <- function() {
  contour(x = Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1),
          drawlabels = F, col = "black", add = T)
}
plot(Andes_SSP_anomalies_base_std_all_median_indep_contrasted, col = pal_grn_red, addfun = contour_fun)

# Crop maps and remove pixels not computed for indices maps
Andes_SSP_anomalies_base_std_all_median_indep_contrasted <- crop(SSP_anomalies_base_std_all_median_indep_contrasted, Andes_ext)
Andes_SSP_anomalies_base_std_all_median_indep_contrasted <- mask(x = Andes_SSP_anomalies_base_std_all_median_indep_contrasted, mask = Andes_DEM)

crop_mask_shp <- readRDS(file = paste0("./input_data/Env_data/crop_mask_shp_", res, ".rds"))

# Load and crop Indices maps
tot.sp.richness_Jaccard.80 <- readRDS(file = paste0("./outputs/Indices_maps/tot.sp.richness_Jaccard.80.rds"))
MPD.raster_Jaccard.80_contrasted <- readRDS(file = "./outputs/Indices_Maps/MPD.raster_Jaccard.80_contrasted.rds")

Andes.tot.sp.richness_Jaccard.80 <- crop(tot.sp.richness_Jaccard.80, Andes_ext)
Andes.MPD.raster_Jaccard.80_contrasted <- crop(MPD.raster_Jaccard.80_contrasted, Andes_ext)

# Save plot
pdf(file = paste0("./maps/anomalies/Andes_zoom_",period,".pdf"), height = 12, width = 12)

par(mfrow=(c(2,2)))

image(Andes.tot.sp.richness_Jaccard.80, col = pal_bl_red_Mannion, legend = T,
      ylab = "", xlab = "", main = "Species richness",
      cex.axis = 1.4, cex.main = 1.6)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
        drawlabels = F, col = "black", add = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(Andes.tot.sp.richness_Jaccard.80, locs = seq(20, 120, 20), cex.axis = 1.1, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.5, 4.5))
graphics::text(x = -86.7, y = 6, font = 2, cex = 1.1, label = "Species")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 1.1, text.font = 2)
graphics::text(x = -86.7, y = -10.5, font = 2, cex = 1.1, label = "Elevation")

image(Andes.MPD.raster_Jaccard.80_contrasted, col = pal_bl_red_Mannion, legend = T,
      ylab = "", xlab = "", main = "Mean pairwise Phylogenetic Distance",
      cex.axis = 1.4, cex.main = 1.6)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
        drawlabels = F, col = "black", add = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(Andes.MPD.raster_Jaccard.80_contrasted, locs = seq(38, 44, 2), cex.axis = 1.1, ramp = pal_bl_red_Mannion, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.5, 4.5))
graphics::text(x = -86.5, y = 6.5, font = 2, cex = 1.1, label = "Evolutionary\nTime (My)")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 1.1, text.font = 2)
graphics::text(x = -86.7, y = -10.5, font = 2, cex = 1.1, label = "Elevation")

image(Andes_SSP_anomalies_base_std_all_median_indep_contrasted[[1]], col = pal_grn_red, legend = T,
      ylab = "", xlab = "", main = "SSP2 - RCP4.5",
      cex.axis = 1.4, cex.main = 1.6)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
        drawlabels = F, col = "black", add = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(Andes_SSP_anomalies_base_std_all_median_indep_contrasted[[1]], locs = seq(0.45, 0.7, 0.05), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.5, 4.5))
graphics::text(x = -86.5, y = 6.5, font = 2, cex = 1.1, label = "Climate\nanomalies")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 1.1, text.font = 2)
graphics::text(x = -86.7, y = -10.5, font = 2, cex = 1.1, label = "Elevation")

image(Andes_SSP_anomalies_base_std_all_median_indep_contrasted[[2]], col = pal_grn_red, legend = T,
      ylab = "", xlab = "", main = "SSP5 - RCP8.5",
      cex.axis = 1.4, cex.main = 1.6)
contour(x = Close_Andes_DEM, levels = c(500, 2500), lwd = 1.5, lty = c(2, 1) ,
        drawlabels = F, col = "black", add = T)
plot(crop_mask_shp, lwd = 1.2, border = "grey20", add = T)
addRasterLegend(Andes_SSP_anomalies_base_std_all_median_indep_contrasted[[2]], locs = seq(0.6, 1.1, 0.1), cex.axis = 1.1, ramp = pal_grn_red, ncolors = 200, border = T, location = c(-88.22, -87.2, -8.5, 4.5))
graphics::text(x = -86.5, y = 6.5, font = 2, cex = 1.1, label = "Climate\nanomalies")
legend(x = "bottomleft", title = "", legend = c("2500m", "500m"), lty = c(1, 2) , lwd = 2, bty = "n", cex = 1.1, text.font = 2)
graphics::text(x = -86.7, y = -10.5, font = 2, cex = 1.1, label = "Elevation")

dev.off()

par(mfrow=(c(1,1)))

### Use tmap ?
# Issue = Legend is not pretty, but could probably be improved

isolines <- rasterToContour(x = Andes_DEM, levels = c(500, 2500)) # Values to generate isolines

library(tmap)

tm_shape(Andes_SSP_anomalies_base_std_all_median_indep_contrasted[[1]]) +
  tm_raster(title = "Anomalies", palette = pal_grn_red, n = 5, style = "cont", legend.reverse = T) +
  tm_legend(legend.position = c("left", "bottom"), legend.frame = T) +
  tm_shape(isolines) + tm_iso(group = "Isolines")

?tm_raster
?tm_legend

# To do a multiple plot with tmap
library(grid)
grid.newpage()
pushViewport(viewport(layout=grid.layout(1,2)))
print(plot1, vp=viewport(layout.pos.col = 1))
print(plot1, vp=viewport(layout.pos.col = 2))

