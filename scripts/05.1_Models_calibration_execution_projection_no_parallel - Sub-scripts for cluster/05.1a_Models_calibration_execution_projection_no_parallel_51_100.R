

##### Script 5: Models Calibration, Execution, and Projection ######

##### Version not parallelized (script for parallelization in comments)

# Remove environment
rm(list = ls())

# setwd("D:/Mael/R_projects/ithomiini_current/")

# ### Set up the cluster ####
# 
# library(foreach)
# library(doParallel)
# 
# # Detect the number of threads
# numCores <- detectCores()
# # registerDoParallel(numCores)  # To set the nb of threads to use when using %dopar%
# 
# cl <- makeCluster(numCores) # Generate the virtual cluster from the number of threads
# registerDoParallel(cl)      # To set the virtual cluster when using %dopar%
# 
# # getDoParWorkers() # Pour vérifier le nombre de coeurs enregistrés


### Load stuff useful for all units (outside the loop) ####

### Select Env data resolution 
# res <- "5"
res <- "15"

# Load environmental stack
library(raster)
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))

# Load Full unit list
load(file = paste0("./input_data/list.models.RData"))

# Generate new column in the summary table for units to keep tracks of models results

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
 index_model_to_compute <- c(51:100)

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

# Set seed to ensure repetability of random PsA draws
set.seed(158340)

source("./functions/progcombine.R")

k <- 1 # Index to keep track of number of iterations when not parallelized, and run on limited nb of units
for (i in index_model_to_compute) { # Unparallelized version
# # Parallelized version
# temp <- foreach (i = index_model_to_compute, .combine = progcombine_rbind(nreps = length(index_model_to_compute)), .init = NULL,
#                  .inorder = F, # Improve performance. Anyway, using the index stored in the final object, we can retrieve the order afterwards
#                  .errorhandling = "remove", # If an error occurs during one task/instance, do not stop, just do not provide an output for .combine
#                  .packages = c("biomod2", "raster", "reshape2")) %dopar% { 
  
  # i <- 4
  # i <- 6
  
  # Load unit/OMU name
  unit <- as.character(list.models$Tag.model[i]) 
  model_type <- list.models$initial_model_type[i]
  
  if (model_type == "rasterized") {
    
    cat(paste("\n", Sys.time(),"------ NO MODELING for", unit, "because of low N.obs \n")) # Warning message for OMU with no SDM
    
    # Save NA values in the object to retrieve at each iteration: Parallelized version
    # Start_Time <- End_Time <- Model_ID <- Models_success <- Mean_TSS_full <- Models_select_TSS_0.5 <- Mean_TSS_select_TSS_0.5 <- NA
    
  } else {
    
    cat(paste0("\n", Sys.time()," ------ Starting process for ", unit, " = Unit N°",k,"/",length(index_model_to_compute)," ------\n")) # Affichage le temps au début de la procédure
    
    # Unparallelized version
     list.models$Start_Time[list.models$Tag.model == unit] <- as.character(Sys.time()) # Store date and time of starting of last modeling process
    
    # Parallelized version
    # Start_Time <- as.character(Sys.time())
    
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
                                                expl.var = unit.env,            # Env stack
                                                resp.name = unit,              # Unit/OMU name
                                                PA.nb.rep = ncol(PA.table),    # Nb of PsA sets
                                                PA.nb.absences = sum(unit.points.PsA@data, na.rm = T), # Nb of occurrences = Nb of PsA per set
                                                PA.strategy = "user.defined",  # To allow to provide a custom PA_table
                                                PA.table = PA.table)            # Provide the custom PA_table
    
    save(formated.input.data, file = paste0(initial.wd,"/models/",res,"/",unit,"/formated.input.data.RData"))
    saveRDS(formated.input.data, file = paste0(initial.wd,"/models/",res,"/",unit,"/formated.input.data.rds"))
    
    # Plot allow to check for the location of data in each PA set
    # x11()
    # plot(formated.input.data) 
    
    ################# 2/ # Models parametrization and execution ######################
    
    # Starting warning
    cat(paste0("\n", Sys.time(),"------ Modeling for ", unit, " = Unit N°",k,"/",length(index_model_to_compute)," ------\n"))
    
    # Check what is going on with the do.full or DataSplit = 100 pour l'évaluation par la suite et la gestion du nom des modèles
    # Faire un if pour 2 setup différents !
    
    # DataSplitTable erase NbRunEva, DataSplit et do.full
    
    # Version complete
    
    if (model_type == "complete") { # Modeling for OMU with enough obs (N.obs ≥ 30) # 3 PsA sets # Spatial CV with 3 folds 
      
      # Load CV array for biomod2, only for OMU with enough obs (N.obs ≥ 30)
      calib_array <- readRDS(file = paste0("../../input_data/Species_data/",res,"/CV_blocks/calib_arrays/calib_array_", unit,".rds"))
      
      model.runs <- BIOMOD_Modeling(data = formated.input.data,
                                    models =  c('RF','GBM','ANN'), # Algorithm: RF = Random Forest, GBM = Gradient Boosting Model, ANN = Artifical Neural Network
                                    DataSplitTable = calib_array, # Array describing CV blocks for calibration/evaluation
                                    VarImport = 3, # Nb of data permutation to evaluate variable importance
                                    models.eval.meth = c("TSS"), # Evaluation metric = True Skill Statistic
                                    do.full.models = FALSE, # No calibration with 100% data
                                    SaveObj = T, # Save everything generated during the modeling process
                                    Prevalence = NULL, # Weights P/A. Poids relatifs d'une absence par rapport à une présence. Si 0.5 => poids répartis également selon le nb de P/A => somme des poids des P et A est identiques. Lorsque qu'utilisation de pseudo-absences, PR et PA équivalentes par défaut => Prévalence = NULL (= 0.5 par défaut). Si Prévalence > 0.5 => more weights for PR. If Prévalence < 0.5 => more weights for Absences.                                
                                    rescal.all.models = FALSE, # No rescaling of output needed with the last version
                                    modeling.id = , # Empty = generate a new ID each time
                                    models.options = BIOMOD_ModelingOptions( # All default options for algorithm...
                                      GBM = list(n.minobsinnode = 2)) # ... except minimum nb of obs per leaf in GBM (n = 2)
      )
      
    } else {  # Model type = "restricted" # 10 PsA sets # Calibration and evaluation on 100% of data
      
      model.runs <- BIOMOD_Modeling(data = formated.input.data,
                                    models =  c('RF','GBM','ANN'), # Algorithm: RF = Random Forest, GBM = Gradient Boosting Model, ANN = Artifical Neural Network
                                    DataSplit = 100, # Proportion of data to use for calibration
                                    NbRunEval = 1, # Nb of evaluation run = 1 because we use 100% of data for calibration. No CV possible for low sample size (N < 30)
                                    VarImport = 3, # Nb of data permutation to evaluate variable importance
                                    models.eval.meth = c("TSS"), # Evaluation metric = True Skill Statistic
                                    do.full.models = TRUE, # Compute model with all data. RUN = "Full"
                                    SaveObj = T, # Save everything generated during the modeling process
                                    Prevalence = NULL, # Weights P/A. Poids relatifs d'une absence par rapport à une présence. Si 0.5 => poids répartis également selon le nb de P/A => somme des poids des P et A est identiques. Lorsque qu'utilisation de pseudo-absences, PR et PA équivalentes par défaut => Prévalence = NULL (= 0.5 par défaut). Si Prévalence > 0.5 => more weights for PR. If Prévalence < 0.5 => more weights for Absences.                                
                                    rescal.all.models = FALSE, # No rescaling of output needed with the last version
                                    modeling.id = , # Empty = generate a new ID each time
                                    models.options = BIOMOD_ModelingOptions( # All default options for algorithm...
                                      GBM = list(n.minobsinnode = 2, bag.fraction = 0.7)) # ... except minimum nb of obs per leaf (= 2), and fraction of obs to select at each iteration to fit a new tree (= 0.7) in GBM (n = 2)
      )
      
    }
    
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
    
    # # Store infos from modeling results: Parallelized version
    # Model_ID <- model.runs@modeling.id # Keep model ID
    # Models_success <- length(get_built_models(model.runs)) # Save nb of converging models
    
    
    # Project everything, but do not do the Ensemble. We need to compute Jaccard and select proper models that pass our evaluatio check prior to do the Ensembel manually (also allow to compute sd)
    
    
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
      
      # # Store info on Biomod2 Evaluation with TSS: Parallelized version
      # Mean_TSS_full <- round(mean(eval.runs$value, na.rm = T),3) # Mean value for TSS among all converged models
      # Models_select_TSS_0.5 <- sum(eval.runs$value >= 0.5, na.rm = T) # Nb of model with TSS ≥ 0.5
      # if (length(eval.runs$value[eval.runs$value >= 0.5]) > 0) {
      #   Mean_TSS_select_TSS_0.5 <- round(mean(eval.runs$value[eval.runs$value >= 0.5], na.rm = T),3) # Mean value for TSS among all models with TSS ≥ 0.5
      # } else {
      #   Mean_TSS_select_TSS_0.5 <- NA # To avoid error when there are no models with TSS ≥ 0.5
      # }
      
      
      ########## 3/ Models projection ##############
      
      cat(paste0("\n", Sys.time()," ------ Projection of all models for ", unit, " = Unit N°",k,"/",length(index_model_to_compute)," ------\n"))
      
      # Projection of each "model" = RUN * PA * Algo
      projection.runs <- BIOMOD_Projection(modeling.output = model.runs, # Output of BIOMOD_Modeling
                                           new.env = stack(envData), # Environnmental stack to use for the projection
                                           proj.name = "Current", # Name of the projection based on name of the climatic model
                                           selected.models = 'all', # Use all models
                                           binary.meth = NULL, # No thresholding for binary maps. We will do it outide of Biomod2
                                           omit.na = T, # Do not predict pixel with incomplete env data
                                           on_0_1000 = T, # Predictions store on 0 to 1000 scale for LOT OF memory saving
                                           compress = TRUE, # Compress output files
                                           build.clamping.mask = TRUE, # To build a mask showing which predictions are made out of the environmentalrange out of the data used for calibraiton => i.e., show the area of extrapolation
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
    
    # Store date and time of the end of modeling process: parallelized version
    # End_Time <- as.character(Sys.time())
    
    # print(warnings()) # To control potential errors
      
  }    
  
  # For unparallized version: Save summary table in specific file at each iteration to not lose info if crashing
  save(list.models, file = paste0(initial.wd,"/models/",res,"/0_List.models_outputs/list.models_",min(index_model_to_compute),"_to_",max(index_model_to_compute),".RData")) 
  saveRDS(list.models, file = paste0(initial.wd,"/models/",res,"/0_List.models_outputs/list.models_",min(index_model_to_compute),"_to_",max(index_model_to_compute),".rds"))
  
   
  # Set wd back to initial working directory
  setwd(initial.wd)
  
  # # Variables to retrieve for each instance in parallelized computing
  # return(list(index = i, Start_Time = Start_Time, End_Time = End_Time, Model_ID = Model_ID, 
  #             Models_success = Models_success, Mean_TSS_full = Mean_TSS_full, 
  #             Models_select_TSS_0.5 = Models_select_TSS_0.5, Mean_TSS_select_TSS_0.5 = Mean_TSS_select_TSS_0.5))

  cat(paste0("\n", Sys.time()," ------ Process over for ", unit, " = Unit N°",k,"/",length(index_model_to_compute)," ------\n")) # Affichage le temps au début de la procédure
  
  k <- k + 1
  
}

setwd(initial.wd) # Pour reset le working directory dans le PC


# Unparallelized version
save(list.models, file = paste0(initial.wd,"/models/",res,"/0_List.models_outputs/list.models_",min(index_model_to_compute),"_to_",max(index_model_to_compute),".RData")) 
saveRDS(list.models, file = paste0(initial.wd,"/models/",res,"/0_List.models_outputs/list.models_",min(index_model_to_compute),"_to_",max(index_model_to_compute),".rds"))
# save(list.models, file = "./input_data/list.models.RData")


# Parallelized version

# ### Manage outputs from the parallelized runs and save them in the summary table: list.models
# 
# # For whatever reason, the reshaping of the output as a row.bind df does not work, so need to do it manually afterwards
# library(dplyr)
# library(tidyr)
# 
# temp <- temp %>% 
#   dplyr::as_tibble() %>%  # Get the list arrange in a df
#   tidyr::unnest() %>%     # Change columns from list to vector of the proper type
#   dplyr::arrange(index)   # Sort in index order
# 
# # Save the info in case the addition into list.models fails or running on the cluster (not a good idea to reand/write the same list.model object in all instances)
# saveRDS(temp, file = paste0("./models/",res,"/0_temp_outputs/temp_",min(temp$index, na.rm = T),"_to_",max(temp$index, na.rm = T),".rds"))

# # Actualize list.model with full temp table
# library(tidyverse)
#
# list.models %>% 
#   select(-intersect(names(.), names(temp))) %>% 
#   bind_cols(., temp)
# 
# # Insert results in the summary table unit by unit
# for (i in seq_along(temp)) {
#   list.models$Start_Time[temp$index[i]] <- temp$Start_Time[i]
#   list.models$End_Time[temp$index[i]] <- temp$End_Time[i]
#   list.models$Model_ID[temp$index[i]] <- temp$Model_ID[i]
#   list.models$Models_success[temp$index[i]] <- temp$Models_success[i]
#   list.models$Mean_TSS_full[temp$index[i]] <- temp$Mean_TSS_full[i]
#   list.models$Models_select_TSS_0.5[temp$index[i]] <- temp$Models_select_TSS_0.5[i]
#   list.models$Mean_TSS_select_TSS_0.5[temp$index[i]] <- temp$Mean_TSS_select_TSS_0.5[i]
# }
# 
# View(list.models)

# save(list.models, file = "./input_data/list.models.RData")

### Close cluster
# # Close clusters generated with registerDoParallel(numCores)
# # stopImplicitCluster()
# # Close clusters generated with makeCluster 
# stopCluster() 


# ##### Explore model projection outputs #####
# 
# modeling.wd <- paste0(initial.wd, "/models/")
# 
# # Choice of unit
# load(file = paste0(initial.wd,"./input_data/list.models.RData"))
# i <- 8
# unit <- list.models$Tag.model[i]
# 
# # Choice of resolution
# res <- "15"
# 
# 
# # Loading occurrence data
# load(paste0(initial.wd,"/input_data/Species_data/",res,"/Spatial_Points_Objects/occurrences_", unit,".RData")) # 
# 
# 
# 
# # Loading projection stack
# proj_stack <- readRDS(file = paste0(modeling.wd,res,"/",unit,"/proj_Current/proj_Current_",unit,".rds"))
# 
# plot(proj_stack)
# points(unit.points, pch = 16)
# 
# 
# 
# # Loading Clamping mask
# Clamping_mask <- readRDS(file = paste0(initial.wd,"/models/",res,"/",unit,"/proj_Current/proj_Current_",unit,"_ClampingMask.rds"))
# 
# plot(Clamping_mask)
# points(unit.points, pch = 16)
