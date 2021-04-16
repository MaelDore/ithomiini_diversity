
##### Script 10: Ensemble #####

# Generate stacks of all submodels for each OMU/unit
# Compute EM: median, incertitude (sd), binary, and committee averaging
# Update boxplot of evaluations with only the selected submodels, and the EM evaluation

##### 

# Inputs:
   # Summary table for OMU models (list.models)
   # Submodel maps from Script 05.1
   # List of selected models for ensemble from Script 08
   # Optimized thresholds for binarization from Script 07

# Outputs:
   # Multiple pages PDF with all selected submodel maps per OMU
   # Ensemble maps of each OMU: median, incertitude (sd), binary, and committee averaging
   # Thresholds and optimized evaluation metric for Ensemble
   # Updated boxplots of evaluation metrics per OMU, and globally, with only the selected submodels, and the EM evaluation

#####


#### Load stuff to prepare the loop

# Effacer l'environnement
rm(list = ls())

# setwd("D:/Mael/R_projects/ithomiini_current/")

# Load libraries

library(raster)
library(biomod2)
library(tidyverse)

# Select Env data resolution 
# res <- "5"
res <- "15"

# Load Summary table for OMU models
load(file = paste0("./input_data/list.models.RData"))
modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]

# Load summary table for sub_models
list.submodels <- readRDS(file = "./input_data/list.submodels.rds")

# Choose the eval metric for selection
valid_submodels_names_Jaccard <- list.submodels$model_name[list.submodels$Ensemble_OK_Jaccard] # Jaccard
valid_submodels_names_TSS <- list.submodels$model_name[list.submodels$Ensemble_OK_TSS] # TSS

# Load the cut-off data
all.units_cutoffs <- readRDS(file = "./controls/evaluations/all.units_cutoffs.rds")

# load eval data to update boxplots with EM metric and only selected submodels
gg_jaccard <- readRDS(file = "./controls/evaluations/Full_df_Jaccard.rds")
gg_TSS <- readRDS(file = "./controls/evaluations/Full_df_TSS.rds")

# Set up directory
initial.wd <- getwd()
setwd(paste0(initial.wd,"/models/",res,"/"))

# Create object to keep track of cutoff research
jaccard.test.list_EM <- TSS.test.list_EM <- NULL

############## Modif between scripts ######################
 index_model_to_compute <- c(1:563)
# index_model_to_compute <- c(1:200)
# index_model_to_compute <- c(201:400)
# index_model_to_compute <- c(401:563)
###########################################################

### Loop for all modeled OMU/unit
# for (i in 1:nrow(modeled_OMU)) 
for (i in index_model_to_compute) 
{
  # i <- 1
  
  unit <- as.character(modeled_OMU$Tag.model[i])
  sp <- as.character(modeled_OMU$Sp_full[i])
  
  cat(paste0("\n", Sys.time()," ------ Starts for ", unit, " = Unit N°",i,"/",nrow(modeled_OMU)," ------\n"))
  
  # Create directory to store stacks by OMU
  if(!dir.exists(paste0("../../outputs/By_OMU/",unit))) { # Test if the folder exists already or not
    dir.create(paste0("../../outputs/By_OMU/",unit), recursive = T) # Create folder if absent
  }
  
  # Create directory to store stacks by species
  if(!dir.exists(paste0("../../outputs/By_Species/",sp))) { # Test if the folder exists already or not
    dir.create(paste0("../../outputs/By_Species/",sp), recursive = T) # Create folder if absent
  }
  
  # Create directory to store maps by OMU
  if(!dir.exists(paste0("../../maps/By_OMU/",unit))) { # Test if the folder exists already or not
    dir.create(paste0("../../maps/By_OMU/",unit), recursive = T) # Create folder if absent
  }
  
  # Create directory to store maps by species
  if(!dir.exists(paste0("../../maps/By_Species/",sp))) { # Test if the folder exists already or not
    dir.create(paste0("../../maps/By_Species/",sp), recursive = T) # Create folder if absent
  }
  
  ### 1/ Generate stacks of selected submodels ####
  
  ### 1.1/ Continous stacks ####
  
  # Load stacks of all submodels projections
  full_stack_cont <- readRDS(file = paste0(initial.wd,"/models/",res,"/",unit,"/proj_Current/proj_Current_",unit,".rds"))
  
  Select_Jaccard_stack_cont <- full_stack_cont[[which(names(full_stack_cont) %in% valid_submodels_names_Jaccard)]]/1000
  Select_TSS_stack_cont <- full_stack_cont[[which(names(full_stack_cont) %in% valid_submodels_names_TSS)]]/1000
  
  # Save the .RData/rds
  save(Select_Jaccard_stack_cont, file = paste0("../../outputs/By_OMU/",unit,"/Select_Jaccard_stack_cont.RData"), version = "2")
  saveRDS(Select_Jaccard_stack_cont, file = paste0("../../outputs/By_OMU/",unit,"/Select_Jaccard_stack_cont.rds"), version = "2")
  save(Select_TSS_stack_cont, file = paste0("../../outputs/By_OMU/",unit,"/Select_TSS_stack_cont.RData"), version = "2")
  saveRDS(Select_TSS_stack_cont, file = paste0("../../outputs/By_OMU/",unit,"/Select_TSS_stack_cont.rds"), version = "2")
  
  # Get the occurrences coordinates
  formated.input.data <- readRDS(file = paste0(initial.wd,"/models/",res,"/",unit,"/formated.input.data.rds"))
  coor_presences <- formated.input.data@coord[which(formated.input.data@data.species == 1), ]
  
  # Plot the maps for Jaccard in a multiple pages pdf
  pdf(file = paste0("../../maps/By_OMU/",unit,"/Select_Jaccard_stack_cont_",unit,".pdf"), height = 8, width = 8)
  for (j in names(Select_Jaccard_stack_cont))
  {
    plot(Select_Jaccard_stack_cont[[j]])
    points(coor_presences, pch = 16)
    points(coor_presences, pch = 16, cex = 0.6, col = "red")
    title(main = j)
  }
  dev.off()
  # Copy in Species folder
  # sp <- as.character(modeled_OMU$Sp_full[i])
  file.copy(from = paste0("../../maps/By_OMU/",unit,"/Select_TSS_stack_cont_",unit,".pdf"), to = paste0("../../maps/By_Species/",sp,"/Select_TSS_stack_cont_",unit,".pdf"), overwrite = T)
  
  # Plot the maps for TSS in a multiple pages pdf
  pdf(file = paste0("../../maps/By_OMU/",unit,"/Select_TSS_stack_cont_",unit,".pdf"), height = 8, width = 8)
  for (j in names(Select_TSS_stack_cont))
  {
    plot(Select_TSS_stack_cont[[j]])
    points(coor_presences, pch = 16, add.plot = T)
    points(coor_presences, pch = 16, cex = 0.6, col = "red", add.plot = T)
    title(main = j)
  }
  dev.off()
  # Copy in Species folder
  # sp <- as.character(modeled_OMU$Sp_full[i])
  file.copy(from = paste0("../../maps/By_OMU/",unit,"/Select_TSS_stack_cont_",unit,".pdf"), to = paste0("../../maps/By_Species/",sp,"/Select_TSS_stack_cont_",unit,".pdf"), overwrite = T)
  
  ### 1.2/ Binary stacks ####
  
  unit_cutoffs <- all.units_cutoffs[unit, , , , ]
  
  # For Jaccard selection
  Select_Jaccard_stack_bin <- Select_Jaccard_stack_cont
  for (j in names(Select_Jaccard_stack_cont)) 
  {
    pa <- strsplit(j, split = "_")[[1]][2]
    cv <- strsplit(j, split = "_")[[1]][3]
    algo <- strsplit(j, split = "_")[[1]][4]
    
    # Extract threshold
    cutoff <- unit_cutoffs[algo, cv, pa, "Jaccard"]
    
    # Generate binary map and store it in the stack
    Select_Jaccard_stack_bin[[j]] <- Select_Jaccard_stack_cont[[j]] > (cutoff/1000)
  }
  
  # For TSS selection
  Select_TSS_stack_bin <- Select_TSS_stack_cont
  for (j in names(Select_TSS_stack_cont)) 
  {
    pa <- strsplit(j, split = "_")[[1]][2]
    cv <- strsplit(j, split = "_")[[1]][3]
    algo <- strsplit(j, split = "_")[[1]][4]
    
    # Extract threshold
    cutoff <- unit_cutoffs[algo, cv, pa, "TSS"]
    
    # Generate binary map and store it in the stack
    Select_TSS_stack_bin[[j]] <- Select_TSS_stack_cont[[j]] > (cutoff/1000)
  }
  
  # Save the .RData/rds
  save(Select_Jaccard_stack_bin, file = paste0("../../outputs/By_OMU/",unit,"/Select_Jaccard_stack_bin.RData"), version = "2")
  saveRDS(Select_Jaccard_stack_bin, file = paste0("../../outputs/By_OMU/",unit,"/Select_Jaccard_stack_bin.rds"), version = "2")
  save(Select_TSS_stack_bin, file = paste0("../../outputs/By_OMU/",unit,"/Select_TSS_stack_bin.RData"), version = "2")
  saveRDS(Select_TSS_stack_bin, file = paste0("../../outputs/By_OMU/",unit,"/Select_TSS_stack_bin.rds"), version = "2")
  
  # Plot the maps for Jaccard in a multiple pages pdf
  pdf(file = paste0("../../maps/By_OMU/",unit,"/Select_Jaccard_stack_bin_",unit,".pdf"), height = 8, width = 8)
  for (j in names(Select_Jaccard_stack_bin))
  {
    plot(Select_Jaccard_stack_bin[[j]])
    points(coor_presences, pch = 16, add.plot = T)
    points(coor_presences, pch = 16, cex = 0.6, col = "red", add.plot = T)
    title(main = j)
  }
  dev.off()
  # Copy in Species folder
  # sp <- as.character(modeled_OMU$Sp_full[i])
  file.copy(from = paste0("../../maps/By_OMU/",unit,"/Select_TSS_stack_bin_",unit,".pdf"), to = paste0("../../maps/By_Species/",sp,"/Select_TSS_stack_bin_",unit,".pdf"), overwrite = T)
  
  # Plot the maps for TSS in a multiple pages pdf
  pdf(file = paste0("../../maps/By_OMU/",unit,"/Select_TSS_stack_bin_",unit,".pdf"), height = 8, width = 8)
  for (j in names(Select_TSS_stack_bin))
  {
    plot(Select_TSS_stack_bin[[j]])
    points(coor_presences, pch = 16, add.plot = T)
    points(coor_presences, pch = 16, cex = 0.6, col = "red", add.plot = T)
    title(main = j)
  }
  dev.off()
  # Copy in Species folder
  # sp <- as.character(modeled_OMU$Sp_full[i])
  file.copy(from = paste0("../../maps/By_OMU/",unit,"/Select_TSS_stack_bin_",unit,".pdf"), to = paste0("../../maps/By_Species/",sp,"/Select_TSS_stack_bin_",unit,".pdf"), overwrite = T)
  
  
  ### 2/ Ensemble
  
  cat(paste0("\n", Sys.time()," ------ Ensemble computation for ", unit, " = Unit N°",i,"/",nrow(modeled_OMU)," ------\n"))
  
  ### 2.1/ Continuous Ensemble ####
  
  Ensemble_Jaccard_median <- calc(x = Select_Jaccard_stack_cont, fun = median, na.rm = T)
  Ensemble_TSS_median <- calc(x = Select_TSS_stack_cont, fun = median, na.rm = T)
  
  save(Ensemble_Jaccard_median, file = paste0("../../outputs/By_OMU/",unit,"/Ensemble_Jaccard_median.RData"), version = "2")
  saveRDS(Ensemble_Jaccard_median, file = paste0("../../outputs/By_OMU/",unit,"/Ensemble_Jaccard_median.rds"), version = "2")
  save(Ensemble_TSS_median, file = paste0("../../outputs/By_OMU/",unit,"/Ensemble_TSS_median.RData"), version = "2")
  saveRDS(Ensemble_TSS_median, file = paste0("../../outputs/By_OMU/",unit,"/Ensemble_TSS_median.rds"), version = "2")
  
  ### 2.2/ Incertitude on Ensemble (sd) ####
  
  Ensemble_Jaccard_sd <- calc(x = Select_Jaccard_stack_cont, fun = sd, na.rm = T)
  Ensemble_TSS_sd <- calc(x = Select_TSS_stack_cont, fun = sd, na.rm = T)
  
  save(Ensemble_Jaccard_sd, file = paste0("../../outputs/By_OMU/",unit,"/Ensemble_Jaccard_sd.RData"), version = "2")
  saveRDS(Ensemble_Jaccard_sd, file = paste0("../../outputs/By_OMU/",unit,"/Ensemble_Jaccard_sd.rds"), version = "2")
  save(Ensemble_TSS_sd, file = paste0("../../outputs/By_OMU/",unit,"/Ensemble_TSS_sd.RData"), version = "2")
  saveRDS(Ensemble_TSS_sd, file = paste0("../../outputs/By_OMU/",unit,"/Ensemble_TSS_sd.rds"), version = "2")
  
  ### 2.3/ Committee Averaging for binary map ####
  
  # Jaccard
  Ensemble_Jaccard_CA <- sum(Select_Jaccard_stack_bin)                                        # Sum all binary maps to get the nb of votes
  Ensemble_Jaccard_CA[Ensemble_Jaccard_CA < (nlayers(Select_Jaccard_stack_bin) / 2)] <- 0     # If less then half votes for presence, ensemble displays an absence
  Ensemble_Jaccard_CA[Ensemble_Jaccard_CA >= (nlayers(Select_Jaccard_stack_bin) / 2)] <- 1    # If half or more votes for presence, ensemble displays a presence
  # plot(Ensemble_Jaccard_CA)
  
  # TSS
  Ensemble_TSS_CA <- sum(Select_TSS_stack_bin)                                    # Sum all binary maps to get the nb of votes
  Ensemble_TSS_CA[Ensemble_TSS_CA < (nlayers(Select_TSS_stack_bin) / 2)] <- 0     # If less then half votes for presence, ensemble displays an absence
  Ensemble_TSS_CA[Ensemble_TSS_CA >= (nlayers(Select_TSS_stack_bin) / 2)] <- 1    # If half or more votes for presence, ensemble displays a presence
  # plot(Ensemble_TSS_CA)
  
  save(Ensemble_Jaccard_CA, file = paste0("../../outputs/By_OMU/",unit,"/Ensemble_Jaccard_CA.RData"), version = "2")
  saveRDS(Ensemble_Jaccard_CA, file = paste0("../../outputs/By_OMU/",unit,"/Ensemble_Jaccard_CA.rds"), version = "2")
  save(Ensemble_TSS_CA, file = paste0("../../outputs/By_OMU/",unit,"/Ensemble_TSS_CA.RData"), version = "2")
  saveRDS(Ensemble_TSS_CA, file = paste0("../../outputs/By_OMU/",unit,"/Ensemble_TSS_CA.rds"), version = "2")
  
  ### 2.4/ Binary map based on evaluation of all observations fitted in all selected models ####
  
  # Note this is biased because evaluation is done on observations used to calibrate some of the models, and nb of PsA depend on model type (more PsA for "restricted")
  
  # Retrieve input data
  obs.data <- formated.input.data@data.species
  obs.data[is.na(obs.data)] <- 0   # Change pseudo-absences (NA) in real absences (0) to be able to compute confusion matrix
  
  # Retrieve coordinates
  coor_all <- formated.input.data@coord
  
  ## Jaccard
  
  # Extract predictions from the EM
  pred.data <- raster::extract(Ensemble_Jaccard_median*1000, coor_all)
  
  jaccard.test <- NULL   # Object to store the cutoff research                                
  for(cutoff in seq(0, 1000, by = 1))                    
  {
    pred.pa <- pred.data              # Predictions of the model for the test set
    pred.pa[pred.pa < cutoff] <- 0    # We assigne a 0 (Absence) for all predicts under the cutoff 
    pred.pa[pred.pa >= cutoff] <- 1   # We assign a 1 (presence) for all predicts equal or above the cutoff
    
    # We compare predicts and inputs = confusion matrix
    TP <- sum(obs.data == 1 & pred.pa == 1, na.rm = T)   # True positives
    FN <- sum(obs.data == 1 & pred.pa == 0, na.rm = T)   # False negatives
    FP <- sum(obs.data == 0 & pred.pa == 1, na.rm = T)   # False positives
    
    # Compute Jaccard index
    jaccard <- TP / (TP + FP + FN)
    
    # Store all infos regarding the cut.off test for TSS in a df
    jaccard.test <- rbind.data.frame(jaccard.test,
                                     data.frame(cutoff = cutoff,
                                                TP = TP,
                                                FN = FN,
                                                FP = FP,
                                                jaccard = jaccard))
  }
  # Store the test track
  jaccard.test.list_EM[[unit]] <- jaccard.test
  
  # Extract the cutoff as the value that maximize the evaluation
  EM.cutoff_Jaccard <- round(mean(jaccard.test$cutoff[which(jaccard.test$jaccard == max(jaccard.test$jaccard))])) # Take the mean value because several cutoffs can maximize the evaluation metric
  EM.metric_Jaccard <- round(jaccard.test$jaccard[round(which(jaccard.test$cutoff == EM.cutoff_Jaccard))],3)
  # Store it in list.models
  list.models$EM.cutoff_Jaccard[list.models$Tag.model == unit] <- EM.cutoff_Jaccard
  list.models$EM.metric_Jaccard[list.models$Tag.model == unit] <- EM.metric_Jaccard
  
  # Binariazation of the continous EM
  Ensemble_Jaccard_bin <- BinaryTransformation(Ensemble_Jaccard_median, EM.cutoff_Jaccard/1000)
  # plot(Ensemble_Jaccard_bin)
  
  ## TSS
  
  # Extract predictions from the EM
  pred.data <- raster::extract(Ensemble_TSS_median*1000, coor_all)
  
  TSS.test <- NULL   # Object to store the cutoff research                                
  for(cutoff in seq(0, 1000, by = 1))                    
  {
    pred.pa <- pred.data              # Predictions of the model for the test set
    pred.pa[pred.pa < cutoff] <- 0    # We assigne a 0 (Absence) for all predicts under the cutoff 
    pred.pa[pred.pa >= cutoff] <- 1   # We assign a 1 (presence) for all predicts equal or above the cutoff
    
    # We compare predicts and inputs = confusion matrix
    TP <- sum(obs.data == 1 & pred.pa == 1, na.rm = T)   # True positives
    FN <- sum(obs.data == 1 & pred.pa == 0, na.rm = T)   # False negatives
    FP <- sum(obs.data == 0 & pred.pa == 1, na.rm = T)   # False positives
    TN <- sum(obs.data == 0 & pred.pa == 0, na.rm = T)   # True negatives
    
    sensitivity <- round(TP / (TP + FN), 3) # Ability to detect presence among presences points
    specificity <- round(TN / (TN + FP), 3)  # Ability to detect absences among absences points
    
    TSS <- sensitivity + specificity - 1  # Compute TSS
    
    # Store all infos regarding the cut.off test for TSS in a df
    TSS.test <- rbind.data.frame(TSS.test,                  
                                 data.frame(cutoff = cutoff,
                                            sensitivity = sensitivity,
                                            specificity = specificity,
                                            TSS = TSS))
  }
  # Store the test track
  TSS.test.list_EM[[unit]] <- TSS.test
  
  # Extract the cutoff as the value that maximize the evaluation
  EM.cutoff_TSS <- round(mean(TSS.test$cutoff[which(TSS.test$TSS == max(TSS.test$TSS))])) # Take the mean value because several cutoffs can maximize the evaluation metric
  EM.metric_TSS <- round(TSS.test$TSS[round(which(TSS.test$cutoff == EM.cutoff_TSS))],3)
  # Store it in list.models
  list.models$EM.cutoff_TSS[list.models$Tag.model == unit] <- EM.cutoff_TSS
  list.models$EM.metric_TSS[list.models$Tag.model == unit] <- EM.metric_TSS
  
  # Binariazation of the continous EM
  Ensemble_TSS_bin <- BinaryTransformation(Ensemble_TSS_median, EM.cutoff_TSS/1000)
  # plot(Ensemble_Jaccard_bin)
  
  ## Save binary EMs
  save(Ensemble_Jaccard_bin, file = paste0("../../outputs/By_OMU/",unit,"/Ensemble_Jaccard_bin.RData"), version = "2")
  saveRDS(Ensemble_Jaccard_bin, file = paste0("../../outputs/By_OMU/",unit,"/Ensemble_Jaccard_bin.rds"), version = "2")
  save(Ensemble_TSS_bin, file = paste0("../../outputs/By_OMU/",unit,"/Ensemble_TSS_bin.RData"), version = "2")
  saveRDS(Ensemble_TSS_bin, file = paste0("../../outputs/By_OMU/",unit,"/Ensemble_TSS_bin.rds"), version = "2")
  
  ### 2.5/ Plot all EMs on one page ####
  
  # Jaccard
  all_EM_Jaccard <- stack(Ensemble_Jaccard_median, Ensemble_Jaccard_sd, Ensemble_Jaccard_bin, Ensemble_Jaccard_CA)
  names(all_EM_Jaccard) <- c("Continuous", "Incertitude (sd)", "Binary", "Committee Averaging")
  
  pdf(file = paste0("../../maps/By_OMU/",unit,"/Ensemble_maps_Jaccard_",unit,".pdf"), height = 6, width = 7)
  internal_margins <- par()$mar
  external_margins <- par()$oma
  par(mar = c(5.1,6.1,9.1,3.1))
  par(oma = c(1,2,4,3))
  plot(all_EM_Jaccard)
  title(main = unit, outer = T)
  par(mar = internal_margins, oma = external_margins)
  dev.off()
  # Copy in species folder
  # sp <- as.character(modeled_OMU$Sp_full[i])
  file.copy(from = paste0("../../maps/By_OMU/",unit,"/Ensemble_maps_Jaccard_",unit,".pdf"), to = paste0("../../maps/By_Species/",sp,"/Ensemble_maps_Jaccard_",unit,".pdf"), overwrite = T)
  
  save(all_EM_Jaccard, file = paste0("../../outputs/By_OMU/",unit,"/all_EM_Jaccard.RData"), version = "2")
  saveRDS(all_EM_Jaccard, file = paste0("../../outputs/By_OMU/",unit,"/all_EM_Jaccard.rds"), version = "2")
  
  # TSS
  all_EM_TSS <- stack(Ensemble_TSS_median, Ensemble_TSS_sd, Ensemble_TSS_bin, Ensemble_TSS_CA)
  names(all_EM_TSS) <- c("Continuous", "Incertitude (sd)", "Binary", "Committee Averaging")
  
  pdf(file = paste0("../../maps/By_OMU/",unit,"/Ensemble_maps_TSS_",unit,".pdf"), height = 6, width = 7)
  internal_margins <- par()$mar
  external_margins <- par()$oma
  par(mar = c(5.1,6.1,9.1,3.1))
  par(oma = c(1,2,4,3))
  plot(all_EM_TSS)
  title(main = unit, outer = T)
  par(mar = internal_margins, oma = external_margins)
  dev.off()
  # Copy in species folder
  # sp <- as.character(modeled_OMU$Sp_full[i])
  file.copy(from = paste0("../../maps/By_OMU/",unit,"/Ensemble_maps_TSS_",unit,".pdf"), to = paste0("../../maps/By_Species/",sp,"/Ensemble_maps_TSS_",unit,".pdf"), overwrite = T)
  
  save(all_EM_TSS, file = paste0("../../outputs/By_OMU/",unit,"/all_EM_TSS.RData"), version = "2")
  saveRDS(all_EM_TSS, file = paste0("../../outputs/By_OMU/",unit,"/all_EM_TSS.rds"), version = "2")
  
  ### 3/ Update boxplots for evaluation with only selected submodels and Eval metric value ####
  
  cat(paste0("\n", Sys.time()," ------ Evaluation boxplots for ", unit, " = Unit N°",i,"/",nrow(modeled_OMU)," ------\n"))
  
  unit_model_type <- modeled_OMU$initial_model_type[i]
  
  ## Jaccard
  
  gg_jaccard_unit <- gg_jaccard[gg_jaccard$unit == unit, ]
  
  # Update the df
  gg_jaccard_unit <- gg_jaccard_unit %>% 
    mutate(model_name = paste0(unit, "_", PA_sets, "_", CV_runs, "_", Algo)) %>%  # Add a model_name column
    filter(model_name %in% valid_submodels_names_Jaccard) %>%  # Filter only selected submodels
    add_row(unit = unit, Algo = "EM", value = EM.metric_Jaccard) # Add the EM value
  
  # Order the factor levels
  gg_jaccard_unit$Algo <- factor(gg_jaccard_unit$Algo, levels = c("RF", "GBM", "ANN", "EM"))
  
  # Plot
  pdf(file = paste0("../../controls/evaluations/By_OMU/boxplots/Jaccard/Postselect_Eval_Jaccard_boxplot_",unit,".pdf"), height = 6, width = 10)
  g <- ggplot(gg_jaccard_unit, aes(x = Algo, y = value, fill = Algo)) +
    geom_boxplot(show.legend = F) +
    ylim(c(min(gg_jaccard_unit$value)*0.9, 1)) +
    geom_hline(yintercept = c(0.6, 0.95)[c("complete", "restricted") == unit_model_type], col = "red", lwd = 1, lty = 2) +
    labs(title = paste0("Evaluation of selected models for ",unit," with Jaccard indices"),
         subtitle = paste0("Model type = ",unit_model_type),
         y = "Jaccard indices",
         x = "Algorithm or Ensemble")
  print(g)
  dev.off()
  
  ## TSS
  
  gg_TSS_unit <- gg_TSS[gg_TSS$unit == unit, ]
  
  # Update the df
  gg_TSS_unit <- gg_TSS_unit %>% 
    mutate(model_name = paste0(unit, "_", PA_sets, "_", CV_runs, "_", Algo)) %>%  # Add a model_name column
    filter(model_name %in% valid_submodels_names_TSS) %>%  # Filter only selected submodels
    add_row(unit = unit, Algo = "EM", value = EM.metric_TSS) # Add the EM value
  
  # Order the factor levels
  gg_TSS_unit$Algo <- factor(gg_TSS_unit$Algo, levels = c("RF", "GBM", "ANN", "EM"))
  
  # Plot
  pdf(file = paste0("../../controls/evaluations/By_OMU/boxplots/TSS/Postselect_Eval_TSS_boxplot_",unit,".pdf"), height = 6, width = 10)
  g <- ggplot(gg_TSS_unit, aes(x = Algo, y = value, fill = Algo)) +
    geom_boxplot(show.legend = F) +
    ylim(c(min(gg_TSS_unit$value)*0.9, 1)) +
    geom_hline(yintercept = c(0.6, 0.95)[c("complete", "restricted") == unit_model_type], col = "red", lwd = 1, lty = 2) +
    labs(title = paste0("Evaluation of selected models for ",unit," with TSS indices"),
         subtitle = paste0("Model type = ",unit_model_type),
         y = "TSS indices",
         x = "Algorithm or Ensemble")
  print(g)
  dev.off()
 
  cat(paste0("\n", Sys.time()," ------ Done for ", unit, " = Unit N°",i,"/",nrow(modeled_OMU)," ------\n"))
  
}

# save(Jaccard.test.list_EM, file = paste0("../../controls/evaluations/Jaccard.test.list_EM.RData"), version = "2")
# saveRDS(Jaccard.test.list_EM, file = paste0("../../controls/evaluations/Jaccard.test.list_EM.rds"), version = "2")
# save(TSS.test.list_EM, file = paste0("../../controls/evaluations/TSS.test.list_EM.RData"), version = "2")
# saveRDS(TSS.test.list_EM, file = paste0("../../controls/evaluations/TSS.test.list_EM.rds"), version = "2")

save(Jaccard.test.list_EM, file = paste0("../../controls/evaluations/Jaccard.test.list_EM_",first(index_model_to_compute),"_",last(index_model_to_compute),".RData"), version = "2")
saveRDS(Jaccard.test.list_EM, file = paste0("../../controls/evaluations/Jaccard.test.list_EM_",first(index_model_to_compute),"_",last(index_model_to_compute),".rds"), version = "2")
save(TSS.test.list_EM, file = paste0("../../controls/evaluations/TSS.test.list_EM_",first(index_model_to_compute),"_",last(index_model_to_compute),".RData"), version = "2")
saveRDS(TSS.test.list_EM, file = paste0("../../controls/evaluations/TSS.test.list_EM_",first(index_model_to_compute),"_",last(index_model_to_compute),".rds"), version = "2")


# setwd("D:/Mael/R_projects/ithomiini_current/")
setwd(initial.wd) # Back to initial wd

# Save Summary table for OMU models
# save(list.models, file = paste0("./input_data/list.models.RData"), version = "2")
# saveRDS(list.models, file = paste0("./input_data/list.models.rds"), version = "2")



### 4/ Update of global boxplots for evaluation metrics ####

library(tidyverse)
library(ggplot2)


## Jaccard

gg_jaccard <- readRDS(file = "./controls/evaluations/Full_df_Jaccard.rds")

# Update the df  # Add EM values from the list.models
gg_jaccard <- gg_jaccard %>% 
  mutate(model_name = paste0(unit, "_", PA_sets, "_", CV_runs, "_", Algo)) %>%  # Add a model_name column
  filter(model_name %in% valid_submodels_names_Jaccard) %>%  # Filter only selected submodels
  add_row(unit = list.models$Tag.model, Algo = "EM", CV_runs = "All", PA_sets = "All", value = list.models$EM.metric_Jaccard, model_type = list.models$initial_model_type, model_name = paste0("EM_",unit)) %>%  # Add the EM value
  na.omit()
  
# Order the factor levels
gg_jaccard$Algo <- factor(gg_jaccard$Algo, levels = c("RF", "GBM", "ANN", "EM"))

pdf(file = paste0("./controls/evaluations/All_units_Postselect_Eval_Jaccard_boxplot.pdf"), height = 6, width = 10) 
g <- ggplot(gg_jaccard, aes(x = Algo, y = value, fill = model_type)) +
       geom_boxplot(show.legend = T) +
       # ylim(c(0, 1)) +
       geom_hline(yintercept = 0.6, col = "red", lwd = 1, lty = 2) +
       geom_hline(yintercept = 0.95, col = "dodgerblue", lwd = 1, lty = 2) +
       labs(title = "Global evaluation of all models after selection with Jaccard indices",
            y = "Jaccard indices",
            x = "Algorithm")
print(g)
dev.off()

# Update the df  # Add EM values from the list.models
gg_jaccard <- gg_jaccard %>% 
  mutate(model_name = paste0(unit, "_", PA_sets, "_", CV_runs, "_", Algo)) %>%  # Add a model_name column
  filter(model_name %in% valid_submodels_names_Jaccard) %>%  # Filter only selected submodels
  add_row(unit = list.models$Tag.model, Algo = "EM", CV_runs = "All", PA_sets = "All", value = list.models$EM.metric_Jaccard, model_type = list.models$initial_model_type, model_name = paste0("EM_",unit)) %>%  # Add the EM value
  na.omit()

# Order the factor levels
gg_jaccard$Algo <- factor(gg_jaccard$Algo, levels = c("RF", "GBM", "ANN", "EM"))

pdf(file = paste0("./controls/evaluations/All_units_Postselect_Eval_Jaccard_boxplot.pdf"), height = 6, width = 10) 
g <- ggplot(gg_jaccard, aes(x = Algo, y = value, fill = model_type)) +
  geom_boxplot(show.legend = T) +
  # ylim(c(0, 1)) +
  geom_hline(yintercept = 0.6, col = "red", lwd = 1, lty = 2) +
  geom_hline(yintercept = 0.95, col = "dodgerblue", lwd = 1, lty = 2) +
  labs(title = "Global evaluation of all models after selection with Jaccard indices",
       y = "Jaccard indices",
       x = "Algorithm")
print(g)
dev.off()



## TSS

gg_TSS <- readRDS(file = "./controls/evaluations/Full_df_TSS.rds")

# Update the df  # Add EM values from the list.models
gg_TSS <- gg_TSS %>% 
  mutate(model_name = paste0(unit, "_", PA_sets, "_", CV_runs, "_", Algo)) %>%  # Add a model_name column
  filter(model_name %in% valid_submodels_names_TSS) %>%  # Filter only selected submodels
  add_row(unit = list.models$Tag.model, Algo = "EM", CV_runs = "All", PA_sets = "All", value = list.models$EM.metric_TSS, model_type = list.models$initial_model_type, model_name = paste0("EM_",unit)) %>%  # Add the EM value
  na.omit()

# Order the factor levels
gg_TSS$Algo <- factor(gg_TSS$Algo, levels = c("RF", "GBM", "ANN", "EM"))

pdf(file = paste0("./controls/evaluations/All_units_Postselect_Eval_TSS_boxplot.pdf"), height = 6, width = 10) 
g <- ggplot(gg_TSS, aes(x = Algo, y = value, fill = model_type)) +
  geom_boxplot(show.legend = T) +
  # ylim(c(0, 1)) +
  geom_hline(yintercept = 0.6, col = "red", lwd = 1, lty = 2) +
  geom_hline(yintercept = 0.95, col = "dodgerblue", lwd = 1, lty = 2) +
  labs(title = "Global evaluation of all models after selection with TSS indices",
       y = "TSS indices",
       x = "Algorithm")
print(g)
dev.off()

# Update the df  # Add EM values from the list.models
gg_TSS <- gg_TSS %>% 
  mutate(model_name = paste0(unit, "_", PA_sets, "_", CV_runs, "_", Algo)) %>%  # Add a model_name column
  filter(model_name %in% valid_submodels_names_TSS) %>%  # Filter only selected submodels
  add_row(unit = list.models$Tag.model, Algo = "EM", CV_runs = "All", PA_sets = "All", value = list.models$EM.metric_TSS, model_type = list.models$initial_model_type, model_name = paste0("EM_",unit)) %>%  # Add the EM value
  na.omit()

# Order the factor levels
gg_TSS$Algo <- factor(gg_TSS$Algo, levels = c("RF", "GBM", "ANN", "EM"))

pdf(file = paste0("./controls/evaluations/All_units_Postselect_Eval_TSS_boxplot.pdf"), height = 6, width = 10) 
g <- ggplot(gg_TSS, aes(x = Algo, y = value, fill = model_type)) +
  geom_boxplot(show.legend = T) +
  # ylim(c(0, 1)) +
  geom_hline(yintercept = 0.6, col = "red", lwd = 1, lty = 2) +
  geom_hline(yintercept = 0.95, col = "dodgerblue", lwd = 1, lty = 2) +
  labs(title = "Global evaluation of all models after selection with TSS indices",
       y = "TSS indices",
       x = "Algorithm")
print(g)
dev.off()


### 5/ Build a global response plot for all OMU using only submodels retained for ensemble ####

# Generate df to store data for global response plot
resp.df_all_units_Jaccard <- resp.df_all_units_TSS <- tibble()

# Load Summary table for OMU models to extract list of unit
load(file = paste0("./input_data/list.models.RData"))

# Extract only OMU with a valid EM
list.models_valid <- list.models[list.models$initial_model_type %in% c("complete", "restricted"),]

# Load summary table for sub_models to retrieve the list of valid submodel to include in the ensemble
list.submodels <- readRDS(file = "./input_data/list.submodels.rds")

valid_submodels_names_Jaccard <- list.submodels$model_name[list.submodels$Ensemble_OK_Jaccard] # Jaccard
valid_submodels_names_TSS <- list.submodels$model_name[list.submodels$Ensemble_OK_TSS] # TSS

# To choose the eval metric for selection
evaluation_metric_list <- c("Jaccard", "TSS")

# i <- 1

# Loop per unit
for (i in seq_along(list.models_valid$Tag.model))
{
  unit <- list.models_valid$Tag.model[i]
  
  # Get df storing all info to build response plot for one OMU
  resp.df <- readRDS(file = paste0("./controls/response_plots/By_OMU/",unit,"/resp.df_",unit,".rds"))
  
  # For each evaluation metric
  for(j in seq_along(evaluation_metric_list))
  {
    evaluation_metric <- evaluation_metric_list[j]
    
    # Filter the resp.df
    if(evaluation_metric == "Jaccard")
    {
      resp.df_filtered <- resp.df[which(resp.df$Model %in% valid_submodels_names_Jaccard), ] # Jaccard
    } else {
      resp.df_filtered <- resp.df[which(resp.df$Model %in% valid_submodels_names_TSS), ] # TSS
    }
    
    # Extract median value to built df for all OMU EM
    resp.df_unit_EM <- resp.df_filtered %>% 
      group_by(Variable, Var.value) %>%  # Regroup rows by predictor (Variable) and predictor value (Var.value)
      summarize(unit = unit,
                median_pred = median(Response, na.rm = TRUE),  # Compute median prediction = EM
                Range_min = first(Range_min), # Keep range info
                Range_max = first(Range_max)) %>% # Keep range info
      ungroup()
    
    
    # Store the data into the global df for Jaccard and TSS
    if(evaluation_metric == "Jaccard")
    {
      resp.df_all_units_Jaccard <- rbind(resp.df_all_units_Jaccard, resp.df_unit_EM)
      saveRDS(resp.df_all_units_Jaccard, file = paste0("./controls/response_plots/resp.df_all_units_",evaluation_metric,".rds"))
    } else {
      resp.df_all_units_TSS <- rbind(resp.df_all_units_TSS, resp.df_unit_EM)
      saveRDS(resp.df_all_units_TSS, file = paste0("./controls/response_plots/resp.df_all_units_",evaluation_metric,".rds"))
      
    }
        
    # Plot only the EM response curves for this unit
    
    # p_unit <- ggplot(resp.df_unit_EM, aes(x = Var.value, y = median_pred)) + 
    #   geom_line(alpha = 1) + # Plot lines for each model
    #   facet_wrap(~Variable, scales = "free_x") +  # make a plot per variable
    #   theme_bw() + 
    #   ylim(0, 1) +
    #   xlab("Variable value") +
    #   labs(title = paste0("Post-selection response plots for ", unit, "\n", evaluation_metric, " threshold")) +
    #   geom_vline(aes(xintercept = Range_min), lwd = 1, col = "red") +
    #   geom_vline(aes(xintercept = Range_max), lwd = 1, col = "red")
    # 
    # print(p_unit)
    
    # 5.1/ Update response plots for this unit/OMU. ####
    
    # Post-selection (Models are discarded on the basis of evaluation or weird response curve) 
    
    p_unit <- ggplot(resp.df_filtered, aes(x = Var.value, y = Response)) + 
      geom_line(alpha = 0.2, aes(group = Model)) + # Plot lines for each model
      stat_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) + # Generate GAM smoothed curve for each mini-plot
      facet_wrap(~Variable, scales = "free_x") +  # make a plot per variable
      theme_bw() + 
      ylim(0, 1) +
      xlab("Variable value") +
      labs(title = paste0("Post-selection response plots for ", unit, "\n", evaluation_metric, " threshold")) +
      geom_vline(aes(xintercept = Range_min), lwd = 1, col = "red") +
      geom_vline(aes(xintercept = Range_max), lwd = 1, col = "red")
    
    # Save plot for this OMU
    pdf(file = paste0("./controls/response_plots/By_OMU/",unit,"/Post_select_response_plots_",evaluation_metric, "_", unit,".pdf"), height = 6, width = 10)
    print(p_unit)
    dev.off()
    
  }
  
  # Check advancement
  if(i %% 10 == 0)
  {
    cat(paste0("\n", Sys.time()," ------ ", unit, " = Unit N°",i,"/",nrow(list.models_valid)," ------\n")) 
  }
  
}

# 5.2/ Global response plot for all OMUs ####

evaluation_metric_list <- c("Jaccard", "TSS")

for(j in seq_along(evaluation_metric_list))
{
  evaluation_metric <- evaluation_metric_list[j]
  
  # Load data
  resp.df_all_units <- readRDS(file = paste0("./controls/response_plots/resp.df_all_units_",evaluation_metric,".rds"))
  
  # Update range min and max for all OMUs
  resp.df_all_units_ranges <- resp.df_all_units %>% 
    group_by(Variable) %>% 
    summarize(Range_min = min(Range_min), # Keep range info
              Range_max = max(Range_max)) %>% # Keep range info
    ungroup() %>% 
    left_join(x = resp.df_all_units[, !(names(resp.df_all_units) %in% c("Range_min", "Range_max"))], y = ., by = "Variable")
  
  # Generate plot
  p_all_units <- ggplot(resp.df_all_units_ranges, aes(x = Var.value, y = median_pred)) + 
    geom_line(alpha = 0.2, aes(group = unit)) + # Plot lines for each model
    stat_smooth(method = 'gam', formula = y ~ s(x, bs = "cs")) + # Generate GAM smoothed curve for each mini-plot
    facet_wrap(~Variable, scales = "free_x") +  # make a plot per variable
    theme_bw() + 
    ylim(0, 1) +
    xlab("Variable value") +
    labs(title = paste0("Post-selection response plots for all modeled OMUs\n", evaluation_metric, " threshold")) +
    geom_vline(aes(xintercept = Range_min), lwd = 1, col = "red") +
    geom_vline(aes(xintercept = Range_max), lwd = 1, col = "red")
  
  # Save response plot for all OMU
  pdf(file = paste0("./controls/response_plots/All_OMUs_Post_select_response_plots_",evaluation_metric,".pdf"), height = 6, width = 10)
  print(p_all_units)
  dev.off()
}











