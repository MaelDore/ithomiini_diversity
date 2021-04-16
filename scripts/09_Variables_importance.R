
##### Script 9: Variables Importance #####

# Extract variable importance from Biomod2 results
# Plot boxplot per variable for each OMU and globally

#####

# Inputs:
   # Results of Biomod2 from Script 05.1

# Outputs:
   # Boxplot of variable importance per OMU, and globally
   # List of array storing variable importance for each submodel


#####


#### 1/ Plots per OMU/unit ####

# Effacer l'environnement
rm(list = ls())

# setwd("D:/Mael/R_projects/ithomiini_current/")

# Load libraries

library(biomod2)
library(reshape2)
library(ggplot2)

# Select Env data resolution 
# res <- "5"
res <- "15"

# Load Summary table for OMU models
load(file = paste0("./input_data/list.models.RData"))
modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]

# Load summary table for sub_models
list.submodels <- readRDS(file = "./input_data/list.submodels.rds")

# Choose the eval metric for selection
valid_submodels_names <- list.submodels$model_name[list.submodels$Ensemble_OK_Jaccard] # Jaccard
# valid_submodels_names <- list.submodels$model_name[list.submodels$Ensemble_OK_TSS] # TSS

# Set up directory
initial.wd <- getwd()
setwd(paste0(initial.wd,"/models/",res,"/"))

# Lists to store df of all units
var.imp_full_list <- var.imp_ensemble_list <- NULL

### Loop for all modeled OMU/unit

for (i in 1:nrow(modeled_OMU)) 
{
  # i <- 1
  
  unit <- as.character(modeled_OMU$Tag.model[i])
  sp <- as.character(modeled_OMU$Sp_full[i])
  
  cat(paste0("\n", Sys.time()," ------ Starts for ", unit, " = Unit N°",i,"/",nrow(modeled_OMU)," ------\n"))
  
  # Create directory to store variable importance by OMU
  if(!dir.exists(paste0("../../controls/variable_importance/By_OMU/",unit))) { # Test if the folder exists already or not
    dir.create(paste0("../../controls/variable_importance/By_OMU/",unit), recursive = T) # Create folder if absent
  }
  
  # Create directory to store variable importance by Species
  if(!dir.exists(paste0("../../controls/variable_importance/By_species/",sp))) { # Test if the folder exists already or not
    dir.create(paste0("../../controls/variable_importance/By_species/",sp), recursive = T) # Create folder if absent
  }
  
  # Retrieve calibration lines to know which data point to use as test set for each sub_model
  model.runs <- readRDS(file = paste0(unit,"/model.runs.rds"))
  
  # Extract variable importance
  var.imp_array <- get_variables_importance(model.runs)
  # Rename variable
  dimnames(var.imp_array)[[1]] <- c("Tmean", "Tvar", "Htot", "Hvar", "Elevation", "Forests")
  
  # Retrieve variable importance under df format
  Var.imp_unit <- reshape2::melt(var.imp_array) 
  colnames(Var.imp_unit) <- c("Variable", "Algo", "CV_runs", "PA_runs", "Variable.importance")
  Var.imp_unit$model_name <- paste0(unit,"_",Var.imp_unit$PA_runs,"_",Var.imp_unit$CV_runs,"_",Var.imp_unit$Algo)
  Var.imp_unit <- na.omit(Var.imp_unit)
  
  # Extract only submodel that are considered valid for ensemble
  Var.imp_unit_selected <- Var.imp_unit[Var.imp_unit$model_name %in% valid_submodels_names,]
  
  # Order variable factor level as increasing median importance
  Var.imp_unit$Variable <- reorder(Var.imp_unit$Variable, Var.imp_unit$Variable.importance, median, na.rm = TRUE)
  Var.imp_unit_selected$Variable <- reorder(Var.imp_unit_selected$Variable, Var.imp_unit_selected$Variable.importance, median, na.rm = TRUE)
  # Decreasing order
  Var.imp_unit$Variable <- factor(Var.imp_unit$Variable, levels = rev(levels(Var.imp_unit$Variable)))
  Var.imp_unit_selected$Variable <- factor(Var.imp_unit_selected$Variable, levels = rev(levels(Var.imp_unit_selected$Variable)))
  
  # Save the complete df
  save(Var.imp_unit, file = paste0(initial.wd, "/controls/variable_importance/By_OMU/",unit,"/Var.imp_full_",unit,".RData"))
  saveRDS(Var.imp_unit, file = paste0(initial.wd, "/controls/variable_importance/By_OMU/",unit,"/Var.imp_full_",unit,".rds"))
  var.imp_full_list[[unit]] <- Var.imp_unit
  
  # Save the df for selected models for Ensemble
  save(Var.imp_unit_selected, file = paste0(initial.wd, "/controls/variable_importance/By_OMU/",unit,"/Var.imp_selected_",unit,".RData"))
  saveRDS(Var.imp_unit_selected, file = paste0(initial.wd, "/controls/variable_importance/By_OMU/",unit,"/Var.imp_selected_",unit,".rds"))
  var.imp_ensemble_list[[unit]] <- Var.imp_unit_selected
  
  # Plot for All models
  p_unit <- ggplot(Var.imp_unit, aes(x = Variable, y = Variable.importance, fill = Variable)) +
    geom_boxplot() + 
    labs(title = paste0("All models for ",unit),
         y = "Variable importance")
    # theme(axis.text.x = element_text(angle = -90))
  
  # Save plots for all models
  pdf(file = paste0("../../controls/variable_importance/By_OMU/",unit,"/All_models_var.imp_",unit,".pdf"), height = 7, width = 8)
  print(p_unit)
  dev.off()
  # Copy in species folder
  file.copy(from = paste0("../../controls/variable_importance/By_OMU/",unit,"/All_models_var.imp_",unit,".pdf"), to = paste0("../../controls/variable_importance/By_species/",sp,"/All_models_var.imp_",unit,".pdf"), overwrite = T)
  
  # Plot for models selected for Ensemble
  p_unit_selected <- ggplot(Var.imp_unit_selected, aes(x = Variable, y = Variable.importance, fill = Variable)) +
    geom_boxplot() + 
    labs(title = paste0("Ensemble models for ",unit),
         y = "Variable importance")
  # theme(axis.text.x = element_text(angle = -90))
  
  # Save plots for all models
  pdf( file = paste0("../../controls/variable_importance/By_OMU/",unit,"/Ensemble_models_var.imp_",unit,".pdf"), height = 7, width = 8)
  print(p_unit_selected)
  dev.off()
  # Copy in species folder
  file.copy(from = paste0("../../controls/variable_importance/By_OMU/",unit,"/All_models_var.imp_",unit,".pdf"), to = paste0("../../controls/variable_importance/By_species/",sp,"/All_models_var.imp_",unit,".pdf"), overwrite = T)
  
  cat(paste0("\n", Sys.time()," ------ Done for ", unit, " = Unit N°",i,"/",nrow(modeled_OMU)," ------\n"))
  
}

save(var.imp_ensemble_list, file = paste0("../../controls/variable_importance/var.imp_ensemble_list.RData"), version = "2")
saveRDS(var.imp_ensemble_list, file = paste0("../../controls/variable_importance/var.imp_ensemble_list.RData"), version = "2")

save(var.imp_ensemble_list, file = paste0("../../controls/variable_importance/var.imp_ensemble_list.RData"), version = "2")
saveRDS(var.imp_full_list, file = paste0("../../controls/variable_importance/var.imp_full_list.RData"), version = "2")

# setwd("D:/Mael/R_projects/ithomiini_current/")
setwd(initial.wd) # Back to initial wd


#### 2/ Plot for all units confounded ####

var.imp_ensemble_list <- readRDS(file = paste0("./controls/variable_importance/var.imp_ensemble_list.RData"))
var.imp_full_list <- readRDS(file = paste0("./controls/variable_importance/var.imp_full_list.RData"))

### Generate a df from the lists
library(tidyverse)
library(stringr)

var.imp_full_df <- reduce(.x = var.imp_full_list, .f = rbind)
var.imp_ensemble_df <- reduce(.x = var.imp_ensemble_list, .f = rbind)

# Add the unit column
var.imp_full_df <- var.imp_full_df %>% 
  mutate(unit = stringr::str_split(string = var.imp_full_df$model_name, pattern = "_", n = Inf, simplify = T)[,1])
var.imp_ensemble_df <- var.imp_ensemble_df %>% 
  mutate(unit = stringr::str_split(string = var.imp_ensemble_df$model_name, pattern = "_", n = Inf, simplify = T)[,1])

# Save all.units df
save(var.imp_full_df, file = paste0("../../controls/variable_importance/var.imp_full_df.RData"), version = "2")
saveRDS(var.imp_full_df, file = paste0("../../controls/variable_importance/var.imp_full_df.RData"), version = "2")
save(var.imp_ensemble_df, file = paste0("../../controls/variable_importance/var.imp_ensemble_df.RData"), version = "2")
saveRDS(var.imp_ensemble_df, file = paste0("../../controls/variable_importance/var.imp_ensemble_df.RData"), version = "2")


# Extract mean per unit to avoid weighting by nb of submodels
var.imp_full_df_unit <- var.imp_full_df %>% 
  group_by(unit, Variable) %>% 
  summarize(mean.var.imp = mean(Variable.importance, na.rm = T))

var.imp_ensemble_df_unit <- var.imp_ensemble_df %>% 
  group_by(unit, Variable) %>% 
  summarize(mean.var.imp = mean(Variable.importance, na.rm = T))


# Order variable factor level as increasing median importance
var.imp_full_df_unit$Variable <- reorder(var.imp_full_df_unit$Variable, var.imp_full_df_unit$mean.var.imp, median, na.rm = TRUE)
var.imp_full_df_unit$Variable <- factor(var.imp_full_df_unit$Variable, levels = rev(levels(var.imp_full_df_unit$Variable))) # Decreasing order

# Order variable factor level as increasing median importance
var.imp_ensemble_df_unit$Variable <- reorder(var.imp_ensemble_df_unit$Variable, var.imp_ensemble_df_unit$mean.var.imp, median, na.rm = TRUE)
var.imp_ensemble_df_unit$Variable <- factor(var.imp_ensemble_df_unit$Variable, levels = rev(levels(var.imp_ensemble_df_unit$Variable))) # Decreasing order

# Plot for All models
p_full <- ggplot(var.imp_full_df_unit, aes(x = Variable, y = mean.var.imp, fill = Variable)) +
  geom_boxplot() + 
  labs(title = paste0("All units confounded - before ensemble model selection"),
       y = "Variable importance")
# theme(axis.text.x = element_text(angle = -90))

# Save plots for all models
pdf(file = paste0("./controls/variable_importance/All_units_full_var.imp.pdf"), height = 7, width = 8)
print(p_full)
dev.off()


# Plot for Ensemble models
p_ensemble <- ggplot(var.imp_ensemble_df_unit, aes(x = Variable, y = mean.var.imp, fill = Variable)) +
  geom_boxplot() + 
  labs(title = paste0("All units confounded - only ensemble models"),
       y = "Variable importance")
# theme(axis.text.x = element_text(angle = -90))

# Save plots for all models
pdf(file = paste0("./controls/variable_importance/All_units_ensemble_var.imp.pdf"), height = 7, width = 8)
print(p_ensemble)
dev.off()

