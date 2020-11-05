
##### Script 8: Selection of suitable sub_models to use in the ensemble #####

# Choice of threshold to discard submodels with limited quality based on evaluation metrics

# Plot submodels response curves with warnings to visually check response curve credibility and discard them or not

# Register information on evaluation metrics and nb of valid sub_models at OMU level

# Actualize response plots with only the selected submodels


#####

# Inputs:
   # Evaluation data from Script 07
   # Summary table for Response Plots with warnings associated, from Script 06 (Response_plot_summary_table)

# Outputs:
   # Evaluation metric thresholds for Jaccard & TSS to discard submodel from ensemble, for complete and restricted models
   # Histogram of nb of valid submodels per OMU depedning on the threshold chosen
   # List of valid submodels for each OMU

#####


# Effacer l'environnement
rm(list = ls())

# Load libraries

library(biomod2)
library(tidyverse)

### 1/ Selection based on the evaluation metric

### Make post-selection boxplot. Global and per_units. Not now ! Wait for the ensemble to be made and add the EM eval in the plot.


### 1.1/ Select the proper threshold for Jaccard index based on histogram of number of models selected per OMU ####

gg_jaccard <- readRDS(file = "./controls/evaluations/Full_df_Jaccard.rds")
gg_TSS <- readRDS(file = "./controls/evaluations/Full_df_TSS.rds")

# To add distinction of model_type
# Load Summary table for OMU models
load(file = paste0("./input_data/list.models.RData"))
modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]
model_type_patch <- dplyr::select(modeled_OMU, unit = Tag.model, model_type = initial_model_type)

# Global
Eval_per_unit_Jaccard <- gg_jaccard %>% 
  # mutate(model_name = paste0(unit, "_", PA_sets, "_", CV_runs, "_", Algo)) %>% 
  mutate(select_0.5 = (value > 0.5), select_0.6 = (value > 0.6), select_0.7 = (value > 0.7), select_0.9 = (value > 0.9), select_0.95 = (value > 0.95), select_0.975 = (value > 0.975)) %>% 
  group_by(unit) %>%
  summarize(mean_jaccard = mean(value, na.rm = T), N_select_0.5 = sum(select_0.5), N_select_0.6 = sum(select_0.6), N_select_0.7 = sum(select_0.7), N_select_0.9 = sum(select_0.9), N_select_0.95 = sum(select_0.95), N_select_0.975 = sum(select_0.975)) %>% 
  left_join(model_type_patch, by = "unit")

ggplot(Eval_per_unit_Jaccard, aes(x = N_select_0.5, fill = model_type)) +
  geom_histogram(position = position_dodge())

table(Eval_per_unit_Jaccard$N_select_0.6)
ggplot(Eval_per_unit_Jaccard, aes(x = N_select_0.6, fill = model_type)) +
  geom_histogram(position = position_dodge())

table(Eval_per_unit_Jaccard$N_select_0.7)
ggplot(Eval_per_unit_Jaccard, aes(x = N_select_0.7, fill = model_type)) +
  geom_histogram(position = position_dodge())

# For complete models with independant test set (and spatial CV)
Eval_per_unit_Jaccard_complete <- Eval_per_unit_Jaccard %>% 
  filter(model_type == "complete")

table(Eval_per_unit_Jaccard_complete$N_select_0.5)
ggplot(Eval_per_unit_Jaccard_complete, aes(x = N_select_0.5, fill = model_type)) +
  geom_histogram()

### Selected threshold = 0.6
table(Eval_per_unit_Jaccard_complete$N_select_0.6)
pdf(file = paste0("./controls/evaluations/Histo_Eval_complete_Jaccard_0.6.pdf"), height = 6, width = 10)
ggplot(Eval_per_unit_Jaccard_complete, aes(x = N_select_0.6, fill = model_type)) +
  geom_histogram() +
  xlim(c(0,30)) +
  xlab("Nb of valid sub_models") +
  ylab("Nb of OMUs")
dev.off()

table(Eval_per_unit_Jaccard_complete$N_select_0.7)
ggplot(Eval_per_unit_Jaccard_complete, aes(x = N_select_0.7, fill = model_type)) +
  geom_histogram()

# For restricted models with test set = learning set. No independance !
Eval_per_unit_Jaccard_restricted <- Eval_per_unit_Jaccard %>% 
  filter(model_type == "restricted")

table(Eval_per_unit_Jaccard_restricted$N_select_0.9)
ggplot(Eval_per_unit_Jaccard_restricted, aes(x = N_select_0.9, fill = model_type)) +
  geom_histogram() +
  scale_fill_manual(values = "#00BFC4")

### Selected threshold = 0.95
table(Eval_per_unit_Jaccard_restricted$N_select_0.95)
pdf(file = paste0("./controls/evaluations/Histo_Eval_restricted_Jaccard_0.95.pdf"), height = 6, width = 10)
ggplot(Eval_per_unit_Jaccard_restricted, aes(x = N_select_0.95, fill = model_type)) +
  geom_histogram() +
  scale_fill_manual(values = "#00BFC4") +
  xlim(c(0,30)) +
  xlab("Nb of valid sub_models") +
  ylab("Nb of OMUs")
dev.off()

table(Eval_per_unit_Jaccard_restricted$N_select_0.975)
ggplot(Eval_per_unit_Jaccard_restricted, aes(x = N_select_0.975, fill = model_type)) +
  geom_histogram() +
  scale_fill_manual(values = "#00BFC4")

ggplot(Eval_per_unit_Jaccard_restricted, aes(x = mean_jaccard, fill = model_type)) +
  geom_histogram() +
  scale_fill_manual(values = "#00BFC4")

### 1.2/ Select the proper threshold for TSS index based on histogram of number of models selected per OMU ####

# Global
Eval_per_unit_TSS <- gg_TSS %>% 
  # mutate(model_name = paste0(unit, "_", PA_sets, "_", CV_runs, "_", Algo)) %>% 
  mutate(select_0.5 = (value > 0.5), select_0.6 = (value > 0.6), select_0.7 = (value > 0.7), select_0.9 = (value > 0.9), select_0.95 = (value > 0.95), select_0.975 = (value > 0.975)) %>% 
  group_by(unit) %>%
  summarize(mean_TSS = mean(value, na.rm = T), N_select_0.5 = sum(select_0.5), N_select_0.6 = sum(select_0.6), N_select_0.7 = sum(select_0.7), N_select_0.9 = sum(select_0.9), N_select_0.95 = sum(select_0.95), N_select_0.975 = sum(select_0.975)) %>% 
  left_join(model_type_patch, by = "unit")

table(Eval_per_unit_TSS$N_select_0.5)
ggplot(Eval_per_unit_TSS, aes(x = N_select_0.5, fill = model_type)) +
  geom_histogram(position = position_dodge())

table(Eval_per_unit_TSS$N_select_0.6)
ggplot(Eval_per_unit_TSS, aes(x = N_select_0.6, fill = model_type)) +
  geom_histogram(position = position_dodge())

table(Eval_per_unit_TSS$N_select_0.7)
ggplot(Eval_per_unit_TSS, aes(x = N_select_0.7, fill = model_type)) +
  geom_histogram(position = position_dodge())

# For complete models with independant test set (and spatial CV)
Eval_per_unit_TSS_complete <- Eval_per_unit_TSS %>% 
  filter(model_type == "complete")

### Selected threshold = 0.5
table(Eval_per_unit_TSS_complete$N_select_0.5)
pdf(file = paste0("./controls/evaluations/Histo_Eval_complete_TSS_0.5.pdf"), height = 6, width = 10)
ggplot(Eval_per_unit_TSS_complete, aes(x = N_select_0.5, fill = model_type)) +
  geom_histogram() +
  xlim(c(0,30)) +
  xlab("Nb of valid sub_models") +
  ylab("Nb of OMUs")
dev.off()

table(Eval_per_unit_TSS_complete$N_select_0.6)
ggplot(Eval_per_unit_TSS_complete, aes(x = N_select_0.6, fill = model_type)) +
  geom_histogram()

table(Eval_per_unit_TSS_complete$N_select_0.7)
ggplot(Eval_per_unit_TSS_complete, aes(x = N_select_0.7, fill = model_type)) +
  geom_histogram()

# For restricted models with test set = learning set. No independance !
Eval_per_unit_TSS_restricted <- Eval_per_unit_TSS %>% 
  filter(model_type == "restricted")

table(Eval_per_unit_TSS_restricted$N_select_0.9)
ggplot(Eval_per_unit_TSS_restricted, aes(x = N_select_0.9, fill = model_type)) +
  geom_histogram() +
  scale_fill_manual(values = "#00BFC4")

### Selected threshold = 0.95
table(Eval_per_unit_TSS_restricted$N_select_0.95)
pdf(file = paste0("./controls/evaluations/Histo_Eval_restricted_TSS_0.95.pdf"), height = 6, width = 10)
ggplot(Eval_per_unit_TSS_restricted, aes(x = N_select_0.95, fill = model_type)) +
  geom_histogram() +
  scale_fill_manual(values = "#00BFC4") +
  xlim(c(0,30)) +
  xlab("Nb of valid sub_models") +
  ylab("Nb of OMUs")
dev.off()

table(Eval_per_unit_TSS_restricted$N_select_0.975)
ggplot(Eval_per_unit_TSS_restricted, aes(x = N_select_0.975, fill = model_type)) +
  geom_histogram() +
  scale_fill_manual(values = "#00BFC4")

ggplot(Eval_per_unit_TSS_restricted, aes(x = mean_TSS, fill = model_type)) +
  geom_histogram() +
  scale_fill_manual(values = "#00BFC4")

### 1.3/ Registered valid sub_models ####

library(tidyverse)

gg_jaccard <- readRDS(file = "./controls/evaluations/Full_df_Jaccard.rds")
gg_TSS <- readRDS(file = "./controls/evaluations/Full_df_TSS.rds")

gg_jaccard <- gg_jaccard %>% 
  mutate(model_name = paste0(unit, "_", PA_sets, "_", CV_runs, "_", Algo))
gg_TSS <- gg_TSS %>% 
  mutate(model_name = paste0(unit, "_", PA_sets, "_", CV_runs, "_", Algo))

# Create df to summarize results of evaluation of all sub_models
Sub_models_evaluations <- left_join(x = gg_jaccard, gg_TSS, by = "model_name") %>% 
  select(model_name = model_name, unit = unit.x, PA_sets = PA_sets.x, CV_runs = CV_runs.x, model_type = model_type.x, Algo = Algo.x, Jaccard = value.x, TSS = value.y) %>% 
  mutate(Eval_Jaccard = ((model_type == "complete") & (Jaccard > 0.6)) | ((model_type == "restricted") & (Jaccard > 0.95)),
         Eval_TSS = ((model_type == "complete") & (TSS > 0.5)) | ((model_type == "restricted") & (TSS > 0.95)))

# plot(Sub_models_evaluations$Jaccard ~ Sub_models_evaluations$TSS)

save(Sub_models_evaluations, file = "./controls/evaluations/Sub_models_evaluations.RData", version = "2")
saveRDS(Sub_models_evaluations, file = "./controls/evaluations/Sub_models_evaluations.rds", version = "2")


### 2/ Additionnal control based on response plots ecological credibility ####

# Load Summary table for response plots
Response_plot_summary_table <- readRDS(file = "./controls/response_plots/Response_plot_summary_table.rds")

library(tidyverse)

# Count number of warnings per sub_models
Warning_summary_per_submodel <- Response_plot_summary_table %>% 
  group_by(model_name) %>%  # Group by sub_model
  summarize(n_warnings = sum(warning, na.rm = TRUE))

table(Warning_summary_per_submodel$n_warnings)
hist(Warning_summary_per_submodel$n_warnings)

# Count number of warnings per OMU
# Warning_summary_per_unit <- Response_plot_summary_table %>% 
#   group_by(unit) %>%  # Group features by continent
#   summarize(n_warnings = sum(warning, na.rm = TRUE))
# 
# table(Warning_summary_per_unit$n_warnings)
# hist(Warning_summary_per_unit$n_warnings)

# Merge info with evaluation to create list.submodels and check only submodels with at least two warnings that passed the evaluation
Sub_models_evaluations <- readRDS(file = "./controls/evaluations/Sub_models_evaluations.rds")
Response_plot_summary_table <- readRDS(file = "./controls/response_plots/Response_plot_summary_table.rds")
list.submodels <- unique(Response_plot_summary_table[, c("ID_sub_model", "model_name", "ID_unit", "unit")])
list.submodels <- list.submodels %>% 
  left_join(y = Sub_models_evaluations, by = c("model_name", "unit")) %>% 
  left_join(y = Warning_summary_per_submodel, by = c("model_name"))


# Select models to check with at least two warnings for ecological credibility, and that pass the Eval cut-off

list.submodels <- list.submodels %>% 
  mutate(to_check = (n_warnings >= 2) & (Eval_Jaccard|Eval_TSS))

table(list.submodels$to_check)

submodels_to_check <- list.submodels$model_name[list.submodels$to_check]
unit_for_submodels_to_check <- list.submodels$unit[list.submodels$to_check]

# Generate a pdf of response plots of models to check
pdf(file = paste0("./controls/response_plots/Models_to_check.pdf"), height = 6, width = 10) 
for (i in 1:length(submodels_to_check)) {
  unit <- unit_for_submodels_to_check[i]
  model_name <- submodels_to_check[i]
  p_sub_model <- readRDS(file = paste0("./controls/response_plots/By_OMU/",unit,"/By_model/Data/Facetted_response_plot_",model_name,".rds"))
  p_sub_model <- p_sub_model + labs(caption = paste0("N°",i)) # Add caption with the submodel index
  print(p_sub_model)
}
dev.off()

# Check manually if models are valid or not before ensemble and register model to discard in the list.submodel file

# Index of submodels to discard after visual check of ecological credibility
submodels_to_discard <- c(1:2, 5:7, 9:10, 21:24, 26, 33, 35:36, 40:41, 46, 48, 53, 57, 60, 62, 76)
submodels_to_discard <- submodels_to_check[submodels_to_discard]

# Record submodels that can be used for Ensemble based on eval metric
list.submodels$Ensemble_OK_Jaccard <- list.submodels$Eval_Jaccard
list.submodels$Ensemble_OK_TSS <- list.submodels$Eval_TSS

# Remove the ones not passing the ecological crediblity check
list.submodels$Ensemble_OK_Jaccard[which(list.submodels$model_name %in% submodels_to_discard)] <- F
list.submodels$Ensemble_OK_TSS[which(list.submodels$model_name %in% submodels_to_discard)] <- F

# save(list.submodels, file = "./input_data/list.submodels.RData", version = "2")
# saveRDS(list.submodels, file = "./input_data/list.submodels.rds", version = "2")


### 3/ Register information on evaluation metrics and nb of valid sub_models at OMU level ####

# list.submodels <- readRDS(file = "./input_data/list.submodels.rds")

list.models_patch_full <- list.submodels %>% 
  group_by(unit) %>% 
  summarize(mean_TSS = round(mean(TSS), 3), mean_Jaccard = round(mean(Jaccard), 3))

list.models_patch_Jaccard <- list.submodels %>% 
  filter(Ensemble_OK_Jaccard == T) %>%
  group_by(unit) %>% 
  summarize(mean_TSS_Select_Jaccard = round(mean(TSS), 3), N_models_Select_Jaccard = n())

list.models_patch_TSS <- list.submodels %>% 
  filter(Ensemble_OK_TSS == T) %>%
  group_by(unit) %>% 
  summarize(mean_TSS_Select_TSS = round(mean(TSS), 3), N_models_Select_TSS = n())

# Add infos to the initial list.models object
list.models <- rename(list.models, unit = Tag.model)
list.models <- list.models %>% 
  left_join(list.models_patch_full, by = "unit") %>% 
  left_join(list.models_patch_Jaccard, by = "unit") %>% 
  left_join(list.models_patch_TSS, by = "unit")
list.models <- rename(list.models, Tag.model = unit)

# Remove old stuff from biomod2 evaluation
# list.models <- select(list.models, -c("Mean_TSS_full", "Models_select_TSS_0.5", "Mean_TSS_select_TSS_0.5"))

# save(list.models, file = "./input_data/list.models.RData", version = "2")
# saveRDS(list.models, file = "./input_data/list.models.rds", version = "2")



##### 4/ Actualize response plots at Unit/OMU level to display only valid submodels used in the EM #####

# Effacer l'environnement
rm(list = ls())

# Chose your evaluation metric
eval_metric <- "Jaccard"
# eval_metric <- "TSS"

# Load list.models
load(file = paste0("./input_data/list.models.RData"))
modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]

# Load submodels summary table
list.submodels <-  readRDS(file = "./input_data/list.submodels.rds")

### Portion to modify for cluster (Not parallelized)
# index_model_to_compute <- 1:nrow(modeled_OMU)
 index_model_to_compute <- 303:nrow(modeled_OMU)
# index_model_to_compute <- 505
###

for (i in index_model_to_compute) 
{
  # for (i in 1:nrow(modeled_OMU)) {
  # i <- 1
  
  # OMU/Unit level infos
  
  unit <- as.character(modeled_OMU$Tag.model[i])
  sp <- as.character(modeled_OMU$Sp_full[i])
  
  cat(paste0("\n", Sys.time()," ------ Initialised for ", unit, " = Unit N°",i,"/",nrow(modeled_OMU)," ------\n"))
  
  resp.df <- readRDS(file = paste0("./controls/response_plots/By_OMU/",unit,"/resp.df_",unit,".rds"))
  
  # Extract only selected models for the Ensemble
  selected_submodels <- list.submodels$model_name[(list.submodels$unit == unit) & (eval(parse(text = paste0("list.submodels$Ensemble_OK_",eval_metric))))]
  resp.df_selected <-  resp.df[resp.df$Model %in% selected_submodels, ]
  
  ### 4.1/ Global response plots for this unit/OMU. Post-selection (only models used for Ensemble) ####
  
  p_unit_postselect <- ggplot(resp.df_selected, aes(x = Var.value, y = Response)) + 
    geom_line(alpha = 0.2, aes(group = Model)) + # Plot lines for each model
    stat_smooth() + # Generate GAM smoothed curve for each mini-plot
    facet_wrap(~Variable, scales = "free_x") +  # make a plot per variable
    theme_bw() + 
    ylim(0, 1) +
    xlab("Variable value") +
    labs(title = paste0("Pre-selection response plots for ",unit)) +
    geom_vline(aes(xintercept = Range_min), lwd = 1, col = "red") +
    geom_vline(aes(xintercept = Range_max), lwd = 1, col = "red")
  
  # Save plot for the OMU
  pdf(file = paste0("./controls/response_plots/By_OMU/",unit,"/Post_select_response_plots_",unit,".pdf"), height = 6, width = 10)
  print(p_unit_postselect)
  dev.off()
  # Save a copy in the Species folder
  # sp <- as.character(modeled_OMU$Sp_full[i])
  file.copy(from = paste0("../../controls/response_plots/By_OMU/",unit,"/Post_select_response_plots_",unit,".pdf"), to = paste0("../../controls/response_plots/By_species/",sp,"/Pre_select_response_plots_unit_",unit,".pdf"), overwrite = T)
  
  # Save plot in R format
  save(p_unit_postselect, file = paste0("./controls/response_plots/By_OMU/",unit,"/Post_select_response_plots_",unit,".Rdata"), version = "2")
  saveRDS(p_unit_postselect, file = paste0("./controls/response_plots/By_OMU/",unit,"/Post_select_response_plots_",unit,".rds"), version = "2")
  
  ### 4.2/ Make a multiple pages pdf with all selected submodels, and the global post_selected plot in the first page ####
  
  # Save plot for the OMU
  pdf(file = paste0("./controls/response_plots/By_OMU/",unit,"/All_models_selected_",eval_metric,"_response_plots_",unit,".pdf"), height = 6, width = 10)
  print(p_unit_postselect)
  for (j in 1:length(selected_submodels)) {
    model_name <- selected_submodels[j]
    p_sub_model <- readRDS(file = paste0("./controls/response_plots/By_OMU/",unit,"/By_model/Data/Facetted_response_plot_",model_name,".rds"))
    print(p_sub_model)
  }
  dev.off()
  # Save a copy in the Species folder
  # sp <- as.character(modeled_OMU$Sp_full[i])
  file.copy(from = paste0("./controls/response_plots/By_OMU/",unit,"/All_models_selected_",eval_metric,"_response_plots_",unit,".pdf"), to = paste0("./controls/response_plots/By_species/",sp,"/All_models_selected_",eval_metric,"_response_plots_",unit,".pdf"), overwrite = T)
  
  cat(paste0("\n", Sys.time()," ------ All plots generated for ", unit, " = Unit N°",i,"/",nrow(modeled_OMU)," ------\n"))
  
}



 

