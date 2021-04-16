

##### Script 7: Manual evaluation of sub.models #####

# Evaluate the quality of all submodels using optimized Jaccard indices and TSS
# Estimate the threshold that maximize those evaluation metrics

#####

# Inputs:

   # Summary table for OMU models (list.models)
   # Model results from Script 05.1

# Outputs:

   # Jaccard and TSS optimized values, with their associated threshold, for each submodel
   # PDF boxplots of evaluation metric between algorithm types
       # Globally (all OMU confounded)
       # Per OMU
  
#####  




# Effacer l'environnement
rm(list = ls())

# setwd("D:/Mael/R_projects/ithomiini_current/")

# Load libraries

library(biomod2)

# Select Env data resolution 
# res <- "5"
res <- "15"

### 1/ Create empty arrays to gather the evaluation data ####

# Load exemple OMUs/units to use as template to create the empty array that will gather the evaluation data for all units

# Load Summary table for OMU models
load(file = paste0("./input_data/list.models.RData"))
modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]

# Extact one OMU of each modeling type to encompass the two types of modeling structure
unit_complete <- as.character(modeled_OMU$Tag.model[which.max(modeled_OMU$initial_model_type == "complete")])
unit_restricted <- as.character(modeled_OMU$Tag.model[which.max(modeled_OMU$initial_model_type == "restricted")])

model.runs_complete <- readRDS(file = paste0("./models/",res,"/",unit_complete,"/model.runs.rds"))
model.runs_restricted <- readRDS(file = paste0("./models/",res,"/",unit_restricted,"/model.runs.rds"))

# Get all options for each modeling level
Algo <- union(dimnames(model.runs_complete@models.evaluation@val)[[3]], dimnames(model.runs_restricted@models.evaluation@val)[[3]])
CV_runs <- union(dimnames(model.runs_complete@models.evaluation@val)[[4]], dimnames(model.runs_restricted@models.evaluation@val)[[4]])
PA_sets <- union(dimnames(model.runs_complete@models.evaluation@val)[[5]], dimnames(model.runs_restricted@models.evaluation@val)[[5]])
eval_stats <- c("Jaccard", "TSS")

# Create an array to store evaluation for all units
all.units_evaluation <- array(dim = c(nrow(modeled_OMU),  # Nb of units/OMUs
                                   length(Algo),          # Nb of algorithms = 3
                                   length(CV_runs),       # Nb of CV runs = 3 + 1 for the "Full" runs of restricted OMU
                                   length(PA_sets),       # Nb of PsA sets = 3 or 10 => 10 max
                                   length(eval_stats)),   # Nb of evaluation stats = 2 
                              dimnames = list(unit = modeled_OMU$Tag.model,   # Unit = OMU
                                              Algo = Algo,                    # Algo names
                                              CV_runs = CV_runs,              # CV runs names + "Full"
                                              PA_sets = PA_sets,              # PA set names 
                                              eval_stats = eval_stats))       # Names of evaluation stats = Jaccard, TSS      
# Create a similar array to store the cut-off value that maximize the evaluation stat for each sub_model. Usefull for binarization
all.units_cutoffs <- all.units_evaluation

# Save empty files
save(all.units_evaluation, file = "./controls/evaluations/all.units_evaluation.RData", version = 2)
saveRDS(all.units_evaluation, file = "./controls/evaluations/all.units_evaluation.rds", version = 2)
save(all.units_cutoffs, file = "./controls/evaluations/all.units_cutoffs.RData", version = 2)
saveRDS(all.units_cutoffs, file = "./controls/evaluations/all.units_cutoffs.rds", version = 2)



### 2/ Loop to evaluate each sub_model with Jaccard and TSS ####

# Effacer l'environnement
rm(list = ls())

# setwd("D:/Mael/R_projects/ithomiini_current/")

# Load libraries
library(biomod2)

# Select Env data resolution 
# res <- "5"
res <- "15"

# Load Summary table for OMU models
load(file = paste0("./input_data/list.models.RData"))
modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]

# Load summary files to fill
all.units_evaluation <- readRDS(file = "./controls/evaluations/all.units_evaluation.rds")
all.units_cutoffs <- readRDS(file = "./controls/evaluations/all.units_cutoffs.rds")

# Set up directory
initial.wd <- getwd()
setwd(paste0(initial.wd,"/models/",res,"/"))

for (i in 1:nrow(modeled_OMU)) 
{
  
  unit <- as.character(modeled_OMU$Tag.model[i])
  model_type <- modeled_OMU$initial_model_type[i]
  
  cat(paste0("\n", Sys.time()," ------ Evaluation starts for ", unit, " = Unit N°",i,"/",nrow(modeled_OMU)," ------\n"))
  
  # Retrieve calibration lines to know which data point to use as test set for each sub_model
  model.runs <- readRDS(file = paste0(unit,"/model.runs.rds"))
  
  # Generate empty list to store results of cut-off tests
  jaccard.test.list <- TSS.test.list <- NULL
  
  # Retrieve foramted data, especially the infos on P/PsA
  formated.input.data <- readRDS(file = paste0(unit,"/formated.input.data.rds"))
  # load(model.runs@formated.input.data@link) ; formated.input.data <- data  # Retrieve formated.input.data for infos on PA runs
  input_data <- formated.input.data@data.species # Presence = 1 ; Absence = 0 ; Pseudo-absence = NA
  
  # Retrieve calib.lines. TRUE = Points used for calibration (70%) ; FALSE = Points used for evaluation (30%) ; Points not used for this submodel
  load(model.runs@calib.lines@link)          
  
  ### Explore calib.lines ; 3D array => calib.lines[data, CV_run, PA_set]
  # dim(calib.lines)
  # dimnames(calib.lines)
  # length(which(calib.lines[ , 1, 1]))         # Points for calibration (70% or 100%)
  # length(which(!calib.lines[ , 1, 1]))        # Points for evaluation (30% or 0%)
  # length(which(is.na(calib.lines[ , 1, 1])))  # Others PsA not used in this sub_model
  
  # Load models.prediction ; 4D array saving predicts => models.prediction[points, Algo, CV_run, PA_set]
  load(model.runs@models.prediction@link)
  # Predicts go from 0 to 1000
  
  # dimnames(models.prediction)
  
  get_built_models(model.runs)
  models.prediction
  
  pa.runs <- dimnames(models.prediction)[[4]] # Si on a pas généré des PsA car on avait déjà des absences, on attribue le nom que BioMOd attribue automatiquement à ce type de run, soit "AllData"
  cv.runs <- dimnames(models.prediction)[[3]]
  algos <- dimnames(models.prediction)[[2]]
  
  
  
  for (pa in pa.runs) # Loop for PA_set
  {
    # pa <- pa.runs[1]
    
    for (cv in cv.runs)  # Loop for CV_run
    {
      # cv <- cv.runs[1]

      if (model_type == "complete") # Model type = complete. Evaluation is done on 30% of data, not used for calibration
      { 
        # Retrieve input values of data use for evaluation
        obs.data <- input_data[which(!calib.lines[,                     # We take the "F" lines = lines for test set
                                                  paste0("_", cv),      # Need to add "_" to the names to fit the names in calib.lines
                                                  paste0("_", pa)])]    
        obs.data[is.na(obs.data)] <- 0 # /!\ Transformation of pseudo-absences in absences needed to generate the confusion matrix
        
        # Retrieve lines used in the model. # T = calibration/learning set. # F = evaluation/test set
        cur.eval <- calib.lines[, 
                                paste0("_", cv), 
                                paste0("_", pa)][which(!is.na(calib.lines[,   # Get only the non-NA lines = the ones used for these models
                                                                          paste0("_", cv), 
                                                                          paste0("_", pa)]))]
        
      } else { # Model type = restricted. Evaluation is done on 100% of data, already used for calibration (because N is too low for CV)
        
        # Retrieve input values of data use for evaluation
        obs.data <- input_data[which(calib.lines[,                      # We take the "T" lines = lines for calibration set but that we will use for evaluation, because anyway, we already used all of them. (No lines put aside for evaluation because of low N)
                                                  paste0("_", cv),      # Need to add "_" to the names to fit the names in calib.lines
                                                  paste0("_", pa)])]    
        obs.data[is.na(obs.data)] <- 0 # /!\ Transformation of pseudo-absences in absences needed to generate the confusion matrix
        
        # Retrieve lines used in the model. # T = calibration/learning set. # F = evaluation/test set
        # Since we already used all of them for calibration, we turn them into F to use them for evaluation too
        cur.eval <- !calib.lines[, 
                                paste0("_", cv), 
                                paste0("_", pa)][which(!is.na(calib.lines[,   # Get only the non-NA lines = the ones used for these models
                                                                          paste0("_", cv), 
                                                                          paste0("_", pa)]))]
      }
      
      for (algo in algos) # Loop for algorithms
      {
        # algo <- algos[1]
        
        # We use line dedicated for evaluation : !cur.eval = lines for evaluation in cur.eval (F)
        # models.prediction = Predicts for each sub_model
        cur.preds <- models.prediction[, algo, cv, pa][!cur.eval]  # We extract only the ones for evaluation
        
        if (!all(is.na(cur.preds))) # Sub_models that did NOT converge, do not have predictions, and do not need to be evaluated
        {
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
          
          jaccard.test.list[[paste0(unit, "_", pa, "_", cv, "_", algo)]] <- jaccard.test
          TSS.test.list[[paste0(unit, "_", pa, "_", cv, "_", algo)]] <- TSS.test
          
          # Find the best cut-off = threshold that maximize the stat
          # We take a mean because several values can lead to equal stat
          all.units_cutoffs[unit, algo, cv, pa, "Jaccard"] <- mean(jaccard.test$cutoff[which(jaccard.test$jaccard == max(jaccard.test$jaccard))])
          all.units_cutoffs[unit, algo, cv, pa, "TSS"] <- mean(TSS.test$cutoff[which(TSS.test$TSS == max(TSS.test$TSS))])
          
          # Save the eval stat value. Will be used to discard bad models prior ensemble, and eventually for binarization
          all.units_evaluation[unit, algo, cv, pa, "Jaccard"] <- jaccard.test$jaccard[
            which(jaccard.test$cutoff == round(all.units_cutoffs[unit, algo, cv, pa, "Jaccard"]))] # We round the mean to get an existing cut-off value
          
          all.units_evaluation[unit, algo, cv, pa, "TSS"] <- TSS.test$TSS[
            which(TSS.test$cutoff == round(all.units_cutoffs[unit, algo, cv, pa, "TSS"]))] # We round the mean to get an existing cut-off value
        
        }
      }
    }
  }
  
  # Save unit level arrays with all the eval and cut-off to avoid losing everything if the loop stops!
  array_cutoffs <- all.units_cutoffs[unit, , , , ]
  array_evals <- all.units_evaluation[unit, , , , ]
  
  save(array_evals, file = paste0("../../controls/evaluations/By_OMU/eval_results/array_evals_",unit,".RData"), version = "2")
  saveRDS(array_evals, file = paste0("../../controls/evaluations/By_OMU/eval_results/array_evals_",unit,".rds"), version = "2")
  save(array_cutoffs, file = paste0("../../controls/evaluations/By_OMU/eval_results/array_cutoffs_",unit,".RData"), version = "2")
  saveRDS(array_cutoffs, file = paste0("../../controls/evaluations/By_OMU/eval_results/array_cutoffs_",unit,".rds"), version = "2")
  
  # Store cut-off test infos for this unit, for reproductibility
  save(jaccard.test.list, file = paste0("../../controls/evaluations/By_OMU/test_data/jaccard.tests_",unit,".RData"), version = "2")
  saveRDS(jaccard.test.list, file = paste0("../../controls/evaluations/By_OMU/test_data/jaccard.tests_",unit,".rds"), version = "2")
  save(TSS.test.list, file = paste0("../../controls/evaluations/By_OMU/test_data/TSS.tests_",unit,".RData"), version = "2")
  saveRDS(TSS.test.list, file = paste0("../../controls/evaluations/By_OMU/test_data/TSS.tests_",unit,".rds"), version = "2")
  
  cat(paste0("\n", Sys.time()," ------ Evaluation done for ", unit, " = Unit N°",i,"/",nrow(modeled_OMU)," ------\n"))
  
}

# setwd("D:/Mael/R_projects/ithomiini_current/")
setwd(initial.wd) # Back to initial wd

# Save cut-off data
save(all.units_cutoffs, file = "./controls/evaluations/all.units_cutoffs.RData", version = "2")
saveRDS(all.units_cutoffs, file = "./controls/evaluations/all.units_cutoffs.rds", version = "2")

# Save evaluations
save(all.units_evaluation, file = "./controls/evaluations/all.units_evaluation.RData", version = "2")
saveRDS(all.units_evaluation, file = "./controls/evaluations/all.units_evaluation.rds", version = "2")


### Explore results
# i <- 2
# unit <- as.character(modeled_OMU$Tag.model[i])
# 
# load(file = paste0("../../controls/evaluations/By_OMU/test_data/TSS.tests_",unit,".RData"))
# load(file = paste0("../../controls/evaluations/By_OMU/eval_results/array_evals_",unit,".RData"))
###


##### 3/ Graphs #####

# Effacer l'environnement
rm(list = ls())

# setwd("D:/Mael/R_projects/ithomiini_current/")

library(ggplot2)
library(tidyverse)


# Load evaluation arrays
all.units_evaluation <- readRDS(file = "./controls/evaluations/all.units_evaluation.rds")

# Turn into a df usable by ggplot
gg_jaccard <- reshape2::melt(all.units_evaluation[ , , , , "Jaccard"])
gg_TSS <- reshape2::melt(all.units_evaluation[ , , , , "TSS"])

# Remove empty lines
gg_jaccard <- gg_jaccard[complete.cases(gg_jaccard), ]
gg_TSS <- gg_TSS[complete.cases(gg_TSS), ]

# Add distinction of model_type
# Load Summary table for OMU models
load(file = paste0("./input_data/list.models.RData"))
modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]
model_type_patch <- dplyr::select(modeled_OMU, unit = Tag.model, model_type = initial_model_type)

gg_jaccard <- dplyr::left_join(x = gg_jaccard, y = model_type_patch, by = "unit")
gg_TSS <- dplyr::left_join(x = gg_TSS, y = model_type_patch, by = "unit")

# Save df of sub_model evaluations
save(gg_jaccard, file = "./controls/evaluations/Full_df_Jaccard.Rdata", version = "2")
saveRDS(gg_jaccard, file = "./controls/evaluations/Full_df_Jaccard.rds", version = "2")
save(gg_TSS, file = "./controls/evaluations/Full_df_TSS.Rdata", version = "2")
saveRDS(gg_TSS, file = "./controls/evaluations/Full_df_TSS.rds", version = "2")


### 3.1/ Global plots ####

# Load df of sub_model evaluations
gg_jaccard <- readRDS(file = "./controls/evaluations/Full_df_Jaccard.rds")
gg_TSS <- readRDS(file = "./controls/evaluations/Full_df_TSS.rds")

# Jaccard index interpretation: % of predict area intersecting observed area. 

pdf(file = paste0("./controls/evaluations/All_units_Preselect_Eval_Jaccard_boxplot.pdf"), height = 6, width = 10) 
ggplot(gg_jaccard, aes(x = Algo, y = value, fill = model_type)) +
  geom_boxplot(show.legend = T) +
  # ylim(c(0, 1)) +
  geom_hline(yintercept = 0.6, col = "red", lwd = 1, lty = 2) +
  geom_hline(yintercept = 0.95, col = "dodgerblue", lwd = 1, lty = 2) +
  labs(title = "Global evaluation of all models with Jaccard indices",
       y = "Jaccard indices",
       x = "Algorithm")
dev.off()

pdf(file = paste0("./controls/evaluations/All_units_Preselect_Eval_TSS_boxplot.pdf"), height = 6, width = 10)
ggplot(gg_TSS, aes(x = Algo, y = value, fill = model_type)) +
  geom_boxplot(show.legend = T) +
  ylim(c(0, 1)) +
  geom_hline(yintercept = 0.5, col = "red", lwd = 1, lty = 2) +
  geom_hline(yintercept = 0.95, col = "dodgerblue", lwd = 1, lty = 2) +
  labs(title = "Global evaluation of all models with TSS indices",
       y = "TSS indices",
       x = "Algorithm")
dev.off()

# Enhanced plot for Supplementaries
pdf(file = paste0("./supplementaries/Eval_Jaccard_boxplot.pdf"), height = 6, width = 10) 
gg_jaccard_plot <- ggplot(gg_jaccard, aes(x = Algo, y = value, fill = model_type)) +
  geom_boxplot(show.legend = T) +
  ylim(c(0.4, 1)) +
  geom_hline(yintercept = 0.6, col = "red", lwd = 1, lty = 2) +
  geom_hline(yintercept = 0.95, col = "dodgerblue", lwd = 1, lty = 2) +
  labs(y = "Jaccard index",
       x = "Algorithm") +
  # ggthemes::geom_rangeframe(data = data.frame(y = c(0.4, 1)), aes(y = y), sides = "l", size = 1.4) +
  scale_fill_manual(name = " Model type",
                      labels = c("Complete", "Restricted"),
                      values = c("#F8766D", "#00BFC4")) +

  theme(panel.background = element_rect(fill = "white", color = "black", size = 1.0),
        legend.position = c(0.095, 0.12),
        legend.background = element_rect(fill = "white", color = "white", size = 10),
        legend.key = element_rect(fill = "white", colour = "white"),
        legend.title = element_text(size = 14, vjust = 2, face = "bold"),
        legend.text = element_text(size = 12, face = "bold"),
        # margin = margin(t = 0, unit = "pt")),
        # legend.key.size = unit(1.5, 'lines'),
        # legend.spacing.y = unit(1,"cm"),
        # axis.line.y = element_line(size = 1, color = "black", linetype = 1),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(size = 1.2),
        axis.ticks.length.y = unit(5, "pt"),
        axis.text.y = element_text(size = 16, face = "bold", color = "black", margin = margin(r = 5)),
        axis.text.x = element_text(size = 18, face = "bold", color = "black", margin = margin(t = 10)),
        axis.title = element_text(size = 18, face = "bold"),
        # axis.title.x = element_text(margin = margin(t = 10, b = 5)),
        axis.title.x = element_blank(),
        axis.title.y = element_text(margin = margin(l = 5, r = 15))) +
  
  annotate(geom = "text", x = 1.18, y = 0.492, size = 5, fontface = "bold", 
           label = "Threshold", color = "black") +
  annotate(geom = "segment", x = 1.0, xend = 1.4, y = 0.445, yend = 0.445,
           size = 1.3, linetype = 2, color = "red") +
  annotate(geom = "segment", x = 1.0, xend = 1.4, y = 0.417, yend = 0.417,
           size = 1.3, linetype = 2, color = "dodgerblue")
print(gg_jaccard_plot)
dev.off()


### 3.2/ Plot per OMU/unit ####

# Load Summary table for OMU models
load(file = paste0("./input_data/list.models.RData"))
modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]

for (i in 1:nrow(modeled_OMU))
{
  # i <- 1
  
  unit <- as.character(modeled_OMU$Tag.model[i])
  unit_model_type <- modeled_OMU$initial_model_type[i]
  
  gg_jaccard_unit <- gg_jaccard[gg_jaccard$unit == unit, ]
  
  pdf(file = paste0("./controls/evaluations/By_OMU/boxplots/Jaccard/Preselect_Eval_Jaccard_boxplot_",unit,".pdf"), height = 6, width = 10)
  g <- ggplot(gg_jaccard_unit, aes(x = Algo, y = value, fill = Algo)) +
    geom_boxplot(show.legend = F) +
    ylim(c(min(gg_jaccard_unit$value)*0.9, 1)) +
    geom_hline(yintercept = c(0.6, 0.95)[c("complete", "restricted") == unit_model_type], col = "red", lwd = 1, lty = 2) +
    labs(title = paste0("Evaluation of models for ",unit," with Jaccard indices"),
         subtitle = paste0("Model type = ",unit_model_type),
         y = "Jaccard indices",
         x = "Algorithm")
  print(g)
  dev.off()
  
  gg_TSS_unit <- gg_TSS[gg_TSS$unit == unit, ]
  
  pdf(file = paste0("./controls/evaluations/By_OMU/boxplots/TSS/Preselect_Eval_TSS_boxplot_",unit,".pdf"), height = 6, width = 10)
  g <- ggplot(gg_TSS_unit, aes(x = Algo, y = value, fill = Algo)) +
    geom_boxplot(show.legend = F) +
    ylim(c(min(gg_TSS_unit$value)*0.9, 1)) +
    geom_hline(yintercept = c(0.5, 0.95)[c("complete", "restricted") == unit_model_type], col = "red", lwd = 1, lty = 2) +
    labs(title = paste0("Evaluation of models for ",unit," with TSS indices"),
         subtitle = paste0("Model type = ",unit_model_type),
         y = "TSS indices",
         x = "Algorithm")
  print(g)
  dev.off()
  
  if (i %% 10 == 0) {
    print(i)
  }
  
}

### 4/ Quick stats ####

# Load df of sub_model evaluations


gg_jaccard <- readRDS(file = "./controls/evaluations/Full_df_Jaccard.rds")

# For RF models

RF_complete <- subset(gg_jaccard, subset = (Algo == "RF") & (model_type == "complete"))$value
round(sum(RF_complete > 0.6) / length(RF_complete) * 100, 1)

RF_restricted <- subset(gg_jaccard, subset = (Algo == "RF") & (model_type == "restricted"))$value
round(sum(RF_restricted > 0.95) / length(RF_restricted) * 100, 1)

(sum(RF_complete > 0.6) + sum(RF_restricted > 0.95)) / (length(RF_complete) + length(RF_restricted))

# For GBM models

GBM_complete <- subset(gg_jaccard, subset = (Algo == "GBM") & (model_type == "complete"))$value
round(sum(GBM_complete > 0.6) / length(GBM_complete) * 100, 1)

GBM_restricted <- subset(gg_jaccard, subset = (Algo == "GBM") & (model_type == "restricted"))$value
round(sum(GBM_restricted > 0.95) / length(GBM_restricted) * 100, 1)

(sum(GBM_complete > 0.6) + sum(GBM_restricted > 0.95)) / (length(GBM_complete) + length(GBM_restricted))

# For ANN models

ANN_complete <- subset(gg_jaccard, subset = (Algo == "ANN") & (model_type == "complete"))$value
round(sum(ANN_complete > 0.6) / length(ANN_complete) * 100, 1)

ANN_restricted <- subset(gg_jaccard, subset = (Algo == "ANN") & (model_type == "restricted"))$value
round(sum(ANN_restricted > 0.95) / length(ANN_restricted) * 100, 1)

(sum(ANN_complete > 0.6) + sum(ANN_restricted > 0.95)) / (length(ANN_complete) + length(ANN_restricted))

