

##### Script 6: Response plots & response check #####

# Effacer l'environnement
rm(list = ls())

# setwd("D:/Mael/R_projects/ithomiini_current/")

##### 1/ Compute and plot response curve per variable, per model (i.e., OMU*algo*PA*RUN) #####

# ### 1.1/ Generate summary table per variable, per model (i.e., OMU*algo*PA*RUN) ####
# 
# library(tidyverse)
# 
# ### Compute the length of the summary table as the sum of nb of model success * nb of variables used in the model, for each modeled OMU
# 
# # Load Summary table for OMU models
# load(file = paste0("./input_data/list.models.RData"))
# 
# # Select Env data resolution 
# # res <- "5"
# res <- "15"
# # Load environmental stack to control nb of variables used for modeling
# envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))
# 
# N.sub_models <- sum(list.models$Models_success, na.rm = T)
# N.vars <- nlayers(envData)
# N.resp.plots <- N.sub_models*N.vars
# 
# # Generate ID for each response plot
# ID_plot <- 1:(N.resp.plots)
# 
# # Generate the tibble
# Response_plot_summary_table <- tibble(ID_plot = ID_plot)
# 
# Response_plot_summary_table <- Response_plot_summary_table %>% 
#   add_column(ID_sub_model = integer(length = N.resp.plots),
#              model_name = character(length = N.resp.plots),
#              ID_unit = integer(length = N.resp.plots),
#              unit = character(length = N.resp.plots),
#              variable = character(length = N.resp.plots),
#              path_to_plot = character(length = N.resp.plots),
#              # plot = list(length = N.resp.plots),
#              # model.fit_x2 = list(length = N.resp.plots),
#              trend = logical(length = N.resp.plots),
#              beta_x2 = double(length = N.resp.plots),
#              p_value = double(length = N.resp.plots),
#              quad_issue = logical(length = N.resp.plots),
#              multimod = logical(length = N.resp.plots),
#              warning = logical(length = N.resp.plots))
# 
# # Fill the summary table with OMU model infos
# modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]
# 
# library(biomod2)
# 
# vars <- names(envData)
# 
# l <- m <- 1
# for (i in 1:nrow(modeled_OMU)) {
#   
#   N.sub_models <- modeled_OMU$Models_success[i]
#   
#   # OMU/Unit level infos
#   
#   ID_unit <- modeled_OMU$ID_unit[i]
#   unit <- as.character(modeled_OMU$Tag.model[i])
#   
#   model.runs <- readRDS(file = paste0("./models/",res,"/",unit,"/model.runs.rds"))
#   sub.models_names <- get_built_models(model.runs)
#   
#   for (j in 1:N.sub_models) {
#     
#     # Sub_model level
#     ID_sub_model <- m
#     model_name <- sub.models_names[j]
#     
#     m <- m+1
#     
#     for (k in 1:N.vars) {
#       
#       Response_plot_summary_table$variable[l] <- vars[k]
#       
#       Response_plot_summary_table$ID_unit[l] <- ID_unit
#       Response_plot_summary_table$unit[l] <- unit
#       Response_plot_summary_table$model_name[l] <- model_name
#       Response_plot_summary_table$ID_sub_model[l] <- ID_sub_model
#      
#       l <- l+1 
#     }
#   }
#   
#   print(i)
# }
# 
# # save(Response_plot_summary_table, file = "./controls/response_plots/Response_plot_summary_table.RData", version = "2")
# # saveRDS(Response_plot_summary_table, file = "./controls/response_plots/Response_plot_summary_table.rds", version = "2")

### 1.2/ Compute and save response plots ####

rm(list = ls())

# Load librairies
library(biomod2)
library(raster)
library(ggplot2)
# library(MuMIn)
library(QuantPsyc)
library(tidyverse)
library(diptest)
library(grid)


# ?biomod2::response.plot2
# Load modified version of response plot function from Biomod2. This version use the mean of occurrences as fixed values instead of the mean of all points, including pseudo-absences
source("./functions/response.plot_modified.R", local = TRUE) # Function = response.plot3

# Load Summary table and list.models
load(file = paste0("./input_data/list.models.RData"))
modeled_OMU <- list.models[!is.na(list.models$Model_ID), ]

# Select Env data resolution 
# res <- "5"
res <- "15"
# Load environmental stack to control nb of variables used for modeling
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))

# Load summary table for repsonse plots
Response_plot_summary_table <- readRDS(file = "./controls/response_plots/Response_plot_summary_table.rds")


initial.wd <- getwd()
setwd(paste0(initial.wd,"/models/",res,"/"))

### Portion to modify for cluster (Not parallelized)
index_model_to_compute <- 1:nrow(modeled_OMU)
index_model_to_compute <- 501:563
# index_model_to_compute <- 505
###

### Initialize line count
# l <- m <- 1

l <- which.max(Response_plot_summary_table$ID_unit == modeled_OMU$ID_unit[index_model_to_compute[1]])
# m <- Response_plot_summary_table$ID_sub_model[l]

for (i in index_model_to_compute) {
# for (i in 1:nrow(modeled_OMU)) {
  # i <- 505
  
  N.sub_models <- modeled_OMU$Models_success[i]
  
  # OMU/Unit level infos
  
  unit <- as.character(modeled_OMU$Tag.model[i])
  sp <- as.character(modeled_OMU$Sp_full[i])
  
  cat(paste0("\n", Sys.time()," ------ Initialised for ", unit, " = Unit N°",i,"/",nrow(modeled_OMU)," ------\n"))
  
  # Create directory to store model outputs by species
  if(!dir.exists(paste0("../../controls/response_plots/By_species/",sp))) { # Test if the folder exists already or not
    dir.create(paste0("../../controls/response_plots/By_species/",sp), recursive = T) # Create folder if absent
  }
  
  # Create directories to store model outputs by Unit/OMU, and by sub_models for each units/OMU
  if(!dir.exists(paste0("../../controls/response_plots/By_OMU/",unit,"/By_model/Data"))) { # Test if the folder exists already or not
    dir.create(paste0("../../controls/response_plots/By_OMU/",unit,"/By_model/Data"), recursive = T) # Create folder if absent
  }
  
  # Create directories to store model outputs by Unit/OMU, and by variable for each units/OMU
  if(!dir.exists(paste0("../../controls/response_plots/By_OMU/",unit,"/By_model/By_variable"))) { # Test if the folder exists already or not
    dir.create(paste0("../../controls/response_plots/By_OMU/",unit,"/By_model/By_variable"), recursive = T) # Create folder if absent
  }
  
  ### 1.2.1/ Retrieve data for response curves ####
  
  # Retrieve Env Stack
  load(file = paste0("../../input_data/Species_data/",res,"/Env_Stacks/Env_stack_",unit,".RData"))
  
  # Retrieve input data
  formated.input.data <- readRDS(file = paste0(unit,"/formated.input.data.rds"))
  
  # Retrieve calibrated models
  model.runs <- readRDS(file = paste0(unit,"/model.runs.rds"))
  
  sub.models_names <- get_built_models(model.runs)
  
  # Get variables name
  cur.vars <- model.runs@expl.var.names
  
  # Load individual models
  
  models.to.plot <- BIOMOD_LoadModels(model.runs)
  
  # Retrieve info for each plot using the in-built biomod2 function
  resp <- response.plot3(models = models.to.plot,           # List of names of models to include in the plot
                         Data = unit.env[[cur.vars]],        # Raster stack of the env variables
                         fixed.var.metric = "sp.mean",      # Mean, only on presence points
                         show.variables = cur.vars,         # Display all variables plots
                         run.data = formated.input.data,    # Initial data of models
                         plot = F)                          # Store the plot, do not display it
  
  # Remove loaded models
  rm(list = ls()[grep(x = ls(), pattern = unit)])
  
  ### Explore information extracted to plot our response plots
  # str(resp) 
  # dimnames(resp) # Array with 4D : resp[Each tested value, Value & Prediction, Variable, Sub_model (Algo*RUN*PA)]
  
  # Rename env. variables
  dimnames(resp)[[3]] <- vars <- c("Tmean", "Tvar", "Htot", "Hvar", "Elevation", "Forests")
  
  # resp[,(x,y),Var,Model_name]
  
  # Melt all infos in a df suitable for ggplot2
  resp.df <- cbind(reshape2::melt(resp[, "Var", , ]), # Retrieve tested values (x-axis)
                   reshape2::melt(resp[, "Pred", , ])$value) # Retrieve associated predictions (y-axis)
  colnames(resp.df) <- c("Index", "Variable", "Model", "Var.value", "Response") # Rename columns
  
  # Change units
  resp.df$Var.value[resp.df$Variable == "Tmean"] <- resp.df$Var.value[resp.df$Variable == "Tmean"] / 10   # Tmean in °C
  resp.df$Var.value[resp.df$Variable == "Tvar"] <- resp.df$Var.value[resp.df$Variable == "Tvar"] / 1000  # Tvar in °C
  resp.df$Var.value[resp.df$Variable == "Htot"] <- resp.df$Var.value[resp.df$Variable == "Htot"] / 1000   # Htot in % water mass/air mass
  resp.df$Var.value[resp.df$Variable == "Hvar"] <- resp.df$Var.value[resp.df$Variable == "Hvar"] / 1000   # Hvar in % of variation (CV * 100)
  
  ### Step to filter out non_selected models to built unit response curve for the EM (after selection)
  #
  ###
  
  # Compute range of each variable, among presence data
  ranges <- as.data.frame(apply(X = formated.input.data@data.env.var[!is.na(formated.input.data@data.species), ], MARGIN = 2, FUN = range))
  names(ranges) <- c("Tmean", "Tvar", "Htot", "Hvar", "Elevation", "Forests")
  
  # Change units
  ranges$Tmean <- ranges$Tmean / 10   # Tmean in °C
  ranges$Tvar <- ranges$Tvar / 1000  # Tvar in °C
  ranges$Htot <- ranges$Htot / 1000   # Htot in % water mass/air mass
  ranges$Hvar <- ranges$Hvar / 1000   # Hvar in % of variation (CV * 100)
  
  # Rearrange range data and merge them with the resp.df
  ranges_df <- ranges %>% 
    tidyr::pivot_longer(cols = names(ranges),        # List of the name of the "columns" that are actually factor levels
                        names_to = "Variable",       # name of the new column that is the real variable for these factor level
                        values_to = c("Ranges")) %>%  # name of the new columns that will get the values in the cells of the previous dataset
    add_column(Range_type = c(rep("Range_min", nrow(.)/2), rep("Range_max", nrow(.)/2))) %>%  # Add column for type of range
    tidyr::pivot_wider(names_from = "Range_type", # List of the name of the "columns" that are contains actually name of real variables 
                       values_from = "Ranges") # name of the columns that contains the values of the new columns to create from "names_from"
  resp.df <- resp.df %>% 
    left_join(y = ranges_df, by = "Variable")
  
  save(resp.df, file = paste0("../../controls/response_plots/By_OMU/",unit,"/resp.df_",unit,".RData"), version = "2")
  saveRDS(resp.df, file = paste0("../../controls/response_plots/By_OMU/",unit,"/resp.df_",unit,".rds"), version = "2")
  
  ### 1.2.2/ Global response plots for this unit/OMU. Pre-selection (no models are discarded on the basis of evaluation or weird response curve) ####
  
  p_unit <- ggplot(resp.df, aes(x = Var.value, y = Response)) + 
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
  pdf(file = paste0("../../controls/response_plots/By_OMU/",unit,"/Pre_select_response_plots_",unit,".pdf"), height = 6, width = 10)
  print(p_unit)
  dev.off()
  # Save a copy in the Species folder
  # sp <- as.character(modeled_OMU$Sp_full[i])
  file.copy(from = paste0("../../controls/response_plots/By_OMU/",unit,"/Pre_select_response_plots_",unit,".pdf"), to = paste0("../../controls/response_plots/By_species/",sp,"/Pre_select_response_plots_unit_",unit,".pdf"), overwrite = T)
  
  # Save plot in R format
  save(p_unit, file = paste0("../../controls/response_plots/By_OMU/",unit,"/Pre_select_response_plots_",unit,".Rdata"), version = "2")
  saveRDS(p_unit, file = paste0("../../controls/response_plots/By_OMU/",unit,"/Pre_select_response_plots_",unit,".rds"), version = "2")
  
  ### 2/ By sub_model, to check ecological credibility of the response curve ####
  for (j in 1:N.sub_models) {
    # j <- 26
    
    # Sub_model level
    # ID_sub_model <- m
    model_name <- Response_plot_summary_table$model_name[l]
    
    resp.df_sub.model <- resp.df[as.character(resp.df$Model) == model_name, ]
    # save(resp.df_sub.model, file = paste0("../../controls/response_plots/By_OMU/",unit,"/Data_for_Model_response_plots_",model_name,".RData"), version = "2") 
    
    ### 2.1/ Global facetted response plots for one sub-model (one mini-plot per variable) ####
    
    p_sub_model <- ggplot(resp.df_sub.model, aes(x = Var.value, y = Response)) + 
      geom_line(alpha = 0.8, lwd = 1.1, aes(group = Model)) + # Plot lines for each model
      stat_smooth() + # Generate GAM smoothed curve for each mini-plot
      facet_wrap(~Variable, scales = "free_x") +  # make a plot per variable
      theme_bw() + 
      ylim(0, 1) +
      xlab("Variable value") +
      labs(title = paste0("Response plots for ",unit),
           subtitle = paste0("Model = ",model_name)) +
      geom_vline(aes(xintercept = Range_min), lwd = 1, col = "red") +
      geom_vline(aes(xintercept = Range_max), lwd = 1, col = "red")
    
    # # Save plot for the sub-model
    # pdf(file = paste0("../../controls/response_plots/By_OMU/",unit,"/Model_response_plots_",model_name,".pdf"), height = 6, width = 10)
    # print(p_sub_model)
    # dev.off()
    
    # Save the plot in R format
    save(p_sub_model, file = paste0("../../controls/response_plots/By_OMU/",unit,"/By_model/Data/Facetted_response_plot_",model_name,".Rdata"), version = "2")
    saveRDS(p_sub_model, file = paste0("../../controls/response_plots/By_OMU/",unit,"/By_model/Data/Facetted_response_plot_",model_name,".rds"), version = "2")
    
    
    # m <- m+1
    
    ### 2.2/ Plot at variable level, for more details ####
    
    # Create directories to store model outputs by Unit/OMU, and by models for each units/OMU, for unique Variable plots
    if(!dir.exists(paste0("../../controls/response_plots/By_OMU/",unit,"/By_model/By_variable/",model_name))) { # Test if the folder exists already or not
      dir.create(paste0("../../controls/response_plots/By_OMU/",unit,"/By_model/By_variable/",model_name), recursive = T) # Create folder if absent
    }
    
    for (k in 1:length(vars)) {
      # k <- 4
      
      Response_plot_summary_table$variable[l] <- vars[k]
      
      resp.df_sub.model_var <- resp.df_sub.model[resp.df_sub.model$Variable == vars[k], ]
      # # For Model fit and test for multimodality restricited to the range of presence data
      # min <- resp.df_sub.model_var$Range_min[1] ; max <- resp.df_sub.model_var$Range_max[1]
      # resp.df_sub.model_var_in_range <- resp.df_sub.model_var[(resp.df_sub.model_var$Var.value > min) & (resp.df_sub.model_var$Var.value < max), ]
      # Var_values <- resp.df_sub.model_var_in_range$Var.value
      # Predicts <- resp.df_sub.model_var_in_range$Response
      Var_values <- resp.df_sub.model_var$Var.value
      Predicts <- resp.df_sub.model_var$Response
      
      # Check linear trend
      model.fit_x <- lm(Predicts ~ Var_values)
      summary_model.fit_x <- summary(model.fit_x)
      if(nrow(coefficients(summary_model.fit_x) > 1)) { # Check if there is a coefficient to extract
        p_value_x <- round(coefficients(summary_model.fit_x)[2,4],3)
        trend <- p_value_x < 0.001
      } else { # If not (i.e., predict is flat), provide values manually
        p_value_x <- 1
        trend <- F
      }
      
      
      # Check for quadratic effect (a positive one would be weird ecologically)
      model.fit_x2 <- lm(Predicts ~ Var_values + I(Var_values^2))
      summary_model.fit_x2 <- summary(model.fit_x2)
      if(nrow(coefficients(summary_model.fit_x2) > 1)) { # Check if there is a coefficient to extract
        # summary_model.fit_x2_std <- MuMIn::std.coef(x = model.fit_x2, partial.sd = F) # Not available for R 3.4.2 on cluster
        summary_model.fit_x2_std <- QuantPsyc::lm.beta(model.fit_x2)
        
        lm_pred_x2 <- predict(model.fit_x2)
        # beta_x2 <- round(summary_model.fit_x2_std[3,1],2)
        beta_x2 <- round(summary_model.fit_x2_std[2],2)
        p_value_x2 <- round(coefficients(summary_model.fit_x2)[3,4],3)
      } else { # If not (i.e., predict is flat), provide values manually
        lm_pred_x2 <- rep(coefficients(summary_model.fit_x2)[1,1], lenght(Predicts)) # Provide Intercept
        beta_x2 <- 0
        p_value_x2 <- 1
      }
      
      # Issue only for significant positive quadratic fit, and when no strong linear trend is detected ! # Issue only for p < 0.001.
      quad_issue <- (beta_x2 > 0) & (p_value_x2 < 0.001) & (p_value_x > 0.001)
      
      
      # Hartigan's dip test for unimodality
      # ?diptest::dip.test
      
      # dip_test <- diptest::dip.test(x =  Predicts)
      # print(dip_test)
      
      # Detect multimodal distribution if p.value < 0.05. But too sensitive so we set the warning for p < 0.01
      multimod <- diptest::dip.test(x =  Predicts)$p.value < 0.01
      
      # Warning only if quatratic issue with positive quadratic fit, or no quadratic fit or trend but multimodal
      warning <- quad_issue | ((p_value_x2 >= 0.05) & multimod & !trend)
      
      ### Plot
      
      # Create text label
      legend_coefs <- grobTree(textGrob(paste0("Beta-coef = ", beta_x2,
                                               "\np = ",p_value_x2,
                                               "\nTrend = ", trend,
                                               "\nQuad issue = ",as.character(quad_issue),
                                               "\nMultimodal = ",as.character(multimod),
                                               "\nWarning = ",as.character(warning)),
                                x=0.03,  y=0.80, hjust=0,
                                gp=gpar(fontsize=13, fontface = "bold")))
      
      # Need to add copy of last point to generate a vector of the same length than the aesthetics of the plot when plotting on limited range (does nothing if plotting on full range)
      x_lines <- c(Var_values, rep(last(Var_values), times = nrow(resp.df_sub.model_var) - length(Var_values)))
      y_lines <- c(lm_pred_x2, rep(last(lm_pred_x2), times = nrow(resp.df_sub.model_var) - length(lm_pred_x2)))
      
      p_var <- ggplot(resp.df_sub.model_var, aes(x = Var.value, y = Response)) +
        geom_line(lwd = 1.5) +
        # geom_line(aes(x = Var.value, y = lm_pred_x2), col = "dodgerblue", lwd = 1.5) +
        geom_line(aes(x = x_lines, y = y_lines), col = "dodgerblue", lwd = 1.5) +
        labs(title = paste0("Response plots for ",unit),
             subtitle = paste0("Model = ",model_name),
             x = vars[k]) +
        ylim(0, 1) +
        geom_vline(aes(xintercept = Range_min), lwd = 1.2, col = "red") +
        geom_vline(aes(xintercept = Range_max), lwd = 1.2, col = "red") +
        annotation_custom(legend_coefs)
      
      # Save plot for this variable
      pdf(file = paste0("../../controls/response_plots/By_OMU/",unit,"/By_model/By_variable/",model_name,"/Variable_response_plots_",model_name,"_",vars[k],".pdf"), height = 6, width = 8)
      print(p_var)
      dev.off()
      
      # Save in R Format, including the customed x_lines and y_lines for the quadratic predict on limited range
      save(p_var, x_lines, y_lines, file = paste0("../../controls/response_plots/By_OMU/",unit,"/By_model/Data/Variable_response_plots_",model_name,"_",vars[k],".RData"), version = "2")
      rm(y_lines, x_lines)
      
      # Save infos at variable level
      Response_plot_summary_table$path_to_plot[l] <- paste0("./controls/response_plots/By_OMU/",unit,"/By_model/By_variable/",model_name,"Variable_response_plots_",model_name,"_",vars[k],".pdf")
      
      # Response_plot_summary_table$model.fit_x2[[l]] <- model.fit_x2
      Response_plot_summary_table$trend[l] <- trend
      Response_plot_summary_table$beta_x2[l] <- beta_x2
      Response_plot_summary_table$p_value[l] <- p_value_x2
      Response_plot_summary_table$multimod[l] <- multimod
      Response_plot_summary_table$warning[l] <- warning
      
      l <- l+1 
    }
    
    ### 2.3/ Make a multiple pages pdf with the facetted plot, followed by the detailed plot for each variable of the sub_model ####
    pdf(file = paste0("../../controls/response_plots/By_OMU/",unit,"/By_model/Model_response_plots_",model_name,".pdf"), height = 6, width = 10)
    print(p_sub_model)
    for (k in 1:length(vars)) {
      load(file = paste0("../../controls/response_plots/By_OMU/",unit,"/By_model/Data/Variable_response_plots_",model_name,"_",vars[k],".RData"))
      print(p_var)
    }
    dev.off()
    
  }
  
  ### 2.4/ Make a multiple pages pdf with each sub_model of the OMU, and the Pre_select in last page ####
  # Save plot for the OMU
  pdf(file = paste0("../../controls/response_plots/By_OMU/",unit,"/All_models_response_plots_",unit,".pdf"), height = 6, width = 10)
  for (j in 1:N.sub_models) {
    model_name <- sub.models_names[j]
    p_sub_model <- readRDS(file = paste0("../../controls/response_plots/By_OMU/",unit,"/By_model/Data/Facetted_response_plot_",model_name,".rds"))
    print(p_sub_model)
  }
  print(p_unit)
  dev.off()
  # Save a copy in the Species folder
  # sp <- as.character(modeled_OMU$Sp_full[i])
  file.copy(from = paste0("../../controls/response_plots/By_OMU/",unit,"/All_models_response_plots_",unit,".pdf"), to = paste0("../../controls/response_plots/By_species/",sp,"/All_models_response_plots_",unit,".pdf"), overwrite = T)
  
  ####### Later, add the post-selection summary plot in first page, or better, the EM response plot !
  
  cat(paste0("\n", Sys.time()," ------ All plots generated for ", unit, " = Unit N°",i,"/",nrow(modeled_OMU)," ------\n"))
}

# setwd("D:/Mael/R_projects/ithomiini_current/")
setwd(initial.wd) # Retour au wd de R_codes

save(Response_plot_summary_table, file = paste0("./controls/response_plots/Response_plot_summary_table_",min(index_model_to_compute),"_to_",max(index_model_to_compute),".RData"), version = "2")
saveRDS(Response_plot_summary_table, file = paste0("./controls/response_plots/Response_plot_summary_table_",min(index_model_to_compute),"_to_",max(index_model_to_compute),".rds"), version = "2")
# save(Response_plot_summary_table, file = "./controls/response_plots/Response_plot_summary_table.RData", version = "2")
# saveRDS(Response_plot_summary_table, file = "./controls/response_plots/Response_plot_summary_table.rds", version = "2")



