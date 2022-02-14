##### Functions to compute residual maps from Biodiversity index maps #####

# Author: Maël Doré
# Contact: mael.dore@gmail.com

### Contents ###

# 1/ Map residuals between two raster mapswith indices computed per communities/pixels
# 2/ Map residuals between two raster maps with indices computed per regions

### For all function, can specify the statistial model to use to compute the relationship between the two indices
  # Chose the model type with the "method" argument
    # "LM" = Linear model on every pixels. Can only model linear relationship such as Y = aX + b. Do not account for spatial autocorrelation.
    # "LM_quadratic" = Linear model that includes a quadratic effect such as Y = aX + bX² + c. Do not account for spatial autocorrelation.
    # "GAM" = Generalized Additive Model that provides a smooth fit to the data (default). Do not account for spatial autocorrelation.
    # "GLS" = Generalized Least Square model with spatial covariance structure (e.g., Devictor et al., 2010). Explicitely account for spatial autocorrelation.
       # Can chose the type of spatial autocorrelation structure to fit to the data as in nlme::corSpatial() with the argument "corSpatial".


##### 1/ Function to compute and map residuals between two Biodiversity index maps #####

# Can randomly subsample data to limit spatial autocorrelation issues with the "subsample_size" argument


# Inputs = Two biodiversity index Raster Layers

# Outputs = Raster Layer of residuals of Y ~ X
#           Fitted model (optional: "include_model" = T)

# Prints = Scatterplot of Y ~ X with predicts and R² (Optional: "scatterplot" = T)
#          Variogram and correlogram of GLS models (Optional: "plot_variogram_correlogram" = T)


### To improve
  # Add a number of runs to aggregate to account for the spatial heterogeneity of each subsample


map_residuals <- function (y_index, x_index, 
                           method = "GAM", # Either LM for linear relationship, LM_quadratic to add quadratic effects, or GLS to account for spatial autocorrelation, or GAM to allow non-linear relationship.
                           corSpatial = "corExp", # To choose the type of spatial autocorrelation structure as in nlme::corSpatial()
                           subsample_size = NA,  # Number of sites to subsample
                           subsample_runs = 10,  # Number of random sample runs to aggregate
                           include_model = F, # To include the fitted model in the output
                           scatterplot = T,  # To display the scatter plot with regression line and R²
                           plot_variogram_correlogram = F) # To display the variogram and correlogram in case of spatial GLS
  
{
  y_values <- getValues(y_index)
  x_values <- getValues(x_index)
  indices_df_full <- data.frame(y = y_values, x = x_values)
  indices_df <- indices_df_full[apply(X = indices_df_full, MARGIN = 1, FUN = sum, na.rm = T) > 0, ] # Keep only non-null communities
  
  # See what is in common between methods and move it after, outside if
  
  if (method == "LM") 
  {
    model <- lm(data = indices_df, formula = y ~ x)
    
    # Extract R²
    R2 <- round(summary(model)$adj.r.squared, 3)
    
    # Print output
    cat(paste0("\nLM for ", names(y_index), " ~ ", names(x_index), ":\n\n"))
    summary(model)
  }
  
  if (method == "LM_quadratic") 
  {
    model <- lm(data = indices_df, formula = y ~ x + I(x^2))
    
    # Extract R²
    R2 <- round(summary(model)$adj.r.squared, 3)
    
    # Print output
    cat(paste0("\nLM with quadratic effects for ", names(y_index), " ~ ", names(x_index), ":\n\n"))
    summary(model)
  }
  
  if (method == "GAM")
  {
    model <- mgcv::gam(data = indices_df, formula = y ~ s(x))
    
    # Extract R²
    R2 <- round(summary(model)$r.sq, 3)
    
    # Print output
    cat(paste0("\nGAM for ", names(y_index), " ~ ", names(x_index), ":\n\n"))
    summary(model)
  }
  
  if (method == "GLS")
  {
    # ?nlme::gls
    # ?nlme::corSpatial
    # 
    # ?nlme::corExp # Exponential spatial correlation structure
    # ?nlme::corLin # Linear spatial correlation structure
    # ?nlme::corGaus # Gaussian spatial correlation structure = Sinusoidal shape
    # ?nlme::corSpher # Spherical spatial correlation structure = Exponential structure with distance threshold for plateau
    # ?nlme::corRatio # Rational quadratic spatial correlation structure = Quadratic shape
    
    # Extract coordinates
    coords <- raster::coordinates(y_index)
    
    set.seed(seed = 1)
    
    ### Add a loop to do this N times and compute the mean coefficients to avois issue of random sampling ###
    
    ### Extract the sampled sites
    
    sample_indices <- sample(x = 1:nrow(indices_df), size = subsample_size)
    indices_df <- indices_df[sample_indices, ]
    sampled_coords <- coords[sample_indices, ]
    
    ### Compute correlation matrix depending on the type of Spatial correlation structure chosen
    if(corSpatial == "corExp") { corr_matrix <- nlme::corExp(form = ~ sampled_coords[, 1] + sampled_coords[, 2], nugget = TRUE) }
    if(corSpatial == "corLin") { corr_matrix <- nlme::corLin(form = ~ sampled_coords[, 1] + sampled_coords[, 2], nugget = TRUE) }
    if(corSpatial == "corGaus") { corr_matrix <- nlme::corGaus(form = ~ sampled_coords[, 1] + sampled_coords[, 2], nugget = TRUE) }
    if(corSpatial == "corSpher") { corr_matrix <- nlme::corGaus(form = ~ sampled_coords[, 1] + sampled_coords[, 2], nugget = TRUE) }
    if(corSpatial == "corRatio") { corr_matrix <- nlme::corRatio(form = ~ sampled_coords[, 1] + sampled_coords[, 2], nugget = TRUE) }
    
    ### Compute GLS models
    gls_null_model <- nlme::gls(data = indices_df, model = y ~ x + I(x^2))
    model <- nlme::gls(data = indices_df, model = y ~ x + I(x^2), correlation = corr_matrix)
    
    # Extract residuals for null model and true model
    resids_null <- residuals(gls_null_model)
    resids <- residuals(model)
    
    # Extract R²
    R2 <- round(performance::r2(model)[[1]], 3)
    
    # Print output
    cat(paste0("\nGLS for ", names(y_index), " ~ ", names(x_index), " with spatial autocorrelation structure:\n\n"))
    summary(model)
    
    ### Plot variogram and correlogram to check fit of the model (only one time)
    if (plot_variogram_correlogram)
    {
      ### Find a way to get the variogram on a plot grid! Recreate them from output ? ###
      
      par(mfrow = c(2,2))
      
      # Check the fit of the model to the semi-variogram for Pearson residuals
      gls_variogram <- nlme::Variogram(model, form = ~ sampled_coords[, 1] + sampled_coords[, 2], resType = "pearson")
      plot(gls_variogram, smooth = TRUE, ylim = c(0, 1.2), main = "Semi-variogram for Pearson residuals")
      
      # Check for the absence of trends on the semi-variogram for normalized residuals
      gls_variogram_norm <- nlme::Variogram(model, form = ~ sampled_coords[, 1] + sampled_coords[, 2], resType = "normalized")
      plot(gls_variogram_norm, smooth = TRUE, ylim = c(0, 1.2), main = "Semi-variogram for Normalized residuals")
      
      # Define bins size
      sp_dist <- geosphere::distm(sampled_coords, fun = geosphere::distHaversine)/1000
      bins <- round(x= 1 + 3.3*log(x = max(sp_dist), base = 2), digits = 0) # Struges's rule
      # bins <- 20
      increment <- round(max(sp_dist)/bins,0)
      
      # ?ncf::correlog
      
      # Generate correlogram of residuals for the GLS model
      gls_correlogram <- ncf::correlog(x = sampled_coords[, 1], y = sampled_coords[, 2], z = resids_null, latlon = T, increment = increment, resamp = 100, na.rm = T, quiet = T)
      # Plot correlogram for GLS with spatial autocorrelation structure
      plot(gls_correlogram$mean.of.class[1:10], gls_correlogram$correlation[1:10], type = "b", col = c("black","red")[(gls_correlogram$p > 0.05) + 1], pch = 16, cex = 1, lwd = 1.5, main = "Correlogram for GLS model", xlab = "Distance (km)", ylab = "Moran's I", cex.lab = 1.2, cex.axis = 1) ; abline(h = 0)
      
      # Generate correlogram of residuals for the null model
      gls_correlogram_null <- ncf::correlog(x = sampled_coords[, 1], y = sampled_coords[, 2], z = resids, latlon = T, increment = increment, resamp = 100, na.rm = T, quiet = T)
      # Plot correlogram for GLS without spatial autocorrelation structure
      plot(gls_correlogram_null$mean.of.class[1:10], gls_correlogram_null$correlation[1:10], type = "b", col = c("black","red")[(gls_correlogram_null$p > 0.05) + 1], pch = 16, cex = 1, lwd = 1.5, main = "Correlogram for Null model", xlab = "Distance (km)", ylab = "Moran's I", cex.lab = 1.2, cex.axis = 1) ; abline(h=0)
      
      par(mfrow = c(1,1))
    }
  }
  
  # Plot indices values with model predictions and GOF
  if (scatterplot)
  {
    # indices_df <- indices_df[sample(x = 1:nrow(indices_df), size = 100), ]
    
    # Extract residuals for the color scheme
    resids <- residuals(model)
    indices_df <- indices_df %>% 
      mutate(Residuals = resids)
    
    # Compute predict for the regression line
    predicts <- predict(object = model, newdata = data.frame(x = seq(from = min(x_values, na.rm = T), to = max(x_values, na.rm = T), length.out = 1000)))
    predicts_df <- data.frame(x = seq(from = min(x_values, na.rm = T), to = max(x_values, na.rm = T), length.out = 1000), y = predicts)
    
    # Extract AICc
    AICc <- round(MuMIn::AICc(model), 1)
    
    # Plot
    g <- ggplot(data = indices_df) +
      geom_point(data = indices_df, aes(x = x, y = y, color = Residuals), alpha = 0.7, show.legend = T) +
      geom_line(data = predicts_df, aes(x = x, y = y), size = 1.5, col = "red", linetype = "solid") +
      
      # Set limits to display origin
      xlim(0, max(x_values, na.rm = T)) +
      ylim(0, max(y_values, na.rm = T)) +
      
      # Display R²
      annotate("text", x = 0.12*max(x_values, na.rm = T), y = 0.95*max(y_values, na.rm = T), label = paste0("R² = ", R2), size = 5, fontface = 2) +
      
      # Display AICc
      annotate("text", x = 0.12*max(x_values, na.rm = T), y = 0.88*max(y_values, na.rm = T), label = paste0("AICc = ", AICc), size = 5, fontface = 2) +
      
      # Add legend for predict
      geom_segment(x = 0.00*max(x_values, na.rm = T), y = 0.81*max(y_values, na.rm = T), xend = 0.10*max(x_values, na.rm = T), yend = 0.81*max(y_values, na.rm = T), color = "#FF000080", alpha = 0.7, size = 1.5) +
      annotate("text", x = 0.25*max(x_values, na.rm = T), y = 0.81*max(y_values, na.rm = T), label = paste0(method, " predict"), size = 5, fontface = 2) +
      
      # Set axis labels
      ylab(label = names(y_index)) +
      xlab(label = names(x_index)) +
      
      # Set title
      ggtitle(paste0(names(y_index), " ~ ", names(x_index))) +
      
      # Arrange axes limits
      ggthemes::geom_rangeframe(data = indices_df, aes(x = x, y = y), size = 1.4) +
      
      # Set color gradient
      scale_color_gradient2(low = scales::muted("blue"), high = scales::muted("red")) +
      
      # Adjust 
      theme(panel.background = element_rect(fill = "white"),
            plot.title = element_text(face = 2, size = 18, hjust = 0.5, margin = margin(t = 5, b = 10)),
            
            legend.position = c(0.85, 0.3),
            legend.background = element_rect(fill = NA, colour = NA),
            legend.title = element_text(size = 14, vjust = 3, face = "bold"),
            legend.text = element_text(size = 12, face = "bold"),
            # # margin = margin(t = 0, unit = "pt")),
            # # legend.key.size = unit(1.5, 'lines'),
            # # legend.spacing.y = unit(1,"cm"),
            
            axis.line = element_line(size = 1.2),
            axis.ticks = element_line(size = 1.2),
            axis.ticks.length = unit(10, "pt"),
            axis.text = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            axis.title.x = element_text(margin = margin(t = 10, b = 5)),
            axis.title.y = element_text(margin = margin(l = 5, r = 10)))
    
    print(g)
  }
  
  ### Predict values for all cells
  predicts_values <- predict(object = model, newdata = indices_df_full)
  
  ### Put predicts in a Raster Layer
  predicts_raster <- y_index
  predicts_raster@data@values <- c(predicts_values)
  
  ### Compute residuals
  residuals_raster <- y_index - predicts_raster
  
  ### Export
  if (include_model) 
  { 
    return(list(residuals_raster, model)) 
  } else {
    return(residuals_raster)
  }
}


### 2/ Function to compute residuals within regions ####

# Inputs = Two biodiversity index Raster Layers computed per regions
#          List of SpatialPolygons for each region 

# Outputs = Raster Layer of residuals of Y ~ X
#           Fitted model (optional: "include_model" = T)

# Prints = Scatterplot of Y ~ X with predicts and R² (Optional: "scatterplot" = T)
#          Variogram and correlogram of GLS models (Optional: "plot_variogram_correlogram" = T)


map_residuals_within_regions <- function (y_index, x_index, 
                                          list_sp_regions, # Provide region Spatial Polygon to define borders
                                          method = "GAM", # Either LM for linear relationship, LM_quadratic to add quadratic effects, or GLS to account for spatial autocorrelation, or GAM to allow non-linear relationship.
                                          corSpatial = "corExp", # To choose the type of spatial autocorrelation structure as in nlme::corSpatial() for GLS
                                          include_model = F, # To include the fitted model in the output
                                          scatterplot = T,  # To display the scatter plot with regression line and R²
                                          plot_variogram_correlogram = F) # To display the variogram and correlogram in case of spatial GLS
  
{
  y_values <- sapply(X = list_sp_regions, FUN = raster::extract, x = y_index)
  y_values_regions <- round(unlist(lapply(X = y_values, FUN = median, na.rm = T)),3)
  
  x_values <- sapply(X = list_sp_regions, FUN = raster::extract, x = x_index)
  x_values_regions <- round(unlist(lapply(X = x_values, FUN = median, na.rm = T)),3)
  
  indices_df <- data.frame(y = y_values_regions, x = x_values_regions)
  
  # See what is in common between methods and move it after, outside if
  
  if (method == "LM") 
  {
    model <- lm(data = indices_df, formula = y ~ x)
    
    # Extract R²
    R2 <- round(summary(model)$adj.r.squared, 3)
    
    # Print output
    cat(paste0("\nLM for ", names(y_index), " ~ ", names(x_index), ":\n\n"))
    summary(model)
  }
  
  if (method == "LM_quadratic") 
  {
    model <- lm(data = indices_df, formula = y ~ x + I(x^2))
    
    # Extract R²
    R2 <- round(summary(model)$adj.r.squared, 3)
    
    # Print output
    cat(paste0("\nLM with quadratic effects for ", names(y_index), " ~ ", names(x_index), ":\n\n"))
    summary(model)
  }
  
  if (method == "GAM")
  {
    model <- mgcv::gam(data = indices_df, formula = y ~ s(x))
    
    # Extract R²
    R2 <- round(summary(model)$r.sq, 3)
    
    # Print output
    cat(paste0("\nGAM for ", names(y_index), " ~ ", names(x_index), ":\n\n"))
    summary(model)
  }
  
  if (method == "GLS")
  {
    # ?nlme::gls
    # ?nlme::corSpatial
    # 
    # ?nlme::corExp # Exponential spatial correlation structure
    # ?nlme::corLin # Linear spatial correlation structure
    # ?nlme::corGaus # Gaussian spatial correlation structure = Sinusoidal shape
    # ?nlme::corSpher # Spherical spatial correlation structure = Exponential structure with distance threshold for plateau
    # ?nlme::corRatio # Rational quadratic spatial correlation structure = Quadratic shape
    
    # Extract coordinates of the centroid of regions
    all_sp_polygons <- do.call(raster::bind, list_sp_regions)
    centroids_coords <- rgeos::gCentroid(spgeom = all_sp_polygons, byid=TRUE)@coords
    
    set.seed(seed = 1)
    
    ### Compute correlation matrix depending on the type of Spatial correlation structure chosen
    if(corSpatial == "corExp") { corr_matrix <- nlme::corExp(form = ~ centroids_coords[, 1] + centroids_coords[, 2], nugget = TRUE) }
    if(corSpatial == "corLin") { corr_matrix <- nlme::corLin(form = ~ centroids_coords[, 1] + centroids_coords[, 2], nugget = TRUE) }
    if(corSpatial == "corGaus") { corr_matrix <- nlme::corGaus(form = ~ centroids_coords[, 1] + centroids_coords[, 2], nugget = TRUE) }
    if(corSpatial == "corSpher") { corr_matrix <- nlme::corGaus(form = ~ centroids_coords[, 1] + centroids_coords[, 2], nugget = TRUE) }
    if(corSpatial == "corRatio") { corr_matrix <- nlme::corRatio(form = ~ centroids_coords[, 1] + centroids_coords[, 2], nugget = TRUE) }
    
    ### Compute GLS models
    gls_null_model <- nlme::gls(data = indices_df, model = y ~ x + I(x^2))
    model <- nlme::gls(data = indices_df, model = y ~ x + I(x^2), correlation = corr_matrix)
    
    # Extract residuals for null model and true model
    resids_null <- residuals(gls_null_model)
    resids <- residuals(model)
    
    # Extract R²
    R2 <- round(performance::r2(model)[[1]], 3)
    
    # Print output
    cat(paste0("\nGLS for ", names(y_index), " ~ ", names(x_index), " with spatial autocorrelation structure:\n\n"))
    summary(model)
    
    ### Plot variogram and correlogram to check fit of the model (only one time)
    if (plot_variogram_correlogram)
    {
      ### Find a way to get the variogram on a plot grid! Recreate them from output ? ###
      
      par(mfrow = c(2,2))
      
      # Check the fit of the model to the semi-variogram for Pearson residuals
      gls_variogram <- nlme::Variogram(model, form = ~ centroids_coords[, 1] + centroids_coords[, 2], resType = "pearson")
      plot(gls_variogram, smooth = TRUE, ylim = c(0, 1.2), main = "Semi-variogram for Pearson residuals")
      
      # Check for the absence of trends on the semi-variogram for normalized residuals
      gls_variogram_norm <- nlme::Variogram(model, form = ~ centroids_coords[, 1] + centroids_coords[, 2], resType = "normalized")
      plot(gls_variogram_norm, smooth = TRUE, ylim = c(0, 1.2), main = "Semi-variogram for Normalized residuals")
      
      # Define bins size
      sp_dist <- geosphere::distm(centroids_coords, fun = geosphere::distHaversine)/1000
      bins <- round(x= 1 + 3.3*log(x = max(sp_dist), base = 2), digits = 0) # Struges's rule
      # bins <- 20
      increment <- round(max(sp_dist)/bins,0)
      
      # ?ncf::correlog
      
      # Generate correlogram of residuals for the GLS model
      gls_correlogram <- ncf::correlog(x = centroids_coords[, 1], y = centroids_coords[, 2], z = resids_null, latlon = T, increment = increment, resamp = 100, na.rm = T, quiet = T)
      # Plot correlogram for GLS with spatial autocorrelation structure
      plot(gls_correlogram$mean.of.class[1:10], gls_correlogram$correlation[1:10], type = "b", col = c("black","red")[(gls_correlogram$p > 0.05) + 1], pch = 16, cex = 1, lwd = 1.5, main = "Correlogram for GLS model", xlab = "Distance (km)", ylab = "Moran's I", cex.lab = 1.2, cex.axis = 1) ; abline(h = 0)
      
      # Generate correlogram of residuals for the null model
      gls_correlogram_null <- ncf::correlog(x = centroids_coords[, 1], y = centroids_coords[, 2], z = resids, latlon = T, increment = increment, resamp = 100, na.rm = T, quiet = T)
      # Plot correlogram for GLS without spatial autocorrelation structure
      plot(gls_correlogram_null$mean.of.class[1:10], gls_correlogram_null$correlation[1:10], type = "b", col = c("black","red")[(gls_correlogram_null$p > 0.05) + 1], pch = 16, cex = 1, lwd = 1.5, main = "Correlogram for Null model", xlab = "Distance (km)", ylab = "Moran's I", cex.lab = 1.2, cex.axis = 1) ; abline(h=0)
      
      par(mfrow = c(1,1))
    }
  }
  
  # Plot indices values with model predictions and GOF
  if (scatterplot)
  {
    # indices_df <- indices_df[sample(x = 1:nrow(indices_df), size = 100), ]
    
    # Extract residuals for the color scheme
    resids <- residuals(model)
    indices_df <- indices_df %>% 
      mutate(Residuals = resids)
    
    # Compute predict for the regression line
    predicts <- predict(object = model, newdata = data.frame(x = seq(from = min(x_values_regions, na.rm = T), to = max(x_values_regions, na.rm = T), length.out = 1000)))
    predicts_df <- data.frame(x = seq(from = min(x_values_regions, na.rm = T), to = max(x_values_regions, na.rm = T), length.out = 1000), y = predicts)
    
    # Extract AICc
    AICc <- round(MuMIn::AICc(model), 1)
    
    # Plot
    g <- ggplot(data = indices_df) +
      geom_point(data = indices_df, aes(x = x, y = y, color = Residuals), alpha = 0.7, show.legend = T) +
      geom_line(data = predicts_df, aes(x = x, y = y), size = 1.5, col = "red", linetype = "solid") +
      
      # Set limits to display origin
      xlim(0, max(x_values_regions, na.rm = T)) +
      ylim(0, max(y_values_regions, na.rm = T)) +
      
      # Display R²
      annotate("text", x = 0.12*max(x_values_regions, na.rm = T), y = 0.95*max(y_values_regions, na.rm = T), label = paste0("R² = ", R2), size = 5, fontface = 2) +
      
      # Display AICc
      annotate("text", x = 0.12*max(x_values_regions, na.rm = T), y = 0.88*max(y_values_regions, na.rm = T), label = paste0("AICc = ", AICc), size = 5, fontface = 2) +
      
      # Add legend for predict
      geom_segment(x = 0.00*max(x_values_regions, na.rm = T), y = 0.81*max(y_values_regions, na.rm = T), xend = 0.10*max(x_values_regions, na.rm = T), yend = 0.81*max(y_values_regions, na.rm = T), color = "#FF000080", alpha = 0.7, size = 1.5) +
      annotate("text", x = 0.25*max(x_values_regions, na.rm = T), y = 0.81*max(y_values_regions, na.rm = T), label = paste0(method, " predict"), size = 5, fontface = 2) +
      
      # Set axis labels
      ylab(label = names(y_index)) +
      xlab(label = names(x_index)) +
      
      # Set title
      ggtitle(paste0(names(y_index), " ~ ", names(x_index))) +
      
      # Arrange axes limits
      ggthemes::geom_rangeframe(data = indices_df, aes(x = x, y = y), size = 1.4) +
      
      # Set color gradient
      scale_color_gradient2(low = scales::muted("blue"), high = scales::muted("red")) +
      
      # Adjust 
      theme(panel.background = element_rect(fill = "white"),
            plot.title = element_text(face = 2, size = 18, hjust = 0.5, margin = margin(t = 5, b = 10)),
            
            legend.position = c(0.85, 0.3),
            legend.background = element_rect(fill = NA, colour = NA),
            legend.title = element_text(size = 14, vjust = 3, face = "bold"),
            legend.text = element_text(size = 12, face = "bold"),
            # # margin = margin(t = 0, unit = "pt")),
            # # legend.key.size = unit(1.5, 'lines'),
            # # legend.spacing.y = unit(1,"cm"),
            
            axis.line = element_line(size = 1.2),
            axis.ticks = element_line(size = 1.2),
            axis.ticks.length = unit(10, "pt"),
            axis.text = element_text(size = 14, face = "bold"),
            axis.title = element_text(size = 16, face = "bold"),
            axis.title.x = element_text(margin = margin(t = 10, b = 5)),
            axis.title.y = element_text(margin = margin(l = 5, r = 10)))
    
    print(g)
  }
  
  ### Predict values for all regions
  predicts_values <- predict(object = model, newdata = indices_df)
  
  ### Compute residuals per regions
  residuals_regions <- y_values_regions - predicts_values
  
  ### Turn residuals into a Raster Layer
  
  # Generate a mask for terrestrial areas to use as background
  continental_mask <- (y_index < Inf) - 1
  
  # Initiate stack per region
  region_stack <- stack()
  # Loop per region
  for (i in 1:length(list_sp_regions))
  {
    # i <- 1
    
    region_stack <- addLayer(region_stack, rasterize(x = list_sp_regions[[i]], 
                                                     y = continental_mask, # Provide the grid to fill with CRS, bbox and resolution
                                                     field = residuals_regions[i], # How to fill non empty cells. With the value of a variable in the df of the sp_obj, or directly with a fixed value ?
                                                     background = NA)) # Value to use to fill empty cells)
  }
  # plot(region_stack) 
  
  # Aggregate all regions in one layer
  all_regions_raster <- calc(x = region_stack, fun = median, na.rm = T)
  # Add null values for continental borders
  residuals_raster <- continental_mask
  residuals_raster@data@values[!is.na(all_regions_raster[])] <- all_regions_raster[!is.na(all_regions_raster[])]
  
  # plot(all_regions_raster)
  # plot(residuals_raster)
  
  ### Export
  if (include_model) 
  { 
    return(list(residuals_raster, model)) 
  } else {
    return(residuals_raster)
  }
}

