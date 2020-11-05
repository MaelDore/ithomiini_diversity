
####### Script 1: Prepare environmental stacks #######

# Import environmental data, set extent and resolution of study area, select variables, and stack them

# 1/ Import climatic data
# 2/ Select climatic variables based on correlation, clustering, and biological signifiance
# 3/ Import elevation data
# 4/ Import forest cover data
# 5/ Synchronize NA and stack all environnemental layers 
# 6/ Check for correlation in the final stack

### Inputs:

   # Climatic data from MERRAclim (Vega et al., 2017)
   # Elevation data from SRTM (Farr et al., 2007)
   # Forest cover data from Landsat Vegetation Continuous Fields (Sexton et al., 2013)

### Outputs:

   # Final environmental stack
   # Climatic variable selection assessment: HAC dendrogram, PCA correlation circle, RDA results
   # Correlation matrix for final stack

#####



# Remove environment
rm(list = ls())

# Load libraries
library(raster)


###### 1/ Import Climatic data ######

# 1.1/ Resolution = 5min ####

# Load a Bioclim layer to use as mask to define terrestrial boundaries
BioClim_mask <- stack(raster("./input_data/Bioclim/bio_5m_esri/bio_1"))

# plot(BioClim_mask)

# Load MERRA Clim stack
climate.vars <- paste0("bio", 1:19)
envData <- NA
for (j in climate.vars) {
  
  if(j == climate.vars[1]) { # Create the stack
    
    envData <- stack(raster(paste0("./input_data/MERRA_Clim/5m_mean_00s/5m_mean_00s_", j,".tif")))
  
  } else { # Add layer iteratively
     
    envData <- addLayer(envData, raster(paste0("./input_data/MERRA_Clim/5m_mean_00s/5m_mean_00s_", j,".tif")))
    
  }
}
names(envData) <- paste0("bio", 1:19) # Create layer name

envData <- resample(x = envData, y = BioClim_mask, method = "bilinear") # Ajust extent and resolution to match the one of Bioclim
Merra_clim_5 <- mask(envData, BioClim_mask) # Apply Bioclim mask to keep only terrestrial data

# x11()
# plot(Merra_clim_5)
# plot(Merra_clim_5[[1]])

save(Merra_clim_5, file = "./input_data/MERRA_Clim/MERRAClim_5m_full.RData")

# Crop to Neotropics

xmin = -120 ; xmax = -30 ; ymin = -37 ; ymax = 28
e <- extent(c(xmin,xmax,ymin,ymax))

Merra_clim_5_neotrop <- crop(Merra_clim_5, e)

# x11()
# plot(Merra_clim_5_neotrop[[1]])

save(Merra_clim_5_neotrop, file = "./input_data/MERRA_Clim/MERRAClim_5m_cropped.RData")
saveRDS(Merra_clim_5_neotrop, file = "./input_data/MERRA_Clim/MERRAClim_5m_cropped.rds")

# 1.2/ Downgrade resolution to final modeling resolution = 15min ####

Merra_clim_15_neotrop <- aggregate(x = Merra_clim_5_neotrop, fact=3, fun=mean, expand=TRUE, na.rm = F)

# x11()
# plot(Merra_clim_15_neotrop[[1]])

save(Merra_clim_15_neotrop, file = "./input_data/MERRA_Clim/MERRAClim_15m_cropped.RData")
saveRDS(Merra_clim_15_neotrop, file = "./input_data/MERRA_Clim/MERRAClim_15m_cropped.rds")



##### 2/ Climatic variable selection #####

# Effacer l'environnement
rm(list = ls())

# Choose resolution for the analysis
res <- "5"
res <- "15"

clim_stack <- readRDS(file = paste0("./input_data/MERRA_Clim/MERRAClim_", res,"m_cropped.rds"))


### 2.1/ HAC clustering on environmental variables ####

# library(virtualspecies)
# removeCollinearity(raster.stack = clim_stack, method = "spearman", multicollinearity.cutoff = 0.7, sample.points = T, nb.points = 10000, plot = T)

source(file = "./functions/raster_multicollinearity.R")

pdf(file = paste0("./controls/variable_selection/HAC_clim_", res ,".pdf"), width = 8, height = 6)

# raster_multicollinearity(raster.stack = clim_stack, method = "spearman", multicollinearity.cutoff = 0.7, sample.points = T, nb.points = 10000, plot = T)
raster_multicollinearity(raster.stack = clim_stack, method = "pearson", multicollinearity.cutoff = 0.7, sample.points = T, nb.points = 10000, plot = T)

dev.off()

### 2.2/ PCA and RDA on environmental variables ####

library(ade4)
library(vegan)

# Need to subsample for only 10000 pixels
data_clim <- data.frame(rep(NA,10000))

index.sample <- sample(x = which(!is.na(clim_stack[[1]]@data@values)), size = 10000) # Get index of cells which are not NA

for (i in 1:(length(names(clim_stack)))) {
  data_clim[,i] <- clim_stack[[i]]@data@values[index.sample] # Retrieve 10000 values from each layer
}
names(data_clim) <- names(clim_stack)

### PCA on all climatic variables

PCA_full <-  dudi.pca(df = data_clim, center = T, nf = 3, scannf = F) # Generate the PCA

# Compute % variance explained by each PC axis
var_perc <- round(PCA_full$eig/sum(PCA_full$eig)*100, 1) 
cum_var_perc <- round(cumsum(var_perc), 1) # Cumulated explained variance
PCA_full$co # Correlation: variables ~ PC axis

# Plot correlation circle of PCA with results of RDA

pdf(file = paste0("./controls/variable_selection/PCA_All_clim_var_", res, ".pdf"), width = 8, height = 6)

original_ext_margins <- par()$oma
original_int_margins <- par()$mar
par(oma = c(0,2,0,3), mar = c(5.1,7,4.1,7), xpd=NA)
s.corcircle(PCA_full$co, xax=1, yax=2, clabel = 1.5)
text(x = c(-1, 0), labels = paste0("PC1: ", var_perc[1], "%"), font = 2, cex = 1.3)
text(x = -1.2 , y = -0.6, adj = c(0, 0), labels = paste0("PC2: ", var_perc[2], "%\nSum: ", cum_var_perc[2], "%"), font = 2, cex = 1.3, srt = 90)
par(mar = original_int_margins, oma = original_ext_margins, xpd = F)

dev.off()

### Compute RDA for 5 variables (BIO1, BIO4, BIO12, BIO15) + BIO3 = Isothermality, to test representativity

vars_index <- c(1,3,4,12,15) # Indices of climatic variable to keep

Y5 <- data_clim[,vars_index] # Explanatory variables = set of selected variables
X5 <- data_clim[,-vars_index] # Responses variables = other climatic variables

RDA5 <- rda(X = X5, Y = Y5 , scale = T) 
res_RDA <- round(summary(RDA5)$cont[[1]][2:3,]*100, 1) ; res_RDA # To see proportions of constrained and unconstrained inertia
r2_RDA <- round(RsquareAdj(RDA5)$adj.r.squared, 3) ; r2_RDA

# Test significiance of the RDA with permutation on the pseudo F-ratio
F_test_RDA <- anova.cca(object = RDA5, permutations = 1000) ; F_test_RDA 
F_value <- round(F_test_RDA$F[1], 0)
p_value <- round(F_test_RDA$`Pr(>F)`[1], 3)

### PCA on 5 selected bioclim var.
PCA_5vars <-  dudi.pca(df = data_clim[,c(1,3,4,12,15)], center = T, nf = 3, scannf = F) # Generate the PCA

# Compute % variance explained by each PC axis.
var_perc <- round(PCA_5vars$eig/sum(PCA_5vars$eig)*100, 1) 
cum_var_perc <- round(cumsum(var_perc),1) # Cumulated explained variance
PCA_5vars$co # Correlation: variables ~ PC axis

pdf(file = paste0("./controls/variable_selection/PCA_Five_clim_var_", res, ".pdf"), width = 8, height = 6)

original_int_margins <- par()$mar
original_ext_margins <- par()$oma
par(oma = c(0,2,2,3), mar = c(5.1,7,4.1,7), xpd=NA)
s.corcircle(PCA_5vars$co, xax=1, yax=2, clabel = 1.5)
text(x = 1.2 , y = -0.8, labels = paste0("PC1: ", var_perc[1], "%"), font = 2, cex = 1.3)
text(x = -1.4 , y = -0.6, font = 2, cex = 1.3, srt = 90,
     labels = paste0("PC2: ", var_perc[2], "%\nSum: ", cum_var_perc[2], "%"))
text(x = 0 , y = 1.1, font = 2, cex = 1.3,
     labels = paste0("RDA: ", res_RDA[2, length(vars_index)], "%, F = ", F_value, ", p = ", p_value))
par(mar = original_int_margins, oma = original_ext_margins, xpd = F)

dev.off()


### Compute RDA for 4 variables (BIO1, BIO4, BIO12, BIO15), to test representativity

vars_index <- c(1,4,12,15) # Indices of climatic variable to keep

Y4 <- data_clim[,vars_index] # Explanatory variables = set of selected variables
X4 <- data_clim[,-vars_index] # Responses variables = other climatic variables

RDA4 <- rda(X = X4, Y = Y4 , scale = T) 
res_RDA <- round(summary(RDA4)$cont[[1]][2:3,]*100, 1) ; res_RDA # To see proportions of constrained and unconstrained inertia
r2_RDA <- round(RsquareAdj(RDA4)$adj.r.squared, 3) ; r2_RDA

# Test significiance of the RDA with permutation on the pseudo F-ratio
F_test_RDA <- anova.cca(object = RDA4, permutations = 1000) ; F_test_RDA 
F_value <- round(F_test_RDA$F[1], 0)
p_value <- round(F_test_RDA$`Pr(>F)`[1], 3)

### PCA on 4 selected bioclim var.
PCA_4vars <-  dudi.pca(df = data_clim[,c(1,4,12,15)], center = T, nf = 3, scannf = F) # Generate the PCA

# Compute % variance explained by each PC axis.
var_perc <- round(PCA_4vars$eig/sum(PCA_4vars$eig)*100, 1) 
cum_var_perc <- round(cumsum(var_perc),1) # Cumulated explained variance
PCA_4vars$co # Correlation: variables ~ PC axis

pdf(file = paste0("./controls/variable_selection/PCA_Four_clim_var_", res, ".pdf"), width = 8, height = 6)

original_int_margins <- par()$mar
original_ext_margins <- par()$oma
par(oma = c(0,2,2,3), mar = c(5.1,7,4.1,7), xpd=NA)
s.corcircle(PCA_4vars$co, xax=1, yax=2, clabel = 1.5)
text(x = 1.2 , y = -0.8, labels = paste0("PC1: ", var_perc[1], "%"), font = 2, cex = 1.3)
text(x = -1.4 , y = -0.6, font = 2, cex = 1.3, srt = 90,
     labels = paste0("PC2: ", var_perc[2], "%\nSum: ", cum_var_perc[2], "%"))
text(x = 0 , y = 1.1, font = 2, cex = 1.3,
     labels = paste0("RDA: ", res_RDA[2, length(vars_index)], "%, F = ", F_value, ", p = ", p_value))
par(mar = original_int_margins, oma = original_ext_margins, xpd = F)

dev.off()




##### 3/ Import Elevation data #####

# Chose resolution
res <- "5"
res <- "15"

# Original script to import SRTM in Google Earth Engine already cropped for Neotropics
Elevation <- raster(x = paste0("./input_data/SRTM/SRTM_5min.tif"))

plot(Elevation)  # Already cropped for Neotropics
Elevation@extent # Not exactly the same than the env. layers

# Define proper extend
xmin = -120 ; xmax = -30 ; ymin = -37 ; ymax = 28
e <- extent(c(xmin,xmax,ymin,ymax))

## Remove Ocean pixels

# Load BioClim layer to apply as mask
BioClim_mask <- stack(raster("./input_data/Bioclim/bio_5m_esri/bio_1"))

# Crop properly the mask
BioClim_mask <- crop(BioClim_mask, e)
plot(BioClim_mask)

Elevation <- resample(x = Elevation, y = BioClim_mask, method = "bilinear")  # Ajust extent and resolution to match the one of Bioclim
Elevation <- mask(Elevation, BioClim_mask) # Apply Bioclim mask to keep only terrestrial data

names(Elevation) <- "Elevation"

# Aggregate for 15min res
if (res == "15") {
  Elevation <- aggregate(x = Elevation, fact=3, fun=mean, expand=TRUE, na.rm = F)
}

plot(Elevation)

# Save results
save(Elevation, file = paste0("./input_data/SRTM/Elevation_", res, ".RData"))
saveRDS(Elevation, file = paste0("./input_data/SRTM/Elevation_", res, ".rds"))





##### 4/ Import Forest cover data #####

# Chose resolution
res <- "5"
res <- "15"


# Original script to import VCF data in Google Earth Engine already cropped for Neotropics

Forests <- raster(x = "./input_data/VCF/VCF_5min.tif")

plot(Forests)  # Already cropped for Neotropics
Forests@extent # Not exactly the same than the env. layers

# Define proper extend
xmin = -120 ; xmax = -30 ; ymin = -37 ; ymax = 28
e <- extent(c(xmin,xmax,ymin,ymax))

## Remove Ocean pixels

# Load BioClim layer to apply as mask
BioClim_mask <- stack(raster("./input_data/Bioclim/bio_5m_esri/bio_1"))

# Crop properly the mask
BioClim_mask <- crop(BioClim_mask, e)
plot(BioClim_mask)

Forests <- resample(x = Forests, y = BioClim_mask, method = "bilinear")  # Ajust extent and resolution to match the one of Bioclim
Forests <- mask(Forests, BioClim_mask) # Apply Bioclim mask to keep only terrestrial data

names(Forests) <- "Forests"

# Aggregate for 15min res
if (res == "15") {
  Forests <- aggregate(x = Forests, fact=3, fun=mean, expand=TRUE, na.rm = F)
}

plot(Forests)

# Save results
save(Forests, file = paste0("./input_data/VCF/Forests_", res, ".RData"))
saveRDS(Forests, file = paste0("./input_data/VCF/Forests_", res, ".rds"))





##### 5/ Final Environmental stack #####

# Chose resolution
res <- "5"
res <- "15"

# Load stacks and layers
Clim_stack <- readRDS(file = paste0("./input_data/MERRA_Clim/MERRAClim_", res, "m_cropped.rds"))
Elevation <- readRDS(file = paste0("./input_data/SRTM/Elevation_", res, ".rds"))
Forests <- readRDS(file = paste0("./input_data/VCF/Forests_", res, ".rds"))

library(raster)
Full_env <- addLayer(Clim_stack, Elevation, Forests)

# Check synchronization of layers data
temp <- Full_env[[1]] ; temp@data@values <- 0 # Create an empty raster
temp@data@values[is.na(Full_env[["Forests"]]@data@values) & !is.na(Full_env[[1]]@data@values)] <- 1 # Show pixels with climatic data but no Forests data
plot(temp) # Need to synchronize coastlines!

# Synchronization of all layers

Full_env_synchro <- virtualspecies::synchroniseNA(Full_env)
plot(Full_env_synchro)

save(Full_env_synchro, file = paste0("./input_data/Env_data/Full_env_synchro_", res, ".RData"))
saveRDS(Full_env_synchro, file = paste0("./input_data/Env_data/Full_env_synchro_", res, ".rds"))

# Synchronization of NA using only selected variables

names(Full_env)

Select_env <- Full_env[[c(paste0("bio", c(1, 4, 12, 15)), "Elevation", "Forests")]]
Select_env <- virtualspecies::synchroniseNA(Select_env)
plot(Select_env)

save(Select_env, file = paste0("./input_data/Env_data/Select_env_", res, ".RData"))
saveRDS(Select_env, file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))





##### 6/ Check for correlation in the final environmental stack #####

# Chose resolution
res <- "5"
res <- "15"

envData <- readRDS(Select_env, file = paste0("./input_data/Env_data/Select_env_", res, ".rds"))

# Need to subsample for only 10000 pixels
data_full <- data.frame(rep(NA,10000))

index.sample <- sample(x = which(!is.na(envData[[1]]@data@values)), size = 10000) # Get index of cells which are not NA

for (i in 1:(length(names(envData)))) {
  data_full[,i] <- envData[[i]]@data@values[index.sample] # Retrieve 10000 values from each layer
}
names(data_full) <- names(envData)

# Check correlation, including Elevation and Forests
cor(data_full)

# Rarity::corPlot(df = data_full, method = "spearman")

library(corrgram)

pdf(file = paste0("./controls/variable_selection/Cor_controls_", res ,".pdf"), width = 6, height = 6)

corrgram::corrgram(x = data_full, type = "data", order = T, method = "spearman",
                   lower.panel = panel.pie,  upper.panel = panel.conf, diag.panel = panel.density,
                   col.regions = colorRampPalette(c("navy", "royalblue", "white", "salmon", "red"))) # pour plotter les correlations entre variables de manière stylée

dev.off()

# Check VIF, including Elevation and Forests
usdm::vif(data_full)
