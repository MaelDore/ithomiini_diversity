### Scripts for several methods to detect outliers in species distribution data ###

# Clean environment
rm(list = ls())

### Load data ###
load(file = "./input_data/Ithomiini_final.RData")
load(file = "./input_data/Env_data/Select_env_15.RData")

# Extract env data for all ithomiines
quick_data <- Ithomiini_final[complete.cases(raster::extract(x = Select_env, y = Ithomiini_final[, c("Longitude", "Latitude")])), c("Longitude", "Latitude")]

quick_data_env <- raster::extract(x = Select_env, y = quick_data)

# Summary df for outlierness values
outlier_summary_outlierness <- data.frame(OutlierMahdist = rep(NA, nrow(quick_data)))

# Summary df for outlier status (binary as TRUE or FALSE)
outlier_summary_binary <- data.frame(OutlierMahdist = rep(NA, nrow(quick_data)))

### Method 1: Robust Mahalanobis distance

# Robust because estimators for the centroid and the covariance matrix are not sensitive to clusters of outliers 
# because they are estimated on a subsample of the data within the minimum ellipsoid (typically 50%) 
# Take only the points within the center of the data to evaluate the centroid and the shape of the data cloud. 
# Significance evaluated against a Khi² distribution. 

library(rrcovHD)

?rrcovHD::OutlierMahdist
mahalanaobis_test_OutlierMahdist <- rrcovHD::OutlierMahdist(quick_data_env) 

outlier_summary_outlierness$OutlierMahdist <- getDistance(mahalanaobis_test_OutlierMahdist)
# getOutliers(mahalanaobis_test_OutlierMahdist)
outlier_summary_binary$OutlierMahdist <- !as.logical(mahalanaobis_test_OutlierMahdist@flag)

### Method 2: Robust Mahalanobis + Classic Mahalanobis

library(mvoutlier)

?mvoutlier::dd.plot
mvoutlier::dd.plot()

data(humus)
mahalanaobis_test_dd.plot <- dd.plot(log(humus[,c("Co","Cu","Ni")]))

color.plot(log(humus[,c("Co","Cu")]))

plot(x = mahalanaobis_test_dd.plot$md.cla, y = mahalanaobis_test_dd.plot$md.rob, pch = c(16, 1)[mahalanaobis_test_dd.plot$outliers+1])

mahalanaobis_test_dd.plot <- mvoutlier::dd.plot(quick_data_env)
plot(x = mahalanaobis_test_dd.plot$md.cla, y = mahalanaobis_test_dd.plot$md.rob, pch = 16, col = c("limegreen", "red")[mahalanaobis_test_dd.plot$outliers+1])

# Extract robust Mahalanobis distances as outlierness
outlier_summary_outlierness$dd.plot <- mahalanaobis_test_dd.plot$md.rob

# Extract outliers
mahalanaobis_test_dd.plot$outliers # Based on robust Mahalanobis + classic Mahalanobis test against Khi² distribution ???
outlier_summary_binary$dd.plot <- mahalanaobis_test_dd.plot$outliers

plot(x = mahalanaobis_test_dd.plot$md.cla, y = mahalanaobis_test_dd.plot$md.rob, pch = 16, col = c("limegreen", "red")[outlier_summary_binary$OutlierMahdist+1])

### Methods based on PCA

# First apply a PCA, then compute distances in reduced space with specific weighting scheme to highlight dimensions where the outliers stick out

?rrcovHD::OutlierPCDist
?rrcovHD::OutlierPCOut
?rrcovHD::OutlierSign1
?rrcovHD::OutlierSign2

?rrcovHD::getOutliers
?rrcovHD::getDistance

### PCOut
mahalanaobis_test_OutlierPCOut <- rrcovHD::OutlierPCOut(quick_data_env)

outlier_summary_outlierness$OutlierPCOut <- getDistance(mahalanaobis_test_OutlierPCOut)
# getOutliers(mahalanaobis_test_OutlierPCOut)
outlier_summary_binary$OutlierPCOut <- !as.logical(mahalanaobis_test_OutlierPCOut@flag)

### Sign2
mahalanaobis_test_OutlierSign2 <- rrcovHD::OutlierSign2(quick_data_env)

outlier_summary_outlierness$OutlierSign2 <- getDistance(mahalanaobis_test_OutlierSign2)
# getOutliers(mahalanaobis_test_OutlierSign2)
outlier_summary_binary$OutlierSign2 <- !as.logical(mahalanaobis_test_OutlierSign2@flag)


### Method based on local point density = Local Outlier Factor

library(Rlof)

?Rlof::lof

# Applied in geographical space but can also be applied in environmental space
mahalanaobis_test_lof <- Rlof::lof(data = quick_data, k = 15) # Too big to work on all data... so subsample...

# Need to decide a nb of neighboring points... k = 5

sample_quick_data <- quick_data[sample(x = 1:nrow(quick_data), size = 100), ]
mahalanaobis_test_lof <- Rlof::lof(data = sample_quick_data, k = 5)
mahalanaobis_test_lof <- Rlof::lof(data = sample_quick_data, k = 5, cores = 5)

# Plot results
plot(sample_quick_data, cex = mahalanaobis_test_lof)

# Need to decide a threshold...
hist(mahalanaobis_test_lof)
threshold <- 2

mahalanaobis_test_lof_outliers <- mahalanaobis_test_lof > threshold

### Methods based on classification trees to detect the originality of a point as it not associated with other points

## Isotree based on isolation

# Compute how much random cut are needed to isolate an observation in a classification tree within a RF (random subset of variables and observations for each tree) 

library(isotree)
?isotree::isolation.forest()

mahalanaobis_test_isotree <- isotree::isolation.forest(df = quick_data_env, output_score = T) # See many parameters

outlier_summary_outlierness$isotree <- mahalanaobis_test_isotree$scores

# Need to decide a threshold...
hist(mahalanaobis_test_isotree$scores)
threshold <- quantile(x = mahalanaobis_test_isotree$scores, prob = 0.95) # Quantile = 95%

mahalanaobis_test_isotree_outliers <- mahalanaobis_test_isotree$scores > threshold
outlier_summary_binary$isotree <- mahalanaobis_test_isotree_outliers

## Outlier_RF based on proximity

# Compute how often an observation fall with others in the terminal leaves, as pairwise proximity.
# Then find outliers as observations with the least proximity towards other points => compute standardized outlierness as n / sum(squared proximity), then normalized with median and MAD.

library(randomForest)

?randomForest::randomForest()
?randomForest::outlier.randomForest()

# Use integer instead of double to save memory usage?
quick_data_env <- apply(X = apply(X = quick_data_env, MARGIN = 2, FUN = round, digits = 1)*10, MARGIN = 2, FUN = as.integer)

# Subset to try on smaller sample
sample_quick_data_env <- quick_data_env[sample(x = 1:nrow(quick_data_env), size = 20000), ] # Crash at 25000 observations...

# Run unsupervised classification RF while computing pairwise proximity between observations
mahalanaobis_test_RF_proximity <- randomForest::randomForest(x = sample_quick_data_env, proximity = T, do.trace = 10)
mahalanaobis_test_RF_proximity <- randomForest::randomForest(x = sample_quick_data_env, proximity = T, do.trace = T)

# Run unsupervised classification RF while computing pairwise proximity between observations
mahalanaobis_test_RF_proximity <- randomForest::randomForest(x = quick_data_env, proximity = T, do.trace = 10)

mahalanaobis_test_RF_outlierness <- randomForest::outlier(x = mahalanaobis_test_RF_proximity)

outlier_summary_outlierness$Outlier_RF <- mahalanaobis_test_RF_outlierness

# Need to decide a threshold...
hist(mahalanaobis_test_RF_outlierness)
threshold <- quantile(x = mahalanaobis_test_RF_outlierness, prob = 0.95) # Quantile = 95%

mahalanaobis_test_RF_outliers <- mahalanaobis_test_RF_outlierness > threshold
outlier_summary_binary$Outlier_RF <- mahalanaobis_test_RF_outliers


### Save and compare all ###

saveRDS(outlier_summary_outlierness, file = "./outputs/outlier_summary_outlierness.rds")
saveRDS(outlier_summary_binary, file = "./outputs/outlier_summary_binary.rds")

outlier_summary_outlierness <- readRDS(file = "./outputs/outlier_summary_outlierness.rds")
outlier_summary_binary <- readRDS(file = "./outputs/outlier_summary_binary.rds")

cor(outlier_summary_outlierness)
cor(outlier_summary_binary)
