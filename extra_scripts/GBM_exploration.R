### Test GBM ###

library(gbm)

# Load environmental data
envData <- readRDS(file = paste0("./input_data/Env_data/Select_env_15.rds"))

# Get indices of cells with data
available_cells <- which(!apply(X = envData@data@values, MARGIN = 1, FUN = anyNA))

# Extract data
set.seed(seed = 648654)
subset <- sample(x = available_cells, size = 100)
xy <- xyFromCell(object = envData, cell = subset)
env_data_for_model <- raster::extract(envData, SpatialPoints(xy))

PA_data <- (rank(subset) <= 50)*1

# Create random PA data
PA_data <- sample(c(rep(x = 0, times = 50), rep(x = 1, times = 50)))

# GBM modeling

test <- gbm(formula = PA_data ~ ., data = data.frame(cbind(PA_data, env_data_for_model)), 
            distribution = "bernoulli", 
            n.trees = 5000, interaction.depth = 5, 
            n.minobsinnode = 5, shrinkage = 0.001, 
            bag.fraction = 0.5, train.fraction = 1, 
            cv.folds = 3, verbose = F)

str(test)

# Check for overfit
which.min(test$cv.error) # Retrieve the iteration that minimizes the error
plot(1:5000, test$train.error, pch = 16, ylab = "Loss function\nBernoulli deviance", xlab = "Iteration", col = "dodgerblue")
points(1:5000, test$cv.error, pch = 16, col = "limegreen")
abline(v = which.min(test$cv.error), lwd = 3, lty = 2, col = "red")
gbm.perf(test, method="cv", overlay = T)

# Check tree for specific iteration
pretty.gbm.tree(object = test, i.tree = 1)
pretty.gbm.tree(object = test, i.tree = which.min(test$cv.error))


# Check final predicts from CV-folds algo and global algo
gbm$fit # fitted values on the scale of the predictors at the last iteration (log_odds for Bernoulli, log for poisson) ; need to be back-transform (invlogit for Bernoulli, exp for Poisson)
gbm$cv.fitted # fitted values on the scale of the predictors at the last iteration when predicted on the CV-models calibrated with data from other CV-folds (used to compute de cv.errors at each iterations)

arm::invlogit(test$fit) # Fitted values in log-odds scale for bernoulli. Need to apply invlogit to retrieve probabilities
plot(test$cv.fitted, test$fit)
plot(arm::invlogit(test$cv.fitted), arm::invlogit(test$fit))


# If you want to retrieve the predicts for the iteration that minimizes the cv.error, you need to use gbm.pred !
?predict.gbm()
pred_last <- predict.gbm(object = test, newdata = data.frame(env_data_for_model), n.trees = 5000, type = "link")
pred_cv <- predict.gbm(object = test, newdata = data.frame(env_data_for_model), n.trees = which.min(test$cv.error), type = "link")


plot(pred_last, test$fit) # gbm$fit = predictions for the last iteration
plot(pred_last, pred_cv) # Predictions for the last iteration vs. predictions from the iteration that minimizes the cv.error
plot(pred_cv, test$cv.fitted) # gbm$cv.fitted = predictions for the last iteration using models calibrated on data from other CV-fold
                              # They are not the prediction from the proper iteration that minimizes the cv.error !

