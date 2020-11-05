# Get mean environmental values per OMU and store it in list.models

load(file = "./input_data/list.models.RData")

list.models$Forests <- list.models$Elevation <- list.models$bio15 <- list.models$bio12 <- list.models$bio4 <- list.models$bio1 <- NA

for (i in 1:nrow(list.models))  # Per unit/OMU
{
  # i <- 7
  
  unit <- list.models$Tag.model[i] # Get OMU name
  
  load(file = paste0("./input_data/Species_data/15/Env_Stacks/Env_stack_",unit,".RData"))       # Load raster of environmental data
  load(file = paste0("./input_data/Species_data/15/Spatial_Points_Objects/occurrences_", unit,".RData"))  # Load sp file of occurrences points
  
  unit_env_data <- extract(x = unit.env, y = unit.points)
  
  unit_mean_env <- apply(X = unit_env_data, MARGIN = 2, FUN = mean)
  
  list.models[i, names(unit_mean_env)] <- unit_mean_env
  
  if (i %% 100 == 0) {print(i)}
}

save(list.models, file = "./input_data/list.models.RData")
saveRDS(list.models, file = "./input_data/list.models.rds")
