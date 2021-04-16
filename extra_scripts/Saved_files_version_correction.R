##### Script to correct saving version of objects in the entire Project folder #####

rm(list = ls())

##### Get a list of path to access saved files #####

all_save_files <- list.files(pattern = "(\\.RData)$|(\\.rds)$", ignore.case = T, all.files = T, recursive = T, full.names = T)
Rdata_files <- list.files(pattern = "(\\.RData)$", ignore.case = T, all.files = T, recursive = T, full.names = T)
rds_files <- list.files(pattern = "(\\.rds)$", ignore.case = T, all.files = T, recursive = T, full.names = T)

# Remove files from the project library to save time and avoid errors
Rdata_files <- Rdata_files[!stringr::str_detect(Rdata_files, pattern = "/packrat/")]
rds_files <- rds_files[!stringr::str_detect(rds_files, pattern = "/packrat/")]

##### Convert from version 3 to version 2 #####

### For .Rdata ###

for (custom_unique_index in seq_along(Rdata_files)) { # Need to use a "custom_unique_index" to avoid error when the loaded object name is "i" !
  
  before <- ls() # Store env. objects list before loading new .Rdata
  load(file = Rdata_files[custom_unique_index]) # Load the object (no need to specify the version)
  after <- ls()  # Store env. objects list after loading new .Rdata
  loaded <- setdiff(after, before)  # Find loaded objects names
  loaded <- loaded[loaded != "before"]
  
  save(list = loaded, file = Rdata_files[custom_unique_index], version = "2") # Save under version 2
  rm(list = loaded)     # Remove from env.
  
  print(custom_unique_index)
}




### For .rds ###

for (custom_unique_index in seq_along(rds_files)) { # Need to use a "custom_unique_index" to avoid error when the loaded object name is "i" !
  
  file_to_save <- readRDS(file = rds_files[custom_unique_index]) # Load the object 
  
  saveRDS(file_to_save, file = rds_files[custom_unique_index], version = "2") # Save under version 2
  rm(file_to_save)     # Remove from env.
  
  print(custom_unique_index)
}


##### Convert from version 2 to version 3 #####

### For .Rdata ###

for (i in seq_along(Rdata_files)) {
  
  before <- ls() # Store env. objects list before loading new .Rdata
  load(file = Rdata_files[i]) # Load the object (no need to specify the version)
  after <- ls()  # Store env. objects list after loading new .Rdata
  loaded <- setdiff(after, before)  # Find loaded objects names
  loaded <- loaded[loaded != "before"]
  
  save(list = loaded, file = Rdata_files[i]) # Save under version 3 (the default one for R 3.6.2)
  rm(list = loaded)     # Remove from env.
}

### For .rds ###

for (i in seq_along(rds_files)) {
  
  file_to_save <- readRDS(file = rds_files[i]) # Load the object (no need to specify the version)
  
  saveRDS(file_to_save, file = rds_files[i]) # Save under version 3 (the default one for R 3.6.2)
  rm(file_to_save)     # Remove from env.
}



