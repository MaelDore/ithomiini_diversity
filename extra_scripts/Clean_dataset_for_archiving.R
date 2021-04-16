##### Generate a clean dataset for Zenodo archiving ######

library(tidyverse)

load(file = "./input_data/Ithomiini_final.RData")

# Extract list of curators for acknowledgments
table(Ithomiini_final$FirstOfregistro_contacto)[order(table(Ithomiini_final$FirstOfregistro_contacto), decreasing = T)]

### Remove useless columns
Ithomiini_Zenodo_archive <- Ithomiini_final %>% 
  select(ID_obs, Genus, Species, Sub.species, "Latitude" = Latitude_raster, "Longitude" = Longitude_raster, M.mimicry, F.mimicry, country, Country)

### Clean the country variable

# Merge in one variable, and fill NA from Andre's dataset with 'Brazil'
for (i in seq_along(Ithomiini_Zenodo_archive$ID_obs)) 
{
  if(is.na(Ithomiini_Zenodo_archive$country[i]))
  {
    Ithomiini_Zenodo_archive$final_country[i] <- as.character(Ithomiini_Zenodo_archive$Country[i])
  }
  if(is.na(Ithomiini_Zenodo_archive$Country[i]))
  {
    Ithomiini_Zenodo_archive$final_country[i] <- as.character(Ithomiini_Zenodo_archive$country[i])
  }
  if(is.na(Ithomiini_Zenodo_archive$final_country[i]))
  {
    Ithomiini_Zenodo_archive$final_country[i] <- "Brazil"
  }
}

# Convert to factor
Ithomiini_Zenodo_archive$final_country <- as.factor(Ithomiini_Zenodo_archive$final_country)
table(Ithomiini_Zenodo_archive$final_country) # Check errors

# Correct errors manually
Ithomiini_Zenodo_archive$final_country[which(Ithomiini_Zenodo_archive$final_country == "PERU")] <- "Peru"
Ithomiini_Zenodo_archive$final_country[which(Ithomiini_Zenodo_archive$final_country == "Peru ")] <- "Peru"
Ithomiini_Zenodo_archive$final_country[which(Ithomiini_Zenodo_archive$final_country == "Surinam")] <- "Suriname"
Ithomiini_Zenodo_archive$final_country[which(Ithomiini_Zenodo_archive$final_country == "Trinidad")] <- "Trinidad & Tobago"

# Relevel factors
Ithomiini_Zenodo_archive$final_country <- factor(Ithomiini_Zenodo_archive$final_country)

# Remove useless columns
Ithomiini_Zenodo_archive <- Ithomiini_Zenodo_archive %>% 
  select(ID_obs, Genus, Species, Sub.species, Latitude, Longitude, M.mimicry, F.mimicry, Country = "final_country")

### Clean sp name
Ithomiini_Zenodo_archive$Species <- str_replace(string = Ithomiini_Zenodo_archive$Species, pattern = "AbsPhylo", replacement = "")

### Clean mimicry ring name (add the "-")
Ithomiini_Zenodo_archive$M.mimicry <- as.character(Ithomiini_Zenodo_archive$M.mimicry)
Ithomiini_Zenodo_archive$M.mimicry[which(Ithomiini_Zenodo_archive$M.mimicry == "BANJANAM")] <- "BANJANA-M"
Ithomiini_Zenodo_archive$M.mimicry[which(Ithomiini_Zenodo_archive$M.mimicry == "THABENAF")] <- "THABENA-F"
Ithomiini_Zenodo_archive$M.mimicry <- as.factor(Ithomiini_Zenodo_archive$M.mimicry)
table(Ithomiini_Zenodo_archive$M.mimicry)

Ithomiini_Zenodo_archive$F.mimicry <- as.character(Ithomiini_Zenodo_archive$F.mimicry)
Ithomiini_Zenodo_archive$F.mimicry[which(Ithomiini_Zenodo_archive$F.mimicry == "BANJANAM")] <- "BANJANA-M"
Ithomiini_Zenodo_archive$F.mimicry[which(Ithomiini_Zenodo_archive$F.mimicry == "THABENAF")] <- "THABENA-F"
Ithomiini_Zenodo_archive$F.mimicry <- as.factor(Ithomiini_Zenodo_archive$F.mimicry)
table(Ithomiini_Zenodo_archive$F.mimicry)

### Round coordinates
Ithomiini_Zenodo_archive$Latitude <- round(Ithomiini_Zenodo_archive$Latitude, 3)
Ithomiini_Zenodo_archive$Longitude <- round(Ithomiini_Zenodo_archive$Longitude, 3)

### Update taxonomy

taxonomy_update <- xlsx::read.xlsx(file = "./input_data/Taxonomy_update_2021.xlsx", sheetIndex = 1, stringsAsFactors = F)

Ithomiini_Zenodo_archive[, c("Genus", "Species", "Sub.species")] <- apply(X = Ithomiini_Zenodo_archive[, c("Genus", "Species", "Sub.species")], MARGIN = 2, FUN = as.character)

# Update subspecies

taxonomy_update_subspecies <- taxonomy_update[!is.na(taxonomy_update$Sub.species_old), ]

for(i in 1:nrow(taxonomy_update_subspecies))
{
  match_indices <- which((Ithomiini_Zenodo_archive$Genus %in% taxonomy_update_subspecies$Genus_old[i]) & 
                           (Ithomiini_Zenodo_archive$Species %in% taxonomy_update_subspecies$Species_old[i]) &
                           (Ithomiini_Zenodo_archive$Sub.species %in% taxonomy_update_subspecies$Sub.species_old[i]))
  Ithomiini_Zenodo_archive$Sub.species[match_indices] <- taxonomy_update_subspecies$Sub.species_new[i]
}

# Update species

for(i in 1:nrow(taxonomy_update))
{
  match_indices <- which((Ithomiini_Zenodo_archive$Genus %in% taxonomy_update$Genus_old[i]) & (Ithomiini_Zenodo_archive$Species %in% taxonomy_update$Species_old[i]))
  Ithomiini_Zenodo_archive$Genus[match_indices] <- taxonomy_update$Genus_new[i]
  Ithomiini_Zenodo_archive$Species[match_indices] <- taxonomy_update$Species_new[i]
  
}

Ithomiini_Zenodo_archive$Genus <- as.factor(Ithomiini_Zenodo_archive$Genus)
Ithomiini_Zenodo_archive$Species <- as.factor(Ithomiini_Zenodo_archive$Species)
Ithomiini_Zenodo_archive$Sub.species <- as.factor(Ithomiini_Zenodo_archive$Sub.species)

# Update species name in the list.sp file
load(file = paste0("./input_data/list.sp.RData"))

list.sp$Updated_name <- F
list.sp$Genus_new <- as.character(list.sp$Genus)
list.sp$Species_new <- as.character(list.sp$Species)

for(i in 1:nrow(taxonomy_update))
{
  match_indices <- which((list.sp$Genus %in% taxonomy_update$Genus_old[i]) & (list.sp$Species %in% taxonomy_update$Species_old[i]))
  list.sp$Updated_name[match_indices] <- T
  list.sp$Genus_new[match_indices] <- taxonomy_update$Genus_new[i]
  list.sp$Species_new[match_indices] <- taxonomy_update$Species_new[i]
  
}

list.sp$sp_full_new <- paste(as.character(list.sp$Genus_new), as.character(list.sp$Species_new), sep = ".")
list.sp$Genus_new <- as.factor(list.sp$Genus_new)
list.sp$Species_new <- as.factor(list.sp$Species_new)

# save(list.sp, file = paste0("./input_data/list.sp.RData"))

# Update species name in the phylogeny
phylo.Ithomiini <- readRDS(file = "./input_data/Phylogenies/Final_phylogeny.rds")

replacement_indices <- match(list.sp$Sp_full[list.sp$Updated_name], phylo.Ithomiini$tip.label)
replacement_names <- list.sp$sp_full_new[list.sp$Updated_name]

# Remove NA from species with updated names absent form the phylogeny
na_cleaning <- !is.na(replacement_indices)
replacement_indices <- replacement_indices[na_cleaning]
replacement_names <- replacement_names[na_cleaning]

phylo.Ithomiini$tip.label[replacement_indices] <- replacement_names

# saveRDS(phylo.Ithomiini, file = "./input_data/Phylogenies/Final_phylogeny_with_updated_names.rds")


### Save with subspecies
library("xlsx")
write.xlsx(Ithomiini_Zenodo_archive, file = "./input_data/Ithomiini_records_with_ssp.xlsx",
           sheetName = "Ithomiini records", append = FALSE, row.names = FALSE)
write.csv2(Ithomiini_Zenodo_archive, file = "./input_data/Ithomiini_records_with_ssp.csv", row.names = FALSE)
saveRDS(Ithomiini_Zenodo_archive, file = "./input_data/Ithomiini_records_with_ssp.rds")

# Ithomiini_Zenodo_archive <- readRDS(file = "./input_data/Ithomiini_records_with_ssp.rds")

# Remove duplicates ?
# duplicated(x = select(Ithomiini_Zenodo_archive, -"ID_obs"))
# Ithomiini_Zenodo_archive_no_duplicates <- Ithomiini_Zenodo_archive[!(duplicated(x = select(Ithomiini_Zenodo_archive, -"ID_obs"))), ]


### Remove subspecies
Ithomiini_Zenodo_archive <- select(Ithomiini_Zenodo_archive, -"Sub.species")

### Save without subspecies
library("xlsx")
write.xlsx(Ithomiini_Zenodo_archive, file = "./input_data/Ithomiini_records.xlsx",
           sheetName = "Ithomiini records", append = FALSE, row.names = FALSE)
write.csv2(Ithomiini_Zenodo_archive, file = "./input_data/Ithomiini_records.csv", row.names = FALSE)
saveRDS(Ithomiini_Zenodo_archive, file = "./input_data/Ithomiini_records.rds")

# Remove duplicates ?
# duplicated(x = select(Ithomiini_Zenodo_archive, -"ID_obs"))
# Ithomiini_Zenodo_archive_no_duplicates <- Ithomiini_Zenodo_archive[!(duplicated(x = select(Ithomiini_Zenodo_archive, -"ID_obs"))), ]

### Save without subspecies and duplicates
library("xlsx")
write.xlsx(Ithomiini_Zenodo_archive_no_duplicates, file = "./input_data/Ithomiini_records_no_duplicates.xlsx",
           sheetName = "Ithomiini records", append = FALSE, row.names = FALSE)
write.csv2(Ithomiini_Zenodo_archive_no_duplicates, file = "./input_data/Ithomiini_records_no_duplicates.csv", row.names = FALSE)
saveRDS(Ithomiini_Zenodo_archive_no_duplicates, file = "./input_data/Ithomiini_records_no_duplicates.rds")




### Export current taxonomy

# Generate final table
Species_list <- as_tibble(unique(Ithomiini_Zenodo_archive[, c("Genus", "Species")]))
Species_list$N_rings <- NA
Mimicry_rings <- as_tibble(matrix(nrow = nrow(Species_list), ncol = 8))
names(Mimicry_rings) <- paste0("Mimicry_ring_", 1:8)
Species_list <- cbind(Species_list, Mimicry_rings)

# Retrieve list of mimicry ring per species
for(i in 1: nrow(Species_list))
{
  match_indices <- which((Ithomiini_Zenodo_archive$Genus %in% Species_list$Genus[i]) & 
                          (Ithomiini_Zenodo_archive$Species %in% Species_list$Species[i]))
  mimicry_rings <- unique(c(as.character(Ithomiini_Zenodo_archive$M.mimicry[match_indices]), as.character(Ithomiini_Zenodo_archive$F.mimicry[match_indices])))
  N_rings <- length(mimicry_rings)
  Species_list[i, 4:(4 + N_rings - 1)] <- mimicry_rings
  Species_list$N_rings[i] <- N_rings
}

# Save
write.xlsx(Species_list, file = "./input_data/Ithomiini_taxonomy.xlsx",
           sheetName = "Ithomiini taxonomy", append = FALSE, row.names = FALSE, showNA = FALSE)
write.csv2(Species_list, file = "./input_data/Ithomiini_taxonomy.csv", row.names = FALSE, na = "")
saveRDS(Species_list, file = "./input_data/Ithomiini_taxonomy.rds")
       

load(file = "./input_data/list.models.RData")

table(list.models$Mimicry.model)[order(table(list.models$Mimicry.model), decreasing = T)]

list.models$Sp_full[list.models$Mimicry.model == "ILLINISSA"]
