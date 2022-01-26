### Function to compute Fair Proportions (FP) ###

# Author: Maël Doré
# Contact: mael.dore@gmail.com

### 2 levels

# Species FP from a phylogeny
# Total and mean Community FP from raster stack of SDM continuous outputs
   # Total = sum of species FP
   # Mean = standardized by species richness

# At species level = Sum of all modified branches length from the root to the species. The most original species have the highest scores
# At community level = Sum of FP values of each species present in the community. The community with the most original species has the highest scores
# Can be standardized by the number of species to get the mean species FP index

# Species FP can be computed for the full phylogeny provided, or the pruned version with only the species presents in the raster Stack

### 1/ Compute Species FP from a phylogeny ####

# Inputs = Phylogeny of the clade

# Output = Dataframe with species Fair-proportions

compute_species_FP <- function(phylo)
{
  # Generate table that store information on branches
  phylo_branches_df = as.data.frame(matrix(NA, nrow(phylo$edge), ncol = 5)) 
  names(phylo_branches_df) <- c("Starting_node","Ending_node","Branch_Length","Nb_desc_sp","FP")
  
  # Retrieve starting and ending nodes for each branches 
  phylo_branches_df[, c("Starting_node","Ending_node")] = phylo$edge 
  # Retrieve edge length
  phylo_branches_df$Branch_Length = phylo$edge.length 
  # Retrieve number of descending species/leaves (terminal nods)
  for (i in 1:nrow(phylo$edge)) 
  { 
    phylo_branches_df$Nb_desc_sp[i] = length(geiger::tips(phylo, phylo_branches_df$Ending_node[i]))
  }
  # Divide branch length by nb of descendant species for FP
  phylo_branches_df$FP = phylo_branches_df$Branch_Length/phylo_branches_df$Nb_desc_sp 
  
  # Generate table to store species Fair-Proportions
  species_FP_df <- data.frame(Taxon = phylo$tip.label, FP = NA)
  
  # Compute species Fair-Proportions
  for (i in 1:length(phylo$tip.label))
  {
    # Retrieve species name
    sp <- phylo$tip.label[i] 
    
    # Get indices of nods on the path from the root to the species terminal nod
    sp_path <- ape::nodepath(phy = phylo, 
                             from = phylo_branches_df$Starting_node[which.max(phylo_branches_df$Nb_desc_sp)], # Root
                             to = which(phylo$tip.label == sp)) # Tip
    
    # Compute FP index for the species by summing FP values of the branches on the path to the species tip in the phylogeny 
    species_FP_df$FP[species_FP_df$Taxon == sp] <- sum(phylo_branches_df[which(phylo_branches_df$Ending_node %in% sp_path), "FP"]) 
  }
  
  return(species_FP_df)
}


#### 2/ Function that matches species list in stack and phylogeny and clean them ####

# Inputs = Raster Stack of species/OMU probability of presence
#          Phylogeny of the clade

# Outputs = Stack with only species form the phylogeny
#           Pruned phylogeny with only species in the raster Stack
# Print = Lists of species removed from stack because absent from the phylogeny, and respectively.

match_stack_and_phylo <- function (proba_stack, phylo)
{
  # List species in the stack
  species_list <- names(proba_stack)
  
  # Remove species not in the phylogeny from the stack
  proba_stack_for_phylo <- subset(proba_stack, which(species_list %in% phylo$tip.label))
  
  # Print species remove from stack because they are not in the phylogeny
  not_in_phylo <- species_list[!(species_list %in% phylo$tip.label)]
  if (length(not_in_phylo) > 0)
  {
    cat(paste0("\n", length(not_in_phylo), " species removed from stack because they are absent from the phylogeny:\n\n"))
    print(not_in_phylo)
    cat("\n")
  }
  
  # Prune the phylogeny to keep only species in the stack
  pruned_tree <- ape::keep.tip(phylo, names(proba_stack_for_phylo))
  
  # Print species remove from phylogeny because they are not in the stack
  not_in_stack <- phylo$tip.label[!(phylo$tip.label %in% names(proba_stack_for_phylo))]
  if (length(not_in_stack) > 0)
  {
    cat(paste0("\n", length(not_in_stack), " species removed from the phylogeny because they are absent from the Raster Stack:\n\n"))
    print(not_in_stack)
    cat("\n")
  }
  
  # Return cleaned stack and phylogeny in a list
  cleaned_stack_and_phylo <- list(proba_stack_for_phylo, pruned_tree)
  return(cleaned_stack_and_phylo)
}


### 3/ Compute Total and mean Community FP from raster stack of SDM continuous outputs ####

# Inputs = Raster Stack of species/OMU probability of presence
#          Phylogeny of the clade

# Output = Raster Layer of Community FP. Mean or Total.
# Print = Lists of species removed from stack because absent from the phylogeny, and respectively.

# Choose between "total" FP and "mean" FP with the "index" argument to select which index to compute
# Choose between "full" and "pruned" with the "phylo_type" argument to select which phylogeny to use to compute species FP

compute_community_FP <- function(proba_stack, phylo, phylo_type = "full", index = "total")
{
  # Match raster Stack and Phylogeny species lists
  clean_stack_and_phylo <- match_stack_and_phylo(proba_stack, phylo)
  proba_stack_for_phylo <- clean_stack_and_phylo[[1]]
  
  # Get matrix of probabilities of presence per species (columns) per community (rows)
  proba_mat <- getValues(proba_stack_for_phylo)
  
  # By default, compote the species FP on the full phylogeny
  if (phylo_type == "full") 
  {
    # Compute the species FP based on the full phylogeny (because the pruned tree may exhibit sampling bias)
    species_FP_df <- compute_species_FP(phylo = phylo)
    # Keep only the species present in the stack
    species_FP_df <- species_FP_df[species_FP_df$Taxon %in% names(proba_stack_for_phylo), ]
  }
  
  # If requested, use the pruned phylogeny with only the species from the stack to compute the species FP
  if (phylo_type == "pruned")
  {
    pruned_tree <- clean_stack_and_phylo[[2]]
    species_FP_df <- compute_species_FP(phylo = pruned_tree)
  }
  
  # Make a loop per community to compute community FP
  all_FP <- rep(NA, nrow(proba_mat))
  for (k in 1:nrow(proba_mat)) 
  {
    proba_com <- proba_mat[k,] # Extract the probabilities of presence for the kth community
    
    # Compute FP only if no NA is present in the community
    if (any(is.na(proba_com))) 
    {
      all_FP[k] <- NA # if NA present, FP = NA
    } else {
      
      if(index == "total") # Compute FP as the weighted sum of the FP index of species, weighted by their probability of presence in the community
      {
        FP <- sum(species_FP_df$FP * proba_com) 
      } else { # Compute mean FP as the weighted mean of the FP index of species, weighted by their probability of presence in the community
        FP <- sum(species_FP_df$FP * proba_com) / sum(proba_com)
      }
      all_FP[k] <- FP
    }
    
    # Show k every 1000 iterations and save a back-up file
    if (k %% 1000 == 0) {
      cat(paste0(Sys.time(), " - ", k," on ",nrow(proba_mat),"\n"))
      # save(all_FP, file = "./outputs/Indices_maps/backup_FP.RData")
    }
  }
  
  # Write community FP in a raster
  FP_raster <- proba_stack[[1]]
  FP_raster@data@values <- all_FP
  
  # Repair issue with max value
  FP_raster@data@max <- max(FP_raster[], na.rm = T)
  
  return(FP_raster)
}
