### Function to compute Faith's Phylogenetic Diversity from raster stack of SDM continuous outputs ###

# Author: Maël Doré
# Contact: mael.dore@gmail.com


#### 1/ Function that matches species list in stack and phylogeny and clean them ####

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


#### 2/ Function to compute Faith's Phylogenetic Diversity from raster stack of SDM continuous outputs ####

# Inputs = Raster Stack of species/OMU probability of presence
#          Phylogeny of the clade

# Output = Raster Layer of Faith's PD
# Print = Lists of species removed from stack because absent from the phylogeny, and respectively.

compute_Faith_PD <- function(proba_stack, phylo)
{
  # Create function to aggregate probabilities to higher hierarchical level (aggregate pixel, or go up on a phylogenetic tree)
  aggreg_prob <- function(x, na.rm) 
  { 
    y <- 1 - prod(1 - x) 
    return(y) # Output
  }
  
  # Match raster Stack and Phylogeny species lists
  clean_stack_and_phylo <- match_stack_and_phylo(proba_stack, phylo)
  proba_stack_for_phylo <- clean_stack_and_phylo[[1]]
  pruned_tree <- clean_stack_and_phylo[[2]]
  
  # Get matrix of probabilities of presence per species (columns) per community (rows)
  proba_mat <- getValues(proba_stack_for_phylo)
  
  # Final nb of taxonomic units (species)
  n_leaves <- length(pruned_tree$tip.label)
  
  # Generate matrix to store info on edges/branches from the pruned phylogeny
  branches = as.data.frame(matrix(NA, nrow(pruned_tree$edge), ncol = 4)) 
  names(branches) <- c("Starting nod", "Ending nod", "Length", "Proba_presence")
  branches[, 1:2] = pruned_tree$edge # Retrieve starting and ending node 
  branches[, 3] = round(pruned_tree$edge.length, 4) # Retrieve edge length
  
  # Make a loop per community
  all_PD <- NA
  for (k in 1:nrow(proba_mat)) {
    proba_com <- proba_mat[k,] # Extract the proba for the k community
    
    # Compute PD only if no NA is present in the community
    if (any(is.na(proba_com))) {
      
      all_PD[k] <- NA # if NA present, PD = NA
      
    } else {
      
      # Compute probability of presence of each edge/branch (i.e., a species descending from this branch) in this community
      for (i in 1:nrow(branches)) {
        descending_species = geiger::tips(pruned_tree, branches[i, 2]) # Retrieve the set of species descending from branch i
        index <- which(colnames(proba_mat) %in% descending_species) # Find position of the species in the stack/matrix
        prob_edge <- aggreg_prob(proba_com[index]) # Compute probability of presence of the edge
        branches[i, 4] <- prob_edge # Store info
      }
      PD <- round(sum(branches[,3] * branches[,4]), 4) # Compute PD as the weighted sum of the edge length (weighted by the probability of presence of this edge in the community)
      all_PD[k] <- PD
    }
    
    # Show k every 1000 iterations and save a backup
    if (k %% 1000 == 0) {
      cat(paste0(Sys.time(), " - ", k," on ",nrow(proba_mat),"\n"))
      # save(all_PD, file = "./outputs/Indices_Maps/PD_backup.RData")
    }
  }
  
  # Write community PD in a raster
  PD_raster <- proba_stack[[1]]
  PD_raster@data@values <- all_PD
  
  # Repair issue with max value
  PD_raster@data@max <- max(PD_raster[], na.rm = T)
  
  return(PD_raster)
  
}
