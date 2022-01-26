### Functions to compute MPD from raster stack of SDM continuous outputs ###

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


#### 2/ Function to compute MPD from raster stack of SDM continuous outputs ####

# Inputs = Raster Stack of species/OMU probability of presence
#          Phylogeny of the clade

# Output = Raster Layer of Mean pairwise Phylogenetic Distances
# Print = Lists of species removed from stack because absent from the phylogeny, and respectively.


compute_MPD <- function(proba_stack, phylo)
{
  # Match raster Stack and Phylogeny species lists
  clean_stack_and_phylo <- match_stack_and_phylo(proba_stack, phylo)
  proba_stack_for_phylo <- clean_stack_and_phylo[[1]]
  pruned_tree <- clean_stack_and_phylo[[2]]
  
  # Get matrix of probabilities of presence per species (columns) per community (rows)
  proba_mat <- getValues(proba_stack_for_phylo)
  
  # Final nb of taxonomic units (species)
  n_leaves <- length(pruned_tree$tip.label)
  
  # Compute the patristic phylogenetic distances
  phylo_dist_mat <- ape::cophenetic.phylo(x = pruned_tree)
  
  # Loop to compute indices per community
  all_MPD <- NA
  for (k in 1:nrow(proba_mat)) {
    proba_com <- proba_mat[k,] # Extract the probabilities for the k community
    
    # Compute MPD only if no NA is present in the community
    if (any(is.na(proba_com))) {
      
      all_MPD[k] <- NA # if NA present, MPD = NA
      
    } else {
      
      # Compute a matrix of probability of existence of pairs in the community
      proba_pairs_mat <- matrix(data = NA, ncol = n_leaves, nrow = n_leaves)
      for (i in 1:n_leaves) {
        for (j in 1:n_leaves) {
          pair_probas <- c(proba_com[i], proba_com[j])
          proba_pairs_mat[i,j] <- prod(pair_probas) # Probability of presence of a pair = product of probability of presence of each of the species of the pair
        }
      }
      
      # Multiply them to weight phylogenetic distance by probability of presence of the pair of species
      weighted_phylo_dist_mat <- proba_pairs_mat * phylo_dist_mat
      
      # Extract unique pair values and get the mean
      total_phylo_dist <- sum(weighted_phylo_dist_mat[lower.tri(weighted_phylo_dist_mat, diag = F)]) # Extract weighted pair values and sum them
      MPD <- total_phylo_dist/sum(proba_pairs_mat[lower.tri(proba_pairs_mat, diag = F)]) # Divide by the sum of the weights to compute the weighted mean
      all_MPD[k] <- MPD # Store the MDP of each community in the final vector
      
    }
    # Show k every 1000 iterations and save a back-up file
    if (k %% 1000 == 0) 
    {
      cat(paste0(Sys.time(), " - ", k," on ",nrow(proba_mat),"\n"))
      # save(all_MPD, file = "./outputs/Indices_Maps/MPD_backup.RData")
    }
  }
  
  # Generate a mask for terrestrial areas
  continental_mask <- (readAll(calc(proba_stack, fun = sum)) >= 0) - 1
  
  # Write community MPD in a raster with null values as background for terrestrial areas
  MPD_raster <- continental_mask
  MPD_raster@data@values[!is.na(all_MPD)] <- all_MPD[!is.na(all_MPD)]
  
  # For a raster with only communities with data
  # MPD_raster <- proba_stack[[1]]
  # MPD_raster@data@values <- all_MPD
  
  # Repair issue with max value
  MPD_raster@data@max <- max(MPD_raster[], na.rm = T)
  
  return(MPD_raster)
}
