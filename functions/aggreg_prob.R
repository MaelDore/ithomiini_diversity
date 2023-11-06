#### Function to compute probability of presence at larger spatial or hierarchical scale ####


aggreg_prob = function(x, na.rm) { 
  y <- 1-prod(1-x) # Probability of presence of any unit = probability of presence of at least element/in at least one cell = opposite of probability of absence of all units/within all cells
  return(y) # Output
}
