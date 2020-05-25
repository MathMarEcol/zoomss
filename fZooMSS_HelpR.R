# This script contains a list of helper functions which can be used to analyse and plot the ZooMSS output

# Summarise the biomass for each gridcell by species
fZooMSS_GlobalBiomass = function(res,w) {
  bioms <- matrix(NA, nrow = length(res), ncol = length(res[[1]][,1])) # Preallocate Fish Proportions
  for(j in 1:length(res)){
    bioms[j,] = apply(sweep(res[[j]], 2, w, '*'), 1, sum) # Get the reference/control data
  }
  return(bioms)
}
