# This script contains a list of helper functions which can be used to analyse and plot the ZooMSS output

# Summarise the biomass for each gridcell by species
fZooMSS_GlobalBiomass = function(res,w) {
  bioms <- matrix(NA, nrow = length(res), ncol = length(res[[1]][,1])) # Preallocate Fish Proportions
  for(j in 1:length(res)){
    bioms[j,] = apply(sweep(res[[j]], 2, w, '*'), 1, sum) # Get the reference/control data
  }
  return(bioms)
}

# Function to calculate the mean of the last 50 % of the model
fZooMSS_AveOutput = function(x){
  ave_x <- colMeans(x[(ceiling(0.5*(dim(x)[1])):dim(x)[1]),,], dims = 1)
  return(ave_x)
}

# Remove nonsense attributes if we are working for speed and memory efficiency.
untibble <- function (tibble) {
  data.frame(unclass(tibble), check.names = FALSE, stringsAsFactors = FALSE)
}  ## escape the nonsense