# This script contains a list of helper functions which can be used to analyse and plot the ZooMSS output


# Sum ZooMSS output across size bins
fZooMSS_SumSize = function(list_in) {
  out <- map(list_in, function(x) apply(x, 1, sum)) # Sum ZooMSS output across the size bins
  return(out)
}

# Sum ZooMSS output across species bins
fZooMSS_SumSpecies = function(list_in) {
  out <- map(list_in, function(x) apply(x, 2, sum)) # Sum ZooMSS output across the species bins
  return(out)
}

# Summarise the biomass for each grid-cell by species
fZooMSS_SpeciesBiomass = function(res, w) {
  Biomass <- map(res, function(x) apply(sweep(x, 2, w, '*'), 1, sum))
  return(Biomass)
}

# Sum ZooMSS output across all species and sizes
fZooMSS_SumAll = function(list_in) {
  out <- unlist(map(list_in, function(x) sum(x))) # Sum ZooMSS output across the species bins
  return(out)
}

# Convert Abundance to Biomass for all species and weight classes
fZooMSS_Biomass <- function(res, w) {
  Biomass <- map(res, function(x) sweep(x, 2, w, '*')) # Biomass in grams
  return(Biomass)
}

# Convert Abundance to Carbon Biomass for all species and weight classes
fZooMSS_CarbonBiomass <- function(res, w, carbon) {
  Biomass <- map(res, function(x) sweep(x, 2, w, '*'))  # Biomass in grams (WW)
  Biomass <- map(Biomass, function(x) sweep(x, 1, carbon, '*')) # Now convert to Carbon
  return(Biomass)
}

# Summarise the biomass for each grid-cell by size-class
fZooMSS_SizeBiomass = function(res,w) {
  Biomass <- map(res, function(x) apply(sweep(x, 2, w, '*'), 2, sum))
  return(Biomass)
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





