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
fZooMSS_SpeciesBiomass = function(res, vmdl) {
  # if (dim(res[[1]])[2] != length(mdl$param$w)){print("error")}
  Biomass <- map(res, function(x) apply(sweep(x, 2, vmdl$param$w, '*'), 1, sum))
  return(Biomass)
}

# Sum ZooMSS output across all species and sizes
fZooMSS_SumAll = function(list_in) {
  out <- unlist(map(list_in, function(x) sum(x))) # Sum ZooMSS output across the species bins
  return(out)
}

# Convert Abundance to Biomass for all species and weight classes
fZooMSS_Biomass <- function(res, vmdl) {
  if (dim(res[[1]])[2] != length(vmdl$param$w)){print("error")}
  Biomass <- map(res, function(x) sweep(x, 2, vmdl$param$w, '*')) # Biomass in grams
  return(Biomass)
}

# Convert Abundance to Carbon Biomass for all species and weight classes
fZooMSS_CarbonBiomass <- function(res, vmdl) {
  if (dim(res[[1]])[2] != length(vmdl$param$w)){print("error")}
  Biomass <- map(res, function(x) sweep(x, 2, vmdl$param$w, '*'))  # Biomass in grams (WW)
  Biomass <- map(Biomass, function(x) sweep(x, 1, vmdl$param$Groups$Carbon, '*')) # Now convert to Carbon
  return(Biomass)
}

# Convert Abundance to Carbon Biomass for all species and weight classes
fZooMSS_SpeciesCarbonBiomass <- function(res, vmdl) {
  if (dim(res[[1]])[2] != length(vmdl$param$w)){print("error")}
  Biomass <- map(res, function(x) sweep(x, 2, vmdl$param$w, '*'))  # Biomass in grams (WW)
  Biomass <- map(Biomass, function(x) sweep(x, 1, vmdl$param$Groups$Carbon, '*')) # Now convert to Carbon
  Biomass <- fZooMSS_SumSize(Biomass)

    return(Biomass)
}

# Summarise the biomass for each grid-cell by size-class
fZooMSS_SizeBiomass = function(res,w) {
  if (dim(res[[1]])[2] != length(w)){print("error")}
  Biomass <- map(res, function(x) apply(sweep(x, 2, w, '*'), 2, sum))
  return(Biomass)
}

# Sum ZooMSS output across size bins
fZooMSS_ExtractSizeRange = function(list_in, minb, maxb) {
  out <- map(list_in, function(x) x[,minb:maxb] )
  return(out)
}


# Function to calculate the mean of the last 50 % of the model
fZooMSS_AveOutput = function(x, prop = 0.5){
  ave_x <- colMeans(x[(ceiling(dim(x)[1] - prop*(dim(x)[1])):dim(x)[1]),,], dims = 1)
  return(ave_x)
}

# Remove nonsense attributes if we are working for speed and memory efficiency.
untibble <- function (tibble) {
  data.frame(unclass(tibble), check.names = FALSE, stringsAsFactors = FALSE)
}  ## escape the nonsense


# At the moment you need to subset the zoomss data by species or size in order to use this function.
# I will rewrite sometime to include other variables but at the moment its only for 2D data
fZooMSS_Convert2Tibble <- function(li, vmdl){
  df <- as_tibble(matrix(unlist(li), nrow=length(li), byrow=T), .name_repair = "unique") %>%
    rename_with(~vmdl$param$Groups$Species)
    return(df)
}

fZooMSS_AddEnviro <- function(Zoo, venviro){
  df <- Zoo %>%
    mutate(cellID = 1:n()) %>% # Create a cellID
    left_join(dplyr::select(venviro, cellID, chlo, sst, phyto_int, phyto_slope, phyto_max, Lat, Lon, geometry), by = "cellID") %>%
    rename(SST = sst, Chl = chlo) %>%
    mutate(Chl_log10 = log10(Chl))
  return(df)
}


# Return the diet matrix as a long tibble
fZooMSS_MakeDietTibble <- function(mat, mdl){
  suppressMessages(
    out <- as_tibble(mat, .name_repair = "unique") %>%
      rename_with(~c("Phyto_Small", "Phyto_Med", "Phyto_Large", mdl$param$Groups$Species)) %>%
      mutate(Predator = mdl$param$Groups$Species) %>%
      pivot_longer(cols = Phyto_Small:Fish_Large, names_to = "Prey", values_to = "Diet")
  )
  return(out)
}


PPMR_plot = function(dat){

  min_size = min(dat$model$param$Groups$W0) # smallest size class
  max_size = max(dat$model$param$Groups$Wmax) # largest size class
  w = 10^(seq(from = min_size, to = max_size, 0.1)) # all size classes

  # Calculate PPMR (beta) table, where dim1 = group, dim2 = body size with
  # value being PPMR for that body size (this is not realised PPMR - not
  # emergent from diet but calculated from m-values and Wirtz, 2012 equation)
  D.z = 2*(3*(w)*1e12/(4*pi))^(1/3) # convert body mass g to ESD (um)
  zoo_m = dat$model$param$Groups$PPMRscale # pull out PPMR scaling values from parameter table
  betas =  log10(t(sapply(zoo_m, function(x){(exp(0.02*log(D.z)^2 - x + 1.832))^3}))) # Convert m to betas, using Wirtz 2012 equation
  betas = betas[-which(is.na(dat$model$param$Groups$PPMRscale)),] # remove fish rows

  ## Modify beta matrix for larvaceans and salps - all size classes for these groups feed on same prey, so log10PPMR increases by 0.1 for each 0.1 log10 size interval
  betas[which(dat$model$param$Groups$Species=="Larvaceans"),45:75] <- betas[which(dat$model$param$Groups$Species=="Larvaceans"),44] + seq(0.1,3.1,0.1) # Larvaceans (index 44 in w vector is smallest size class, 75 is maximum size class)
  betas[which(dat$model$param$Groups$Species=="Salps"),61:121] <- betas[which(dat$model$param$Groups$Species=="Salps"),61] + seq(0.1,6.1,0.1) # Larvaceans (index 61 in w vector is smallest size class, 121 is maximum size class

  # Calculate ave abundances across oligo/eutro grid squares, then calculate ave
  # biomass and proportion of total zoo biomass that is from each group size class

  ave = matrix(0, nrow = dim(dat$model$param$Groups)[1], ncol = length(w))
  for(i in 1:length(dat$abundances)){
    ave = ave + dat$abundances[[i]]/length(dat$abundances)
  }
  ave_biom = sweep(ave, 2, w, "*") # Calculate oligo biomass for zoo groups
  ave_biom = ave_biom[-which(is.na(dat$model$param$Groups$PPMRscale)),] # remove rows for fish
  beta_props = ave_biom/sum(ave_biom) # Calculate fraction of zoo biomass in each group, in each size class

  out <- list()
  out[[1]] <- betas
  out[[2]] <- beta_props
  names(out) <- c("betas", "beta_props")

  temp <- density(betas, weights = beta_props)

  out <- tibble("x" = temp$x, "y" = temp$y, "mn_beta" = sum(beta_props*betas))

  spbeta_props = ave_biom/rowSums(ave_biom) # Species specific proportions
  spPPMR <- tibble("Species" = as.factor(dat$model$param$Groups$Species[-which(is.na(dat$model$param$Groups$PPMRscale))]), "Betas" = rowSums(spbeta_props*betas), "y" = NA) # Get species-specific PPMR

  for (s in 1:length(spPPMR$Species)){
    spPPMR$y[s] <- out$y[which.min(abs(out$x - spPPMR$Betas[s]))]
  }

  spPPMR <- spPPMR %>%
    mutate(y = y * 0) %>%
    bind_rows(spPPMR)

  out2 <- list()
  out2[[1]] <- out
  out2[[2]] <- spPPMR

  return(out2)
}



## Function to calculate slope intercept and maximum size of phytoplankton spectrum, for zooplankton
## resolved size spectrum model (Heneghan et al. in prep).

# Last updated 17th September 2020

fZooMSS_CalculatePhytoParam = function(df){ # chlo is chlorophyll concentration in mg m^-3

  ## Calculate pico, nano, micro phytoplankton proportions of total chlorophyll
  ## BREWIN ET AL., 2015
  pico <- (0.13*(1-exp(-0.8/0.13*df$chlo)))/df$chlo
  nano <- (0.77*(1-exp(-0.94/0.77*df$chlo)))/df$chlo - pico
  micro <- (df$chlo - 0.77*(1-exp(-0.94/0.77*df$chlo)))/df$chlo

  ## Convert total chlorophyll to g m^-3 total wet weight - biomass
  ## Allocate total chlorophyll to the three size classes
  c_chl <- ((df$chlo^0.89)*(10^1.79))/df$chlo # chlo:carbon ratio, from Mara??on et al. 2014
  tot_biom_c <- c_chl*df$chlo/1000 # (convert to grams carbon)
  tot_biom <- tot_biom_c*(1/0.1) # convert to grams wet weight, assuming 0.1 C:ww

  # Break up total biom into pico, nano and micro
  df$pico_biom <- pico*tot_biom
  df$nano_biom <- nano*tot_biom
  df$micro_biom <- micro*tot_biom

  ## Find abundances at boundaries of pico, nano size ranges, by analytically
  ## solving integral of N = aw^b

  w_0 <- -14.5 # log minimum size of picophytoplankton
  w_1 <- -11.5 # log minimum size of nanophytoplankton (max size of pico also)
  w_2 <- -8.5 # log minimum size of macrophytoplankton (max size of nano also)

  df$phyto_slope <- (log10(df$pico_biom) - log10(df$nano_biom) - w_1 + w_2)/(w_1 - w_2)  # Calculate slope
  df$phyto_int <- log10(df$pico_biom*(df$phyto_slope+1)/((10^(w_1))^(df$phyto_slope+1) - (10^(w_0))^(df$phyto_slope+1))) # Calculate intercept

  ## Calculate maximum size
  df$phyto_max <- 0.1*round((-8.4 + 2*micro)/0.1) # Maximum size depends on the proportion of micro
  max_phyto <- rep(-7, length(df$chlo)) # Set -7 to be the max possible size for phyto
  df$phyto_max <- pmin(max_phyto, df$phyto_max)

  return(df)
}

