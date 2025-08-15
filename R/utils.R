#' ZooMSS Utility Functions for Analysis and Post-Processing
#'
#' @title Collection of helper functions for analyzing ZooMSS model outputs
#' @description This file contains utility functions for processing, analyzing, and
#'   transforming ZooMSS model outputs for visualization and interpretation.
#' @details The utility functions in this file provide tools for:
#'   - Converting between abundance and biomass
#'   - Aggregating results across size classes or functional groups
#'   - Calculating ecological metrics (trophic levels, PPMR)
#'   - Processing environmental data for model input
#'   - Data format conversions for analysis workflows
#'
#'   These functions are essential for the ZooMSS analysis pipeline and help users
#'   work with model outputs in different formats depending on their research needs.

#' Sum ZooMSS Output Across Size Bins
#'
#' @title Aggregate ZooMSS abundances across all size classes
#' @description Sums abundance values across all size classes for each functional group,
#'   providing total abundance per group.
#' @details This function collapses the size dimension of ZooMSS output by summing
#'   across all size classes. Useful for analyzing total abundance patterns without
#'   size structure detail.
#'
#' @param list_in List of abundance matrices (typically from multiple spatial cells)
#'
#' @return List of vectors with total abundance per functional group
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model
#' results <- zoomss_model(input_params, Groups, SaveTimeSteps = FALSE)
#' total_abundances <- zSumSize(results$abundances)
#' }
#'
zSumSize = function(list_in) {
  out <- purrr::map(list_in, function(x) apply(x, 1, sum)) # Sum ZooMSS output across the size bins
  return(out)
}

#' Sum ZooMSS Output Across Functional Groups
#'
#' @title Aggregate ZooMSS abundances across all functional groups
#' @description Sums abundance values across all functional groups for each size class,
#'   providing total abundance per size class.
#' @details This function collapses the functional group dimension by summing across
#'   all groups for each size class. Useful for analyzing community size spectrum
#'   patterns without functional group detail.
#'
#' @param list_in List of abundance matrices (typically from multiple spatial cells)
#'
#' @return List of vectors with total abundance per size class
#' @export
#'
zSumSpecies = function(list_in) {
  out <- purrr::map(list_in, function(x) apply(x, 2, sum)) # Sum ZooMSS output across the species bins
  return(out)
}

#' Calculate Species-Level Biomass
#'
#' @title Sum biomass across size classes for each functional group
#' @description Converts abundance to biomass and sums across all size classes
#'   to provide total biomass per functional group per spatial cell.
#' @details This function combines abundance-to-biomass conversion with size-class
#'   aggregation in one step, providing species-level biomass summaries useful
#'   for spatial analyses and ecological comparisons.
#'
#' @param res List of abundance matrices from ZooMSS output
#' @param vmdl ZooMSS model object containing weight vector (param$w)
#'
#' @return List of vectors with total biomass per functional group (grams wet weight)
#' @export
#'
zSpeciesBiomass = function(res, vmdl) {
  # if (dim(res[[1]])[2] != length(mdl$param$w)){print("error")}
  Biomass <- purrr::map(res, function(x) apply(sweep(x, 2, vmdl$param$w, '*'), 1, sum))
  return(Biomass)
}

#' Sum All ZooMSS Output
#'
#' @title Sum abundances across all groups and size classes
#' @description Calculates total abundance across all functional groups and size classes,
#'   providing a single abundance value per spatial cell.
#' @details This function provides the most aggregated view of ZooMSS output by
#'   summing across both functional groups and size classes. Useful for comparing
#'   total community abundance between locations or time periods.
#'
#' @param list_in List of abundance matrices (typically from multiple spatial cells)
#'
#' @return Vector of total abundance values (one per spatial cell)
#' @export
#'
zSumAll = function(list_in) {
  out <- unlist(purrr::map(list_in, function(x) sum(x))) # Sum ZooMSS output across the species bins
  return(out)
}

#' Convert Abundance to Biomass
#'
#' @title Convert ZooMSS abundance matrices to biomass by multiplying by body weights
#' @description Converts abundance data to wet weight biomass by multiplying abundances
#'   by the corresponding body weights for each size class.
#' @details This function transforms abundance matrices to biomass by applying the
#'   weight vector across size classes. Essential for analyses requiring biomass
#'   units rather than abundance counts.
#'
#' @param res List of abundance matrices from ZooMSS output
#' @param vmdl ZooMSS model object containing weight vector (param$w)
#'
#' @return List of biomass matrices in grams wet weight
#' @export
#'
zBiomass <- function(res, vmdl) {
  # if (dim(res[[1]])[2] != length(vmdl$param$w)){print("error")}
  Biomass <- purrr::map(res, function(x) sweep(x, 2, vmdl$param$w, '*')) # Biomass in grams
  return(Biomass)
}

#' Convert Abundance to Carbon Biomass
#'
#' @title Convert ZooMSS abundances to carbon biomass across all size classes
#' @description Converts abundance data to carbon biomass by multiplying by body weights
#'   and then by carbon content factors for each functional group.
#' @details This function performs a two-step conversion:
#'   1. Abundance to wet weight biomass (using body weights)
#'   2. Wet weight to carbon biomass (using group-specific carbon content)
#'
#'   Carbon biomass is essential for biogeochemical analyses and comparisons
#'   with field data that are often reported in carbon units.
#'
#' @param res List of abundance matrices from ZooMSS output
#' @param vmdl ZooMSS model object containing weight vector and carbon content factors
#'
#' @return List of carbon biomass matrices (grams carbon)
#' @export
#'
zCarbonBiomass <- function(res, vmdl) {
  if (dim(res[[1]])[2] != length(vmdl$param$w)){print("error")}
  Biomass <- purrr::map(res, function(x) sweep(x, 2, vmdl$param$w, '*'))  # Biomass in grams (WW)
  Biomass <- purrr::map(Biomass, function(x) sweep(x, 1, vmdl$param$Groups$Carbon, '*')) # Now convert to Carbon
  return(Biomass)
}

#' Convert Abundance to Species-Level Carbon Biomass
#'
#' @title Convert abundances to carbon biomass and sum across size classes
#' @description Converts abundance data to carbon biomass and then sums across all
#'   size classes to provide total carbon biomass per functional group.
#' @details This function combines carbon biomass conversion with size-class aggregation:
#'   1. Converts abundance to wet weight biomass
#'   2. Converts to carbon biomass using group-specific factors
#'   3. Sums across all size classes for each functional group
#'
#'   Provides species-level carbon biomass useful for ecological stoichiometry
#'   and biogeochemical cycle analyses.
#'
#' @param res List of abundance matrices from ZooMSS output
#' @param vmdl ZooMSS model object containing weight vector and carbon content factors
#'
#' @return List of vectors with total carbon biomass per functional group (grams carbon)
#' @export
#'
zSpeciesCarbonBiomass <- function(res, vmdl) {
  if (dim(res[[1]])[2] != length(vmdl$param$w)){print("error")}
  Biomass <- purrr::map(res, function(x) sweep(x, 2, vmdl$param$w, '*'))  # Biomass in grams (WW)
  Biomass <- purrr::map(Biomass, function(x) sweep(x, 1, vmdl$param$Groups$Carbon, '*')) # Now convert to Carbon
  Biomass <- zSumSize(Biomass)

  return(Biomass)
}

#' Calculate Size-Class Biomass
#'
#' @title Sum biomass across functional groups for each size class
#' @description Converts abundance to biomass and sums across all functional groups
#'   to provide total biomass per size class per spatial cell.
#' @details This function provides size-class-level biomass by summing across
#'   functional groups. Useful for analyzing community size structure and
#'   comparing size spectrum patterns between locations.
#'
#' @param res List of abundance matrices from ZooMSS output
#' @param w Vector of body weights for each size class (grams)
#'
#' @return List of vectors with total biomass per size class (grams wet weight)
#' @export
#'
zSizeBiomass = function(res,w) {
  if (dim(res[[1]])[2] != length(w)){print("error")}
  Biomass <- purrr::map(res, function(x) apply(sweep(x, 2, w, '*'), 2, sum))
  return(Biomass)
}

#' Extract Size Range from ZooMSS Output
#'
#' @title Extract specific size class range from abundance matrices
#' @description Subsets ZooMSS output to include only specified size classes,
#'   useful for focusing analysis on particular size ranges.
#' @details This function extracts a subset of size classes from the full
#'   ZooMSS output matrices. Useful for analyzing specific size ranges
#'   (e.g., microzooplankton, mesozooplankton) or excluding boundary effects
#'   from model analysis.
#'
#' @param list_in List of abundance matrices from ZooMSS output
#' @param minb Minimum size class index to extract
#' @param maxb Maximum size class index to extract
#'
#' @return List of abundance matrices with only specified size classes
#' @export
#'
zExtractSizeRange = function(list_in, minb, maxb) {
  out <- purrr::map(list_in, function(x) x[,minb:maxb] )
  return(out)
}


#' Calculate Average Output from Model Time Series
#'
#' @title Calculate mean of final portion of ZooMSS time series
#' @description Calculates the mean of the final portion (default 50%) of a time series
#'   to obtain equilibrium values after model spin-up period.
#' @details This function removes the initial transient period from time series data
#'   and calculates the mean of the remaining portion, providing representative
#'   steady-state values. Essential for obtaining equilibrium abundances, growth rates,
#'   and other model outputs after the model has reached dynamic equilibrium.
#'
#' @param x 3D array with dimensions (time, groups, size_classes)
#' @param prop Proportion of final time series to average (default: 0.5)
#'
#' @return 2D array with averaged values (groups x size_classes)
#' @export
#'
zAveOutput = function(x, prop = 0.5){
  ave_x <- colMeans(x[(ceiling(dim(x)[1] - prop*(dim(x)[1])):dim(x)[1]),,], dims = 1)
  return(ave_x)
}

#' Remove Tibble Attributes
#'
#' @title Convert tibble to data frame for efficiency
#' @description Removes tibble attributes and converts to a plain data frame
#'   for improved speed and memory efficiency in computational workflows.
#' @details This utility function strips tibble-specific attributes that can
#'   slow down operations in tight computational loops. Used internally by
#'   ZooMSS for performance optimization when working with large datasets.
#'
#' @param tibble A tibble object to convert
#'
#' @return Plain data frame without tibble attributes
#' @export
#'
untibble <- function (tibble) {
  data.frame(unclass(tibble), check.names = FALSE, stringsAsFactors = FALSE)
}  ## escape the nonsense


#' Convert List to Tibble Format
#'
#' @title Convert ZooMSS list output to tibble with species names
#' @description Converts ZooMSS list output to a tibble format with proper column
#'   names based on functional group species names.
#' @details This function converts aggregated ZooMSS output (typically from
#'   zSumSize or similar functions) into a tibble format suitable for
#'   analysis and visualization. Currently designed for 2D data (species x cells).
#'
#' @param li List of vectors/matrices from ZooMSS aggregation functions
#' @param vmdl ZooMSS model object containing species names (param$Groups$Species)
#'
#' @return Tibble with columns named by species and rows representing spatial cells
#' @export
#'
zConvert2Tibble <- function(li, vmdl){
  df <- tibble::as_tibble(matrix(unlist(li), nrow=length(li), byrow=T), .name_repair = "unique") %>%
    dplyr::rename_with(~vmdl$param$Groups$Species)
  return(df)
}


#' Create Diet Matrix in Long Tibble Format
#'
#' @title Convert diet matrix to long format for analysis and visualization
#' @description Converts ZooMSS diet matrix from wide format to long (tidy) format
#'   with predator-prey relationships clearly defined.
#' @details This function transforms diet matrices into a long format suitable for
#'   analysis and visualization of feeding relationships. The resulting tibble
#'   contains predator-prey pairs with diet fraction values, making it easy to
#'   analyze trophic interactions and create food web visualizations.
#'
#' @param mat Diet matrix from ZooMSS output (predators x prey)
#' @param mdl ZooMSS model object containing species names for labeling
#'
#' @return Long tibble with columns: Predator, Prey, Diet
#' @export
#'
zMakeDietTibble <- function(mat, mdl){
  suppressMessages(
    out <- tibble::as_tibble(mat, .name_repair = "unique") %>%
      dplyr::rename_with(~c("Phyto_Small", "Phyto_Med", "Phyto_Large", mdl$param$Groups$Species)) %>%
      dplyr::mutate(Predator = mdl$param$Groups$Species) %>%
      tidyr::pivot_longer(cols = .data$Phyto_Small:.data$Fish_Large, names_to = "Prey", values_to = "Diet")
  )
  return(out)
}


#' Calculate PPMR Data for Plotting
#'
#' @title Calculate predator-prey mass ratio data for visualization
#' @description Calculates predator-prey mass ratio (PPMR) values and biomass weightings
#'   for creating PPMR distribution plots in ZooMSS analysis.
#' @details This function computes theoretical and realized PPMR patterns by:
#'   - Calculating size-dependent PPMR values using Wirtz 2012 equations
#'   - Weighting by biomass to show community-level patterns
#'   - Computing species-specific PPMR values
#'   - Handling special cases for filter feeders (larvaceans, salps)
#'
#'   This is a helper function primarily used by zPlot_PPMR for visualization.
#'   PPMR analysis provides insights into food web structure and predation patterns.
#'
#' @param dat ZooMSS results object containing abundances and model parameters
#'
#' @return List containing PPMR density data and species-specific values for plotting
#' @export
#'
zExtract_PPMR = function(dat){

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

  # Check for non-finite values and handle edge cases
  total_biom = sum(ave_biom)
  if (!is.finite(total_biom) || total_biom == 0) {
    # If total biomass is zero or non-finite, create uniform weights
    beta_props = matrix(1/length(ave_biom), nrow = nrow(ave_biom), ncol = ncol(ave_biom))
    warning("Non-finite or zero total biomass detected. Using uniform weights for density calculation.")
  } else {
    beta_props = ave_biom/total_biom # Calculate fraction of zoo biomass in each group, in each size class
  }

  # Ensure beta_props values are finite for density function
  beta_props[!is.finite(beta_props)] <- 0

  out <- list()
  out[[1]] <- betas
  out[[2]] <- beta_props
  names(out) <- c("betas", "beta_props")

  temp <- stats::density(betas, weights = beta_props)

  out <- tibble::tibble("x" = temp$x, "y" = temp$y, "mn_beta" = sum(beta_props*betas))

  # Calculate species-specific proportions with safety checks
  row_sums <- rowSums(ave_biom)
  spbeta_props = ave_biom
  for(i in 1:nrow(ave_biom)) {
    if(is.finite(row_sums[i]) && row_sums[i] > 0) {
      spbeta_props[i,] = ave_biom[i,] / row_sums[i]
    } else {
      spbeta_props[i,] = 1/ncol(ave_biom)  # uniform distribution if row sum is invalid
    }
  }
  spbeta_props[!is.finite(spbeta_props)] <- 0  # ensure all values are finite
  spPPMR <- tibble::tibble("Species" = as.factor(dat$model$param$Groups$Species[-which(is.na(dat$model$param$Groups$PPMRscale))]), "Betas" = rowSums(spbeta_props*betas), "y" = NA) # Get species-specific PPMR

  for (s in 1:length(spPPMR$Species)){
    spPPMR$y[s] <- out$y[which.min(abs(out$x - spPPMR$Betas[s]))]
  }

  spPPMR <- spPPMR %>%
    dplyr::mutate(y = .data$y * 0) %>%
    dplyr::bind_rows(spPPMR)

  out2 <- list()
  out2[[1]] <- out
  out2[[2]] <- spPPMR

  return(out2)
}



#' Calculate Phytoplankton Size Spectrum Parameters
#'
#' @title Calculate phytoplankton abundance spectrum from chlorophyll data
#' @description Converts chlorophyll concentration data to phytoplankton size spectrum
#'   parameters (slope, intercept, maximum size) using established oceanographic relationships.
#' @details This function implements the Brewin et al. (2015) algorithm to partition
#'   chlorophyll among picophytoplankton, nanophytoplankton, and microphytoplankton size
#'   classes, then calculates:
#'   - Size spectrum slope and intercept parameters
#'   - Maximum phytoplankton size based on micro proportion
#'   - Biomass estimates for each size class
#'
#'   These parameters drive the dynamic phytoplankton spectrum in ZooMSS that serves
#'   as the base of the food web. The function can work with either chlorophyll-only
#'   data (using empirical relationships) or direct phytoplankton biomass measurements.
#'
#' @param df Data frame containing chlorophyll data (chl column in mg/m^3) and
#'   optionally phytoplankton biomass (phy column in g/m^3)
#'
#' @return Data frame with added columns:
#'   \itemize{
#'     \item phyto_slope: Power law slope for phytoplankton size spectrum
#'     \item phyto_int: Log10 intercept for phytoplankton abundance
#'     \item phyto_max: Maximum phytoplankton size (log10 grams)
#'     \item pico_biom, nano_biom, micro_biom: Biomass in each size class
#'   }
#' @export
#'
#' @references
#' Brewin, R.J.W., et al. (2015). A three-component model of phytoplankton size class
#' for the Atlantic Ocean. Ecological Modelling, 306, 90-101.
#'
#' Maranon, E., et al. (2014). Resource supply overrides temperature as a controlling
#' factor of marine phytoplankton growth. PLoS ONE, 9(6), e99312.
#'
zCalculatePhytoParam <- function(df){ # chl is chlorophyll concentration in mg m^-3, phy is mean euphotic zone phyto in g wet weight m-3

  ## Calculate pico, nano, micro phytoplankton proportions of total chlorophyll
  ## BREWIN ET AL., 2015
  pico <- (0.13*(1-exp(-0.8/0.13*df$chl)))/df$chl
  nano <- (0.77*(1-exp(-0.94/0.77*df$chl)))/df$chl - pico
  micro <- (df$chl - 0.77*(1-exp(-0.94/0.77*df$chl)))/df$chl

  if("phy" %in% colnames(df)){
    tot_biom <- df$phy
  } else {
    ## Convert total chlorophyll to g m^-3 total wet weight - biomass
    ## Allocate total chlorophyll to the three size classes
    c_chl <- ((df$chl^0.89)*(10^1.79))/df$chl # chl:carbon ratio, from Maranon et al. 2014
    tot_biom_c <- c_chl*df$chl/1000 # (convert to grams carbon)
    tot_biom <- tot_biom_c*(1/0.1) # convert to grams wet weight, assuming 0.1 C:ww
  }

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
  max_phyto <- rep(-7, length(df$chl)) # Set -7 to be the max possible size for phyto
  df$phyto_max <- pmin(max_phyto, df$phyto_max)

  return(df)
}


#' Calculate Trophic Levels from Diet Matrix
#'
#' @title Compute trophic levels for functional groups using diet composition
#' @description Calculates trophic levels for each functional group based on their
#'   diet composition using an iterative Gauss-Seidel algorithm.
#' @details This function computes trophic levels by:
#'   - Starting with phytoplankton at trophic level 1.0
#'   - Initializing all other groups at trophic level 2.0
#'   - Iteratively updating trophic levels based on weighted diet composition
#'   - Continuing until convergence (difference < 0.01) or maximum iterations (100)
#'
#'   Trophic level calculation follows: TL = 1 + sum(diet_fraction_i * TL_prey_i)
#'
#'   This provides a quantitative measure of each group's position in the food web
#'   and is useful for analyzing ecosystem structure and energy transfer efficiency.
#'
#' @param diet_matrix 12x15 matrix where rows are predators (functional groups) and
#'   columns are prey (first 3 columns are phytoplankton size classes, remaining 12 are
#'   zooplankton/fish groups). Values represent diet fractions.
#'
#' @return Vector of trophic levels for each functional group (length 12)
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model
#' results <- zoomss_model(input_params, Groups, SaveTimeSteps = FALSE)
#' trophic_levels <- zTrophicLevel(results$diets)
#'
#' # View trophic levels by group
#' names(trophic_levels) <- results$model$param$Groups$Species
#' print(trophic_levels)
#' }
#'
zTrophicLevel <- function(diet_matrix){

  phyto_tl <- 1 # Phyto TL is 1
  start_dynam_tl <- rep(2,12) # Start TL - start at 2 for all zoo and fish groups

  # To be truly generic I need to fix these up with testgroup references #TODO
  curr_phyto_diet <- rowSums(diet_matrix[,1:3]) # Current phyto diet
  curr_dynam_diet <- diet_matrix[,4:15] # Current heterotroph diet

  total_diet <- curr_phyto_diet + rowSums(curr_dynam_diet) # Total consumption, in grams wet weight, by pred group

  curr_phyto_frac <- curr_phyto_diet/total_diet # Fraction of diet from phyto, by pred group and pred sizes
  curr_dynam_frac <- sweep(curr_dynam_diet, 1, total_diet, '/') # Fraction of diet from each prey group for each pred group

  n <- 1
  eps_diff <- 1

  while(eps_diff > 0.01 & n < 100){ # Gauss-Siedel iterative loop to calculate trophic levels, stops when converged or number of loops reaches 100
    n <- n + 1
    eps <- start_dynam_tl[10]

    calc_dynam_tl = sweep(curr_dynam_frac, 2, start_dynam_tl, '*')
    calc_dynam_tl[which(is.nan(calc_dynam_tl) == TRUE)] = 0 # Get rid of nans - these are entrys where there is no biomass for a given group
    #calc_dynam_tl[which(calc_dynam_tl == Inf)] = 0 # Get rid of infinite values, occurs with asymptotic size bins, because there is no biomass to have a diet in those bins
    start_dynam_tl = 1+phyto_tl*curr_phyto_frac + rowSums(calc_dynam_tl) # Update trophic level matrix

    eps_diff <- abs(eps - start_dynam_tl[10])
  } # End Gauss-Siedel loop

  return(start_dynam_tl)

} # End trophic level function
