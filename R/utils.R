#' Sum ZooMSS Output Across Size Bins
#'
#' @title Aggregate ZooMSS abundances across all size classes
#' @description Sums abundance values across all size classes for each functional group,
#'   providing total abundance per group.
#' @details This function collapses the size dimension of ZooMSS output by summing
#'   across all size classes. Useful for analyzing total abundance patterns without
#'   size structure detail.
#' @param x 3D array outptut from ZooMSS model
#' @param method Character string specifying aggregation method: "sum" (default) or "mean".
#' @return List of vectors with total abundance per functional group
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model
#' results <- zoomss_model(input_params, Groups)
#' total_abundances <- reduceSize(results$abundances)
#' }
#'
reduceSize = function(x, method = "sum") {

  assertthat::assert_that(is.array(x))
  assertthat::assert_that(method %in% c("sum", "mean"))

  apply(x, c(1, 2), match.fun(method)) # Sum or mean ZooMSS output across the size bins
}


#' Aggregate ZooMSS abundances across all species
#'
#' @title Aggregate ZooMSS abundances across all species
#' @description Aggregates abundance values across all species bins for each functional group and size class using the specified method.
#' @details This function collapses the species dimension by applying the specified method (sum or mean) across all species bins.
#' @param x 3D array outptut from ZooMSS model
#' @param method Character string specifying aggregation method: "sum" (default) or "mean".
#' @return Array with species dimension reduced using the specified method
#' @export
reduceSpecies = function(x, method = "sum") {
  assertthat::assert_that(is.array(x))
  assertthat::assert_that(method %in% c("sum", "mean"))

  apply(x, c(1, 3), match.fun(method))
}


#' Aggregate abundances across all groups and size classes
#'
#' @title Aggregate abundances across all groups and size classes
#' @description Calculates total abundance across all functional groups and size classes using the specified method.
#' @details This function provides the most aggregated view of ZooMSS output by applying the method across both functional groups and size classes.
#' @param x 3D array outptut from ZooMSS model
#' @param method Character string specifying aggregation method: "sum" (default) or "mean".
#' @return Vector of total abundance values (one per spatial cell)
#' @export
reduceAll = function(x, method = "sum") {
  assertthat::assert_that(is.array(x))
  assertthat::assert_that(method %in% c("sum", "mean"))

  apply(x, 1, match.fun(method))
}

#' Convert Abundance to Biomass
#'
#' @title Convert ZooMSS abundance matrices to biomass by multiplying by body weights
#' @description Converts abundance data to wet weight biomass by multiplying abundances
#'   by the corresponding body weights for each size class. Optionally converts to carbon biomass.
#' @details This function transforms abundance matrices to biomass by applying the
#'   weight vector across size classes. Essential for analyses requiring biomass
#'   units rather than abundance counts. Works with 3D arrays (time, groups, size_classes).
#'   Can convert to either wet weight or carbon biomass units.
#'
#' @param mdl ZooMSS model object containing abundance array (N) and weight vector (param$w)
#' @param units Character string specifying biomass units: "ww" for wet weight (default) or "carbon" for carbon biomass
#'
#' @return 3D array of biomass values with same dimensions as N. Units depend on the units parameter:
#'   \itemize{
#'     \item "ww": grams wet weight
#'     \item "carbon": grams carbon
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' # Run ZooMSS model
#' results <- zoomss_model(input_params, Groups)
#'
#' # Convert abundances to wet weight biomass
#' biomass_ww <- getBiomass(results, units = "ww")
#'
#' # Convert abundances to carbon biomass
#' biomass_carbon <- getBiomass(results, units = "carbon")
#' }
#'
getBiomass <- function(mdl, units = "ww") {

  # Validate units parameter
  if (!units %in% c("ww", "carbon")) {
    stop("units must be either 'ww' (wet weight) or 'carbon'")
  }

  # Check that N exists in the model
  if (!"abundance" %in% names(mdl)) {
    stop("Abundance array 'N' not found in model output")
  }

  # Check that weights exist in model parameters
  if (is.null(mdl$param$w)) {
    stop("Weight vector 'w' not found in mdl$param$w")
  }

  # Get abundance array and weights
  N <- mdl$abundance
  w <- mdl$param$w

  # Check dimensions match
  if (dim(N)[3] != length(w)) {
    stop("Size dimension of abundance (", dim(N)[3], ") does not match length of weight vector (", length(w), ")")
  }

  # Convert abundance to wet weight biomass by multiplying by weights across size dimension (3rd dimension)
  Biomass <- sweep(N, 3, w, '*')  # Multiply each size class by its corresponding weight

  # Convert to carbon biomass if requested
  if (units == "carbon") {
    # Check that carbon content factors exist
    if (is.null(mdl$param$Groups$Carbon)) {
      stop("Carbon content factors 'Groups$Carbon' not found in mdl$param$Groups$Carbon")
    }

    # Check dimensions match
    if (dim(N)[2] != length(mdl$param$Groups$Carbon)) {
      stop("Group dimension of N (", dim(N)[2], ") does not match length of carbon vector (", length(mdl$param$Groups$Carbon), ")")
    }

    # Convert to carbon biomass by multiplying by carbon content across group dimension (2nd dimension)
    Biomass <- sweep(Biomass, 2, mdl$param$Groups$Carbon, '*')
  }

  return(Biomass)
}


#' Extract Size Range from ZooMSS Output
#'
#' @title Extract specific size class range from model variable
#' @description Subsets ZooMSS model output to include only specified size range,
#'   useful for focusing analysis on particular size ranges.
#' @details This function extracts a subset of size classes from the specified
#'   ZooMSS model variable. Useful for analyzing specific size ranges
#'   (e.g., microzooplankton, mesozooplankton) or excluding boundary effects
#'   from model analysis. The function converts log10 size values to size class
#'   indices automatically.
#'
#' @param mdl ZooMSS model results object containing model parameters and output arrays
#' @param var Character string specifying which variable to extract (e.g., "N", "Z", "Growth")
#' @param min_size Minimum size (log10 grams) to extract
#' @param max_size Maximum size (log10 grams) to extract
#'
#' @return Array with same dimensions as original variable but subsetted size range
#' @export
#'
#' @examples
#' \dontrun{
#' # Run ZooMSS model
#' results <- zoomss_model(input_params, Groups)
#'
#' # Extract mesozooplankton size range from abundance data
#' meso_abundance <- extractSizeRange(results, "N", min_size = -8, max_size = -5)
#'
#' # Extract microzooplankton size range from growth data
#' micro_growth <- extractSizeRange(results, "Growth", min_size = -10, max_size = -8)
#' }
#'
extractSizeRange = function(mdl, var, min_size, max_size) {

  # Extract the specified variable from model output
  if (!var %in% names(mdl)) {
    stop("Variable '", var, "' not found in model output. Available variables: ",
         paste(names(mdl), collapse = ", "))
  }
  x <- mdl[[var]]

  # Get size grid from model parameters
  if (is.null(mdl$param$w_log10)) {
    stop("Cannot find size grid (w_log10) in mdl$param$w_log10")
  }
  w_log10 <- mdl$param$w_log10

  # Find indices corresponding to size range
  min_idx <- which.min(abs(w_log10 - min_size))
  max_idx <- which.min(abs(w_log10 - max_size))

  # Ensure min_idx <= max_idx
  if (min_idx > max_idx) {
    temp <- min_idx
    min_idx <- max_idx
    max_idx <- temp
  }

  # Validate indices are within bounds
  n_size_classes <- length(w_log10)
  min_idx <- max(1, min_idx)
  max_idx <- min(n_size_classes, max_idx)

  # Provide feedback about the extraction
  actual_min_size <- w_log10[min_idx]
  actual_max_size <- w_log10[max_idx]
  n_extracted <- max_idx - min_idx + 1
  cat("Extracting size range", actual_min_size, "to", actual_max_size,
      "log10(g) (", n_extracted, "size classes) from variable", var, "\n")

  # Extract the size range from the variable
  # Handle different array dimensions
  if (length(dim(x)) == 2) {
    # 2D array: (groups, size_classes)
    out <- x[, min_idx:max_idx, drop = FALSE]
  } else if (length(dim(x)) == 3) {
    # 3D array: (time, groups, size_classes)
    out <- x[, , min_idx:max_idx, drop = FALSE]
  } else {
    stop("Variable '", var, "' must be a 2D or 3D array")
  }

  return(out)
}


#' Calculate Average Output from Model Time Series
#'
#' @title Calculate mean of final portion of ZooMSS time series
#' @description Calculates the mean of the final n years of a time series
#'   to obtain equilibrium values after model spin-up period.
#' @details This function removes the initial transient period from time series data
#'   and calculates the mean of the final n years, providing representative
#'   steady-state values. Essential for obtaining equilibrium abundances, growth rates,
#'   and other model outputs after the model has reached dynamic equilibrium.
#'
#' @param mdl ZooMSS model results object containing model parameters and output arrays
#' @param var Character string specifying which variable to extract and average (e.g., "N", "Growth", "Mort")
#' @param n_years Number of years from the end of the time series to average (default: 10)
#'
#' @return 2D array with averaged values (groups x size_classes)
#' @export
#'
#' @examples
#' \dontrun{
#' # Run ZooMSS model
#' results <- zoomss_model(input_params, Groups)
#'
#' # Average final 3 years of abundance data
#' avg_abundance <- averageTimeSeries(results, "N", n_years = 3)
#'
#' # Average final 10 years of growth data (default)
#' avg_growth <- averageTimeSeries(results, "gg")
#' }
#'
averageTimeSeries = function(mdl, var, n_years = 10){

  # Extract dt from model parameters
  if (is.null(mdl$param$dt)) {
    stop("Cannot find dt (time step) in mdl$param$dt")
  }
  dt <- mdl$param$dt

  # Extract isave parameter (how often results were saved)
  if (is.null(mdl$param$isave)) {
    stop("Cannot find isave (save interval) in mdl$param$isave")
  }
  isave <- mdl$param$isave

  # Extract the specified variable from model output
  if (!var %in% names(mdl)) {
    stop("Variable '", var, "' not found in model output. Available variables: ",
         "N", "Z", "gg")
  }
  x <- mdl[[var]]

  # Check that we have a 3D array
  if (length(dim(x)) != 3) {
    stop("Variable '", var, "' must be a 3D array with dimensions (time, groups, size_classes)")
  }

  # Calculate number of saved time steps corresponding to n_years
  # Each saved time step represents isave * dt years
  dt_saved <- isave * dt  # Years represented by each saved time step
  n_time_steps <- round(n_years / dt_saved)

  # Get total number of time steps in the array
  total_time_steps <- dim(x)[1]

  # Ensure we don't try to average more time steps than available
  n_time_steps <- min(n_time_steps, total_time_steps)

  # Calculate start index for averaging (from end)
  start_idx <- max(1, total_time_steps - n_time_steps + 1)

  # Provide helpful feedback about what's being averaged
  actual_years <- n_time_steps * dt_saved
  cat("Averaging final ", actual_years, " years (", n_time_steps, " saved time steps with isave = ", isave, ") of ",
      var, " from ", total_time_steps, " total saved time steps.\n", sep = "")

  # Calculate the average
  ave_x <- colMeans(x[start_idx:total_time_steps,,], dims = 1)
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
#' @param tibble A tibble or data frame object to convert
#'
#' @return Plain data frame without tibble attributes
#'
#' @noRd
#'
untibble <- function (tibble) {
  data.frame(unclass(tibble), check.names = FALSE, stringsAsFactors = FALSE)
}  ## escape the nonsense




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
#'   - Processing either time-averaged abundances (2D) or full time series (3D)
#'
#'   The function dynamically determines size class ranges for larvaceans and salps
#'   based on their W0 and Wmax values. For 3D abundance arrays, it calculates
#'   PPMR for each time step separately.
#'
#'   This is a helper function primarily used by plotPPMR for visualization.
#'   PPMR analysis provides insights into food web structure and predation patterns.
#'
#' @param mdl ZooMSS results object containing abundance data (N) and model parameters
#'
#' @return For 2D input: List containing PPMR density data and species-specific values for plotting
#'         For 3D input: Array where first dimension is time, containing PPMR results for each timestep
#' @export
#'
extractPPMR = function(mdl){

  # min_size = min(mdl$param$Groups$W0) # smallest size class
  # max_size = max(mdl$param$Groups$Wmax) # largest size class

  w = mdl$param$w # all size classes
  w_log10 = log10(mdl$param$w) # log10 size classes

  # Calculate PPMR (beta) table, where dim1 = group, dim2 = body size with
  # value being PPMR for that body size (this is not realised PPMR - not
  # emergent from diet but calculated from m-values and Wirtz, 2012 equation)
  D.z = 2*(3*(w)*1e12/(4*pi))^(1/3) # convert body mass g to ESD (um)
  zoo_m = mdl$param$Groups$PPMRscale # pull out PPMR scaling values from parameter table
  betas =  log10(t(sapply(zoo_m, function(x){(exp(0.02*log(D.z)^2 - x + 1.832))^3}))) # Convert m to betas, using Wirtz 2012 equation
  betas = betas[-which(is.na(mdl$param$Groups$PPMRscale)),] # remove fish rows

  ## Dynamically determine size class indices for special groups
  # Find larvaceans indices
  larv_idx <- which(mdl$param$Groups$Species=="Larvaceans")
  if(length(larv_idx) > 0 && !is.na(mdl$param$Groups$PPMRscale[larv_idx])) {
    larv_w0 <- mdl$param$Groups$W0[larv_idx]
    larv_wmax <- mdl$param$Groups$Wmax[larv_idx]
    larv_start_idx <- which.min(abs(w_log10 - larv_w0))
    larv_end_idx <- which.min(abs(w_log10 - larv_wmax))
    larv_beta_idx <- which(mdl$param$Groups$Species[-which(is.na(mdl$param$Groups$PPMRscale))]=="Larvaceans")

    if(length(larv_beta_idx) > 0 && larv_start_idx <= larv_end_idx) {
      # Set PPMR increases by 0.1 for each 0.1 log10 size interval
      n_larv_sizes <- larv_end_idx - larv_start_idx + 1
      betas[larv_beta_idx, larv_start_idx:larv_end_idx] <- betas[larv_beta_idx, larv_start_idx] +
        seq(0, (n_larv_sizes-1)*0.1, 0.1)
    }
  }

  # Find salps indices
  salp_idx <- which(mdl$param$Groups$Species=="Salps")
  if(length(salp_idx) > 0 && !is.na(mdl$param$Groups$PPMRscale[salp_idx])) {
    salp_w0 <- mdl$param$Groups$W0[salp_idx]
    salp_wmax <- mdl$param$Groups$Wmax[salp_idx]
    salp_start_idx <- which.min(abs(w_log10 - salp_w0))
    salp_end_idx <- which.min(abs(w_log10 - salp_wmax))
    salp_beta_idx <- which(mdl$param$Groups$Species[-which(is.na(mdl$param$Groups$PPMRscale))]=="Salps")

    if(length(salp_beta_idx) > 0 && salp_start_idx <= salp_end_idx) {
      # Set PPMR increases by 0.1 for each 0.1 log10 size interval
      n_salp_sizes <- salp_end_idx - salp_start_idx + 1
      betas[salp_beta_idx, salp_start_idx:salp_end_idx] <- betas[salp_beta_idx, salp_start_idx] +
        seq(0, (n_salp_sizes-1)*0.1, 0.1)
    }
  }

  # Handle 3D abundance array (time, groups, size_classes)
  n_timesteps <- dim(mdl$abundance)[1]

  # Initialize output array
  results_array <- array(list(), dim = c(n_timesteps))

  for(t in 1:n_timesteps) {
    # Extract abundance for this timestep
    N_t <- mdl$abundance[t, , ]

    # Calculate biomass for this timestep
    ave_biom = sweep(N_t, 2, w, "*") # Calculate biomass for each group and size
    ave_biom = ave_biom[-which(is.na(mdl$param$Groups$PPMRscale)),] # remove rows for fish

    # Check for non-finite values and handle edge cases
    total_biom = sum(ave_biom)
    if (!is.finite(total_biom) || total_biom == 0) {
      # If total biomass is zero or non-finite, create uniform weights
      beta_props = matrix(1/length(ave_biom), nrow = nrow(ave_biom), ncol = ncol(ave_biom))
      warning("Non-finite or zero total biomass detected at timestep ", t, ". Using uniform weights for density calculation.")
    } else {
      beta_props = ave_biom/total_biom # Calculate fraction of zoo biomass in each group, in each size class
    }

    # Ensure beta_props values are finite for density function
    beta_props[!is.finite(beta_props)] <- 0

    # Calculate density with bandwidth selection using weights
    betas_vec <- as.vector(betas)
    beta_props_vec <- as.vector(beta_props)
    non_zero_weights <- beta_props_vec > 0
    unique_betas <- unique(betas_vec[non_zero_weights])

    if (length(unique_betas) < 2) {
      # Not enough unique values for automatic bandwidth selection
      bw <- if (length(unique_betas) == 1) 0.1 else diff(range(unique_betas))/10
      temp <- suppressWarnings(stats::density(betas_vec, weights = beta_props_vec, bw = bw))
    } else {
      temp <- suppressWarnings(stats::density(betas_vec, weights = beta_props_vec))
    }
    out <- tibble::tibble("x" = temp$x, "y" = temp$y, "mn_beta" = sum(beta_props*betas))

    # Calculate species-specific proportions with safety checks
    row_sums <- rowSums(ave_biom)
    spbeta_props = ave_biom
    for(i in seq_len(nrow(ave_biom))) {
      if(is.finite(row_sums[i]) && row_sums[i] > 0) {
        spbeta_props[i,] = ave_biom[i,] / row_sums[i]
      } else {
        spbeta_props[i,] = 1/ncol(ave_biom)  # uniform distribution if row sum is invalid
      }
    }
    spbeta_props[!is.finite(spbeta_props)] <- 0  # ensure all values are finite
    spPPMR <- tibble::tibble("Species" = as.factor(mdl$param$Groups$Species[-which(is.na(mdl$param$Groups$PPMRscale))]), "Betas" = rowSums(spbeta_props*betas), "y" = NA) # Get species-specific PPMR

    for (s in seq_along(spPPMR$Species)){
      spPPMR$y[s] <- out$y[which.min(abs(out$x - spPPMR$Betas[s]))]
    }

    spPPMR <- spPPMR %>%
      dplyr::mutate(y = .data$y * 0) %>%
      dplyr::bind_rows(spPPMR)

    # Store results for this timestep
    results_array[[t]] <- list("ppmr_density" = out, "species_ppmr" = spPPMR)
  }

  return(results_array)
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
#' @param dat Data frame containing chlorophyll data (chl column in mg/m^3) and
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
calculatePhytoParam <- function(dat){

  ## Calculate pico, nano, micro phytoplankton proportions of total chlorophyll
  ## BREWIN ET AL., 2015
  pico <- (0.13*(1-exp(-0.8/0.13*dat$chl)))/dat$chl
  nano <- (0.77*(1-exp(-0.94/0.77*dat$chl)))/dat$chl - pico
  micro <- (dat$chl - 0.77*(1-exp(-0.94/0.77*dat$chl)))/dat$chl

  if("phy" %in% colnames(dat)){
    tot_biom <- dat$phy
  } else {
    ## Convert total chlorophyll to g m^-3 total wet weight - biomass
    ## Allocate total chlorophyll to the three size classes
    c_chl <- ((dat$chl^0.89)*(10^1.79))/dat$chl # chl:carbon ratio, from Maranon et al. 2014
    tot_biom_c <- c_chl*dat$chl/1000 # (convert to grams carbon)
    tot_biom <- tot_biom_c*(1/0.1) # convert to grams wet weight, assuming 0.1 C:ww
  }

  # Break up total biom into pico, nano and micro
  dat$pico_biom <- pico*tot_biom
  dat$nano_biom <- nano*tot_biom
  dat$micro_biom <- micro*tot_biom

  ## Find abundances at boundaries of pico, nano size ranges, by analytically
  ## solving integral of N = aw^b

  w_0 <- -14.5 # log minimum size of picophytoplankton
  w_1 <- -11.5 # log minimum size of nanophytoplankton (max size of pico also)
  w_2 <- -8.5 # log minimum size of macrophytoplankton (max size of nano also)

  dat$phyto_slope <- (log10(dat$pico_biom) - log10(dat$nano_biom) - w_1 + w_2)/(w_1 - w_2)  # Calculate slope
  dat$phyto_int <- log10(dat$pico_biom*(dat$phyto_slope+1)/((10^(w_1))^(dat$phyto_slope+1) - (10^(w_0))^(dat$phyto_slope+1))) # Calculate intercept

  ## Calculate maximum size
  dat$phyto_max <- 0.1*round((-8.4 + 2*micro)/0.1) # Maximum size depends on the proportion of micro
  max_phyto <- rep(-7, length(dat$chl)) # Set -7 to be the max possible size for phyto
  dat$phyto_max <- pmin(max_phyto, dat$phyto_max)

  return(dat)
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
#'   - Processing 3D diet arrays with time series data
#'
#'   Trophic level calculation follows: TL = 1 + sum(diet_fraction_i * TL_prey_i)
#'
#'   The function calculates trophic levels for each time step separately and
#'   dynamically determines the number of groups from the diet matrix dimensions.
#'
#'   This provides a quantitative measure of each group's position in the food web
#'   and is useful for analyzing ecosystem structure and energy transfer efficiency.
#'
#' @param mdl ZooMSS model results object containing 3D diet data (mdl$diet).
#'   Dimensions are time, groups, prey_items where columns 1:3 are always
#'   phytoplankton size classes and remaining columns are zooplankton/fish groups.
#'
#' @return Matrix where rows are time steps and columns are functional groups
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model with 3D time series
#' results <- zoomss_model(input_params, Groups)
#' trophic_levels <- extractTrophicLevels(results)  # Returns matrix (time x groups)
#'
#' # View trophic levels by group for final time step
#' final_tl <- trophic_levels[nrow(trophic_levels), ]
#' names(final_tl) <- results$param$Groups$Species
#' print(final_tl)
#' }
#'
extractTrophicLevels <- function(mdl){

  phyto_tl <- 1 # Phyto TL is 1

  # Extract 3D diet array (time, groups, prey)
  if(!"diet" %in% names(mdl)) {
    stop("Diet data (mdl$diet) not found in model output")
  }

  diet_array <- mdl$diet

  if(length(dim(diet_array)) != 3) {
    stop("Diet array must be 3D with dimensions [time, groups, prey_items]")
  }

  n_timesteps <- dim(diet_array)[1]
  n_groups <- dim(diet_array)[2]
  n_prey_items <- dim(diet_array)[3]

  # Dynamically determine number of zooplankton/fish groups
  n_dynamic_groups <- n_prey_items - 3  # Subtract 3 phytoplankton columns

  # Initialize output matrix (time x groups)
  trophic_levels <- matrix(NA, nrow = n_timesteps, ncol = n_groups)

  for(t in 1:n_timesteps) {
    # Extract diet matrix for this timestep
    diet_matrix <- diet_array[t,,]

    # Calculate trophic levels for this timestep
    start_dynam_tl <- rep(2, n_groups) # Start TL at 2 for all zoo and fish groups

    # Current phyto diet (columns 1:3)
    curr_phyto_diet <- rowSums(diet_matrix[,1:3])

    # Current heterotroph diet (columns 4 onwards)
    if(n_dynamic_groups > 0) {
      curr_dynam_diet <- diet_matrix[,4:(3+n_dynamic_groups)]
      if(n_dynamic_groups == 1) {
        # Handle case where there's only one dynamic group (ensure it's a matrix)
        curr_dynam_diet <- matrix(curr_dynam_diet, ncol = 1)
      }
    } else {
      # If no dynamic groups, create empty matrix
      curr_dynam_diet <- matrix(0, nrow = n_groups, ncol = 0)
    }

    # Calculate total diet
    total_diet <- curr_phyto_diet + rowSums(curr_dynam_diet)

    # Avoid division by zero
    valid_diet <- total_diet > 0

    if(any(valid_diet)) {
      # Calculate diet fractions only for groups with valid diets
      curr_phyto_frac <- rep(0, n_groups)
      curr_phyto_frac[valid_diet] <- curr_phyto_diet[valid_diet] / total_diet[valid_diet]

      if(n_dynamic_groups > 0) {
        curr_dynam_frac <- matrix(0, nrow = n_groups, ncol = n_dynamic_groups)
        curr_dynam_frac[valid_diet,] <- sweep(curr_dynam_diet[valid_diet,, drop = FALSE], 1, total_diet[valid_diet], '/')
      } else {
        curr_dynam_frac <- matrix(0, nrow = n_groups, ncol = 0)
      }

      # Gauss-Seidel iterative loop
      n_iter <- 1
      eps_diff <- 1

      while(eps_diff > 0.01 && n_iter < 100) {
        n_iter <- n_iter + 1
        eps <- start_dynam_tl[min(n_groups, 10)]  # Use a representative group for convergence check

        if(n_dynamic_groups > 0) {
          calc_dynam_tl <- sweep(curr_dynam_frac, 2, start_dynam_tl[1:min(n_dynamic_groups, n_groups)], '*')
          calc_dynam_tl[is.nan(calc_dynam_tl)] <- 0 # Remove NaNs
          start_dynam_tl <- 1 + phyto_tl * curr_phyto_frac + rowSums(calc_dynam_tl)
        } else {
          start_dynam_tl <- 1 + phyto_tl * curr_phyto_frac
        }

        eps_diff <- abs(eps - start_dynam_tl[min(n_groups, 10)])
      }
    }

    trophic_levels[t,] <- start_dynam_tl
  }

  return(trophic_levels)
}
