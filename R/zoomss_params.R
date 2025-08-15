#' Set Up ZooMSS Model Parameters
#'
#' @title Initialize and validate ZooMSS model parameters
#' @description Sets up the complete parameter list for ZooMSS model runs, including 
#'   functional group parameters, model dimensions, and environmental forcing data.
#' @details This function creates a comprehensive parameter object that contains:
#'   
#'   **Static Parameters (fixed across time steps):**
#'   - Model dimensions (number of groups, size classes, time steps)
#'   - Biological parameters (growth efficiency, mortality rates)
#'   - Size class definitions and ranges for each functional group
#'   - Phytoplankton size spectrum parameters
#'   
#'   **Dynamic Parameters (calculated from environmental data):**  
#'   - Phytoplankton abundance time series based on chlorophyll
#'   - Temperature effects on metabolism for zooplankton and fish
#'   - Environmental forcing validation and interpolation
#'   
#'   The function validates that environmental time series data covers the full
#'   simulation period and pre-calculates time-varying parameters to optimize
#'   model performance during the main simulation loop.
#'
#' @param Groups Data frame containing functional group definitions with columns:
#'   Species, Type, W0 (log min size), Wmax (log max size), and various biological parameters
#' @param input_params Data frame with model parameters including:
#'   tmax (max years), dt (time step), isave (save frequency), and environmental time series
#'   (time_step, sst, chlo) for dynamic forcing
#'
#' @return List containing comprehensive model parameters:
#'   \itemize{
#'     \item Groups: Functional group definitions
#'     \item ngrps: Number of functional groups
#'     \item ngrid: Number of size classes
#'     \item w: Size class weights (g)
#'     \item tmax, dt, isave: Temporal parameters
#'     \item zoo_grps, fish_grps: Indices for different organism types
#'     \item phyto_int_ts, phyto_slope_ts: Time series of phytoplankton parameters
#'     \item temp_eff_zoo_ts, temp_eff_fish_ts: Time series of temperature effects
#'     \item Additional biological and physical parameters
#'   }
#' @export
#'
#' @examples
#' \dontrun{
#' # Load functional groups
#' data(Groups)
#' 
#' # Create environmental time series
#' env_data <- zCreateSimpleTimeSeries(1000, 0.01)
#' 
#' # Set up input parameters with environmental forcing
#' input_params <- data.frame(
#'   tmax = 10,
#'   dt = 0.01,
#'   isave = 10, 
#'   time_step = 1:1000,
#'   sst = env_data$sst,
#'   chlo = env_data$chlo
#' )
#' 
#' # Generate parameter list
#' params <- zoomss_params(Groups, input_params)
#' }
#'
zoomss_params <- function(Groups, input_params){

  source("zXtras.R")
  
  param <- list(
    Groups = Groups, # Read in functional group specific parameters from file
    ngrps = dim(Groups)[1], # no. of Groups
    dx = 0.1, # log10 weight step
    day = 12, # day length (hours of each day in sun)
    tmax = input_params$tmax[1], # max years - use scalar value
    dt = input_params$dt[1], # timestep - use scalar value
    w0 = 10^(min(Groups$W0)),		# minimum size class
    wMax = 10^(max(Groups$Wmax)),# maximum size class
    gge_base = 0.25, # baseline gross growth efficiency
    ZSpre = 1, # senescence mortality prefactor
    ZSexp = 0.3, # senescence mortality exponent
    w0_phyto = 10^(-14.5), # minimum phytoplankton size class (1um)
    # wMax_phyto will be calculated from time series
    zoo_grps = which(Groups$Type == "Zooplankton"), # Which rows are zooplankton
    fish_grps = which(Groups$Type == "Fish"), # Which rows are fish
    num_zoo = sum(Groups$Type == "Zooplankton"), # How many zooplankton
    num_fish = sum(Groups$Type == "Fish"), # How many fish
    cc_phyto = 0.1, # Carbon content of phytoplankton size classes
    isave = input_params$isave[1] # how often to save results every 'isave' time steps - use scalar value
  )

  ## Add additional parameters which are based on the parameter set
  param2 <- list(
    nsave  = max(1, ceiling(param$tmax/(param$dt*param$isave))), # Number of time slots to save - use ceiling to ensure enough space
    ntime = ceiling(param$tmax / param$dt) # Total number of time steps
  )

  # Process provided time series
  n_time_steps <- nrow(input_params)
  
  # Validate that environmental data covers the full simulation
  required_time_steps <- ceiling(param$tmax / param$dt)
  
  # Ensure these are scalar values for comparison
  n_time_steps_scalar <- n_time_steps[1]
  required_time_steps_scalar <- required_time_steps[1]
  
  if (n_time_steps_scalar < required_time_steps_scalar) {
    stop("Environmental time series too short! Need ", required_time_steps_scalar, 
         " timesteps (", param$tmax, " years with dt=", param$dt, 
         ") but only have ", n_time_steps_scalar, " timesteps (",
         round(n_time_steps_scalar * param$dt, 2), " years)")
  }
  
  cat("✅ Environmental validation: ", n_time_steps_scalar, " timesteps cover ", 
      round(n_time_steps_scalar * param$dt[1], 1), " years (need ", param$tmax[1], " years)\n")
  
  # Pre-calculate phytoplankton parameters for each time step
  
  # Initialize arrays to store time-varying parameters (optimize memory allocation)
  param2$phyto_int_ts <- numeric(n_time_steps)
  param2$phyto_slope_ts <- numeric(n_time_steps)
  param2$phyto_max_ts <- numeric(n_time_steps)
  param2$wMax_phyto_ts <- numeric(n_time_steps)
  
  # Pre-calculate temperature effects for each time step (more efficient than in-loop)
  param2$temp_eff_zoo_ts <- matrix(NA, nrow = n_time_steps, ncol = param$num_zoo)
  param2$temp_eff_fish_ts <- matrix(NA, nrow = n_time_steps, ncol = param$num_fish)
  
  # Check if phytoplankton parameters are already in input_params (expanded format)
  if (nrow(input_params) == n_time_steps && all(c("phyto_int", "phyto_slope", "phyto_max") %in% names(input_params))) {
    cat("✅ Using pre-calculated phytoplankton parameters from input_params\n")
    
    # Use pre-calculated values (most efficient path)
    param2$phyto_int_ts <- input_params$phyto_int
    param2$phyto_slope_ts <- input_params$phyto_slope
    param2$phyto_max_ts <- input_params$phyto_max
    param2$wMax_phyto_ts <- 10^input_params$phyto_max
    
    # Calculate temperature effects for each time step - consistent temperature formula throughout model
    for (i in 1:n_time_steps) {
      sst_i <- input_params$sst[i]
      temp_factor <- 2^((sst_i - 30)/10)
      param2$temp_eff_zoo_ts[i, ] <- temp_factor
      param2$temp_eff_fish_ts[i, ] <- temp_factor
    }
    
  } else {
    cat("⚠️  Calculating phytoplankton parameters from environmental time series\n")
    
    # Calculate phytoplankton parameters for each time step (fallback method)
    for (i in 1:n_time_steps) {
      sst_i <- input_params$sst[i]
      chlo_i <- input_params$chlo[i]
      
      # Calculate phytoplankton parameters using the existing function
      temp_df <- data.frame(cellID = 1, sst = sst_i, chlo = chlo_i)
      phyto_params <- zCalculatePhytoParam(temp_df)
      
      param2$phyto_int_ts[i] <- phyto_params$phyto_int
      param2$phyto_slope_ts[i] <- phyto_params$phyto_slope
      param2$phyto_max_ts[i] <- phyto_params$phyto_max
      param2$wMax_phyto_ts[i] <- 10^phyto_params$phyto_max
      
      # Temperature effects (same calculation throughout model)
      temp_factor <- 2^((sst_i - 30)/10)
      param2$temp_eff_zoo_ts[i, ] <- temp_factor
      param2$temp_eff_fish_ts[i, ] <- temp_factor
    }
  }
  
  # Store the maximum wMax_phyto from time series for grid creation
  param2$wMax_phyto <- max(param2$wMax_phyto_ts)
  
  # Add final parameters that depend on the complete parameter set (calculate only once)
  param2$w_log10 <- round(seq(from = min(Groups$W0), to = max(Groups$Wmax), param$dx), digits = 2) # Set up log10 weight grid
  param2$w <- 10^(seq(from = min(Groups$W0), to = max(Groups$Wmax), param$dx)) # Set up weight grid
  param2$w_phyto <- 10^(seq(from = log10(param$w0_phyto), to = log10(param2$wMax_phyto), param$dx)) # Set up phytoplankton size classes
  param2$ngrid <- length(param2$w) # total number of size classes for zoo and fish
  param2$ngridPP <- length(param2$w_phyto) # total number of size classes for phyto

  # Final parameter combination
  # Exclude time series vectors from input_params since they're now stored as _ts arrays in param2
  input_params_filtered <- input_params[!names(input_params) %in% c("tmax", "dt", "isave", "time_step", "sst", "chlo", "phyto_int", "phyto_slope", "phyto_max")]
  
  param_final <- c(input_params_filtered, param, param2)
  return(param_final)
}
