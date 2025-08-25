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
#'   time (time vector in years), sst (sea surface temperature), and chl (chlorophyll).
#'   The time vector can start at any value and the model automatically calculates dt (time step) and tmax (maximum time).
#' @param isave Save frequency in time steps (default: 50)
#'
#' @return List containing comprehensive model parameters:
#'   \itemize{
#'     \item Groups: Functional group definitions
#'     \item ngrps: Number of functional groups
#'     \item ngrid: Number of size classes
#'     \item w: Size class weights (g)
#'     \item tmax, dt, isave: Temporal parameters
#'     \item zoo_grps, fish_grps: Indices for different organism types
#'     \item phyto_int, phyto_slope: Time series of phytoplankton parameters
#'     \item temp_eff_zoo, temp_eff_fish: Time series of temperature effects
#'     \item Additional biological and physical parameters
#'   }
#'
#' @examples
#' \dontrun{
#' # Load functional groups
#' data(Groups)
#'
#' # Create environmental time series
#' env_data <- createEnviroData(10, 0.01)
#' input_params <- createInputParams(env_data$time, env_data$sst, env_data$chl)
#'
#' # Generate parameter list
#' params <- zoomss_params(Groups, input_params, isave = 50)
#' }
#'
#' @noRd
#'
zoomss_params <- function(Groups, input_params, isave){

  # Calculate dt and tmax from time column in input_params
  time_values <- input_params$time
  dt_calc <- time_values[2] - time_values[1]  # Use first time step difference
  tmax_calc <- max(time_values)  # Maximum time value (not duration)

  # Check if time steps are uniform - ERROR if not consistent
  if (length(time_values) > 2) {
    dt_diffs <- diff(time_values)
    max_diff <- max(abs(dt_diffs - dt_calc))
    if (max_diff > dt_calc * 0.001) {  # Allow only 0.1% variation
      stop("Time steps are not uniform in input_params. Maximum deviation: ",
           round(max_diff, 6), " (", round(100 * max_diff / dt_calc, 2), "% of dt). ",
           "ZooMSS requires uniform time steps for accurate results.")
    }
  }

  param <- list(
    Groups = Groups, # Read in functional group specific parameters from file
    ngrps = dim(Groups)[1], # no. of Groups
    dx = 0.1, # log10 weight step
    day = 12, # day length (hours of each day in sun)
    tmax = tmax_calc, # max years - calculated from time
    dt = dt_calc, # timestep - calculated from time
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
    isave = isave # how often to save results every 'isave' time steps
  )

  # Process provided time series
  n_time_steps <- nrow(input_params)

  ## Add additional parameters which are based on the parameter set
  param2 <- list(
    nsave  = max(1, floor(n_time_steps / param$isave)), # Number of time slots to save (regular intervals only, minimum 1)
    ntime = n_time_steps, # Total number of time steps to simulate
    itimemax = n_time_steps # Number of iterations for time loop (one for each time point)
  )
  # Pre-calculate phytoplankton parameters for each time step

  # Initialize arrays to store time-varying parameters (optimize memory allocation)
  param2$phyto_int <- numeric(n_time_steps)
  param2$phyto_slope <- numeric(n_time_steps)
  param2$phyto_max <- numeric(n_time_steps)
  param2$wMax_phyto <- numeric(n_time_steps)

  # Pre-calculate temperature effects for each time step (more efficient than in-loop)
  param2$temp_eff_zoo <- matrix(NA, nrow = n_time_steps, ncol = param$num_zoo)
  param2$temp_eff_fish <- matrix(NA, nrow = n_time_steps, ncol = param$num_fish)

  # Check if phytoplankton parameters are already in input_params (expanded format)
  if (nrow(input_params) == n_time_steps && all(c("phyto_int", "phyto_slope", "phyto_max") %in% names(input_params))) {
    cat("Using pre-calculated phytoplankton parameters from input_params\n")

    # Use pre-calculated values (most efficient path)
    param2$phyto_int <- input_params$phyto_int
    param2$phyto_slope <- input_params$phyto_slope
    param2$phyto_max <- input_params$phyto_max
    param2$wMax_phyto <- 10^input_params$phyto_max

    # Calculate temperature effects for each time step - consistent temperature formula throughout model
    temp_factor <- 2^((input_params$sst - 30)/10)
    param2$temp_eff_zoo[,] <- temp_factor # preallocated above for all groups
    param2$temp_eff_fish[,] <- temp_factor # preallocated above for all groups

  } else {
    cat("Calculating phytoplankton parameters from environmental time series\n")

    # Calculate phytoplankton parameters for each time step

    # Calculate phytoplankton parameters using the existing function
    phyto_params <- calculatePhytoParam(input_params)

    param2$phyto_int <- phyto_params$phyto_int
    param2$phyto_slope <- phyto_params$phyto_slope
    param2$phyto_max <- phyto_params$phyto_max
    param2$wMax_phyto <- 10^phyto_params$phyto_max

    # Temperature effects (same calculation throughout model)
    temp_factor <- 2^((input_params$sst - 30)/10)
    param2$temp_eff_zoo[,] <- temp_factor # preallocated above for all groups
    param2$temp_eff_fish[,] <- temp_factor # preallocated above for all groups
  }

  # Store the maximum wMax_phyto from time series for grid creation
  param2$wMax_phyto <- max(param2$wMax_phyto)

  # Add final parameters that depend on the complete parameter set (calculate only once)
  param2$w_log10 <- round(seq(from = min(Groups$W0), to = max(Groups$Wmax), param$dx), digits = 2) # Set up log10 weight grid
  param2$w <- 10^(seq(from = min(Groups$W0), to = max(Groups$Wmax), param$dx)) # Set up weight grid
  param2$w_phyto <- 10^(seq(from = log10(param$w0_phyto), to = log10(param2$wMax_phyto), param$dx)) # Set up phytoplankton size classes
  param2$ngrid <- length(param2$w) # total number of size classes for zoo and fish
  param2$ngridPP <- length(param2$w_phyto) # total number of size classes for phyto

  # Final parameter combination
  # Exclude time series vectors from input_params since they're now stored as _ts arrays in param2
  input_params_filtered <- input_params[!names(input_params) %in% c("tmax", "dt", "isave", "time_step", "phyto_int", "phyto_slope", "phyto_max")]

  param_final <- c(input_params_filtered, param, param2)
  return(param_final)
}
