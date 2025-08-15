#' Create ZooMSS Input Parameters Object
#'
#' @title Create input parameters data frame for ZooMSS model runs
#' @description Creates a properly formatted input parameters data frame for ZooMSS model
#'   simulations, combining temporal parameters with environmental time series data.
#' @details This function combines environmental time series (SST and chlorophyll) with
#'   model temporal parameters to create the input_params object required by zoomss_model().
#'   The function performs validation checks using assertthat to ensure:
#'   - All input vectors are numeric and of equal length
#'   - SST values are within reasonable ocean range (-2 to 35 deg C)
#'   - Chlorophyll values are positive and within typical range (0 to 50 mg/m^3)
#'   - Temporal parameters are positive and reasonable
#'
#'   The resulting data frame includes time_step indices, environmental data, and
#'   model parameters needed for ZooMSS simulations.
#'
#' @param time Numeric vector of time values (any units, used for time_step sequence)
#' @param sst Numeric vector of sea surface temperature values in deg C
#' @param chl Numeric vector of chlorophyll concentration values in mg/m^3
#' @param cellID Optional numeric vector of cell identifiers for spatial data (default: NULL)
#' @param dt Time step size in years (default: 0.01)
#' @param tmax Maximum simulation time in years (default: 250)
#' @param isave Save frequency in time steps (default: 100)
#'
#' @return Data frame with columns: time_step, sst, chlo, dt, tmax, isave, and cellID (if provided)
#' @export
#'
#' @examples
#' \dontrun{
#' # Create simple environmental time series
#' time_vec <- 1:100
#' sst_vec <- 15 + 3*sin(2*pi*time_vec/50)
#' chl_vec <- 0.5 + 0.2*cos(2*pi*time_vec/50)
#'
#' # Create input parameters object
#' input_params <- zcreateInputs(time_vec, sst_vec, chl_vec,
#'                                      dt = 0.01, tmax = 10, isave = 50)
#'
#' # Use with ZooMSS model
#' results <- zoomss_model(input_params, Groups, SaveTimeSteps = TRUE)
#' }
#'
zcreateInputs <- function(time, sst, chl,
                                 cellID = NULL,
                                 dt = 0.01,
                                 tmax = 250,
                                 isave = 100) {

  # Load assertthat package for validation
  if (!requireNamespace("assertthat", quietly = TRUE)) {
    stop("assertthat package required for input validation")
  }

  # Validate input data types and structure
  assertthat::assert_that(is.numeric(time), msg = "time must be numeric")
  assertthat::assert_that(is.numeric(sst), msg = "sst must be numeric")
  assertthat::assert_that(is.numeric(chl), msg = "chl must be numeric")

  # Validate equal lengths
  assertthat::assert_that(length(time) == length(sst),
                         msg = "time and sst must have the same length")
  assertthat::assert_that(length(time) == length(chl),
                         msg = "time and chl must have the same length")

  # Validate cellID if provided
  if (!is.null(cellID)) {
    assertthat::assert_that(is.numeric(cellID), msg = "cellID must be numeric")
    assertthat::assert_that(length(cellID) == length(time),
                           msg = "cellID must have the same length as time")
  }

  # Validate temporal parameters
  assertthat::assert_that(is.numeric(dt) && length(dt) == 1 && dt > 0,
                         msg = "dt must be a positive number")
  assertthat::assert_that(is.numeric(tmax) && length(tmax) == 1 && tmax > 0,
                         msg = "tmax must be a positive number")
  assertthat::assert_that(is.numeric(isave) && length(isave) == 1 && isave > 0,
                         msg = "isave must be a positive number")

  # Validate environmental data ranges
  assertthat::assert_that(all(!is.na(sst)), msg = "sst cannot contain NA values")
  assertthat::assert_that(all(!is.na(chl)), msg = "chl cannot contain NA values")
  assertthat::assert_that(all(sst >= -2 & sst <= 35),
                         msg = "sst values must be within ocean range (-2 to 35 deg C)")
  assertthat::assert_that(all(chl >= 0 & chl <= 50),
                         msg = "chl values must be within range (0 to 50 mg/m^3)")

  # Create formatted data frame
  if (is.null(cellID)) {
    formatted_data <- data.frame(
      time_step = seq_along(time),
      sst = sst,
      chlo = chl,
      dt = dt,
      tmax = tmax,
      isave = isave
    )
  } else {
    formatted_data <- data.frame(
      time_step = seq_along(time),
      sst = sst,
      chlo = chl,
      cellID = cellID,
      dt = dt,
      tmax = tmax,
      isave = isave
    )
  }

  # Provide summary information
  cat("ZooMSS input parameters created:\n")
  cat("- Time steps:", nrow(formatted_data), "\n")
  cat("- SST range:", round(min(formatted_data$sst), 1), "to",
      round(max(formatted_data$sst), 1), "deg C\n")
  cat("- Chlorophyll range:", round(min(formatted_data$chlo), 2), "to",
      round(max(formatted_data$chlo), 2), "mg/m^3\n")
  cat("- Model parameters: dt =", dt, "years, tmax =", tmax, "years, isave =", isave, "steps\n")

  return(formatted_data)
}


#' Create Simple Environmental Time Series for Testing
#'
#' @title Generate synthetic environmental data for ZooMSS testing
#' @description Creates simple synthetic environmental time series with optional seasonal
#'   variation for testing ZooMSS model runs when real environmental data is not available.
#' @details This function generates synthetic sea surface temperature and chlorophyll
#'   time series that can be used for testing ZooMSS model behavior. The function can
#'   create either static environmental conditions or seasonal cycles with sinusoidal
#'   variation. This is particularly useful for:
#'   - Testing model sensitivity to environmental forcing
#'   - Creating idealized scenarios for model exploration
#'   - Generating data when real environmental data is unavailable
#'
#'   The seasonal option creates SST and chlorophyll cycles that are out of phase,
#'   mimicking typical ocean patterns where chlorophyll peaks when SST is lower.
#'
#' @param n_time_steps Number of time steps to generate
#' @param dt Time step size in years
#' @param base_sst Base sea surface temperature in deg C (default: 15)
#' @param base_chlo Base chlorophyll concentration in mg/m^3 (default: 0.5)
#' @param seasonal Logical, whether to add seasonal variation (default: TRUE)
#' @param sst_amplitude Amplitude of SST seasonal variations in deg C (default: 3)
#' @param chlo_amplitude Amplitude of chlorophyll seasonal variations in mg/m^3 (default: 0.2)
#'
#' @return Data frame with columns: time, sst, chlo
#' @export
#'
#' @examples
#' # Create seasonal environmental data
#' env_data <- zCreateSimpleTimeSeries(
#'   n_time_steps = 1000,
#'   dt = 0.01,
#'   seasonal = TRUE
#' )
#'
#' # Create static environmental conditions
#' static_data <- zCreateSimpleTimeSeries(
#'   n_time_steps = 500,
#'   dt = 0.01,
#'   seasonal = FALSE,
#'   base_sst = 20,
#'   base_chlo = 1.0
#' )
#'
zCreateSimpleTimeSeries <- function(n_time_steps, dt,
                                          base_sst = 15, base_chlo = 0.5,
                                          seasonal = TRUE,
                                          sst_amplitude = 3, chlo_amplitude = 0.2) {

  # Create time vector
  time_years <- seq(0, (n_time_steps-1) * dt, by = dt)

  if (seasonal) {
    # Seasonal patterns (peaks in different seasons)
    sst_values <- base_sst + sst_amplitude * sin(2 * pi * time_years)  # Annual cycle
    chlo_values <- base_chlo + chlo_amplitude * sin(2 * pi * time_years + pi)  # Inverse to SST
  } else {
    # Static values
    sst_values <- rep(base_sst, n_time_steps)
    chlo_values <- rep(base_chlo, n_time_steps)
  }

  return(data.frame(
    time = time_years,
    sst = sst_values,
    chlo = chlo_values
  ))
}
