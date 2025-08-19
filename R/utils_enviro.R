#' Create ZooMSS Input Parameters Object
#'
#' @title Create input parameters data frame for ZooMSS model runs
#' @description Creates a properly formatted input parameters data frame for ZooMSS model
#'   simulations, combining temporal parameters with environmental time series data.
#' @details This function combines environmental time series (SST and chlorophyll) with
#'   time data to create the input_params object required by zoomss_model().
#'   The function performs validation checks using assertthat to ensure:
#'   - All input vectors are numeric and of equal length
#'   - SST values are within reasonable ocean range (-2 to 35 deg C)
#'   - Chlorophyll values are positive and within typical range (0 to 50 mg/m^3)
#'   - Time values are increasing and reasonable
#'
#'
#' @param time Numeric vector of time values in years (must be increasing and uniform, can start at any value)
#' @param sst Numeric vector of sea surface temperature values in deg C
#' @param chl Numeric vector of chlorophyll concentration values in mg/m^3
#' @param cellID Optional numeric vector of cell identifiers for spatial data (default: NULL)
#' @param isave Save frequency in time steps (default: 100)
#'
#' @return Data frame with columns: time, time_step, sst, chl, isave, and cellID (if provided)
#' @export
#'
#' @examples
#' \dontrun{
#' # Create simple environmental time series
#' time_vec <- seq(0, 10, 0.01)  # 10 years with 0.01 year time steps
#' sst_vec <- 15 + 3*sin(2*pi*time_vec/1)  # annual cycle
#' chl_vec <- 0.5 + 0.2*cos(2*pi*time_vec/1)  # annual cycle
#'
#' # Create input parameters object
#' input_params <- zCreateInputs(time_vec, sst_vec, chl_vec, isave = 50)
#'
#' # Use with ZooMSS model
#' results <- zoomss_model(input_params, Groups, SaveTimeSteps = TRUE)
#' }
#'
zCreateInputs <- function(time,
                          sst,
                          chl,
                          cellID = NULL,
                          isave = 10) {

  # Load assertthat package for validation
  if (!requireNamespace("assertthat", quietly = TRUE)) {
    stop("assertthat package required for input validation")
  }

  # Validate input data types and structure
  assertthat::assert_that(is.numeric(time), msg = "time must be numeric")
  assertthat::assert_that(is.numeric(sst), msg = "sst must be numeric")
  assertthat::assert_that(is.numeric(chl), msg = "chl must be numeric")

  # Validate equal lengths if length of sst and chl > 1
  if (length(sst) > 1 && length(chl) > 1){
    assertthat::assert_that(length(time) == length(sst),
                            msg = "time and sst must have the same length")
    assertthat::assert_that(length(time) == length(chl),
                            msg = "time and chl must have the same length")
  }

  # Validate cellID if provided
  if (!is.null(cellID)) {
    assertthat::assert_that(is.numeric(cellID), msg = "cellID must be numeric")
    assertthat::assert_that(length(cellID) == length(time),
                            msg = "cellID must have the same length as time")
  }

  # Validate time vector properties
  assertthat::assert_that(length(time) > 1, msg = "time must have at least 2 values")
  assertthat::assert_that(all(!is.na(time)), msg = "time cannot contain NA values")
  assertthat::assert_that(all(diff(time) > 0), msg = "time values must be increasing")

  # Calculate dt and tmax from time vector
  dt_values <- diff(time)
  dt <- dt_values[1]  # Use first time step as dt

  # Check if time steps are uniform - ERROR if not consistent
  max_dt_diff <- max(abs(dt_values - dt))
  if (max_dt_diff > dt * 0.001) {  # Allow only 0.1% variation (much stricter)
    stop("Time steps are not uniform. Maximum deviation: ", round(max_dt_diff, 6),
         " (", round(100 * max_dt_diff / dt, 2), "% of dt). ",
         "ZooMSS requires uniform time steps for accurate results.")
  }

  tmax <- max(time)  # Maximum time value (not duration)

  # Validate temporal parameters
  assertthat::assert_that(dt > 0, msg = "calculated dt must be positive")
  # Note: tmax can be any value (positive, negative, or zero) as it's the final time point
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
      time = time,
      time_step = seq_along(time),
      sst = sst,
      chl = chl,
      isave = isave
    )
  } else {
    formatted_data <- data.frame(
      time = time,
      time_step = seq_along(time),
      sst = sst,
      chl = chl,
      cellID = cellID,
      isave = isave
    )
  }

  # Provide summary information
  n_time_points <- nrow(formatted_data)
  n_time_steps <- n_time_points - 1

  cat("ZooMSS input parameters created:\n")
  cat("- Time points:", n_time_points, "(time values provided)\n")
  cat("- Time steps:", n_time_steps, "(intervals to simulate)\n")
  cat("- Time range:", round(min(formatted_data$time), 3), "to",
      round(max(formatted_data$time), 3), "years\n")
  cat("- SST range:", round(min(formatted_data$sst), 1), "to",
      round(max(formatted_data$sst), 1), "deg C\n")
  cat("- Chlorophyll range:", round(min(formatted_data$chl), 2), "to",
      round(max(formatted_data$chl), 2), "mg/m^3\n")
  cat("- Model parameters: dt =", round(dt, 4), "years, tmax =", round(tmax, 3), "years, isave =", isave, "steps\n")

  # Helpful reminder about time vector interpretation
  if (length(time) > 1 && all(diff(time) == 1) && min(time) %% 1 == 0 && max(time) %% 1 == 0) {
    cat("- Note: Time vector", min(time), ":", max(time), "creates", n_time_steps,
        "time steps (intervals) from", length(time), "time points.\n")
  }

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
#' @param base_chl Base chlorophyll concentration in mg/m^3 (default: 0.5)
#' @param seasonal Logical, whether to add seasonal variation (default: TRUE)
#' @param sst_amplitude Amplitude of SST seasonal variations in deg C (default: 3)
#' @param chl_amplitude Amplitude of chlorophyll seasonal variations in mg/m^3 (default: 0.2)
#'
#' @return Data frame with columns: time, sst, chl
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
#'   base_chl = 1.0
#' )
#'
zCreateSimpleTimeSeries <- function(n_time_steps,
                                    dt,
                                    base_sst = 15,
                                    base_chl = 0.5,
                                    seasonal = TRUE,
                                    sst_amplitude = 3,
                                    chl_amplitude = 0.2) {

  # Create time vector
  time_years <- seq(0, (n_time_steps-1) * dt, by = dt)

  if (seasonal) {
    # Seasonal patterns (peaks in different seasons)
    sst_values <- base_sst + sst_amplitude * sin(2 * pi * time_years)  # Annual cycle
    chl_values <- base_chl + chl_amplitude * sin(2 * pi * time_years + pi)  # Inverse to SST
  } else {
    # Static values
    sst_values <- rep(base_sst, n_time_steps)
    chl_values <- rep(base_chl, n_time_steps)
  }

  return(data.frame(
    time = time_years,
    sst = sst_values,
    chl = chl_values
  ))
}
