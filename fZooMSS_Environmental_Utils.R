## Utility functions for processing environmental data for dynamic ZooMSS
## These functions help users work with real environmental data

#' Format real environmental data for ZooMSS
#'
#' @param env_data Data frame with environmental data
#' @param time_col Name of time column
#' @param sst_col Name of SST column
#' @param chlo_col Name of chlorophyll column
#' @param dt Time step for model (years)
#'
#' @return Formatted data frame ready for ZooMSS
#'
fZooMSS_FormatRealData <- function(env_data, time_col = "time", sst_col = "sst", chlo_col = "chlo", dt = 0.01) {

  # Standardize column names
  formatted_data <- data.frame(
    time_step = seq_len(nrow(env_data)),
    time = env_data[[time_col]],
    sst = env_data[[sst_col]],
    chlo = env_data[[chlo_col]]
  )

  # Validate data
  if (any(is.na(formatted_data$sst))) {
    warning("SST data contains NA values")
  }
  if (any(is.na(formatted_data$chlo))) {
    warning("Chlorophyll data contains NA values")
  }

  # Check for reasonable ranges
  if (any(formatted_data$sst < -2 | formatted_data$sst > 35)) {
    warning("SST values outside typical ocean range (-2 to 35°C)")
  }
  if (any(formatted_data$chlo < 0 | formatted_data$chlo > 50)) {
    warning("Chlorophyll values outside typical range (0 to 50 mg/m³)")
  }

  cat("✅ Environmental data formatted:\n")
  cat("- Time steps:", nrow(formatted_data), "\n")
  cat("- SST range:", round(min(formatted_data$sst, na.rm=TRUE), 1), "to",
      round(max(formatted_data$sst, na.rm=TRUE), 1), "°C\n")
  cat("- Chlorophyll range:", round(min(formatted_data$chlo, na.rm=TRUE), 2), "to",
      round(max(formatted_data$chlo, na.rm=TRUE), 2), "mg/m³\n")

  return(formatted_data)
}

#' Plot environmental time series
#'
#' @param env_data Environmental data frame with time, sst, chlo columns
#'
#' @return ggplot object
#'
fZooMSS_PlotEnvironment <- function(env_data) {

  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 package required for plotting")
  }
  if (!requireNamespace("tidyr", quietly = TRUE)) {
    stop("tidyr package required for plotting")
  }

  # Convert to long format for plotting
  env_long <- tidyr::pivot_longer(env_data,
                                  cols = c("sst", "chlo"),
                                  names_to = "variable",
                                  values_to = "value")

  # Create separate y-axes for SST and chlorophyll
  p1 <- ggplot2::ggplot(data = subset(env_long, variable == "sst"),
                        ggplot2::aes(x = time_step*dt, y = value)) +
    ggplot2::geom_line(color = "red", linewidth = 1) +
    ggplot2::labs(y = "SST (°C)", title = "Environmental Forcing") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())

  p2 <- ggplot2::ggplot(data = subset(env_long, variable == "chlo"),
                        ggplot2::aes(x = time_step*dt, y = value)) +
    ggplot2::geom_line(color = "green", linewidth = 1) +
    ggplot2::labs(x = "Time (years)", y = "Chlorophyll (mg/m³)") +
    ggplot2::theme_bw()

  # Combine plots
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(patchwork::wrap_plots(p1, p2, ncol = 1))
  } else {
    cat("Install gridExtra package to combine plots\n")
    return(list(sst_plot = p1, chlo_plot = p2))
  }
}

#' Create simple environmental time series for testing
#'
#' @param n_time_steps Number of time steps
#' @param dt Time step size in years
#' @param base_sst Base sea surface temperature (°C)
#' @param base_chlo Base chlorophyll concentration (mg/m³)
#' @param seasonal Logical, add seasonal variation?
#' @param sst_amplitude Amplitude of SST variations (°C)
#' @param chlo_amplitude Amplitude of chlorophyll variations (mg/m³)
#'
#' @return Data frame with time, sst, and chlo columns
#'
fZooMSS_CreateSimpleTimeSeries <- function(n_time_steps, dt,
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
