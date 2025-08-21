#' Plot Predator-Prey Mass Ratio (PPMR)
#'
#' @title Visualize predator-prey mass ratio patterns in ZooMSS results
#' @description Creates a plot showing the distribution of predator-prey mass ratios (PPMR)
#'   across functional groups, providing insights into the trophic structure of the ecosystem.
#' @details This function calculates and visualizes PPMR patterns by:
#'   - Computing theoretical PPMR values for each functional group and size class
#'   - Weighting by biomass to show realized community patterns
#'   - Creating a density plot of PPMR distribution across the community
#'   - Overlaying species-specific PPMR values as points
#'
#'   PPMR is a key ecological metric that describes the size relationship between
#'   predators and their prey, providing insight into food web structure and
#'   energy transfer efficiency in marine ecosystems.
#'
#' @param mdl ZooMSS results object containing model outputs and parameters
#' @param idx The time index to plot
#'
#' @return ggplot object showing PPMR distribution with species-specific overlays
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model
#' results <- zoomss_model(input_params, Groups)
#' ppmr_plot <- plotPPMR(results)
#' print(ppmr_plot)
#' }
#'
plotPPMR <- function(mdl, idx){

  out <- extractPPMR(mdl)

  out <- out[[idx]] # Subset to the required timestep

  ggplot2::ggplot() +
    ggplot2::geom_line(data = out[[2]], mapping = ggplot2::aes(x = .data$Betas, y = .data$y, colour = .data$Species), linewidth = 1) +
    ggplot2::geom_line(data = out[[1]], mapping = ggplot2::aes(x = .data$x, y = .data$y), linewidth = 1.2) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
    ggplot2::labs(x = expression('log' [10] * PPMR),
         y = "Zoop. Biomass Proportion", subtitle = "PPMR") +
    ggplot2::geom_vline(data = out[[1]], mapping = ggplot2::aes(xintercept = .data$mn_beta), colour = 'black') +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_colour_manual(values = c("Flagellates" = mdl$param$Groups$PlotColour[mdl$param$Groups$Species=="Flagellates"],
                                   "Ciliates" = mdl$param$Groups$PlotColour[mdl$param$Groups$Species=="Ciliates"],
                                   "Larvaceans" = mdl$param$Groups$PlotColour[mdl$param$Groups$Species=="Larvaceans"],
                                   "Salps" = mdl$param$Groups$PlotColour[mdl$param$Groups$Species=="Salps"],
                                   "Jellyfish" = mdl$param$Groups$PlotColour[mdl$param$Groups$Species=="Jellyfish"],
                                   "CarnCopepods" = mdl$param$Groups$PlotColour[mdl$param$Groups$Species=="CarnCopepods"],
                                   "Chaetognaths" = mdl$param$Groups$PlotColour[mdl$param$Groups$Species=="Chaetognaths"],
                                   "Euphausiids" = mdl$param$Groups$PlotColour[mdl$param$Groups$Species=="Euphausiids"],
                                   "OmniCopepods" = mdl$param$Groups$PlotColour[mdl$param$Groups$Species=="OmniCopepods"]))
}

#' Plot Size Spectra for ZooMSS Results
#'
#' @title Visualize abundance size spectra across functional groups
#' @description Creates a log-log plot of abundance versus body size for all functional groups,
#'   showing the classic size spectrum pattern in marine ecosystems.
#' @details This function visualizes the abundance size spectrum by:
#'   - Converting abundance data to long format with body weights
#'   - Filtering out zero abundances to focus on active size classes
#'   - Creating log-log plots colored by functional group
#'   - Using species-specific colors defined in the Groups parameter table
#'
#'   Size spectra are fundamental patterns in marine ecology, typically showing
#'   declining abundance with increasing body size. This visualization helps
#'   assess model realism and identify dominant size classes within each
#'   functional group.
#'
#' @param mdl ZooMSS results object containing model outputs and parameters
#' @param n_years The number of years (from the end) over which to average the size spectra
#'
#' @return ggplot object showing log abundance vs log body weight by species
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model
#' results <- zoomss_model(input_params, Groups)
#' size_plot <- plotSizeSpectra(results)
#' print(size_plot)
#' }
#'
plotSizeSpectra <- function(mdl, n_years) {

  species <- averageTimeSeries(mdl, "N", n_years = n_years)

  rownames(species) <- mdl$param$Groups$Species
  species <- tibble::as_tibble(t(species))

  species <- species %>%
    tibble::add_column("Weight" = mdl$param$w) %>%
    tidyr::pivot_longer(-"Weight", names_to = "Species", values_to = "Abundance") %>%
    dplyr::filter(.data$Abundance > 0) %>%
    dplyr::mutate(Species = factor(.data$Species, levels = mdl$param$Groups$Species))

  gg <- ggplot2::ggplot(data = species, mapping = ggplot2::aes(x = log10(.data$Weight), y = log10(.data$Abundance), colour = .data$Species)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = mdl$param$Groups$PlotColour) +
    ggplot2::theme_bw() +
    ggplot2::labs(subtitle = "Abundance Spectrum")

  return(gg)
}

#' Plot Time Series Data for ZooMSS Results
#'
#' @title Unified function to visualize time series changes for different metrics
#' @description Creates time series plots showing how abundance, biomass, mortality, or growth
#'   rates of functional groups change throughout the ZooMSS simulation period.
#' @details This function creates time series visualizations by:
#'   - **Abundance**: Summing abundances across size classes, log-transformed y-axis
#'   - **Biomass**: Calculating biomass (abundance × weight), with optional stacking and proportional scaling
#'   - **Mortality**: Averaging predation mortality rates across size classes
#'   - **Growth**: Averaging growth rates across size classes, log-transformed y-axis
#'
#'   All plots use species-specific colors and filter out zero values. Time series plots help identify:
#'   - Equilibration time for model runs
#'   - Seasonal or cyclical patterns in ecological metrics
#'   - Relative patterns between functional groups
#'   - Model stability and convergence behavior
#'
#' @param mdl ZooMSS results object containing model outputs with time series data
#' @param by Character string specifying the metric to plot. Options: "abundance", "biomass", "mortality", "growth" (default: "abundance")
#' @param stacked Logical, whether to create stacked area plot for biomass (default: FALSE, only applies to biomass)
#' @param proportional Logical, whether to show proportions for biomass (default: FALSE, only applies to biomass)
#' @param species Character vector of species names to include. If NULL, all species included (default: NULL, applies to all metrics)
#'
#' @return ggplot object showing the requested time series by species
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model
#' results <- zoomss_model(input_params, Groups)
#'
#' # Plot different metrics
#' abundance_plot <- plotTimeSeries(results, by = "abundance")
#' biomass_plot <- plotTimeSeries(results, by = "biomass")
#' mortality_plot <- plotTimeSeries(results, by = "mortality")
#' growth_plot <- plotTimeSeries(results, by = "growth")
#'
#' # Biomass with special options
#' stacked_plot <- plotTimeSeries(results, by = "biomass", stacked = TRUE)
#' prop_plot <- plotTimeSeries(results, by = "biomass", proportional = TRUE)
#'
#' # Focus on specific species (works for all metrics)
#' copepod_plot <- plotTimeSeries(results, by = "biomass",
#'                               species = c("OmniCopepods", "CarnCopepods"))
#' abundance_copepods <- plotTimeSeries(results, by = "abundance",
#'                                     species = c("OmniCopepods", "CarnCopepods"))
#' mortality_copepods <- plotTimeSeries(results, by = "mortality",
#'                                     species = c("OmniCopepods", "CarnCopepods"))
#' growth_copepods <- plotTimeSeries(results, by = "growth",
#'                                  species = c("OmniCopepods", "CarnCopepods"))
#' }
#'
plotTimeSeries <- function(mdl, by = "abundance", stacked = FALSE, proportional = FALSE, species = NULL) {

  # Validate inputs
  by <- match.arg(by, choices = c("abundance", "biomass", "mortality", "growth"))

  if (!("N" %in% names(mdl))) {
    stop("Time series data not available. Model may not have been run correctly.")
  }

  # Calculate data based on requested metric
  if (by == "abundance") {
    # Sum abundances across size classes for each group
    data_matrix <- rowSums(mdl$N, dims = 2)
    value_name <- "Abundance"
    y_label <- "log10(Abundance)"
    subtitle <- "Abundance"
    log_transform <- TRUE

  } else if (by == "biomass") {
    # Calculate biomass (abundance × weight) and sum across size classes
    data_matrix <- rowSums(sweep(mdl$N, 3, mdl$param$w, "*"), dims = 2)
    value_name <- "Biomass"
    y_label <- "Biomass (mg C/m^3)"
    subtitle <- "Total Biomass Time Series"
    log_transform <- FALSE

  } else if (by == "mortality") {
    # Average predation mortality rates across size classes
    data_matrix <- rowSums(mdl$Z, dims = 2) / length(mdl$param$w)
    value_name <- "Mortality"
    y_label <- "Mortality Rate"
    subtitle <- "Mortality Rate"
    log_transform <- FALSE

  } else if (by == "growth") {
    # Average growth rates across size classes
    data_matrix <- rowSums(mdl$gg, dims = 2) / length(mdl$param$w)
    value_name <- "Growth"
    y_label <- "log10(Growth Rate)"
    subtitle <- "Growth Rate"
    log_transform <- TRUE
  }

  # Set up data frame
  colnames(data_matrix) <- mdl$param$Groups$Species
  data_df <- tibble::as_tibble(data_matrix)
  data_df$Time <- mdl$time

  # Filter species if specified (applies to all plot types)
  if (!is.null(species)) {
    missing_species <- species[!species %in% mdl$param$Groups$Species]
    if (length(missing_species) > 0) {
      warning("Species not found in data: ", paste(missing_species, collapse = ", "))
    }
    valid_species <- species[species %in% mdl$param$Groups$Species]
    if (length(valid_species) == 0) {
      stop("No valid species specified. Available species: ", paste(mdl$param$Groups$Species, collapse = ", "))
    }
    data_df <- data_df[, c("Time", valid_species)]
  }

  # Convert to long format
  data_long <- data_df %>%
    tidyr::pivot_longer(-"Time", names_to = "Species", values_to = value_name) %>%
    dplyr::filter(!!rlang::sym(value_name) > 0) %>%
    dplyr::mutate(Species = factor(.data$Species, levels = mdl$param$Groups$Species))

  # Handle proportional scaling for biomass
  if (by == "biomass" && proportional && (stacked || length(unique(data_long$Species)) > 1)) {
    data_long <- data_long %>%
      dplyr::group_by(.data$Time) %>%
      dplyr::mutate(!!rlang::sym(value_name) := !!rlang::sym(value_name) / sum(!!rlang::sym(value_name), na.rm = TRUE)) %>%
      dplyr::ungroup()
    y_label <- "Proportion"
    subtitle <- "Proportional Biomass Time Series"
  }

  # Get colors for plotting
  if (!is.null(species)) {
    species_indices <- match(intersect(species, mdl$param$Groups$Species), mdl$param$Groups$Species)
    plot_colors <- mdl$param$Groups$PlotColour[species_indices]
    names(plot_colors) <- mdl$param$Groups$Species[species_indices]
  } else {
    plot_colors <- mdl$param$Groups$PlotColour
    names(plot_colors) <- mdl$param$Groups$Species
  }

  # Create plot based on options
  if (by == "biomass" && (stacked || proportional)) {
    # Stacked area plot for biomass
    if (proportional) subtitle <- "Proportional Biomass Time Series"
    else if (stacked) subtitle <- "Stacked Biomass Time Series"

    gg <- ggplot2::ggplot(data = data_long,
                          mapping = ggplot2::aes(x = .data$Time,
                                                 y = !!rlang::sym(value_name),
                                                 fill = .data$Species)) +
      ggplot2::geom_area(position = "stack", alpha = 0.7) +
      ggplot2::scale_fill_manual(values = plot_colors) +
      ggplot2::theme_bw() +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::labs(subtitle = subtitle, y = y_label) +
      ggplot2::xlab("Time (Years)")
  } else {
    # Line plot (standard for all types)
    y_aes <- if (log_transform) ggplot2::aes(y = log10(!!rlang::sym(value_name))) else ggplot2::aes(y = !!rlang::sym(value_name))

    gg <- ggplot2::ggplot(data = data_long, mapping = ggplot2::aes(x = .data$Time, colour = .data$Species)) +
      y_aes +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 1.2) +
      ggplot2::scale_color_manual(values = plot_colors) +
      ggplot2::theme_bw() +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::labs(subtitle = subtitle, y = y_label) +
      ggplot2::xlab("Time (Years)")
  }

  return(gg)
}






#' Plot Environmental Time Series
#'
#' @title Plot environmental forcing data
#' @description Creates plots of sea surface temperature and chlorophyll time series
#'   for visualizing environmental forcing data used in ZooMSS model runs.
#' @details This function creates two separate plots with different y-axes scales:
#'   - SST plot (red line) with temperature in deg C
#'   - Chlorophyll plot (green line) with concentration in mg/m^3
#'
#'   The plots can be combined using the patchwork package if available, otherwise
#'   separate plots are returned as a list. This helps users visualize the
#'   environmental forcing that drives ZooMSS model dynamics.
#'
#' @param env_data Environmental data frame with time, sst, chl columns
#'
#' @return ggplot object (if patchwork available) or list of two ggplot objects
#' @export
#'
#' @examples
#' # Create sample data and plot
#' env_data <- data.frame(
#'   time = 1:100,
#'   dt = 0.01,
#'   sst = 15 + 3*sin(2*pi*(1:100)/50),
#'   chl = 0.5 + 0.2*cos(2*pi*(1:100)/50)
#' )
#' plots <- plotEnvironment(env_data)
#'
plotEnvironment <- function(env_data) {

  # Convert to long format for plotting
  env_long <- tidyr::pivot_longer(env_data,
                                  cols = c("sst", "chl"),
                                  names_to = "variable",
                                  values_to = "value")

  # Create separate y-axes for SST and chlorophyll
  p1 <- ggplot2::ggplot(data = dplyr::filter(env_long, .data$variable == "sst"),
                        ggplot2::aes(x = .data$time, y = .data$value)) +
    ggplot2::geom_line(color = "red", linewidth = 1) +
    ggplot2::labs(y = "SST (deg C)", title = "Environmental Forcing") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title.x = ggplot2::element_blank())

  p2 <- ggplot2::ggplot(data = dplyr::filter(env_long, .data$variable == "chl"),
                        ggplot2::aes(x = .data$time, y = .data$value)) +
    ggplot2::geom_line(color = "green", linewidth = 1) +
    ggplot2::labs(x = "Time (years)", y = "Chlorophyll (mg/m^3)") +
    ggplot2::theme_bw()

  # Combine plots
  if (requireNamespace("patchwork", quietly = TRUE)) {
    return(patchwork::wrap_plots(p1, p2, ncol = 1))
  } else {
    cat("Install patchwork package to combine plots\n")
    return(list(sst_plot = p1, chl_plot = p2))
  }
}

