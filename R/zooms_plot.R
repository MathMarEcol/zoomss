#' Plot Predator-Prey Mass Ratio (PPMR) Distribution
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
#' @param dat ZooMSS results object containing model outputs and parameters
#'
#' @return ggplot object showing PPMR distribution with species-specific overlays
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model
#' results <- fZooMSS_Model(input_params, Groups, SaveTimeSteps = FALSE)
#' ppmr_plot <- fZooMSS_Plot_PPMR(results)
#' print(ppmr_plot)
#' }
#'
fZooMSS_Plot_PPMR <- function(dat){

  out <- PPMR_plot(dat)

  gg <- ggplot2::ggplot() +
    ggplot2::geom_line(data = out[[2]], mapping = ggplot2::aes(x = Betas, y = y, colour = Species), linewidth = 1) +
    ggplot2::geom_line(data = out[[1]], mapping = ggplot2::aes(x = x, y = y), linewidth = 1.2) +
    ggplot2::theme_bw() +
    ggplot2::theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
    ggplot2::labs(x = expression('log' [10] * PPMR),
         y = "Zoop. Biomass Proportion", subtitle = "PPMR") +
    ggplot2::geom_vline(data = out[[1]], mapping = ggplot2::aes(xintercept = mn_beta), colour = 'black') +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::scale_colour_manual(values = c("Flagellates" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Flagellates"],
                                   "Ciliates" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Ciliates"],
                                   "Larvaceans" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Larvaceans"],
                                   "Salps" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Salps"],
                                   "Jellyfish" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Jellyfish"],
                                   "CarnCopepods" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="CarnCopepods"],
                                   "Chaetognaths" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Chaetognaths"],
                                   "Euphausiids" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Euphausiids"],
                                   "OmniCopepods" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="OmniCopepods"]))
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
#' @param dat ZooMSS results object containing model outputs and parameters
#'
#' @return ggplot object showing log abundance vs log body weight by species
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model
#' results <- fZooMSS_Model(input_params, Groups, SaveTimeSteps = FALSE)
#' size_plot <- fZooMSS_Plot_SizeSpectra(results)
#' print(size_plot)
#' }
#'
fZooMSS_Plot_SizeSpectra <- function(dat) {
  species <- dat$abundances

  rownames(species) <- dat$model$param$Groups$Species
  species <- tibble::as_tibble(t(species))

  species <- species %>%
    tibble::add_column("Weight" = dat$model$param$w) %>%
    tidyr::pivot_longer(-Weight, names_to = "Species", values_to = "Abundance") %>%
    dplyr::filter(Abundance > 0) %>%
    dplyr::mutate(Species = factor(Species, levels = dat$model$param$Groups$Species))

  gg <- ggplot2::ggplot(data = species, mapping = ggplot2::aes(x = log10(Weight), y = log10(Abundance), colour = Species)) +
    ggplot2::geom_line() +
    ggplot2::geom_point() +
    ggplot2::scale_color_manual(values = dat$model$param$Groups$PlotColour) +
    ggplot2::theme_bw() +
    ggplot2::labs(subtitle = "Abundance Spectrum")

  return(gg)
}

#' Plot Abundance Time Series
#'
#' @title Visualize abundance changes over time for each functional group
#' @description Creates time series plots showing how total abundance of each functional
#'   group changes throughout the ZooMSS simulation period.
#' @details This function creates time series visualization by:
#'   - Summing abundances across all size classes for each functional group
#'   - Converting to long format for ggplot visualization
#'   - Plotting log-transformed abundance over time
#'   - Using species-specific colors and filtering out zero abundances
#'
#'   Time series plots help identify:
#'   - Equilibration time for model runs
#'   - Seasonal or cyclical patterns in abundance
#'   - Relative abundance patterns between functional groups
#'   - Model stability and convergence behavior
#'
#' @param dat ZooMSS results object containing model outputs with time series data
#'
#' @return ggplot object showing abundance time series by species
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model with SaveTimeSteps = TRUE
#' results <- fZooMSS_Model(input_params, Groups, SaveTimeSteps = TRUE)
#' time_plot <- fZooMSS_Plot_AbundTimeSeries(results)
#' print(time_plot)
#' }
#'
fZooMSS_Plot_AbundTimeSeries <- function(dat){
  tspecies <- rowSums(dat$model$N, dims = 2)
  colnames(tspecies) <- dat$model$param$Groups$Species
  tspecies <- tibble::as_tibble(tspecies)
  tspecies$Time <- seq(dat$model$param$dt * dat$model$param$isave,
                       dat$model$param$tmax,
                       dat$model$param$dt * dat$model$param$isave)
  tspecies <- tspecies %>%
    tidyr::pivot_longer(-Time, names_to = "Species", values_to = "Abundance") %>%
    dplyr::filter(Abundance > 0) %>%
    dplyr::mutate(Species = factor(Species, levels = dat$model$param$Groups$Species))

  gg <- ggplot2::ggplot(data = tspecies, mapping = ggplot2::aes(x = Time, y = log10(Abundance), colour = Species)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_color_manual(values = dat$model$param$Groups$PlotColour) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(subtitle = "Abundance") +
    ggplot2::xlab("Time (Years)")

  return(gg)
}

#' Plot Growth Rate Time Series
#'
#' @title Visualize growth rate changes over time for each functional group
#' @description Creates time series plots showing how average growth rates of each functional
#'   group change throughout the ZooMSS simulation period.
#' @details This function creates growth rate time series by:
#'   - Averaging growth rates across all size classes for each functional group
#'   - Converting to long format for ggplot visualization
#'   - Plotting log-transformed growth rates over time
#'   - Using species-specific colors and filtering out zero values
#'
#'   Growth rate time series help assess:
#'   - Environmental effects on organism growth
#'   - Seasonal patterns in productivity
#'   - Differences in growth potential between functional groups
#'   - Model response to changing environmental conditions
#'
#' @param dat ZooMSS results object containing model outputs with time series data
#'
#' @return ggplot object showing growth rate time series by species
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model with SaveTimeSteps = TRUE
#' results <- fZooMSS_Model(input_params, Groups, SaveTimeSteps = TRUE)
#' growth_plot <- fZooMSS_Plot_GrowthTimeSeries(results)
#' print(growth_plot)
#' }
#'
fZooMSS_Plot_GrowthTimeSeries <- function(dat){
  gr <- rowSums(dat$model$gg, dims = 2) / length(dat$model$param$w)
  colnames(gr) <- dat$model$param$Groups$Species
  gr <- tibble::as_tibble(gr)
  gr$Time <- seq(dat$model$param$dt * dat$model$param$isave,
                 dat$model$param$tmax,
                 dat$model$param$dt * dat$model$param$isave)
  gr <- gr %>%
    tidyr::pivot_longer(-Time, names_to = "Species", values_to = "Growth") %>%
    dplyr::filter(Growth > 0) %>%
    dplyr::mutate(Species = factor(Species, levels = dat$model$param$Groups$Species))

  gg <- ggplot2::ggplot(data = gr, mapping = ggplot2::aes(x = Time, y = log10(Growth), colour = Species)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_color_manual(values = dat$model$param$Groups$PlotColour) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(subtitle = "Growth Rate") +
    ggplot2::xlab("Time (Years)")

  return(gg)
}

#' Plot Predation Mortality Time Series
#'
#' @title Visualize predation mortality changes over time for each functional group
#' @description Creates time series plots showing how average predation mortality rates of each
#'   functional group change throughout the ZooMSS simulation period.
#' @details This function creates predation mortality time series by:
#'   - Averaging predation mortality rates across all size classes for each functional group
#'   - Converting to long format for ggplot visualization
#'   - Plotting mortality rates over time without log transformation
#'   - Using species-specific colors and filtering out zero values
#'
#'   Predation mortality time series help assess:
#'   - Predation pressure on different functional groups over time
#'   - Seasonal or temporal patterns in predation intensity
#'   - Relative vulnerability of functional groups to predation
#'   - Model dynamics and predator-prey interactions
#'
#' @param dat ZooMSS results object containing model outputs with time series data
#'
#' @return ggplot object showing predation mortality time series by species
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model with SaveTimeSteps = TRUE
#' results <- fZooMSS_Model(input_params, Groups, SaveTimeSteps = TRUE)
#' mortality_plot <- fZooMSS_Plot_PredTimeSeries(results)
#' print(mortality_plot)
#' }
#'
fZooMSS_Plot_PredTimeSeries <- function(dat){

  Z <- rowSums(dat$model$Z,dims = 2) / length(dat$model$param$w)
  colnames(Z) <- dat$model$param$Groups$Species
  Z <- tibble::as_tibble(Z)
  Z$Time <- seq(dat$model$param$dt * dat$model$param$isave,
                dat$model$param$tmax,
                dat$model$param$dt * dat$model$param$isave)
  Z <- Z %>%
    tidyr::pivot_longer(-Time, names_to = "Species", values_to = "Mortality") %>%
    dplyr::filter(Mortality > 0) %>%
    dplyr::mutate(Species = factor(Species, levels = dat$model$param$Groups$Species))

  gg <- ggplot2::ggplot(data = Z, mapping = ggplot2::aes(x = Time, y = Mortality, colour = Species)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_color_manual(values = dat$model$param$Groups$PlotColour) +
    ggplot2::theme_bw() +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::labs(subtitle = "Mortality Rate") +
    ggplot2::xlab("Time (Years)")

  return(gg)
}

#' Plot Biomass Time Series
#'
#' @title Visualize biomass changes over time with multiple display options
#' @description Creates flexible time series plots showing how total biomass of functional
#'   groups changes throughout the ZooMSS simulation, with options for line plots,
#'   stacked area plots, and proportional displays.
#' @details This function creates biomass time series visualization with multiple options:
#'   - **Line plots**: Individual species biomass trajectories over time
#'   - **Stacked plots**: Cumulative biomass showing total ecosystem biomass
#'   - **Proportional plots**: Relative biomass contributions (0-1 scale)
#'   - **Species filtering**: Focus on specific functional groups
#'
#'   The function calculates biomass by multiplying abundance by body weights and
#'   summing across size classes for each functional group. Different plot types help
#'   visualize different aspects of ecosystem dynamics:
#'   - Line plots show individual group patterns and relative magnitudes
#'   - Stacked plots show total ecosystem biomass and contributions
#'   - Proportional plots highlight shifts in community composition
#'
#' @param dat ZooMSS results object containing model outputs with time series data
#' @param stacked Logical, whether to create stacked area plot instead of line plot (default: FALSE)
#' @param proportional Logical, whether to show proportions instead of absolute values (default: FALSE)
#' @param species Character vector of species names to include in plot. If NULL, all species included (default: NULL)
#'
#' @return ggplot object showing biomass time series by species
#' @export
#'
#' @examples
#' \dontrun{
#' # After running ZooMSS model with SaveTimeSteps = TRUE
#' results <- fZooMSS_Model(input_params, Groups, SaveTimeSteps = TRUE)
#'
#' # Basic line plot of all species
#' biomass_plot <- fZooMSS_Plot_BiomassTimeSeries(results)
#'
#' # Stacked area plot showing total biomass
#' stacked_plot <- fZooMSS_Plot_BiomassTimeSeries(results, stacked = TRUE)
#'
#' # Proportional plot showing relative contributions
#' prop_plot <- fZooMSS_Plot_BiomassTimeSeries(results, proportional = TRUE)
#'
#' # Focus on specific groups
#' copepod_plot <- fZooMSS_Plot_BiomassTimeSeries(results,
#'                                               species = c("OmniCopepods", "CarnCopepods"))
#' }
#'
fZooMSS_Plot_BiomassTimeSeries <- function(dat, stacked = FALSE, proportional = FALSE, species = NULL){
  if (!("N" %in% names(dat$model))) {
    stop("Abundance data not available. Make sure SaveTimeSteps=TRUE when running the model.")
  }

  # Calculate biomass from abundance and weights
  # dat$model$N dims: [time, groups, sizes]
  # dat$model$param$w: weights for each size class
  # Result: sum across size classes for each group at each time step -> [time, groups]
  biomass <- rowSums(sweep(dat$model$N, 3, dat$model$param$w, "*"), dims = 2)
  time_steps <- seq_len(nrow(biomass))
  time_years <- time_steps * dat$model$param$dt * dat$model$param$isave

  # Create long format dataframe
  biomass_df <- as.data.frame(biomass)
  colnames(biomass_df) <- dat$model$param$Groups$Species
  biomass_df$Time <- time_years

  # Filter species if specified
  if (!is.null(species)) {
    # Check if specified species exist
    missing_species <- species[!species %in% dat$model$param$Groups$Species]
    if (length(missing_species) > 0) {
      warning("Species not found in data: ", paste(missing_species, collapse = ", "))
    }
    # Keep only specified species that exist
    valid_species <- species[species %in% dat$model$param$Groups$Species]
    if (length(valid_species) == 0) {
      stop("No valid species specified. Available species: ", paste(dat$model$param$Groups$Species, collapse = ", "))
    }
    # Select only specified species columns plus Time
    biomass_df <- biomass_df[, c("Time", valid_species)]
  }

  # Convert to long format
  biomass_long <- biomass_df %>%
    tidyr::pivot_longer(-Time, names_to = "Species", values_to = "Biomass") %>%
    dplyr::mutate(Species = factor(Species, levels = dat$model$param$Groups$Species))

  # Calculate proportions if needed for proportional stacked plot
  if (proportional && (stacked || length(unique(biomass_long$Species)) > 1)) {
    biomass_long <- biomass_long %>%
      dplyr::group_by(Time) %>%
      dplyr::mutate(Biomass = Biomass / sum(Biomass, na.rm = TRUE)) %>%
      dplyr::ungroup()
  }

  # Get colors for selected species
  if (!is.null(species)) {
    # Get indices of selected species
    species_indices <- match(intersect(species, dat$model$param$Groups$Species), dat$model$param$Groups$Species)
    plot_colors <- dat$model$param$Groups$PlotColour[species_indices]
    names(plot_colors) <- dat$model$param$Groups$Species[species_indices]
  } else {
    plot_colors <- dat$model$param$Groups$PlotColour
    names(plot_colors) <- dat$model$param$Groups$Species
  }

  # Create plot based on options
  if (stacked || proportional) {
    # Stacked area plot (absolute or proportional)
    y_label <- if (proportional) "Proportion" else "Biomass (mg C/m³)"
    subtitle <- if (proportional) "Proportional Biomass Time Series" else "Stacked Biomass Time Series"

    gg <- ggplot2::ggplot(data = biomass_long, mapping = ggplot2::aes(x = Time, y = Biomass, fill = Species)) +
      ggplot2::geom_area(position = "stack", alpha = 0.7) +
      ggplot2::scale_fill_manual(values = plot_colors) +
      ggplot2::theme_bw() +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::labs(subtitle = subtitle, y = y_label) +
      ggplot2::xlab("Time (Years)")
  } else {
    # Line plot (original)
    gg <- ggplot2::ggplot(data = biomass_long, mapping = ggplot2::aes(x = Time, y = Biomass, colour = Species)) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 1.2) +
      ggplot2::scale_color_manual(values = plot_colors) +
      ggplot2::theme_bw() +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::scale_y_continuous(expand = c(0, 0)) +
      ggplot2::labs(subtitle = "Total Biomass Time Series", y = "Biomass (mg C/m³)") +
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
#'   - SST plot (red line) with temperature in °C
#'   - Chlorophyll plot (green line) with concentration in mg/m³
#'
#'   The plots can be combined using the patchwork package if available, otherwise
#'   separate plots are returned as a list. This helps users visualize the
#'   environmental forcing that drives ZooMSS model dynamics.
#'
#' @param env_data Environmental data frame with time, sst, chlo columns
#'
#' @return ggplot object (if patchwork available) or list of two ggplot objects
#' @export
#'
#' @examples
#' # Create sample data and plot
#' env_data <- data.frame(
#'   time_step = 1:100,
#'   dt = 0.01,
#'   sst = 15 + 3*sin(2*pi*(1:100)/50),
#'   chlo = 0.5 + 0.2*cos(2*pi*(1:100)/50)
#' )
#' plots <- fZooMSS_PlotEnvironment(env_data)
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

