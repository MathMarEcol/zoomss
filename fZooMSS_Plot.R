library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# Plot PPMR
fZooMSS_Plot_PPMR <- function(dat){

  out <- PPMR_plot(dat)

  gg <- ggplot() +
    geom_line(data = out[[2]], mapping = aes(x = Betas, y = y, colour = Species), linewidth = 1) +
    geom_line(data = out[[1]], mapping = aes(x = x, y = y), linewidth = 1.2) +
    theme_bw() +
    theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
    labs(x = expression('log' [10] * PPMR),
         y = "Zoop. Biomass Proportion", subtitle = "PPMR") +
    geom_vline(data = out[[1]], mapping = aes(xintercept = mn_beta), colour = 'black') +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_colour_manual(values = c("Flagellates" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Flagellates"],
                                   "Ciliates" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Ciliates"],
                                   "Larvaceans" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Larvaceans"],
                                   "Salps" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Salps"],
                                   "Jellyfish" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Jellyfish"],
                                   "CarnCopepods" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="CarnCopepods"],
                                   "Chaetognaths" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Chaetognaths"],
                                   "Euphausiids" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="Euphausiids"],
                                   "OmniCopepods" = dat$model$param$Groups$PlotColour[dat$model$param$Groups$Species=="OmniCopepods"]))
}

# Plot Size Spectra
fZooMSS_Plot_SizeSpectra <- function(dat) {
  species <- dat$abundances

  rownames(species) <- dat$model$param$Groups$Species
  species <- as_tibble(t(species))

  species <- species %>%
    add_column("Weight" = dat$model$param$w) %>%
    pivot_longer(-Weight, names_to = "Species", values_to = "Abundance") %>%
    filter(Abundance > 0) %>%
    mutate(Species = factor(Species, levels = dat$model$param$Groups$Species))

  gg <- ggplot(data = species, mapping = aes(x = log10(Weight), y = log10(Abundance), colour = Species)) +
    geom_line() +
    geom_point() +
    scale_color_manual(values = dat$model$param$Groups$PlotColour) +
    theme_bw() +
    labs(subtitle = "Abundance Spectrum")

  return(gg)
}

# Plot abundance by time
fZooMSS_Plot_AbundTimeSeries <- function(dat){
  tspecies <- rowSums(dat$model$N, dims = 2)
  colnames(tspecies) <- dat$model$param$Groups$Species
  tspecies <- as_tibble(tspecies)
  tspecies$Time <- seq(dat$model$param$dt * dat$model$param$isave,
                       dat$model$param$tmax,
                       dat$model$param$dt * dat$model$param$isave)
  tspecies <- tspecies %>%
    pivot_longer(-Time, names_to = "Species", values_to = "Abundance") %>%
    filter(Abundance > 0) %>%
    mutate(Species = factor(Species, levels = dat$model$param$Groups$Species))

  gg <- ggplot(data = tspecies, mapping = aes(x = Time, y = log10(Abundance), colour = Species)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    scale_color_manual(values = dat$model$param$Groups$PlotColour) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(subtitle = "Abundance") +
    xlab("Time (Years)")

  return(gg)
}

# Plot growth by time
fZooMSS_Plot_GrowthTimeSeries <- function(dat){
  gr <- rowSums(dat$model$gg, dims = 2) / length(dat$model$param$w)
  colnames(gr) <- dat$model$param$Groups$Species
  gr <- as_tibble(gr)
  gr$Time <- seq(dat$model$param$dt * dat$model$param$isave,
                 dat$model$param$tmax,
                 dat$model$param$dt * dat$model$param$isave)
  gr <- gr %>%
    pivot_longer(-Time, names_to = "Species", values_to = "Growth") %>%
    filter(Growth > 0) %>%
    mutate(Species = factor(Species, levels = dat$model$param$Groups$Species))

  gg <- ggplot(data = gr, mapping = aes(x = Time, y = log10(Growth), colour = Species)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    scale_color_manual(values = dat$model$param$Groups$PlotColour) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(subtitle = "Growth Rate") +
    xlab("Time (Years)")

  return(gg)
}


# Plot predation by time
fZooMSS_Plot_PredTimeSeries <- function(dat){

  Z <- rowSums(dat$model$Z,dims = 2) / length(dat$model$param$w)
  colnames(Z) <- dat$model$param$Groups$Species
  Z <- as_tibble(Z)
  Z$Time <- seq(dat$model$param$dt * dat$model$param$isave,
                dat$model$param$tmax,
                dat$model$param$dt * dat$model$param$isave)
  Z <- Z %>%
    pivot_longer(-Time, names_to = "Species", values_to = "Mortality") %>%
    filter(Mortality > 0) %>%
    mutate(Species = factor(Species, levels = dat$model$param$Groups$Species))

  gg <- ggplot(data = Z, mapping = aes(x = Time, y = Mortality, colour = Species)) +
    geom_line(linewidth = 1) +
    geom_point(size = 1.2) +
    scale_color_manual(values = dat$model$param$Groups$PlotColour) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(subtitle = "Mortality Rate") +
    xlab("Time (Years)")

  return(gg)
}

# Plot Biomass Time Series
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
    pivot_longer(-Time, names_to = "Species", values_to = "Biomass") %>%
    mutate(Species = factor(Species, levels = dat$model$param$Groups$Species))
  
  # Calculate proportions if needed for proportional stacked plot
  if (proportional && (stacked || length(unique(biomass_long$Species)) > 1)) {
    biomass_long <- biomass_long %>%
      group_by(Time) %>%
      mutate(Biomass = Biomass / sum(Biomass, na.rm = TRUE)) %>%
      ungroup()
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
    
    gg <- ggplot(data = biomass_long, mapping = aes(x = Time, y = Biomass, fill = Species)) +
      geom_area(position = "stack", alpha = 0.7) +
      scale_fill_manual(values = plot_colors) +
      theme_bw() +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(subtitle = subtitle, y = y_label) +
      xlab("Time (Years)")
  } else {
    # Line plot (original)
    gg <- ggplot(data = biomass_long, mapping = aes(x = Time, y = Biomass, colour = Species)) +
      geom_line(linewidth = 1) +
      geom_point(size = 1.2) +
      scale_color_manual(values = plot_colors) +
      theme_bw() +
      scale_x_continuous(expand = c(0, 0)) +
      scale_y_continuous(expand = c(0, 0)) +
      labs(subtitle = "Total Biomass Time Series", y = "Biomass (mg C/m³)") +
      xlab("Time (Years)")
  }
  
  return(gg)
}
