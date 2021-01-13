ZooMSS_Plot_SizeSpectra <- function(dat) {
library(tidyverse)
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