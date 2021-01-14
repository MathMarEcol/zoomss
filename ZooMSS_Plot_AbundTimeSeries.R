# Plot abundance by time
ZooMSS_Plot_AbundTimeSeries <- function(dat){
  library(tidyverse)
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
    geom_line(size = 1) +
    geom_point(size = 1.2) +
    scale_color_manual(values = dat$model$param$Groups$PlotColour) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(subtitle = "Abundance") +
    xlab("Time (Years)")

  return(gg)
}