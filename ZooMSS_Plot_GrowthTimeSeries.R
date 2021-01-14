# Plot growth by time

ZooMSS_Plot_GrowthTimeSeries <- function(dat){
  library(tidyverse)
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

  gg <- ggplot(data = gg, mapping = aes(x = Time, y = log10(Growth), colour = Species)) +
    geom_line(size = 1) +
    geom_point(size = 1.2) +
    scale_color_manual(values = dat$model$param$Groups$PlotColour) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(subtitle = "Growth Rate") +
    xlab("Time (Years)")

  return(gg)
}
