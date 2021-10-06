# Plot predation by time
ZooMSS_Plot_PredTimeSeries <- function(dat){
  library(tidyverse)
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
    geom_line(size = 1) +
    geom_point(size = 1.2) +
    scale_color_manual(values = dat$model$param$Groups$PlotColour) +
    theme_bw() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(subtitle = "Mortality Rate") +
    xlab("Time (Years)")

  return(gg)
}
