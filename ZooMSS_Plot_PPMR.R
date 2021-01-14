ZooMSS_Plot_PPMR <- function(dat){
  library(tidyverse)
  source("fZooMSS_Xtras.R")

  out <- PPMR_plot(dat)

  gg <- ggplot() +
    geom_line(data = out[[2]], mapping = aes(x = Betas, y = y, colour = Species), size = 1) +
    geom_line(data = out[[1]], mapping = aes(x = x, y = y), size = 1.2) +
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

  return(gg)
}
