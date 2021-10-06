fZooMSS_Plot_RuntimeAveraging <- function(filename, tit = ""){

  dat <- read_rds(filename)

  # Function to calculate the mean of the last pd % of the model
  pd <- seq(0.1,0.9, 0.05)
  for (p in 1:length(pd)){
    if (p == 1){
      out <- fZooMSS_AveOutput(dat$model$N, pd[p]) %>%
        list() %>%
        fZooMSS_SpeciesBiomass(dat$model) %>%
        fZooMSS_Convert2Tibble(dat$model)
    }

    if (p != 1){
      out[p,] <- fZooMSS_AveOutput(dat$model$N, pd[p]) %>%
        list() %>%
        fZooMSS_SpeciesBiomass(dat$model) %>%
        fZooMSS_Convert2Tibble(dat$model)
    }
  }

  out <- out %>%
    mutate(pd = pd) %>%
    pivot_longer(cols = Flagellates:Fish_Large, names_to = "Species", values_to = "Biomass") %>%
    mutate(Species = factor(Species, levels = dat$model$param$Groups$Species))

  gg1 <- ggplot(data = out, aes(x = pd, y = Biomass, colour = Species)) +
    geom_line() +
    scale_colour_manual(values = dat$model$param$Groups$PlotColour) +
    ggtitle(tit)


  # Now check actual run time with 50 % averaging.
  rt <- seq(100,1000, 50)

  for (r in 1:length(rt)){
    if (r == 1){
      out <- fZooMSS_AveOutput(dat$model$N[1:rt[r],,], 0.5) %>%
        list() %>%
        fZooMSS_SpeciesBiomass(dat$model) %>%
        fZooMSS_Convert2Tibble(dat$model)
    }

    if (r != 1){
      out[r,] <- fZooMSS_AveOutput(dat$model$N[1:rt[r],,], 0.5) %>%
        list() %>%
        fZooMSS_SpeciesBiomass(dat$model) %>%
        fZooMSS_Convert2Tibble(dat$model)
    }
  }

  out <- out %>%
    mutate(rt = rt) %>%
    pivot_longer(cols = Flagellates:Fish_Large, names_to = "Species", values_to = "Biomass") %>%
    mutate(Species = factor(Species, levels = dat$model$param$Groups$Species))

  gg2 <- ggplot(data = out, aes(x = rt, y = Biomass, colour = Species)) +
    geom_line() +
    scale_colour_manual(values = dat$model$param$Groups$PlotColour) +
    ggtitle(tit)

  return(list(gg1, gg2))

}

