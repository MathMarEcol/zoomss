## Function to calculate slope intercept and maximum size of phytoplankton spectrum, for zooplankton
## resolved size spectrum model (Heneghan et al. in prep).

# Last updated March 2020

fZooMSS_CalculatePhytoParam = function(chlo){ # chlo is chlorophyll concentration in mg m^-3
  
  ## Calculate pico, nano, micro phytoplankton proportions of total chlorophyll
  ## BREWIN ET AL., 2015
  pico <- (0.13*(1-exp(-0.8/0.13*chlo)))/chlo
  nano <- (0.77*(1-exp(-0.94/0.77*chlo)))/chlo - pico
  micro <- (chlo - 0.77*(1-exp(-0.94/0.77*chlo)))/chlo
  
  ## Convert total chlorophyll to g m^-3 total wet weight - biomass
  ## Allocate total chlorophyll to the three size classes
  c_chl <- ((chlo^0.89)*(10^1.79))/chlo # chlo:carbon ratio, from Mara??on et al. 2014
  tot_biom_c <- c_chl*chlo/1000 # (convert to grams carbon)
  tot_biom <- tot_biom_c*(1/0.1) # convert to grams wet weight, assuming 0.1 C:ww
  
  # Break up total biom into pico, nano and micro
  pico_biom <- pico*tot_biom
  nano_biom <- nano*tot_biom
  micro_biom <- micro*tot_biom

  ## Find abundances at boundaries of pico, nano size ranges, by analytically
  ## solving integral of N = aw^b

  w_0 <- -14.5 # log minimum size of picophytoplankton
  w_1 <- -11.5 # log minimum size of nanophytoplankton (max size of pico also)
  w_2 <- -8.5 # log minimum size of macrophytoplankton (max size of nano also)
    
  b <- (log10(pico_biom) - log10(nano_biom) - w_1 + w_2)/(w_1 - w_2)  # Calculate slope
  a <- pico_biom*(b+1)/((10^(w_1))^(b+1) - (10^(w_0))^(b+1)) # Calculate intercept
  
  ## Calculate maximum size
  w_max_phyto <- 0.1*round((-8.4 + 2*micro)/0.1) # Maximum size depends on the proportion of micro
  w_max_phyto <- min(-7, w_max_phyto)
  
  out <- data.frame("a" = a, "b" = b, "phyto_max" = w_max_phyto)
  return(out)
}