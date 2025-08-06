


library(tidyverse)
library(patchwork)

source("fZooMSS_Plot.R")


load("20250806_chl0_01_0001.RData")
out1 <- out

load("20250806_chl4_0_0001.RData")
out2 <- out

load("20250806_linear_0001.RData")
out3 <- out

load("20250806_Seasonal_0001.RData")
out4 <- out


fZooMSS_Plot_AbundTimeSeries(out1) /
  fZooMSS_Plot_AbundTimeSeries(out2) /
  fZooMSS_Plot_AbundTimeSeries(out3) /
  fZooMSS_Plot_AbundTimeSeries(out4)


(fZooMSS_Plot_BiomassTimeSeries(out1, proportional = TRUE) /
  fZooMSS_Plot_BiomassTimeSeries(out2, proportional = TRUE) /
  fZooMSS_Plot_BiomassTimeSeries(out3, proportional = TRUE) /
  fZooMSS_Plot_BiomassTimeSeries(out4, proportional = TRUE)) +
  plot_layout(guides = 'collect')

spp <- c("Larvaceans", "Salps", "Jellyfish", "OmniCopepods", "Euphausiids", "CarnCopepods", "Chaetognaths")
(fZooMSS_Plot_BiomassTimeSeries(out1, proportional = TRUE, species = spp) /
    fZooMSS_Plot_BiomassTimeSeries(out2, proportional = TRUE, species = spp) /
    fZooMSS_Plot_BiomassTimeSeries(out3, proportional = TRUE, species = spp) /
    fZooMSS_Plot_BiomassTimeSeries(out4, proportional = TRUE, species = spp)) +
  plot_layout(guides = 'collect')

