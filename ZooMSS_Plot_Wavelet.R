fWavelet = function(filename, species){

  library(tidyverse)
  library(WaveletComp)

  # Load the data
  dat = read_rds(filename)

  tspecies <- as_tibble(rowSums(dat$model$N, dims = 2))
  colnames(tspecies) <- dat$model$param$Groups$Species

  # -------------------------------------------------------------------------
  wt <- analyze.wavelet(tspecies, species,
                           loess.span = 0.75, dt = 1, dj = 1/20,
                           lowerPeriod = 2, upperPeriod = 333,
                           make.pval = T, n.sim = 100)

  wt.image(wt, color.key = "quantile", n.levels = 250,
           legend.params = list(lab = "Wavelet Power Spectrum2",
                                mar = 4.7, label.digits = 2, n.ticks=10), periodlab = "period")

  wt.avg(wt, my.series = 1, exponent = 1,
         show.siglvl = TRUE, siglvl = c(0.05, 0.1),
         sigcol = c("red", "blue"), sigpch = 20, sigcex = 1,
         minimum.level = NULL, maximum.level = NULL,
         label.avg.axis = TRUE,
         averagelab = NULL, averagetck = 0.02, averagetcl = 0.5,
         spec.avg.axis = list(at = NULL, labels = TRUE,
                              las = 1, hadj = NA, padj = NA),
         label.period.axis = TRUE,
         periodlab = NULL, periodtck = 0.02, periodtcl = 0.5,
         spec.period.axis = list(at = NULL, labels = TRUE,
                                 las = 1, hadj = NA, padj = NA),
         show.legend = TRUE, legend.coords = "bottomright",
         main = NULL,
         lwd = 1, col = 1,
         lwd.axis = 1,
         verbose = FALSE)
return(wt)

}