# library(tidyverse)
# library(patchwork)
# 
# # Load all files
# fi <- list.files("RawOutput/", full.names = TRUE)
# 
# # Create empty df to store data
# df <- tibble("No" = rep(NA, length(fi)), 
#              "SST" = rep(NA, length(fi)), 
#              "Chl" = rep(NA, length(fi)), 
#              "RunTime" = rep(NA, length(fi)))
# 
# for (i in 1:length(fi)){
# 
#   dat <- read_rds(fi[i])
# 
#   df$No[i] <- i
#   df$SST[i] <- dat$model$param$sst
#   df$Chl[i] <- dat$model$param$chlo
#   df$RunTime[i] <- as.numeric(dat$model$model_runtime[3])/60
#   rm(dat)
# 
# }
# # Plot
# gg1 <- ggplot(data = df, aes(x = log10(Chl), y = RunTime)) + 
#   geom_point() + 
#   geom_smooth(method = "loess") + 
#   theme_bw() + 
#   ylab("Model Run Time (mins)")
# 
# gg2 <- ggplot(data = df, aes(x = SST, y = RunTime)) + 
#   geom_point() + 
#   geom_smooth(method = "loess") + 
#   theme_bw() + 
#   ylab("Model Run Time (mins)")
# 
# gg1 / gg2
# ggsave("Examine_Time.png", dpi = 300)


## Now lets have a look at how the output changes in 50 year increments

library(tidyverse)
library(patchwork)
fi <- list.files("~/GitHub/ZooMSS_Bryce/RawOutput/", full.names = TRUE)
tm <- seq(50,500,50)

# Create empty df to store data
out <- tibble("SST" = numeric(),
             "Chl" = numeric(),
             "Run" = numeric(),
             "Time" = numeric(),
             "Abund" = numeric(),
             "Growth" = numeric(),
             "Pred" = numeric())

for (i in 1:length(fi)){
  dat <- read_rds(fi[i])
  for (j in 1:length(tm)){
    
    sp <- dat$model$param$Groups$Species

    df <- tibble("SST" = rep(NA, length(sp)),
                 "Chl" = rep(NA, length(sp)),
                 "Run" = rep(NA, length(sp)),
                 "Time" = rep(NA, length(sp)),
                 "Abund" = rep(NA, length(sp)),
                 "Growth" = rep(NA, length(sp)),
                 "Pred" = rep(NA, length(sp)))
    
    
    
    
    for (k in 1:length(sp)){

      df$Run[k] <- i
      df$Time[k] <- tm[j]
      df$SST[k] <- dat$model$param$sst
      df$Chl[k] <- dat$model$param$chlo
      df$Species[k] <- sp[k]
      df$Abund[k] <- sum(colMeans(dat$model$N[(ceiling(0.5*tm[j]):tm[j]),k,], dims = 1))
      df$Growth[k] <- sum(colMeans(dat$model$gg[(ceiling(0.5*tm[j]):tm[j]),k,], dims = 1))
      df$Pred[k] <- sum(colMeans(dat$model$M2[(ceiling(0.5*tm[j]):tm[j]),k,], dims = 1))
    }
    
    out <- bind_rows(out, df)
    rm(df)
    
  }
}


x11(width = 8, height = 12)
out %>% 
  # filter(Species == "OmniCopepods") %>% 
  mutate(Run = as.factor(Run),
         Species = factor(Species, levels=sp),
         ChlSST = factor(paste0(Chl, "Chl-", SST, "SST"))) %>% 
  ggplot(aes(x = Time, y = log10(Abund), colour = ChlSST)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~ Species, ncol = 2) + 
  guides(colour=guide_legend(ncol=1))

ggsave("Figures/ModelRunTimes.png", dpi = 300)


