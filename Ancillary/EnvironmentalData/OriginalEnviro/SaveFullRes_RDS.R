
library(tidyverse)
source("PhytoplanktonModelAnalytical.R")

enviro <- read_csv("enviro_5d_3.csv")

enviro_orig <- enviro %>%
  filter(is.na(chlo)==FALSE) %>%
  filter(sst >= 1)

enviro_cold <- enviro %>%
  filter(is.na(chlo)==FALSE) %>%
  filter(sst < 1)

# Rejoin them here so that the old rows are on the top and I know exactly which ones to run moving forward.
enviro <- bind_rows(enviro_orig, enviro_cold)

enviro <- enviro %>% 
  select(-c("a", "b", "phyto_max"))

p_save <- tibble("a" = NA, "b" = NA, "phyto_max" = NA)
for (p in 1:length(enviro$chlo)){
  out <- phyto_info(enviro$chlo[p])
  out$a <- log10(out$a)
  p_save[p,] <- out  
}

enviro <- bind_cols(enviro, p_save)
enviro$dt <- 0.01

write_rds(enviro,paste0("envirofull_",str_replace_all(ymd(today()),"-",""),".RDS"))
