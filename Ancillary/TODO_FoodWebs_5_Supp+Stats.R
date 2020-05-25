library(tidyverse)
library(patchwork)

source("FoodWebs_0a_Functions.R")
source("FoodWebs_0b_Initialise.R")

# diets <- read_rds("/Users/jason/Dropbox/Multi-Zoo Size Spectrum Model/_LatestModel/20200212_Control_Full_UNSW/Output/diets_20200212_Control_Full_UNSW.RDS")
# res <- read_rds("/Users/jason/Dropbox/Multi-Zoo Size Spectrum Model/_LatestModel/20200212_Control_Full_UNSW/Output/res_20200212_Control_Full_UNSW.RDS")
# model <- read_rds("/Users/jason/Dropbox/Multi-Zoo Size Spectrum Model/_LatestModel/20200212_Control_Full_UNSW/Output/ModelParameters.RDS")
# enviro_data <- read_rds("/Users/jason/Dropbox/Multi-Zoo Size Spectrum Model/_LatestModel/20200212_Control_Full_UNSW/envirofull_20200209.RDS")

#### TROPHIC LEVEL ####
# TrophLev <- as.data.frame(TrophicLevel_Calc(diets, TestGroups))
# TrophLev$Chl <- enviro_data$chlo

enviro200 <- read_rds("/Users/jason/Dropbox/Multi-Zoo Size Spectrum Model/_LatestModel/20200209_GroupExperiments_UNSW/NoLarvaceans/enviro200_20200209.RDS")

diets_noLarv <- read_rds("/Users/jason/Dropbox/Multi-Zoo Size Spectrum Model/_LatestModel/20200209_GroupExperiments_UNSW/NoLarvaceans/Output/diets_NoLarvaceans.RDS")
TrophLev_noLarv <- as.data.frame(TrophicLevel_Calc(diets_noLarv, TestGroups))
TrophLev_noLarv$Chl <- enviro200$chlo

diets_noSalps <- read_rds("/Users/jason/Dropbox/Multi-Zoo Size Spectrum Model/_LatestModel/20200209_GroupExperiments_UNSW/NoSalps/Output/diets_NoSalps.RDS")
TrophLev_noSalps <- as.data.frame(TrophicLevel_Calc(diets_noSalps, TestGroups))
TrophLev_noSalps$Chl <- enviro200$chlo

diets_all <- read_rds("/Users/jason/Dropbox/Multi-Zoo Size Spectrum Model/_LatestModel/20200204_Control_200cells_UNSW/Output/diets_20200204_Control_200cells_UNSW.RDS")
TrophLev_all <- as.data.frame(TrophicLevel_Calc(diets_all, TestGroups))
TrophLev_all$Chl <- enviro200$chlo

## I need to do a delta plot of the TL change. Wait for my full model runs to finish.

larv_diff = TrophLev_all$Fish_Small - TrophLev_noLarv$Fish_Small
salps_diff = TrophLev_all$Fish_Small - TrophLev_noSalps$Fish_Small
# fish_props = ((fish_bioms-control_fish_bioms)/control_fish_bioms)*100

TrophLev <- tibble(Chl = enviro200$chlo, diff_NoLarv = larv_diff, diff_NoSalps = salps_diff)

gg_troph <- ggplot(data = TrophLev, aes(x = log10(Chl), y = diff_NoLarv)) +
  geom_point(alpha = 1, size = 1, colour = TestGroups$Colour[TestGroups$species=="Larvaceans"], show.legend = FALSE) +
  geom_smooth(colour = TestGroups$Colour[TestGroups$species=="Larvaceans"], fill = TestGroups$Colour[TestGroups$species=="Larvaceans"], method = "loess", span = 1, se = TRUE, show.legend = TRUE) +
  geom_point(data = TrophLev, aes(x = log10(Chl), y = diff_NoSalps),
             alpha = 1, size = 1, colour = TestGroups$Colour[TestGroups$species=="Salps"], show.legend = FALSE) +
  geom_smooth(data = TrophLev, aes(x = log10(Chl), y = diff_NoSalps),
              colour = TestGroups$Colour[TestGroups$species=="Salps"], fill = TestGroups$Colour[TestGroups$species=="Salps"], method = "loess", span = 1, se = TRUE, show.legend = TRUE) +
  ylab("Trophic Level") +
  scale_fill_manual(values = TestGroups$Colour[TestGroups$Feeding=="FilterFeeder"],
                    labels = TestGroups$CommonName[TestGroups$Feeding=="FilterFeeder"],
                    name = element_blank()) +
  scale_x_continuous(limits = c(-1.3001, 0.5001), expand = c(0, 0)) +
  theme_bw(base_size = 14) +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),
        legend.position="right") +
  xlab("log10(Chlorophyll)") +
  ylab("Planktivorous Fish Trophic Level")

graphics.off()
x11(width = 10, height = 5)
gg_troph
ggsave(filename="Figures/FoodWebs_Supp_TrophicLevel_Gelatinous.png", dpi = 400)

