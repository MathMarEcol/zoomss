library(tidyverse)
library(patchwork)
library(raster)
library(ggplotify)
library(circlize)

source("FoodWebs_0a_Functions.R")
source("FoodWebs_0b_Initialise.R")

data_dir <- "~/Nextcloud/MME2Work/ZooMSS/_LatestModel/20200212_Control_Full_UNSW/"
diets <- read_rds(paste0(data_dir, "Output/diets_20200212_Control_Full_UNSW.RDS"))
res <- read_rds(paste0(data_dir, "Output/res_20200212_Control_Full_UNSW.RDS"))
model <- read_rds(paste0(data_dir, "Output/ModelParameters.RDS"))
enviro_data <- read_rds(paste0(data_dir, "envirofull_20200312.RDS"))


# data_dir <- "~/Nextcloud/MME2Work/ZooMSS/_LatestModel/20200428_NoDiffusion/"
# diets <- read_rds(paste0(data_dir, "Output/diets_20200428_NoDiffusion.RDS"))
# res <- read_rds(paste0(data_dir, "Output/res_20200428_NoDiffusion.RDS"))
# model <- read_rds(paste0(data_dir, "Output/ModelParameters.RDS"))
# enviro_data <- read_rds(paste0(data_dir, "envirofull_20200317.RDS"))

enviros_oligo = which(enviro_data$chlo <= 0.1) # Oligo grid squares are where chlo is <= 0.1 mg m-3
enviros_eutro = which(enviro_data$chlo >= 1) # Eutro grid squares are where chlo is > 1 mg m-3

low_Chl <- "#e5f5f9" # "#bcdcbc"# "#E5FFCC"
high_Chl <- "#99d8c9" # "#228b22" #"#4ea24e"# "#2e8b57"# "#009900"

#### CHORD DIAGRAMS ####
colr <- TestGroups$Colour[c(3,8,4,6,5,7,9,10)]
model_color <- c("darkgreen", colr) # Add green for phytoplankton
txt <- 0.9

diets_mat_oligo = Chord_Prepare(diets, enviros_oligo, TestGroups)
diets_mat_oligo <- diets_mat_oligo[-c(2:3), -c(2:3)]
diets_mat_oligo <- diets_mat_oligo[c(1,2,7,3,5,4,6,8,9),c(1,2,7,3,5,4,6,8,9)]

diets_mat_eutro = Chord_Prepare(diets, enviros_eutro, TestGroups)
diets_mat_eutro <- diets_mat_eutro[-c(2:3), -c(2:3)]
diets_mat_eutro <- diets_mat_eutro[c(1,2,7,3,5,4,6,8,9),c(1,2,7,3,5,4,6,8,9)]

## ##  bottom, left, top and right margins
gg_fw <- as.ggplot(function(){
  circos.clear()
  par(mar = c(0, 1, 0, 2), xpd = NA, mfrow = c(1,2)) #oma = c(0,0,0,0)
  # circos.par(start.degree = 90)
  chordDiagram(diets_mat_oligo,
               directional = 1,
               grid.col = model_color,
               direction.type = c("diffHeight"),
               link.arr.type = "triangle",
               diffHeight = -uh(2, "mm"),
               annotationTrack = c("grid"),
               annotationTrackHeight = c(0.1, 0.01),
               transparency = 0.2,
               # preAllocateTracks = 1
  )
  text(0.9, -1, "Phytoplankton", cex = txt)
  text(-1.3, 0.3, "Larvaceans", cex = txt)
  text(-0.6, 0.98, "Salps", cex = txt)
  text(0, 1.3, "Omnivorous\n Copepods", cex = txt)
  text(0.7, 1, "Euphausiids", cex = txt)
  text(1.5, 0.57, "Carnivorous\n Copepods", cex = txt, adj = c(1,0))
  text(1.3, 0.6, "Chaetognaths", cex = txt)
  text(1.2, 0.6, "Jellyfish", cex = txt)
  text(1.1, 0.25, "Fish", cex = txt)
  # text(-1,1,"E)",cex = 1.2)
  # }) +
  #   labs(tag = "E)") +
  #   theme(plot.tag.position = c(-0.01, 1),
  #         plot.tag = element_text(size = 14, face = "bold"))
  #
  # gg_fw_eutro <- as.ggplot(function(){
  #   par(mar = c(0, 0, 0, 0), oma = c(0,0,0,0), xpd = TRUE)
  chordDiagram(diets_mat_eutro,
               directional = 1,
               grid.col = model_color,
               direction.type = c("diffHeight"),
               link.arr.type = "triangle",
               diffHeight = -uh(2, "mm"),
               annotationTrack = c("grid"),
               annotationTrackHeight = c(0.1, 0.01),
               transparency = 0.2,
               # preAllocateTracks = 1
  )
  text(0.6, -1.1, "Phytoplankton", cex = txt)
  text(-1.8, 0.9, "Omnivorous\n Copepods", cex = txt, adj = c(-1,0))
  text(0.45, 1, "Euphausiids", cex = txt)
  text(1.2, 0.3, "Fish", cex = txt)
  arrows(x0 = -0.85, y0 = -0.4, x1 = -0.7, y1 = -0.4, length = 0, lwd = 2)
  text(-1.2, -0.25, "Larvaceans\n Salps", cex = txt, pos = 1)
  arrows(x0 = 0.7, y0 = 0.55, x1 = 0.6, y1 =0.5, length = 0, lwd = 2)
  text(1.5, 1, "Carnivorous Copepods\n Chaetognaths\n Jellyfish", cex = txt,  pos = 2)
  # text(-1,1,"F)", cex = 1.2)
})

#  +
#   labs(tag = "F)") +
#   theme(plot.tag.position = c(-0.01, 1),
#         plot.tag = element_text(size = 14, face = "bold"))



gg_fw_eutro <- as.ggplot(function(){
  par(mar = c(0, 0, 0, 0), oma = c(0,0,0,0), xpd = TRUE)
  chordDiagram(diets_mat_eutro,
               directional = 1,
               grid.col = model_color,
               direction.type = c("diffHeight"),
               link.arr.type = "triangle",
               diffHeight = -uh(2, "mm"),
               annotationTrack = c("grid"),
               annotationTrackHeight = c(0.1, 0.01),
               transparency = 0.2)}
) +
  labs(tag = "F)") +
  theme(plot.tag.position = c(-0.01, 1),
        plot.tag = element_text(size = 14, face = "bold"))






#### PHYTO BIOMASS PROPORTION ####
PhytoBioProp <- PhytoBio_Proportion() # Get the Biomass Proportion Data

gg_numPbio <- ggplot(PhytoBioProp, aes(x = Chl, y = Proportion, fill = Taxa)) +
  scale_x_continuous(limits = c(-1.3001, 0.5001), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.001,1.001), expand = c(0, 0)) +
  labs(fill = "Taxa") +
  ylab("Biomass Proportion") +
  geom_area(alpha = 0.8) +
  scale_fill_manual(values = c("PicoPhytoplankton" = "#c7e9c0",
                               "NanoPhytoplankton" = "#a1d99b",
                               "MicroPhytoplankton" = "#74c476",
                               "Flagellates" = "#238b45",
                               "Ciliates" = "#00441b")) + # Greens: #e5f5e0 #c7e9c0 #a1d99b #74c476 #41ab5d #238b45 #006d2c #00441b
  theme_bw(base_size = 12) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        # plot.margin = margin(0,0,1,10),
        legend.key.size = unit(5, "mm"),
        panel.grid = element_line(color = "darkgrey")) +
  # annotate("text", x = -1.3, y = 0.9, label = "A)") +
  labs(tag = "A)") +
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 14, face = "bold"))


#### MESOZOOP BIOMASS PROPORTION ####
BiomassProp <- Biomass_Proportion() # Get the Biomass Proportion Data

gg_numbio <- ggplot(BiomassProp, aes(x = Chl, y = Proportion, fill = Taxa)) +
  scale_x_continuous(limits = c(-1.3001, 0.5001), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.001,1.001), expand = c(0, 0)) +
  labs(fill = "Taxa") +
  ylab("Biomass Proportion") +
  xlab(expression(paste("log"[10], "(Chlorophyll ",italic(a)," mg m"^-3, ")"))) +
  geom_area(alpha = 0.8) +
  scale_fill_manual(values = c("Larvaceans" = TestGroups$Colour[TestGroups$species=="Larvaceans"],
                               "Salps" = TestGroups$Colour[TestGroups$species=="Salps"],
                               "Jellyfish" = TestGroups$Colour[TestGroups$species=="Jellyfish"],
                               "CarnCopepods" = TestGroups$Colour[TestGroups$species=="CarnCopepods"],
                               "Chaetognaths" = TestGroups$Colour[TestGroups$species=="Chaetognaths"],
                               "Euphausiids" = TestGroups$Colour[TestGroups$species=="Euphausiids"],
                               "OmniCopepods" = TestGroups$Colour[TestGroups$species=="OmniCopepods"])) +
  theme_bw(base_size = 12) +
  theme(legend.title = element_blank(),
        plot.margin=grid::unit(c(0,0,1,0), "mm"),
        legend.key.size = unit(5, "mm"),
        # panel.background = element_rect(fill = NA),
        # panel.ontop = TRUE,
        panel.grid = element_line(color = "darkgrey")) +
  labs(tag = "B)") +
  theme(plot.tag.position = c(0, 1),
        plot.tag = element_text(size = 14, face = "bold"))


#### FISH DIET ####
diet_df <- tibble(Diet = numeric(), Chl = numeric(), Taxa = character())
Chl <- log10(enviro_data$chlo)
bins <- seq(-1.3, 0.5, by = 0.1)

for (j in bins){
  m <- (Chl > j-0.05 & Chl <= j+0.05)
  diets_mat <- Diet_Prepare(diets, m)
  temp <- tibble(Diet = as.vector(diets_mat[c(1,4,5,6,7,8,9,10,11),"PlanktivorousFish"]), Chl = j, Taxa = rownames(diets_mat[c(1, 4:11),]))
  diet_df <- bind_rows(diet_df, temp)
}

sum_diet_df <- diet_df %>%
  group_by(Chl) %>%
  summarise(SumDiet = sum(Diet, na.rm=T)) %>%
  ungroup()

diet_df <- left_join(diet_df, sum_diet_df, by = "Chl")

diet_df <- diet_df %>%
  mutate(DietProp = Diet/SumDiet,
         Taxa = as.factor(Taxa),
         Taxa = factor(Taxa, levels = c("Phytoplankton", "OmniCopepods", "Euphausiids", "Jellyfish", "Chaetognaths", "CarnCopepods", "Salps", "Larvaceans", "PlanktivorousFish")))

gg_diet <- ggplot(diet_df, aes(x = Chl, y = DietProp, fill = Taxa)) +
  scale_x_continuous(limits = c(-1.3001, 0.5001), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.001,1.001), expand = c(0, 0)) +
  # labs(subtitle = "Diet Proportion") +
  labs(fill = "Taxa") +
  ylab("Planktivorous Fish\n Diet") +
  geom_area(alpha = 0.8) +
  theme_bw(base_size = 12) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        legend.title = element_blank(),
        plot.margin=grid::unit(c(0,0,1,0), "mm"),
        legend.key.size = unit(5, "mm"),
        # panel.background = element_rect(fill = NA),
        # panel.ontop = TRUE,
        panel.grid = element_line(color = "darkgrey")) +
  scale_fill_manual(values = c("Phytoplankton" = "darkgreen",
                               "Larvaceans" = TestGroups$Colour[TestGroups$species=="Larvaceans"],
                               "Salps" = TestGroups$Colour[TestGroups$species=="Salps"],
                               "Jellyfish" = TestGroups$Colour[TestGroups$species=="Jellyfish"],
                               "CarnCopepods" = TestGroups$Colour[TestGroups$species=="CarnCopepods"],
                               "Chaetognaths" = TestGroups$Colour[TestGroups$species=="Chaetognaths"],
                               "Euphausiids" = TestGroups$Colour[TestGroups$species=="Euphausiids"],
                               "OmniCopepods" = TestGroups$Colour[TestGroups$species=="OmniCopepods"],
                               "PlanktivorousFish" = TestGroups$Colour[TestGroups$species=="Fish_Small"])) +
  labs(tag = "G)") +
  theme(plot.tag.position = c(0.0, 1),
        plot.tag = element_text(size = 14, face = "bold"))




#### TROPHIC LEVEL ####
TrophLev <- as.data.frame(TrophicLevel_Calc(diets, TestGroups))
TrophLev$Chl <- enviro_data$chlo

gg_troph <- ggplot() +
  geom_rect(aes(xmin = -Inf, xmax = -1, ymin = 2, ymax = 4.6), fill = low_Chl, alpha = 0.5) +
  geom_rect(aes(xmin = 0.1, xmax = Inf, ymin = 2, ymax = 4.6), fill = high_Chl, alpha = 0.5) +
  geom_point(data = TrophLev, aes(x = log10(Chl), y = Fish_Small), alpha = 1, size = 1, colour = "black", show.legend = FALSE) +
  geom_smooth(data = TrophLev, aes(x = log10(Chl), y = Fish_Small), colour = "blue", fill = "black", method = "loess", span = 1, se = TRUE, show.legend = TRUE) +
  # labs(subtitle = "Planktivorous Fish Trophic Level") +
  xlab(expression(paste("log"[10], "(Chlorophyll ",italic(a)," mg m"^-3, ")"))) +
  ylab("Food Chain Length") +
  scale_x_continuous(limits = c(-1.3001, 0.5001), expand = c(0, 0)) +
  scale_y_continuous(limits = c(2, 4.6), expand = c(0, 0)) +
  theme_bw(base_size = 12) +
  theme(
    # panel.ontop = TRUE,
    # panel.background = element_blank(),
    # panel.grid.major = element_line(colour = "black"),
    # axis.title.x = element_blank(),
    #     axis.text.x = element_blank(),
    legend.key.size = unit(5, "mm"),
    plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  labs(tag = "H)") +
  theme(plot.tag.position = c(0.0, 1),
        plot.tag = element_text(size = 14, face = "bold"))


# #### FISH:PHYTO RATIOS ####
# BiomassRat <- as.data.frame(biom_ratios(res, TestGroups, enviro_data))
# BiomassRat$Chl <- enviro_data$chlo
#
# gg_ratio <- ggplot(data = BiomassRat, aes(x = log10(Chl), y = BiomassRat$FishPhyto)) +
#   geom_point(alpha = 1, size = 1, colour = "black", show.legend = FALSE) +
#   geom_smooth(colour = "blue", fill = "black", method = "loess", span = 1, se = TRUE, show.legend = TRUE) +
#   xlab("log10(Chlorophyll)") +
#   ylab("Fish:Phyto") +
#   scale_x_continuous(limits = c(-1.3001, 0.5001), expand = c(0, 0)) +
#   theme_bw(base_size = 12)


#### PPMR PLOTS ####
out_list_oligo <- PPMR_plot(TestGroups, res, enviros_oligo)
out_PPMR_oligo <- out_list_oligo[[1]]
sp_PPMR_oligo <- out_list_oligo[[2]]
sp_PPMR_oligo <- sp_PPMR_oligo %>%
  mutate(Species = as.factor(Species),
  Species = factor(Species, levels = c("Flagellates", "Ciliates", "OmniCopepods", "Euphausiids", "Jellyfish", "Chaetognaths", "CarnCopepods", "Salps", "Larvaceans")))

spPPMR2 <- sp_PPMR_oligo
spPPMR2$y <- spPPMR2$y * 0
sp_PPMR_oligo <- bind_rows(sp_PPMR_oligo, spPPMR2)


out_list_eutro <- PPMR_plot(TestGroups, res, enviros_eutro)
out_PPMR_eutro <- out_list_eutro[[1]]
sp_PPMR_eutro <- out_list_eutro[[2]]
sp_PPMR_eutro <- sp_PPMR_eutro %>%
  mutate(Species = as.factor(Species),
         Species = factor(Species, levels = c("Flagellates", "Ciliates", "OmniCopepods", "Euphausiids", "Jellyfish", "Chaetognaths", "CarnCopepods", "Salps", "Larvaceans")))

spPPMR2 <- sp_PPMR_eutro
spPPMR2$y <- spPPMR2$y * 0
sp_PPMR_eutro <- bind_rows(sp_PPMR_eutro, spPPMR2)


x1 <- -1.5
x2 <- 14
y1 <- -0.001
y2 <- 0.31

gg_ppmr_o <- ggplot() +
  geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = low_Chl, alpha = 0.5) +
  geom_line(data = sp_PPMR_oligo, mapping = aes(x = Betas, y = y, colour = Species), size = 1) +
  geom_line(data = out_PPMR_oligo, mapping = aes(x = x, y = y), size = 1.2) +
  theme_bw() +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  scale_x_continuous(limits = c(x1, x2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(y1, y2), expand = c(0, 0)) +
  labs(x = expression('log' [10] * PPMR),
       y = "Biomass Proportion") +
  geom_vline(data = out_PPMR_oligo, mapping = aes(xintercept = mn_beta), colour = 'black') +
  theme_bw(base_size = 12) +
  geom_text(data = data.frame(x = x2-0.5, y = 0.29, label = "Oligotrophic"), aes(x = x, y = y, label = label), fontface = "bold", hjust = "right") +
  labs(tag = "C)") +
  theme(plot.tag.position = c(0.0, 1.03),
        plot.tag = element_text(size = 14, face = "bold"),
        legend.position = "none") +
  scale_colour_manual(values = c("Flagellates" = TestGroups$Colour[TestGroups$species=="Flagellates"],
                                 "Ciliates" = TestGroups$Colour[TestGroups$species=="Ciliates"],
                                 "Larvaceans" = TestGroups$Colour[TestGroups$species=="Larvaceans"],
                                 "Salps" = TestGroups$Colour[TestGroups$species=="Salps"],
                                 "Jellyfish" = TestGroups$Colour[TestGroups$species=="Jellyfish"],
                                 "CarnCopepods" = TestGroups$Colour[TestGroups$species=="CarnCopepods"],
                                 "Chaetognaths" = TestGroups$Colour[TestGroups$species=="Chaetognaths"],
                                 "Euphausiids" = TestGroups$Colour[TestGroups$species=="Euphausiids"],
                                 "OmniCopepods" = TestGroups$Colour[TestGroups$species=="OmniCopepods"]))

gg_ppmr_e <- ggplot() +
  geom_rect(aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2), fill = high_Chl, alpha = 0.5) +
  geom_line(data = sp_PPMR_eutro, mapping = aes(x = Betas, y = y, colour = Species), size = 1) +
  geom_line(data = out_PPMR_eutro, mapping = aes(x = x, y = y), size = 1.2) +
  theme_bw() +
  theme(plot.margin=grid::unit(c(0,0,0,0), "mm"),
        axis.text.y = element_blank(),
        legend.title = element_blank()) +
  scale_x_continuous(limits = c(x1, x2), expand = c(0, 0)) +
  scale_y_continuous(limits = c(y1, y2), expand = c(0, 0)) +
  geom_text(data = data.frame(x = x2-0.5, y = 0.29, label = "Eutrophic"), aes(x = x, y = y, label = label), fontface = "bold", hjust = "right") +
  labs(x = expression('log' [10] * PPMR), y = element_blank()) +
  geom_vline(data = out_PPMR_eutro, mapping = aes(xintercept = mn_beta), colour = 'black') +
  # theme_bw(base_size = 12) +
  labs(tag = "D)") +
  theme(plot.tag.position = c(0, 1.03),
        plot.tag = element_text(size = 14, face = "bold")) +
  scale_colour_manual(values = c("Flagellates" = TestGroups$Colour[TestGroups$species=="Flagellates"],
                                 "Ciliates" = TestGroups$Colour[TestGroups$species=="Ciliates"],
                                 "Larvaceans" = TestGroups$Colour[TestGroups$species=="Larvaceans"],
                                 "Salps" = TestGroups$Colour[TestGroups$species=="Salps"],
                                 "Jellyfish" = TestGroups$Colour[TestGroups$species=="Jellyfish"],
                                 "CarnCopepods" = TestGroups$Colour[TestGroups$species=="CarnCopepods"],
                                 "Chaetognaths" = TestGroups$Colour[TestGroups$species=="Chaetognaths"],
                                 "Euphausiids" = TestGroups$Colour[TestGroups$species=="Euphausiids"],
                                 "OmniCopepods" = TestGroups$Colour[TestGroups$species=="OmniCopepods"]),
                      guide = guide_legend(keyheight = grid::unit(5, "mm")))


#### Do PPMR/Biomass plots
graphics.off()
x11(height = 14, width = 8)
layout <- "
AA
AA
BB
BB
CC
CC
DD
DD
DD
DD
EE
EE
FF
FF
"
gg_numPbio + gg_numbio + (gg_ppmr_o + gg_ppmr_e + plot_layout(guides = 'collect')) + plot_spacer() +
     gg_diet  + gg_troph + plot_layout(design = layout)
# gg_numPbio + gg_numbio + gg_ppmr_o + gg_ppmr_e +
#   ((wrap_elements(panel = gg_fw_oligo, clip = FALSE) + wrap_elements(panel = gg_fw_eutro, clip = TRUE)) /
#      gg_diet + plot_layout(guides = 'collect'))  + gg_troph + plot_layout(design = layout)
ggsave("Figures/FoodWebs_3.png", dpi = 400)




print(
  paste("Mean PPMR for oligo is: ",round((10^out_PPMR_oligo$mn_beta[1])/1e3,3), "thousand")
)
print(
  paste("Mean PPMR for eutro is: ",round((10^out_PPMR_eutro$mn_beta[1])/1e3,3), "thousand")
)

print(
  paste("Mean TL for oligo is: ",mean(TrophLev$Fish_Small[enviros_oligo]))
)

print(
  paste("Mean TL for eutro is: ",mean(TrophLev$Fish_Small[enviros_eutro]))
)


diets <- read_rds("~/Nextcloud/MME2Work/ZooMSS/_LatestModel/20200304_GroupExperiments_Full_UNSW/NoLarvaceans/Output/diets_NoLarvaceans.RDS")
TrophLevNoLarv <- as.data.frame(TrophicLevel_Calc(diets, TestGroups))
TrophLevNoLarv$Chl <- enviro_data$chlo

print(
  paste("Mean TL for No Larvaceans in oligo is: ",mean(TrophLevNoLarv$Fish_Small[enviros_oligo]))
)



