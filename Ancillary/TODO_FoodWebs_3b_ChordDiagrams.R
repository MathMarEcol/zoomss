library(circlize)
library(tidyverse)

source("FoodWebs_0a_Functions.R")
source("FoodWebs_0b_Initialise.R")

# data_dir <- "~/Nextcloud/MME2Work/ZooMSS/_LatestModel/20200212_Control_Full_UNSW/"
# diets <- read_rds(paste0(data_dir, "Output/diets_20200212_Control_Full_UNSW.RDS"))
# res <- read_rds(paste0(data_dir, "Output/res_20200212_Control_Full_UNSW.RDS"))
# model <- read_rds(paste0(data_dir, "Output/ModelParameters.RDS"))
# enviro_data <- read_rds(paste0(data_dir, "envirofull_20200312.RDS"))

data_dir <- "~/Nextcloud/MME2Work/ZooMSS/_LatestModel/20200428_NoDiffusion/"
diets <- read_rds(paste0(data_dir, "Output/diets_20200428_NoDiffusion.RDS"))
res <- read_rds(paste0(data_dir, "Output/res_20200428_NoDiffusion.RDS"))
model <- read_rds(paste0(data_dir, "Output/ModelParameters.RDS"))
enviro_data <- read_rds(paste0(data_dir, "envirofull_20200317.RDS"))


enviros_oligo = which(enviro_data$chlo <= 0.1) # Oligo grid squares are where chlo is <= 0.1 mg m-3
enviros_eutro = which(enviro_data$chlo >= 1) # Eutro grid squares are where chlo is > 1 mg m-3

#### CHORD DIARGRAMS ####

### DO Chords without Hetero
diets_mat_oligo = Chord_Prepare(diets, enviros_oligo, TestGroups)
diets_mat_oligo <- diets_mat_oligo[-c(2:3), -c(2:3)]
diets_mat_oligo <- diets_mat_oligo[c(1,2,7,3,5,4,6,8,9),c(1,2,7,3,5,4,6,8,9)]

diets_mat_eutro = Chord_Prepare(diets, enviros_eutro, TestGroups)
diets_mat_eutro <- diets_mat_eutro[-c(2:3), -c(2:3)]
diets_mat_eutro <- diets_mat_eutro[c(1,2,7,3,5,4,6,8,9),c(1,2,7,3,5,4,6,8,9)]

enviros_new <- 1429:1638
diets_mat_new = Chord_Prepare(diets, enviros_new, TestGroups)
diets_mat_new <- diets_mat_new[-c(2:3), -c(2:3)]
diets_mat_new <- diets_mat_new[c(1,2,7,3,5,4,6,8,9),c(1,2,7,3,5,4,6,8,9)]

colr <- TestGroups$Colour[c(3,8,4,6,5,7,9,10)]

model_color = c("darkgreen", colr) # Add green for phytoplankton

x11(width = 11, height = 5)
 ##  bottom, left, top and right margins
par(mar = c(0, 3, 0, 3), mfrow = c(1,2), xpd = NA)

### OLIGOTROPHIC FOOD WEB
chordDiagram(diets_mat_oligo,
             directional = 1,
             grid.col = model_color,
             direction.type = c("diffHeight"),
             link.arr.type = "triangle",
             diffHeight = -uh(2, "mm"),
             annotationTrack = c("grid"),
             annotationTrackHeight = c(0.1, 0.01),
             transparency = 0.2)

txt <- 0.9
text(0.6, -1, "Phytoplankton", cex = txt)
text(-1.1, 0.45, "Larvaceans", cex = txt)
text(-0.55, 0.95, "Salps", cex = txt)
text(0, 1.1, "Omnivorous\n Copepods", cex = txt)
text(0.5, 1, "Euphausiids", cex = txt)
text(1.02, 0.76, "Carnivorous\n Copepods", cex = txt, adj = c(1,0))
text(1.01, 0.68, "Chaetognaths", cex = txt)
text(0.98, 0.6, "Jellyfish", cex = txt)
text(1.2, 0.25, "Planktivorous\n Fish", cex = txt)

text(-1,1,"E)",cex = 1.5, font = 2)

### EUTROPHIC FOOD WEB
# circos.par(gap.after = c(rep(5, nrow(diets_mat_eutro))))
chordDiagram(diets_mat_eutro,
             directional = 1,
             grid.col = model_color,
             direction.type = c("diffHeight"),
             link.arr.type = "triangle",
             diffHeight = -uh(2, "mm"),
             annotationTrack = c("grid"),
             annotationTrackHeight = c(0.1, 0.01),
             transparency = 0.2)


text(0.6, -1, "Phytoplankton", cex = txt)
text(-1.4, 0.75, "Omnivorous\n Copepods", cex = txt, adj = c(-1,0))
text(0.45, 1, "Euphausiids", cex = txt)
text(1.15, 0.35, "Planktivorous\n Fish", cex = txt)

arrows(x0 = -1.1, y0 = -0.6, x1 = -0.9, y1 = -0.5, length = 0, lwd = 2)
text(-1.2, -0.6, "Larvaceans\n Salps", cex = txt, pos = 1)

arrows(x0 = 0.79, y0 = 0.65, x1 = 0.85, y1 = 0.7, length = 0, lwd = 2)
text(0.98, 0.9, "Carnivorous Copepods\n Chaetognaths\n Jellyfish", cex = txt,  pos = 1)

text(-1,1,"F)",cex = 1.5, font = 2)

dev.print(pdf,'Figures/FoodWebs_3ef_ChordDiagram.pdf')
circos.clear()


#
#
#
#
#
#
#
#
#
# ### Do Chords with Hetero
# diets_mat_oligo = Chord_Prepare(diets, enviros_oligo, TestGroups)
# diets_mat_eutro = Chord_Prepare(diets, enviros_eutro, TestGroups)
#
# model_color = c("darkgreen", TestGroups$Colour[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)]) # Add green for phytoplankton
#
# x11(width = 10, height = 5)
# par(mfrow = c(1,2))
#
# ### OLIGOTROPHIC FOOD WEB
# chordDiagram(diets_mat_oligo,
#              directional = 1,
#              grid.col = model_color,
#              direction.type = c("diffHeight"),
#              link.arr.type = "triangle",
#              diffHeight = uh(0.8, "mm"),
#              annotationTrack = c("grid"),
#              annotationTrackHeight = c(0.1, 0.01),
#              transparency = 0.1,
#              preAllocateTracks = 1)
#
# for(si in get.all.sector.index()) {
#   if(get.cell.meta.data("xlim", sector.index = si)[2] > 0.02){
#     if (get.cell.meta.data("xlim", sector.index = si)[2] > 0.05) {
#       prop <- seq(0,1,0.2)} else{prop <- seq(0,1,0.5)}
#     circos.axis(h = "top",
#                 labels.cex = 0.5,
#                 sector.index = si,
#                 track.index = 2,
#                 labels.facing = "outside",
#                 major.at = (get.cell.meta.data("xlim", sector.index = si)[2] *prop),
#                 labels = (get.cell.meta.data("xlim", sector.index = si)[2] *prop)/(get.cell.meta.data("xlim", sector.index = si)[2])*100)
#
#     xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
#     ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
#     circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
#                 facing = "bending.inside", niceFacing = TRUE, col = "white", cex = 0.8)
#   }
# }
#
# ### EUTROPHIC FOOD WEB
# chordDiagram(diets_mat_eutro,
#              directional = 1,
#              grid.col = model_color,
#              direction.type = c("diffHeight"),
#              link.arr.type = "triangle",
#              diffHeight = uh(0.8, "mm"),
#              annotationTrack = c("grid"),
#              annotationTrackHeight = c(0.1, 0.01),
#              transparency = 0.1,
#              preAllocateTracks = 1)
#
# for(si in get.all.sector.index()) {
#   if(get.cell.meta.data("xlim", sector.index = si)[2] > 5){
#     if (get.cell.meta.data("xlim", sector.index = si)[2] > 2) {
#       prop <- seq(0,1,0.2)} else{prop <- seq(0,1,0.5)}
#     circos.axis(h = "top",
#                 labels.cex = 0.5,
#                 sector.index = si,
#                 track.index = 2,
#                 labels.facing = "outside",
#                 major.at = (get.cell.meta.data("xlim", sector.index = si)[2] *prop),
#                 labels = (get.cell.meta.data("xlim", sector.index = si)[2] *prop)/(get.cell.meta.data("xlim", sector.index = si)[2])*100)
#
#     xlim = get.cell.meta.data("xlim", sector.index = si, track.index = 2)
#     ylim = get.cell.meta.data("ylim", sector.index = si, track.index = 2)
#     circos.text(mean(xlim), mean(ylim), si, sector.index = si, track.index = 2,
#                 facing = "bending.inside", niceFacing = TRUE, col = "white", cex = 0.8)
#   }
# }
#
# dev.print(pdf,'Figures/FoodWebs_Supp_ChordDiagram.pdf')
# circos.clear()
#
