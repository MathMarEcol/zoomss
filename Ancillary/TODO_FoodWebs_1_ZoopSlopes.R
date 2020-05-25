## script to take the outputs and produce summaries of abundance and biomass by species and weight class

library(tidyverse)
library(patchwork)

library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(raster)

source("FoodWebs_0a_Functions.R")
source("FoodWebs_0b_Initialise.R")

enviro_data <- read_rds("~/Nextcloud/MME2Work/ZooMSS/_LatestModel/20200221_OneZoo_Full_UNSW/envirofull_20200312.RDS") # Load in environmental data
res <- read_rds("~/Nextcloud/MME2Work/ZooMSS/_LatestModel/20200212_Control_Full_UNSW/Output/res_20200212_Control_Full_UNSW.RDS")

txt_size <- 4.5

# Original
oligo <- c(0,0.1)
eutro <- c(1,100)

Bins <- read_csv("~/Dropbox/MissingLink/Output/LOPC_Bins_250.csv", col_names = FALSE)
Limits <- read_csv("~/Dropbox/MissingLink/Output/LOPC_Limits_250.csv", col_names = FALSE)

### Plot data
plot_params <- tibble(dx = 0.1, # log10 weight step
                      w0 = 10^(min(TestGroups$W0)), # minimum dynamic size class
                      wMax = 10^(max(TestGroups$Wmax)), # maximum dynamic size class
                      MinSize = log10(Limits$X1[1]), # Min zoo size
                      MaxSize = log10(Limits$X1[length(Limits$X1)]), # Max zoo size
                      XMin = -2.5, # Plot x-min
                      XMax = 2.7, # Plot x-max
                      YMin = -4,#-1.2, # Plot y-min
                      YMax = 4) # Plot y-max
# Set weight bins
w = 10^(seq(from = log10(plot_params$w0), to =  log10(plot_params$wMax), plot_params$dx))

### Load OPC Data
LOPC <- read_csv("~/Dropbox/MissingLink/Output/MasterDatabase_v20200308.csv",
                 col_types = cols(Voyage = "f", VolMethod = "f", User = "f", Unit = "f",
                                  Season = "f",Model = "f", Region = "f", LonghurstBiome = "f",
                                  LonghurstProvince = "f",OceanName = "f", Orientation = "f"))

LOPC <- LOPC %>%
  filter(is.na(Chl)==0) %>%
  filter(Depth < 200)


## NOW DO MODEL VS OBSERVATIONS #

## Data fits
NBSSdata_list_oligo <- NBSS_data(LOPC, Bins, oligo)
NBSSdata_oligo <- NBSSdata_list_oligo[[1]]
NBSSdata_full_oligo <- NBSSdata_list_oligo[[2]]
NBSSdata_full_oligo$Enviro <- "Oligotrophic"

NBSSdata_list_eutro <- NBSS_data(LOPC, Bins, eutro)
NBSSdata_eutro <- NBSSdata_list_eutro[[1]]
NBSSdata_full_eutro <- NBSSdata_list_eutro[[2]]
NBSSdata_full_eutro$Enviro <- "Eutrophic"

NBSSdata_all <- bind_rows(NBSSdata_full_eutro, NBSSdata_full_oligo)
NBSSdata_all$Enviro <- as.factor(NBSSdata_all$Enviro)

modelfit_LOPC <- lm(log10(NB) ~ log10(Bins) + Enviro, data = NBSSdata_all) # No interactions
summary(modelfit_LOPC)

modelfit_LOPC_int <- lm(log10(NB) ~ log10(Bins) * Enviro, data = NBSSdata_all) # Add in the interaction
summary(modelfit_LOPC_int)
anova(modelfit_LOPC_int, modelfit_LOPC)

plot(modelfit_LOPC_int)


m1<- glm(NB ~ log10(Bins) + Enviro, family = Gamma(link = "log"), data = NBSSdata_all) # Add in the interaction
m2 <- glm(NB ~ log10(Bins) * Enviro, family = Gamma(link = "log"), data = NBSSdata_all) # Add in the interaction
anova(m1, m2, test = "F")
plot(m2)



## ZooMSS fits

NBSS_ZooMSS_oligo <- NBSS_ZooMSS(res, oligo, enviro_data, w, plot_params) # Compile Oligotrophic Model Output
modelfit_ZooMSS_oligo <- lm(log10(NB) ~ log10(Bins), data = NBSS_ZooMSS_oligo) # Do NBSS Fit

NBSS_ZooMSS_eutro <- NBSS_ZooMSS(res, eutro, enviro_data, w, plot_params) # Compile Eutrophic Model Output
modelfit_ZooMSS_eutro <- lm(log10(NB) ~ log10(Bins), data = NBSS_ZooMSS_eutro) # Do NBSS Fit

NBSS_ZooMSS_oligo$Enviro <- "Oligotrophic" # Add factor for bound df below.
NBSS_ZooMSS_eutro$Enviro <- "Eutrophic"

NBSS_ZooMSS_all <- bind_rows(NBSS_ZooMSS_eutro, NBSS_ZooMSS_oligo)
NBSS_ZooMSS_all$Enviro <- as.factor(NBSS_ZooMSS_all$Enviro)

##
gg_dat <- ggplot() +
  scale_x_continuous(limits = c(plot_params$XMin, plot_params$XMax), expand = c(0, 0)) +
  # scale_y_continuous(limits = c(plot_params$YMin, plot_params$YMax), expand = c(0, 0)) +
  geom_point(data = NBSSdata_oligo, mapping = aes(x = log10(Bins), y = log10(NB)), colour = "red3") +
  annotate("segment", x = plot_params$MinSize, xend = plot_params$MaxSize,
           y = modelfit_LOPC_oligo$coefficients[1] + modelfit_LOPC_oligo$coefficients[2]*plot_params$MinSize,
           yend = modelfit_LOPC_oligo$coefficients[1] + modelfit_LOPC_oligo$coefficients[2]*plot_params$MaxSize,
           colour = "red3", size = 1.5) +
  ylab(expression(paste("log"[10], "(Normalised biomass (" * m^-3 * "))"))) +
  xlab(expression(paste("log"[10], "(Zooplankton size class (mg))"))) +
  annotate("text", x = -1, y = 0.5,
           label = "Oligotrophic",
           colour = "red3", size = txt_size, hjust = "centre") +
  annotate("text", x = -1, y = 0,
           label = paste("Slope == ",sprintf("%.2f", modelfit_LOPC_oligo$coefficients[2])),
           parse = TRUE, colour = "red3", size = txt_size, hjust = "centre") +
  geom_point(data = NBSSdata_eutro, mapping = aes(x = log10(Bins), y = log10(NB)), colour = "blue3") +
  annotate("segment", x = plot_params$MinSize, xend = plot_params$MaxSize,
           y = modelfit_LOPC_eutro$coefficients[1] + modelfit_LOPC_eutro$coefficients[2]*plot_params$MinSize,
           yend = modelfit_LOPC_eutro$coefficients[1] + modelfit_LOPC_eutro$coefficients[2]*plot_params$MaxSize,
           colour = "blue3", size = 1.5) +
  annotate("text", x = 0.5, y = 3.5,
           label = "Eutrophic",
           colour = "blue3", size = txt_size, hjust = "centre") +
  annotate("text", x = 0.5, y = 3,
           label = paste("Slope == ",deparse(sprintf("%.2f", round(modelfit_LOPC_eutro$coefficients[2], digits = 2)))),
           parse = TRUE, colour = "blue3", size = txt_size, hjust = "centre") +
  theme_bw() +
  ggtitle(label = "",subtitle = "D) Data") +
  theme(plot.title = element_text(size = 18),
        plot.margin = unit(c(0,0,0,0), "mm"))

gg_mod <- ggplot() +
  scale_x_continuous(limits = c(plot_params$XMin, plot_params$XMax), expand = c(0, 0)) +
  # scale_y_continuous(limits = c(plot_params$YMin, plot_params$YMax), expand = c(0, 0)) +
  geom_point(data = NBSS_ZooMSS_oligo, mapping = aes(x = log10(Bins), y = log10(NB)), colour = "red3") +
  annotate("segment", x = plot_params$MinSize, xend = plot_params$MaxSize,
           y = modelfit_ZooMSS_oligo$coefficients[1] + modelfit_ZooMSS_oligo$coefficients[2]*plot_params$MinSize,
           yend = modelfit_ZooMSS_oligo$coefficients[1] + modelfit_ZooMSS_oligo$coefficients[2]*plot_params$MaxSize,
           colour = "red3", size = 1.5) +
  ylab(expression(paste("log"[10], "(Normalised biomass (" * m^-3 * "))"))) +
  xlab(expression(paste("log"[10], "(Zooplankton size class (mg))"))) +
  annotate("text", x = -0.5, y = -0.5,
           label = "Oligotrophic",
           colour = "red3", size = txt_size, hjust = "centre") +
  annotate("text", x = -0.5, y = -1,
           label = paste("Slope == ",sprintf("%.2f", round(modelfit_ZooMSS_oligo$coefficients[2], digits = 2))),
           parse = TRUE, colour = "red3", size = txt_size, hjust = "centre") +
  geom_point(data = NBSS_ZooMSS_eutro, mapping = aes(x = log10(Bins), y = log10(NB)), colour = "blue3") +
  annotate("segment", x = plot_params$MinSize, xend = plot_params$MaxSize,
           y = modelfit_ZooMSS_eutro$coefficients[1] + modelfit_ZooMSS_eutro$coefficients[2]*plot_params$MinSize,
           yend = modelfit_ZooMSS_eutro$coefficients[1] + modelfit_ZooMSS_eutro$coefficients[2]*plot_params$MaxSize,
           colour = "blue3", size = 1.5) +
  annotate("text", x = 0.5, y = 2.5,
           label = "Eutrophic",
           colour = "blue3", size = txt_size, hjust = "centre") +
  annotate("text", x = 0.5, y = 2,
           label = paste("Slope == ",sprintf("%.2f", round(modelfit_ZooMSS_eutro$coefficients[2], digits = 2))),
           parse = TRUE, colour = "blue3", size = txt_size, hjust = "centre") +
  theme_bw() +
  ggtitle(label = "",subtitle = "E) Model") +
  theme(plot.title = element_text(size = 18),
        plot.margin = unit(c(0,0,0,0), "mm"))


world <- ne_countries(scale = "medium", returnclass = "sf")
world_sf <- st_transform(world, crs = 54009) # Convert to Mollweide

theme_opts <- list(theme(
  # panel.grid.minor = element_blank(),
  # panel.grid.major = element_blank(),
  panel.background = element_blank(),
  # panel.border = element_blank(),
  # plot.background = element_rect(fill = NA),
  plot.margin = unit(c(0,0,0,0), "mm"),
  # axis.line = element_blank(),
  # axis.text.x = element_blank(),
  # axis.text.y = element_blank(),
  # axis.ticks = element_blank(),
  # axis.title.x = element_blank(),
  # axis.title.y = element_blank(),
  # legend.title = element_text(size = 6),
  # legend.text = element_text(size = 6),
  # legend.position = "right",
  # legend.direction = "horizontal",
  # legend.background = element_rect(fill = NA),
  # legend.key.height = unit(9, "mm"),
  # legend.key.width = unit(4, "mm"),
  # legend.position = c(0.5, -0.05),
))

dat <- LOPC %>%
  dplyr::select(c(Lon, Lat))
sub_sf <- st_as_sf(dat, coords = c('Lon', 'Lat'), crs = 4326) # Convert to sf
sub_sf <- st_transform(sub_sf, crs = 54009) # Alter the CRS
#
# nb <- readRDS("../GlobalZoopAbundModels/ModelOutput/GlobalLayers/Num_Brick_Biomass.rds")
# Chl_poly <- rasterToPolygons(nb$Chl, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) # Convert to polygon which is better for plotting
# Chl <- st_as_sf(Chl_poly, coords = c("x", "y"))
# st_crs(Chl) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
# Chl <- st_transform(Chl, crs = 54009) # Alter the CRS

enviro_data <- read_rds("~/GitHub/ZooMSS/envirofull_20200317.RDS")
Chl_data <- enviro_data %>%
  dplyr::select(Lon, Lat, chlo) %>%
  rename("Lon" = Lon, "Lat" = Lat, "Chl" = chlo)

Chl_raster <- rasterFromXYZ(Chl_data)  #Convert first two columns as lon-lat and third as value
crs(Chl_raster) <- CRS('+init=EPSG:4326')
Chl_poly <- rasterToPolygons(Chl_raster, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) # Convert to polygon which is better for plotting
Chl_sf <- st_as_sf(Chl_poly, crs = 4326) # Convert to sf
st_crs(Chl_sf) <- 4326 # Set base projection
Chl_sf_moll <- st_transform(Chl_sf, crs = 54009) # Alter the CRS to mollweide


gg_map <- ggplot() +
  geom_sf(data = Chl_sf_moll, mapping = aes(fill = log10(Chl)), colour = NA) +
  geom_sf(data = sub_sf, size = 0.4, colour = "purple", fill = "purple", inherit.aes = T) +
  geom_sf(data = world_sf, size = 0.05, fill = "grey50") +
  ggtitle("C) Data Distribution") +
  scale_fill_gradient(name = element_blank(),
                      limits = c(-1, 0.5),
                      low = "#e5f5e0", #"#c7e9c0",
                      high = "#00441b",
                      na.value = "grey50",
                      guide = guide_colourbar(title = expression(paste("log"[10], "(Chlorophyll ",italic(a),"\n mg m"^-3, ")")),
                                              title.position = "right",
                                              title.vjust = 0.5,
                                              title.hjust = 0.5,
                                              title.theme = element_text(angle = 270, size = 10),
                                              barheight = 10),
                      aesthetics = "fill",
                      oob = scales::squish) +
  theme(plot.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "grey50", size = 0.2),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank(),
        plot.margin = unit(c(0,0,0,0), "mm"))

graphics.off()
x11(height = 7, width = 8)

layout <- "
AA
AA
AA
AA
BC
BC
BC
"

gg_map + gg_dat + gg_mod + plot_layout(design = layout)
ggsave("Figures/FoodWebs_1_Slopes.png", dpi = 400)



## Do the statistics

load("NBSS_Data.RData")

modelfit_LOPC <- lm(log10(NB) ~ log10(Bins) + Enviro, data = NBSSdata_all) # No interactions
summary(modelfit_LOPC)

modelfit_LOPC_int <- lm(log10(NB) ~ log10(Bins) * Enviro, data = NBSSdata_all) # Add in the interaction
summary(modelfit_LOPC_int)

anova(modelfit_LOPC_int) # ? ANCOVA


modelfit_ZooMSS <- lm(log10(NB) ~ log10(Bins) + Enviro, data = NBSS_ZooMSS_all) # No interaction
summary(modelfit_ZooMSS)

modelfit_ZooMSS_int <- lm(log10(NB) ~ log10(Bins) * Enviro, data = NBSS_ZooMSS_all) # With interaction
summary(modelfit_ZooMSS_int)

anova(modelfit_ZooMSS_int) # ? ANCOVA




