
# Choose environmental data to use
enviro_data <- readRDS("envirofull_20200312.RDS") # full set of 5x5 degree grids

quantile(enviro_data$chlo, c(0.1, 0.9))
quantile(enviro_data$chlo, c(0.2, 0.8))
quantile(enviro_data$chlo, c(0.25, 0.75))





# mollCRS <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84")
# latlonCRS <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# rasChl <- raster("~/GitHub/FoodWebEfficiency/A20021852019334.L3m_CU_CHL_chlor_a_9km.nc",
#                  varname = "chlor_a", nrows=4320, ncols=8640, xmn=-180, xmx=180, ymn=-90, ymx=90)
rasChl <- read_rds("~/GitHub/FoodWebEfficiency/SatelliteDataAnalysis/Chl5Deg.rds")
rasChl_poly <- rasterToPolygons(rasChl, fun=NULL, n=4, na.rm=TRUE, digits=2, dissolve=FALSE) # Convert to polygon which is better for plotting
rasChl_sf <- st_as_sf(rasChl_poly, crs = 4326) # Convert to sf
rasChl_sf_moll <- st_transform(rasChl_sf, crs = 54009) # Alter the CRS to mollweide
#

old <- read_rds("/Users/jason/Dropbox/Multi-Zoo Size Spectrum Model/_LatestModel/20200212_Control_Full_UNSW/envirofull_20200209.RDS")
rasChl <- old %>%
  dplyr::select(x, y, chlo) %>%
  rename("Chl" = chlo)

rasChl_raster <- rasterFromXYZ(rasChl)  #Convert first two columns as lon-lat and third as value
crs(rasChl_raster) <- CRS('+init=EPSG:4326')
rasChl_poly <- rasterToPolygons(rasChl_raster, fun=NULL, n=4, na.rm=TRUE, digits=2, dissolve=FALSE) # Convert to polygon which is better for plotting
rasChl_sf <- st_as_sf(rasChl_poly, crs = 4326) # Convert to sf
st_crs(rasChl_sf) <- 4326 # Set base projection
rasChl_sf_moll <- st_transform(rasChl_sf, crs = 54009) # Alter the CRS to mollweide

# across oligotrophic (which we define as areas where chlorophyll a < 0.1 mg m-3)
# to eutrophic (which we define as areas where chlorophyll a > 1 mg m-3) waters

minChl <- -1
maxChl <- 0
minChl <- -1
maxChl <- -0.3
rasChl_sf_moll <- rasChl_sf_moll %>%
  mutate(Chl2 = case_when(Chl > 10^maxChl ~ "High",
                          Chl < 10^(minChl) ~ "Low",
                          Chl > 10^(minChl) & Chl < 10^maxChl ~ "Mid"))


dat <- tibble(Chl = rasChl_sf_moll$Chl2)

rasChl_sf_moll$Chl2 <- factor(rasChl_sf_moll$Chl2, levels = c("High", "Mid", "Low"))

# What bin limits should I use for
ggplot(rasChl_sf_moll, aes(log10(Chl))) +
  geom_histogram(binwidth = 0.1)

rasChl_Summary <- dat %>%
  mutate(Chl2 = as.factor(Chl)) %>%
  group_by(Chl) %>%
  summarise(n = n())



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


gg_map <- ggplot() +
  geom_sf(data = rasChl_sf_moll, aes(fill = Chl2)) +
  geom_sf(data = world_sf, size = 0.05, fill = "grey50") +
  # geom_sf(data = sub_sf, size = 0.4, colour = "purple", fill = "purple", inherit.aes = T) +
  ggtitle("C) Data Distribution") +
  # scale_fill_gradient(name = element_blank(),
  #                     limits = c(log10(0.1), log10(1)),
  #                     low = "#ece2f0",
  #                     high = "darkgreen",
  #                     na.value = "grey50",
  #                     guide = "colourbar",
  #                     aesthetics = "fill",
  #                     oob = scales::squish) +
  theme_opts +
  theme(panel.grid.major = element_line(colour = "grey50", size = 0.2))


