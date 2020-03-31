# This code loads and plots the data from the old and new satellite datasets to compare.
#
# I have discovered some differences between my code and Ryans code,
#   1. I use newer data
#   2. I use 4 km data, not 9 km data
#   3. I use 'resample' with bilinear methot from raster, while Ryan uses 'aggregate' with a mean.
#
# Written by Jason D. Everett (January 2020)

library(tidyverse)
library(patchwork)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(raster)

world <- ne_countries(scale = "medium", returnclass = "sf")
world_sf <- st_transform(world, crs = 54009) # Convert to Mollweide

## Old environmental data sets
enviro_data <- read_rds("envirofull_20200309.RDS")

Chl_data <- enviro_data %>%
  dplyr::select(x, y, chlo) %>%
  rename("Lon" = x, "Lat" = y, "Chl" = chlo)

Chl_raster <- rasterFromXYZ(Chl_data)  #Convert first two columns as lon-lat and third as value
crs(Chl_raster) <- CRS('+init=EPSG:4326')
Chl_poly <- rasterToPolygons(Chl_raster, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) # Convert to polygon which is better for plotting
Chl_sf <- st_as_sf(Chl_poly, crs = 4326) # Convert to sf
st_crs(Chl_sf) <- 4326 # Set base projection
Chl_sf_moll <- st_transform(Chl_sf, crs = 54009) # Alter the CRS to mollweide


SST_data <- enviro_data %>%
  dplyr::select(x, y, sst) %>%
  rename("Lon" = x, "Lat" = y, "SST" = sst)

SST_raster <- rasterFromXYZ(SST_data)  #Convert first two columns as lon-lat and third as value
crs(SST_raster) <- CRS('+init=EPSG:4326')
SST_poly <- rasterToPolygons(SST_raster, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) # Convert to polygon which is better for plotting
SST_sf <- st_as_sf(SST_poly, crs = 4326) # Convert to sf
st_crs(SST_sf) <- 4326 # Set base projection
SST_sf_moll <- st_transform(SST_sf, crs = 54009) # Alter the CRS to mollweide




gg_SST <- ggplot() +
  geom_sf(data = SST_sf_moll, aes(fill = SST), color = NA) +
  geom_sf(data = world_sf, size = 0.05, fill = "grey20") +
  scale_fill_distiller(palette = "Spectral",
                       name = element_blank(),
                       limits = c(1, 30),
                       na.value = "grey50",
                       guide = "colourbar",
                       aesthetics = "fill",
                       oob = scales::squish) +
  theme(plot.background = element_rect(fill = NA),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()) +
  ggtitle(label = "",subtitle = "Satellite SST")


gg_Chl <- ggplot() +
  geom_sf(data = Chl_sf_moll, aes(fill = log10(Chl)), color = NA) +
  geom_sf(data = world_sf, size = 0.05, fill = "grey20") +
  scale_fill_distiller(palette = "Spectral",
                       name = element_blank(),
                       limits = c(-1, 1),
                       na.value = "grey50",
                       guide = "colourbar",
                       aesthetics = "fill",
                       oob = scales::squish) +
  theme(plot.background = element_rect(fill = NA),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()) +
  ggtitle(label = "",subtitle = "Satellite Chl")



x11(height = 5, width = 10)
gg_SST + gg_Chl + plot_annotation(tag_levels = "A", tag_suffix = ")")
ggsave(filename="ZooMSS_SatelliteData.png", dpi = 400)



