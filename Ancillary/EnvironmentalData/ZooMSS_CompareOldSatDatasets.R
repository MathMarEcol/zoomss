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
old <- read_csv('/Users/jason/Dropbox/Multi-Zoo Size Spectrum Model/_LatestModel/20191106full/enviro_5d_3.csv')
old <- old %>%
  dplyr::select(x, y, chlo) %>%
  rename("Lon" = x, "Lat" = y, "Chl" = chlo)

old_raster <- rasterFromXYZ(old)  #Convert first two columns as lon-lat and third as value
crs(old_raster) <- CRS('+init=EPSG:4326')
old_poly <- rasterToPolygons(old_raster, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) # Convert to polygon which is better for plotting
old_sf <- st_as_sf(old_poly, crs = 4326) # Convert to sf
st_crs(old_sf) <- 4326 # Set base projection
old_sf_moll <- st_transform(old_sf, crs = 54009) # Alter the CRS to mollweide



## Now Do SST
## Old environmental data sets
oldSST <- read_csv('/Users/jason/Dropbox/Multi-Zoo Size Spectrum Model/_LatestModel/20191106full/enviro_5d_3.csv')
oldSST <- oldSST %>%
  dplyr::select(x, y, sst) %>%
  rename("Lon" = x, "Lat" = y, "SST" = sst)

oldSST_raster <- rasterFromXYZ(oldSST)  #Convert first two columns as lon-lat and third as value
crs(oldSST_raster) <- CRS('+init=EPSG:4326')
oldSST_poly <- rasterToPolygons(oldSST_raster, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) # Convert to polygon which is better for plotting
oldSST_sf <- st_as_sf(oldSST_poly, crs = 4326) # Convert to sf
st_crs(oldSST_sf) <- 4326 # Set base projection
oldSST_sf_moll <- st_transform(oldSST_sf, crs = 54009) # Alter the CRS to mollweide



##### ##### Create the new resolution to match the SSM ##### #####
# 5 x 5 degree grid squares globally
lon <- seq(-177.5,177.5,5)
lat <- seq(-87.5,87.5,5)

# 1 x 1 degree grid squares globally
# lon <- seq(-179.5,179.5,1)
# lat <- seq(-89.5,89.5,1)

# 0.1 x 0.1 degree
# lon <- seq(-179.5,179.5,0.1)
# lat <- seq(-89.5,89.5,0.1)

Coord_ref <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
lonlat <- expand.grid(Long = lon, Lat = lat)
coordinates(lonlat) <- ~Long+Lat
proj4string(lonlat) <- CRS('+init=EPSG:4326')
gridded(lonlat) = TRUE
grd_lonlat <- raster(lonlat)

# There seems to be no cumulative SST dataset. Not sure why? For the moment, I just average the monthly ones.

##### ##### Now do SST ##### #####
SST_files <- list.files(path = "~/GitHub/GlobalZoopAbundModels/EnvironmentalData/nc",
                        pattern = "MC_SST", recursive = FALSE, full.names = TRUE)

for(i in 1:12){
  rasSST <- raster(SST_files[i], varname = "sst", nrows=4320, ncols=8640, xmn=-180, xmx=180, ymn=-90, ymx=90)
  names(rasSST) <- "SST"
  if (i == 1){
    SSTbrick <- brick(rasSST)} else{
      SSTbrick <- addLayer(SSTbrick,rasSST)
    }
}

rasSST  <- calc(SSTbrick, fun = function(x){sum(x, na.rm = TRUE)/12})
names(rasSST) <- "SST"
saveRDS(rasSST, file = "SST_raster.rds")

SST5Deg <- resample(rasSST, grd_lonlat, method='bilinear')
names(SST5Deg) <- "SST"
saveRDS(SST5Deg, file = "SST5Deg.rds")


##### ##### Now do Chl ##### #####

# rasChl <- raster("/Users/jason/GitHub/GlobalZoopAbundModels/EnvironmentalData/nc/A20021852019304.L3m_CU_CHL_chlor_a_4km.nc",
#                  varname = "chlor_a", nrows=4320, ncols=8640, xmn=-180, xmx=180, ymn=-90, ymx=90)

Chl_files <- list.files(path = "~/GitHub/GlobalZoopAbundModels/EnvironmentalData/nc",
                        pattern = "MC_CHL", recursive = FALSE, full.names = TRUE)

for(i in 1:12){
  rasChl <- raster(Chl_files[i], varname = "chlor_a", nrows=4320, ncols=8640, xmn=-180, xmx=180, ymn=-90, ymx=90)
  names(rasChl) <- "Chl"
  if (i == 1){
    Chlbrick <- brick(rasChl)} else{
      Chlbrick <- addLayer(Chlbrick,rasChl)
    }
}

rasChl  <- calc(Chlbrick, fun = function(x){sum(x, na.rm = TRUE)/12})
names(rasChl) <- "Chl"
saveRDS(rasChl, file = "Chl_raster.rds")

rasChl[rasChl>10] <- 10

Chl5Deg <- resample(rasChl, grd_lonlat, method='bilinear')
names(Chl5Deg) <- "Chl"
saveRDS(Chl5Deg, file = "Chl5Deg.rds")



## New environmental data sets
new_raster <- Chl5Deg
new_poly <- rasterToPolygons(new_raster, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) # Convert to polygon which is better for plotting
new_sf <- st_as_sf(new_poly, crs = 4326) # Convert to sf
new_sf_moll <- st_transform(new_sf, crs = 54009) # Alter the CRS to mollweide

## Difference
# Get the difference between the rasters
diff_raster <- ((new_raster - old_raster)/old_raster)*100
names(diff_raster) <- "Chl_diff"
diff_poly <- rasterToPolygons(diff_raster, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) # Convert to polygon which is better for plotting
diff_sf <- st_as_sf(diff_poly) # Convert to sf
diff_sf_moll <- st_transform(diff_sf, crs = 54009) # Alter the CRS to mollweide

gg_oldChl <- ggplot() +
  geom_sf(data = old_sf_moll, aes(fill = log10(Chl)), color = NA) +
  geom_sf(data = world_sf, size = 0.05, fill = "grey20") +
  scale_fill_gradient(name = element_blank(),
                      limits = c(-1, 0.5),
                      low = "white",
                      high = "darkgreen",
                      na.value = "grey50",
                      guide = "colourbar",
                      aesthetics = "fill",
                      oob = scales::squish) +
  theme(plot.background = element_rect(fill = NA),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()) +
  ggtitle(label = "",subtitle = "Original Satellite (Chl. a)")

gg_newChl <- ggplot() +
  geom_sf(data = new_sf_moll, aes(fill = log10(Chl)), color = NA) +
  geom_sf(data = world_sf, size = 0.05, fill = "grey20") +
  scale_fill_gradient(name = element_blank(),
                      limits = c(-1, 0.5),
                      low = "white",
                      high = "darkgreen",
                      na.value = "grey50",
                      guide = "colourbar",
                      aesthetics = "fill",
                      oob = scales::squish) +
  theme(plot.background = element_rect(fill = NA),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()) +
  ggtitle(label = "",subtitle = "New Satellite (Chl. a)")

gg_diffChl <- ggplot() +
  geom_sf(data = diff_sf_moll, aes(fill = Chl_diff), color = NA) +
  geom_sf(data = world_sf, size = 0.05, fill = "grey20") +
  scale_fill_gradient2(name = element_blank(),
                       limits = c(-99, 99),
                       low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0,
                       na.value = "grey50",
                       guide = "colourbar",
                       aesthetics = "fill",
                       oob = scales::squish) +
  theme(plot.background = element_rect(fill = NA),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()) +
  ggtitle(label = "",subtitle = "Difference (%)")




newSST_raster <- SST5Deg
newSST_poly <- rasterToPolygons(newSST_raster, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) # Convert to polygon which is better for plotting
newSST_sf <- st_as_sf(newSST_poly, crs = 4326) # Convert to sf
newSST_sf_moll <- st_transform(newSST_sf, crs = 54009) # Alter the CRS to mollweide

## Difference
# Get the difference between the rasters
diff_raster <- ((newSST_raster - oldSST_raster)/oldSST_raster)*100
names(diff_raster) <- "SST_diff"
diff_poly <- rasterToPolygons(diff_raster, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) # Convert to polygon which is better for plotting
diff_sf <- st_as_sf(diff_poly) # Convert to sf
diff_sf_moll <- st_transform(diff_sf, crs = 54009) # Alter the CRS to mollweide

gg_oldSST <- ggplot() +
  geom_sf(data = oldSST_sf_moll, aes(fill = SST), colour = NA) +
  geom_sf(data = world_sf, size = 0.05, fill = "grey20") +
  scale_fill_gradient(name = element_blank(),
                      limits = c(-1, 30),
                      low = "white",
                      high = "darkgreen",
                      na.value = "grey50",
                      guide = "colourbar",
                      aesthetics = "fill",
                      oob = scales::squish) +
  theme(plot.background = element_rect(fill = NA),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()) +
  ggtitle(label = "",subtitle = "Original Satellite (SST)")

gg_newSST <- ggplot() +
  geom_sf(data = newSST_sf_moll, aes(fill = SST), color = NA) +
  geom_sf(data = world_sf, size = 0.05, fill = "grey20") +
  scale_fill_gradient(name = element_blank(),
                      limits = c(-1, 30),
                      low = "white",
                      high = "darkgreen",
                      na.value = "grey50",
                      guide = "colourbar",
                      aesthetics = "fill",
                      oob = scales::squish) +
  theme(plot.background = element_rect(fill = NA),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()) +
  ggtitle(label = "",subtitle = "New Satellite (SST)")

gg_diffSST <- ggplot() +
  geom_sf(data = diff_sf_moll, aes(fill = SST_diff), color = NA) +
  geom_sf(data = world_sf, size = 0.05, fill = "grey20") +
  scale_fill_gradient2(name = element_blank(),
                       limits = c(-99, 99),
                       low = "blue",
                       mid = "white",
                       high = "red",
                       midpoint = 0,
                       na.value = "grey50",
                       guide = "colourbar",
                       aesthetics = "fill",
                       oob = scales::squish) +
  theme(plot.background = element_rect(fill = NA),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()) +
  ggtitle(label = "",subtitle = "Difference (%)")





x11(height = 10, width = 10)
(gg_oldSST + gg_oldChl) / (gg_newSST + gg_newChl) / (gg_diffSST + gg_diffChl) / plot_annotation(tag_levels = "A", tag_suffix = ")")
ggsave(filename="Figures/CompareSatelliteData.png", dpi = 400)







