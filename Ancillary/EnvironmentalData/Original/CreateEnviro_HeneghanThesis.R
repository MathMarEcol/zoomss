rm(list = ls())

library(ncdf4) # For reading/manipulating netCDFs
library(raster)
library(sp)

### IMPORT BATHYMETRY - USE AS LAND MASK, ASSIGN ALL LAND VALUES AS NA
bathy_f <- list.files(pattern = "GEBCO", recursive = TRUE)
bathy_file <- read.csv(bathy_f, header = FALSE)
bathy_filee <- as.vector(t(as.matrix(bathy_file)))
bathy_store = matrix(NA, nrow = length(bathy_filee), ncol = 3)
bathy_store[,c(1,2)] = as.matrix(expand.grid("Long" = seq(-179.5, 179.5, 1), "Lat" = seq(89.5, -89.5, -1)))
bathy_store[,3] <- bathy_filee
colnames(bathy_store) <- c("Long", "Lat", "BATHY")
bathy_store <- as.data.frame(bathy_store)
bathy_store[c(bathy_store$BATHY >= 0), "BATHY"] <- NA
#land <- c(bathy_filee > 0) # which cells are land

### Import chlorophyll and sst data
setwd("C:/Users/1539391/Desktop/enviro_5d_create")

sst_files <- list.files(pattern = "MC_SST", recursive = TRUE)
chlo_files <- list.files(pattern = "MC_CHL", recursive = TRUE)

sst_files <- sst_files[c(7:12, 1:6)] # Reorder files so Jan-June go first
chlo_files <- chlo_files[c(7:12, 1:6)] 

# IMPORT NETCDF AND EXTRACT GEOGRAPHIC REFERENCES
nc.sst <- nc_open(sst_files[1])

sst_lat <- ncvar_get(nc.sst, "lat")
sst_lon <- ncvar_get(nc.sst, "lon")
sst_store <- chlo_store <- matrix(NA, nrow = length(sst_lat)*length(sst_lon), ncol = 14)
sst_store[,c(1,2)] <- chlo_store[,c(1,2)] <- as.matrix(expand.grid(lon = sst_lon, lat = sst_lat))

### OPEN NETCDFS
sst_nc_all <- lapply(sst_files, nc_open)
chlo_nc_all <- lapply(chlo_files, nc_open)

## EXTRACT SST AND CHLOROPHYLL
sst_nc_vall <- lapply(sst_nc_all, ncvar_get, varid = "sst4")
chlo_nc_vall <- lapply(chlo_nc_all, ncvar_get, varid = "chl_ocx")

## PUT MONTHLY SST AND CHLO IN RESPECTIVE MATRICES, LABEL COLUMNS AND CONVERT TO DATA FRAME
for(i in 1:12){
  sst_store[,c(i+2)] <- as.vector((unlist(sst_nc_vall[[i]])))
  chlo_store[,c(i+2)] <- as.vector((unlist(chlo_nc_vall[[i]])))
}

colnames(sst_store) <- colnames(chlo_store) <- c("Long","Lat","January", "February", "March", "April", "May", 
                              "June", "July", "August", "September", "October", "November", "December")

sst_store <- as.data.frame(sst_store)
chlo_store <- as.data.frame(chlo_store)

####### CALCULATE ANNUAL AVERAGES
# Here, we're assuming grid squares that are missing coverage in certain months are 0 for chlo and temp (no/little light in those months at those latitudes)
sst_store$annual_ave <- rowSums(sst_store[,3:14], na.rm = TRUE)/12
chlo_store$annual_ave <- rowSums(chlo_store[,3:14], na.rm = TRUE)/12 

###### AGGREGATE TO 5x5 GRID SQUARES, THIS TAKES A WHILE
sst_store <- sst_store[,c(1,2,15)] # we only need annual ave
chlo_store <- chlo_store[,c(1,2,15)]

# 5x5 degree grid squares globally
lon <- seq(-177.5,177.5,5)
lat <- seq(-87.5,87.5,5)
lonlat <- expand.grid(Long = lon, Lat = lat)

# Give coordinates to lonlat
coordinates(lonlat) <- ~Long+Lat
crs(lonlat) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Convert lonlat to grid
grd_lonlat <- points2grid(lonlat)
grd_lonlat <- SpatialGrid(grd_lonlat, proj4string = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

### Aggregate sst to 5x5 lonlat grid
coordinates(sst_store) <- ~ Long + Lat
crs(sst_store) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
sst_store5 <- raster(aggregate(sst_store, grd_lonlat, mean, na.rm = TRUE))
sst_store5d <- as.data.frame(sst_store5, xy = TRUE)

### Aggregate CHLO to 5x5 lonlat grid
coordinates(chlo_store) <- ~ Long + Lat
crs(chlo_store) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
chlo_store5 <- raster(aggregate(chlo_store, grd_lonlat, mean, na.rm = TRUE))
chlo_store5d <- as.data.frame(chlo_store5, xy = TRUE)

### CREATE LANDMASK: Aggregate BATHY to 5x5 lonlat grid
coordinates(bathy_store) <- ~ Long + Lat
crs(bathy_store) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
bathy_store5 <- raster(aggregate(bathy_store, grd_lonlat, mean, na.rm = TRUE))
bathy_store5d <- as.data.frame(bathy_store5, xy = TRUE)


# Apply land masks to sst and chlo 5x5 global grid
chlo_store5d[which(is.na(bathy_store5d$BATHY) == TRUE), 'annual_ave'] <- NA
sst_store5d[which(is.na(bathy_store5d$BATHY) == TRUE), 'annual_ave'] <- NA


# Create and save enviro_data csv
enviro_data <- data.frame('x' = sst_store5d$x, 'y' = sst_store5d$y, 'sst' = sst_store5d$annual_ave, 'chlo' = chlo_store5d$annual_ave)
write.csv(enviro_data, 'enviro_data.csv', row.names = FALSE)
