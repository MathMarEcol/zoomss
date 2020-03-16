library(tidyverse)
library(ncdf4)
library(raster)

res = "five"

# 5 x 5 degree grid squares globally
if(str_detect(res, "five")){
  lon <- seq(-177.5,177.5,5)
  lat <- seq(-87.5,87.5,5)
}

# 1 x 1 degree grid squares globally
if(str_detect(res, "one")){
  lon <- seq(-179.5,179.5,1)
  lat <- seq(-89.5,89.5,1)
}

Coord_ref <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
lonlat <- expand.grid(Long = lon, Lat = lat)
coordinates(lonlat) <- ~Long+Lat
proj4string(lonlat) <- CRS(Coord_ref)
gridded(lonlat) = TRUE
grd_lonlat <- raster(lonlat)

##### Bathy #####
bathy_file <- paste0("Data",.Platform$file.sep,"GEBCO_2014_2D.nc")
rasBathy <- raster(bathy_file, varname = "elevation")
Bathy <- resample(rasBathy, grd_lonlat, method='bilinear')
Bathy <- reclassify(Bathy, c(0,Inf,NA))
Bathy <- abs(Bathy)
names(Bathy) <- "Bathy"
saveRDS(Bathy, file = paste0("Data",.Platform$file.sep,"Bathy_raster_",res,"Deg.rds"))

##### SST #####
SST_files <- list.files(path = "Data", pattern = "_sst4_", recursive = FALSE, full.names = TRUE)

for(i in 1:12){
  rasSST <- raster(SST_files[i], varname = "sst4", nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90)
  SSTagg <- resample(rasSST, grd_lonlat, method='bilinear')
  names(SSTagg) <- "SST"
  if (i == 1){
    SSTbrick <- brick(SSTagg)} else{
      SSTbrick <- addLayer(SSTbrick,SSTagg)
    }
}

SST_mn  <- calc(SSTbrick, fun = mean, na.rm = TRUE)
names(SST_mn) <- "SST"
saveRDS(SSTbrick_mn, file = paste0("SST_raster_00Mean_",res,"Deg.rds"))

#### Chlorophyll ####

Chl_files <- paste0("Data",.Platform$file.sep,"A20021852019334.L3m_CU_CHL_chl_ocx_9km.nc")
rasChl <- raster(Chl_files, varname = "chl_ocx", nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90)
Chl_mn <- resample(rasChl, grd_lonlat, method='bilinear')
names(Chl_mn) <- "Chl"
saveRDS(Chl_mn, file = paste0("Chl_raster_00Mean_",res,"Deg.rds"))


# Now convert to data-frame and save as .RDS
brick <- brick(Bathy)
brick <- addLayer(brick, SST_mn, Chl_mn)

df <- as.data.frame(brick, xy=TRUE)

df <- df %>% 
  drop_na



