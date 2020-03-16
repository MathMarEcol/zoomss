

## Lets compare the NASA CU with my arithmatic average of the monthly data

## I don't expect them to match because I assume NASA will use the daily raw data to create the climatology,
## and so I am averaging an average which will give a different answer.

ChlAve <- readRDS("/Users/jason/GitHub/GlobalZoopAbundModels/EnvironmentalData/five/Chl_raster_00Mean_fiveDeg.rds")
ChlAve_poly <- rasterToPolygons(ChlAve, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) # Convert to polygon which is better for plotting
ChlAve_sf <- st_as_sf(ChlAve_poly, crs = 4326) # Convert to sf
ChlAve_sf_moll <- st_transform(ChlAve_sf, crs = 54009) # Alter the CRS to mollweide


ChlNASA <- raster("/Users/jason/GitHub/GlobalZoopAbundModels/EnvironmentalData/nc/A20021852019304.L3m_CU_CHL_chlor_a_4km.nc",
                            varname = "chlor_a", nrows=180, ncols=360, xmn=-180, xmx=180, ymn=-90, ymx=90)
ChlNASA5Deg <- resample(ChlNASA, grd_lonlat, method='bilinear')
names(ChlNASA5Deg) <- "Chl"
## New environmental data sets
# new_raster <- readRDS('/Users/jason/GitHub/GlobalZoopAbundModels/EnvironmentalData/five/Chl_raster_00Mean_fiveDeg.rds')
ChlNASA5Deg_poly <- rasterToPolygons(ChlNASA5Deg, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) # Convert to polygon which is better for plotting
ChlNASA5Deg_sf <- st_as_sf(ChlNASA5Deg_poly, crs = 4326) # Convert to sf
ChlNASA5Deg_sf_moll <- st_transform(ChlNASA5Deg_sf, crs = 54009) # Alter the CRS to mollweide

## Difference
# Get the difference between the rasters
diff_raster <- ((ChlAve - old_raster)/ChlNASA5Deg)*100
names(diff_raster) <- "Chl_diff"
diff_poly <- rasterToPolygons(diff_raster, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) # Convert to polygon which is better for plotting
diff_sf <- st_as_sf(diff_poly) # Convert to sf
diff_sf_moll <- st_transform(diff_sf, crs = 54009) # Alter the CRS to mollweide



gg_jase <- ggplot() +
  geom_sf(data = ChlAve_sf_moll, aes(fill = log10(Chl)), color = NA) +
  geom_sf(data = world_sf, size = 0.05, fill = "grey20") +
  scale_fill_gradient(name = element_blank(),
                      limits = c(-1, 0.5),
                      low = "white",
                      high = "green",
                      na.value = "grey50",
                      guide = "colourbar",
                      aesthetics = "fill",
                      oob = scales::squish) +
  theme(plot.background = element_rect(fill = NA),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()) +
  ggtitle(label = "",subtitle = "Jason Average")

gg_old <- ggplot() +
  geom_sf(data = old_sf_moll, aes(fill = log10(Chl)), color = NA) +
  geom_sf(data = world_sf, size = 0.05, fill = "grey20") +
  scale_fill_gradient(name = element_blank(),
                      limits = c(-1, 0.5),
                      low = "white",
                      high = "green",
                      na.value = "grey50",
                      guide = "colourbar",
                      aesthetics = "fill",
                      oob = scales::squish) +
  theme(plot.background = element_rect(fill = NA),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()) +
  ggtitle(label = "",subtitle = "Ryans Data")

gg_diff <- ggplot() +
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
  ggtitle(label = "",subtitle = "Diff")



x11(height = 10, width = 7)
(gg_jase / gg_old) / (gg_diff) + plot_annotation(tag_levels = "A", tag_suffix = ")")
ggsave(filename="Figures/CompareSatelliteData_CompareAveraging.png", dpi = 400)










## Lets compare the NASA CU with my arithmatic average of the monthly data

## I don't expect them to match because I assume NASA will use the daily raw data to create the climatology,
## and so I am averaging an average which will give a different answer.

ChlAve <- readRDS("/Users/jason/GitHub/GlobalZoopAbundModels/EnvironmentalData/five/Chl_raster_00Mean_fiveDeg.rds")
ChlAve_poly <- rasterToPolygons(ChlAve, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) # Convert to polygon which is better for plotting
ChlAve_sf <- st_as_sf(ChlAve_poly, crs = 4326) # Convert to sf
ChlAve_sf_moll <- st_transform(ChlAve_sf, crs = 54009) # Alter the CRS to mollweide


## Difference
# Get the difference between the rasters
diff_raster <- ((ChlAve - Chl5Deg)/Chl5Deg)*100
names(diff_raster) <- "Chl_diff"
diff_poly <- rasterToPolygons(diff_raster, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) # Convert to polygon which is better for plotting
diff_sf <- st_as_sf(diff_poly) # Convert to sf
diff_sf_moll <- st_transform(diff_sf, crs = 54009) # Alter the CRS to mollweide



gg_jase <- ggplot() +
  geom_sf(data = ChlAve_sf_moll, aes(fill = log10(Chl)), color = NA) +
  geom_sf(data = world_sf, size = 0.05, fill = "grey20") +
  scale_fill_gradient(name = element_blank(),
                      limits = c(-1, 0.5),
                      low = "white",
                      high = "green",
                      na.value = "grey50",
                      guide = "colourbar",
                      aesthetics = "fill",
                      oob = scales::squish) +
  theme(plot.background = element_rect(fill = NA),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()) +
  ggtitle(label = "",subtitle = "Jason Average")

gg_nasa <- ggplot() +
  geom_sf(data = Chl5Deg_sf_moll, aes(fill = log10(Chl)), color = NA) +
  geom_sf(data = world_sf, size = 0.05, fill = "grey20") +
  scale_fill_gradient(name = element_blank(),
                      limits = c(-1, 0.5),
                      low = "white",
                      high = "green",
                      na.value = "grey50",
                      guide = "colourbar",
                      aesthetics = "fill",
                      oob = scales::squish) +
  theme(plot.background = element_rect(fill = NA),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        panel.border = element_blank()) +
  ggtitle(label = "",subtitle = "Old Satellite")

gg_diff <- ggplot() +
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
  ggtitle(label = "",subtitle = "Diff")



x11(height = 10, width = 7)
(gg_jase / gg_nasa) / (gg_diff) + plot_annotation(tag_levels = "A", tag_suffix = ")")
ggsave(filename="Figures/CompareSatelliteData_Averaging.png", dpi = 400)