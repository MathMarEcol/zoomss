library(tidyverse)
library(patchwork)
library(raster)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)



df_raster <- rasterFromXYZ(df)  #Convert first two columns as lon-lat and third as value
df_poly <- rasterToPolygons(df_raster, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) # Convert to polyg



biomass_brick <- readRDS("../GlobalZoopAbundModels/ModelOutput/GlobalLayers/Num_Brick_Biomass.rds")
BB <- dropLayer(biomass_brick, c("Fish_Small", "Fish_Med", "Fish_Large", "SST", "Chl"))

BB_sum <- calc(BB, fun = sum) # Sum the data for the proportion
names(BB_sum) <- "ZoopSum"

BBProp <- BB/BB_sum # Calculate Proportion
names(BBProp) <- names(BB) # Copy layer names across

BBPropPoly <- rasterToPolygons(BBProp, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) # Convert to polygon which is better for plotting
BBProp_sf <- st_as_sf(BBPropPoly, crs = 4326) # Convert to sf

st_crs(BBProp_sf) <- 4326 # Set base projection
BBProp_sf_moll <- st_transform(BBProp_sf, crs = 54009) # Alter the CRS to mollweide
BBProp_sf_moll <- BBProp_sf_moll %>%
  dplyr::select(c("Larvaceans", "CarnCopepods",
                  "Salps", "Chaetognaths",
                  "OmniCopepods", "Jellyfish",
                  "Euphausiids", "geometry" ))

BB_Group_sum <- BB_sum
BB_Group_sum <- addLayer(BB_Group_sum, BB$Larvaceans + BB$Salps)
BB_Group_sum <- addLayer(BB_Group_sum, BB$CarnCopepods + BB$Jellyfish + BB$Chaetognaths)
BB_Group_sum <- addLayer(BB_Group_sum, BB$OmniCopepods + BB$Euphausiids)
names(BB_Group_sum) <- c("ZoopSum", "FilterFeeders", "Carnivores", "Omnivores")
BB_GroupProp <- BB_Group_sum/BB_Group_sum$ZoopSum # Calculate Proportion
names(BB_GroupProp) <- names(BB_Group_sum) # Copy layer names across

BB_GroupPropPoly <- rasterToPolygons(BB_GroupProp, fun=NULL, n=4, na.rm=TRUE, digits=12, dissolve=FALSE) # Convert to polygon which is better for plotting
BB_GroupProp_sf <- st_as_sf(BB_GroupPropPoly, crs = 4326) # Convert to sf

st_crs(BB_GroupProp_sf) <- 4326 # Set base projection
BB_GroupProp_sf_moll <- st_transform(BB_GroupProp_sf, crs = 54009) # Alter the CRS to mollweide

rm(biomass_brick, BB, BB_sum, BBProp, BBPropPoly, BBProp_sf)

world <- ne_countries(scale = "medium", returnclass = "sf")
world_sf <- st_transform(world, crs = 54009) # Convert to Mollweide

theme_opts <- list(theme(panel.grid.minor = element_blank(),
                         panel.grid.major = element_line(colour = "grey50", size = 0.2),
                         panel.background = element_blank(),
                         panel.border = element_blank(),
                         plot.background = element_rect(fill = NA),
                         plot.title = element_text(hjust = 0.5),
                         plot.margin = unit(c(0,0,0,0), "mm"),
                         axis.line = element_blank(),
                         axis.text.x = element_blank(),
                         axis.text.y = element_blank(),
                         axis.ticks = element_blank(),
                         axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         # legend.title = element_text(size = 6),
                         # legend.text = element_text(size = 6),
                         legend.position = "right",
                         # legend.direction = "horizontal",
                         # legend.background = element_rect(fill = NA),
                         legend.key.height = unit(9, "mm"),
                         legend.key.width = unit(4, "mm"),
                         # legend.position = c(0.5, -0.05),
))

colr <- list()
colr[[1]] <- c("#eff3ff", "#c6dbef", "#9ecae1", "#6baed6", "#3182bd", "#08519c") # Blues
colr[[2]] <- c("#fee5d9", "#fcbba1", "#fc9272", "#fb6a4a", "#de2d26", "#a50f15") # Reds
colr[[3]] <- c("#f5f5f5", "#f6e8c3", "#dfc27d", "#bf812d", "#8c510a", "#543005") # Browns
#c("#edf8fb", "#ccece6", "#99d8c9", "#66c2a4", "#2ca25f", "#006d2c") # Greens

myplots <- list()
for(j in names(BBProp_sf_moll)){
  if(j != "geometry"){
    myplots[[j]] <- ggplot() +
      geom_sf(data = BBProp_sf_moll, aes_string(fill = j), color = NA) +
      geom_sf(data = world_sf, size = 0.05, fill = "grey20") +
      scale_fill_gradient2(limits = c(as.numeric(quantile(as.vector(dplyr::select(as_tibble(BBProp_sf_moll), !! j)[[1]]), 0.05)),
                                      as.numeric(quantile(as.vector(dplyr::select(as_tibble(BBProp_sf_moll), !! j)[[1]]), 0.95))),
                           low = "white",
                           high = TestGroups$Colour[TestGroups$species==j],
                           space = "Lab",
                           na.value = "grey50",
                           aesthetics = "fill",
                           oob = scales::squish,
                           guide = guide_colourbar(title = "Proportion of\nMesozooplankton",
                                                   title.position = "right",
                                                   title.hjust = 0.5,
                                                   title.theme = element_text(angle = 270, size = 10))) +
      theme_opts +
      scale_alpha(range = c(-0, 0.5)) +
      (if(j=="OmniCopepods"){ggtitle("Omnivorous Copepods")}
       else{if(j=="CarnCopepods"){ggtitle("Carnivorous Copepods")}
         else{if(j=="Euphausiids"){ggtitle("Krill")}
           ggtitle(j)}})
  }
}

x11(height = 10, width = 10)
wrap_plots(myplots) + plot_layout(ncol = 2) + plot_annotation(tag_levels = "A", tag_suffix = ")")
ggsave("Figures/FoodWebs_Supp_BiomassProp.png", dpi = 500)

# Now plot by feeding group

BB_GroupProp_sf_moll <- BB_GroupProp_sf_moll %>%
  dplyr::select(Omnivores, Carnivores, FilterFeeders, everything()) # Reorder

jj <- 0

myplots <- list()
for(j in names(BB_GroupProp_sf_moll)){
  if(j != "geometry" & j != "ZoopSum"){
    jj <- jj + 1
    myplots[[j]] <- ggplot() +
      geom_sf(data = BB_GroupProp_sf_moll, aes_string(fill = j), color = NA) +
      geom_sf(data = world_sf, size = 0.05, fill = "grey20") +
      scale_fill_gradientn(name =  "Proportion",
                           colours = colr[[jj]], # BuGn from RColorBRewer
                           limits = c(as.numeric(quantile(as.vector(dplyr::select(as_tibble(BB_GroupProp_sf_moll), !! j)[[1]]), 0.05)),
                                      as.numeric(quantile(as.vector(dplyr::select(as_tibble(BB_GroupProp_sf_moll), !! j)[[1]]), 0.95))
                           ),
                           oob = scales::squish,
                           guide = guide_colourbar(title = "Proportion of\nMesozooplankton",
                                                   title.position = "right",
                                                   title.hjust = 0.5,
                                                   title.theme = element_text(angle = 270, size = 10))) +
      theme_opts +
      scale_alpha(range = c(-0, 0.5)) +
      (if(j=="FilterFeeders"){ggtitle("Filter Feeders")} else{ggtitle(j)})
  }
}

x11(height = 7, width = 5)
wrap_plots(myplots) + plot_layout(ncol = 1) + plot_annotation(tag_levels = "A", tag_suffix = ")")
ggsave("Figures/FoodWebs_2_BiomassProp_ByGroup.png", dpi = 500)


