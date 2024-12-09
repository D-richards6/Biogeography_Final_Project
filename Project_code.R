library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)
library(rmapshaper)
library(rnaturalearth)
library(rnaturalearthdata)
library(classInt)
library(gridExtra)
library(eHOF)
library(devtools)
install_github("ropensci/rnaturalearthhires")
library(rnaturalearthhires)

setwd("D:\\Biogeography\\Data Files")
ecoregions <- st_read("D:\\Biogeography\\Data Files\\Ecoregions")
head(ecoregions)

unique(ecoregions$NA_L2NAME)

appalachia <- ecoregions %>%
  filter(NA_L2NAME == "OZARK/OUACHITA-APPALACHIAN FORESTS")
head(appalachia)
head(appalachia$geometry)

# Converting from polygon to coordinate system

# Creating Lat & Lon coordinates for the Appalachia only:
appalachia_wgs84 <- st_transform(appalachia, crs = 4326)
appalachia_cords <- st_coordinates(appalachia_wgs84)
head(appalachia_cords)
app_cords_df <- as.data.frame(appalachia_cords)
colnames(app_cords_df) <- c("X", "Y", "Shape 1", "Shape 2")
plot(app_cords_df)

# Only the area of Appalachia to be mapped, clean up study area map
fin_appalachia <- app_cords_df %>%
  filter(X >= -89.5 & X <= -73) %>%
  select(-`Shape 2`)
colnames(fin_appalachia) <- c("LON","LAT","Shape")
fin_appalachia <- fin_appalachia %>%
  select(LAT, LON, Shape)

# Adding the state shape files as background
usa <- ne_states(country = "United States of America", returnclass = "sf")
usa31 <- usa %>%
  filter(longitude >= -95 & longitude <= -65)

iv_data <- readRDS("D:\\Biogeography\\Data Files\\df_master.rds")
head(iv_data)

top10 <- iv_data %>%
  group_by(common_name) %>%
  summarise(Frequency = n()) %>%
  arrange(desc(Frequency)) %>%
  head(10) %>%
  pull(common_name)
print(top10)

# Isolating only the study area from the IV_data file. Keep all columns but only have it in Appalachia:
iv_sf <- st_as_sf(iv_data, coords = c("LON","LAT"), crs = 4326)
bbox <- st_bbox(fin_appalachia_sf)
subset_data <- st_crop(iv_sf, bbox)

temp_coords <- st_coordinates(subset_data)
coords_df <- as.data.frame(temp_coords)

names(subset_data)
subset_data <- subset_data %>%
  mutate(
    LAT = coords_df$Y,
    LON = coords_df$X
  )

# Finding species richness for the refined study area
species_richness2 <- subset_data %>%
  group_by(GRIDID, LAT, LON) %>%
  summarise(richness = n_distinct(SPCD))

richness_classes2 <- classIntervals(species_richness2$richness, n = 5, style = "jenks")
species_rich_class2 <- species_richness2 %>%
  mutate(JenksClass = cut(richness, breaks = richness_classes2$brks, labels = 1:5, include.lowest = TRUE))
ggplot()+
  geom_point(data = red_map, aes(x = LON, y = LAT, colour = IV), size = 1)

ggplot() +
  #geom_sf(data = usa31, fill = "lightblue", color = "black") +
  geom_point(data = species_rich_class2, aes(x = LON, y = LAT, color = JenksClass), size = 2) +
  scale_color_manual(values = c("#ffffb2","#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026"), name = "Species Richness")+
  geom_point(data = fin_appalachia, aes(x = LON, y = LAT), size = 0.5, color = "blue")


red_maple <- tree %>%
  filter(COMMON_NAME == "red maple")
latitude_values <- red_maple$LAT
longitude_values <- red_maple$LON

lat_analysis <- latitude_values
lat_range <- range(latitude_values, na.rm = TRUE)
lat_min <- min(latitude_values, na.rm = TRUE)
lat_max <- max(latitude_values, na.rm = TRUE)
lat_midpoint <- (lat_min + lat_max) / 2 
latitudinal_range <- lat_range[2] - lat_range[1]

lon_range <- range(longitude_values, na.rm = TRUE)
lon_min <- min(longitude_values, na.rm = TRUE)
lon_max <- max(longitude_values, na.rm = TRUE)
lon_midpoint <- (lon_min + lon_max) / 2

red_maple_center <- c(LAT = lat_midpoint, LON = lon_midpoint)
red_maple_center <- as.data.frame(red_maple_center)
red_maple_center <- red_maple_center %>%
  mutate(
    LAT = lat_midpoint,
    LON = lon_midpoint
  )

find_midpoint <- function(species_name, fia_data) {
  target_species <- tree %>%
    filter(COMMON_NAME == species_name)
  
  if(nrow(target_species) == 0) {
    return(data.frame(Species = species_name, Latitude_Midpoint = NA, Longitude_Midpoint = NA))
  }
  
  target_lat_values <- target_species$LAT
  target_lon_values <- target_species$LON
  
  target_lat_values <- target_lat_values[!is.na(target_lat_values)]
  target_lon_values <- target_lon_values[!is.na(target_lon_values)]
  
  if(length(target_lat_values) == 0 | length(target_lon_values) == 0) {
    return(data.frame(Species = species_name, Latitude_Midpoint = NA, Longitude_Midpoint = NA))
  }
  
  lat_min <- min(target_lat_values, na.rm = TRUE)
  lat_max <- max(target_lat_values, na.rm = TRUE)
  lon_min <- min(target_lon_values, na.rm = TRUE)
  lon_max <- max(target_lon_values, na.rm = TRUE)
  
  target_lat_midpoint <- (lat_min + lat_max) / 2
  target_lon_midpoint <- (lon_min + lon_max) / 2
  
  midpoint_df <- data.frame(
    Species = species_name,
    LAT_midpoint = target_lat_midpoint,
    LON_midpoint = target_lon_midpoint
  )
  return(midpoint_df)
}

find_midpoint("boxelder", tree)
species_list <- c("red maple","black cherry","American elm","white oak","northern red oak","green ash","white ash","blackgum","sugar maple","black oak")
species_midpoint_df <- bind_rows(lapply(tree_list, function(species_name) find_midpoint(species_name, fia_data)))
find_midpoint(tree_list, tree)

species_midpoints_list <- lapply(tree_list, function(species) find_midpoint(species, fia_data))
midpoints_df <- bind_rows(species_midpoints_list)
print(midpoints_df)

tree_list <- unique(na.omit(tree$COMMON_NAME))
tree_list

ggplot() +
  geom_point(data = iv_data, aes(x = LON, y = LAT, color = IV), size = 1) +
  geom_point(data = midpoints_df, aes(x = LON_midpoint, y = LAT_midpoint), size = 3, color = "hotpink") +
  geom_point(data = fin_appalachia, aes(x = LON, y = LAT), size = 0.5, color = "red")





