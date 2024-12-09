# Required packages used in this project:
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

# Where the data is being sourced from and the original files before augmentation:
setwd("D:\\Biogeography\\Data Files")
ecoregions <- st_read("D:\\Biogeography\\Data Files\\Ecoregions")
tree <- readRDS("D:\\Biogeography\\Data Files\\FIA_tree_master1.RDS")
iv_data <- readRDS("D:\\Biogeography\\Data Files\\df_master.rds")
# iv_data will need to be converted to a shape file for later subset of study area
iv_sf <- st_as_sf(iv_data, coords = c("LON","LAT"), crs = 4326)

# Manipulating the ecoregion data to display a target region (Appalachia):
appalachia <- ecoregions %>%
  filter(NA_L2NAME == "OZARK/OUACHITA-APPALACHIAN FORESTS")

appalachia_wgs84 <- st_transform(appalachia, crs = 4326)
appalachia_cords <- st_coordinates(appalachia_wgs84)
app_cords_df <- as.data.frame(appalachia_cords)
colnames(app_cords_df) <- c("X", "Y", "Shape 1", "Shape 2")
fin_appalachia <- app_cords_df %>%
  filter(X >= -89.5 & X <= -73) %>%
  select(-`Shape 2`)
colnames(fin_appalachia) <- c("LON","LAT","Shape")
fin_appalachia <- fin_appalachia %>%
  select(LAT, LON, Shape)

# Now the cleaned area of interest needs to be used to subset the FIA database:
fin_appalachia_sf <- st_as_sf(fin_appalachia, coords = c("LON","LAT"), crs = 4326)

bbox <- st_bbox(fin_appalachia_sf)
subset_data <- st_crop(iv_sf, bbox)

# The subset data only has polygons, it doesn't have coordinates so we will add them back:
temp_coords <- st_coordinates(subset_data)
coords_df <- as.data.frame(temp_coords)

subset_data <- subset_data %>%
  mutate(
    LAT = coords_df$Y,
    LON = coords_df$X )

# Calculating species richness within the subset data frame:
species_richness <- subset_data %>%
  group_by(GRIDID, LAT, LON) %>%
  summarise(richness = n_distinct(SPCD))

richness_classes <- classIntervals(species_richness$richness, n = 5, style = "jenks")
species_rich_class <- species_richness %>%
  mutate(JenksClass = cut(richness, breaks = richness_classes$brks, labels = 1:5, include.lowest = TRUE))

# Plot of current data:
ggplot() +
  geom_point(data = species_rich_class, aes(x = LON, y = LAT, color = JenksClass), size = 2) +
  scale_color_manual(values = c("#ffffb2","#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026"), name = "Species Richness")+
  geom_point(data = fin_appalachia, aes(x = LON, y = LAT), size = 0.5, color = "blue") +
  geom_point(data = red_maple_center, aes(x = LON, y = LAT), size = 8, color = 'hotpink')
