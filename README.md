### Required packages:
```{r, warning=FALSE, message=FALSE}
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
library(rnaturalearthhires)
```
### Where the data is being sourced from and the original files before augmentation:
```{r}
setwd("D:\\Biogeography\\Data Files")
ecoregions <- st_read("D:\\Biogeography\\Data Files\\Ecoregions")
tree <- readRDS("D:\\Biogeography\\Data Files\\FIA_tree_master1.RDS")
iv_data <- readRDS("D:\\Biogeography\\Data Files\\df_master.rds")
### iv_data will need to be converted to a shape file for later subset of study area
iv_sf <- st_as_sf(iv_data, coords = c("LON","LAT"), crs = 4326)
```
### Manipulating the ecoregion data to display a target region (Appalachia):
```{r}
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
```
### Now the cleaned area of interest needs to be used to subset the FIA database:
```{r}
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
```
### Calculating species richness within the subset data frame:
```{r}
species_richness <- subset_data %>%
  group_by(GRIDID, LAT, LON) %>%
  summarise(richness = n_distinct(SPCD))

richness_classes <- classIntervals(species_richness$richness, n = 5, style = "jenks")
species_rich_class <- species_richness %>%
  mutate(JenksClass = cut(richness, breaks = richness_classes$brks, labels = 1:5, include.lowest = TRUE))
```
### Species Richness Plot of the Appalachia Province with normalized classes:
```{r}
ggplot() +
  geom_point(data = species_rich_class, aes(x = LON, y = LAT, color = JenksClass), size = 2) +
  scale_color_manual(values = c("#ffffb2","#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026"), name = "Species Richness")+
  geom_point(data = fin_appalachia, aes(x = LON, y = LAT), size = 0.5, color = "blue") +
  labs(title = "Species Richness of the Appalachia Province",
       x = "Longitude",
       y = "Latitude" )
```
![Appalachia Species Richness](https://github.com/user-attachments/assets/65339738-2ded-426b-8b06-548c6a5b9bcd)

# Rapoport Effect
### To calculate the Rapoport Index for a given tree species you need both the latitudinal range & latitudinal midpoint of that species habitat.
```{r}
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
```
### This was how I first calculated the Rapoport Index, but this code is messy with parts that weren't even used in the end. Let's clean it up and also make it a function:
```{r}
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
  lat_range <- range(target_lat_values, na.rm = TRUE)
  lon_min <- min(target_lon_values, na.rm = TRUE)
  lon_max <- max(target_lon_values, na.rm = TRUE)
  lon_range <- range(target_lon_values, na.rm = TRUE)
  
  target_lat_midpoint <- (lat_min + lat_max) / 2
  target_lon_midpoint <- (lon_min + lon_max) / 2
  
  midpoint_df <- data.frame(
    Species = species_name,
    LAT_midpoint = target_lat_midpoint,
    LON_midpoint = target_lon_midpoint,
    Rapoport_Index = lat_range / lat_midpoint
  )
  return(midpoint_df)
}
```
### Now that it is a function, might as well find the rapoport index for every tree species in the FIA dataset. I added all of the common names into a list and then used lapply to run each common name, output the value and then store all of them together in a master data frame.
```{r}
tree_list <- unique(na.omit(tree$COMMON_NAME))
species_midpoints_list <- lapply(tree_list, function(species) find_midpoint(species, fia_data))
midpoints_df <- bind_rows(species_midpoints_list)
print(midpoints_df)
