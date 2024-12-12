# Introduction:
The goal of this project is to investigate the well documented Rapoport Effect, but on a smaller scale of observation than is typical. I am very curious if within a specific ecoregion of the US the rapoport effect can be observed. The ecoregion I choose for this project is the Appalachia Province, which is well know for its abundant biological diversity.
There is an abundance of tree species found within the Appalachia Province, and the fact that it is geographically centered in the Eastern United States leads me to believe many tree species midpoints will fall within this area. For those factors I wanted to investigate how present the Rapoport Effect may be in this biological hotspot.

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
All of the tree characterstic data such as species richness and latitudinal range are sourced from the FIA data base.
The shape & extent of my target ecoregion was sourced from the Environmental Protection Agency. (https://www.epa.gov/eco-research/ecoregions-north-america)
The ecoregion data is level 2, which is detailed enough to isolate Appalachia but not so detailed to break it into subregions.
```{r}
setwd("D:\\Biogeography\\Data Files")
ecoregions <- st_read("D:\\Biogeography\\Data Files\\Ecoregions")
tree <- readRDS("D:\\Biogeography\\Data Files\\FIA_tree_master1.RDS")
iv_data <- readRDS("D:\\Biogeography\\Data Files\\df_master.rds")
### iv_data will need to be converted to a shape file for later subset of study area
iv_sf <- st_as_sf(iv_data, coords = c("LON","LAT"), crs = 4326)
```
### Manipulating the ecoregion data to display a target region (Appalachia):
This step proved to be critical to the success of the project. The ecoregion data are shapefiles with the regions themselves expressed as polygons.
The FIA database is a projected coordinate system in Latitude & Longitude. Significant resources and time were devoted to projected the ecoregion data, extracting Lat & Lon,
and finally making a structured data frame.
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
```
|Species | LAT_midpoint | LON_midpoint | Rapoport_Index |
|--------|--------------|--------------|----------------|
| boxelder | 38.80832 | -82.46690 | 1.3137777 |
| baldcypress | 32.81618 | -84.85485 | 1.0658709 |
| American Sycamore | 36.64056 | -84.26546 | 1.1652512 |
| American Elm | 38.16106 | -82.17432 | 1.3228237 |
| red maple | 34.24173 | -81.38487 | 1.3042623 |
```{r}
rap_all <- classIntervals(midpoints_df$Rapoport_Index, n = 5, style = "jenks")
rap_all_class <- midpoints_df %>%
  mutate(JenksClass = cut(Rapoport_Index, breaks = rap_all$brks, labels = 1:5, include.lowest = TRUE))

ggplot() +
  geom_point(data = iv_data, aes(x = LON, y = LAT), size = 1, color = "lightblue") +
  geom_point(data = rap_all_class, aes(x = LON_midpoint, y = LAT_midpoint, color = JenksClass), size = 2) +
  geom_point(data = fin_appalachia, aes(x = LON, y = LAT), size = 0.5, color = "red") +
  scale_color_manual(values = c("#ffffb2","#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026"), name = "Rapoport Index") +
  labs(title = "Eastern US Rapoport Index",
       x = "Longitude",
       y = "Latitude")
```
![Eastern US Big map](https://github.com/user-attachments/assets/bc812105-5e6f-495b-b001-a4e3064d79b5)

# Final Map
### Now that the Rapoport Index is calculated and I have isolated my study area to just the Appalachia Province, I can map them together:
```{r}
usa <- ne_states(country = "United States of America", returnclass = "sf")
usa31 <- usa %>%
  filter(longitude >= -95 & longitude <= -65)

usa31_sf <- st_as_sf(usa31, coords = c("longitude","latitude"), crs = 4326)
states_bbox <- st_crop(usa31_sf, bbox)

ggplot()+
  geom_sf(data = states_bbox, fill = "lightblue", color = "black") +
  geom_point(data = fin_appalachia, aes(x = LON, y = LAT), size = 0.5, color = "#006d2c") +
  geom_point(data = rap_classes, aes(x = LON_midpoint, y = LAT_midpoint, color = JenksClass), size = 2) +
  scale_color_manual(values = c("#ffffb2","#fed976","#feb24c","#fd8d3c","#f03b20","#bd0026"), name = "Rapoport Index") +
  labs(title = "Appalachia Province & Rapoport Effect",
       x = "Longitude",
       y = "Latitude")
```
![Area_wRapoport](https://github.com/user-attachments/assets/be72f9da-a836-440d-bc9e-9d2a9cc525e7)

# Extra
### I've also added the graph of Latitudinal Midpoint compared to Rapoport Index:
``` {r}
ggplot(midpoints_df, aes(x = LAT_midpoint, y = Rapoport_Index)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "blue") +
  labs(title = "Rapoport Effect",
       x = "Midpoint Latitude",
       y = "Rapoport Index") +
  theme_minimal()
```
![Rapoport Effect](https://github.com/user-attachments/assets/86fd4de5-5129-4d32-b648-5c0d2c3b658a)
