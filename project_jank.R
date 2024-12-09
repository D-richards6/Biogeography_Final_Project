install.packages("raster")
library(raster)
setwd("D:\\Biogeography\\Data Files")
tree <- readRDS("D:\\Biogeography\\Data Files\\FIA_tree_master1.RDS")
iv_data <- readRDS("D:\\Biogeography\\Data Files\\df_master.rds")

fin_appalachia_sf <- st_as_sf(fin_appalachia, coords = c("LON","LAT"), crs = 4326)
plot(fin_appalachia_sf)

hull_index <- chull(fin_appalachia$LAT, fin_appalachia$LON)
ordered <- fin_appalachia[c(hull_index, hull_index[1]), ]

fin_appalachia_poly <- st_cast(fin_appalachia_sf, "POLYGON")

convex_hull <- st_convex_hull(fin_appalachia_sf)
convex_hull_polygon <- st_cast(convex_hull, "POLYGON")

plot(appalachia)
fia_sf <- st_as_sf(tree, coords = c("LON","LAT"), crs = 4326)
iv_sf <- st_as_sf(iv_data, coords = c("LON","LAT"), crs = 4326)

appalachia_sf <- st_as_sf(appalachia, )
intersection_test <- st_intersection(fia_sf, fin_appalachia_sf)

bbox <- st_bbox(fin_appalachia_sf)
subset_data <- st_crop(iv_sf, bbox)

plot(st_geometry(subset_data))

ggplot() +
  geom_point(data = subset_richness, aes(x = LON, y = LAT, color = richness), size = 1) +
  geom_point(data = fin_appalachia, aes(x = LON, y = LAT), size = 0.5, color = "red")

temp_coords <- st_coordinates(subset_data)
coords_df <- as.data.frame(temp_coords)

names(subset_data)
subset_data <- subset_data %>%
  mutate(
    LAT = coords_df$Y,
    LON = coords_df$X
  )
col_num <- 15:16
subset_data <- subset_data[, -col_num]

subset_richness <- subset_data%>%
  group_by(GRIDID, LAT, LON) %>%
  summarise(richness = n_distinct(SPCD))
  