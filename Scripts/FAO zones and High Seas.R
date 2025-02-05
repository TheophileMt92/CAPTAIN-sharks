library(sf)
library(raster)
library(dplyr)
library(tidyr)
library(here)
library(ggplot2)

FAO=st_read(here::here("Data", "fao", "World_Fao_Zones.shp"))
plot(st_geometry(FAO))

High_seas=st_read(here::here("Data", "World_High_Seas_v2_20241010", "High_Seas_v2.shp"))
plot(st_geometry(High_seas))

EEZ_world=st_read(here::here("Data", "World_EEZ_v12_20231025", "eez_boundaries_v12.shp"))
plot(st_geometry(EEZ_world))

# Ensure both layers are in the same CRS
High_seas <- st_transform(High_seas, st_crs(FAO))

# Create a template raster
# Adjust the resolution as needed. Lower numbers mean higher resolution but more processing time
template_raster <- raster(extent(st_bbox(FAO)), resolution = 0.1)

# Rasterize FAO zones
fao_raster <- rasterize(FAO, template_raster, field = "zone", background = NA)

# Rasterize High Seas (1 for High Seas, 0 for others)
highseas_raster <- rasterize(High_seas, template_raster, field = 1, background = 0)

# Stack the rasters
stacked_raster <- stack(fao_raster, highseas_raster)
names(stacked_raster) <- c("fao_zone", "high_seas")

# Create categories, merging all high seas zones
raster_df <- raster_df %>%
  mutate(category = case_when(
    high_seas == 1 ~ "High Seas",
    is.na(fao_zone) & high_seas == 0 ~ "Outside FAO zones",
    !is.na(fao_zone) & high_seas == 0 ~ paste("FAO zone", fao_zone)
  ))

# Convert category to factor and assign numeric values
raster_df$category_num <- as.numeric(factor(raster_df$category))

# Create a new raster from the dataframe
result_raster <- rasterFromXYZ(raster_df[, c("x", "y", "category_num")])

# Plot the raster
plot_data <- as.data.frame(result_raster, xy = TRUE)
colnames(plot_data)[3] <- "value"

ggplot() +
  geom_raster(data = plot_data, aes(x = x, y = y, fill = factor(value))) +
  scale_fill_viridis_d(name = "Category", 
                       labels = levels(factor(raster_df$category))) +
  coord_equal() +
  theme_minimal() +
  labs(title = "FAO Zones and High Seas",
       x = "Longitude",
       y = "Latitude")

# Save the plot
ggsave(here("Data", "FAO_HighSeas_map.png"), width = 12, height = 8, dpi = 300)

# If you want to save the raster for further use
writeRaster(result_raster, here("Data", "FAO_HighSeas_categories.tif"), overwrite=TRUE)