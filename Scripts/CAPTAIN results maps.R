# Load required libraries
library(readr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(gridExtra)
library(biscale)
library(colorspace)
library(grid)

# EDGE ----

# Buget 0.1
# Read the RDS file from the Data folder
data <- readRDS("Data/EDGE_full_results_averaged_budget0.1_replicates100.rds")

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define the McBryde-Thomas 2 projection
mcbryde_thomas_2 <- "+proj=mbt_s"

# Transform the data to sf object and project
data_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = mcbryde_thomas_2)

# Project the world map
world_projected <- st_transform(world, crs = mcbryde_thomas_2)

# Create the outline from the projected world map
outline <- st_union(world_projected) %>% st_boundary()

# Create the globe bounding box
globe_bbox <- rbind(c(-180, -90), c(-180, 90), 
                    c(180, 90), c(180, -90), c(-180, -90))

# Create the globe border
globe_border <- st_polygon(list(globe_bbox)) %>%
  st_sfc(crs = 4326) %>%
  st_sf(data.frame(rgn = 'globe', geom = .)) %>%
  smoothr::densify(max_distance = 0.5) %>%
  st_transform(crs = mcbryde_thomas_2)

# Create the plot
ggplot() +
  # Add points for each planning unit, colored by priority
  geom_sf(data = data_sf, aes(color = Priority), size = 0.5, alpha = 0.7) +
  # Add the world map with a light gray fill
  geom_sf(data = world_projected, fill = "lightgrey", color = "lightgrey", size = 0.1) +
  # Add the black outline around the globe
  geom_sf(data = globe_border, fill = NA, color = "black", size = 0.5) +  # Add the globe border
  # Use a white to yellow to blue color gradient
  scale_color_gradientn(
    colors = c("white", "yellow", "darkblue"),
    values = c(0, 0.5, 1),
    name = "Priority",
    guide = guide_colorbar(barwidth = 20, barheight = 0.5, 
                           title.position = "top", title.hjust = 0.5)
  ) +
  # Add labels and title
  labs(title = "Global Distribution of Conservation Priorities",
       subtitle = "Index: EDGE2, Budget: 0.1, Replicates: 100",
       x = "Longitude", y = "Latitude") +
  # Adjust theme
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.margin = margin(t = 20, r = 0, b = 0, l = 0),
    legend.title = element_text(margin = margin(b = 10)),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank()
  )

# Save the plot
ggsave("priority_world_map_mcbryde_outline_EDGE2_0.1_100reps.png", width = 15, height = 10, dpi = 300)

# Buget 0.3
# Read the RDS file from the Data folder
data <- readRDS("Data/EDGE_full_results_averaged_budget0.3_replicates30.rds")

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define the McBryde-Thomas 2 projection
mcbryde_thomas_2 <- "+proj=mbt_s"

# Transform the data to sf object and project
data_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = mcbryde_thomas_2)

# Project the world map
world_projected <- st_transform(world, crs = mcbryde_thomas_2)

# Create the outline from the projected world map
outline <- st_union(world_projected) %>% st_boundary()

# Create the globe bounding box
globe_bbox <- rbind(c(-180, -90), c(-180, 90), 
                    c(180, 90), c(180, -90), c(-180, -90))

# Create the globe border
globe_border <- st_polygon(list(globe_bbox)) %>%
  st_sfc(crs = 4326) %>%
  st_sf(data.frame(rgn = 'globe', geom = .)) %>%
  smoothr::densify(max_distance = 0.5) %>%
  st_transform(crs = mcbryde_thomas_2)

# Create the plot
ggplot() +
  # Add points for each planning unit, colored by priority
  geom_sf(data = data_sf, aes(color = Priority), size = 0.5, alpha = 0.7) +
  # Add the world map with a light gray fill
  geom_sf(data = world_projected, fill = "lightgrey", color = "lightgrey", size = 0.1) +
  # Add the black outline around the globe
  geom_sf(data = globe_border, fill = NA, color = "black", size = 0.5) +  # Add the globe border
  # Use a white to yellow to blue color gradient
  scale_color_gradientn(
    colors = c("white", "yellow", "darkblue"),
    values = c(0, 0.5, 1),
    name = "Priority",
    guide = guide_colorbar(barwidth = 20, barheight = 0.5, 
                           title.position = "top", title.hjust = 0.5)
  ) +
  # Add labels and title
  labs(title = "Global Distribution of Conservation Priorities",
       subtitle = "Index: EDGE2, Budget: 0.3 Replicates: 30",
       x = "Longitude", y = "Latitude") +
  # Adjust theme
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.margin = margin(t = 20, r = 0, b = 0, l = 0),
    legend.title = element_text(margin = margin(b = 10)),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank()
  )

# Save the plot
ggsave("priority_world_map_mcbryde_outline_EDGE2_0.3_30reps.png", width = 15, height = 10, dpi = 300)

# FUSE ---- 
# Read the RDS file from the Data folder
data <- readRDS("Data/FUSE_full_results_averaged_budget0.1_replicates100.rds")

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define the McBryde-Thomas 2 projection
mcbryde_thomas_2 <- "+proj=mbt_s"

# Transform the data to sf object and project
data_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = mcbryde_thomas_2)

# Project the world map
world_projected <- st_transform(world, crs = mcbryde_thomas_2)

# Create the outline from the projected world map
outline <- st_union(world_projected) %>% st_boundary()

# Create the plot
ggplot() +
  # Add points for each planning unit, colored by priority
  geom_sf(data = data_sf, aes(color = Priority), size = 0.5, alpha = 0.7) +
  # Add the world map with a light gray fill
  geom_sf(data = world_projected, fill = "lightgrey", color = "lightgrey", size = 0.1) +
  # Add the black outline around the globe
  geom_sf(data = globe_border, fill = NA, color = "black", size = 0.5) +  # Add the globe border
  # Use a white to yellow to blue color gradient
  scale_color_gradientn(
    colors = c("white", "yellow", "darkblue"),
    values = c(0, 0.5, 1),
    name = "Priority",
    guide = guide_colorbar(barwidth = 20, barheight = 0.5, 
                           title.position = "top", title.hjust = 0.5)
  ) +
  # Add labels and title
  labs(title = "Global Distribution of Conservation Priorities",
       subtitle = "Index: FUSE, Budget: 0.1, Replicates: 100",
       x = NULL, y = NULL) +
  # Adjust theme
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.margin = margin(t = 20, r = 0, b = 0, l = 0),
    legend.title = element_text(margin = margin(b = 10)),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank()
  )

# Budget 0.3
data <- readRDS("Data/FUSE_full_results_averaged_budget0.3_replicates56.rds")

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define the McBryde-Thomas 2 projection
mcbryde_thomas_2 <- "+proj=mbt_s"

# Transform the data to sf object and project
data_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = mcbryde_thomas_2)

# Project the world map
world_projected <- st_transform(world, crs = mcbryde_thomas_2)

# Create the outline from the projected world map
outline <- st_union(world_projected) %>% st_boundary()

# Create the plot
ggplot() +
  # Add points for each planning unit, colored by priority
  geom_sf(data = data_sf, aes(color = Priority), size = 0.5, alpha = 0.7) +
  # Add the world map with a light gray fill
  geom_sf(data = world_projected, fill = "lightgrey", color = "lightgrey", size = 0.1) +
  # Add the black outline around the globe
  geom_sf(data = globe_border, fill = NA, color = "black", size = 0.5) +  # Add the globe border
  # Use a white to yellow to blue color gradient
  scale_color_gradientn(
    colors = c("white", "yellow", "darkblue"),
    values = c(0, 0.5, 1),
    name = "Priority",
    guide = guide_colorbar(barwidth = 20, barheight = 0.5, 
                           title.position = "top", title.hjust = 0.5)
  ) +
  # Add labels and title
  labs(title = "Global Distribution of Conservation Priorities",
       subtitle = "Index: FUSE, Budget: 0.3, Replicates: 56",
       x = NULL, y = NULL) +
  # Adjust theme
  theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    legend.box = "vertical",
    legend.margin = margin(t = 20, r = 0, b = 0, l = 0),
    legend.title = element_text(margin = margin(b = 10)),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    panel.grid = element_blank()
  )

# RGB map ----
# Budget 0.1
# Read the RDS files from the Data folder
edge_data <- readRDS("Data/EDGE_full_results_averaged_budget0.1_replicates100.rds")
fuse_data <- readRDS("Data/FUSE_full_results_averaged_budget0.1_replicates100.rds")

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define the McBryde-Thomas 2 projection
mcbryde_thomas_2 <- "+proj=mbt_s"

# Combine EDGE and FUSE data
combined_data <- edge_data %>%
  rename(EDGE_Priority = Priority) %>%
  left_join(fuse_data %>% rename(FUSE_Priority = Priority),
            by = c("Longitude", "Latitude"))

# Normalize priorities to 0-1 range
combined_data <- combined_data %>%
  mutate(
    EDGE_Priority_Norm = (EDGE_Priority - min(EDGE_Priority)) / (max(EDGE_Priority) - min(EDGE_Priority)),
    FUSE_Priority_Norm = (FUSE_Priority - min(FUSE_Priority)) / (max(FUSE_Priority) - min(FUSE_Priority))
  )

# Transform the data to sf object and project
data_sf <- st_as_sf(combined_data, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = mcbryde_thomas_2)

# Project the world map
world_projected <- st_transform(world, crs = mcbryde_thomas_2)

# Create the color palette
map_pal_raw <- bi_pal(pal = 'PurpleOr', dim = 4, preview = FALSE)

map_pal_mtx <- matrix(map_pal_raw, nrow = 4, ncol = 4)
map_pal_mtx[3, ] <- colorspace::lighten(map_pal_mtx[3, ], .1)
map_pal_mtx[2, ] <- colorspace::lighten(map_pal_mtx[2, ], .2)
map_pal_mtx[1, ] <- colorspace::lighten(map_pal_mtx[1, ], .3)
map_pal_mtx[ , 3] <- colorspace::lighten(map_pal_mtx[ , 3], .1)
map_pal_mtx[ , 2] <- colorspace::lighten(map_pal_mtx[ , 2], .2)
map_pal_mtx[ , 1] <- colorspace::lighten(map_pal_mtx[ , 1], .3)
map_pal_mtx[1, 1] <- '#ffffee'

map_pal <- as.vector(map_pal_mtx) %>% setNames(names(map_pal_raw))

# Create a custom color mapping function
get_color <- function(edge, fuse) {
  edge_class <- cut(edge, breaks = c(-Inf, 0.25, 0.5, 0.75, Inf), labels = 1:4)
  fuse_class <- cut(fuse, breaks = c(-Inf, 0.25, 0.5, 0.75, Inf), labels = 1:4)
  return(map_pal[(as.numeric(fuse_class)-1)*4 + as.numeric(edge_class)])
}

# Apply the function to get colors
data_sf$new_color <- mapply(get_color, data_sf$EDGE_Priority_Norm, data_sf$FUSE_Priority_Norm)
unique(data_sf$new_color)

# Create the main plot
main_plot <- ggplot() +
  geom_sf(data = data_sf, aes(color = new_color), size = 0.1, alpha = 1) +
  geom_sf(data = world_projected, fill = "lightgray", color = "lightgray") +
  geom_sf(data = globe_border, fill = NA, color = "grey70", size = 0.5) +  # Add the globe border
  scale_color_identity() +
  coord_sf() +
  theme_minimal() +
  labs(title = "Bivariate Map of EDGE and FUSE Priorities",
       x = "Longitude", y = "Latitude") +
  theme(plot.title = element_text(hjust = 0.5))

# Create the legend
legend_plot <- bi_legend(pal = map_pal, dim = 4,
                         xlab = 'EDGE',
                         ylab = 'FUSE')

# Define the layout matrix
layout <- rbind(c(1, 1),
                c(2, 3))

# Combine the main plot and legend using grid.arrange
combined_plot <- grid.arrange(
  main_plot, 
  legend_plot,
  rectGrob(gp = gpar(col = "white")),  # Blank space
  layout_matrix = layout,
  widths = c(0.15, 0.85),  # Legend takes 15% width, main plot and blank space take 85%
  heights = c(0.8, 0.2)  # Main plot takes 80% of height, legend row takes 20%
)

# Display the combined plot
print(combined_plot)

# Save the combined plot
ggsave("bivariate_priority_map_with_legend_100reps_budget01.png", combined_plot, width = 15, height = 12, dpi = 300)

# Budget 0.1
# Read the RDS files from the Data folder
edge_data <- readRDS("Data/EDGE_full_results_averaged_budget0.3_replicates30.rds")
fuse_data <- readRDS("Data/FUSE_full_results_averaged_budget0.3_replicates56.rds")

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Define the McBryde-Thomas 2 projection
mcbryde_thomas_2 <- "+proj=mbt_s"

# Combine EDGE and FUSE data
combined_data <- edge_data %>%
  rename(EDGE_Priority = Priority) %>%
  left_join(fuse_data %>% rename(FUSE_Priority = Priority),
            by = c("Longitude", "Latitude"))

# Normalize priorities to 0-1 range
combined_data <- combined_data %>%
  mutate(
    EDGE_Priority_Norm = (EDGE_Priority - min(EDGE_Priority)) / (max(EDGE_Priority) - min(EDGE_Priority)),
    FUSE_Priority_Norm = (FUSE_Priority - min(FUSE_Priority)) / (max(FUSE_Priority) - min(FUSE_Priority))
  )

# Transform the data to sf object and project
data_sf <- st_as_sf(combined_data, coords = c("Longitude", "Latitude"), crs = 4326) %>%
  st_transform(crs = mcbryde_thomas_2)

# Project the world map
world_projected <- st_transform(world, crs = mcbryde_thomas_2)

# Create the color palette
map_pal_raw <- bi_pal(pal = 'PurpleOr', dim = 4, preview = FALSE)

map_pal_mtx <- matrix(map_pal_raw, nrow = 4, ncol = 4)
map_pal_mtx[3, ] <- colorspace::lighten(map_pal_mtx[3, ], .1)
map_pal_mtx[2, ] <- colorspace::lighten(map_pal_mtx[2, ], .2)
map_pal_mtx[1, ] <- colorspace::lighten(map_pal_mtx[1, ], .3)
map_pal_mtx[ , 3] <- colorspace::lighten(map_pal_mtx[ , 3], .1)
map_pal_mtx[ , 2] <- colorspace::lighten(map_pal_mtx[ , 2], .2)
map_pal_mtx[ , 1] <- colorspace::lighten(map_pal_mtx[ , 1], .3)
map_pal_mtx[1, 1] <- '#ffffee'

map_pal <- as.vector(map_pal_mtx) %>% setNames(names(map_pal_raw))

# Create a custom color mapping function
get_color <- function(edge, fuse) {
  edge_class <- cut(edge, breaks = c(-Inf, 0.25, 0.5, 0.75, Inf), labels = 1:4)
  fuse_class <- cut(fuse, breaks = c(-Inf, 0.25, 0.5, 0.75, Inf), labels = 1:4)
  return(map_pal[(as.numeric(fuse_class)-1)*4 + as.numeric(edge_class)])
}

# Apply the function to get colors
data_sf$new_color <- mapply(get_color, data_sf$EDGE_Priority_Norm, data_sf$FUSE_Priority_Norm)
unique(data_sf$new_color)

# Create the main plot
main_plot <- ggplot() +
  geom_sf(data = data_sf, aes(color = new_color), size = 0.1, alpha = 1) +
  geom_sf(data = world_projected, fill = "lightgray", color = "lightgray") +
  geom_sf(data = globe_border, fill = NA, color = "grey70", size = 0.5) +  # Add the globe border
  scale_color_identity() +
  coord_sf() +
  theme_minimal() +
  labs(title = "Bivariate Map of EDGE and FUSE Priorities",
       x = "Longitude", y = "Latitude") +
  theme(plot.title = element_text(hjust = 0.5))

# Create the legend
legend_plot <- bi_legend(pal = map_pal, dim = 4,
                         xlab = 'EDGE',
                         ylab = 'FUSE')

# Define the layout matrix
layout <- rbind(c(1, 1),
                c(2, 3))

# Combine the main plot and legend using grid.arrange
combined_plot <- grid.arrange(
  main_plot, 
  legend_plot,
  rectGrob(gp = gpar(col = "white")),  # Blank space
  layout_matrix = layout,
  widths = c(0.15, 0.85),  # Legend takes 15% width, main plot and blank space take 85%
  heights = c(0.8, 0.2)  # Main plot takes 80% of height, legend row takes 20%
)

# Display the combined plot
print(combined_plot)

# Save the combined plot
ggsave("bivariate_priority_map_with_legend_50reps_budget03.png", combined_plot, width = 15, height = 12, dpi = 300)

# Create a custom palette
#custom_pal <- c("#ffffee", "#F4CEB4", "#FAAC7B", "#FE7C4B", 
#                "#D0C5E1", "#D0AAAB", "#CF896F", "#CA5F3A", 
#                "#AF9CD1", "#A7819A", "#9D625D", "#913D23", 
#                "#866EB6", "#775380", "#663746", "#551601")



# Protect fraction summary ----

#EDGE 0.1 
prot_frac=readRDS("Data/protect_fraction_summary_EDGE_0.1.rds")

sp <- fromJSON(here("Data", "shark_conservation_metrics_no_freshwater.json"))

Species = sp$EDGE2$info$Species
EDGE2 = sp$EDGE2$info$EDGE2

prot_frac$"Species"=Species
prot_frac$"EDGE2"=EDGE2

hist(prot_frac$Mean_Protect_Fraction)
hist(as.numeric(prot_frac$EDGE2))

str(prot_frac)

plot(prot_frac$EDGE2,prot_frac$Mean_Protect_Fraction)
