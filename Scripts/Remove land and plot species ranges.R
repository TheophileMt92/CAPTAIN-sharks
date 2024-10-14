library(dplyr)
library(tidyr)
library(sf)
library(raster)
library(ggplot2)

# Remove land from the dataset and identifty which species this corresponds too ----
# Loading the data
load("Data/Raw/CAPTAIN_Occurance_Grid.RData")

puvsp <- grids %>%
  dplyr::mutate(across(4:ncol(.), ~ifelse(. == 0, NA, .))) %>% # Replace all 0 values with NA in every column from the 4th to the last column
  dplyr::mutate(id = dplyr::row_number(), .before = "biomes") %>% #Adding a new column called "pu" to the dataframe, filled with sequential numbers (1, 2, 3)
  dplyr::select(-c(biomes)) %>% #Remove biomes columns
  dplyr::filter(if_any(4:ncol(.), ~!is.na(.)))

# Load the bathymetry raster
Bathy <- raster(here::here("Data", "bathymetry-0.1deg-adjusted.tif"))

# Extract raster to dataframe at its original resolution
bathy_df <- as.data.frame(Bathy, xy = TRUE)
colnames(bathy_df) <- c("lon", "lat", "depth")

# Aggregate to 0.5 degree resolution
bathy_0.5deg <- bathy_df %>%
  mutate(
    lon = floor(lon * 2) / 2,
    lat = floor(lat * 2) / 2
  ) %>%
  group_by(lon, lat) %>%
  summarise(depth = mean(depth, na.rm = TRUE)) %>%
  ungroup()

# Ensure puvsp coordinates are rounded to 0.5 degrees
puvsp <- puvsp %>%
  mutate(
    lon = round(lon * 2) / 2,
    lat = round(lat * 2) / 2
  )

# Filter puvsp to keep only marine points (depth > 0)
puvsp_marine <- puvsp %>%
  inner_join(bathy_0.5deg %>% filter(depth < 0), by = c("lon", "lat"))

# Now let's plot the extent of the remaining dataset
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  geom_sf(data = world, fill = "lightgray", color = "darkgray") +
  geom_tile(data = bathy_0.5deg %>% filter(depth > 0), aes(x = lon, y = lat), fill = "lightblue", alpha = 0.3) +
  geom_point(data = puvsp_marine, aes(x = lon, y = lat), color = "blue", size = 0.5, alpha = 0.5) +
  coord_sf(crs = st_crs(4326)) +
  scale_x_continuous(limits = c(-180, 180), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-90, 90), expand = c(0, 0)) +
  labs(x = "Longitude", y = "Latitude", 
       title = "Global Extent of Marine Species Dataset",
       subtitle = paste("Total points:", nrow(puvsp_marine))) +
  theme_minimal() +
  theme(panel.grid.major = element_line(color = "white", size = 0.2),
        panel.grid.minor = element_blank())

# Identify species with only NAs
species_with_only_NAs <- puvsp_marine %>%
  select(-id, -lon, -lat, -depth) %>%  # Remove non-species columns
  summarise(across(everything(), ~all(is.na(.)))) %>%
  pivot_longer(everything(), names_to = "species", values_to = "all_NA") %>%
  filter(all_NA == TRUE) %>%
  pull(species)

# Print the species names
print(species_with_only_NAs)

# Count how many species have only NAs
n_species_only_NAs <- length(species_with_only_NAs)
print(paste("Number of species with only NAs:", n_species_only_NAs))

# Calculate the percentage of species with only NAs
total_species <- ncol(puvsp_marine) - 4  # Subtract 4 for id, lon, lat, and depth columns
percentage_species_only_NAs <- (n_species_only_NAs / total_species) * 100
print(paste("Percentage of species with only NAs:", round(percentage_species_only_NAs, 2), "%"))

# Plot the ranges of the removed species ----
library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)
library(rnaturalearth)

# List of species to plot
species_to_plot <- c("Bythaelurus clevai", "Eridacnis barbouri", "Eusphyra blochii", 
                     "Galeus sauteri", "Glaucostegus cemiculus", "Heterodontus portusjacksoni", 
                     "Myliobatis aquila", "Pseudoginglymostoma brevicaudatum", "Pseudoraja fischeri", 
                     "Pseudotriakis microdon", "Pteroplatytrygon violacea", "Raja maderensis", 
                     "Raja microocellata", "Rhinobatos horkelii", "Rhinobatos productus", 
                     "Sinobatis melanosoma")

# Filter the grids dataset for these species and create a long format dataframe
species_data <- grids %>%
  select(lon, lat, all_of(species_to_plot)) %>%
  pivot_longer(cols = -c(lon, lat), names_to = "species", values_to = "presence") %>%
  filter(presence != 0)  # Filter out cells where the species is not present

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Create the plot
p <- ggplot() +
  geom_sf(data = world, fill = "lightgray", color = "darkgray") +
  geom_tile(data = species_data, aes(x = lon, y = lat, fill = species), alpha = 0.6) +
  facet_wrap(~ species, ncol = 4) +
  coord_sf(crs = st_crs(4326)) +
  scale_x_continuous(limits = c(-180, 180), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-90, 90), expand = c(0, 0)) +
  labs(x = "Longitude", y = "Latitude", 
       title = "Species Range Extents",
       subtitle = "Based on grids dataset") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_line(color = "white", size = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# Save the plot
ggsave(here::here("species_range_extents_grids.pdf"), plot = p, width = 15, height = 15)

# Plot the ranges of 16 random species ----
#Randomly selected 16 species 
library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)
library(rnaturalearth)
library(here)

# Get the list of all species in the dataset
all_species <- colnames(grids)[!(colnames(grids) %in% c("lon", "lat", "biomes"))]

# Randomly select 16 species
set.seed(123)  # for reproducibility
species_to_plot <- sample(all_species, 16)

# Filter the grids dataset for these species and create a long format dataframe
species_data <- grids %>%
  select(lon, lat, all_of(species_to_plot)) %>%
  pivot_longer(cols = -c(lon, lat), names_to = "species", values_to = "presence") %>%
  filter(presence != 0)  # Filter out cells where the species is not present

# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")

# Create the plot
p <- ggplot() +
  geom_sf(data = world, fill = "lightgray", color = "darkgray") +
  geom_tile(data = species_data, aes(x = lon, y = lat, fill = species), alpha = 0.6) +
  facet_wrap(~ species, ncol = 4) +
  coord_sf(crs = st_crs(4326)) +
  scale_x_continuous(limits = c(-180, 180), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-90, 90), expand = c(0, 0)) +
  labs(x = "Longitude", y = "Latitude", 
       title = "Distribution of 16 Randomly Selected Species",
       subtitle = "Based on grids dataset") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_line(color = "white", size = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# Save the plot as PNG
ggsave(here::here("random_16_species_distribution_grids2.png"), plot = p, width = 15, height = 15, dpi = 300, bg = "white")

# Save the plot as PDF
ggsave(here("random_16_species_distribution.pdf"), plot = p, width = 15, height = 15)

# Plot the ranges of the "freshwater species" in grids_biomes_no.syn.Rdata ---- 
#
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(rnaturalearth)

#data <- readRDS("Data/Map_PA_biomes_0_5.RDS")

load(here::here("Data", "grids_biomes_no.syn.Rdata"))
head(grids)

# List of species to plot
species_to_plot <- c("Bythaelurus clevai", "Eridacnis barbouri", "Eusphyra blochii", 
                     "Galeus sauteri", "Glaucostegus cemiculus", "Heterodontus portusjacksoni", 
                     "Myliobatis aquila", "Pseudoginglymostoma brevicaudatum", "Pseudoraja fischeri", 
                     "Pseudotriakis microdon", "Pteroplatytrygon violacea", "Raja maderensis", 
                     "Raja microocellata", "Rhinobatos horkelii", "Rhinobatos productus", 
                     "Sinobatis melanosoma")

# Filter the grids dataset to include the species and create a long format dataframe
species_data <- grids %>%
  select(lon, lat, all_of(species_to_plot)) %>%
  pivot_longer(cols = -c(lon, lat), names_to = "species", values_to = "presence") %>%
  filter(presence != 0)  # Filter out cells where the species is not present

# Get world map data using rnaturalearth package
world <- ne_countries(scale = "medium", returnclass = "sf")

# Create the plot
p <- ggplot() +
  geom_sf(data = world, fill = "lightgray", color = "darkgray") +
  geom_tile(data = species_data, aes(x = lon, y = lat, fill = species), alpha = 0.6) +
  facet_wrap(~ species, ncol = 4) +
  coord_sf(crs = st_crs(4326)) +
  scale_x_continuous(limits = c(-180, 180), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-90, 90), expand = c(0, 0)) +
  labs(x = "Longitude", y = "Latitude", 
       title = "Species Range Extents",
       subtitle = "Based on grids dataset") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_line(color = "white", size = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))

# Show the plot
print(p)

# Plot and save species ranges comparing Map_PA_biomes_0_5.RDS and grids_biomes_no.syn.Rdata ----
# Load necessary libraries
library(ggplot2)
library(rnaturalearth)
library(sf)

# Load the dataset
data <- readRDS("Data/Map_PA_biomes_0_5.RDS")

# Select the first species to plot ("Glaucostegus cemiculus")
first_species <- "Myliobatis aquila"

# Manually select the lon, lat, and the first species columns
species_data <- data[, c("X_ENTIER", "Y_ENTIER", first_species)]

# Rename the columns for clarity
colnames(species_data) <- c("lon", "lat", "presence")

# Filter out rows where the species is absent (i.e., presence is 0)
species_data <- species_data[species_data$presence != 0, ]

# Get world map data using rnaturalearth package
world <- ne_countries(scale = "medium", returnclass = "sf")

# Create the plot for the first species
p <- ggplot() +
  geom_sf(data = world, fill = "lightgray", color = "darkgray") +
  geom_tile(data = species_data, aes(x = lon, y = lat), fill = "blue", alpha = 0.6) +
  coord_sf(crs = st_crs(4326)) +
  scale_x_continuous(limits = c(-180, 180), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-90, 90), expand = c(0, 0)) +
  labs(x = "Longitude", y = "Latitude", 
       title = paste("Distribution of", first_species),
       subtitle = "Based on grids dataset") +
  theme_minimal()

print(p)

# Save the plot as an image file
output_filename <- paste0("Distribution of ", first_species, " from Map_PA_biomes_0_5.RDS.pdf")
ggsave(output_filename, plot = p, device = "pdf", width = 10, height = 6)

# Load necessary libraries
library(ggplot2)
library(rnaturalearth)
library(sf)
library(here)

# Load the dataset
load(here::here("Data", "grids_biomes_no.syn.Rdata"))

# Assuming the dataset is named 'grids' after loading
# Select the species to plot ("Glaucostegus cemiculus")
first_species <- "Myliobatis aquila"

# Manually select the lon, lat, and the first species columns
species_data <- grids[, c("lon", "lat", first_species)]

# Rename the columns for clarity
colnames(species_data) <- c("lon", "lat", "presence")

# Filter out rows where the species is absent (i.e., presence is 0)
species_data <- species_data[species_data$presence != 0, ]

# Get world map data using rnaturalearth package
world <- ne_countries(scale = "medium", returnclass = "sf")

# Create the plot for the species
p <- ggplot() +
  geom_sf(data = world, fill = "lightgray", color = "darkgray") +
  geom_tile(data = species_data, aes(x = lon, y = lat), fill = "blue", alpha = 0.6) +
  coord_sf(crs = st_crs(4326)) +
  scale_x_continuous(limits = c(-180, 180), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-90, 90), expand = c(0, 0)) +
  labs(x = "Longitude", y = "Latitude", 
       title = paste("Distribution of", first_species),
       subtitle = "Based on grids dataset") +
  theme_minimal()

# Show the plot
print(p)

# Save the plot as a PDF file
output_filename <- paste0("Distribution of ", first_species, " from grids_biomes_no.syn.Rdata.pdf")
ggsave(output_filename, plot = p, device = "pdf", width = 10, height = 6)

head(data)
