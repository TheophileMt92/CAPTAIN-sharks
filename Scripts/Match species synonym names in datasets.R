# Load necessary libraries
library(dplyr)
library(tidyr)
library(sf)
library(raster)
library(ggplot2)
library(rnaturalearth)
library(here)

# Change species names with synonyms ----
# Load spatial data
distr <- as_tibble(readRDS("Data/Map_PA_biomes_0_5.RDS"))

# Load data frame with synonyms
load('Data/iucn.names_traits.names.RData') # distr.names

# 1. Identify rows where Species and new.names don't match
mismatches <- distr.names %>%
  filter(Species != new.names)

# 2. Print the mismatched rows
cat("Species that don't match new.names:\n")
print(mismatches)

# Now compare Species names in distr.names by species names in distr
# 1. Extract species names from distr (starting from the 6th column onwards)
species_in_distr <- colnames(distr)[6:ncol(distr)]

# 2. Extract species names from distr.names
species_in_distr_names <- distr.names$Species

# 3. Find species in distr that are not in distr.names$Species
species_not_in_distr_names <- setdiff(species_in_distr, species_in_distr_names)

# 4. Print species in distr that are not in distr.names$Species
cat("Species in distr that are not in distr.names$Species:\n")
print(species_not_in_distr_names)

# Now exclude these species from distr
distr_cleaned <- distr %>%
  dplyr::select(-all_of(species_not_in_distr_names))

# Now, before renaming columns, let's check if the species names in distr_cleaned 
# match the Species column in distr.names

# 1. Extract species names from distr_cleaned (starting from the 6th column onwards)
species_in_distr_cleaned <- colnames(distr_cleaned)[6:ncol(distr_cleaned)]

# 2. Check which species in distr_cleaned are not in distr.names$Species
species_not_in_distr_names_cleaned <- setdiff(species_in_distr_cleaned, species_in_distr_names)

# 3. Print mismatches
cat("Species in distr_cleaned that are not in distr.names$Species:\n")
print(species_not_in_distr_names_cleaned)

# If there are no mismatches, proceed with renaming the columns using the new names
if (length(species_not_in_distr_names_cleaned) == 0) {
  
  # =========================================
  # Now proceed with renaming columns
  # =========================================
  
  # 1. Rename `X_ENTIER` and `Y_ENTIER` to `lon` and `lat`
  distr_cleaned <- distr_cleaned %>%
    rename(lon = X_ENTIER, lat = Y_ENTIER)
  
  # 2. Remove `Region`, `Biome`, and `Biome1` as these are not species
  distr_cleaned <- distr_cleaned %>%
    dplyr::select(-Region, -Biome, -Biome1)
  
  # 3. Now, rename species columns using `new.names`
  # Extract species columns (after removing the first 2, which are lon, lat)
  species_columns <- colnames(distr_cleaned)[3:ncol(distr_cleaned)]
  
  # Extract the corresponding new names (synonyms) from distr.names
  new_species_names <- distr.names$new.names
  
  # Check length consistency before renaming
  if (length(species_columns) != length(new_species_names)) {
    stop("The number of species columns in distr_cleaned does not match the number of new species names.")
  }

  # Rename the species columns
  colnames(distr_cleaned)[3:ncol(distr_cleaned)] <- new_species_names
  
  # 4. Print the cleaned and renamed tibble to check
  cat("Final column names in distr_cleaned:\n")
  print(colnames(distr_cleaned))
  
} else {
  cat("There are species in distr_cleaned that do not match the Species column in distr.names. Please review the mismatches before renaming.\n")
}

which(duplicated(names(distr_cleaned)))

library(dplyr)

# Function to combine duplicate columns and ensure presence/absence format
combine_duplicates <- function(df) {
  # Separate the first two columns
  first_two_cols <- df[, 1:2]
  
  # Process the remaining columns
  remaining_cols <- df[, -(1:2)]
  unique_names <- unique(names(remaining_cols))
  
  combined_cols <- lapply(unique_names, function(name) {
    cols <- which(names(remaining_cols) == name)
    if (length(cols) > 1) {
      # Combine duplicate columns and convert to presence/absence
      as.integer(rowSums(remaining_cols[, cols, drop = FALSE]) > 0)
    } else {
      # For single columns, ensure presence/absence format
      as.integer(remaining_cols[[name]] > 0)
    }
  })
  
  names(combined_cols) <- unique_names
  
  # Combine the first two columns with the processed columns
  bind_cols(first_two_cols, as_tibble(combined_cols))
}

# Apply the function to combine duplicate columns
distr_combined <- combine_duplicates(distr_cleaned)

# Check for any remaining duplicates
duplicate_cols <- which(duplicated(names(distr_combined)))

if(length(duplicate_cols) > 0) {
  cat("There are still duplicate column names after combining:\n")
  print(names(distr_combined)[duplicate_cols])
} else {
  cat("All duplicate columns have been successfully combined.\n")
}

# Print the dimensions of the new tibble
cat("\nDimensions of the combined tibble:\n")
print(dim(distr_combined))

# Print the first few column names
cat("\nFirst few column names of the combined tibble:\n")
print(head(names(distr_combined), 10))

# Print the head of the combined tibble
cat("\nHead of the combined tibble:\n")
print(head(distr_combined))

# Check the range of values in the combined tibble (excluding first two columns)
cat("\nRange of values in the combined tibble (excluding first two columns):\n")
print(range(as.matrix(distr_combined[, -(1:2)])))

# Plot species ranges ----
# Load necessary libraries
library(ggplot2)
library(rnaturalearth)
library(sf)
library(here)
library(dplyr)

# Assuming distr_combined is already loaded and contains the data

# Function to plot species distribution
plot_species_distribution <- function(data, species_name) {
  # Select the species data
  species_data <- data %>%
    dplyr::select(lon, lat, !!species_name) %>%
    rename(presence = !!species_name) %>%
    filter(presence != 0)
  
  # Get world map data
  world <- ne_countries(scale = "medium", returnclass = "sf")
  
  # Create the plot
  p <- ggplot() +
    geom_sf(data = world, fill = "lightgray", color = "darkgray") +
    geom_tile(data = species_data, aes(x = lon, y = lat), fill = "blue", alpha = 0.6) +
    coord_sf(crs = st_crs(4326)) +
    scale_x_continuous(limits = c(-180, 180), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-90, 90), expand = c(0, 0)) +
    labs(x = "Longitude", y = "Latitude", 
         title = paste("Distribution of", species_name),
         subtitle = "Based on distr_combined dataset") +
    theme_minimal()
  
  return(p)
}

# Get all species names (assuming first two columns are lon and lat)
species_names <- names(distr_combined)[3:ncol(distr_combined)]

# Plot the first species as an example
first_species <- species_names[1091]
p <- plot_species_distribution(distr_combined, first_species)
print(p)

# Now remove freshwater areas from the dataset ----
puvsp <- distr_combined %>%
  dplyr::mutate(across(3:ncol(.), ~ifelse(. == 0, NA, .))) %>% # Replace all 0 values with NA in every column from the 3rd to the last column
  dplyr::mutate(id = dplyr::row_number(), .before = 1) %>% # Adding a new column called "id" to the dataframe, filled with sequential numbers
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

# Filter puvsp to keep only marine points (depth < 0)
puvsp_marine <- puvsp %>%
  inner_join(bathy_0.5deg %>% filter(depth < 0), by = c("lon", "lat")) %>%
  dplyr::select(-last_col())  # This removes the last column (which should be 'depth')

# Now let's plot the extent of the remaining dataset
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")

ggplot() +
  #geom_sf(data = world, fill = "lightgray", color = "darkgray") +
  geom_tile(data = bathy_0.5deg %>% filter(depth < 0), aes(x = lon, y = lat), fill = "lightblue", alpha = 0.3) +
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
  dplyr::select(-id, -lon, -lat, -depth) %>%  # Remove non-species columns
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
species_to_plot <- species_with_only_NAs

# Filter the grids dataset for these species and create a long format dataframe
species_data <- distr_combined %>%
  dplyr::select(lon, lat, all_of(species_to_plot)) %>%
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
       title = "Freshwater Species Range Extents",
       subtitle = "Based on distr_combined dataset") +
  theme_minimal() +
  theme(legend.position = "none",
        panel.grid.major = element_line(color = "white", size = 0.2),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 6),
        axis.title = element_text(size = 8),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA))
p
