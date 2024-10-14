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

# Define a list of possible biome names
biome_names <- c("Tundra", "Taiga", "Temperate Deciduous Forest", "Tropical Rainforest", 
                 "Grassland", "Savanna", "Desert", "Alpine", "Aquatic")

# Generate random biome values
random_biomes <- sample(biome_names, nrow(distr_combined), replace = TRUE)

# Add the new Biome column at position 3
library(tibble)

distr_combined <- distr_combined %>%
  add_column(Biome = random_biomes, .after = 2)

# Print the dimensions of the updated tibble
cat("\nDimensions of the updated tibble:\n")
print(dim(distr_combined))

# Print the first few column names
cat("\nFirst few column names of the updated tibble:\n")
print(head(names(distr_combined), 10))

# Print the head of the updated tibble
cat("\nHead of the updated tibble:\n")
print(head(distr_combined))

grids=distr_combined

save(grids, file="grids_biomes_no.syn_corrected_04102024.Rdata")

summary(grids)

load("Data/grids_biomes_no.syn.Rdata")

str(distr_combined)
str(grids)
