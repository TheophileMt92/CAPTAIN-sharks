# Loading in the required libraries ----
library(tidyverse)
library(readxl)
library(data.table)
library(raster)
library(rasterVis)
library(sf)
library(vegan)
library(dplyr)
library(tidyr)
library(here)

# First three databases: pu_vsp, pu_fishing, pu_coords ----
# Loading the data
load(here::here("Data", "puvsp_marine.Rdata"))

# GRID data transformation

# Define the freshwater species to be removed
puvsp <- puvsp_marine %>%
  dplyr::mutate(across(4:ncol(.), ~ifelse(. == 0, NA, .))) %>% # Replace all 0 values with NA in every column from the 4th to the last column
  dplyr::mutate(id = dplyr::row_number()) %>% #Adding a new column called "id" to the dataframe, filled with sequential numbers (1, 2, 3)
  dplyr::select(-c(lon, lat)) %>% #Remove coordinate columns
  dplyr::filter(if_any(4:ncol(.), ~!is.na(.))) %>% # Remove rows where all species columns are NA
  tidyr::pivot_longer(cols = -id, 
                      names_to = "sp", 
                      values_to = "amount", 
                      values_drop_na = TRUE) %>% #Melting, removing NAS and renaming variables  
  dplyr::mutate(sp = as.numeric(factor(sp))) %>% #Sp to numeric  
  dplyr::select(sp, id, amount) #Three variables of interest selected

# puvsp: Species ID, PU ID and Species Presence
#Not saving this dataframe since it needs harmonising with the shark_conservation metrics 

# Prepare coordinates
shark.coordinates <- puvsp_marine %>%
  dplyr::transmute(
    FID = dplyr::row_number() - 1,
    id = dplyr::row_number(),
    Coord_x = lon,
    Coord_y = lat
  )

# pu_coords: id, lat, lon 
write.csv(shark.coordinates[,-1], here::here("Data", "My dataframes","pu_coords.csv"), row.names=FALSE)

#-------------- Fishing hours dataframe
library(tidyverse)
# Load the data
load("Data/Raw/predicted_fishing_everywhere_05Deg.Rdata")
head(prediction_data_env)
summary(prediction_data_env)

# Create the dataframe with pu and fishing hours
Fishing_preds <- prediction_data_env %>%
  dplyr::rename(Coord_x = 1, Coord_y = 2)

# Join the two dataframes together and process
pu_fishing <- shark.coordinates %>%
  dplyr::select(id, Coord_x, Coord_y) %>%
  dplyr::left_join(Fishing_preds, by = c("Coord_x", "Coord_y")) %>%
  dplyr::select(id, predicted_fishing_hours) %>%
  dplyr::mutate(
    # Replace NA values with 0.1
    predicted_fishing_hours = dplyr::coalesce(predicted_fishing_hours, 0.1),
    
    # Apply log1p (log(1 + x)) transformation
    predicted_fishing_hours = log1p(predicted_fishing_hours)
  )

summary(pu_fishing)

# pu_fishing: PU ID and log1p(fishing hours)
colnames(pu_fishing)[2]="cost"
write.csv(pu_fishing, here::here("Data", "My dataframes","pu_fishing_everywhere.csv"), row.names=FALSE)

# Create Species Biodiversity Index Values lists: shark_conservation_metrics.rds ----
# Read data
shark.info <- readRDS("Data/Raw/CAPTAIN_SharkInformation.rds")
shark.EDGE <- readRDS("Data/Raw/sharks_EDGE_final.rds")
shark.FUSE <- read.csv("Data/Raw/shark_sp_metrics.csv")
shark.IUCN <- read.csv("Data/Raw/species_updated.iucn_no.syn.csv")

# Function to process and join metrics
process_and_join_metric <- function(info_data, metric_data, metric_name) {
  metric_data %>%
    dplyr::select(Species, !!metric_name) %>%
    dplyr::mutate(
      Species = gsub("_", " ", Species),
      !!metric_name := vegan::decostand(!!sym(metric_name), "range")
    ) %>%
    {
      missing_species <- setdiff(info_data$Species, .$Species)
      message(paste("Missing species in", metric_name, ":", length(missing_species)))
      .
    } %>%
    dplyr::right_join(info_data, by = "Species")
}

# Process and join all metrics
shark.info <- shark.info %>%
  process_and_join_metric(shark.EDGE, "ED") %>%
  process_and_join_metric(shark.EDGE, "EDGE_P100") %>%
  process_and_join_metric(shark.FUSE, "FUSE") %>%
  process_and_join_metric(shark.FUSE, "EDGE2") %>%
  process_and_join_metric(shark.FUSE, "FUn") %>%
  process_and_join_metric(shark.FUSE, "FSp")

# Reorder columns
shark.info <- shark.info %>%
  dplyr::relocate(c(ED, EDGE_P100, EDGE2, FUn, FSp, FUSE), .after = IUCN)

# Display the first few rows of the final dataset
print(head(shark.info))

# Transforming IUCN status

# NE (Not Evaluated) and DD (Data Deficient) turn into NA
# Turn the others from 0 to 4 (LC, NT, VU, EN, CR)
shark.info <- shark.info %>%
  mutate(IUCN = case_when(
    IUCN %in% c("NE", "DD") ~ NA_real_,
    IUCN == "LC" ~ 0,
    IUCN == "NT" ~ 1,
    IUCN == "VU" ~ 2,
    IUCN == "EN" ~ 3,
    IUCN == "CR" ~ 4,
    TRUE ~ as.numeric(as.character(IUCN))
  )) %>%
  mutate(IUCN = decostand(IUCN, "range", na.rm = TRUE))

# Get the species names from puvsp_marine (columns 4 to the last)
species_in_puvsp_marine <- names(puvsp_marine)[4:ncol(puvsp_marine)]

# Filter shark.info to keep only the species found in puvsp_marine
shark.info <- shark.info %>%
  filter(Species %in% species_in_puvsp_marine) %>%
  mutate(SpeciesID = row_number())

# Print the number of species retained
print(paste("Number of species retained in shark.info:", nrow(shark.info)))

# Optional: Print the first few rows of the filtered shark.info
print(head(shark.info))

process_dataset <- function(shark.info, puvsp, metric) {
  # Drop NAs and extract species IDs
  info <- drop_na(shark.info, {{metric}})
  species_ids <- info$SpeciesID
  
  # Create transformed dataset using the metric values as amount
  transformed <- info %>%
    dplyr::select(SpeciesID, {{metric}}) %>%
    dplyr::rename(sp = SpeciesID, amount = {{metric}}) %>%
    # Scale the amount values so the maximum is 1
    dplyr::mutate(amount = amount / max(amount, na.rm = TRUE))
  
  return(list(info = info, species = species_ids, transformed = transformed))
}

# Process each dataset
metrics <- c("EDGE_P100", "FUSE", "EDGE2", "IUCN", "ED", "FSp", "FUn")
results <- lapply(metrics, function(metric) {
  process_dataset(shark.info, puvsp, metric)
})

names(results) <- metrics

# Define the file name
file_name <- here::here("Data", "My dataframes","shark_conservation_metrics_no_freshwater.rds")

# Save the results list
saveRDS(results, file = file_name)
jsonlite::write_json(results, here::here("Data", "My dataframes","shark_conservation_metrics_no_freshwater.json"))

# Access results like this in R:
# EDGE.info <- results$EDGE_P100$info
# EDGE.species <- results$EDGE_P100$species
EDGE.transformed <- results$EDGE_P100$transformed

# Access results like this in Python:
# EDGE_info = results['EDGE_P100']['info']
# EDGE_species = results['EDGE_P100']['species']
# EDGE_transformed = results['EDGE_P100']['transformed']



# Add the MPA data ----
# Function to read and filter MPA data
read_and_filter_mpa <- function(file_number) {
  file_path <- here("Data", "Raw", "WDPA_WDOECM_Sep2024_Public_marine_shp", 
                    paste0("WDPA_WDOECM_Sep2024_Public_marine_shp_", file_number),
                    "WDPA_WDOECM_Sep2024_Public_marine_shp-polygons.shp")
  
  tryCatch({
    mpa <- st_read(file_path)
    mpa_filtered <- mpa %>% filter(IUCN_CAT %in% c("Ia", "Ib", "II", "III"))
    list(all = mpa, filtered = mpa_filtered)
  }, error = function(e) {
    warning(paste("Error reading file", file_path, ":", e$message))
    NULL
  })
}

# Read and filter MPA data
mpa_data <- map(0:2, read_and_filter_mpa)

# Combine MPA data
mpa_NT <- do.call(rbind, map(mpa_data, "filtered"))
mpa_all <- do.call(rbind, map(mpa_data, "all"))

# Reset row names
row.names(mpa_NT) <- NULL
row.names(mpa_all) <- NULL

# Create raster layout for MPA data
rast <- raster(xmn=-180.25, xmx=180.75, ymn=-90.25, ymx=90.75, resolution=c(0.5,0.5))

# Rasterize MPA data
mpa_NT_raster <- rasterize(st_sf(mpa_NT$geometry), rast)
mpa_ALL_raster <- rasterize(st_sf(mpa_all$geometry), rast)

# Modified process_raster function
process_raster <- function(raster_data) {
  df <- as.data.frame(raster_data, xy=TRUE) %>%
    dplyr::mutate(layer = tidyr::replace_na(layer, 0),
                  layer = dplyr::if_else(layer > 0, 1, 0),
                  id = dplyr::row_number()) %>%
    dplyr::select(id, status = layer)
  return(df)
}

# Apply the function
mpa_NT_df_raster <- process_raster(mpa_NT_raster)
mpa_ALL_df_raster <- process_raster(mpa_ALL_raster)

# Join with pu_fishing data
pu_fishing_notake <- dplyr::left_join(pu_fishing, mpa_NT_df_raster, by = "id") %>%
  dplyr::select(id, status, cost)  # Swap columns 2 and 3

pu_fishing_allMPAs <- dplyr::left_join(pu_fishing, mpa_ALL_df_raster, by = "id") %>%
  dplyr::select(id, status, cost)  # Swap columns 2 and 3

# Save results
readr::write_csv(pu_fishing_notake, here("Data", "My dataframes","pu_fishing_everywhere_notake.csv"))
readr::write_csv(pu_fishing_allMPAs, here("Data", "My dataframes","pu_fishing_everywhere_allMPAs.csv"))

# Create binary rasters from the existing rasters
mpa_NT_binary <- mpa_NT_raster
mpa_NT_binary[!is.na(mpa_NT_binary)] <- 1
mpa_NT_binary[is.na(mpa_NT_binary)] <- 0

mpa_ALL_binary <- mpa_ALL_raster
mpa_ALL_binary[!is.na(mpa_ALL_binary)] <- 1
mpa_ALL_binary[is.na(mpa_ALL_binary)] <- 0

plot(mpa_NT_binary)
plot(mpa_ALL_binary)

# Save the binary rasters
writeRaster(mpa_NT_binary, filename=here::here("mpa_NT_binary.tif"), format="GTiff", overwrite=TRUE)
writeRaster(mpa_ALL_binary, filename=here::here("mpa_ALL_binary.tif"), format="GTiff", overwrite=TRUE)

# Harmonise datasets ----

#You have these datasets:
# pu_coords
# pu_fishing_notake and pu_fishing_allMPAs
# puvsp

load(here::here("Data", "puvsp_marine.Rdata"))

pu_coords <- read.csv("Data/My dataframes/pu_coords.csv")
head(pu_coords)

#puvsp <- read.csv("Data/My dataframes/puvsp.csv")
summary(puvsp)
head(puvsp)

pu_fishing_notake <- read.csv("Data/My dataframes/pu_fishing_everywhere_notake.csv")
head(pu_fishing_notake)

pu_fishing_allMPAs <- read.csv("Data/My dataframes/pu_fishing_everywhere_allMPAs.csv")
head(pu_fishing_allMPAs)

file_name <- here::here("Data", "My dataframes", "shark_conservation_metrics_no_freshwater.rds")
results <- readRDS(file = file_name)
summary(results$EDGE2$transformed)

EDGE2 <- results$EDGE2$transformed
summary(EDGE2)

#Harmonise the number of species in all datasets: puvsp and shark_conservation_metric
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(here)

# Create a mapping of species names to IDs for results
results_species_mapping <- shark.EDGE %>%
  dplyr::mutate(
    latin_name = stringr::str_replace_all(Species, "_", " "),
    sp = dplyr::row_number()
  ) %>%
  dplyr::select(sp, latin_name)

# Print the first few rows to verify
print(head(results_species_mapping))

# Create a mapping of species names to IDs for puvsp
puvsp_species_mapping <- tibble(
  latin_name = colnames(puvsp_marine)[4:ncol(puvsp_marine)],
  sp = 1:length(latin_name)
)

# Modify your existing puvsp transformation to include this mapping
puvsp <- puvsp_marine %>%
  dplyr::mutate(dplyr::across(4:ncol(.), ~ifelse(. == 0, NA, .))) %>%
  dplyr::mutate(id = dplyr::row_number()) %>%
  dplyr::select(-c(lon, lat)) %>%
  tidyr::pivot_longer(cols = -id, 
                      names_to = "latin_name", 
                      values_to = "amount", 
                      values_drop_na = TRUE) %>%
  dplyr::left_join(puvsp_species_mapping, by = "latin_name") %>%
  dplyr::select(sp, id, amount)

compare_and_reconcile_species <- function(results, puvsp, results_mapping, puvsp_mapping) {
  # Get all unique species IDs from results
  results_species <- unique(unlist(purrr::map(results, ~unique(.$transformed$sp))))
  
  # Get all unique species IDs from puvsp
  puvsp_species <- unique(puvsp$sp)
  
  # Create a comparison dataframe
  comparison <- dplyr::full_join(
    results_mapping %>% dplyr::filter(sp %in% results_species),
    puvsp_mapping %>% dplyr::filter(sp %in% puvsp_species),
    by = "latin_name", suffix = c("_results", "_puvsp")
  )
  
  # Find common species
  common_species <- comparison %>%
    dplyr::filter(!is.na(sp_results) & !is.na(sp_puvsp)) %>%
    dplyr::pull(sp_results)
  
  # Check for mismatches
  mismatches <- comparison %>%
    dplyr::filter(is.na(sp_results) | is.na(sp_puvsp))
  
  # Print summary
  cat("Total species in results:", length(results_species), "\n")
  cat("Total species in puvsp:", length(puvsp_species), "\n")
  cat("Common species:", length(common_species), "\n")
  cat("Mismatched species:", nrow(mismatches), "\n")
  
  # Print mismatches if any
  if (nrow(mismatches) > 0) {
    print(mismatches)
  }
  
  # Create a mapping between results sp and puvsp sp
  sp_mapping <- comparison %>%
    dplyr::filter(!is.na(sp_results) & !is.na(sp_puvsp)) %>%
    dplyr::select(sp_results, sp_puvsp)
  
  # Filter results and puvsp to include only common species
  results_filtered <- purrr::map(results, function(metric) {
    metric$transformed <- metric$transformed %>% 
      dplyr::filter(sp %in% common_species) %>%
      dplyr::left_join(sp_mapping, by = c("sp" = "sp_results")) %>%
      dplyr::mutate(sp = sp_puvsp) %>%
      dplyr::select(-sp_puvsp)
    metric$species <- intersect(metric$species, common_species)
    metric
  })
  
  puvsp_filtered <- puvsp %>% dplyr::filter(sp %in% sp_mapping$sp_puvsp)
  
  list(results = results_filtered, puvsp = puvsp_filtered, common_species = common_species, sp_mapping = sp_mapping)
}

# Use the function
reconciled <- compare_and_reconcile_species(results, puvsp, results_species_mapping, puvsp_species_mapping)

# Recreate the species mappings
results_species_mapping <- shark.EDGE %>%
  dplyr::mutate(
    latin_name = stringr::str_replace_all(Species, "_", " "),
    sp = dplyr::row_number()
  ) %>%
  dplyr::select(sp, latin_name)

puvsp_species_mapping <- tibble(
  latin_name = colnames(puvsp_marine)[4:ncol(puvsp_marine)],
  sp = 1:length(latin_name)
)

# Create the comparison dataframe
comparison <- dplyr::full_join(
  results_species_mapping,
  puvsp_species_mapping,
  by = "latin_name", suffix = c("_results", "_puvsp")
)

# Get the common species names
common_species_names <- comparison %>%
  dplyr::filter(!is.na(sp_results) & !is.na(sp_puvsp)) %>%
  dplyr::pull(latin_name)

# GRID data transformation
puvsp <- puvsp_marine %>%
  # Select only the columns we want (common species and the first three columns)
  dplyr::select(1:3, dplyr::all_of(common_species_names)) %>%
  dplyr::mutate(across(4:ncol(.), ~ifelse(. == 0, NA, .))) %>% # Replace all 0 values with NA in every column from the 4th to the last column
  dplyr::mutate(id = dplyr::row_number()) %>% #Adding a new column called "id" to the dataframe, filled with sequential numbers
  dplyr::select(-c(lon, lat)) %>% #Remove three columns
  tidyr::pivot_longer(cols = -id, 
                      names_to = "sp", 
                      values_to = "amount", 
                      values_drop_na = TRUE) %>% #Melting, removing NAs and renaming variables  
  dplyr::mutate(sp = as.numeric(factor(sp))) %>% #Sp to numeric  
  dplyr::select(sp, id, amount) #Three variables of interest selected

# Print the first few rows and summary to verify
print(head(puvsp))
print(summary(puvsp))
print(length(unique(puvsp$sp)))

write.csv(puvsp, here::here("Data", "My dataframes","puvsp_no_freshwater_harmonised.csv"), row.names=FALSE)

# Prepare coordinates
shark.coordinates <- puvsp_marine %>%
  dplyr::select(1:3, dplyr::all_of(common_species_names)) %>%
  dplyr::transmute(
    FID = dplyr::row_number() - 1,
    id = dplyr::row_number(),
    Coord_x = lon,
    Coord_y = lat
  )

# pu_coords: id, lat, lon 
write.csv(shark.coordinates[,-1], here::here("Data", "My dataframes","pu_coords_no_freshwater_harmonised.csv"), row.names=FALSE)

# Check all datasets ----
pu_fishing_data <- read.csv(here::here("Data", 
                                       "Theos datasets", 
                                       "Budget 30% all MPAs", 
                                       "pu_fishing_everywhere_allMPAs.csv"))
summary(pu_fishing_data)

puvsp <- read.csv(here::here("Data", 
                                       "Theos datasets", 
                                       "Budget 30% all MPAs", 
                                       "puvsp_no_freshwater_harmonised.csv"))
summary(puvsp)

pu_coords <- read.csv(here::here("Data", 
                             "Theos datasets", 
                             "Budget 30% all MPAs", 
                             "pu_coords_no_freshwater_harmonised.csv"))
summary(pu_coords)

shark_conservation_metrics <- jsonlite::fromJSON(here::here("Data", "My dataframes", "shark_conservation_metrics_no_freshwater.json"))

# View the structure of the data
str(shark_conservation_metrics)
shark_conservation_metrics$EDGE2$transformed
summary(as.data.frame(shark_conservation_metrics$EDGE2$transformed))
