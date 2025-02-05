library(sf)
library(terra)

# First database, species grid: pu_vsp ----
# Loading the data
load(here::here("Data", "puvsp_marine.Rdata"))

High_seas=st_read(here::here("Data", "World_High_Seas_v2_20241010", "High_Seas_v2.shp"))
#plot(st_geometry(High_seas))

# Create a template raster with the same resolution as your species data (0.5 degrees)
template_raster <- rast(ext(High_seas), resolution = 0.5)

# Rasterize the high seas polygon
highseas_raster <- rasterize(High_seas, template_raster, field = 1, background = 0)

# Plot to verify
plot(highseas_raster, 
     main = "High Seas (1) vs Continental Waters (0)",
     col = c("lightblue", "darkblue"))

# Extract values from the raster for each point in your species data
highseas_values <- extract(highseas_raster, puvsp_marine[, c("lon", "lat")])

# Create two separate datasets
highseas_species <- puvsp_marine[highseas_values[,2] == 1, ]
continental_species <- puvsp_marine[highseas_values[,2] == 0, ]

# Function to process the datasets
process_dataset <- function(df) {
  df %>%
    dplyr::mutate(across(4:ncol(.), ~ifelse(. == 0, NA, .))) %>%
    dplyr::mutate(id = dplyr::row_number()) %>%
    dplyr::select(-c(lon, lat)) %>%
    dplyr::filter(if_any(4:ncol(.), ~!is.na(.))) %>%
    tidyr::pivot_longer(cols = -id,
                        names_to = "sp",
                        values_to = "amount",
                        values_drop_na = TRUE) %>%
    dplyr::mutate(sp = as.numeric(factor(sp))) %>%
    dplyr::select(sp, id, amount)
}

# Apply the transformation to both datasets
highseas_puvsp <- process_dataset(highseas_species)
continental_puvsp <- process_dataset(continental_species)

#puvsp: Species ID, PU ID and Species Presence
#Not saving this dataframe since it needs harmonising with the shark_conservation metrics 

# Prepare coordinates ----

# Function to process coordinates
prepare_coordinates <- function(df, output_name) {
  coords <- df %>%
    dplyr::transmute(
      FID = dplyr::row_number() - 1,
      id = dplyr::row_number(),
      Coord_x = lon,
      Coord_y = lat
    )
  
  # Write the coordinates file (without FID)
  write.csv(coords[,-1], 
            here::here("Data", "My dataframes", paste0(output_name, "_pu_coords.csv")), 
            row.names = FALSE)
  
  return(coords)
}

# Process both datasets
highseas_coordinates <- prepare_coordinates(highseas_species, "highseas")
continental_coordinates <- prepare_coordinates(continental_species, "continental")

# Fishing datasets ---- 
# First create separate fishing predictions dataframes
load("Data/Raw/predicted_fishing_everywhere_05Deg.Rdata")

Fishing_preds <- prediction_data_env %>%
  dplyr::rename(Coord_x = 1, Coord_y = 2)

# Process for high seas
highseas_fishing <- highseas_coordinates %>%
  dplyr::select(id, Coord_x, Coord_y) %>%
  dplyr::left_join(Fishing_preds, by = c("Coord_x", "Coord_y")) %>%
  dplyr::select(id, predicted_fishing_hours) %>%
  dplyr::mutate(
    predicted_fishing_hours = dplyr::coalesce(predicted_fishing_hours, 0.1),
    predicted_fishing_hours = log1p(predicted_fishing_hours)
  )

# Calculate mean for high seas
highseas_mean_fishing <- mean(highseas_fishing$predicted_fishing_hours, na.rm = TRUE)

# Replace with mean and prepare for export
highseas_fishing_final <- highseas_fishing %>%
  dplyr::mutate(predicted_fishing_hours = highseas_mean_fishing) %>%
  dplyr::rename(cost = predicted_fishing_hours)

# Process for continental waters
continental_fishing <- continental_coordinates %>%
  dplyr::select(id, Coord_x, Coord_y) %>%
  dplyr::left_join(Fishing_preds, by = c("Coord_x", "Coord_y")) %>%
  dplyr::select(id, predicted_fishing_hours) %>%
  dplyr::mutate(
    predicted_fishing_hours = dplyr::coalesce(predicted_fishing_hours, 0.1),
    predicted_fishing_hours = log1p(predicted_fishing_hours)
  )

# Calculate mean for continental waters
continental_mean_fishing <- mean(continental_fishing$predicted_fishing_hours, na.rm = TRUE)

# Replace with mean and prepare for export
continental_fishing_final <- continental_fishing %>%
  dplyr::mutate(predicted_fishing_hours = continental_mean_fishing) %>%
  dplyr::rename(cost = predicted_fishing_hours)

# Save the files
write.csv(highseas_fishing_final, 
          here::here("Data", "My dataframes", "highseas_pu_fishing.csv"), 
          row.names = FALSE)

write.csv(continental_fishing_final, 
          here::here("Data", "My dataframes", "continental_pu_fishing.csv"), 
          row.names = FALSE)

# Print summary statistics
cat("High seas mean fishing hours (log1p):", highseas_mean_fishing, "\n")
cat("Continental waters mean fishing hours (log1p):", continental_mean_fishing, "\n")

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

# Transforming IUCN status
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
  mutate(IUCN = vegan::decostand(IUCN, "range", na.rm = TRUE))

# Filter to species with complete metrics
shark.info_complete <- shark.info %>%
  drop_na(ED, EDGE_P100, FUSE, EDGE2, FUn, FSp)

cat("Species with complete metrics:", nrow(shark.info_complete), "\n")

# Filter PUVSP data to match species in metrics
filter_puvsp_to_metrics <- function(puvsp_data, shark_info_complete) {
  # Get species columns (assuming columns 4+ are species)
  species_cols <- names(puvsp_data)[4:ncol(puvsp_data)]
  
  # Keep only species that are in shark_info_complete
  species_to_keep <- species_cols[species_cols %in% shark_info_complete$Species]
  
  # Select these species plus the coordinate columns
  filtered_puvsp <- puvsp_data %>%
    select(id, lon, lat, all_of(species_to_keep))
  
  cat("Original species count:", length(species_cols), "\n")
  cat("Species after filtering:", length(species_to_keep), "\n")
  
  return(filtered_puvsp)
}

# Filter both PUVSP datasets
highseas_species <- filter_puvsp_to_metrics(highseas_species, shark.info_complete)
continental_species <- filter_puvsp_to_metrics(continental_species, shark.info_complete)

# Then use shark.info_complete in process_zone_metrics
# Process for high seas
highseas_results <- process_zone_metrics(shark.info_complete, 
                                         highseas_species, 
                                         "highseas")

# Process for continental waters
continental_results <- process_zone_metrics(shark.info_complete, 
                                            continental_species, 
                                            "continental")
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
process_raster <- function(raster_data, coords_df) {
  # Convert coordinates to cell numbers
  cells <- cellFromXY(raster_data, 
                      xy = data.frame(x = coords_df$Coord_x, 
                                      y = coords_df$Coord_y))
  
  # Extract values for just those cells
  values <- raster::extract(raster_data, cells)
  
  # Create dataframe with the correct IDs
  df <- data.frame(
    id = coords_df$id,
    status = ifelse(is.na(values), 0, 1)
  )
  
  return(df)
}

# Process separately for high seas and continental waters
mpa_NT_df_highseas <- process_raster(mpa_NT_raster, highseas_coordinates)
mpa_NT_df_continental <- process_raster(mpa_NT_raster, continental_coordinates)

mpa_ALL_df_highseas <- process_raster(mpa_ALL_raster, highseas_coordinates)
mpa_ALL_df_continental <- process_raster(mpa_ALL_raster, continental_coordinates)

# Process for high seas
highseas_fishing_notake <- dplyr::left_join(highseas_fishing_final, mpa_NT_df_highseas, by = "id") %>%
  dplyr::select(id, status, cost)

highseas_fishing_allMPAs <- dplyr::left_join(highseas_fishing_final, mpa_ALL_df_highseas, by = "id") %>%
  dplyr::select(id, status, cost)

# Process for continental waters
continental_fishing_notake <- dplyr::left_join(continental_fishing_final, mpa_NT_df_continental, by = "id") %>%
  dplyr::select(id, status, cost)

continental_fishing_allMPAs <- dplyr::left_join(continental_fishing_final, mpa_ALL_df_continental, by = "id") %>%
  dplyr::select(id, status, cost)

# Save results for high seas
readr::write_csv(highseas_fishing_notake, 
                 here("Data", "My dataframes", "highseas_fishing_notake.csv"))
readr::write_csv(highseas_fishing_allMPAs, 
                 here("Data", "My dataframes", "highseas_fishing_allMPAs.csv"))

# Save results for continental waters 
readr::write_csv(continental_fishing_notake, 
                 here("Data", "My dataframes", "continental_fishing_notake.csv"))
readr::write_csv(continental_fishing_allMPAs, 
                 here("Data", "My dataframes", "continental_fishing_allMPAs.csv"))

# Create binary rasters from the existing rasters
mpa_NT_binary <- mpa_NT_raster
mpa_NT_binary[!is.na(mpa_NT_binary)] <- 1
mpa_NT_binary[is.na(mpa_NT_binary)] <- 0

mpa_ALL_binary <- mpa_ALL_raster
mpa_ALL_binary[!is.na(mpa_ALL_binary)] <- 1
mpa_ALL_binary[is.na(mpa_ALL_binary)] <- 0

# Plot the binary rasters
plot(mpa_NT_binary, main="No-Take MPAs")
plot(mpa_ALL_binary, main="All MPAs")

# Create masks using coordinates
create_zone_mask <- function(rast_template, coords_df) {
  # Create empty raster with same properties as template
  zone_mask <- rast_template
  zone_mask[] <- 0
  
  # Convert coordinates to cell numbers in the raster
  cells <- terra::cellFromXY(zone_mask, 
                             xy = data.frame(x = coords_df$Coord_x, 
                                             y = coords_df$Coord_y))
  
  # Set those cells to 1
  zone_mask[cells] <- 1
  
  return(zone_mask)
}

# Create masks for each zone using the coordinates
highseas_mask <- create_zone_mask(rast, highseas_coordinates)
continental_mask <- create_zone_mask(rast, continental_coordinates)

# Create binary rasters for high seas
highseas_NT_binary <- mpa_NT_raster * highseas_mask
highseas_ALL_binary <- mpa_ALL_raster * highseas_mask

# Create binary rasters for continental waters
continental_NT_binary <- mpa_NT_raster * continental_mask
continental_ALL_binary <- mpa_ALL_raster * continental_mask

# Convert to binary (0/1)
highseas_NT_binary[highseas_NT_binary > 0] <- 1
highseas_NT_binary[is.na(highseas_NT_binary)] <- 0

highseas_ALL_binary[highseas_ALL_binary > 0] <- 1
highseas_ALL_binary[is.na(highseas_ALL_binary)] <- 0

continental_NT_binary[continental_NT_binary > 0] <- 1
continental_NT_binary[is.na(continental_NT_binary)] <- 0

continental_ALL_binary[continental_ALL_binary > 0] <- 1
continental_ALL_binary[is.na(continental_ALL_binary)] <- 0

# Plot to verify
par(mfrow=c(2,2))
plot(highseas_NT_binary, main="High Seas No-Take MPAs")
plot(highseas_ALL_binary, main="High Seas All MPAs")
plot(continental_NT_binary, main="Continental No-Take MPAs")
plot(continental_ALL_binary, main="Continental All MPAs")
par(mfrow=c(1,1))

# Save the binary rasters
writeRaster(highseas_NT_binary, 
            filename=here::here("Data", "My dataframes", "highseas_mpa_NT_binary.tif"), 
            format="GTiff", 
            overwrite=TRUE)
writeRaster(highseas_ALL_binary, 
            filename=here::here("Data", "My dataframes", "highseas_mpa_ALL_binary.tif"), 
            format="GTiff", 
            overwrite=TRUE)

writeRaster(continental_NT_binary, 
            filename=here::here("Data", "My dataframes", "continental_mpa_NT_binary.tif"), 
            format="GTiff", 
            overwrite=TRUE)
writeRaster(continental_ALL_binary, 
            filename=here::here("Data", "My dataframes", "continental_mpa_ALL_binary.tif"), 
            format="GTiff", 
            overwrite=TRUE)

# Harmonise the number of species in all datasets: puvsp and shark_conservation_metric ----
library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(here)

# Load the existing files
highseas_metrics <- readRDS(here::here("Data", "My dataframes", "highseas_shark_conservation_metrics.rds"))
continental_metrics <- readRDS(here::here("Data", "My dataframes", "continental_shark_conservation_metrics.rds"))

# Function to harmonize metrics with PUVSP species
harmonize_data <- function(metrics_data, puvsp_data, zone_name) {
  if(zone_name == "highseas") {
    # For highseas: filter metrics to match PUVSP species
    puvsp_species_ids <- unique(puvsp_data$sp)
    
    harmonized_metrics <- lapply(metrics_data, function(metric) {
      metric$transformed <- metric$transformed %>%
        filter(sp %in% puvsp_species_ids)
      
      metric$species <- metric$species[metric$species %in% puvsp_species_ids]
      metric$info <- metric$info %>%
        filter(SpeciesID %in% puvsp_species_ids)
      
      return(metric)
    })
    
    harmonized_puvsp <- puvsp_data  # No change needed for highseas PUVSP
    
  } else {
    # For continental: filter PUVSP to match metrics species
    metrics_species_ids <- unique(metrics_data$EDGE_P100$transformed$sp)
    
    harmonized_metrics <- metrics_data  # No change needed for continental metrics
    harmonized_puvsp <- puvsp_data %>%
      filter(sp %in% metrics_species_ids)
  }
  
  # Save harmonized metrics
  saveRDS(harmonized_metrics, 
          here::here("Data", "My dataframes", 
                     paste0(zone_name, "_shark_conservation_metrics_harmonised.rds")))
  
  jsonlite::write_json(harmonized_metrics,
                       here::here("Data", "My dataframes", 
                                  paste0(zone_name, "_shark_conservation_metrics_harmonised.json")))
  
  # Save harmonized PUVSP
  write.csv(harmonized_puvsp,
            here::here("Data", "My dataframes", 
                       paste0(zone_name, "_puvsp_harmonised.csv")),
            row.names = FALSE)
  
  # Return summary
  return(list(
    n_species = if(zone_name == "highseas") length(puvsp_species_ids) else length(metrics_species_ids),
    species_ids = if(zone_name == "highseas") puvsp_species_ids else metrics_species_ids
  ))
}

# Harmonize high seas data
highseas_harmonized <- harmonize_data(
  highseas_metrics,
  highseas_puvsp,
  "highseas"
)

# Harmonize continental data
continental_harmonized <- harmonize_data(
  continental_metrics,
  continental_puvsp,
  "continental"
)

# Print summaries
cat("\nHigh Seas Harmonized Summary:\n")
cat("Number of species:", highseas_harmonized$n_species, "\n")

cat("\nContinental Waters Harmonized Summary:\n")
cat("Number of species:", continental_harmonized$n_species, "\n")

# Verify that species match between metrics and PUVSP for each zone
verify_species_match <- function(zone_name) {
  metrics <- readRDS(here::here("Data", "My dataframes", 
                                paste0(zone_name, "_shark_conservation_metrics_harmonised.rds")))
  puvsp <- read.csv(here::here("Data", "My dataframes", 
                               paste0(zone_name, "_puvsp_harmonised.csv")))
  
  metrics_species <- unique(metrics$EDGE_P100$transformed$sp)
  puvsp_species <- unique(puvsp$sp)
  
  cat("\nVerification for", zone_name, ":\n")
  cat("Species in metrics:", length(metrics_species), "\n")
  cat("Species in PUVSP:", length(puvsp_species), "\n")
  cat("Matching species:", length(intersect(metrics_species, puvsp_species)), "\n")
}

verify_species_match("highseas")
verify_species_match("continental")

# Make sure that we've removed the correct species ----
# Function to extract and compare species details
compare_species_details <- function(zone_name) {
  # Load datasets
  metrics <- readRDS(here::here("Data", "My dataframes", 
                                paste0(zone_name, "_shark_conservation_metrics_harmonised.rds")))
  puvsp <- read.csv(here::here("Data", "My dataframes", 
                               paste0(zone_name, "_puvsp_harmonised.csv")))
  
  # Get species from metrics (including names if available)
  metrics_species <- metrics$EDGE_P100$info %>%
    select(SpeciesID, Species) %>%
    distinct()
  
  # Get unique species from PUVSP
  puvsp_species <- data.frame(SpeciesID = unique(puvsp$sp))
  
  # Compare sets
  all_species <- full_join(
    metrics_species %>% mutate(in_metrics = TRUE),
    puvsp_species %>% mutate(in_puvsp = TRUE),
    by = "SpeciesID"
  ) %>%
    mutate(
      in_metrics = ifelse(is.na(in_metrics), FALSE, in_metrics),
      in_puvsp = ifelse(is.na(in_puvsp), FALSE, in_puvsp)
    )
  
  # Save comparison to CSV
  write.csv(all_species,
            here::here("Data", "My dataframes", 
                       paste0(zone_name, "_species_comparison.csv")),
            row.names = FALSE)
  
  # Print summary
  cat("\nSpecies comparison for", zone_name, ":\n")
  cat("Total unique species:", nrow(all_species), "\n")
  cat("Species in both datasets:", sum(all_species$in_metrics & all_species$in_puvsp), "\n")
  cat("Species only in metrics:", sum(all_species$in_metrics & !all_species$in_puvsp), "\n")
  cat("Species only in PUVSP:", sum(!all_species$in_metrics & all_species$in_puvsp), "\n")
  
  # Show examples of mismatches
  cat("\nExamples of mismatches:\n")
  print(head(all_species %>% 
               filter(xor(in_metrics, in_puvsp)) %>%
               arrange(SpeciesID)))
  
  # Return the comparison dataframe
  return(all_species)
}

# Run comparison for both zones
highseas_comparison <- compare_species_details("highseas")
continental_comparison <- compare_species_details("continental")

# Function to check species ID and name correspondence
check_species_names <- function(zone_name) {
  # Load datasets
  metrics <- readRDS(here::here("Data", "My dataframes", 
                                paste0(zone_name, "_shark_conservation_metrics_harmonised.rds")))
  puvsp <- read.csv(here::here("Data", "My dataframes", 
                               paste0(zone_name, "_puvsp_harmonised.csv")))
  
  # Get species info from metrics
  metrics_species <- metrics$EDGE_P100$info %>%
    select(SpeciesID, Species) %>%
    distinct()
  
  # Get unique species IDs from PUVSP
  puvsp_species <- data.frame(SpeciesID = unique(puvsp$sp))
  
  # Compare
  comparison <- full_join(
    metrics_species,
    puvsp_species,
    by = "SpeciesID"
  )
  
  # Print summary
  cat("\nSpecies ID and name check for", zone_name, ":\n")
  cat("Total species IDs:", nrow(comparison), "\n")
  cat("\nExample of species matches (first 10):\n")
  print(head(comparison, 10))
  
  return(comparison)
}

# Run checks for both zones
highseas_names <- check_species_names("highseas")
continental_names <- check_species_names("continental")

# Function to compare specific ID ranges
compare_id_range <- function(zone1_name, zone2_name, id_start, id_end) {
  # Load datasets
  metrics1 <- readRDS(here::here("Data", "My dataframes", 
                                 paste0(zone1_name, "_shark_conservation_metrics_harmonised.rds")))
  metrics2 <- readRDS(here::here("Data", "My dataframes", 
                                 paste0(zone2_name, "_shark_conservation_metrics_harmonised.rds")))
  
  # Extract species info
  species1 <- metrics1$EDGE_P100$info %>%
    select(SpeciesID, Species) %>%
    distinct() %>%
    filter(SpeciesID >= id_start & SpeciesID <= id_end) %>%
    arrange(SpeciesID)
  
  species2 <- metrics2$EDGE_P100$info %>%
    select(SpeciesID, Species) %>%
    distinct() %>%
    filter(SpeciesID >= id_start & SpeciesID <= id_end) %>%
    arrange(SpeciesID)
  
  # Create comparison dataframe
  comparison <- full_join(
    species1 %>% rename(Highseas_Species = Species),
    species2 %>% rename(Continental_Species = Species),
    by = "SpeciesID"
  ) %>%
    mutate(Match = Highseas_Species == Continental_Species)
  
  return(comparison)
}

# Compare species 180-200
comparison <- compare_id_range("highseas", "continental", 180, 200)
print(comparison)

check_species_match <- function(zone_name) {
  # Load datasets
  metrics <- readRDS(here::here("Data", "My dataframes", 
                                paste0(zone_name, "_shark_conservation_metrics_harmonised.rds")))
  puvsp <- read.csv(here::here("Data", "My dataframes", 
                               paste0(zone_name, "_puvsp_harmonised.csv")))
  
  # Get unique species from each dataset
  metrics_species <- metrics$EDGE_P100$info %>%
    select(SpeciesID, Species) %>%
    distinct()
  
  puvsp_species <- data.frame(
    SpeciesID = unique(puvsp$sp)
  )
  
  # Compare presence in both datasets
  comparison <- full_join(
    metrics_species %>% 
      mutate(in_metrics = TRUE),
    puvsp_species %>% 
      mutate(in_puvsp = TRUE),
    by = "SpeciesID"
  ) %>%
    mutate(
      in_metrics = replace_na(in_metrics, FALSE),
      in_puvsp = replace_na(in_puvsp, FALSE)
    )
  
  # Print summary
  cat("\nSpecies match analysis for", zone_name, ":\n")
  cat("Total species in metrics:", sum(comparison$in_metrics), "\n")
  cat("Total species in PUVSP:", sum(comparison$in_puvsp), "\n")
  cat("Species in both datasets:", sum(comparison$in_metrics & comparison$in_puvsp), "\n")
  cat("Species only in metrics:", sum(comparison$in_metrics & !comparison$in_puvsp), "\n")
  cat("Species only in PUVSP:", sum(!comparison$in_metrics & comparison$in_puvsp), "\n\n")
  
  # Show species that are only in one dataset
  cat("First few species only in PUVSP (if any):\n")
  print(head(comparison[!comparison$in_metrics & comparison$in_puvsp, ], 5))
  cat("\nFirst few species only in metrics (if any):\n")
  print(head(comparison[comparison$in_metrics & !comparison$in_puvsp, ], 5))
  
  return(comparison)
}

# Check both zones
highseas_comparison <- check_species_match("highseas")
continental_comparison <- check_species_match("continental")

# Bind conservation metrics to a minimum of 10% and 30% ----
create_budget_scenarios <- function(zone_name) {
  # Load metrics - adjusting path to match your structure
  metrics <- jsonlite::fromJSON(here::here("Data", "My dataframes", 
                                           paste0(zone_name, "_shark_conservation_metrics_harmonised.json")))
  
  # Create 10% scenario
  metrics_10 <- metrics
  for(metric_name in names(metrics_10)) {
    metrics_10[[metric_name]]$transformed$amount <- pmax(metrics_10[[metric_name]]$transformed$amount, 0.1)
  }
  
  # Create 30% scenario
  metrics_30 <- metrics
  for(metric_name in names(metrics_30)) {
    metrics_30[[metric_name]]$transformed$amount <- pmax(metrics_30[[metric_name]]$transformed$amount, 0.3)
  }
  
  # Save the new datasets - adjusting path
  jsonlite::write_json(metrics_10, here::here("Data", "My dataframes", 
                                              paste0(zone_name, "_shark_conservation_metrics_10_harmonised.json")))
  jsonlite::write_json(metrics_30, here::here("Data", "My dataframes", 
                                              paste0(zone_name, "_shark_conservation_metrics_30_harmonised.json")))
}

# Create scenarios for both zones
create_budget_scenarios("continental")
create_budget_scenarios("highseas")

check_metrics <- function(zone_name) {
  metrics_10 <- jsonlite::fromJSON(here::here("Data", "My dataframes", 
                                              paste0(zone_name, "_shark_conservation_metrics_10.json")),
                                   simplifyDataFrame = TRUE)
  metrics_30 <- jsonlite::fromJSON(here::here("Data", "My dataframes", 
                                              paste0(zone_name, "_shark_conservation_metrics_30.json")),
                                   simplifyDataFrame = TRUE)
  
  # Function to get summary statistics for each metric
  summarize_metric <- function(metric_data) {
    # Convert to numeric if necessary and handle NAs
    amount <- as.numeric(metric_data$transformed$amount)
    amount <- amount[!is.na(amount)]
    
    return(data.frame(
      min = min(amount),
      median = median(amount),
      mean = round(mean(amount), 4),
      max = max(amount),
      n_below_threshold = sum(amount < 0.1)
    ))
  }
  
  # Create summary tables
  cat("\nSummary for", zone_name, ":\n")
  
  cat("\n10% Budget Scenario:\n")
  summary_10 <- do.call(rbind, lapply(names(metrics_10), function(metric) {
    stats <- summarize_metric(metrics_10[[metric]])
    cbind(metric = metric, stats)
  }))
  rownames(summary_10) <- NULL
  print(summary_10)
  
  cat("\n30% Budget Scenario:\n")
  summary_30 <- do.call(rbind, lapply(names(metrics_30), function(metric) {
    stats <- summarize_metric(metrics_30[[metric]])
    cbind(metric = metric, stats)
  }))
  rownames(summary_30) <- NULL
  print(summary_30)
}

# Check metrics
check_metrics("continental")
check_metrics("highseas")

# Check all datasets ----
pu_fishing_data <- read.csv(here::here("Data", 
                                       "My Dataframes", 
                                       "Continental", 
                                       "continental_fishing_allMPAs.csv"))
summary(pu_fishing_data)

puvsp <- read.csv(here::here("Data", 
                             "My Dataframes", 
                             "Continental", 
                             "continental_puvsp_harmonised.csv"))
summary(puvsp)
length(unique(puvsp$sp))
pu_coords <- read.csv(here::here("Data", 
                                 "My Dataframes", 
                                 "Continental", 
                                 "continental_coords_harmonised.csv"))
summary(pu_coords)

shark_conservation_metrics <- jsonlite::fromJSON(here::here("Data", "My dataframes", "Continental", 
                                                            "continental_shark_conservation_metrics.json"))

# View the structure of the data
str(shark_conservation_metrics)
length(shark_conservation_metrics$EDGE2$transformed[,1])
summary(as.data.frame(shark_conservation_metrics$EDGE2$transformed))
