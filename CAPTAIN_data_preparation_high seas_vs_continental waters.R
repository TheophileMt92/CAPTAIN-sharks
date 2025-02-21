library(sf)
library(terra)
library(rlang)
library(here)
library(tidyverse)

# First database, species grid: pu_vsp ----
# Loading the data
load(here::here("Data", "puvsp_marine.Rdata"))

# Create global species lookup table
global_species_lookup <- data.frame(
  species_name = names(puvsp_marine)[4:ncol(puvsp_marine)],
  species_id = 1:length(names(puvsp_marine)[4:ncol(puvsp_marine)])
)

High_seas = st_read(here::here("Data", "World_High_Seas_v2_20241010", "High_Seas_v2.shp"))

# Create template raster
template_raster <- rast(ext(High_seas), resolution = 0.5)

# Rasterize high seas polygon
highseas_raster <- rasterize(High_seas, template_raster, field = 1, background = 0)

# Plot to verify
plot(highseas_raster, 
     main = "High Seas (1) vs Continental Waters (0)",
     col = c("lightblue", "darkblue"))

# Convert points to SpatVector for extraction
points_coords <- vect(puvsp_marine, geom = c("lon", "lat"), crs = crs(highseas_raster))

# Extract values
highseas_values <- terra::extract(highseas_raster, points_coords)

# Create two separate datasets
highseas_species <- puvsp_marine[highseas_values[,2] == 1, ]
continental_species <- puvsp_marine[highseas_values[,2] == 0, ]

# Modified function to process the datasets using global lookup
process_dataset <- function(df, species_lookup) {
  # Print original species count
  cat("Original unique species count:", length(unique(names(df)[-(1:3)])), "\n")
  
  df_long <- df %>%
    dplyr::mutate(across(4:ncol(.), ~ifelse(. == 0, NA, .))) %>%
    dplyr::mutate(id = dplyr::row_number()) %>%
    dplyr::select(-c(lon, lat)) %>%
    dplyr::filter(if_any(4:ncol(.), ~!is.na(.))) %>%
    tidyr::pivot_longer(cols = -id,
                        names_to = "species_name",
                        values_to = "amount",
                        values_drop_na = TRUE)
  
  # Print species names before join
  cat("\nSample of species names before join:\n")
  print(head(unique(df_long$species_name)))
  
  # Join with lookup table
  result <- df_long %>%
    dplyr::left_join(species_lookup, by = "species_name") %>%
    dplyr::select(sp = species_id, id, amount)
  
  # Print summary after join
  cat("\nSpecies IDs after join:\n")
  print(summary(result$sp))
  
  return(result)
}

# Apply the transformation to both datasets using the same lookup table
highseas_puvsp <- process_dataset(highseas_species, global_species_lookup)
continental_puvsp <- process_dataset(continental_species, global_species_lookup)

# Save the global species lookup for use in subsequent steps
saveRDS(global_species_lookup, 
        here::here("Data", "My dataframes", "global_species_lookup.rds"))

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
# Load global species lookup created in first part
global_species_lookup <- readRDS(here::here("Data", "My dataframes", "global_species_lookup.rds"))

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

# Modified filter_puvsp_to_metrics function
filter_puvsp_to_metrics <- function(puvsp_data, shark_info_complete, species_lookup) {
  # Get species columns (assuming columns 4+ are species)
  species_cols <- names(puvsp_data)[4:ncol(puvsp_data)]
  
  # Keep only species that are in shark_info_complete and lookup
  species_to_keep <- species_cols[species_cols %in% shark.info_complete$Species]
  
  # Select these species plus the coordinate columns
  filtered_puvsp <- puvsp_data %>%
    select(id, lon, lat, all_of(species_to_keep))
  
  cat("Original species count:", length(species_cols), "\n")
  cat("Species after filtering:", length(species_to_keep), "\n")
  
  return(filtered_puvsp)
}

process_zone_metrics <- function(shark.info, puvsp_zone, zone_name, species_lookup) {
  species_in_zone <- names(puvsp_zone)[4:ncol(puvsp_zone)]
  
  zone_shark_info <- shark.info %>%
    filter(Species %in% species_in_zone) %>%
    left_join(species_lookup, by = c("Species" = "species_name"))
  
  print(paste("Number of species retained in", zone_name, ":", nrow(zone_shark_info)))
  
  # Add diagnostic check for IUCN values
  cat("\nIUCN value check before processing:\n")
  print(summary(zone_shark_info$IUCN))
  
  process_dataset <- function(shark.info, puvsp, metric) {
    # Drop NAs and extract species IDs
    info <- drop_na(shark.info, all_of(metric))
    species_ids <- info$species_id
    
    cat("\nProcessing metric:", metric, "\n")
    cat("Number of species before NA removal:", nrow(shark.info), "\n")
    cat("Number of species after NA removal:", nrow(info), "\n")
    
    transformed <- info %>%
      dplyr::select(species_id, all_of(metric)) %>%
      dplyr::rename(sp = species_id, 
                    amount = !!sym(metric))
    
    # The values should already be normalized from earlier processing
    # but we'll check and normalize if needed
    if(max(transformed$amount, na.rm = TRUE) > 1) {
      transformed <- transformed %>%
        dplyr::mutate(amount = amount / max(amount, na.rm = TRUE))
    }
    
    cat("Summary of transformed values:\n")
    print(summary(transformed$amount))
    
    return(list(
      info = info,
      species = species_ids,
      transformed = transformed,
      species_names = info$Species
    ))
  }
  
  metrics <- c("EDGE_P100", "FUSE", "EDGE2", "IUCN", "ED", "FSp", "FUn")
  results <- lapply(metrics, function(metric) {
    process_dataset(zone_shark_info, puvsp, metric)
  })
  
  names(results) <- metrics
  
  # Save results
  output_data <- list(
    metrics = results,
    species_lookup = species_lookup[species_lookup$species_name %in% species_in_zone, ]
  )
  
  saveRDS(output_data, 
          here::here("Data", "My dataframes", 
                     paste0(zone_name, "_shark_conservation_metrics.rds")))
  
  jsonlite::write_json(output_data,
                       here::here("Data", "My dataframes", 
                                  paste0(zone_name, "_shark_conservation_metrics.json")))
  
  return(results)
}

# Filter both PUVSP datasets
highseas_species <- filter_puvsp_to_metrics(highseas_species, shark.info_complete, global_species_lookup)
continental_species <- filter_puvsp_to_metrics(continental_species, shark.info_complete, global_species_lookup)

#You've lost ten species, find which ones:
# Get the species from shark.info_complete
shark_species <- shark.info_complete$Species

# Get the species from your presence data (excluding the first 3 columns which are id, lon, lat)
presence_species <- names(highseas_species)[4:ncol(highseas_species)]

# Find species that are in shark.info_complete but not in presence data
missing_species <- setdiff(shark_species, presence_species)

# Show the missing species
print("Species with metrics but no presence data:")
print(missing_species)

# Process for high seas
highseas_results <- process_zone_metrics(shark.info_complete, 
                                         highseas_species, 
                                         "highseas",
                                         global_species_lookup)

# Process for continental waters
continental_results <- process_zone_metrics(shark.info_complete, 
                                            continental_species, 
                                            "continental",
                                            global_species_lookup)
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
rast <- raster::raster(xmn=-180.25, xmx=180.75, ymn=-90.25, ymx=90.75, resolution=c(0.5,0.5))

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

# Load the existing files and global species lookup
highseas_data <- readRDS(here::here("Data", "My dataframes", "highseas_shark_conservation_metrics.rds"))
continental_data <- readRDS(here::here("Data", "My dataframes", "continental_shark_conservation_metrics.rds"))
global_species_lookup <- readRDS(here::here("Data", "My dataframes", "global_species_lookup.rds"))

# Extract metrics and lookup from the saved data
highseas_metrics <- highseas_data$metrics
continental_metrics <- continental_data$metrics

# Function to harmonize metrics with PUVSP species
harmonize_data <- function(metrics_data, puvsp_data, zone_name, species_lookup) {
  # Print initial counts
  cat("\nInitial counts for", zone_name, ":\n")
  cat("Metrics species:", length(unique(metrics_data$EDGE_P100$transformed$sp)), "\n")
  cat("PUVSP species:", length(unique(puvsp_data$sp)), "\n")
  
  # Get metrics species IDs
  metrics_species_ids <- unique(metrics_data$EDGE_P100$transformed$sp)
  
  # First, harmonize PUVSP to only include species that exist in metrics
  harmonized_puvsp <- puvsp_data %>%
    filter(sp %in% metrics_species_ids)
  
  # Get the final list of matching species
  final_species_ids <- unique(harmonized_puvsp$sp)
  
  cat("Matching species after PUVSP harmonization:", length(final_species_ids), "\n")
  
  # Then harmonize metrics to match the filtered PUVSP species
  harmonized_metrics <- lapply(metrics_data, function(metric) {
    metric$transformed <- metric$transformed %>%
      filter(sp %in% final_species_ids)
    
    metric$species <- metric$species[metric$species %in% final_species_ids]
    metric$info <- metric$info %>%
      filter(species_id %in% final_species_ids)
    
    # Add species names to transformed data
    metric$transformed <- metric$transformed %>%
      left_join(species_lookup, by = c("sp" = "species_id"))
    
    return(metric)
  })
  
  # Verify final counts
  cat("Final counts:\n")
  cat("Harmonized metrics species:", length(unique(harmonized_metrics$EDGE_P100$transformed$sp)), "\n")
  cat("Harmonized PUVSP species:", length(unique(harmonized_puvsp$sp)), "\n")
  
  # Create harmonized data object
  harmonized_data <- list(
    metrics = harmonized_metrics,
    puvsp = harmonized_puvsp,
    species_lookup = species_lookup %>%
      filter(species_id %in% final_species_ids)
  )
  
  # Save harmonized data
  saveRDS(harmonized_data, 
          here::here("Data", "My dataframes", 
                     paste0(zone_name, "_shark_conservation_metrics_harmonised.rds")))
  
  # Save harmonized PUVSP with species names
  harmonized_puvsp %>%
    left_join(species_lookup, by = c("sp" = "species_id")) %>%
    write.csv(here::here("Data", "My dataframes", 
                         paste0(zone_name, "_puvsp_harmonised.csv")),
              row.names = FALSE)
  
  # Return summary
  species_summary <- species_lookup %>%
    filter(species_id %in% final_species_ids)
  
  return(list(
    n_species = nrow(species_summary),
    species_ids = species_summary$species_id,
    species_names = species_summary$species_name
  ))
}

# Re-run the harmonization
highseas_harmonized <- harmonize_data(
  highseas_metrics,
  highseas_puvsp,
  "highseas",
  global_species_lookup
)

# Harmonize continental data
continental_harmonized <- harmonize_data(
  continental_metrics,
  continental_puvsp,
  "continental",
  global_species_lookup
)

# Print detailed summaries
cat("\nHigh Seas Harmonized Summary:\n")
cat("Number of species:", highseas_harmonized$n_species, "\n")
cat("Species names:\n")
print(data.frame(
  ID = highseas_harmonized$species_ids,
  Species = highseas_harmonized$species_names
))

cat("\nContinental Waters Harmonized Summary:\n")
cat("Number of species:", continental_harmonized$n_species, "\n")
cat("Species names:\n")
print(data.frame(
  ID = continental_harmonized$species_ids,
  Species = continental_harmonized$species_names
))

# Verify that species match between metrics and PUVSP for each zone
verify_species_match <- function(zone_name) {
  # Load harmonized data
  metrics <- readRDS(here::here("Data", "My dataframes", 
                                paste0(zone_name, "_shark_conservation_metrics_harmonised.rds")))
  puvsp <- read.csv(here::here("Data", "My dataframes", 
                               paste0(zone_name, "_puvsp_harmonised.csv")))
  
  # Access species IDs from the harmonized metrics
  metrics_species <- unique(metrics$metrics$EDGE_P100$transformed$sp)
  puvsp_species <- unique(puvsp$sp)
  
  cat("\nVerification for", zone_name, ":\n")
  cat("Species in harmonized metrics:", length(metrics_species), "\n")
  cat("Species in harmonized PUVSP:", length(puvsp_species), "\n")
  cat("Species in metrics (first few):", paste(head(metrics_species), collapse=", "), "\n")
  cat("Species in PUVSP (first few):", paste(head(puvsp_species), collapse=", "), "\n")
  cat("Matching species:", length(intersect(metrics_species, puvsp_species)), "\n")
  
  return(list(
    metrics_species = metrics_species,
    puvsp_species = puvsp_species
  ))
}

# Run verification on harmonized data
verify_species_match("highseas")
verify_species_match("continental")

# Bind conservation metrics to a minimum of 10% and 30% ----
create_budget_scenarios <- function(zone_name) {
  # Load harmonized metrics from RDS file
  metrics <- readRDS(here::here("Data", "My dataframes", 
                                paste0(zone_name, "_shark_conservation_metrics_harmonised.rds")))
  
  # Extract just the metrics part from the loaded data
  metrics <- metrics$metrics
  
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
  
  # Save as RDS files
  saveRDS(metrics_10, here::here("Data", "My dataframes", 
                                 paste0(zone_name, "_shark_conservation_metrics_10_harmonised.rds")))
  saveRDS(metrics_30, here::here("Data", "My dataframes", 
                                 paste0(zone_name, "_shark_conservation_metrics_30_harmonised.rds")))
  
  # Also save as JSON if needed
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

shark_conservation_metrics <- jsonlite::fromJSON(here::here("Data", "My dataframes", #"Continental", 
                                                            "continental_shark_conservation_metrics_10_harmonised.json"))

# View the structure of the data
str(shark_conservation_metrics)
length(shark_conservation_metrics$EDGE2$transformed[,1])
summary(as.data.frame(shark_conservation_metrics$metrics$EDGE2$transformed))

#
pu_fishing_data <- read.csv(here::here("Data", 
                                       "My Dataframes", 
                                       "High Seas", 
                                       "highseas_fishing_allMPAs.csv"))
summary(pu_fishing_data)

puvsp <- read.csv(here::here("Data", 
                             "My Dataframes", 
                             #"High Seas", 
                             "highseas_puvsp_harmonised.csv"))
summary(puvsp)
length(unique(puvsp$sp))
pu_coords <- read.csv(here::here("Data", 
                                 "My Dataframes", 
                                 "High Seas", 
                                 "highseas_pu_coords.csv"))
summary(pu_coords)

shark_conservation_metrics <- jsonlite::fromJSON(here::here("Data", "My dataframes", "High Seas", 
                                                            "highseas_shark_conservation_metrics_30_harmonised.json"))

# View the structure of the data
str(shark_conservation_metrics)
length(shark_conservation_metrics$metrics$EDGE2$transformed[,1])
summary(as.data.frame(shark_conservation_metrics$EDGE2$transformed))