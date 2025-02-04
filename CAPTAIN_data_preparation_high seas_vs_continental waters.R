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
  mutate(IUCN = decostand(IUCN, "range", na.rm = TRUE))

# Function to process datasets for each zone
process_zone_metrics <- function(shark.info, puvsp_zone, zone_name) {
  # Get the species names from zone-specific dataset
  species_in_zone <- names(puvsp_zone)[4:ncol(puvsp_zone)]
  
  # Filter shark.info for zone-specific species
  zone_shark_info <- shark.info %>%
    filter(Species %in% species_in_zone) %>%
    mutate(SpeciesID = row_number())
  
  # Print the number of species retained
  print(paste("Number of species retained in", zone_name, ":", nrow(zone_shark_info)))
  
  # Process each metric for the zone
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
    process_dataset(zone_shark_info, puvsp, metric)
  })
  
  names(results) <- metrics
  
  # Save the results
  saveRDS(results, 
          here::here("Data", "My dataframes", 
                     paste0(zone_name, "_shark_conservation_metrics.rds")))
  
  jsonlite::write_json(results, 
                       here::here("Data", "My dataframes", 
                                  paste0(zone_name, "_shark_conservation_metrics.json")))
  
  return(results)
}

# Process for high seas
highseas_results <- process_zone_metrics(shark.info, 
                                         highseas_species, 
                                         "highseas")

# Process for continental waters
continental_results <- process_zone_metrics(shark.info, 
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

# Function to process already transformed PUVSP data
process_zone_data <- function(puvsp_data, coordinates_data, zone_name) {
  # No need to transform the data, just save it
  write.csv(puvsp_data, 
            here::here("Data", "My dataframes", paste0(zone_name, "_puvsp_harmonised.csv")), 
            row.names=FALSE)
  
  write.csv(coordinates_data, 
            here::here("Data", "My dataframes", paste0(zone_name, "_coords_harmonised.csv")), 
            row.names=FALSE)
  
  # Return summary statistics
  return(list(
    n_species = length(unique(puvsp_data$sp)),
    n_planning_units = nrow(coordinates_data),
    species_ids = unique(puvsp_data$sp)
  ))
}

# Process high seas data
highseas_stats <- process_zone_data(
  highseas_puvsp,
  highseas_coordinates,
  "highseas"
)

# Process continental data
continental_stats <- process_zone_data(
  continental_puvsp,
  continental_coordinates,
  "continental"
)

# Print summaries
cat("\nHigh Seas Summary:\n")
cat("Number of species:", highseas_stats$n_species, "\n")
cat("Number of planning units:", highseas_stats$n_planning_units, "\n")

cat("\nContinental Waters Summary:\n")
cat("Number of species:", continental_stats$n_species, "\n")
cat("Number of planning units:", continental_stats$n_planning_units, "\n")

# Verify common species between zones
common_species_both_zones <- intersect(
  highseas_stats$species_ids,
  continental_stats$species_ids
)
cat("\nSpecies common to both zones:", length(common_species_both_zones), "\n")

# Check final datasets
for (zone in c("highseas", "continental")) {
  cat(paste("\nChecking", zone, "datasets:\n"))
  
  puvsp_file <- read.csv(here::here("Data", "My dataframes", paste0(zone, "_puvsp_harmonised.csv")))
  coords_file <- read.csv(here::here("Data", "My dataframes", paste0(zone, "_coords_harmonised.csv")))
  
  cat("PUVSP dimensions:", dim(puvsp_file), "\n")
  cat("Coordinates dimensions:", dim(coords_file), "\n")
  cat("Unique species in PUVSP:", length(unique(puvsp_file$sp)), "\n")
}
