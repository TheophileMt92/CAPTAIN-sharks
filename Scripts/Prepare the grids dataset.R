# This script does the following:
# 1. assigns biomes/regions to shark species from IUCN data
# 2. Transform grids into occurrences to then calculaye per-species geographic range
# Date: 16 Jul 2021
# Author: Catalina Pimiento

library(readxl)
library(stringr)
library(tidyverse)
library(dplyr)
library(janitor)
library(reshape2)
library (raster)

#########################################################################################################
# REGIONS
#########################################################################################################
# read look up table with all the shark names and synonyms
names.raw<-read_xlsx("~/Dropbox (Smithsonian)/SHARKS/Taxonomy/Lookup_Taxonomy.xlsx", sheet = "Species")

names.simple<- names.raw %>%
  dplyr::select(-Family, -Order, -Superorder)

# load spatial data
distr <- as_tibble(
  readRDS("Data/Map_PA_biomes_0_5.RDS"))

# load("~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/Mat_Pa_CHONDRICHTHYES.rdata")
# distr2<-as_tibble(Mat_Pa_CHONDRICHTHYES)

# load data frame with fixed synonyms
load('Data/iucn.names_traits.names.RData') # distr.names

# replace names in distr with valid names (there will be duplicates)
new.names <- c("lon","lat", distr.names$new.names)
length(new.names)
which(duplicated(new.names)==TRUE)

distr.no.syn<- distr %>%
  rename(lon = X_ENTIER, lat = Y_ENTIER) %>%
  dplyr:: select(-Region, -Biome, -Biome1) %>%
  setNames(new.names)

# If you want to focus on the last few columns:
num_cols <- 5  # Change this to the number of last columns you want to see
last_columns <- tail(names(distr.no.syn), num_cols)

print("\nHead of the tibble with the last few columns:")
print(head(distr.no.syn[last_columns]))

# Try to get the last few columns, with error handling
num_cols <- 50  # Change this to the number of last columns you want to see
total_cols <- ncol(distr.no.syn)

if (total_cols > 0) {
  num_cols <- min(num_cols, total_cols)
  last_columns <- tail(names(distr.no.syn), num_cols)
  
  print(paste("\nLast", num_cols, "column names:"))
  print(last_columns)
  
  if (length(last_columns) > 0) {
    print("\nHead of the tibble with the last few columns:")
    print(head(distr.no.syn[, last_columns, drop = FALSE]))
  } else {
    print("\nNo columns selected.")
  }
} else {
  print("\nThe tibble appears to be empty.")
}

# merge occurrence values from duplicated names (merge junior synonyms)
distr.no.syn <- 
  as_tibble(do.call(cbind,lapply(split(seq_len(ncol(distr.no.syn)),names(distr.no.syn)),
                                 function(x) rowSums(distr.no.syn[x]))))

which(duplicated(names(distr.no.syn))==TRUE)
