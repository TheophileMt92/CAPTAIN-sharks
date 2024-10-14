load(here::here("Data", "Raw", "CAPTAIN_Occurance_Grid.Rdata"))

head(grids[,1:10])

Shark_Info=readRDS(here::here("Data", "Raw", "CAPTAIN_SharkInformation.rds"))
Fishing=readRDS(here::here("Data", "Raw", "CAPTAIN_FishingEffortFabien.rds"))
head(Fishing)

library(readr)
shark_pu <- read_csv(here::here("Data", "Shark_pu.csv"))
head(shark_pu)

shark_pu <- read_csv(here::here("Data", "Shark_puvsp_REST.csv"))
head(shark_pu)

shark_pu <- read_csv(here::here("Data", "Shark_Planning_Units.csv"))
head(shark_pu)

shark_pu_MPA <- read_csv(here::here("Data", "Shark_pu_withMPA.csv"))
head(shark_pu_MPA)

shark_puvsp <- read_csv(here::here("Data", "Shark_puvsp_ALL.csv"))
head(shark_puvsp)

EDGE2_Info <- read_csv(here::here("Data", "EDGE2_Info.csv"))
head(EDGE2_Info)

FUSE_Info <- read_csv(here::here("Data", "FUSE_Info.csv"))
head(FUSE_Info)

Shark_Info_MPA <- read_csv(here::here("Data", "CAPTAIN_SharkInformationWithConservation.csv"))
head(Shark_Info_MPA)
