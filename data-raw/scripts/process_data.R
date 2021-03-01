# ------------------------------------------------------------------------------
#' Proccess data
# ------------------------------------------------------------------------------

library(readr)
library(dplyr)
library(sf)
library(raster)
library(fasterize)
library(data.table)
library(janitor)
library(magrittr)
library(lubridate)
select <- dplyr::select
source("R/utils-data.R")

# Data -------------------------------------

# shapefile
sd_shape <-
  st_read("data-raw/raw/SD_shape/SD_From_HHS/SD_Villages_2002_From_HHS_250m_Smoothed_UTM.shp")
sd_shape <- st_transform(sd_shape, crs = CRS("+init=epsg:32736"))

# district wide census (@ vill level)
sd_pops <- read_csv("data-raw/raw/SerengetiPop.csv")
sd_census <- read_csv(get_latest("data-raw/raw/wisemonkey", "Census"))
case_dt <- read_csv(get_latest("data-raw/raw/wisemonkey", "Animal_Contact_Tracing"))
vacc_dt <- read_csv(get_latest("data-raw/raw/wisemonkey", "Vaccination"))

# Get populations by village (growth & hdr) ----
sd_pops %>%
  clean_names() %>%
  group_by(village_2002) %>%
  summarize(pop_2002 = sum(population_2002, na.rm = TRUE),
            pop_2012 = sum(population_2012, na.rm = TRUE),
            villcode = as.character(villcodes[1])) %>%
  mutate(growth = exp(log(pop_2012/pop_2002)/10),
         row_id = row_number()) %>%
  left_join(clean_names(sd_shape), .,
            by = c("vill_2002" = "village_2002")) -> sd_shapefile

usethis::use_data(sd_shapefile, overwrite = TRUE)

# Clean census data and out neccessary bits -----
sd_census %>%
  clean_names() %>%
  select(village_2002, utm_easting, utm_northing,
         adults, children, adult_dogs = dogs_3_months,
         adult_dogs_vacc = vaccinated_dogs_3_months,
         pups = pups_3_months, pups_vacc = vaccinated_pups_3_months,
         date = actvitity_date) %>%
  filter(!is.na(utm_easting), !is.na(utm_northing)) -> sd_census_data

usethis::use_data(sd_census_data, overwrite = TRUE)

# Case data ----
case_dt %>%
  clean_names() %>%
  select(village_2002, utm_easting, utm_northing, biter_id, species, suspect,
         symptoms_started_known, symptoms_started,
         symptoms_started_accuracy)  %>%
  filter(suspect %in% "Yes", symptoms_started_known,
         species %in% "Domestic dog") -> sd_case_data

usethis::use_data(sd_case_data, overwrite = TRUE)

# Vacc data ----
vacc_dt %>%
  clean_names() %>%
  select(village_2002, start_date, end_date,
         doses, dogs_vaccinated) %>%
  # when dogs_vaccinated not specified assume 95% of doses go to dogs
  mutate(dogs_vaccinated = coalesce(dogs_vaccinated,
                                         round(doses * 0.95))) -> sd_vacc_data

correct_names <- tribble(~village_2002, ~corrected,
                          "Mbirikili", "Bonchugu",
                          "Kebanchabancha" , "Kebanchebanche",
                          "Natta" , "Mbisso",
                          "Natta Mbiso" , "Mbisso",
                          "Nyamisingisi" , "Nyamasingisi")

sd_vacc_data %<>%
  left_join(correct_names) %>%
  mutate(village_2002 = coalesce(corrected, village_2002),
        match_2002 = sd_pops$VILLCODES[match(village_2002,
                                              sd_pops$Village_2002)],
        match_2012 = sd_pops$VILLCODES[match(village_2002,
                                               sd_pops$Village_2012)],
        villcode = coalesce(match_2002, match_2012)) %>%
  select(-contains(c("match", "corrected")))

usethis::use_data(sd_vacc_data, overwrite = TRUE)

# Parameterization defaults (per Mancy et al. 2021) -----
disp <- read.csv("data-raw/raw/parameters_KH/DK_params.csv")
steps <- read.csv("data-raw/raw/parameters_KH/steps.distribution.csv")
serial <- read.csv("data-raw/raw/parameters_KH/SI_params.csv")

param_defaults <-
  list(
    steps_shape = steps$shape, steps_scale = steps$scale,
    disp_meanlog = disp$DK_meanlog, disp_sdlog = disp$DK_sdlog,
    serial_meanlog = serial$SI_ml, serial_sdlog = serial$SI_sdlog,
    detect_alpha = 80.1, detect_beta = 8.9,
    detect_prob = 0.9
  )


usethis::use_data(param_defaults, overwrite = TRUE)

