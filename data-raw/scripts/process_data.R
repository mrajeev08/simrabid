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
source("R/utils-wm.R")

# get latest fun
get_latest <- function(path, pattern) {
  list.files(path, full.names = TRUE)[grep(pattern, list.files(path))][1]
}

# Data -------------------------------------
sd_shape <- st_read("data-raw/raw/SD_shape/Serengeti_villages UTM_region.shp") # shapefile
sd_pops <- read_csv("data-raw/raw/SerengetiPop.csv") # district wide census (@ vill level)
sd_census <- read_csv(get_latest("data-raw/raw/wisemonkey", "Census"))
case_dt <- read_csv(get_latest("data-raw/raw/wisemonkey", "Animal_Contact_Tracing"))
vacc_dt <- read_csv(get_latest("data-raw/raw/wisemonkey", "Vaccination"))

# Get populations by village (growth & hdr) ----
sd_pops %>%
  clean_names(replace = c("Village_2002" = "pop_names")) %>%
  group_by(pop_names) %>%
  summarize(pop_2002 = sum(population_2002, na.rm = TRUE),
            pop_2012 = sum(population_2012, na.rm = TRUE),
            villcode = as.character(villcodes[1])) %>%
  mutate(growth = exp(log(pop_2012/pop_2002)/10),
         row_id = row_number()) %>%
  left_join(clean_names(sd_shape), .) -> sd_shape_pop

# Clean census data and out neccessary bits -----
sd_census %>%
  clean_names() %>%
  select(village_2002, utm_easting, utm_northing,
         adults, children, adult_dogs = dogs_3_months,
         adult_dogs_vacc = vaccinated_dogs_3_months,
         pups = pups_3_months, pups_vacc = vaccinated_pups_3_months) -> sd_census_clean

# Use fasterize to get village ids
r <- raster(sd_shape_pop)
res(r) <- 500 # res in meters (lets do it at 500 to get best idea of locs)
sd_rast <- fasterize(sd_shape_pop, r, field = "row_id")

# Make sure there are no NAs!
sd_census_clean %<>%
  filter(!is.na(utm_easting), !is.na(utm_northing)) %>%
  mutate(villcode = get_ids(x_coord = utm_easting, y_coord = utm_northing,
                             res_m = 500,
                             x_topl = bbox(sd_rast)[1, "min"],
                             y_topl = bbox(sd_rast)[2, "max"],
                             ncol = ncol(sd_rast),
                             nrow = nrow(sd_rast),
                             rast = sd_rast,
                             id_col = sd_shape_pop$villcode))
# ggplot(sd_shape_pop) +
#   geom_sf() +
#   geom_point(data = sd_census_clean,
#              aes(x = utm_easting, y = utm_northing, color = villcode),
#              alpha = 0.5) +
#   guides(color = "none")

# Case data ----
case_dt %>%
  clean_names() %>%
  select(village_2002, utm_easting, utm_northing, biter_id, species, suspect,
         symptoms_started_known, symptoms_started,
         symptoms_started_accuracy)  %>%
  filter(suspect %in% "Yes", symptoms_started_known,
         species %in% "Domestic dog") -> case_data_cleaned

# Vacc data ----
vacc_dt %>%
  clean_names() %>%
  select(village_2002, start_date, end_date,
         doses, dogs_vaccinated) %>%
  # when dogs_vaccinated not specified assume 95% of doses go to dogs
  mutate(dogs_vaccinated = coalesce(dogs_vaccinated,
                                         round(doses * 0.95))) -> vacc_dt_cleaned

correct_names <- tribble(~village_2002, ~corrected,
                          "Mbirikili", "Bonchugu",
                          "Kebanchabancha" , "Kebanchebanche",
                          "Natta" , "Mbisso",
                          "Natta Mbiso" , "Mbisso",
                          "Nyamisingisi" , "Nyamasingisi")

vacc_dt_cleaned %<>%
  left_join(correct_names) %>%
  mutate(village_2002 = coalesce(corrected, village_2002),
        match_2002 = sd_pops$VILLCODES[match(village_2002,
                                              sd_pops$Village_2002)],
         match_2012 = sd_pops$VILLCODES[match(village_2002,
                                               sd_pops$Village_2012)],
         villcode = coalesce(match_2002, match_2012)) %>%
  select(-contains(c("match", "corrected")))

