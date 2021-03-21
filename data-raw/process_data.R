# ------------------------------------------------------------------------------
#' Proccess parameter defaults
# ------------------------------------------------------------------------------

library(readr)
library(dplyr)
library(sf)
library(janitor)
library(here)

# Data -------------------------------------
# shapefile
sd_shape <-
  st_read(here("data-raw/SD_shape/SD_From_HHS/SD_Villages_2002_From_HHS_250m_Smoothed_UTM.shp"))
sd_shape <- st_transform(sd_shape, sf::st_crs(32736))

# district wide census (@ vill level)
sd_pops <- read_csv(here("data-raw/SerengetiPop.csv"))

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

st_write(sd_shapefile, "inst/extdata/sd_shapefile.shp", append = FALSE)

# Parameterization defaults (per Mancy et al. 2021) -----
disp <- read.csv(here("data-raw/parameters/DK_params.csv"))
steps <- read.csv(here("data-raw/parameters/steps.distribution.csv"))
serial <- read.csv(here("data-raw/parameters/SI_params.csv"))

param_defaults <-
  list(
    steps_shape = steps$shape, steps_scale = steps$scale,
    disp_meanlog = disp$DK_meanlog, disp_sdlog = disp$DK_sdlog,
    serial_meanlog = serial$SI_ml, serial_sdlog = serial$SI_sdlog,
    detect_alpha = 80.1, detect_beta = 8.9,
    detect_prob = 0.9,
    max_secondaries = 100
  )

usethis::use_data(param_defaults, overwrite = TRUE)

