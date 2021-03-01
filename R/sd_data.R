#' Title
#'
#' @param sd_shapefile
#' @param res_m
#' @param sd_census_data
#' @param death_rate_annual
#'
#' @return
#' @export
#' @import raster dplyr fasterize
#' @importFrom magrittr %<>%
#'
get_sd_pops <- function(sd_shapefile, res_m, sd_census_data, death_rate_annual) {


  sd_shapefile$id_col <- 1:nrow(sd_shapefile)
  rast <-
    simrabid::setup_space(shapefile = sd_shapefile,
                          resolution = res_m, id_col = "id_col",
                          use_fasterize = TRUE)

  # Match each census data point to the raster cell id
  # within the get_ids fun, it will also match NA cells
  # to the closest non-NA
  sd_census_data %<>%
    mutate(cell_id = cellFromXY(rast, xy = cbind(utm_easting,
                                                 utm_northing)),
           villcode = get_ids(x_coord = utm_easting,
                              y_coord = utm_northing,
                              rast = rast,
                              id_col = sd_shapefile$villcode))

  # update rast with admin ids
  sd_census_data %>%
    distinct(cell_id, villcode) %>%
    mutate(villcode = match(villcode, sd_shapefile$villcode),
           villcode_check = rast[match(cell_id, seq_len(ncell(rast)))]) %>%
    filter(is.na(villcode_check)) -> update

  rast[update$cell_id] <- update$villcode

  # Group by villcode & get human to dog ratios by village
  sd_census_data %>%
    group_by(villcode) %>%
    tidyr::replace_na(list(adults = 0, children = 0,
                           adult_dogs = 0, pups = 0)) %>%
    summarize(human_pop = sum(adults + children, na.rm = TRUE),
              dog_pop = sum(adult_dogs + pups, na.rm = TRUE)) %>%
    mutate(hdr = human_pop/dog_pop) -> sd_hdr


  # Get starting pop + death + birthrates
  sd_hdr %>%
    left_join(sd_shapefile) %>%
    mutate(start_pop = round(pop_2002/hdr),
           birth_rate = death_rate_annual + growth - 1) -> sd_dem

  sd_census_data %>%
    left_join(select(sd_dem, hdr, start_pop, villcode)) %>%
    group_by(cell_id) %>%
    summarize(dogs_total = sum(adult_dogs + pups, na.rm = TRUE),
              start_pop = start_pop[1], villcode = villcode[1]) %>%
    group_by(villcode) %>%
    summarize(cell_id = safe_sample(opts = cell_id,
                                    size = start_pop[1],
                                    prob = dogs_total,
                                    replace = TRUE)) %>%
    group_by(cell_id) %>%
    summarize(start_pop = n()) %>%
    tidyr::complete(cell_id = seq_len(ncell(rast))) %>%
    arrange(cell_id) -> start_pop

  # start pop (in right order)
  pop <- rast
  pop[start_pop$cell_id] <- start_pop$start_pop


  # return updated pop & admin here & also birth rate here optionally
  return(list(start_pop = pop[], birth_rate_annual = sd_dem$birth_rate, rast = rast,
              death_rate_annual = death_rate_annual))

}

#' Title
#'
#' @param sd_vacc_data
#' @param sd_shapefile
#' @param start_date
#' @param date_format
#' @param days_in_step
#' @param rollup
#'
#' @return
#' @export
#' @import dplyr
#'
get_sd_vacc <- function(sd_vacc_data,
                        sd_shapefile,
                        origin_date = "01-Jan-2002",
                        date_fun = lubridate::dmy,
                        units = "weeks",
                        rollup = 4) {

  sd_vacc_data %>%
    mutate(vacc_times = floor(get_timestep(start_date, origin_date,
                                           date_fun, units))) %>%
    group_by(villcode, vacc_times) %>%
    summarize(vacc_est = sum(dogs_vaccinated)) %>%
    mutate(vacc_locs = match(villcode, sd_shapefile$villcode)) -> vacc_data

  if(!is.null(rollup)) {
    vacc_data %>%
      arrange(vacc_locs, vacc_times) %>%
      mutate(window = vacc_times - dplyr::lag(vacc_times, 1),
             group = if_else(window <= rollup & !is.na(window), 0, 1),
             group = cumsum(group)) %>%
      group_by(vacc_locs, group) %>%
      mutate(vacc_times = min(vacc_times)) %>%
      group_by(vacc_locs, vacc_times) %>%
      summarize(vacc_est = sum(vacc_est)) -> vacc_data

  }

  return(as.data.table(vacc_data[, c("vacc_times", "vacc_est", "vacc_locs")]))
}


