get.demdata <- function(shapefile, res_meters = 1000, pop = pop_data,
                        census = census_data) {

  # 1. Get starting pop and births from village censuses
  pop$VILLCODES <- as.factor(pop$VILLCODES)
  shapefile@data %>%
    left_join(pop, by = c("VILLCODE" = "VILLCODES")) %>%
    group_by(Village_2002) %>%
    summarize(pop2002 = sum(Population.2002, na.rm = TRUE),
              pop2012 = sum(Population.2012, na.rm = TRUE),
              villcode = VILLCODE[1]) %>%
    mutate(growth = exp(log(pop2012/pop2002)/10)) %>%
    right_join(shapefile@data, by = c("villcode" = "VILLCODE")) -> shapefile@data

  if (is.numeric(res_meters)) {
    ## 2. Get raster
    r <- raster(shapefile)
    res(r) <- res_meters
    SD_raster <- rasterize(shapefile, r)
    rast_attributes <- SD_raster@data@attributes[[1]]
    shapefile$ID_match <-rast_attributes$ID[match(shapefile$villcode, rast_attributes$villcode)]

    ## 3a. Get proportion of dogs in each cell (from census)
    census_dogs <- rasterize(cbind(census$X, census$Y), r,
                             field = census$dogs + census$pups, fun = sum)
    summarize_dogs <- as.data.frame(cbind(census_dogs@data@values, SD_raster@data@values))
    summarize_dogs %>%
      group_by(village_ID = V2) %>%
      summarize(dogs_total = sum(V1, na.rm = TRUE)) %>%
      right_join(shapefile@data, by = c("village_ID" = "ID_match")) -> shapefile@data
    raster_pop <- rasterize(shapefile, r, field = shapefile$dogs_total)
    prop <- census_dogs/raster_pop

    ## 3b.And also dogs vaccinated
    census_vacc <- rasterize(cbind(census$X, census$Y), r,
                             field = census$vacc_dogs + census$vacc_pups, fun = sum)
    cov <- census_vacc/census_dogs
    ##' add small prob so that allocates equally if vaccinations happen in place where
    ##' no coverage reported in census and can be allocated equally when 'leftovers'
    cov[cov > 1] <- 1 ## max of 1
    cov[is.na(cov) | cov == Inf] <- 0
    cov <- cov + 1e-5

    ## 4. Get HDR and starting dog population from starting vill population
    census_humans <- rasterize(cbind(census$X, census$Y), r,
                               field = census$Adults + census$Children, fun = sum)
    summarize_humans <- as.data.frame(cbind(census_humans@data@values, SD_raster@data@values))
    summarize_humans %>%
      group_by(village_ID = V2) %>%
      summarize(humans_total = sum(V1, na.rm = TRUE)) %>%
      right_join(shapefile@data, by = c("village_ID" = "village_ID")) -> shapefile@data
    shapefile$HDR <- shapefile$humans_total/shapefile$dogs_total
    shapefile$start_pop <- shapefile$pop2002/shapefile$HDR
    HDR <- rasterize(shapefile, r, field = shapefile$HDR)
    growth <- rasterize(shapefile, r, field = shapefile$growth)
    start_pop <- rasterize(shapefile, r, field = shapefile$start_pop) ## just the starting vill pops

    ## 5. Return stacked raster file
    dat <- as.data.table(list(village_ID = SD_raster@data@values,
                              HDR = HDR@data@values, growth = growth@data@values,
                              prop_pop = prop@data@values, start_pop = start_pop@data@values,
                              cov = cov@data@values))
    dat$start_pop <- rbinom(nrow(dat), size = round(dat$start_pop), prob = dat$prop_pop)
    dat$villcode <- rast_attributes$villcode[match(dat$village_ID, rast_attributes$ID)]
    dat$cell_id <- 1:nrow(dat) ## in order to use identity of cell!
    return(dat)
  } else {
    ## For village data get HDR + start_pop
    census %>%
      group_by(Village.2002) %>%
      summarize(dogs_total = sum(dogs + pups, na.rm = TRUE),
                humans_total = sum(Adults + Children, na.rm = TRUE),
                HDR = humans_total/dogs_total) %>%
      right_join(shapefile@data, by = c("Village.2002" = "Village_2002")) %>% # keep all shapefile data
      mutate(start_pop = pop2002/HDR) -> shapefile@data
    return(shapefile@data)
  }
}
