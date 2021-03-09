#' Setup spatial bounds of simulation
#'
#' \code{setup_space} creates the spatial grid of the simulation from a user defined
#' shapefile and resolution. If vaccination data are aggregated at a larger scale,
#' pass the id column (i.e. numeric id of the administrative unit) to assign each grid cell to.
#' Alternatively if vaccination data are aggregated at the grid cell level or vaccination
#' is absent, then pass NULL to `id_col`.
#'
#' To do:
#' - Add examples with sf and sp & different id_col
#' - Tests: resolution too large? if no fasterize, if no id col
#'
#' @param shapefile a shapefile (either an sf or sp object) of the area being simulated
#'   should be projected in UTM coordinates
#' @param resolution numeric, the resolution in meters
#' @param id_col character, the name of the attribute/column to rasterize the data by,
#'   largely for allocating vaccinations
#' @param use_fasterize boolean, if TRUE and class of shapefile is sf will use the
#'   fasterize package (if installed)
#'
#' @importFrom raster rasterize res raster
#' @return raster of simulation bounds rasterized by the id column specify (i.e.
#' each grid cell allocated to a location)
#' @export
#' @keywords setup
#'
setup_space <- function(shapefile, resolution = 1000,
                        id_col,
                        use_fasterize = FALSE) {

  r <- raster(shapefile)
  res(r) <- resolution # res in meters

  if(is.null(id_col)) {

    values(r) <- 1:length(r)
    return(r)

  } else {
    if(use_fasterize && "sf" %in% class(shapefile) && require(fasterize)) {
      rast <- fasterize::fasterize(shapefile, r, field = id_col)
    } else {
      rast <- rasterize(shapefile, r, field = id_col)
    }

    return(rast)
  }

}

#' Set-up simulation
#'
#' \code{setup_sim} takes inputs and generates all needed outputs to be passed
#' to \code{\link{simrabid}}.
#'
#' @param tmax integer, the maximum number of time steps to run the simulation for
#' @param rast raster, from \code{\link{setup_space}} the output raster with each cell allocated
#'   to a location
#' @param death_rate_annual annual death rate, either length 1 or length of the admin units
#' @param birth_rate_annual annual birth rate, either length 1 or length of the admin units
#' @param waning_rate_annual annual waning rate of vaccination immunity,
#'  either length 1 or length of the admin units
#' @param block_fun a function that designates the indices to track/block/consider
#'  out-of-bounds (see \code{\link{block_cells}} for an example)
#' @param params additional params to pass, must include start_pop, a numeric
#'  vector of the starting population in each grid cell of `rast`
#' @param step the number of steps to translate rate to probabilities
#'  (i.e. 52 = weekly, 365 = daily)
#' @param by_admin logical, should mixing be at the admin (or arbitrary) scale?
#'
#' @return a list of objects needed for the simulation
#' @export
#' @import data.table
#' @keywords setup
#'
setup_sim <- function(tmax, rast,
                      death_rate_annual,
                      birth_rate_annual,
                      waning_rate_annual,
                      block_fun = block_cells,
                      params = list(start_pop),
                      step = 52,
                      by_admin = FALSE) {

  # cells to block / not track
  all_inds <- block_fun(rast, params)

  # only tracking populated cells within district (no new pops in new cells)
  cell_ids <- all_inds$track_inds
  coords <- raster::coordinates(rast)

  if(by_admin) {

    # need to group by id and spit back out
    pop_dt <- data.table(admin_id = rast[], start_pop = params$start_pop)
    start_pop <- pop_dt[!is.na(admin_id)][, .(pop = sum(start_pop, na.rm = TRUE)),
                                          by = "admin_id"]$pop
    bins <- max(rast[], na.rm = TRUE)
    loc_ids <- 1:bins

  } else {
    start_pop <- params$start_pop[cell_ids]
    bins <- length(cell_ids)
    loc_ids <- rast[cell_ids]

  }

  if(length(birth_rate_annual) > 1 & !by_admin) {
    if(length(birth_rate_annual) != max(rast[], na.rm = TRUE)) {
      stop("Error, length of birth rates is not equal to 1 or number of admin units")
    }
    birth_rate_annual <- birth_rate_annual[rast[cell_ids]]
  }

  # other dem params
  death_prob <- get_prob(rate = death_rate_annual, step = step) # annual death rate to prob
  birth_prob <- get_prob(rate = birth_rate_annual, step = step) # annual birth rate to prob
  waning_prob <- get_prob(rate = waning_rate_annual, step = step) # annual waning to prob

  # state matrices
  row_ids <- 1:bins
  S_mat <- V_mat <- N_mat <- I_mat <- E_mat <- matrix(0L, nrow = bins, ncol = tmax)


  # Set up I_dt (data table)
  I_dt <- data.table(id = 0L, cell_id = 0L, row_id = 0L,
                     progen_id = 0L, path = 0L, x_coord = 0,
                     y_coord = 0, invalid = TRUE, outbounds = TRUE,
                     t_infected = 0, contact = "N",
                     infected = FALSE, t_infectious = 0)

  I_dt <- I_dt[rep(I_dt[, .I], 1000)]

  empty_dt <-  I_dt[0]

  return(list(row_ids = row_ids, S_mat = S_mat, I_mat = I_mat, E_mat = E_mat,
              V_mat = V_mat, N_mat = N_mat, cell_ids = cell_ids,
              cells_block = all_inds$block_inds,
              cells_out_bounds = all_inds$out_inds,
              I_dt = I_dt, empty_dt = empty_dt,
              loc_ids = loc_ids,
              start_pop = start_pop,
              nrows = nrow(rast),
              ncols = ncol(rast),
              ncells = ncell(rast),
              x_topl = bbox(rast)[1, "min"],
              y_topl = bbox(rast)[2, "max"],
              x_coord = coords[, 1],
              y_coord = coords[, 2],
              admin_ids = rast[], # admin unit ids to aggregate to
              bins = bins,
              res_m = res(rast)[1], tmax = tmax,
              death_prob = death_prob, birth_prob = birth_prob,
              waning_prob = waning_prob))
}

#' Initialize the simulation
#'
#' \code{init} starts up the simulation by drawing the starting state values
#'   and seeds infectious cases
#' @keywords internal
#' @import data.table
#'
init <- function(start_pop, start_vacc, I_seeds, I_dt, cell_ids,
                 admin_ids = NULL, row_ids,
                 bins, x_coord, y_coord,
                 params, incursion_fun, ncells) {

  # Starting pop + sus
  if(length(start_vacc) != 1 & is.integer(start_vacc)) {
    V <- start_vacc
  } else {
    V <- rbinom(n = bins, size = start_pop, prob = start_vacc)
  }

  S <- start_pop - V

  # Seed cases at t0 and create data.table

  if(I_seeds > 0) {
    cell_id_incs <- safe_sample(opts = cell_ids, size = I_seeds, replace = TRUE)

    I_init <- add_incursions(cell_id_incs, cell_ids, ncells,
                             admin_ids, row_ids,
                             x_coord, y_coord, tstep = 1,
                             counter = 0,
                             days_in_step = 7)

    I_dt[I_init$id] <- I_init # this should get updated in global environment
    I <- tabulate(I_dt$row_id, nbins = bins)
  } else {
    I <- rep(0, bins)
  }


  # Summarize and put into I!
  E <- rep(0, bins)

  return(list(S = S, V = V, I = I, E = E, I_dt = I_dt, N = start_pop))
}

#' Doubles the infection line list
#'
#' @keywords internal
#'
double_I <- function(I_dt) {
  I_skeleton <- data.table(id = 0L, cell_id = 0L, row_id = 0L,
                           progen_id = 0L, path = 0L, x_coord = 0,
                           y_coord = 0, invalid = TRUE, outbounds = TRUE,
                           t_infected = 0, contact = "N",
                           infected = FALSE, t_infectious = 0)

  I_dt <- rbind(I_dt, I_skeleton[rep(I_skeleton[, .I], nrow(I_dt))])
  return(I_dt)
}

#
#' Which cells should be considered invalid if movement occurs?
#'
#' This is an example function for desginating certain cells invalid/outofbounds.
#' Can customize this to account for landscape features (i.e. lakes, rivers, mountains)
#' where movement is not possible. Function must have two arguments: the
#' raster input and other params passed as a list. Function must output the
#' indexes to block, the indexes considered out of bounds, and the indices to track.
#'
#' @param rast raster generated using `setup_space` function
#' @param params a list of params, in this case,
#'  integer vector of the population size in each cell of `rast`
#'
#' @return indexes to block, indexes considered out of bounds, and indexes to track
#' @export
#'
block_cells <- function(rast, params = list(start_pop)) {

  in_inds <- which(!is.na(rast[]))
  out_inds <- which(is.na(rast[]))
  no_pop <- which(params$start_pop < 1 | is.na(params$start_pop))
  block_inds <- no_pop[no_pop %in% in_inds]

  # Ones to track: inside district & unblocked
  track_inds <- which(!is.na(rast[] & params$start_pop > 0))

  return(list(block_inds = block_inds, out_inds = out_inds,
              track_inds = track_inds))
}

