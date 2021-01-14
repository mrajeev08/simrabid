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
#' @param start_pop integer vector, the population size in each cell of `rast`
#' @param rast raster, from `setup_space` the output raster with each cell allocated
#'   to a location
#'
#' @return a list of objects needed for the simulation
#' @export
#' @import data.table
#' @keywords setup
#'
setup_sim <- function(tmax, start_pop, rast) {

  inds <- !is.na(rast[]) & !is.na(start_pop)
  cell_ids <- (1:ncell(rast))[inds]
  loc_ids <- rast[inds]
  nlocs <- length(cell_ids)
  coords <- raster::coordinates(rast)[inds, ]
  start_pop <- start_pop[inds]

  # state matrices
  row_ids <- 1:nlocs # don't need to do this function call each time either
  S_mat <- V_mat <- N_mat <- I_mat <- E_mat <- matrix(0L, nrow = nlocs, ncol = tmax)
  rows_pop <- row_ids[start_pop > 0]
  cells_pop <- cell_ids[start_pop > 0]

  # Set up I_dt (data table)
  I_dt <- data.table(id = 0, cell_id = 0, row_id = 0,
                     progen_id = 0, path = 0L, x_coord = 0,
                     y_coord = 0, populated = FALSE, within = FALSE,
                     t_infected = 0, contact = "N",
                     infected = FALSE, t_infectious = 0)

  I_dt <- I_dt[rep(I_dt[, .I], 1000)]

  empty_dt <-  I_dt[0]

  return(list(row_ids = row_ids, S_mat = S_mat, I_mat = I_mat, E_mat = E_mat,
              V_mat = V_mat, N_mat = N_mat, cell_ids = cell_ids,
              cells_pop = cells_pop, rows_pop = rows_pop, I_dt = I_dt,
              empty_dt = empty_dt, nlocs = nlocs, loc_ids = loc_ids,
              start_pop = start_pop,
              nrow = nrow(rast),
              ncol = ncol(rast),
              x_topl = bbox(rast)[1, "min"],
              y_topl = bbox(rast)[2, "max"],
              x_coord = coords[, 1],
              y_coord = coords[, 2], res_m = res(rast)[1], tmax = tmax))
}

#' Initialize the simulation
#'
#' \code{init} starts up the simulation by drawing the starting state values
#'   and seeds infectious cases
#' @keywords internal
#' @import data.table
#'
init <- function(start_pop, start_vacc, I_seeds, I_dt,
                 rows_pop, cell_ids, nlocs,
                 x_coord, y_coord) {

  # Starting pop + sus
  if(length(start_vacc) != 1 & is.integer(start_vacc)) {
    V <- start_vacc
  } else {
    V <- rbinom(n = nlocs, size = start_pop, prob = start_vacc)
  }

  S <- start_pop - V

  # Seed cases at t0 and create data.table
  row_id <- sample(rows_pop, I_seeds, replace = TRUE)

  I_init <- data.table(id = 1:I_seeds, cell_id = cell_ids[row_id],
                      row_id, progen_id = 0, path = 0L,
                      x_coord = x_coord[row_id],
                      y_coord = y_coord[row_id],
                      populated = TRUE, within = TRUE,
                      t_infected = 0, contact = "N",
                      infected = TRUE,
                      t_infectious = 1)

  I_dt[I_init$id] <- I_init # this should get updated in global environment

  # Summarize and put into I!
  I <- tabulate(I_dt$row_id, nbins = nlocs)
  E <- rep(0, nlocs)

  return(list(S = S, V = V, I = I, E = E, I_dt = I_dt, N = start_pop))
}

#' Doubles the infection line list
#'
#' @keywords internal
#'
double_I <- function(I_dt) {
  I_skeleton <- data.table(id = 0, cell_id = 0, row_id = 0,
                           progen_id = 0, path = 0L, x_coord = 0,
                           y_coord = 0, populated = FALSE, within = FALSE,
                           t_infected = 0, contact = "N",
                           infected = FALSE, t_infectious = 0)

  I_dt <- rbind(I_dt, I_skeleton[rep(I_skeleton[, .I], nrow(I_dt))])
  return(I_dt)
}
