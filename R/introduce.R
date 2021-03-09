#' Seed incursions explicitly
#'
#' \code{sim_incursions_hardwired} takes known incursions and seeds them explicitly
#' in time and space.
#'
#' These functions are passed as an argument to the parameter `incursion_fun` in
#' \code{\link{simrabid}} and can be customized. They must include arguments
#' for `cell_ids` and `params` (if you want to change the parameter
#' values or you can fix them within the function).
#'
#' @param cell_ids a dummy variable passed to keep consistent with other
#'  sim incursion function
#' @param params consisting of two numeric vectors:
#'  the cell ids of the known incursions (`cell_ids`) and
#'  the timestep at which the known incursions should be seeded (`tstep`)
#'
#' @inherit sim_incursions_pois return
#' @export
#' @keywords incursions
#'
sim_incursions_hardwired <- function(cell_ids,
                                     params = list(cell_ids, tstep)) {

  current_t <- get("t", envir = parent.frame(1))

  # filter list of empirical incursions
  cell_id <- params$cell_ids[params$tstep %in% current_t]

  return(cell_id)
}

#' Simulate incursions
#'
#' \code{sim_incursions_pois} simulates incursions from a poisson distribution.
#'
#' These functions are passed as an argument to the parameter `incursion_fun` in
#' \code{\link{simrabid}} and can be customized. They must include arguments
#' for `cell_ids` and `params` (if you want to change the parameter
#' values or you can fix them within the function).
#'
#' @param cell_ids the cells to sample incursions to
#' @param params a list with a named parameter `iota`, the number of incursions
#'   on average in each time step
#'
#' @return numeric vector of cell_ids where incursions occured
#'
#' @export
#' @keywords incursions
#'
sim_incursions_pois <- function(cell_ids,
                                params = list(iota = 1)) {

  # number of incursions this week
  n_incs <- rpois(1, params$iota)

  # sample cell ids by bins
  cell_id <- safe_sample(x = cell_ids, size = n_incs, replace = TRUE)

  return(cell_id)
}


#' Add incursions to infectious line list
#'
#' \code{add_incursions} adds incursions generated from `incursion_fun` function
#' to the line list of infections.
#'
#' @param cell_id_incs numeric vector, number of incursions to add to each grid cell, generated
#'   from `incursion_fun`
#' @param cell_ids integer vector, of all possible cell_ids to match between cell & row_ids
#' @param ncells integer, number of grid cells total in the raster of simulated space
#' @param admin_ids integer vector, admin ids corresponding to the row_ids
#' @param row_ids integer vector, row id of state matrix
#' @param x_coord numeric vector of Eastings which correspond to the centroid of grid cells
#' @param y_coord numeric vector of Northings which correspond to the centroid of grid cells
#' @param tstep integer, the time in which the incursions are being added (i.e. the current time step)
#' @param counter integer, the id number to start from to assign ids to the new incursions
#' @param days_in_step integer, the number of days in each tstep, to assign the day
#'   of infectiousness to each incursion
#'
#' @inherit sim_bites return
#' @export
#' @import data.table
#' @keywords incursions internal
#'
add_incursions <- function(cell_id_incs, cell_ids, ncells,
                           admin_ids = NULL, row_ids,
                           x_coord, y_coord, tstep,
                           counter, days_in_step = 7) {

  n_incs <- length(cell_id_incs)

  # date infectious (this tstep just draw the day!)
  t_infectious <- sample(days_in_step, n_incs,
                          replace = TRUE)/days_in_step + tstep
  row_id <- get_rowid(cell_id_incs, cell_ids, admin_ids, row_ids, ncells)

  x_coord <- x_coord[cell_id_incs]
  y_coord <- y_coord[cell_id_incs]

  # incursions have a progenitor id of -1
  data.table(id = counter + 1:n_incs,
             cell_id = cell_id_incs, row_id, progen_id = -1L,
             path = 0L,
             x_coord, y_coord,
             invalid = FALSE, outbounds = FALSE,
             t_infected = 0, contact = "N",
             infected = TRUE, t_infectious)
}

