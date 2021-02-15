#' Simulate incursions
#'
#' \code{sim_incursions_pois} simulates incursions from a poisson distribution.
#'
#' These functions are passed as an argument to the paramter `incursion_fun` in
#' \code{\link{simrabid}} and can be customized. They must include arguments
#' for `nlocs` and `params` (if you want to change the parameter
#' values or you can fix them within the function).
#'
#' @param nlocs the total number of grid cells being simulated
#' @param params a list with a named parameter `iota`, the number of incursions
#'   on average each week
#'
#' @return numeric vector of length `nlocs`, the number of incursions to add to
#' each grid cell being simulated
#' @export
#' @keywords incursions
#'
sim_incursions_pois <- function(nlocs,
                                params = list(iota = 1)) {

  # number of incursions this week
  n_incs <- rpois(1, params$iota)

  # sample cell ids by nlocs
  row_id <- sample(nlocs, n_incs, replace = TRUE)
  incs <- tabulate(row_id, nbins = nlocs)
  return(incs)
}

#' Seed incursions explicitly
#'
#' \code{sim_incursions_hardwired} takes known incursions and seeds them explicitly
#' in time and space.
#'
#' @inheritSection sim_incursions_pois Section details
#'
#' @inherit sim_incursions_pois
#' @param params consisting of two numeric vectors:
#'  the row ids of the known incursions (`row_ids`) and
#'  the timestep at which the known incursions should be seeded (`tsteps`)
#' @param current_tstep the current timestep in the simulation
#'
#' @inheritSection sim_incursions_pois Section return
#' @export
#' @keywords incursions
#'
sim_incursions_hardwired <- function(nlocs,
                                     params = list(row_ids = row_ids_empirical,
                                                   tsteps = tstep_empirical),
                                     current_tstep = t) {

  # filter list of empirical incursions
  row_id <- params$row_ids[params$tsteps %in% current_tstep]
  incs <- tabulate(row_id, nbins = nlocs)
  return(incs)
}

#' Add incursions to infectious line list
#'
#' \code{add_incursions} adds incursions generated from `incursion_fun` function
#' to the line list of infections.
#'
#' @param incs numeric vector, number of incursions to add to each grid cell, generated
#'   from `incursion_fun`
#' @param cell_ids integer vector, the cell_ids which correspond to the row_ids (of same length as `incs`)
#' @param x_coord numeric vector of Eastings which correspond to the centroid of grid cells
#' @param y_coord numeric vector of Northings which correspond to the centroid of grid cells
#' @param tstep integer, the time in which the incursions are being added (i.e. the current time step)
#' @param counter integer, the id number to start from to assign ids to the new incursions
#' @param days_in_step integer, the number of days in each tstep, to assign the day
#'   of infectiousness to each incursion
#'
#' @inheritSection sim_bites Section return
#' @export
#' @import data.table
#' @keywords incursions internal
#'
add_incursions <- function(incs, cell_ids, x_coord, y_coord, tstep,
                           counter, days_in_step = 7) {

  n_incs <- sum(incs)

  # date infectious (this tstep just draw the day!)
  t_infectious <- (sample(1:days_in_step, n_incs,
                          replace = FALSE) + tstep)/days_in_step
  row_id <- rep(which(incs > 0), incs[incs > 0])
  cell_id <- rep(cell_ids[incs > 0], incs[incs > 0])

  # incursions have a progenitor id of -1
  data.table(id = counter + 1:n_incs,
             cell_id, row_id, progen_id = -1L,
             path = 0L,
             x_coord = x_coord[row_id],
             y_coord = y_coord[row_id],
             invalid = FALSE, outbounds = FALSE,
             t_infected = 0, contact = "N",
             infected = TRUE, t_infectious)
}

