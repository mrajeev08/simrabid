#' Simulate transmission outcomes of exposures at individual level
#'
#' \code{sim_trans} simulates the outcome of each exposure/transmission event
#' (i.e. whether the exposure was to a susceptible)
#'
#' This function takes the locations of each exposure/transmission event, tabulates
#' them, and then allocates them to individuals of each state (using sampling). These
#' outcomes can be tracked explicitly (i.e. track if S/E/I/V or M if no individuals left to allocate to
#' in that location) or just whether with a susceptible (S) or not (M) (see `track` argument in
#' \code{\link{simrabid}}.
#'
#' To Do:
#' - Tests: returns should be same length as length of row_id and should be S/E/I/V/M or S/M (no NAs)
#' - rewrite in Rcpp? Will it actually speed this up? Maybe by reducing cost of function calls
#'
#' @param row_id numeric vector (length >= 1) the row ids corresponding to the grid cell of the exposure (one value per exposure)
#' @param S,E,I,V numeric vector of state variables (i.e. # of individuals in
#'   each class in each grid cell) of length nlocs
#' @param bins numeric, either the number of grid cells total OR
#'  the number of admin units, for when simulating at admin rather
#'  than grid cell level
#' @param track boolean, whether to explicitly track the outcome of each exposure or
#'   to only track whether successful or not (i.e. with a suscpetible or no). See details.
#'
#' @return a list of two vectors (contact and infected) of same length as `row_id`.
#'  contact will be NULL if track = FALSE.
#' @keywords transmit internal
#'
sim_trans <- function(row_id, S, E, I, V, bins, track = FALSE) {

  exps <- tabulate(row_id, nbins = bins)
  id_list <- which(exps > 0)
  nexps <- length(row_id)

  # track who exposure was allocated to
  contact <- rep("M", nexps)
  infected <- rep(FALSE, nexps)

  if(!track) {
    # probability that a exposure was with a susceptible (can't have exp with yourself)
    N <- S + E + I + V
    probs <- ifelse(N > 0, S/N, 0)

    exps_out <- rbinom(n = length(exps), size = exps,
                       prob = probs)

    # Number of successful exposures can't be greater than the susceptibles available
    # Issue when N is small
    success <- ifelse(exps_out[id_list] > S[id_list], S[id_list],
                      exps_out[id_list])

    for(i in seq_along(id_list)) {
      if(success[i] > 0) {
        inf_inds <- which(row_id == id_list[i])
        if(length(inf_inds) > success[i]) {
          infected[sample(inf_inds,
                          success[i],
                          replace = FALSE)] <- TRUE
        } else {
          infected[inf_inds] <- TRUE
        }

      }

    }

  } else {

    for(i in seq_along(id_list)) {
      # Sample the states
      ind <- id_list[i]
      # taking out the infectious dog doing the biting (I[i] - 1)
      states <- rep(c("S", "E", "I", "V"),
                    c(S[ind], E[ind],
                      ifelse(I[ind] > 0, I[ind] - 1, I[ind]), V[ind]))

      sample_length <- length(states)

      if(sample_length > 0) {
        rows <- which(row_id == ind)
        sample_length <- ifelse(exps[ind] > sample_length,
                                sample_length, exps[ind])
        rows <- rows[1:sample_length]
        contact[rows] <- sample(states, size = sample_length,
                                replace = FALSE)
        infected[contact == "S" & !is.na(contact)] <- TRUE
      }
    }
  }

  return(list(contact = contact, infected = infected))
}

#' Simulates biting/exposure behavior of currently infected individuals
#'
#' \code{sim_bites} simulates spatially-explicit transmission events per
#' infectious individual.
#'
#' This function takes the currently infectious individuals and the draws of the number of
#' secondary cases that they seed and simulates their movement.
#' Movements can be sequential (i.e. infected individuals moves from origin to
#' first location of exposure, then from this to the second, etc.)
#' or can be kernel based where movements are all seeded from
#' the origin location of the infected individual. Movements outside of the bounds
#' of the simulation or to uninhabitable grid cells can either be accepted (if
#' `leave_bounds` or `allow_invalid` are `TRUE`). If `FALSE`, these movements will not be
#' accepted and will be redrawn `max_tries` times.
#'
#' To Do:
#' - Tests: returns should have same column structure & order as the I_dt template
#'   and rows should = sum of secondaries
#' - rewrite in Rcpp? Will it actually speed this up? Maybe by reducing cost of function calls
#'
#' @param secondaries numeric vector of number of secondary cases for each infectious individual
#' @param ids numeric vector of the ids of infected individuals (progenitor ids)
#' @param x_coords numeric vector, x coordinate (Easting) of infected individuals
#' @param y_coords numeric vector, y coordinate (Northing) of infected individuals
#' @param t_infectious numeric vector, the time (fractional time step) at which each
#'   individual became infectious
#' @param counter integer, the id number to start with for resulting secondary cases
#' @param dispersal_fun function with a single parameter, n, which will draw distances
#'   for each individual to move in meters
#' @param row_ids integer vector, the row ids corresponding to each grid cell in which
#'   the infectious individual is starting from
#' @param cell_ids integer vector, the cell ids corresponding to each grid cell in which
#'   the infectious individual is starting from
#' @param sim_movement the movement function you want to use for simulating movement, options are sim_movement_continuous and sim_movement_prob (pass this through \code{\link{simrabid}})
#' @inheritParams accept
#' @inheritParams sim_movement_continuous
#' @param weights numeric vector, weights to give each grid cell for sampling moves
#'  for use with sim_movement_prob; will be length ncells + 1,
#'  see \code{\link{cell_weights}}) for more details.
#' defaults to NULL
#' @param admin_ids if aggregating by the administrative unit (rather than grid cell), pass the rasterized
#' administrative ids for each cell
#' @param max_tries integer, the maximum number of tries to make before accepting
#'   an invalid movement (i.e. transmission event fails due to either leaving the
#'   bounds of the simulation or moving to an uninhabitable grid cell) if either
#'   or both leave_bounds and allow_invalid are FALSE.
#'
#' @import data.table
#' @return a data.table that corresponds to the columns in I_dt in the simulation.
#' See \code{\link{simrabid}} for full description.
#' @keywords move internal
#'
sim_bites <- function(secondaries, ids = I_now$id,
                      x_coords = I_now$x_coord,
                      y_coords = I_now$y_coord,
                      t_infectious = I_now$t_infectious,
                      counter = max(I_dt$id),
                      sim_movement = sim_movement_continuous,
                      dispersal_fun, res_m,
                      row_ids, cell_ids, cells_block, cells_out_bounds,
                      nrows, ncols, ncells,
                      x_topl, y_topl,
                      weights = NULL, admin_ids = NULL,
                      sequential = TRUE, allow_invalid = TRUE,
                      leave_bounds = TRUE, max_tries = 100,
                      params, ...) {

  # set up
  progen_ids <- rep(ids, secondaries)
  origin_x <- rep(x_coords, secondaries)
  origin_y <- rep(y_coords, secondaries)
  nsim <- length(progen_ids)
  dist_m <- dispersal_fun(nsim, params)
  if(is.null(weights))  angles <- angle_fun(nsim) else angles <- NULL

  if(sequential) {

    # first movement index of each progenitor
    inds <- match(ids, progen_ids) # returns first match of id in progen id
    inds <- inds[!is.na(inds)]
    out <- vector("list", nsim)

    for (i in seq_along(progen_ids)) {
      if (i %in% inds) { # need progenitor coords for 1st movement
        x <- origin_x[i]
        y <- origin_y[i]
        path <- 1
      } else {
        x <- out[[i - 1]]$x_coord
        y <- out[[i - 1]]$y_coord
        path <- path + 1
      }

      out[[i]] <- sim_movement(dist_m = dist_m[i],
                               angle = angles[i],
                               dispersal_fun,
                               x0 = x, y0 = y, x_topl,
                               y_topl, res_m, ncols,
                               nrows, ncells,
                               cells_block, cells_out_bounds, path,
                               leave_bounds, allow_invalid, max_tries,
                               sequential, weights, params)

      }

    out <- rbindlist(out)

    # getting coords of grid centroid here vectorized!

  } else {

  out <- sim_movement(dist_m = dist_m,
                      angle = angles,
                      dispersal_fun,
                      x0 = origin_x, y0 = origin_y, x_topl,
                      y_topl, res_m, ncols, nrows, ncells,
                      cells_block, cells_out_bounds,
                      path = 0,
                      leave_bounds, allow_invalid, max_tries,
                      sequential, weights, params)
  }

  # add in ids
  out$id <- counter + 1:nsim
  out$progen_id <- progen_ids
  out$t_infected <- rep(t_infectious, secondaries) # tstep became infected
  out$contact <- "M"
  out$infected <- FALSE

  # row ids to match to I/E mats (either corresponding to admin unit or cells)
  out$row_id <- get_rowid(cell_id = out$cell_id, cell_ids, admin_ids, row_ids, ncells)

  # Make sure colum order matches that of I_dt (better way to do this)
  setcolorder(out, c('id', 'cell_id', 'row_id', 'progen_id',
                     'path', 'x_coord', 'y_coord', 'invalid',
                         'outbounds', 't_infected'))
  return(out)

}

# want to do this vectorized @ the end
cell_to_admin <- function(cell_id, admin_ids, ncells) {

  cell_id <- ifelse(cell_id %in% 1:ncells, cell_id, NA)

  if(all(is.na(cell_id))) {
    admin_id <- rep(NA, length(cell_id))
  } else {
    admin_id <- admin_ids[cell_id]
  }

  return(admin_id)

}

# getting cellid for use in movement funs (prob/continuous)
get_cellid <- function(x_coord, y_coord, res_m, x_topl, y_topl,
                       ncols, nrows, ncells) {

  col <- ceiling((x_coord - x_topl) / res_m)
  row <- ceiling(-(y_coord - y_topl) / res_m)
  cell_id <- row * ncols - (ncols - col)
  cell_id <- ifelse(cell_id %in% 1:ncells, cell_id, ncells + 1) # means outside district
  return(cell_id)

}

# Getting row ids
get_rowid <- function(cell_id, cell_ids, admin_ids, row_ids, ncells) {

  if(!is.null(admin_ids)) {

    ids <- cell_to_admin(cell_id, admin_ids, ncells)

    if(all(is.na(ids))) {
      row_id <- rep(NA, length(ids))
    } else {
      row_id <- row_ids[ids] # rows indexed by the admin unit
    }
  } else {

    if(all(is.na(cell_id))) {
      row_id <- rep(NA, length(cell_id))
    } else {
      row_id <- row_ids[match(cell_id, cell_ids)]
    }
  }
  return(row_id)
}


