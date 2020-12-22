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
#' @param nlocs numeric, number of grid cells total
#' @param track boolean, whether to explicitly track the outcome of each exposure or
#'   to only track whether successful or not (i.e. with a suscpetible or no). See details.
#'
#' @return a list of two vectors (contact and infected) of same length as `row_id`.
#'  contact will be NULL if track = FALSE.
#' @keywords transmit internal
#'
sim_trans <- function(row_id, S, E, I, V, nlocs, track = FALSE) {

  exps <- tabulate(row_id, nbins = nlocs)
  id_list <- which(exps > 0)
  nexps <- length(row_id)

  # track who exposure was allocated to
  contact <- rep("M", nexps)
  infected <- rep(FALSE, nexps)

  if(!track) {
    # probability that a exposure was with a susceptible
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
#' `leave_bounds` or `allow_empties` are `TRUE`). If `FALSE`, these movements will not be
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
#'   for each individual to move
#' @param row_ids integer vector, the row ids corresponding to each grid cell in which
#'   the infectious individual is starting from
#' @param cell_ids integer vector, the cell ids corresponding to each grid cell in which
#'   the infectious individual is starting from
#' @inheritParams sim_movement
#' @param sequential boolean, if TRUE then movements are sequential, if FALSE, then
#'   movements are kernel based (see Details)
#' @inheritParams accept
#' @param max_tries integer, the maximum number of tries to make before accepting
#'   an invalid movement (i.e. transmission event fails due to either leaving the
#'   bounds of the simulation or moving to an uninhabitable grid cell) if either
#'   or both leave_bounds and allow_empties are FALSE.
#'
#' @import data.table
#' @return a data.table that corresponds to the columns in I_dt in the simulation.
#' See \code{\link{simrabid}} for full description.
#' @keywords transmit internal
#'
sim_bites <- function(secondaries, ids = I_now$id,
                      x_coords = I_now$x_coord, y_coords = I_now$y_coord,
                      t_infectious = I_now$t_infectious,
                      counter = max(I_dt$id),
                      dispersal_fun, res_m,
                      row_ids, cell_ids, cells_pop, nrow, ncol,
                      x_topl, y_topl,
                      sequential = TRUE, allow_empties = TRUE,
                      leave_bounds = TRUE, max_tries = 100) {

  # just use rep(.I situation to replicate?) (see how it does vs. multiple reps)
  progen_ids <- rep(ids, secondaries)
  origin_x <- rep(x_coords, secondaries)
  origin_y <- rep(y_coords, secondaries)

  if(sequential) {

    # first movement index of each progenitor
    inds <- match(ids, progen_ids) # returns first match of id in progen id
    inds <- inds[!is.na(inds)]
    out <- vector("list", length(progen_ids))

    for (i in seq_along(progen_ids)) {
      if (i %in% inds) { # need progenitor coords for 1st movement
        x <- origin_x[i]
        y <- origin_y[i]
        path <- 1
      } else {
        x <- move$x_coord
        y <- move$y_coord
        path <- path + 1
      }

      accept <- 0
      tries <- 0
      while (!accept & tries < max_tries) {

        move <- sim_movement(angle = runif(n = 1, min = 0, max = 2*pi),
                             distance = dispersal_fun(1),
                             x0 = x, y0 = y, x_topl,
                             y_topl, res_m, ncol,
                             nrow, cells_pop, path)
        accept <- accept(leave_bounds, allow_empties, move$within,
                         move$populated)

        tries <- tries + 1

      }
      out[[i]] <- move
    }

    out <- rbindlist(out)

  } else {

    nmoves <- sum(secondaries)

    # get origin from infector if kernel movements
    out <- as.data.table(
      sim_movement(
        angle = runif(n = nmoves, min = 0, max = 2*pi),
        distance = dispersal_fun(nmoves),
        x0 = origin_x, y0 = origin_y, x_topl,
        y_topl, res_m, ncol,
        nrow, cells_pop, path = 0)
      ) # not sequential

    accept <- accept(leave_bounds, allow_empties, within = out$within,
                     populated = out$populated)

    tries <- 0
    while(sum(!accept) > 0 & tries < max_tries) {

      out[!accept, ] <-
        as.data.table(
          sim_movement(angle = runif(n = sum(!accept), min = 0, max = 2*pi),
                     distance = dispersal_fun(sum(!accept)),
                     x0 = origin_x[!accept], y0 = origin_y[!accept],
                     x_topl, y_topl, res_m, ncol, nrow, cells_pop,
                     path = 0L))

      accept <- accept(leave_bounds, allow_empties, within = out$within,
                       populated = out$populated)

      tries <- tries + 1
    }
  }
  # add in ids
  out$id <- counter + 1:length(progen_ids)
  out$progen_id <- progen_ids
  out$t_infected <- rep(t_infectious, secondaries) # tstep became infected
  out$row_id <- row_ids[match(out$cell_id, cell_ids)] # row ids to match to I/E mats
  out$contact <- "M"
  out$infected <- FALSE
  return(out)
}

#' Simulate individual movements
#'
#' Simulates movemement of individuals in continuous space.
#'
#' This simulates movement and uses the top-left coordinates and 1-based indexing
#' of raster cell ids to identify the grid cell moved to by an individual.
#'
#' @param angle numeric vector [0, 360] or value, angle at which to move at
#' @param distance numeric vector or value, distance to move
#' @param x0 numeric vector or value, the x origin of the infected individual
#' @param y0 numeric vector or value, the y origin of the infected individual
#' @param x_topl numeric, the top left x coordinate of the grid on which
#'   movement is being simulated
#' @param y_topl numeric, the top left y coordinate of the grid on which
#'   movement is being simulated
#' @param res_m numeric, the grid cell resolution in meters
#' @param ncol numeric, the number of columns in the grid
#' @param nrow numeric, the number of rows in the grid
#' @param cells_pop integer vector, the cell ids of grid cells that are populated
#'   or generally where movements are valid to
#' @param path integer vector or value, if sequential is FALSE then 0, if TRUE then
#'  the path id (i.e. which step of the movements) to pass through and assign to exposed
#'
#' @return a list with the resulting x and y coordinates, cell ids, whether
#'   the movement was to a populated cell, whether the movement fell within the bounds
#'   of the simulation, and the path id for the exposed.
#' @keywords transmit internal
#'
sim_movement <- function(angle, distance, x0, y0, x_topl,
                         y_topl, res_m, ncol, nrow, cells_pop,
                         path) {

  x_coord <- (sin(angle) * distance * 1000) + x0 # convert to m
  y_coord <- (cos(angle) * distance * 1000) + y0

  # This means can go anywhere
  col <- ceiling((x_coord - x_topl)/res_m)
  row <- ceiling(-(y_coord - y_topl)/res_m)
  cell_id <- row*ncol - (ncol - col)
  populated <- cell_id %in% cells_pop
  within <- row > 0 & row <= nrow & col > 0 & col <= ncol

  # return(data.table(x_coord, y_coord, cell_id, populated, within, path))

  return(list(x_coord = x_coord,
              y_coord = y_coord, cell_id = cell_id,
              populated = populated, within = within, path = path))

}

#' Accept a simulated movement
#'
#' Helper function to determine if a simulated movement is valid.
#' To do: tests (same length as within and boolean no NAs)
#'
#' @param leave_bounds boolean, are movements to outside of the boundaries of the are being
#'   simulated valid
#' @param allow_empties boolean, are movements to empty patches (i.e. with no dogs) valid
#' @param within boolean vector, whether movement falls within bounds
#' @param populated boolean vector, whether movement falls within a populated area
#'
#' @return a boolean vector of length `within`/`populated` corresponding to whether
#'   a movement is accepted as valid
#'
#' @keywords transmit internal
#'
accept <- function(leave_bounds, allow_empties,
                   within, populated) {

  accept <- rep(1, length(within))

  if(!leave_bounds) {
    accept <- ifelse(within, 1, 0)
  }

  if(!allow_empties) {
    accept <- ifelse(populated, 1, 0)
  }

  return(accept)
}

