# These functions should handle zero situations gracefully so you don't have to
# write the conditionals in the IBM part


#' Simulate whether a transmission event was successful (i.e. was with a susceptible)
#'
#' @param cell_id
#' @param S
#' @param N
#' @param nlocs
#' @param track
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' cell_id <- sample(1:1000, 100, replace = FALSE)
#' nlocs <- 1000
#' S <- sample(1000:1500, nlocs, replace = TRUE)
#' N <- sample(2000:3000, nlocs, replace = TRUE)
#' E <- sample(300:400, nlocs, replace = TRUE)
#' I <- sample(100:200, nlocs, replace = TRUE)
#' V <- sample(500:1000, nlocs, replace = TRUE)
#' sim_trans(cell_id = cell_id, nlocs = nlocs,
#'           S = S, N = N, track = FALSE)
#' sim_trans(cell_id = cell_id, nlocs = nlocs,
#'           S = S, N = N, E, I, V, track = TRUE)
# Notes :
# only run this if length(cell_id > 0)
# tests if everything is zero!
# Tests = if exps is zero
# Edge case = if prob S/N is zero or S/E/I/V all zero
# Should be handled before this step?
# use data.table for line list of cases?
sim_trans <- function(exposed, S, N, nlocs, track = FALSE, ...) {

  exps <- tabulate(exposed$row_id, nbins = nlocs)
  id_list <- unique(exposed$row_id)
  nexps <- nrow(exposed)

  if(track == FALSE) {
    # probability that a exposure was with a susceptible
    exps_out <- mapply(rbinom, n = 1, size = exps, prob = S/N)
    outcome <- rep(0, nexps)
    contact <- NULL

    # Number of successful exposures can't be greater than the susceptibles available
    # Issue when N is small
    success <- ifelse(exps_out[id_list] > S[id_list], S[id_list],
                      exps_out[id_list])

    for(i in seq_along(id_list)) {
      if(success[i] > 0) {
        outcome[sample(which(exposed$row_id == id_list[i]),
                       success[i],
                       replace = FALSE)] <- 1
      }
    }
  } else {

    # track who exposure was allocated to
    contact <- rep("", nexps)
    outcome <- rep(0, nexps)

    for(i in seq_along(id_list)) {
      # Sample the states
      ind <- id_list[i]
      # taking out the infectious dog doing the biting (I[i] - 1)
      states <- rep(c("S", "E", "I", "V"),
                    c(S[ind], E[ind], I[ind] - 1, V[ind]))
      contact[which(exposed$row_id == ind)] <- sample(states, size = exps[ind],
                                               replace = FALSE)
      outcome[which(contact == "S")] <- 1
    }
  }
  exposed$outcome <- outcome
  exposed$contact <- contact

  return()
}

#' Simulate biting and movement
#'
#' @param secondaries
#' @param ids
#' @param dispersal_fun
#' @param counter
#' @param cells_pop
#' @param x_topl
#' @param y_topl
#' @param tstep
#' @param sequential
#' @param allow_empties
#' @param leave_district
#' @param max_tries
#' @param res_m
#' @param nrow
#' @param ncol
#' @param ...
sim_bites <- function(secondaries, ids, dispersal_fun, counter, res_m,
                      row_id, cell_id, t_infectious,
                      cells_pop, nrow, ncol, x_topl, y_topl, tstep = t,
                      sequential = TRUE, allow_empties = TRUE,
                      leave_district = TRUE, max_tries = 100,
                      ...) {

  progen_ids <- rep(ids, secondaries)
  origin_x <- rep(x_coord, secondaries)
  origin_y <- rep(y_coord, secondaries)

  if(sequential == TRUE) {

    # first movement index of each progenitor
    inds <- mapply(function(x) which(progen_ids == x)[1], x = ids)
    out <- vector("list", length(progen_ids))

    for (i in seq_along(progen_ids)) {
      if (i %in% inds) { # need progenitor coords for 1st movement
        x <- origin_x[i]
        y <- origin_y[i]
        path_id <- 1
      } else {
        x <- out$x_new[i - 1]
        y <- out$y_new[i - 1]
        path_id <- path_id + 1
      }

      accept <- 0
      tries <- 0
      while (accept == 0 & tries < max_tries) {

        move <- sim_movement(angle = runif(n = 1, min = 0, max = 2*pi),
                             distance = dispersal_fun(1),
                             x0 = x, y0 = y, x_topl,
                             y_topl, res_m, ncol,
                             nrow, cells_pop, path_id)
        accept <- accept(leave_district, allow_empties, move$within,
                         move$populated)

        tries <- tries + 1

      }
      out[[i]] <- move
    }

    out <- rbindlist(out)

  } else {

    nmoves <- sum(secondaries)

    # get origin from infector if kernel movements
    out <- sim_movement(angle = runif(n = nmoves, min = 0, max = 2*pi),
                        distance = dispersal_fun(nmoves),
                        x0 = origin_x, y0 = origin_y, x_topl,
                        y_topl, res_m, ncol,
                        nrow, cells_pop)

    accept <- accept(leave_district, allow_empties, within = out$within,
                     populated = out$populated)

    tries <- 0
    while(sum(accept == 0) > 0 & tries < max_tries) {

      out[accept == 0, ] <-
        sim_movement(angle = runif(n = sum(accept == 0), min = 0, max = 2*pi),
                     distance = dispersal_fun(sum(accept == 0)),
                     x0 = origin_x[accept == 0], y0 = origin_y[accept == 0],
                     x_topl, y_topl, res_m, ncol, nrow, cells_pop,
                     path_id = 0)

      accept <- accept(leave_district, allow_empties, within = out$within,
                       populated = out$populated)

      tries <- tries + 1
    }
  }

  # add in ids
  out$id <- counter + 1:length(progen_ids)
  out$progen_id <- progen_ids
  out$t_infected <- rep(t_infectious, secondaries) # tstep became infected
  out$row_id <- row_id[out$cell_id %in% cell_id] # row ids to match to I/E mats
  return(out)
}

# above should inherit params x_topl, y_topl, res_m, ncol, nrow, cells_pop!
# populated & within are boolean
#' Simulate movement from origin, angle, and distance
#'
#' @param angle
#' @param distance
#' @param x0
#' @param y0
#' @param x_topl
#' @param y_topl
#' @param res_m
#' @param ncol
#' @param nrow
#' @param cells_pop
#'
#' @return
#' @export
#'
#' @examples
sim_movement <- function(angle, distance, x0, y0, x_topl,
                         y_topl, res_m, ncol, nrow, cells_pop,
                         path_id) {

  x_coord <- (sin(angle) * distance * 1000) + x0 # convert to m
  y_coord <- (cos(angle) * distance * 1000) + y0

  # This means can go anywhere
  col <- ceiling((x_new - x_topl)/res_m)
  row <- ceiling(-(y_new - y_topl)/res_m)
  cell_id <- row*ncol - (ncol - col)
  populated <- cell_id %in% cells_pop
  within <- row > 0 & row <= nrow & col > 0 & col <= ncol

  return(data.table(x_coord, y_coord, cell_id, populated, within, path_id))

}

#' accept movement made
#'
#' @param leave_district
#' @param allow_empties
#' @param within
#' @param populated
#'
#' @return
#' @export
#'
#' @examples
#' # inherit params: leave_district allow_empties
accept <- function(leave_district, allow_empties,
                   within, populated) {

  accept <- rep(1, length(within))

  if(leave_district == FALSE) {
    accept <- ifelse(within, 1, 0)
  }

  if(allow_empties == FALSE) {
    accept <- ifelse(populated, 1, 0)
  }

  return(accept)
}

sim_incursions <- function(n_incs, cells_pop, cells_all, x_coord_pop, y_coord_pop,
                           counter, tstep, days_in_step = 7) {

  # sample cell ids by n
  cell_id <- sample(cells_pop, n_incs, replace = TRUE)

  # date infectious (this tstep just draw the day!)
  t_infectious <- (sample(1:days_in_step, n_incs, replace = FALSE) + tstep)/days_in_step

  # incursions have a progenitor id of zero
  data.table(id = counter + 1:n_incursions,
             row_id = row_id[cell_id %in% cells_all],
             cell_id, tstep, progen_id = 0,
             x_coord = x_coord_pop[cell_incursions],
             y_coord = y_coord_pop[cell_incursions],
             populated = TRUE,
             within = TRUE, t_infectious)

}

sim_incursions_hardwired <- function(row_ids_empirical, tstep_empirical,
                                     counter, cell_id, x_coord, y_coord,
                                     tstep, days_in_step = 7) {

  # filter list of empirical incursions
  row_id <- row_ids_empirical[tstep_empirical %in% current_tstep]
  n_incs <- length(cell_id)
  # date infectious (this tstep just draw the day!)
  t_infectious <- (sample(1:days_in_step, n_incs, replace = FALSE) + tstep)/days_in_step

  # incursions have a progenitor id of zero
  data.table(id = counter + 1:n_incs,
             row_id, cell_id = cell_id[row_id], tstep, progen_id = 0,
             x_coord = x_coord[row_id],
             y_coord = y_coord[row_id],
             populated = TRUE,
             within = TRUE, t_infectious)

}




