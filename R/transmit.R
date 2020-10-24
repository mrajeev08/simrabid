# These functions should handle zero situations gracefully so you don't have to
# or write the conditionals in the IBM part? (because then fewer calls)


#' Simulate whether a transmission event was successful (i.e. was with a susceptible)
#'
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

#' Simulate biting and movement
#'
#' @param secondaries
#' @param ids
#' @param dispersal_fun
#' @param counter
#' @param cells_pop
#' @param x_topl
#' @param y_topl
#' @param sequential
#' @param allow_empties
#' @param leave_district
#' @param max_tries
#' @param res_m
#' @param nrow
#' @param ncol
#' @param ...
sim_bites <- function(secondaries, ids = I_now$id,
                      x_coords = I_now$x_coord, y_coords = I_now$y_coord,
                      t_infectious = I_now$t_infectious,
                      counter = max(I_dt$id),
                      dispersal_fun, res_m,
                      row_ids, cell_ids, cells_pop, nrow, ncol,
                      x_topl, y_topl,
                      sequential = TRUE, allow_empties = TRUE,
                      leave_district = TRUE, max_tries = 100) {

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
    out <- as.data.table(
      sim_movement(
        angle = runif(n = nmoves, min = 0, max = 2*pi),
        distance = dispersal_fun(nmoves),
        x0 = origin_x, y0 = origin_y, x_topl,
        y_topl, res_m, ncol,
        nrow, cells_pop, path = 0)
      ) # not sequential

    accept <- accept(leave_district, allow_empties, within = out$within,
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

      accept <- accept(leave_district, allow_empties, within = out$within,
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

  if(!leave_district) {
    accept <- ifelse(within, 1, 0)
  }

  if(!allow_empties) {
    accept <- ifelse(populated, 1, 0)
  }

  return(accept)
}

sim_incursions_pois <- function(nlocs, rows_pop,
                                params = list(iota = 4),
                                steps = 4) {

  # number of incursions this week
  n_incs <- rpois(1, params$iota/steps)

  # sample cell ids by nlocs
  row_id <- sample(rows_pop, n_incs, replace = TRUE)
  incs <- tabulate(row_id, nbins = nlocs)
  return(incs)
}

# pass through args? or params? or ...?
sim_incursions_hardwired <- function(nlocs, rows_pop,
                                     row_ids = args$row_ids_empirical,
                                     tsteps = args$tstep_empirical,
                                     current_tstep = t) {

  # filter list of empirical incursions
  row_id <- row_ids[tsteps %in% current_tstep & row_ids %in% rows_pop]
  incs <- tabulate(row_id, nbins = nlocs)
  return(incs)
}


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
             populated = TRUE, within = TRUE,
             t_infected = 0, contact = "N",
             infected = TRUE, t_infectious)
}

