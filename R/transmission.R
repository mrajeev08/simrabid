#' Simulate whether a transmission event was successful (i.e. with a susceptible)
#'
#' @param cell_id
#' @param S
#' @param N
#' @param max_cells
#' @param track
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' cell_id <- sample(1:1000, 100, replace = FALSE)
#' max_cells <- 1000
#' S <- sample(1000:1500, max_cells, replace = TRUE)
#' N <- sample(2000:3000, max_cells, replace = TRUE)
#' E <- sample(300:400, max_cells, replace = TRUE)
#' I <- sample(100:200, max_cells, replace = TRUE)
#' V <- sample(500:1000, max_cells, replace = TRUE)
#' sim_trans(cell_id = cell_id, max_cells = max_cells,
#'           S = S, N = N, track = FALSE)
#' sim_trans(cell_id = cell_id, max_cells = max_cells,
#'           S = S, N = N, E, I, V, track = TRUE)

sim_trans <- function(row_id, S, N, max_cells, track = FALSE, ...) {

  exps <- tabulate(row_id, nbins = max_cells)
  id_list <- unique(row_id)

  if(track == FALSE) {
    # probability that a exposure was with a susceptible
    exps_out <- mapply(rbinom, n = 1, size = exps, prob = S/N)
    outcome <- rep(0, length(row_id))
    contact <- NULL

    # Number of successful exposures can't be greater than the susceptibles available
    # Issue when N is small
    success <- ifelse(exps_out[id_list] > S[id_list], S[id_list],
                      exps_out[id_list])

    for(i in seq_along(id_list)) {
      if(success[i] > 0) {
        outcome[sample(which(row_id == id_list[i]),
                       success[i],
                       replace = FALSE)] <- 1
      }
    }
  } else {

    # track who exposure was allocated to
    contact <- rep("", length(row_id))
    outcome <- rep(0, length(row_id))

    for(i in seq_along(id_list)) {
      # Sample the states
      ind <- id_list[i]
      states <- rep(c("S", "E", "I", "V"), # taking out the infectious dog doing the biting (I[i] - 1)
                    c(S[ind], E[ind], I[ind] - 1, V[ind]))
      contact[which(row_id == ind)] <- sample(states, size = exps[ind],
                                               replace = FALSE)
      outcome[which(contact == "S")] <- 1
    }
  }

  return(list(outcome = outcome, contact = contact))
}
# Notes :
# only run this if length(cell_id > 0)
# tests if everything is zero!
# Tests = if exps is zero
# Edge case = if prob S/N is zero or S/E/I/V all zero
# Should be handled before this step?
# use data.table for line list of cases?

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
#' @param ...
sim_bites <- function(secondaries, ids, dispersal_fun, counter, res_m,
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
      } else {
        x <- out$x_new[i - 1]
        y <- out$y_new[i - 1]
      }

      accept <- 0
      tries <- 0
      while (accept == 0 & tries < max_tries) {

        move <- sim_movement(angle = runif(n = 1, min = 0, max = 2*pi),
                            distance = dispersal_fun(1),
                            x0 = x, y0 = y, x_topl,
                            y_topl, res_m, ncol,
                            nrow, cells_pop)
        accept <- accept(leave_district, allow_empties, move$within,
                         move$populated)

        tries <- tries + 1
      }
      out[i] <- move
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
                     x_topl, y_topl, res_m, ncol, nrow, cells_pop)

      accept <- accept(leave_district, allow_empties, within = out$within,
                       populated = out$populated)

      tries <- tries + 1
    }
  }

  # add in ids
  out$id <- counter + 1:length(progen_ids)
  out$progen_id <- progen_ids

  return(out)
}

# above should inherit params x_topl, y_topl, res_m, ncol, nrow, cells_pop!
# populated & within are boolean
sim_movement <- function(angle, distance, x0, y0, x_topl,
                         y_topl, res_m, ncol, nrow, cells_pop) {

  x_new <- (sin(angle) * distance * 1000) + x0 # convert to m
  y_new <- (cos(angle) * distance * 1000) + y0

  # This means can go anywhere
  col <- ceiling((x_new - x_topl)/res_m)
  row <- ceiling(-(y_new - y_topl)/res_m)
  cell <- row*ncol - (ncol - col)
  populated <- cell %in% cells_pop
  within <- row > 0 & row <= nrow & col > 0 & col <= ncol

  return(data.table(x_new = x_new, y_new = y_new, cell = cell,
                    populated = populated, within = within))

}

# inherit params: leave_district allow_empties
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



