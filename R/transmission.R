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

sim_trans <- function(cell_id, S, N, max_cells, track = FALSE, ...) {

  exps <- tabulate(cell_id, nbins = max_cells)
  id_list <- unique(cell_id)

  if(track == FALSE) {
    # probability that a exposure was with a susceptible
    exps_out <- mapply(rbinom, n = 1, size = exps, prob = S/N)
    outcome <- rep(0, length(cell_id))
    contact <- NULL

    # Number of successful exposures can't be greater than the susceptibles available
    # Issue when N is small
    success <- ifelse(exps_out[id_list] > S[id_list], S[id_list],
                      exps_out[id_list])

    for(i in seq_along(id_list)) {
      if(success[i] > 0) {
        outcome[sample(which(cell_id == id_list[i]),
                       success[i],
                       replace = FALSE)] <- 1
      }
    }
  } else {

    # track who exposure was allocated to
    contact <- rep("", length(cell_id))
    outcome <- rep(0, length(cell_id))

    for(i in seq_along(id_list)) {
      # Sample the states
      ind <- id_list[i]
      states <- rep(c("S", "E", "I", "V"), # taking out the infectious dog doing the biting (I[i] - 1)
                    c(S[ind], E[ind], I[ind] - 1, V[ind]))
      contact[which(cell_id == ind)] <- sample(states, size = exps[ind],
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
#' \code{sim_bites} simulates bites and movement on landscape
#'
#' Number of secondary cases is drawn from a _dim, res_m,
#'


sim_bites <- function(secondaries, ids, dispersal_fun, counter,
                      cells_pop, x_topl, y_topl, tstep = t,
                      sequential = TRUE, ...) {

  nmoves <- sum(secondaries)
  progen_ids <- rep(ids, secondaries)
  second_ids <- counter + 1:length(progen_ids)

  origin_x <- rep(x_coord, secondaries)
  origin_y <- rep(y_coord, secondaries)

  distance <- dispersal_fun(nmoves) # length progen_ids
  angle <- runif(n = nmoves, min = 0, max = 2*pi)

  # Convert to a vector and move
  x_moves <- (sin(angle) * distance * 1000)
  y_moves <- (cos(angle) * distance * 1000)

  if(sequential == TRUE) {

    # first movement index of each progenitor
    inds <- mapply(function(x) which(progen_ids == x)[1], x = ids)
    x_new <- y_new <- vector("numeric", length(progen_ids))

    for (i in seq_along(progen_ids)) {

      if (i %in% inds) { # need progenitor coords for 1st movement
        x <- origin_x[i]
        y <- origin_y[i]
      } else {
        x <- last_coords[1]
        y <- last_coords[2]
      }

      # Convert to a vector and move
      x_new[i] <- (sin(angle[i]) * distance[i] * 1000) + x # convert to m
      y_new[i] <- (cos(angle[i]) * distance[i] * 1000) + y

      # This means can go anywhere
      col <- ceiling((x_new - x_topl)/res_m)
      row <- ceiling(-(y_new - y_topl)/res_m)
      cell <- ifelse(row < = row*ncol - (ncol - col)

      within <- ifelse(cell %in% cells_pop, 1, 0)
      last_coords <- c(x_new[i], y_new[i])

    }
  } else {
    # get origin from infector if kernel movements
    x_new <- x_moves + origin_x # convert to m
    y_new <- y_moves + origin_y
  }

}



  sim_trans_constrained <- function()

    if(sequential == FALSE) {

      # Do vectorized
      # get origin from infector if kernel movements
      x_new <- x_moves + origin_x # convert to m
      y_new <- y_moves + origin_y

      ## This means can go anywhere including outside of district
      col <- ceiling((x_new - x_topl)/res_m)
      row <- ceiling(-(y_new - y_topl)/res_m)
      cell <- ifelse(row %in% 1:nrow(grid) & col %in% 1:ncol(grid), grid[row, col], 0)
      ## This bit means only in populated places and doggies don't leave the district

      within <- ifelse(cell %in% cells_pop, 1, 0)

      # Repeat for ones not in district or which went to an empty patch
      if(leave_district == FALSE) {

      }

      if(allow_empties = FALSE) {

      }


    }

  }

  store_coords_all <- vector("list", length(secondaries))

  for (i in 1:length(secondaries)) {
    if (secondaries[i] > 0) {

      store_coords <- vector("list", secondaries[i])

      for (j in 1:secondaries[i]){
        counter <- counter + 1
        if (j == 1) { # need progenitor coords for 1st movement
          origin_x <- x_coord[i]
          origin_y <- y_coord[i]
        } else {
          origin_x <- last_coords[1]
          origin_y <- last_coords[2]
        }
        within <- 0
        while (within == 0) {
          distance <- dispersal_fun(1, shape = dispersalShape,
                               scale = dispersalScale) # Distance to move in km
          angle <- runif(n = 1, min = 0, max = 2*pi)  # Angle to move at

          # Convert to a vector and move
          x_new <- (sin(angle) * distance * 1000) + origin_x # convert to m
          y_new <- (cos(angle) * distance * 1000) + origin_y

          ## This means can go anywhere including outside of district
          col <- ceiling((x_new - x_topl)/res_m)
          row <- ceiling(-(y_new - y_topl)/res_m)
          cell <- ifelse(row %in% 1:nrow(grid) & col %in% 1:ncol(grid), grid[row, col], 0)
          ## This bit means only in populated places and doggies don't leave the district

          within <- ifelse(cell %in% cells_pop, 1, 0)
        }

        store_coords[[j]] <- list(ID = counter, tstep = tstep, x_coord = x_new, y_coord = y_new,
                                  progen_ID = ids[i], path_ID = j, cell_id = cell, sus = NA,
                                  trans = NA_integer_, infectious = 0L, secondaries = NA_integer_)
        last_coords <- c(x_new, y_new)
      }
      store_coords_all[[i]] <- store_coords
    }
  }
  E_coords_now <- rbindlist(lapply(store_coords_all, rbindlist))
  return(list(E_coords_now = E_coords_now, counter = counter))
}
