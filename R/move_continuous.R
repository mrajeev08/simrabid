
#' Simulate individual movements
#'
#' Simulates movemement of individuals in continuous space.
#'
#' This simulates movement and uses the top-left coordinates and 1-based indexing
#' of raster cell ids to identify the grid cell moved to by an individual.
#'
#' @param angle numeric vector [0, 360] or value, angle at which to move at
#' @param dist_m numeric vector or value, distance to move
#' @param x0 numeric vector or value, the x origin of the infected individual
#' @param y0 numeric vector or value, the y origin of the infected individual
#' @param x_topl numeric, the top left x coordinate of the grid on which
#'   movement is being simulated
#' @param y_topl numeric, the top left y coordinate of the grid on which
#'   movement is being simulated
#' @param res_m numeric, the grid cell resolution in meters
#' @param ncol numeric, the number of columns in the grid
#' @param nrow numeric, the number of rows in the grid
#' @param cells_block integer vector, the cell ids of grid cells where movements are valid to;
#'  for use with max_tries
#'  @param cells_block integer vector, the cell ids of grid cells which are outside the bounds/
#'  not covered by the area being simulated (out of bounds)
#' @param path integer vector or value, if sequential is FALSE then 0, if TRUE then
#'  the path id (i.e. which step of the movements) to pass through and assign to exposed
#'  @param leave_bounds boolean, are movements to outside of the boundaries of the are being simulated are valid
#' @param allow_invalid boolean, are movements to empty patches (i.e. with no dogs) valid
#' @param outbounds boolean vector, whether movement falls out of bounds
#' @param invalid boolean vector, whether movement falls within a invalid area
#'
#' @return a list with the resulting x and y coordinates, cell ids, whether
#'   the movement was to a invalid cell, whether the movement fell within the bounds
#'   of the simulation, and the path id for the exposed.
#' @keywords transmit internal
#'
sim_movement_continuous <-
  function(dist_m, angle, dispersal_fun, x0, y0, x_topl,
           y_topl, res_m, ncol, nrow, ncell, cells_block, cells_out_bounds, path,
           leave_bounds, allow_invalid, max_tries, sequential,
           ...) {

    if(sequential) {

      accept <- 0
      tries <- 0

      while (!accept & tries < max_tries) {

        if(tries = 0) {
          dist_now <- dist_m
          angle_now <- angle
        } else {
          dist_now <- dispersal_fun(1)
          angle_now <- angle_fun(1)
        }

        out <- move_continuous(dist_m = dist_now, angle = angle_now,
                               x0, y0, x_topl, y_topl,
                               res_m, ncol, nrow, ncell, cells_block,
                               cells_out_bounds, path = 0)

        accept <- accept(leave_bounds, allow_invalid, outbounds = out$outbounds,
                         invalid = out$invalid)

        tries <- tries + 1
      }

    } else {

      # get origin from infector if kernel movements
      out <- as.data.table(
        move_continuous(dist_m, angle, x0, y0, x_topl, y_topl, res_m, ncol,
                        nrow, cells_block, cells_out_bounds, path = 0)
      )

      accept <- accept(leave_bounds, allow_invalid, outbounds = out$outbounds,
                       invalid = out$invalid)

      tries <- 0
      retry <- sum(!accept)

      while(retry > 0 & tries < max_tries) {
        inds <- !accept

        out[inds, ] <-
          as.data.table(
            move_continuous(dist_m = dispersal_fun(retry),
                            angle = angle_fun(retry),
                            x0 = x0[inds], y0 = y0[inds],
                            x_topl, y_topl, res_m, ncol, nrow, cells_block,
                            cells_out_bounds,
                            path = 0L))

        accept <- accept(leave_bounds, allow_invalid, outbounds = out$outbounds,
                         invalid = out$invalid)

        tries <- tries + 1
      }
    }

    return(out)

  }

move_continuous <- function(dist_m, angle,
                            x0, y0, x_topl,
                            y_topl, res_m, ncol, nrow, ncell,
                            cells_block, cells_out_bounds,
                            path) {

  x_coord <- (sin(angle) * dist_m) + x0 # convert to m
  y_coord <- (cos(angle) * dist_m) + y0

  # Get cell id
  cell_id <- get_cellid(x_coord, y_coord, res_m, x_topl, y_topl, ncol, nrow)
  invalid <- cell_id %in% cells_block
  outbounds <- !(out$cell_id %in% 1:ncell) | cell_id %in% cells_out_bounds

  return(list(x_coord = x_coord,
              y_coord = y_coord, cell_id = cell_id,
              invalid = invalid, outbounds = outbounds, path = path))

}

#' Accept a simulated movement
#'
#' Helper function to determine if a simulated movement is valid.
#' To do: tests (same length as outbounds and boolean no NAs)
#'

#'
#' @return a boolean vector of length `outbounds`/`invalid` corresponding to whether
#'   a movement is accepted as valid
#'
#' @keywords transmit internal
#'
accept <- function(leave_bounds, allow_invalid,
                   outbounds, invalid) {

  accept <- rep(1, length(outbounds))

  if(!leave_bounds) {
    accept <- ifelse(outbounds, 0, 1)
  }

  if(!allow_invalid) {
    accept <- ifelse(invalid, 0, 1)
  }

  return(accept)
}

# Random angle fun
angle_fun <- function(n) {
  runif(n, min = 0, max = 2*pi)
}
