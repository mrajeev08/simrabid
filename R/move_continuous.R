
#' Simulate individual movements
#'
#' Simulates movemement of individuals in continuous space.
#'
#' This simulates movement and uses the top-left coordinates and 1-based indexing
#' of a UTM raster to identify the grid cell moved to by an individual.
#'
#' @param dist_m numeric vector or value, distance to move
#' @param angle numeric vector [0, 360] or value, angle at which to move at
#' @param x0 numeric vector or value, the x origin of the infected individual
#' @param y0 numeric vector or value, the y origin of the infected individual
#' @param x_topl numeric, the top left x coordinate of the grid on which
#'   movement is being simulated
#' @param y_topl numeric, the top left y coordinate of the grid on which
#'   movement is being simulated
#' @param res_m numeric, the grid cell resolution in meters
#' @param ncols numeric, the number of columns in the grid
#' @param nrows numeric, the number of rows in the grid
#' @param ncells numeric, the number of total cells (ncol * nrow)
#' @param cells_block integer vector, the cell ids of grid cells where movements are invalid to,
#'  for use with max_tries
#'  @param cells_out_bounds integer vector, the cell ids of grid cells which are outside the bounds/
#'  not covered by the area being simulated
#' @param path integer vector or value, if sequential is FALSE then 0, if TRUE then
#'  the path id (i.e. which step in the series of movements) to pass through and assign to exposed
#' @param leave_bounds boolean, are movements to outside of the boundaries of the are being simulated are valid
#' @param allow_invalid boolean, are movements to empty patches (i.e. with no dogs or bodies of water, etc.) valid
#' @param sequential boolean, if TRUE then movements are sequential, if FALSE, then
#'   movements are kernel based
#'
#' @return a list (if sequential) or a data.table (if kernel based) with the resulting x and y coordinates, cell ids, whether
#'   the movement was to a invalid cell, whether the movement fell within the bounds
#'   of the simulation, and the path id for the exposed.
#' @keywords move
#' @export
#'
sim_movement_continuous <-
  function(dist_m, angle, dispersal_fun, x0, y0, x_topl,
           y_topl, res_m, ncols, nrows, ncells, cells_block, cells_out_bounds, path,
           leave_bounds, allow_invalid, max_tries, sequential, weights = NULL,
           params) {

    if(sequential) {

      accept <- 0
      tries <- 0

      while (!accept & tries < max_tries) {

        if(tries == 0) {
          dist_now <- dist_m
          angle_now <- angle
        } else {
          dist_now <- dispersal_fun(1, params)
          angle_now <- angle_fun(1)
        }

        out <- move_continuous(dist_m = dist_now, angle = angle_now,
                               x0, y0, x_topl, y_topl,
                               res_m, ncols, nrows, ncells, cells_block,
                               cells_out_bounds, path = 0)

        accept <- valid(leave_bounds, allow_invalid, outbounds = out$outbounds,
                         invalid = out$invalid)

        tries <- tries + 1
      }

    } else {

      # get origin from infector if kernel movements
      out <- as.data.table(
        move_continuous(dist_m, angle, x0, y0, x_topl, y_topl, res_m, ncols,
                        nrows, ncells, cells_block, cells_out_bounds, path = 0)
      )

      accept <- valid(leave_bounds, allow_invalid, outbounds = out$outbounds,
                       invalid = out$invalid)

      tries <- 0
      retry <- sum(!accept)

      while(retry > 0 & tries < max_tries) {
        inds <- !accept

        out[inds, ] <-
          as.data.table(
            move_continuous(dist_m = dispersal_fun(retry, params),
                            angle = angle_fun(retry),
                            x0 = x0[inds], y0 = y0[inds],
                            x_topl, y_topl, res_m, ncols, nrows, ncells,
                            cells_block,
                            cells_out_bounds,
                            path = 0L))

        accept <- valid(leave_bounds, allow_invalid, outbounds = out$outbounds,
                         invalid = out$invalid)

        tries <- tries + 1
      }
    }

    return(out)

  }



#' Movement in continuous space
#'
#'Simulates individual movements
#'
#' @inheritParams sim_movement_continuous
#'
#' @return a list with the resulting x and y coordinates, cell ids, whether
#'   the movement was to a invalid cell, whether the movement fell within the bounds
#'   of the simulation, and the path id for the exposed (can either be length 1 or length of dist_m)
#'
#' @keywords internal
#'
#'
move_continuous <- function(dist_m, angle,
                            x0, y0, x_topl,
                            y_topl, res_m, ncols, nrows, ncells,
                            cells_block, cells_out_bounds,
                            path) {

  x_coord <- (sin(angle) * dist_m) + x0 # convert to m
  y_coord <- (cos(angle) * dist_m) + y0

  # Get cell id
  cell_id <- get_cellid(x_coord, y_coord, res_m, x_topl, y_topl, ncols, nrows, ncells)
  invalid <- cell_id %in% cells_block
  outbounds <- !(cell_id %in% 1:ncells) | cell_id %in% cells_out_bounds

  return(list(x_coord = x_coord,
              y_coord = y_coord, cell_id = cell_id,
              invalid = invalid, outbounds = outbounds, path = path))

}

#' Accept a simulated movement
#'
#' Helper function to determine if a simulated movement is valid.
#' To do: tests (same length as outbounds and boolean no NAs)
#'
#' @inheritParams sim_movement_continuous
#' @param outbounds whether cell id is out-of-bounds of the area being simulated
#' @param invalid whether cell id is considered an ivalid area (i.e. not populated or other barrier to movement)
#'
#' @return a boolean vector of length `outbounds`/`invalid` corresponding to whether
#'   a movement is accepted as valid based on the arguments `leave_bounds` / `allow_invalid`
#'
#' @keywords move internal
#'
valid <- function(leave_bounds, allow_invalid,
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

#' Helper function for generating angles to move at
#'
#' @keywords move internal
#'
angle_fun <- function(n) {
  runif(n, min = 0, max = 2*pi)
}
