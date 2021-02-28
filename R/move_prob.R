#' Simulate movement probabilistically
#'
#' @inheritParams sim_bites
#' @param ...
#'
#' @return a list (if sequential) or a data.table (if kernel based) with the resulting x and y coordinates, cell ids, whether
#'   the movement was to a invalid cell, whether the movement fell within the bounds
#'   of the simulation, and the path id for the exposed.
#' @keywords move
#' @export
#'
sim_movement_prob <-
  function(dist_m, angle = NULL, dispersal_fun, x0, y0, x_topl,
           y_topl, res_m, ncols, nrows, ncells, cells_block, cells_out_bounds, path,
           leave_bounds, allow_invalid, max_tries, sequential, weights = NULL,
           params, ...) {

    if(sequential) {

      accept <- 0
      tries <- 0

      while (!accept & tries < max_tries) {

        if(tries == 0) {
          dist_now <- dist_m
        } else {
          dist_now <- dispersal_fun(1, params)
        }

        out <- movement_prob(dist_m = dist_now, weights, x0, y0, x_topl, y_topl, ncells,
                             res_m, ncols, nrows, path,
                             leave_bounds, cells_block,
                             cells_out_bounds)

        accept <- valid(leave_bounds, allow_invalid, outbounds = out$outbounds,
                         invalid = out$invalid)

        tries <- tries + 1
      }

    } else {

      out <- vector("list", length(dist_m))

      for (i in 1:length(out)) {

        accept <- 0
        tries <- 0

        while (!accept & tries < max_tries) {

          if(tries == 0) {
            dist_now <- dist_m[i]
          } else {
            dist_now <- dispersal_fun(1, params)
          }

          out_i <- movement_prob(dist_m = dist_now, weights, x0 = x0[i],
                                 y0 = y0[i], x_topl, y_topl, ncells,
                                 res_m, ncols, nrows, path,
                                 leave_bounds, cells_block,
                                 cells_out_bounds)

          accept <- valid(leave_bounds, allow_invalid, outbounds = out_i$outbounds,
                           invalid = out_i$invalid)

          tries <- tries + 1
        }

        out[[i]] <- out_i

      }

      out <- rbindlist(out)


    }

    return(out)

}


#' Get cell weights based on covariates or a vector of null probabilities
#'
#' This function generates probabilities for each cell given covariates and
#' parameter estimates in a logistic regression framework.
#' As a default it generates uniform probabilities
#' for each cell (0.5, with invalid and out-of-bounds set to zero if leave_bounds
#' or allow_invalid are false respectively). Covaraites and parameters can be passed,
#' but these must include a covariate length 1 equal to 1L and a parameter estimate
#' of length 1 corresponding to the model intercept. The last weight is for those
#' cell ids that fall outside the range of possible cell ids (i.e. not in 1:ncells).
#'
#' @param covars list of vectors corresponding to covariates for each cell,
#'  one should be the intercept (valued at 1)
#' @param params the parameter corresponding to the effect of each covariate,
#'  one should be length 1 and correspond to the intercept
#' @inheritParams sim_bites
#'
#' @return a vector of weights of length ncells + 1
#'
#' @export
#'
#' @example
#' covars <- list(pop = rpois(100, 100), proximity_to_road = runif(100, 0, 10), intercept = 1)
#' params <- list(pop = 1.2, proximity_to_road = -0.5, intercept = -5)
#' weights <- cell_weights(covars = covars, params = params, ncells = 100,
#'                         leave_bounds = TRUE, allow_invalid = FALSE,
#'                         cells_block = c(1, 35, 20), cells_out_bounds = 90:100)
cell_weights <- function(covars = list(0),
                         params = list(0),
                         ncells,
                         leave_bounds,
                         allow_invalid,
                         cells_block,
                         cells_out_bounds) {

  # convert to value between 0 - 1 (defaults to weights of all = 0.5)
  weights <- plogis(Reduce('+', Map('*', covars, params)))

  if(length(weights) == 1) {
    weights <- rep(weights, ncells)
  }

  if(!allow_invalid) {
    weights[cells_block] <- 0
  }

  if(!leave_bounds) {
    weights[cells_out_bounds] <- 0
    # For those that fall outside of all possible cells, set to zero
    weights <- c(weights, 0)
  } else {
    # If outside movements are allowed
    # For those that fall outside of all possible cells, set to min non zero
    weights <- c(weights, min(weights[weights > 0]))
  }


  if(length(weights) != ncells + 1) {
    stop("Cell weights length should be ncells(rast) + 1!")
  }


  return(weights)

}

#' Subset weights of cells to the candidate cell ids
#'
#' @inheritParams sim_bites
#' Weights should be length ncells + 1
#' @return a vector of weights corresponding to the cell_ids passed through
#' @keywords move internal
#'
#'
get_cellweights <- function(weights, cell_ids) {

  weights <- weights[cell_ids]

  return(weights)

}

#' Movement based on probability
#'
#' @inheritParams sim_bites
#'
#' @return list of length 1 with the chosen with sampled cell id,  x and y coordinates, whether
#'   the movement was to a invalid cell, whether the movement fell within the bounds
#'   of the simulation, and the path id for the exposed.
#'
#' @keywords move internal
#'
movement_prob <- function(dist_m,
                          weights, x0, y0, x_topl, y_topl,
                          ncells,
                          res_m, ncols, nrows, path,
                          leave_bounds, cells_block,
                          cells_out_bounds) {

  # Options for movement at distance x (should handle weights at the top)
  # Should return a list!
  opts <- cells_away(x0, y0, dist_m, res_m, ncols, nrows, x_topl, y_topl,
                     ncells)

  # Get the weights for each cell id
  weights <- get_cellweights(weights, opts$cell_id)

  # Pick the cell to move to based on probabilities
  opt_n <- length(opts$cell_id)

  if(opt_n == 1) { # if only one option, then return that option
    out <- unlist(opts)
  } else {
    if(sum(weights) == 0) { # if all weights == 0, don't do prob
      out <- transpose(opts)[[sample(opt_n, size = 1)]]
    } else {
      out <- transpose(opts)[[sample(opt_n, size = 1, prob = weights)]]
    }
  }

  invalid <- out[1] %in% cells_block
  outbounds <- !(out[1] %in% 1:ncells) | out[1] %in% cells_out_bounds

  return(list(cell_id = out[1], x_coord = out[2], y_coord = out[3],
              invalid = invalid, outbounds = outbounds, path = path))

}

#' Generate cell ids for sampling movement to
#'
#' @inheritParams sim_bites
#'
#' @return a list of the candidate cell_ids and the x and y coords of the movement
#' @keywords move internal
#'
cells_away <- function(x0, y0, dist_m, res_m, ncols, nrows,
                       x_topl, y_topl, ncells) {

  n_away <- floor(dist_m/res_m) * 4
  if(n_away == 0) n_away <- 4

  incr <- (2 * pi) / n_away
  angle <- 0 + 1:n_away * incr

  x_coord <- (sin(angle) * dist_m) + x0 # convert to m
  y_coord <- (cos(angle) * dist_m) + y0

  # Get cell ids
  cell_id <- get_cellid(x_coord, y_coord, res_m,
                        x_topl, y_topl, ncols, nrows, ncells)

  return(list(cell_id = cell_id, x_coord = x_coord, y_coord = y_coord))

}

