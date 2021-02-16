sim_movement_prob <-
  function(dist_m, dispersal_fun, x0, y0, x_topl,
           y_topl, res_m, ncol, nrow, cells_block, cells_out_bounds, path,
           leave_bounds, allow_invalid, max_tries, sequential, weights, ...) {


    if(sequential) {

      accept <- 0
      tries <- 0

      while (!accept & tries < max_tries) {

        if(tries = 0) {
          dist_now <- dist_m
        } else {
          dist_now <- dispersal_fun(1)
        }

        out <- movement_prob(dist_m = dist_now, weights, x0, y0, x_topl, y_topl, ncell,
                             res_m, ncol, nrow, path,
                             leave_bounds, cells_block,
                             cells_out_bounds)

        accept <- accept(leave_bounds, allow_invalid, outbounds = out$outbounds,
                         invalid = out$invalid)

        tries <- tries + 1
      }

    } else {

      out <- vector("list", length(dist_m))

      for (i in 1:length(out)) {

        accept <- 0
        tries <- 0

        while (!accept & tries < max_tries) {

          if(tries = 0) {
            dist_now <- dist_m[i]
          } else {
            dist_now <- dispersal_fun(1)
          }

          out_i <- movement_prob(dist_m = dist_now, weights, x0 = x0[i],
                                 y0 = x0[i], x_topl, y_topl, ncell,
                                 res_m, ncol, nrow, path,
                                 leave_bounds, cells_block,
                                 cells_out_bounds)

          accept <- accept(leave_bounds, allow_invalid, outbounds = out_i$outbounds,
                           invalid = out_i$invalid)

          tries <- tries + 1
        }

        out[[i]] <- out_i

      }

      out <- rbindlist(out)


    }

    return(out)

}



# By angle works! How to do this with populated & within for sampling?
# Don't deal with it here...how to separate populated vs. out of district vs. out of bounds (out of district & out of bounds needs to be combined somehow--or treat it all as one and recommend a buffer of epmty cells on all sides! )
# memoize this fun? across clusters?
cells_away <- function(x0, y0, dist_m, res_m, ncol, nrow,
                       x_topl, y_topl) {

  n_away <- floor(dist_m/res_m) * 4
  if(n_away == 0) n_away <- 4

  incr <- (2 * pi) / n_away
  angle <- 0 + 1:n_away * incr

  x_coord <- (sin(angle) * dist_m) + x0 # convert to m
  y_coord <- (cos(angle) * dist_m) + y0

  # Get cell ids
  cell_id <- get_cellid(x_coord, y_coord, res_m, x_topl, y_topl, ncol, nrow)

  return(list(cell_id = cell_id, x_coord = x_coord, y_coord = y_coord))

}

# only need to do this once for each parameter set!
# one param needs to be intercept, i.e. beta_0 (where covar = 1)
# covars if weighted need to be the length of ncells
# tricky for leave bounds because if true and weighted then you need to
# decide what the proabibility is that doggoes leave the bounds
cell_weights <- function(covars = list(0),
                         params = list(0),
                         ncell,
                         leave_bounds,
                         allow_invalid,
                         cells_block,
                         cells_out_bounds) {

  # convert to value between 0 - 1 (defaults to weights of all = 0.5)
  weights <- plogis(Reduce('+', Map('*', covars, params)))

  if(length(weights) == 1) {
    weights <- rep(weights, ncell)
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


  if(length(weights) != ncell + 1) {
    stop("Cell weights length should be ncell(rast) + 1!")
  }


  return(weights)

}

# get cell weights
# dealing with ones that actually are outside of all cells!
get_cellweights <- function(weights, cell_ids) {

  cell_ids[cell_ids < 0] <- length(weights) # means not inside possible cells
  weights <- weights[cell_ids]

  return(weights)

}

# Movement based on probability
movement_prob <- function(dist_m,
                          weights, x0, y0, x_topl, y_topl,
                          ncell,
                          res_m, ncol, nrow, path,
                          leave_bounds, cells_block,
                          cells_out_bounds) {

  # Options for movement at distance x (should handle weights at the top)
  # Should return a list!
  opts <- cells_away(x0, y0, dist_m, res_m, ncol, nrow, x_topl, y_topl)

  # Get the weights for each cell id
  weights <- get_cellweights(weights, cell_ids)

  # Pick the cell to move to based on probabilities
  opt_n <- length(opts$cell_id)

  if(opt_n == 1) { # if only one option, then return that option
    out <- opts
  } else {
    if(sum(weights) == 0) { # if all weights == 0, don't do prob
      out <- transpose(opts)[[sample(opt_n, size = 1)]]
    } else {
      out <- transpose(opts)[[sample(opt_n, size = 1, prob = weights)]]
    }
  }

  out$invalid <- out$cell_id %in% cells_block
  out$outbounds <- !(out$cell_id %in% 1:ncell) | out$cell_id %in% cells_out_bounds
  out$path <- path

  return(out)

}

