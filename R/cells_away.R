sim_movement_prob <-
  function(dispersal_fun, x0, y0, x_topl,
           y_topl, res_m, ncol, nrow, cells_pop, path,
           leave_bounds, allow_empties, max_tries, sequential, weights,
           ...) {

    # get x_coords from indexing (at end add list element rather than within the loop!
    # this as an if statement into bigger movement fun?)

    # list of one if sequential or as.data.table if not sequential
    return(list(x_coord = x_coord,
                y_coord = y_coord, cell_id = cell_id,
                populated = populated, within = within, path = path))

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

  return(cell_id)

}

# only need to do this once for each parameter set!
# one covar needs to be 1 and one of the params needs to be beta_0
# covars if weighted need to be the length of ncells
cell_weights <- function(covars = list(0),
                         params = list(0)) {

  # convert to value between 0 - 1
  weights <- plogis(Reduce('+', Map('*', covars, params)))

  return(weights)

}

# get cell weights for given cell ids accounting for
get_weights <- function(covars = list(0),
                         params = list(0)) {

  # convert to value between 0 - 1
  weights <- plogis(Reduce('+', Map('*', covars, params)))

  return(weights)

}

movement_prob <- function(dispersal_fun,
                          weights, x0, y0, x_topl, y_topl,
                          res_m, ncol, nrow, cells_pop, cells_in_bounds,
                          sequential,
                          allow_empties, leave_bounds, path) {

  # need to handle length of x0 in here as well and sequential in here

  # draw distance
  dist_m <- dispersal_fun(1)

  opts <- cells_away(x0, y0, dist_m, res_m, ncol, nrow, x_topl, y_topl)
  populated <- opts %in% cells_pop
  within <- opts %in% cells_in_bounds

  accept <- accept(leave_bounds, allow_empties, within = within,
                   populated = populated)
  # choose a cell id based on weights

  # Get cell id
  cell_id <- get_cellid(x_coord, y_coord, res_m, x_topl, y_topl, ncol, nrow)
  populated <- cell_id %in% cells_pop
  within <- row > 0 & row <= nrow & col > 0 & col <= ncol

  return(list(x_coord = x_coord,
              y_coord = y_coord, cell_id = cell_id,
              populated = populated, within = within, path = path))

}





