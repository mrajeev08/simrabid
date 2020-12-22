dispersal_cdf <- function(x, params = list(disp_shape = 1.46,
                                             disp_scale = 16.1)) {
  pgamma(x, shape = params$disp_shape, scale = params$disp_scale)
}


# Do this inside start up?
get_movemat <- function(x_coords = start_up$x_coord, y_coords = start_up$y_coord,
                        cdf = dispersal_cdf,
                        res_m = start_up$res_m,
                        edge_inds = edge_to_row,
                        edge_fun = "min", leave_bounds = TRUE) {

  # distance matrix (will work in Easting/Northing space)
  dist_mat <- floor(as.matrix(dist(cbind(x_coords, y_coords)))/1000)

  # max distance between any two cells
  max_km <- max(dist_mat) + 1
  ints <- seq(0, max_km, by = res_m/1000)

  # first index is 0-1, then 1-2, etc.
  probs <- diff(dispersal_cdf(ints))
  prob_mat <- matrix(probs[dist_mat + 1], ncol = ncol(dist_mat),
                     nrow = nrow(dist_mat)) # symmetric matrix so can be by col
  # fix diagonal
  diag(prob_mat) <- probs[1] # 0 - 1 bin

  # if rabid dogs are allowed to leave SD
  if(leave_bounds) {
    # get the edge distance for each row
    edge_dist <- apply(dist_mat, 1,
                       function(x) get(edge_fun)(x[edge_inds]))
    edge_probs <- 1 - dispersal_cdf(edge_dist + 1) # upper bound beyond that cell

    prob_mat <- cbind(prob_mat, edge_probs) # last col = out of district
  }

  # turn it into a list
  prob_list <- asplit(prob_mat, 1)

  return(prob_list)
}


sim_movement_prob <- function(secondaries, ids, rows, prob_list,
                              x_coords, y_coords, cell_ids,
                              t_infectious,
                              leave_bounds = FALSE, max_rows,
                              max_to_draw = length(prob_list),
                              counter, sequential) {

  progen_ids <- rep(ids, secondaries)
  row_ids <- rep(rows, secondaries)
  out <- rep(0, length(progen_ids))

  if(!sequential) {
    path <- 0
    # Write the below section in Rcpp and see if it's faster...
    for(i in 1:length(row_ids)) {
      out[i] <- sample.int(max_to_draw, 1, prob = prob_list[[row_ids[i]]])
    }
  } else {

    # first movement index of each progenitor
    inds <- match(ids, progen_ids) # returns first match of id in progen id
    inds <- inds[!is.na(inds)]
    path <- rep(0, length(progen_ids))

    # Write the below section in Rcpp and see if it's faster...
    for (i in seq_along(progen_ids)) {
      if (i %in% inds) { # need progenitor coords for 1st movement
        origin <- row_ids[i]
        path[i] <- 1
      } else {
        if(row_moved > max_rows) {
          origin <- row_ids[i] # spit back into disrict if you left it...
        } else {
          origin <- row_moved
        }
        path[i] <- path[i - 1] + 1
      }
      row_moved <- sample.int(max_to_draw, 1, prob = prob_list[[origin]])
      out[i] <- row_moved
    }
  }

  if(leave_bounds) {
    within <- fifelse(out > max_rows, FALSE, TRUE)
  } else {
    within <- TRUE
  }

  out <- data.table(id = counter + 1:length(progen_ids),
                    row_id = out,
                    cell_id = cell_ids[out],
                    x_coord = x_coords[out],
                    y_coord = y_coords[out],
                    progen_ids, path,
                    within,
                    t_infected = rep(t_infectious, secondaries),
                    contact = "M",
                    infected = FALSE)
  return(out)

}



