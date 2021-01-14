cells_away <- function(origin_cell, dist_m, res_m, ncol, nrow, ncells) {

  n_away <- floor(dist_m/res_m)
  origin_row <- ceiling(origin_cell/ncol)
  origin_col <- ifelse(origin_cell %% ncol == 0, ncol, origin_cell %% ncol)

  if(n_away == 0) {
    return(origin_cell)
  } else {
    aways <- data.table(expand.grid(col = -n_away:n_away, row = -n_away:n_away))
    aways <- aways[abs(row) + abs(col) != 0 & pmax(abs(col), abs(row)) == n_away]
    aways <- aways[, c("col", "row") := .(origin_col + col,
                                          origin_row + row)]
    aways <- aways[row <= nrow & row > 0 & col <= ncol & col > 0]
    aways[, cell_id :=  (col - 1) * ncol  + ifelse(row %% ncol == 0,
                                                   ncol,
                                                   row %% ncol)]
    if(nrow(aways) > 0) out <- aways$cell_id else out <- ncells + 1
    return(out)
  }

}

# only need to do this once for each parameter set!
# one covar needs to be 1 and one of the params needs to be beta_0
# covars if non weighted need to be the length of ncells
cell_weights <- function(covars = list(0),
                         params = list(0),
                         cells_pop,
                         allow_empties) {

  weights <- plogis(Reduce('+', Map('*', covars, params)))

  # allow movement into empty patches
  # and allow leaving the bounds of the district

  if(!allow_empties) {
    weights[-cells_pop] <- 0
  }

  # append probability that goes beyond the raster bounds (the minimum)
  # set to zero for out of bounds
  return(c(weights, 0))

}

sim_movement_cells <- function(secondaries, ids = I_now$id,
                                x_coords, y_coords, cell_id,
                                cell_ids, row_ids,
                                cells_pop, weights,
                                t_infectious = I_now$t_infectious,
                                counter = max(I_dt$id),
                                dispersal_fun, res_m,
                                nrow, ncol, ncells,
                                sequential = TRUE) {

  progen_ids <- rep(ids, secondaries)
  cell_id <- rep(cell_id, secondaries)
  xcoords <- x_coords[cell_id]
  ycoords <- y_coords[cell_id]
  out <- rep(0, length(progen_ids))

  if(!sequential) {
    path <- 0
    # Write the below section in Rcpp and see if it's faster...

    for(i in 1:length(cell_id)) {

      distance <- dispersal_fun(1)
      opts <- cells_away(origin_cell = cell_id[i],
                         dist_m = distance * 1000,
                         res_m = res_m, ncol = ncol, nrow = nrow,
                         ncells = ncells)


      if(length(opts) == 1) {
        out[i] <- opts
      } else {
        wts <- weights[opts]

        # set equal opts if somehow all zero
        if(sum(wts) == 0) wts <- rep(1/length(wts), length(wts))

        out[i] <- sample(opts, 1, prob = wts)
      }
    }
  } else {

    # first movement index of each progenitor
    inds <- match(ids, progen_ids) # returns first match of id in progen id
    inds <- inds[!is.na(inds)]
    path <- rep(0, length(progen_ids))

    # Write the below section in Rcpp and see if it's faster...
    for (i in seq_along(progen_ids)) {

      if (i %in% inds) { # need progenitor coords for 1st movement
        origin <- cell_id[i]
        path[i] <- 1
      } else {
        origin <- out[i - 1]
        path[i] <- path[i - 1] + 1
      }

      distance <- dispersal_fun(1)
      opts <- cells_away(origin_cell = origin,
                         dist_m = distance * 1000,
                         res_m = res_m,
                         ncol = ncol, nrow = nrow,
                         ncells = ncells)

      if(length(opts) == 1) {
        out[i] <- opts
      } else {
        wts <- weights[opts]

        # set equal opts if somehow all zero
        if(sum(wts) == 0) wts <- rep(1/length(wts), length(wts))

        out[i] <- sample(opts, 1, prob = wts)
      }
    }
  }

  out <- data.table(id = counter + 1:length(progen_ids),
                    cell_id = out,
                    row_id = row_ids[match(out, cell_ids)],
                    x_coord = x_coords[out],
                    y_coord = y_coords[out],
                    progen_ids, path,
                    within = out %in% 1:ncells,
                    populated = out %in% cells_pop,
                    t_infected = rep(t_infectious, secondaries),
                    contact = "M",
                    infected = FALSE)
  return(out)

}





