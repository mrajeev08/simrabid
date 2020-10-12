
# return a list of the objects in this environment
# then when you start the function you can do a list2env call of setup
# should be faster than setting it up each time you run the simulation
# especially for large matrices
# nlocs <- 10000
# tmax <- 1000
# start_pop <- round(rnorm(10000, 1000))
# I_seed <- 10
# cells_pop <- sample.int(10000, 5000)
# start_vacc <- round(rnorm(10000, 500))

setup <- function(nlocs, tmax, start_pop, cell_ids) {

  # state matrices
  row_id <- 1:nlocs # don't need to do this function call each time either
  S_mat <- V_mat <- N_mat <- matrix(0, nrow = nlocs, ncol = tmax)
  I_mat <- E_mat <- Matrix(matrix(0, nrow = nlocs, ncol = tmax), sparse = TRUE)
  rows_pop <- row_ids[start_pop > 0]
  cells_pop <- cell_ids[start_pop > 0]

  return(list(row_id = row_id, S_mat = S_mat, I_mat = I_mat,
              cells_pop = cells_pop, rows_pop = rows_pop))
}

init <- function(start_pop, start_vacc, rows_pop, cell_id, nlocs,
                 x_coord, y_coord, I_dt) {

  # Starting pop + sus
  if(length(start_vacc) != 1 & is.integer(start_vacc)) {
    V <- start_vacc
  } else {
    V <- rbinom(n = nlocs, size = start_pop, prob = start_vacc)
  }

  S <- start_pop - V

  # Seed cases at t0 and create data.table
  row_id <- sample(rows_pop, I_seeds, replace = FALSE) # which rows are populated

  # build data table for I_dt
  I_dt <- data.table(id = 0, cell_id = NA, row_id = NA,
                     tstep = NA, progen_id = NA, x_coord = NA,
                     y_coord = NA, populated = NA, within = NA,
                     t_infectious = NA)
  I_dt <- I_dt[rep(I_dt[, .I], 1000)]

  I_init <- data.table(id = 1:I_seeds, cell_id = cell_id[row_id],
                     row_id,
                     tstep = 0, progen_id = 0,
                     x_coord = x_coord[row_id],
                     y_coord = y_coord[row_id],
                     populated = TRUE, within = TRUE,
                     t_infectious = 0)

  I_dt[id] <- I_init # this should get updated in global environment

  # Summarize and put into I!
  I <- tabulate(I_dt$row_id, nbins = nlocs)
  E <- rep(0, nlocs)

  return(list(S = S, V = V, I = I, E = E, I_dt = I_dt))
}

# Doubles this
double_I <- function(I_dt) {
  I_skeleton <- data.table(id = 0, cell_id = NA, row_id = NA,
                           tstep = NA, progen_id = NA, x_coord = NA,
                           y_coord = NA, populated = NA, within = NA,
                           t_infectious = NA)
  I_dt <- rbind(I_dt, I_skeleton[rep(I_skeleton[, .I], nrow(I_dt))])
  return(I_dt)
}

# I_dt <- data.table(id = 0, cell_id = NA, row_id = NA,
#                          tstep = NA, progen_id = NA, x_coord = NA,
#                          y_coord = NA, populated = NA, within = NA,
#                          t_infectious = NA)
# double_I(I_dt)
