#' Return a list of named objects from parent environment
#'
#' This function is an example of a summary function that can be
#' passed to `simrabid`. If passing to the summary_fun arg in `simrabid`,
#' it will be called with no arguments,
#' so all arguments must have a default or be empty. This is also a useful
#' helper function for writing your own summary functions, just wrap it in
#' list2env and set envir_num to the parent environment frame number where
#' the objects you're looking for are (see the source for `check_sim` for an example).
#'
#' @param names character vector of object names in parent environment
#'
#' @return a list of objects from simrabid specified in names
#' @export
#'
#' @keywords summarize
use_mget <- function(names = c("I_dt", "I"),
                     envir_num = 1) {
  mget(names, envir = parent.frame(envir_num))
}

#' Return a list of all objects from parent environment
#'
#' This function is an example of a summary function that can be
#' passed to `simrabid`. These functions will be called with no arguments,
#' so all arguments must have a default or be empty.
#'
#' @return a list of all objects from the simrabid environment
#' @export
#'
#' @keywords summarize
#'
return_env <- function() {
  out <- as.list.environment(parent.frame(1))
  return(out)
}


#' Quick summary stats of time series
#'
#' Sanity check and also test different settings in the simulation framework.
#'
#' @param names objects to get from `simrabid` environment
#'
#' @return time series of monthly cases (total, local, and detected), susceptibles,
#'  population size, and coverage.
#'
#' @keywords summarize
#'
check_sim <- function(names = c("I_dt", "S_mat", "N_mat", "tmax", "days_in_step")) {

  # Get the objects you need from the environment above this one
  list2env(use_mget(names, envir_num = 2), envir = environment())

  # Filter to infected (this also means it will no longer point to I_dt)
  I_dt <- I_dt[infected == TRUE]

  # aggregate cols by timestep
  ncols_sum <- floor(30.5 / days_in_step)

  I_total <- tabulate(floor(I_dt$t_infectious), tmax)
  I_detected <- tabulate(floor(I_dt$t_infectious[I_dt$detected == TRUE]), tmax)
  I_local <- tabulate(floor(I_dt$t_infectious[I_dt$progen_id > 0]), tmax)

  # Summarize monthly cases
  I_total <- sum_to_month(I_total, nc = ncols_sum)
  I_local <- sum_to_month(I_local, nc = ncols_sum)
  I_detected <- sum_to_month(I_detected, nc = ncols_sum)

  # Summarize S/N monthly
  S <- inds_to_month(colSums(S_mat), nc = ncols_sum)
  N <- inds_to_month(colSums(N_mat), nc = ncols_sum)
  cov <- S/N

  # distances
  mean_dist_m <- mean(dist_linked(I_dt))
  max_dist_m <- max(dist_linked(I_dt))

  # times
  mean_times_days <- mean(times_linked(I_dt))
  max_times_days <- mean(times_linked(I_dt))

  # data.table with tstep as additional covariate
  data.table(month = 1:length(I_total), I_local, I_total, I_detected,
             S, N, cov, mean_dist_m, max_dist_m, mean_times_days, max_times_days)
}

# summarize every 4 columns
sum_to_month <- function(ts, nc = 4) {

  if(!is.null(dim(ts))) {
    ts <- colSums(ts)
  }
  nv <- length(ts)
  if (nv %% nc)
    ts[ceiling(nv / nc) * nc] <- NA
  colSums(matrix(ts, nc), na.rm = TRUE)

}

# summarize every 4 columns
inds_to_month <- function(ts, nc = 4) {

  maxl <- floor(length(ts) / nc)
  ts <- ts[seq(1, length(ts), by = nc)]
  ts[1:maxl]

}

# average euclidean distance between linked cases
dist_linked <- function(I_dt) {

  coords <- I_dt[, c("x_coord", "y_coord", "id", "progen_id")]
  coords <- coords[coords, on = c("id" = "progen_id")][!is.na(progen_id)]
  coords[, dist_m := sqrt((x_coord - i.x_coord)^2 + (y_coord - i.y_coord)^2)]
  return(coords$dist_m)
}

# average times between linked cases
times_linked <- function(I_dt) {

  times <- I_dt[t_infected > 0][, ts := t_infectious - t_infected]

  return(times$ts)
}




