# Summary stat functions
# To try
# max # of cases
# mean # of cases
# median # of cases
# pop growth overall (at district scale) (so that its okay if timed out)
# acf (1 - 10)
# mean distance between cases at 1 month (Rcpp)
# this mean normalized by mean overall (Rcpp)
# spatial loss function (mismatch between where cases are?)
# ks test statistic for frequency of case counts per month (temporal bit) (so that its okay if timed out)


# I_dt <- data.table(tstep = runif(1e4, 0, 4*24),
#                    x_coord = runif(1e4, 0, 1),
#                    y_coord = runif(1e4, 0, 1))
#

# Incidence stats
# data = inc_hist + breaks (which should be 1 longer than length of inc_hist)
# data = # of cases per cell (order by the cell_id)
# do this with example of different sizes & benchmark
inc_stats <- function(I_dt, tmax, ncells, data) {

  # filter I_dt to successful transmission events & detected cases

  # group by month and count (& keyby sorts it)
  # Complete this so that months with zero cases are accounted for
  month_max <- ceiling(tmax/4)
  I_ts <- I_dt[, .N, keyby = floor(tstep/4)][data.table(floor = 1:month_max)]
  setnafill(I_ts, fill = 0, cols = 'N')

  # return stats on the time series (should not have NAs)
  max_I <- max(I_ts$N)
  median_I <- median(I_ts$N)
  mean_I <- mean(I_ts$N)

  # temporal corr
  acfs <- as.vector(acf(I_ts$N, lag.max = 10, plot = FALSE)$acf)[-1]
  names(acfs) <- paste0("acf_lag", 1:10)

  # ks discrete statistic
  inc_hist <- hist(I_ts$N, breaks = data$breaks)$count
  ks_stat <- max(abs(inc_hist - data$inc_hist))
  hist_ss <- sum((I_cell$N/sum(I_cell$N) - data$cases_by_cell/sum(data$cases_by_cell))^2)
  hist_rmse <- sqrt(mean((I_cell$N/sum(I_cell$N) - data$cases_by_cell/sum(data$cases_by_cell))^2))

  # spatial corr
  normalized <- mean(dist(cbind(I_dt$x_coord, I_dt$y_coord)))
  mean_dist_4wks <- get_mean_dist(t_dt = I_dt[, .(tstep, x_coord, y_coord)],
                             t_window = 8, samp_max = 1e4)
  mean_dist_4wks_norm <- mean_dist/normalized
  mean_dist_8wks <- get_mean_dist(t_dt = I_dt[, .(tstep, x_coord, y_coord)],
                             t_window = 8, samp_max = 1e4)
  mean_dist_8wks_norm <- mean_dist/normalized


  # loss functions
  # will already be ordered by month due to keyby
  # temporal
  temp_rmse <- sqrt(mean((I_ts$N - data$cases_by_month)^2))
  temp_ss <- sum((I_ts$N - data$cases_by_month)^2)

  # spatial
  I_cell <- I_dt[, .N, keyby = cell_id][data.table(cell_id = 1:ncells)][cell_id %in% 1:ncells]
  spat_rmse <- sqrt(mean((I_cell$N - data$cases_by_cell)^2))
  spat_ss <- sum((I_cell$N - data$cases_by_cell)^2)

  spat_loss <- mean(abs(I_cell$N - data$cases_by_cell))

  return(c(list(max_I = max_I, median_I = median_I, mean_I = mean_I,
                ks_stat = ks_stat,
                hist_ss = hist_ss, hist_rmse = hist_rmse,
                mean_dist_4wks = mean_dist_4wks,
                mean_dist_8wks = mean_dist_8wks,
                mean_dist_4wks_norm = mean_dist_4wks_norm,
                mean_dist_8wks_norm = mean_dist_8wks_norm, spat_rmse = spat_rmse,
                spat_ss = spat_ss, temp_rmse = temp_rmse,
                temp_ss = temp_ss), as.list(acfs)))

}


# Fast acf

# Mean distance between cases within one month of each other
# (lag in 1 direction only to avoid mutliple reps?)

# Normalized by mean distance between all cases (dist_mat mean)

# Naive version with coords (beyond 10,000 coords this becomes very slow)
# tstep <- runif(1e4, 0, 4*24)
# x_coord <- runif(1e4, 0, 1)
# y_coord <- runif(1e4, 0, 1)
# mean_dist_window(x_coord, y_coord, tstep, 4)

mean_dist_window <- function(I_dt, t_window = 4,
                             normalize = TRUE) {

  mu_dist <- rep(NA, length(x_coord))

  for(i in seq_len(length(x_coord))) {

    diff_t <- tstep[i] - tstep

    # filter to ones with cases in preceding twindow
    within <- diff_t > 0 & diff_t < t_window

    if(sum(within) > 0 ) {
      mu_dist[i] <- mean((x_coord[i] - x_coord[within])^2 + (y_coord[i] - y_coord[within])^2)
    } else {
      next
    }

  }

  mean_dist <- mean(mu_dist, na.rm = TRUE)

  if(normalize) {
    mean_dist <- mean_dist/mean(dist(cbind(x_coord, y_coord)))
  }

  return(mean_dist)

}


mean_dist_dt<- function(I_dt, t_window = 4,
                             normalize = TRUE) {

  mu_dist <- rep(NA, length(x_coord))

  for(i in seq_len(length(x_coord))) {

    diff_t <- tstep[i] - tstep

    # filter to ones with cases in preceding twindow
    within <- diff_t > 0 & diff_t < t_window

    if(sum(within) > 0 ) {
      mu_dist[i] <- mean((x_coord[i] - x_coord[within])^2 + (y_coord[i] - y_coord[within])^2)
    } else {
      next
    }

  }

  mean_dist <- mean(mu_dist, na.rm = TRUE)

  if(normalize) {
    mean_dist <- mean_dist/mean(dist(cbind(x_coord, y_coord)))
  }

  return(mean_dist)

}

# Get the x & ycoords for up to N lags & then filter to the diff in times and take the means
# create a data.table with the tsteps
# this is the fastest!
# set it to do max 10k samples
dist_fun <- function(x_coord, y_coord, x_coord_to, y_coord_to) {

 mean((x_coord - x_coord_to)^2 + (y_coord - y_coord_to)^2)

}

get_mean_dist <- function(t_dt = I_dt[, .(tstep, x_coord, y_coord)],
                          t_window = 4, samp_max = 1e4) {

  if(nrow(t_dt) > samp_max) {
    t_dt <- t_dt[sample(.N, samp_max)]
  }

  c_dt <- t_dt[ , .(max = tstep,
                      min = tstep - 4,
                      x_coord_to = x_coord,
                      y_coord_to = y_coord)]

  new <- t_dt[c_dt,
                    on = .(tstep < max, tstep >= min),
                    allow.cartesian = TRUE, nomatch = NULL]

  return(dist_fun(new$x_coord, new$y_coord, new$x_coord_to, new$y_coord_to))

}


