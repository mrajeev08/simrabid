library(raster)
library(data.table)
library(microbenchmark)

# distribution of cell ids (should = 1 max)
# sometimes with less than res of cell you'll end up with less than 4 neighbor
# cells because you end up in your own
# it's okay that these are duplicated right? because it means you're in a place where
# it's more likely to stay in your grid cell than leave it...I think.
test_df <- data.frame(x_noise = rnorm(1000, mean = 100, sd = 100), # meters
                      y_noise = rnorm(1000, mean = 100, sd = 100),
                      inds = sample(length(start_up$x_coord), 1000),
                      dist_m = rnorm(1000, mean = 2000, sd = 500))
test <- rep(NA, nrow(test_df))
expected <- rep(NA, nrow(test_df))

for (i in seq_len(nrow(test_df))) {
  tt <- test_df[i, ]
  cells <- cells_away(x0 = start_up$x_coord[tt$inds] + tt$x_noise,
                     y0 = start_up$y_coord[tt$inds] + tt$y_noise,
                     res_m = start_up$res_m,
                     ncol = start_up$ncol,
                     nrow = start_up$nrow,
                     dist_m = tt$dist_m,
                     x_topl = start_up$x_topl,
                     y_topl = start_up$y_topl)
  test[i] <- length(unique(cells))
  exp <- floor(tt$dist_m/start_up$res_m) * 4
  expected[i] <- ifelse(exp == 0, 4, exp)
}

# benchmarking cellfromxy vs. helper function in this
# this much cheaper when calling a bunch (non vectorized)
get_cellid <- function(x_coord = start_up$x_coord,
                       y_coord = start_up$y_coord,
                       res_m = start_up$res_m,
                       ncol = start_up$ncol,
                       nrow = start_up$nrow,
                       x_topl = start_up$x_topl,
                       y_topl = start_up$y_topl) {

  for (i in seq_len(length(x_coord))) {

    col <- ceiling((x_coord[i] - x_topl)/res_m)
    row <- ceiling(-(y_coord[i] - y_topl)/res_m)
    cell_id <- row*ncol - (ncol - col)

  }

  return(cell_id)

}

# when it's sequential (i.e. with prob based)
coords <-  cbind(start_up$x_coord,
                 start_up$y_coord)

get_cellrast <- function(rast_sd = rast, xy = coords) {

  for (i in seq_len(nrow(xy))) {

    cell_id <- cellFromXY(rast_sd, xy[i, ])

  }

  return(cell_id)

}

# testing indexing list vs. data.table
tt <- list(x = rnorm(10), y = rnorm(10), z = rnorm(10))
dt_tt <- as.data.table(tt)

# test and see how this differs between using centroids!
microbenchmark(transpose(tt)[[1]], # fastest! (slightly faster than sapply)
                 dt_tt[1, ], sapply(tt, "[", 2))


# how much does an if else statment add to timing?
# boolean
dist = TRUE
# numeric
num = 13.2

# 145 nanoseconds super cheap don't worry about these!
microbenchmark(if (dist) check = "what",
               if (num < 10) check = "nope")

dist = rbinom(100, size = 1, prob = 0.5)
num = rnorm(100, 10)

microbenchmark(ifelse(dist, "what", "nope"),
               ifelse(num < 10, "what", "nope"))

