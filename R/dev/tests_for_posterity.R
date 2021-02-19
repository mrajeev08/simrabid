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

library(data.table)
library(microbenchmark)

# Splitting a vector using base split vs. data.table -------------------
# Depends on the # of groups rather than the vector size
# So if less than 1000 groups/admin units default to using base & split
# & less than 1e5 cells
vec <- rnorm(100000) # 1e5 cells
group <- rep(c(1:10000), each = 10) # 10k groups

microbenchmark::microbenchmark({
  lapply(split(vec, group), sum)
},
{
  test <- data.table(vec, group)
  test[, .(check = sum(vec, na.rm = TRUE)), by = "group"]
}
)

group <- rep(c(1:100), each = 1000) # 100 groups

microbenchmark::microbenchmark({
  t <- lapply(split(vec, group), sum)
  s <- lapply(split(vec*2, group), sum)
  v <- lapply(split(vec*3, group), sum)
},
{
  test <- data.table(v1 = vec, group, v2 = vec*2, v3 = vec*3)
  test[,lapply(.SD, sum, na.rm = TRUE), by = "group",
       .SDcols = c("v1", "v2", "v3")]
}
)

# Is tabulating over a huge # faster than indexing --------------------

# Is subsetting a vector faster than filtering a data.table
# Also need to compare with a key

# Vectorized rbinom
# with R
vec_rbinom <- Vectorize(rbinom)

sizes <- sample(10:100, 1000, replace = TRUE)
probs <- runif(1000)
microbenchmark::microbenchmark(
  mapply(rbinom, n = 1, size = sizes,
         prob = probs),
  vec_rbinom(n = 1, size = sizes,
             prob = probs)
)

# Testing getting raster from cell id
library(raster)
library(microbenchmark)
rast <- raster(nrows=500, ncols=300, xmn=0, xmx=10)

get_cellid <- function(row, col, nrow, ncol) {
  row*ncol - (ncol - col)
}

row_test <- sample(1:50, 10)
col_test <- sample(1:30, 10)
identical(get_cellid(row_test, col_test, nrow = 500, ncol = 300),
          cellFromRowCol(rast, row_test, col_test))

microbenchmark(get_cellid(row_test, col_test, nrow = 500, ncol = 300),
               cellFromRowCol(rast, row_test, col_test))
# testing logicals
microbenchmark(row_test %in% 1:50,
               row_test > 0 & row_test <= 50,
               1000)

origin_x <- rnorm(100)
origin_y <- rnorm(100)
last_coords <- c(100, 200)
inds <- c(1, 5, 10, 15, 21, 30)
i = 10
microbenchmark(
  {
    if (i %in% inds) { # need progenitor coords for 1st movement
      x <- origin_x[i]
      y <- origin_y[i]
    } else {
      x <- last_coords[1]
      y <- last_coords[2]
    }
  },
  {
    x <- ifelse(i %in% inds, origin_x[i], last_coords[1])
    y <- ifelse(i %in% inds, origin_y[i], last_coords[2])
  }
)

x <- last_coords[1]
microbenchmark(test <- last_coords[1],
               test <- x)


# Rbindlist of data.tables vs. setting index

microbenchmark(
  {
    test <- vector("list", 100)
    for(i in 1:length(test)) {
      test[[i]] <- data.table(x = rnorm(1), y = "zz")
    }
    test <- rbindlist(test) # this is much faster
  },
  {
    test <- data.table(x = rep(0, 100), y = rep("", 100))
    for(i in  1:nrow(test)) {
      test[i] <- data.table(x = rnorm(1), y = "zz")
    }
  }
)

x <- rnorm(100)
y <- rnorm(100)
test <- function() mean(x) # takes x from the global environment (maybe this is the function I want?)
test <- function(x) mean(x)
# takes x from the global environment if
test(x = x)
test(x)
test() # wont find it
test(x = c(1, 2, 3)) # will work as expected

summary_funs <- list(means = function() {x <- 1; mean(x)}, adds = function() x + y)
# functions without parameters (whatever parameters are in the function will be searched for in the global environment; not the safest because if you have something in there thats also named in the global environment you have to be careful!)

summary_out <- lapply(summary_funs, function(x) x())

summary_funs <- list(means = function() mean(x), adds = function() x + y)

# test your summary functions (do you have undeclared variables that are not in the function environment of sim?)

microbenchmark(
  sample(100, size = 10), # this actually calls sample.int
  sample.int(100, size = 10), # this is the fastest by a tad
  sample(1:100, size = 10),
  times = 1000
)

# the higher the number of groups the faster the data.table approach is
# but if less than 1000 groups (i.e. in this case timesteps then data.table wins out)
# switch methods when tmax is less than 1000 or does the switching methods actually cost more overhead than it saves?
# faster to do a matrix?
test <- data.table(t = sample(100, 1e4, replace = TRUE))

microbenchmark(
  {
    for(ind in 1:100) {
      tnow <- test[t == ind]
    }
  },
  {
    ts <- split(test, test$t)
    for(i in 1:100) {
      tnow <- ts[[i]]
    }
  }
)

# Outstanding question
# Is it faster to vacc_mat or vacc_dt (i.e. all locations & vacc at all tsteps
# or just when vaccination happens)
# does this change if nlocs change or based on sparsity?
# sparse matrix?

# which vs. match
cell_id <- c(21, 10, 15) # cell numbers chosen
cells <- c(21, 20, 5, 15, 6, 18, 10) # cell numbers
row_id <- 1:length(cells)
rowfromcell <- function(cell_id, cells, rows) {
  rows[match(cell_id, cells)]
}

microbenchmark(
  row_id[cells %in% cell_id], # this is fastest
  rowfromcell(cell_id, cells, rows = row_id),
  row_id[match(cell_id, cells)])

# To test: whether it is faster to store E & I as data.tables vs. lists
t <- 1:1000

microbenchmark(
  (1:1000)[c(1, 5)],
  {
    t[c(1, 5)]
  }
)

# tabulate vs. data.table N
dt <- data.table(row_id = sample(1000, 10000, replace = TRUE))
microbenchmark(
  dt[, .(total = .N), by = "row_id"],
  tabulate(dt$row_id)
)

# growing vs. rbind!
dt <- rbind(data.table(x = 1, y = 2, z = 1),
            data.table(x = rep(NA, 1000), y = rep(NA, 1000), z = rep(NA, 1000)))
dt_1 <- data.table(x = 1, y = 2, z = 1)
microbenchmark(
  {
    for(i in c(1, 11, 21, 31, 41)) {
      dt[i:(i + 10)] <- data.table(x = 1, y = 2, z = 3)
    }
  },
  {
    for(i in  c(1, 11, 21, 31, 41)) {
      dt_1 <- rbind(dt_1, data.table(x = 1, y = 2, z = 3))
    }
  }
)

# Does Na.rm slow things down?
t1 <- c(sample(1e5, 1e4), NA)
t2 <- sample(1e5, 1e4)
microbenchmark(
  max(t1, na.rm = TRUE),
  max(t2)
)

# list2env vs passing arguments
# list2env adds a 5 second overhead per 1e6 sims
test <- function(start_up) {
  list2env(start_up, envir = environment())
  return(ls())
}

other_test <- function(one = start_up[[1]], two = start_up[[2]],
                       three = start_up[[3]], four = start_up[[4]],
                       five = start_up[[5]], six = start_up[[6]]) {
  ls()
}

list2env(start_up)

third_test <- function(row_ids = row_ids,
                       S_mat = S_mat,
                       I_mat = I_mat,
                       E_mat = E_mat) {
  ls()
}

microbenchmark::microbenchmark(test(start_up), other_test(),
                               third_test())

# control flow outside vs inside function
# so minimal a difference put it inside the functions!

test <- function(x = 0, y = rnorm(1e9), z = rnorm(1e9)) {
  if(sum(x) > 0) {
    out <- x*y*z
  } else {
    out <- 0
  }
  return(out)
}

test2 <- function(x = 0, y = rnorm(1e9), z = rnorm(1e9)) {
  x*y*z
}
x <- 0

microbenchmark({
  if(sum(x) > 0) {
    out <- test2(x = x)
  } else {
    out <- 0
  }
},
test())

# comparing data.table vs. base rep the dt way is slower
test <- data.table(x = rnorm(5000), z = rnorm(5000), d = rnorm(5000),
                   y = rep(c(1, 2, 3, 4, 5), each = 1000))

microbenchmark(
  test[rep(test[, .I], y)][, c("x", "z")],
  {
    x <- rep(test$x, test$y)
    z <- rep(test$z, test$y)
  }
)

# benchmark
ids <- sample(5000, 1000)
secondaries <- rnbinom(1000, size = 1, mu = 1.2)
progen_ids <- rep(ids, secondaries)
out <- c(1, cumsum(rle(progen_ids)$lengths) + 1)
out <- mapply(function(x) which(progen_ids == x)[1], x = ids)

ids <- c(1, 5, 7)
secondaries <- c(2, 3, 1)
progen_ids <- rep(ids, secondaries)
microbenchmark(
  out1 <- as.integer(c(1, cumsum(rle(progen_ids)$lengths) + 1)[1:length(ids)]),
  out2 <- mapply(function(x) which(progen_ids == x)[1], x = ids),
  out3 <- match(ids, progen_ids))
identical(out1, out2, out3)

library(glue)
library(microbenchmark)
library(data.table)

hist(rgamma(1000, shape = 0.8, scale = 1.5))


sum(pgamma(1:10, shape = 0.8, rate = 1.5))

tt <- vector(, 100)
for(i in seq_along(1:100)) {
  tt[i] <- var(rnbinom(10000, size = i, mu = 2.5))
}
plot(tt)
abline(h = 2.5, col = "red")

# Modifying by reference in environment ---
tt <- rep(1, 1000)
tt_now <- 0
ids <- 1:100
name <- "tt"
dtt <- data.table(t = rep(1, 1000))
dtt_now <- data.table(t = rep(0, 100))
microbenchmark(
  eval(parse(text = glue("{name}[ids] <- tt_now"))),
  eval(parse(text = paste0(name, "[ids] <- 0"))),
  eval(parse(text = "tt[ids] <- tt_now")), # fastest
  tt[ids] <- tt_now,
  tt <- f(),
  dtt[ids] <- dtt_now # slowest...(gets faster when larger)
)

# Rcpp ----
library(Rcpp)

cppFunction('LogicalVector repf(int n, NumericVector x) {
  LogicalVector check = rep(false, n);
  int l = x.size();
  for(int i = 0; i < l; ++i) {
    int ind = x[i];
    check[ind] = true;
  }
  return check;
}')

repfb <- function(n, x) {
  check <- rep(FALSE, n)
  check[x] <- TRUE
}

library(microbenchmark)

microbenchmark(
  repf(1e4, sample(1e4, 100)),
  repfb(1e4, sample(1e4, 100))
)

cppFunction('IntegerVector tabulate4(const IntegerVector row_id, double nlocs) {
  IntegerVector counts = rep(0, nlocs);
  for (int i = 0; i < nlocs; i++) {
    if (row_id[i] > 0 && row_id[i] <= nlocs)
      counts[row_id[i] - 1]++;
  }
  return counts;
}')

exe <- sample(1000, 500)
max <- 1000

microbenchmark(tabulate(exe, max), tabulate4(exe, max), times = 1000)


