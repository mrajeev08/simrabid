library(Rcpp)
sourceCpp('src/dist.cpp')
sourceCpp('src/profiler.cpp')

tstep <- runif(1e3, 0, 4*24)
x_coord <- runif(1e3, 0, 1)
y_coord <- runif(1e3, 0, 1)

start_profiler("profile.out")
test <- distmean(x_coord, y_coord, tstep, 4)
stop_profiler()
