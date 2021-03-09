#' Safer implementation of base sample
#'
#' Wraps sample function in base R so that if options are length 1 and an integer,
#' doesn't default to sample.int. See \code{\link{sample}} for arguments. Only
#' with replace = TRUE (as the function will return the option of length 1 with
#' length of size).
#'
#' @export
safe_sample <- function(x, size, prob = NULL) {

  if(length(x) == 1) {
    return(rep(x, size))
  } else {
    return(sample(x, size = size, prob = prob, replace = TRUE))
  }

}

# testing haversine distance fun
coord_haversine_m <- function(lat, long, distance, brng, r = 6378137) {

  # convert to radians
  lat <- lat * pi / 180
  long <- long * pi / 180

  out_lat <- asin(sin(lat) * cos(distance / r) + cos(lat) * sin(distance / r) * cos(brng))

  out_long <- long + atan2(sin(brng) * sin(distance / r) * cos(lat),
                           cos(distance/ r) - sin(lat) * sin(out_lat))

  out_lat <- out_lat * 180 / pi
  out_long <- out_long * 180 / pi

  return(list(x_coord = out_long, y_coord = out_lat))
}

# getting probability from rate
get_prob <- function(rate, step) {
  converted <- (1 + rate)^(1/step) - 1
  return(1 - exp(-converted))
}

# Generation time to the decimal timestep
t_infectious <- function(n, t_infected, days_in_step = 7,
                         serial_fun, params) {

  t_infects <- t_infected + serial_fun(n, params)/days_in_step
  t_infects <- ifelse(floor(t_infects) == floor(t_infected),
                      floor(t_infected) + 1,
                      t_infects)
  return(t_infects)

}
