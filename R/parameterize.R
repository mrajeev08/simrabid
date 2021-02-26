# Parameter distributions
# all functions should have two arguments: n & param list

# secondary function (should take from param list (n & params))
#' Title
#'
#' @param n
#' @param params
#' @param max_secondaries
#'
#' @return
#' @export
#'
#' @examples
nbinom_constrained <- function(n,
                               params = list(R0 = 1.2, k = 1),
                               max_secondaries = 100) {

  secondaries <- rnbinom(n, size = params$k, mu = params$R0)

  # return constrained!
  return(ifelse(secondaries > max_secondaries,
                max_secondaries, secondaries))
}

#' Title
#'
#' @param n
#' @param params
#'
#' @return
#' @export
#'
#' @examples

dispersal_lognorm <- function(n,
                              paramsp) {

  rlnorm(n, meanlog = params$meanlog, sdlog = params$sdlog)

}

# serial interval
#' Title
#'
#' @param n
#' @param params
#'
#' @return
#' @export
#'
#' @examples
serial_lognorm <- function(n,
                           params = param_defaults$serial) {

  rlnorm(n, meanlog = params$meanlog, sdlog = params$sdlog)

}

# serial interval
#' Title
#'
#' @param n
#' @param params
#'
#' @return
#' @export
#'
#' @examples
steps_weibull <- function(n,
                          params = param_defaults$steps) {

  rweibull(n, shape = params$shape, scale = params$scale)

}

# generation time to the decimal timestep
#' Title
#'
#' @param n
#' @param t_infected
#' @param days_in_step
#' @param serial_fun
#'
#' @return
#' @export
#'
#' @examples
t_infectious <- function(n, t_infected, days_in_step = 7,
                         serial_fun) {

  return(t_infected + serial_fun(n)/days_in_step)

}

# Reporting model
#' Title
#'
#' @param I_dt
#' @param params
#'
#' @return
#' @export
#'
#' @examples
binom_detect <- function(I_dt, params = list(detect_prob = 0.9)) {
  I_dt[, detected := rbinom(.N, size = 1, prob = params$detect_prob)]
}

# Reporting model with beta binomial distribution
#' Title
#'
#' @param I_dt
#' @param params
#'
#' @return
#' @export
#'
#' @examples
beta_detect <- function(I_dt,
                              params = list(alpha = 80.1,
                                            beta = 8.9)) {
  probs <- rbeta(1, params$alpha, params$beta)
  I_dt[, detected := rbinom(.N, size = 1, prob = probs)]
}

# Reporting model with beta binomial distribution monthly
#' Title
#'
#' @param I_dt
#' @param tmax
#' @param params
#'
#' @return
#' @export
#'
#' @examples
beta_detect_monthly <- function(I_dt, tmax,
                                      params = list(alpha = 80.1, beta = 8.9)) {
  I_dt[, month := floor(tstep/4)][, detect_prob := rbeta(1, params$alpha,
                                                   params$beta), by = month]
  I_dt[, detected := rbinom(.N, size = 1, prob = detect_prob)]
}


# Getting probability from rate
#' Title
#'
#' @param rate
#' @param step
#'
#' @return
#' @export
#'
#' @examples
get_prob <- function(rate, step) {
  converted <- (1 + rate)^(1/step) - 1
  return(1 - exp(-converted))
}

