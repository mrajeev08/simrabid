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

  secondaries <- rnbinom(n, size = 1/params$k, mu = params$R0)

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
                              params) {

  rlnorm(n, meanlog = params$disp_meanlog, sdlog = params$disp_sdlog)

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
                           params) {

  rlnorm(n, meanlog = params$serial_meanlog, sdlog = params$serial_sdlog)

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
                          params) {

  rweibull(n, shape = params$steps_shape, scale = params$steps_scale)

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
                         serial_fun, params) {

  t_infects <- t_infected + serial_fun(n, params)/days_in_step
  t_infects <- ifelse(floor(t_infects) == floor(t_infected),
                      floor(t_infected) + 1,
                      t_infects)
  return(t_infects)

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
                              params = list(detect_alpha = 80.1,
                                            detect_beta = 8.9)) {
  probs <- rbeta(1, params$detect_alpha, params$detect_beta)
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
                                      params = list(detect_alpha = 80.1,
                                                    detect_beta = 8.9)) {
  I_dt[, month := floor(t_infectious/4)][,
                                         detect_prob := rbeta(1,
                                                              params$detect_alpha,
                                                              params$detect_beta),
                                  by = month]
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

