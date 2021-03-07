#' Secondary case distribution
#'
#' This function draws secondary cases from a negative binomial distribution,
#' but constrained to a maximum (`max_secondaries`),
#' where `R0` is the mean and `k` is the inverse of the dispersion parameter.
#' Passed to the param `secondary_fun` in \code{\link{simrabid}}.
#'
#' @param n number to draw
#' @param params list of parameters
#'
#' @return integer vector, number of secondary cases, of length `n`
#' @export
#'
nbinom_constrained <- function(n,
                               params = list(R0 = 1.2, k = 1,
                                             max_secondaries = 100)) {

  secondaries <- rnbinom(n, size = 1/params$k, mu = params$R0)
  secondaries[secondaries > params$max_secondaries] <- params$max_secondaries

  # return constrained!
  return(secondaries)
}

#' Dispersal kernel distribution
#'
#' This function draws dispersal distances (in meters)
#' from a lognormal distribution, pass `param_defaults` to the parameter list
#' for best-fit parameters from Mancy et al.
#' Can be passed to the `dispersal_fun` in \code{\link{simrabid}} (should
#' set parameter `sequential` to FALSE as well in that case).
#'
#' @inheritParams nbinom_constrained
#'
#' @return numeric vector, dispersal distances of length `n`
#' @export
#'

dispersal_lognorm <- function(n,
                              params) {

  rlnorm(n, meanlog = params$disp_meanlog, sdlog = params$disp_sdlog)

}

#' Serial interval distribution
#'
#' This function draws serial intervals (in days)
#' from a lognormal distribution, pass `param_defaults` to the parameter list
#' for best-fit parameters from Mancy et al.
#' Can be passed to the `serial_fun` in \code{\link{simrabid}}.
#'
#' @inheritParams nbinom_constrained
#'
#' @return numeric vector, serial intervals of length `n`
#' @export
#'
serial_lognorm <- function(n,
                           params) {

  rlnorm(n, meanlog = params$serial_meanlog, sdlog = params$serial_sdlog)

}

#' Step length distribution of rabid animal movements
#'
#' This function draws step lengths (in meters)
#' from a lognormal distribution, pass `param_defaults` to the parameter list
#' for best-fit parameters from Mancy et al.
#' Can be passed to the `dispersal_fun` in \code{\link{simrabid}} (should
#' set parameter `sequential` to TRUE as well in that case).
#'
#' @inheritParams nbinom_constrained
#'
#' @return numeric vector of step lengths of length `n`
#' @export
#'
steps_weibull <- function(n,
                          params) {

  rweibull(n, shape = params$steps_shape, scale = params$steps_scale)

}

#' Observation function: single binomial probability
#'
#' This function simulates the observation proccess given a
#' single detection probability. Defaults to 0.9 when passing p
#' `param_defaults`.
#' Can be passed to the `observe_fun` in \code{\link{simrabid}}.
#'
#' @param I_dt the line list of cases from `simrabid`
#' @inheritParams nbinom_constrained
#'
#' @return The reporting function should modify the line list
#' in place within the simrabid function (adding a `detected` column
#' as a boolean vector).
#'
#' @export
#'
binom_detect <- function(I_dt, params = list(detect_prob = 0.9)) {
  I_dt[, detected := rbinom(.N, size = 1, prob = params$detect_prob)]
}


#' Observation function: with beta binomial distribution of detection probabilities
#'
#' This function simulates the observation proccess drawing a detection
#' probability from a beta binomial distribution for each case. Defaults are in
#' `param_defaults`.
#' Can be passed to the `observe_fun` in \code{\link{simrabid}}.
#'
#' @inheritParams binom_detect
#'
#' @inheritSection binom_detect Section return
#' @export
#'
beta_detect <- function(I_dt,
                              params = list(detect_alpha = 80.1,
                                            detect_beta = 8.9)) {
  probs <- rbeta(1, params$detect_alpha, params$detect_beta)
  I_dt[, detected := rbinom(.N, size = 1, prob = probs)]
}

#' Observation function: with beta binomial distribution of monthly
#' detection probabilities
#'
#' This function simulates the observation proccess drawing a detection
#' probability from a beta binomial distribution for each month.
#' Defaults are in `param_defaults`.
#' Can be passed to the `observe_fun` in \code{\link{simrabid}}.
#'
#' @inheritParams binom_detect
#'
#' @inheritSection binom_detect Section return
#' @export
#'
beta_detect_monthly <- function(I_dt, params = list(detect_alpha = 80.1,
                                                    detect_beta = 8.9)) {

  # Get the number of days in step (to keep compatible with other
  # detect functions)
  list2env(use_mget("days_in_step", envir_num = 2), envir = environment())

  # aggregate cols by timestep
  ncols_month <- floor(30.5 / days_in_step)

  I_dt[, month := floor(t_infectious/ncols_month)][,
                                         detect_prob := rbeta(1,
                                                              params$detect_alpha,
                                                              params$detect_beta),
                                  by = month]
  I_dt[, detected := rbinom(.N, size = 1, prob = detect_prob)]
}


