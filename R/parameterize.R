# Parameter distributions
# think about how flexible you want this to
# all functions should have n & param list
# pass parameters for each into param list in sim.IBM
# if you're fitting these then you want to be able to pass params easily in
# the args of the function (I think this works well? Just need unique names?)
# use best fits from katie as defaults here!
# if you want to change this then just need to change the param list & the functions accordingly

# secondary function (should take from param list (n & params))
negbinom_constrained <- function(n, params = list(R0 = 1.2, k = 1),
                                 max_secondaries = 100) {
  secondaries <- rnbinom(n, size = params$k, mu = params$R0)

  # return constrained!
  return(ifelse(secondaries > max_secondaries,
                max_secondaries, secondaries))
}

# dispersal function
# defaults to Townsend et al. 2013 plos ntds
# units = km?!
dispersal_gamma <- function(n, params = list(disp_shape = 0.8,
                                         disp_scale = 1)) {

  rgamma(n, shape = params$disp_shape, scale = params$disp_scale)

}

# defaults to Townsend et al. 2013 plos ntds
# could also use an empirical sample of generation times per rebecca
gen_gamma <- function(n, params = list(inc_shape = 1.46, inc_scale = 16.1)) {

  rgamma(n, shape = params$inc_shape, scale = params$inc_scale)

}

# generation time to the decimal timestep
t_infectious <- function(n, t_infected, days_in_step = 7,
                         generation_fun) {

  return(t_infected + generation_fun(n)/days_in_step)

}

# Reporting model
binom_detect <- function(I_dt, params = list(detect_prob = 0.9)) {
  I_dt[, reported := rbinom(.N, size = 1, prob = params$detect_prob)]
}

# Reporting model
binom_beta_detect <- function(I_dt,
                              params = list(alpha = 80.1,
                                            beta = 8.9)) {
  probs <- rbeta(1, params$alpha, params$beta)
  I_dt[, reported := rbinom(.N, size = 1, prob = probs)]
}

# Monthly probabilities
binom_beta_detect_monthly <- function(I_dt, tmax,
                                      params = list(alpha = 80.1, beta = 8.9)) {
  I_dt[, month := floor(tstep/4)][, detect_prob := rbeta(1, params$alpha,
                                                   params$beta), by = month]
  I_dt[, reported := rbinom(.N, size = 1, prob = detect_prob)]
}

# beta_pars from mean
est_beta_pars <- function(mean, var) {
  alpha <- ((1 - mean) / var - 1 / mean) * mean ^ 2
  beta <- alpha * (1 / mean - 1)
  return(params = list(alpha = alpha, beta = beta))
}

# Getting probability from rate ------------------------------------------------

# ## example 1: turn annual waning rate of 0.33 to weekly prob
# get.prob(rate = 0.33, step = 52)
# ## example 2: turn annual birth rate of 0.45 to monthly prob
# get.prob(rate = 0.45, step = 12)

get_prob <- function(rate, step) {
  converted <- (1 + rate)^(1/step) - 1
  return(1 - exp(-converted))
}

