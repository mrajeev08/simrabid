# Parameter distributions
# think about how flexible you want this to
# all functions should have n & param list
# pass parameters for each into param list in sim.IBM
# if you're fitting these then you want to be able to pass params easily in
# the args of the function (I think this works well? Just need unique names?)

# if you want to change this then just need to change the param list & the functions accordingly

# secondary function (should take from param list (n & params))
secondary_fun <- function(n, params = c(R0 = 1.2, k = 1,
                                      max_secondaries = 100)) {
  secondaries <- rnbinom(n, size = params$k, mu = params$R0)

  # return constrained!
  return(ifelse(secondaries > params$max_secondaries,
                params$max_secondaries, secondaries))
}

# dispersal function
# defaults to Townsend et al. 2013 plos ntds
dispersal_dist <- function(n, params = c(disp_shape = 1.46,
                                         disp_scale = 16.1)) {

  rgamma(n, shape = params$disp_shape, scale = params$disp_scale)

}

# defaults to Townsend et al. 2013 plos ntds
# could also use an empirical sample of generation times per rebecca
gen_time <- function(n, params = c(inc_shape = 1.46, inc_scale = 16.1),
                     ...) {

  rgamma(n, shape = params$inc_shape, scale = params$inc_scale)

}

# generation time to the decimal timestep
t_infectious <- function(t_infected, days_in_step = 7,
                         generation_fun) {

  return(t_infected + generation_fun(n)/days_in_step)

}
