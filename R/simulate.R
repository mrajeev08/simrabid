
simrabid <- function() {

  # setup (these objects stay the same throughout simulations)
  list2env(setup)

  # init
  list2env(init(start_pop, start_vacc, N, rows_pop, cell_id, nlocs,
                x_coord, y_coord))

  # vaccination steps

  for (t in 1:tmax) {

    # Demography -----
    # Vaccinated class
    V <- V - rbinom(nlocs, size = V, prob = deaths) # die first
    waning <- rbinom(nlocs, size = V, prob = waning)
    V <- V - waning

    # Susceptible class
    S <- S - rbinom(nlocs, size = S, prob = deaths) + waning
    S <- S + rbinom(nlocs, size = S + V, prob = births)

    # components to add = puppy vaccination (in between campaigns)
    # and also immigration, colonization of empty patches?

    # Vaccination ----
    vacc_now <- vacc_dt[tstep == t]
    if (nrow(vacc_now) > 0) {
      nvacc <- sim_vacc(vacc_dt = vacc_now,
                        N = N[, t - 1], S, V, loc_id, prob_revacc, nlocs,
                        row_id,
                        additive, row_probs,
                        vacc_est)
    } else {
      nvacc <- 0
    }

    # balance vaccinated
    S <- S - nvacc
    V <- V + nvacc

    # Transmission ----
    # incursions (to do:pass a function)
    I_dt[max(I_dt$id) + 1:n_incs, ] <- sim_incursions(n_incs, cells_pop, cells_all,
                                             x_coord_pop,
                                             y_coord_pop, counter = max(I_dt$id),
                                             tstep,
                                             days_in_step = 7)


    # exposed -> infectious (those in tstep) (better way to do this?)
    I_now <- I_dt[floor(t_infectious) == t & contact == 1]
    secondaries <- secondary_fun(nrow(I_now), params)

    # infectious contacts (fix parameterization)
    exposed <- sim_bites(secondaries, ids, dispersal_fun,
                         counter = max(I_dt$id), res_m,
                         row_id, cell_id, t_infectious,
                         cells_pop, nrow, ncol, x_topl, y_topl, tstep = t,
                         sequential, allow_empties,
                         leave_district, max_tries)

    # were those contacts with a susceptible? (this updates exposed in the parent env!)
    sim_trans(exposed, S, N, nlocs, track)

    exposed$t_infectious <- t_infectious(t_infected = exposed$t_infected,
                                         days_in_step,
                                         generation_fun)

    # Update I_dt
    if(nrow(I_dt) > max(exposed$id)) {
      I_dt <- double_I(I_dt)
    }

    I_dt[exposed$id] <- exposed

    # Balance ----
    I <- tabulate(I_now$row_id, nbins = nlocs)
    I_loc <- tabulate(I_now[progen_id != 0]$row_id, nbins = nlocs)
    E <- E + tabulate(exposed[outcome == 1]$row_id, nbins = nlocs) - I_loc
    S <- S - E

    # Update matrices ----
    S_mat[, t] <- S
    E_mat[, t] <- E
    I_mat[, t] <- I
    V_mat[, V] <- V
    N_mat[, t] <- S + E + I_loc + V
  }

  # Reporting model (takes I_dt as an argument--modifies in parent env)

  # Summary functions which returns list of objs (or list of lists)
  out <- lapply(summary_funs, function(x) x())
  return(out)
}

