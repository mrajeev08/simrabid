
simrabid <- function(start_up, start_vacc, I_seeds, vacc_dt,
                     params = c(list(R0 = 1.2, k = 1, iota = 4),
                                param_defaults),
                     days_in_step = 7,
                     observe_fun = beta_detect_monthly,
                     serial_fun = serial_lognorm,
                     dispersal_fun = steps_weibull,
                     secondary_fun = nbinom_constrained,
                     incursion_fun = sim_incursions_pois,
                     movement_fun = sim_movement_continuous,
                     sequential = TRUE, allow_invalid = TRUE,
                     leave_bounds = TRUE, max_tries = 100,
                     summary_funs = list(return_env = return_env),
                     track = TRUE,
                     prob_revacc = 0.5,
                     weights = NULL,
                     row_probs = NULL,
                     coverage = TRUE,
                     break_threshold = 0.8) {

  # pass the start_up objects into the function environment
  list2env(start_up, envir = environment())

  # initialize vaccination & infection
  list2env(init(start_pop, start_vacc, I_seeds, I_dt, cell_ids, admin_ids,
                row_ids, bins, nlocs, x_coord, y_coord, params, incursion_fun),
           envir = environment())

  list2env(vacc_dt, envir = environment()) # spit out data table into environment as vectors

  for (t in seq_len(tmax)) {

    # Demography -----
    # Vaccinated class
    V <- V - rbinom(nlocs, size = V, prob = death_prob) # die first
    waning <- rbinom(nlocs, size = V, prob = waning_prob)
    V <- V - waning

    # Susceptible class
    if (sum(is.na(S) | S < 0) > 0) browser()

    S <- S - rbinom(nlocs, size = S, prob = death_prob) + waning
    S <- S + rbinom(nlocs, size = S + V, prob = birth_prob)

    if (sum(is.na(S) | S < 0) > 0) browser()


    # Vaccination ----
    inds <- vacc_times == t

    if (sum(inds) > 0) {
      nvacc <- sim_vacc(vacc_est = vacc_est[inds],
                        vacc_locs = vacc_locs[inds],
                        S, V, N, loc_ids, bins,
                        row_ids, row_probs,
                        coverage)
    } else {
      nvacc <- 0
    }

    # balance vaccinated
    S <- S - nvacc
    V <- V + nvacc

    # Transmission ----
    # incursions
    cell_id_incs <- incursion_fun(cell_ids = cell_ids, params = params)

    if(sum(cell_id_incs) > 0) {
      I_incs <- add_incursions(cell_id_incs, cell_ids, admin_ids,
                               row_ids, x_coord, y_coord,
                               tstep = t, counter = max(I_dt$id),
                               days_in_step)

      I_dt[I_incs$id] <- I_incs

    }

    # exposed -> infectious (those in tstep) (better way to do this?)
    I_now <- I_dt[floor(t_infectious) == t & infected == TRUE]

    # Balance I & E class
    I <- tabulate(I_now$row_id, nbins = bins)
    I_loc <- tabulate(I_now$row_id[I_now$progen_id > 0], nbins = bins)
    E <- E - I_loc

    if(nrow(I_now) > 0) {
      secondaries <- secondary_fun(nrow(I_now), params)
    } else {
      secondaries <- 0
    }

    # infectious contacts
    if(sum(secondaries) > 0) {

      exposed <- sim_bites(secondaries, ids = I_now$id,
                           x_coords = I_now$x_coord, y_coords = I_now$y_coord,
                           t_infectious = I_now$t_infectious,
                           counter = max(I_dt$id),
                           sim_movement = movement_fun,
                           dispersal_fun, res_m,
                           row_ids, cell_ids,
                           cells_block, cells_out_bounds,
                           nrow, ncol, ncell,
                           x_topl, y_topl,
                           weights, admin_ids,
                           sequential, allow_invalid,
                           leave_bounds, max_tries,
                           params)

      # this should only be for ones that were successful (i.e. within & populated)
      exp_inds <- !exposed$invalid & !exposed$outbounds
      out <- sim_trans(row_id = exposed$row_id[exp_inds],
                       S, E, I, V, bins,
                       track)
      exposed$contact[exp_inds] <- out$contact
      exposed$infected[exp_inds] <- out$infected

      # were those contacts with a susceptible?  (better way to do this)
      exposed$t_infectious <- 0
      exposed$t_infectious[exposed$infected] <-
        t_infectious(n = length(exposed$t_infectious[exposed$infected]),
                     t_infected = exposed$t_infected[exposed$infected],
                     days_in_step, serial_fun, params)

      # Make sure colum order matches that of I_dt (better way to do this)
      setcolorder(exposed, c('id', 'cell_id', 'row_id', 'progen_id',
                             'path', 'x_coord', 'y_coord', 'invalid',
                             'outbounds', 't_infected', 'contact', 'infected',
                             't_infectious'))

      # if not tracking outcomes, then only track infected
      if(!track) {
        exposed <- exposed[infected == TRUE]
        if(nrow(exposed) > 0) {
          id_now <- max(I_dt$id)
          exposed$id <- id_now:(id_now + nrow(exposed) - 1)
        }
      }

      # Update I_dt
      if(nrow(exposed) > 0) {
        if(nrow(I_dt) < max(exposed$id)) {
          I_dt <- double_I(I_dt)
        }
      }

      I_dt[exposed$id] <- exposed

    } else {
      exposed <- empty_dt
    }

    # Final balance ----
    E_new <- tabulate(exposed$row_id[exposed$infected == TRUE], nbins = bins)
    E <- E + E_new
    if (sum(is.na(S) | S < E_new) > 0) browser()

    S <- S - E_new
    N <- S + E + I_loc + V

    # Update matrices ----
    S_mat[, t] <- S
    E_mat[, t] <- E
    I_mat[, t] <- I
    V_mat[, t] <- V

    prop_start_pop <- sum(N, na.rm = TRUE)/sum(start_pop, na.rm = TRUE)

    if (prop_start_pop < break_threshold) {
      break
    }

  }

  # Reporting model (adds a reported column to the I data.table) (modifies in place!)
  observe_fun(I_dt, params) # change this so it operates within data.table?

  # Summary functions which returns list of objs (or list of lists)
  out <- lapply(summary_funs, function(x) x())

  return(out)

}

