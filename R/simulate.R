
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
                     weights = NULL,
                     row_probs = NULL,
                     coverage = FALSE,
                     break_threshold = 0.8,
                     by_admin = FALSE) {

  # pass the start_up objects into the function environment
  list2env(start_up, envir = environment())

  # initialize vaccination & infection
  list2env(init(start_pop, start_vacc, I_seeds, I_dt, cell_ids, admin_ids,
                row_ids, bins, nlocs, x_coord, y_coord, params, incursion_fun,
                ncells),
           envir = environment())

  # spit out data table into environment as vectors
  list2env(vacc_dt, envir = environment())

  if(!by_admin) admin_ids <- NULL
  if(by_admin) loc_ids <- row_ids

  for (t in seq_len(tmax)) {

    # Demography -----

    # Vaccinated class
    V <- V - rbinom(bins, size = V, prob = death_prob) # die first
    waning <- rbinom(bins, size = V, prob = waning_prob)
    V <- V - waning

    # Susceptible class
    S <- S - rbinom(bins, size = S, prob = death_prob) + waning
    S <- S + rbinom(bins, size = S + V, prob = birth_prob)

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
      I_incs <- add_incursions(cell_id_incs, cell_ids, ncells, admin_ids,
                               row_ids, x_coord, y_coord,
                               tstep = t, counter = max(I_dt$id),
                               days_in_step)

      if(nrow(I_dt) <= max(I_incs$id)) {
          I_dt <- double_I(I_dt)
      }

      I_dt[I_incs$id] <- I_incs

    }

    # exposed -> infectious (those in tstep)
    I_now <- I_dt[floor(t_infectious) == t & infected]

    # Balance I & E class
    I <- tabulate(I_now$row_id, nbins = bins)
    # only local cases included
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
                           nrows, ncols, ncells,
                           x_topl, y_topl,
                           weights, admin_ids,
                           sequential, allow_invalid,
                           leave_bounds, max_tries,
                           params)

      # Currently handling this in the simulate function (move to sim_trans)
      exp_inds <- which(!exposed$invalid & !exposed$outbounds & !is.na(exposed$row_id))

      if(length(exp_inds) > 0) {

        out <- sim_trans(row_id = exposed$row_id[exp_inds],
                         S, E, I, V, bins, track)

        exposed$contact[exp_inds] <- out$contact
        exposed$infected[exp_inds] <- out$infected

        exposed$t_infectious <-
          t_infectious(n = nrow(exposed),
                       t_infected = exposed$t_infected,
                       days_in_step, serial_fun, params)

      } else {
        exposed <- empty_dt
      }


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
        if(nrow(I_dt) <= max(exposed$id)) {
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

    S <- S - E_new

    N <- S + E + I_loc + V

    # Update matrices ----
    S_mat[, t] <- S
    E_mat[, t] <- E
    I_mat[, t] <- I
    V_mat[, t] <- V
    N_mat[, t] <- N

    prop_start_pop <- sum(N, na.rm = TRUE)/sum(start_pop, na.rm = TRUE)

    if (prop_start_pop < break_threshold) {
      break
    }

  }

  # Filter to incursions + cases seeded locally
  I_dt <- I_dt[id != 0]

  # Observation model (adds a detected column to the I data.table)
  # (modifies I_dt in place)
  observe_fun(I_dt, params)

  # Summary functions which returns list of outputs
  out <- lapply(summary_funs, function(x) x())

  return(out)

}

