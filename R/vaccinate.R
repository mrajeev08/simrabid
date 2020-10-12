
#' Title
#'
#' @param loc_id
#' @param campaign_prob
#' @param coverage
#' @param sim_years
#' @param burn_in_years
#' @param steps_in_year
#' @param sample_tstep
#'
#' @return
#' @export
#'
#' @examples
#' # Simulate annual campaigns
# campaign prob can be either a single val or vector of vals corresponding to
# each loc_id
# cov function example
# coverage = function(n) rbeta(n, shape1 = 10, shape2 = 10)
# alternative sample_tstep (i.e. if you want campaigns to only happen in beginning of year)
# sample_tstep = function(n) sample.int(12, n, replace = true)

sim_campaigns <- function(loc_id, campaign_prob = 0.5,
                          coverage = 0.5,
                          sim_years = 10, burn_in_years = 5,
                          steps_in_year = 52,
                          sample_tstep = function(n) sample.int(52, n, replace = TRUE)) {

  vacc_dt <- data.table(loc_id = rep(loc_id, each = sim_years),
                        years = rep(1:sim_years, length(loc_id)))
  vacc_dt$vaccinated <- rbinom(nrow(vacc_dt), 1, campaign_prob)
  vacc_dt <- vacc_dt[vaccinated == 1]

  if(is.function(coverage)) {
    coverage <- coverage(nrow(vacc_dt))
  }

  vacc_dt[, c("tstep", "vacc") := .(
    sample_tstep(.N) + (years - 1)*steps_in_year + burn_in_years*steps_in_year,
    coverage
  )]

  return(vacc_dt[, c("loc_id", "tstep", "vacc")])
}

#' Simulate vaccination @ scale
#'
#' @param vacc_dt
#' @param N
#' @param S
#' @param V
#' @param loc_id
#' @param prob_revacc
#' @param row_id
#' @param additive
#' @param row_probs
#' @param allocate_by
#' @param vacc_est
#' @return
#' @export
#'
#' @examples
# # Simulate vaccination
# max_cells <- 1000
# vacc_dt <- data.table(loc_id = c(1, 27, 5),
#                       vacc = c(6000, 8000, 5000)) # fake vill vacc
# row_id <- sample.int(1000, max_cells, replace = FALSE)
# loc_id <- sample.int(75, max_cells, replace = TRUE) # village id for each cell
# S <- sample(1000:1500, max_cells, replace = TRUE)
# N <- sample(2000:3000, max_cells, replace = TRUE)
# E <- sample(300:400, max_cells, replace = TRUE)
# I <- sample(100:200, max_cells, replace = TRUE)
# V <- sample(500:1000, max_cells, replace = TRUE)
# prob <- runif(1000)
# prob_revacc = 0.5
# additive = TRUE
# skip control flows for fitting because otherwise it's slower
# Add a control flow argument in top level sim function
# vacc_dt needs to have loc_id, tstep, vacc as variables
# and vacc should be double and between [0, 1] or integer

# Cleaner internals?
# Limit objects created & ifelse statements?
# using row ids to match things up?
# Try matching vacc_new instead of doing vacc_new[1] (benchmark)
# Test & make sure it all works ok

sim_vacc <- function(vacc_dt, N, S, V, loc_id, prob_revacc, nlocs,
                     row_id, additive = TRUE, row_probs = NULL,
                     vacc_est = "coverage") {

  # make & join up data tables
  dem_now <- data.table(S, N, V, loc_id, row_id, row_probs)
  dem_now <- dem_now[loc_id %in% vacc_dt$loc_id]

  # summarize
  loc_dem <- dem_now[, lapply(.SD, sum), by = "loc_id",
                     .SDcols = c("S", "N", "V")]
  vacc_now <- loc_dem[vacc_dt, on = "loc_id"]

  if(vacc_est == "coverage") {
    vacc_now[, vacc := rbinom(length(vacc), size = N, prob = vacc)]
  }

  vacc_now[, new_vacc := vacc - rbinom(length(vacc),
                                        size = V, prob = prob_revacc)]

  # if no new vacc and additive campaigns are not possible
  # then the -new_vacc is the # of vaccinated dogs that went unallocated
  vacc_now[, c("unallocated",
               "new_vacc") := .(fifelse(new_vacc < 0 & additive, 0, -new_vacc),
                                fifelse(new_vacc < 0 & additive, vacc, 0))]

  # if more new vacc than total susceptibles in that vill
  # then the difference is added to unallocated
  # & we also constrain new vacc to be maximum the number of susceptibles available
  vacc_now[, c("unallocated",
               "new_vacc") := .(fifelse(new_vacc > S, unallocated + (new_vacc - S), unallocated),
                                fifelse(new_vacc > S, S, new_vacc))]

  # allocate vaccination (filter to current locations)
  sus_dt <- vacc_now[, c("loc_id", "new_vacc")][dem_now, on = "loc_id"]
  # replicate each row by the number of available susceptibles
  sus_dt <- sus_dt[rep(sus_dt[, .I], S)]

  # sample locations by the number of available susceptibles
  # and optionally weight that sampling by the probability of vacc in that cell
  vacc_ids <- sus_dt[, .(row_id = sample(row_id,
                                          size = new_vacc[1],
                                          replace = FALSE,
                                          prob = row_probs)),
                     by = "loc_id"]

  # update & return V
  V <- V + tabulate(vacc_ids$row_id, nbins = nlocs) # by row id

  return(list(V = V, vacc_now = vacc_now))
}


