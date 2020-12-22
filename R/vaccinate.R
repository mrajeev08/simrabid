#' Simulated vaccination campaigns
#'
#' \code{sim_campaigns} simulates annual vaccination campaigns across locations & time.
#'
#' This function simulates annual vaccination campaigns for a givens et of locations (`locs`)
#' over a certain time frame. Optionally, you can include a burn-in period where
#' no campaigns occur (in addition to the imeframe you specify).
#' You can flexibly pass simulated coverage and the
#' timestep which campaigns occur using custom functions which accept a single
#' parameter (n = the number of values to return). Defaults to using a weekly
#' time step, but this can be changed through the steps_in_year
#' parameter. See examples for more details.
#'
#' To Do:
#' - Make it possible to pass location specific times & cov?
#' - Do we need to be flexible with steps_in_year?
#' - Write catches if sample_tstep & steps_in_year are off
#' - And if locs & probs don't match in length
#' - And if functions passed don't have n as an argument (can they use function env vars though?)
#' - Tests: return should be a data.table with three columns; timestep should be in realm of possible &
#' burn in period should be zero. vacc_cov should be between [0, 1], vacc_locs should be all from locs only; no NAs
#'
#' @param locs numeric or character vector of ids of locations
#' @param campaign_prob either a single numeric probability [0, 1] or vector
#'  of probabilities of the same length as locs which specifies the probability
#'  that a campaign occurs in a given location in a given year
#' @param coverage either a function that returns a coverage estimate [0, 1],
#'  a single numeric value of coverage
#' @param sim_years the number of years to simulate campaigns over
#' @param burn_in_years the number of years to start without any vaccination
#' @param steps_in_year timestep to allocate campaigns across, defaults to weekly
#' @param sample_tstep function to sample the timestep during which a campaign occurs in
#'   a given location
#'
#' @return a data.table with three columns: vacc_times (timestep of campaign),
#'   vacc_cov (coverage estimate), vacc_locs (location of campaign)
#'
#' @export
#' @import data.table
#' @keywords vaccinate
#'
#' @examples
#' ex_campaign <- sim_campaigns(locs = 1:10) # with defaults
#' ex_campaign <- sim_campaigns(locs = 1:10, burn_in_years = 10) # longer burn in
#'
#' # with a beta distribution around coverage (takes n as only parameter)
#' cov_fun <- function(n) rbeta(n, shape1 = 2, shape2 = 2)
#' ex_campaign <- sim_campaigns(locs = 1:10, coverage = cov_fun)
#'
#' # with a different sample_tstep function (takes n as only parameter)
#' only half of year do campaigns occur
#' sample_half <- function(n) sample.int(26, n, replace = TRUE)
#' ex_campaign <- sim_campaigns(locs = 1:10, sample_tstep = sample_half)
#' hist(ex_campaign$vacc_times)
#'
#' # three month gaps between campaigns
#' sample_threemos <- function(n) sample(seq(0, 52, 12), n, replace = TRUE)
#' ex_campaign <- sim_campaigns(locs = 1:10, sample_tstep = sample_threemos)
#' hist(ex_campaign$vacc_times)
#'
#' # simulate on monthly timestep
#' ex_campaign <- sim_campaigns(locs = 1:10, steps_in_year = 12,
#' sample_tstep = function(n) {sample.int(12, n, replace = TRUE)})
#'
sim_campaigns <- function(locs, campaign_prob = 0.5,
                          coverage = 0.5,
                          sim_years = 10, burn_in_years = 5,
                          steps_in_year = 52,
                          sample_tstep = function(n) {
                            sample.int(52, n, replace = TRUE)
                          }) {
  vacc_dt <- data.table(
    vacc_locs = rep(locs, each = sim_years),
    years = rep(1:sim_years, length(locs)),
  )

  if(length(campaign_prob) > 1) {
    vacc_dt$campaign_prob <- rep(campaign_prob, each = sim_years) # vector matched to locs
  }

  vacc_dt[, vaccinated := rbinom(.N, 1, campaign_prob)] # probability that location had a campaign
  vacc_dt <- vacc_dt[vaccinated == TRUE]

  if (is.function(coverage)) {
    vacc_cov <- coverage(nrow(vacc_dt)) # coverage acheived during campaign
  }

  vacc_dt[, c("vacc_times", "vacc_cov") := .(
    sample_tstep(.N) + (years - 1) * steps_in_year + burn_in_years * steps_in_year,
    coverage
  )]

  return(vacc_dt[, c("vacc_times", "vacc_cov", "vacc_locs")])
}

#' Simulate vaccination at individual level
#'
#' \code{sim_vacc} simulates the number of vaccinated individuals at scale
#'
#' Translates coverage estimates from \code{\link{sim_campaigns}}
#' or from user defined data.table of location specific vaccination coverage/numbers
#' filtered to the current timestep into numbers of vaccinated individuals in each
#' grid cell at the scale specified in the simulation. Assumes that all currently vaccinated
#' individuals are revaccinated (so campaigns are not additive). If scale of vaccination is at a coarser
#' resolution (i.e. at village level) allocates to grid cell based on number of susceptibles available.
#' This can also be modified by passing `row_probs`, grid cell level probabilities of vaccination (i.e. based on household surveys or spatial covariates).
#'
#' To Do:
#' - Tests/ catches for lengths? This might slow things down! (maybe should be higher level test/catch)
#' - Tests: return should be same length as nlocs, should be positive integer & non NA
#'
#' @param vacc_cov numeric vector of coverage estimates (either [0, 1] or the number vaccinated)
#' @param vacc_locs vector of numeric or character ids of same length as `vacc_cov`,
#'   the locations corresponding to the coverage estimates
#' @param S,V,N, numeric vector of Susceptibles, Vaccinated, and the Population size
#'   respectively, of same length
#' @param loc_ids vector of numeric or character ids of same length as **S/V/N** that
#'   matches the demographic vectors to the locations from `vacc_locs`
#' @param nlocs numeric, the total number of grid cells being simulated (i.e. at the scale of the simulation,
#'   this can be different than the scale of the vaccination campaigns)
#' @param row_ids numeric vector of same length as **S/V/N** of the row id corresponding to the demographic vectors to
#'   sample which grid cell vaccinations are allocated to
#' @param row_probs numeric vector, optionally can pass the probability of vaccination
#'   for each grid cell, defaults to NULL.
#' @param coverage boolean, defaults to TRUE, if TRUE then coverage estimates [0, 1]
#'   are translated to numbers with rbinom and then allocated, if FALSE then number
#'   vaccinated are allocated (constrained to the available susceptible dogs)
#'
#' @return a numeric vector of length **nlocs** of the number of currently vaccinated
#' individuals in each grid cell
#'
#' @import data.table
#' @keywords internal
#'
sim_vacc <- function(vacc_cov, vacc_locs, S, V, N, loc_ids, nlocs,
                     row_ids, row_probs = NULL, coverage = TRUE) {

  # make & join up data tables
  dem_now <- data.table(S, V, N, loc_ids, row_ids, row_probs)
  dem_now <- dem_now[loc_ids %in% vacc_locs]

  # summarize
  vacc_now <- dem_now[, lapply(.SD, sum),
    by = "loc_ids",
    .SDcols = c("S", "V", "N")
  ]
  vacc_now[, vacc := vacc_cov[match(loc_ids, vacc_locs)]]

  if (coverage) {
    vacc_now[, vacc := rbinom(length(vacc), size = N, prob = vacc)]
  }

  vacc_now[, nvacc := vacc - V]

  # vacc only new individuals (i.e. vacc)
  # & we also constrain vacc to be maximum the number of susceptibles available
  vacc_now[, nvacc := fcase(
    nvacc < 0, 0L,
    nvacc < S & nvacc > 0, as.integer(nvacc),
    nvacc > S, as.integer(S)
  )]

  # allocate vaccination (filter to current locations)
  sus_dt <- vacc_now[, c("loc_ids", "nvacc")][dem_now, on = "loc_ids"]
  # replicate each row by the number of available susceptibles
  sus_dt <- sus_dt[rep(sus_dt[, .I], S)]

  # sample locations by the number of available susceptibles
  # and optionally weight that sampling by the probability of vacc in that cell
  vacc_ids <- sus_dt[, .(row_id = sample(row_ids,
    size = nvacc[1],
    replace = FALSE,
    prob = row_probs
  )),
  by = "loc_ids"
  ]

  # update & return V
  nvacc <- tabulate(vacc_ids$row_id, nbins = nlocs) # by row id

  return(nvacc)
}
