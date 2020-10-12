# Getting probability from rate ------------------------------------------------

# ## example 1: turn annual waning rate of 0.33 to weekly prob
# get.prob(rate = 0.33, step = 52)
# ## example 2: turn annual birth rate of 0.45 to monthly prob
# get.prob(rate = 0.45, step = 12)

get.prob <- function(rate, step) {
  converted <- (1 + rate)^(1/step) - 1
  return(1 - exp(-converted))
}

# get row id from cell id
# (in simulation so that you don't have to have a super large matrix with empty cells!)
rowfromcell <- function(cell_id, cells, rows) {
  rows[cells %in% cell_id]
}
