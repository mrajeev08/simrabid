# Summary function (examples)

# simplest function finds objects in parent environment and returns them in list
# I_dt <- "This is I_dt"
# I <- list(this_list = "This is I", this_too = "This is I, too, a list of lists")
# S <- data.table(x = "This is S")
# # use mget is the basic version, which will return you the objects from the simulation
# # function
# x <- 500

use_mget <- function(names = c("I_dt", "I")) {
  mget(names, envir = parent.frame(3))
}

# returns everything in environment in a named list
return_env <- function() {
  out <- as.list.environment(parent.frame(3))
  out$summary_funs <- NULL
  return(out)
}

# example
# summary_funs <- list(use_mget = use_mget, return_env = return_env)
# test <- function(summary_funs, I_dt = "This is I_dt", I = "This is I",
#                  S = "This is S") {
#   lapply(summary_funs, function(x) x())
# }
#
# test(summary_funs)

# parent.frame 3 means you have to step out of list of functions you're calling it from
# and the function itself
# it won't work if you call the function itself from inside a function
# another_mget <- function(names = "x") {
#   mget(names, envir = parent.frame(3))
# }


# will this work then?
# test <- function(summary_funs, I_dt, I, S, x) {
#   # this is an x not in the global
#   x <- "This is not global x please!"
#   lapply(summary_funs, function(x) x())
# }
#
# # This doesn't work! Has to be passed in a list for the environment to be the right one
# test <- function(use_mget, I_dt, I, S, x) {
#   # this is an x not in the global
#   I <- "This is not global I please!"
#   use_mget()
# }
