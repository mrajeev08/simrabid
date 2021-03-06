% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameterize.R
\name{beta_detect_monthly}
\alias{beta_detect_monthly}
\title{Observation function: with beta binomial distribution of monthly
detection probabilities}
\usage{
beta_detect_monthly(
  I_dt,
  params = list(detect_alpha = 80.1, detect_beta = 8.9)
)
}
\arguments{
\item{I_dt}{the line list of cases from \code{simrabid}}

\item{params}{list of parameters}
}
\value{
The reporting function should modify the line list
in place within the simrabid function (adding a \code{detected} column
as a boolean vector).
}
\description{
This function simulates the observation proccess drawing a detection
probability from a beta binomial distribution for each month.
Defaults are in \code{param_defaults}.
Can be passed to the \code{observe_fun} in \code{\link{simrabid}}.
}
