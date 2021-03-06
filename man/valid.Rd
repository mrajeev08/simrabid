% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/move_continuous.R
\name{valid}
\alias{valid}
\title{Accept a simulated movement}
\usage{
valid(leave_bounds, allow_invalid, outbounds, invalid)
}
\arguments{
\item{leave_bounds}{boolean, are movements to outside of the boundaries of the are being simulated are valid}

\item{allow_invalid}{boolean, are movements to empty patches (i.e. with no dogs or bodies of water, etc.) valid}

\item{outbounds}{whether cell id is out-of-bounds of the area being simulated}

\item{invalid}{whether cell id is considered an ivalid area (i.e. not populated or other barrier to movement)}
}
\value{
a boolean vector of length \code{outbounds}/\code{invalid} corresponding to whether
a movement is accepted as valid based on the arguments \code{leave_bounds} / \code{allow_invalid}
}
\description{
Helper function to determine if a simulated movement is valid.
To do: tests (same length as outbounds and boolean no NAs)
}
\keyword{internal}
\keyword{move}
