% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/setup.R
\name{block_cells}
\alias{block_cells}
\title{Which cells should be considered invalid if movement occurs?}
\usage{
block_cells(rast, params = list(start_pop))
}
\arguments{
\item{rast}{raster generated using \code{setup_space} function}

\item{params}{a list of params, in this case,
integer vector of the population size in each cell of \code{rast}}
}
\value{
indexes to block, indexes considered out of bounds, and indexes to track
}
\description{
This is an example function for desginating certain cells invalid/outofbounds.
Can customize this to account for landscape features (i.e. lakes, rivers, mountains)
where movement is not possible. Function must have two arguments: the
raster input and other params passed as a list. Function must output the
indexes to block, the indexes considered out of bounds, and the indices to track.
}
