% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hrf.R
\name{hrf_library}
\alias{hrf_library}
\alias{gen_hrf_library}
\title{Generate an HRF library from a parameter grid}
\usage{
hrf_library(fun, pgrid, ...)

gen_hrf_library(...)
}
\arguments{
\item{fun}{A function that generates an HRF, given a set of parameters.}

\item{pgrid}{A data frame where each row is a set of parameters.}

\item{...}{Additional arguments passed to \code{fun}.}
}
\value{
A combined HRF object representing the library.
}
\description{
\code{hrf_library} applies a base HRF generating function to each row of a parameter grid.
}
\keyword{internal}
