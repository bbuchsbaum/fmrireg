% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/event_vector.R
\name{.checkEVArgs}
\alias{.checkEVArgs}
\title{Validate event arguments.}
\usage{
.checkEVArgs(name, vals, onsets, blockids, durations = NULL)
}
\arguments{
\item{name}{The name of the event.}

\item{vals}{The event values.}

\item{onsets}{Numeric vector of event onsets.}

\item{blockids}{Numeric vector of block IDs.}

\item{durations}{Numeric vector of durations (or a scalar).}
}
\value{
A list of validated event parameters.
}
\description{
Ensures that the event onsets, values, block IDs, and durations are
consistent. In particular, it asserts that onsets have no NA values and
are strictly increasing within each block.
}
\keyword{internal}
