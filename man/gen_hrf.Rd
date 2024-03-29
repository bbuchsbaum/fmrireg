% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hrf.R
\name{gen_hrf}
\alias{gen_hrf}
\title{Construct an HRF Instance}
\usage{
gen_hrf(
  hrf,
  lag = 0,
  width = 0,
  precision = 0.1,
  summate = TRUE,
  normalize = FALSE,
  name = "gen_hrf",
  span = NULL,
  ...
)
}
\arguments{
\item{hrf}{A function mapping from time to signal.}

\item{lag}{Optional lag in seconds.}

\item{width}{Optional block width in seconds.}

\item{precision}{Sampling precision in seconds.}

\item{summate}{Whether to allow each impulse response function to "add" up (default: TRUE).}

\item{normalize}{Rescale so that the peak of the output is 1 (default: FALSE).}

\item{name}{The assigned name of the generated HRF.}

\item{span}{The span of the HRF (maximum width in seconds after which function reverts to zero).}

\item{...}{Extra parameters for the \code{hrf} function.}
}
\value{
An instance of type \code{HRF} inheriting from \code{function}.
}
\description{
\code{gen_hrf} takes a raw function \code{f(t)} and returns an HRF (Hemodynamic Response Function) instance.
}
\examples{

## Generate an HRF using SPMG1 canonical HRF, a lag of 3, and a width of 2.
grf <- gen_hrf(hrf_spmg1, lag=3, width=2)
grf(0:20)

hg <- purrr::partial(hrf_gaussian, mean=4, sd=1)
grf <- gen_hrf(hg, lag=1, width=2)

vals <- grf(0:20)
}
