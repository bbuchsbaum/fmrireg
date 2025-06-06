% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hrf-functions.R
\name{hrf_toeplitz}
\alias{hrf_toeplitz}
\title{HRF Toeplitz Matrix}
\usage{
hrf_toeplitz(hrf, time, len, sparse = FALSE)
}
\arguments{
\item{hrf}{The hemodynamic response function.}

\item{time}{A numeric vector representing the time points.}

\item{len}{The length of the output Toeplitz matrix.}

\item{sparse}{Logical, if TRUE, the output Toeplitz matrix is returned as a sparse matrix (default: FALSE).}
}
\value{
A Toeplitz matrix for HRF convolution.
}
\description{
Create a Toeplitz matrix for hemodynamic response function (HRF) convolution.
}
