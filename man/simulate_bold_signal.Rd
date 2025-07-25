% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{simulate_bold_signal}
\alias{simulate_bold_signal}
\title{Simulate fMRI Time Series}
\usage{
simulate_bold_signal(
  ncond,
  hrf = fmrihrf::HRF_SPMG1,
  nreps = 12,
  amps = rep(1, ncond),
  isi = c(3, 6),
  ampsd = 0,
  TR = 1.5
)
}
\arguments{
\item{ncond}{The number of conditions to simulate.}

\item{hrf}{The hemodynamic response function to use (default is fmrihrf::HRF_SPMG1).}

\item{nreps}{The number of repetitions per condition (default is 12).}

\item{amps}{A vector of amplitudes for each condition (default is a vector of 1s with length ncond).}

\item{isi}{A vector of length 2 specifying the range of inter-stimulus intervals to sample from (default is c(3, 6) seconds).}

\item{ampsd}{The standard deviation of the amplitudes (default is 0).}

\item{TR}{The repetition time of the fMRI acquisition (default is 1.5 seconds).}
}
\value{
A list with the following components:
\itemize{
\item onset: A vector of the onset times for each trial
\item condition: A vector of condition labels for each trial
\item mat: A matrix containing the simulated fMRI time series:
\itemize{
\item Column 1: Time points (in seconds)
\item Columns 2:(ncond+1): Simulated BOLD responses for each condition
}
}
}
\description{
This function simulates an fMRI time series for multiple experimental conditions with specified parameters.
It generates a realistic event-related design with randomized inter-stimulus intervals and condition orders.
}
\examples{
# Simulate 3 conditions with different amplitudes
sim <- simulate_bold_signal(ncond = 3, amps = c(1, 1.5, 2), TR = 2)

# Plot the simulated time series
matplot(sim$mat[,1], sim$mat[,-1], type = "l", 
        xlab = "Time (s)", ylab = "BOLD Response")

}
