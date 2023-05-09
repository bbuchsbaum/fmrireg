


#' Simulate fMRI Time Series
#'
#' This function simulates an fMRI time series for multiple conditions with specified parameters.
#'
#' @param ncond The number of conditions to simulate.
#' @param nreps The number of repetitions per condition (default is 12).
#' @param amps A vector of amplitudes for each condition (default is a vector of 1s with length ncond).
#' @param isi A vector specifying the range of inter-stimulus intervals to sample from (default is c(3, 6)).
#' @param TR The repetition time of the fMRI acquisition (default is 1.5 seconds).
#'
#' @return A list with two elements:
#'   \itemize{
#'     \item onset: A vector of the onset times for each trial.
#'     \item mat: A matrix containing the simulated fMRI time series, with time in the first column and the simulated responses for each condition in the subsequent columns.
#'   }
#' @keywords internal
sim_ts <- function(ncond, nreps=12, amps=rep(1,ncond), isi=c(3,6), TR=1.5) {
  cond <- letters[1:ncond]
  trials <- sample(rep(cond, nreps))
  isis <- sample(isi[1]:isi[length(isi)], length(trials), replace=TRUE)
  onset <- cumsum(isis)
  
  time <- seq(0, max(onset+12), by=TR)
  ymat <- do.call(cbind, lapply(1:length(cond), function(i) {
    #print(i)
    idx <- which(trials == cond[i])
    reg <- regressor(onset[idx], amplitude=amps[i])
    evaluate(reg, time)
  }))
  
  list(onset=onset, mat=cbind(time, ymat))
}


#' Simulate fMRI Noise using ARMA Model
#'
#' This function simulates fMRI noise using an autoregressive moving average (ARMA) model.
#'
#' @param n The number of time points in the fMRI time series.
#' @param ar A numeric vector containing autoregressive (AR) coefficients (default is c(0.3)).
#' @param ma A numeric vector containing moving average (MA) coefficients (default is c(0.5)).
#' @param sd The standard deviation of the white noise component (default is 1).
#' @param seed An optional seed for reproducibility (default is NULL).
#'
#' @return A numeric vector containing the simulated fMRI noise.
#' @keywords internal
#' @export
simulate_fmri_noise <- function(n, ar = c(0.3), ma = c(0.5), sd = 1, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  noise <- arima.sim(n = n, model = list(ar = ar, ma = ma), sd = sd)
  return(noise)
}