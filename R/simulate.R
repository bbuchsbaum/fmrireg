#' Simulate fMRI Time Series
#'
#' This function simulates an fMRI time series for multiple experimental conditions with specified parameters.
#' It generates a realistic event-related design with randomized inter-stimulus intervals and condition orders.
#'
#' @param ncond The number of conditions to simulate.
#' @param hrf The hemodynamic response function to use (default is HRF_SPMG1).
#' @param nreps The number of repetitions per condition (default is 12).
#' @param amps A vector of amplitudes for each condition (default is a vector of 1s with length ncond).
#' @param apmsd the standard deviation of the amplitudes (default is 0).
#' @param isi A vector of length 2 specifying the range of inter-stimulus intervals to sample from (default is c(3, 6) seconds).
#' @param TR The repetition time of the fMRI acquisition (default is 1.5 seconds).
#'
#' @return A list with the following components:
#'   \itemize{
#'     \item onset: A vector of the onset times for each trial
#'     \item condition: A vector of condition labels for each trial
#'     \item mat: A matrix containing the simulated fMRI time series:
#'       \itemize{
#'         \item Column 1: Time points (in seconds)
#'         \item Columns 2:(ncond+1): Simulated BOLD responses for each condition
#'       }
#'   }
#' 
#' @examples
#' # Simulate 3 conditions with different amplitudes
#' sim <- sim_ts(ncond = 3, amps = c(1, 1.5, 2), TR = 2)
#' 
#' # Plot the simulated time series
#' matplot(sim$mat[,1], sim$mat[,-1], type = "l", 
#'         xlab = "Time (s)", ylab = "BOLD Response")
#' 
#' @export
sim_ts <- function(ncond, hrf=HRF_SPMG1, nreps=12, amps=rep(1,ncond), isi=c(3,6), ampsd=0, TR=1.5) {
  assert_that(length(amps) == ncond, 
              msg = "Length of amplitudes vector must match number of total trials (ncond")
  assert_that(length(isi) == 2 && isi[2] > isi[1], 
              msg = "ISI must be a vector of length 2 with isi[2] > isi[1]")
  assert_that(TR > 0, msg = "TR must be positive")
  
  cond <- letters[1:ncond]
  trials <- sample(rep(cond, nreps))
  isis <- sample(isi[1]:isi[2], length(trials), replace=TRUE)
  onset <- cumsum(isis)
  
  time <- seq(0, max(onset+12), by=TR)
  ymat <- do.call(cbind, lapply(1:length(cond), function(i) {
    idx <- which(trials == cond[i])
    reg <- regressor(onset[idx], hrf, amplitude=rnorm(1, mean=amps[i], sd=ampsd))
    evaluate(reg, time)
  }))
  
  list(onset=onset, condition=trials, mat=cbind(time, ymat))
}

#' Simulate fMRI Noise
#'
#' This function simulates realistic fMRI noise by combining:
#' \itemize{
#'   \item Temporal autocorrelation using an ARMA model
#'   \item Low-frequency drift
#'   \item Physiological noise (cardiac and respiratory)
#' }
#'
#' @param n The number of time points in the fMRI time series
#' @param TR The repetition time in seconds (default is 1.5)
#' @param ar A numeric vector containing autoregressive (AR) coefficients (default is c(0.3))
#' @param ma A numeric vector containing moving average (MA) coefficients (default is c(0.5))
#' @param sd The standard deviation of the white noise component (default is 1)
#' @param drift_freq Frequency of the low-frequency drift in Hz (default is 1/128)
#' @param drift_amplitude Amplitude of the low-frequency drift (default is 2)
#' @param physio Logical; whether to add simulated physiological noise (default is TRUE)
#' @param seed An optional seed for reproducibility (default is NULL)
#'
#' @return A numeric vector containing the simulated fMRI noise
#' 
#' @examples
#' # Simulate noise for a 5-minute scan with TR=2s
#' n_timepoints <- 150  # 5 minutes * 60 seconds / 2s TR
#' noise <- simulate_fmri_noise(n_timepoints, TR = 2)
#' plot(noise, type = "l", xlab = "Time Point", ylab = "Signal")
#' 
#' @export
simulate_fmri_noise <- function(n, TR = 1.5, ar = c(0.3), ma = c(0.5), sd = 1, 
                               drift_freq = 1/128, drift_amplitude = 2,
                               physio = TRUE, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Generate ARMA noise
  noise <- arima.sim(n = n, model = list(ar = ar, ma = ma), sd = sd)
  
  # Add low-frequency drift
  time <- seq(0, (n-1)*TR, by = TR)
  drift <- drift_amplitude * sin(2 * pi * drift_freq * time)
  noise <- noise + drift
  
  # Add physiological noise if requested
  if (physio) {
    # Simulate cardiac (~1.2 Hz) and respiratory (~0.3 Hz) noise
    cardiac <- 0.5 * sin(2 * pi * 1.2 * time)
    respiratory <- 0.8 * sin(2 * pi * 0.3 * time)
    noise <- noise + cardiac + respiratory
  }
  
  return(noise)
}

#' Simulate Complete fMRI Dataset
#'
#' This function simulates a complete fMRI dataset by combining task-related signals
#' with realistic noise. It returns both the clean signals and the noisy data.
#'
#' @param ncond Number of conditions to simulate
#' @param nreps Number of repetitions per condition (default is 12)
#' @param TR Repetition time in seconds (default is 1.5)
#' @param snr Signal-to-noise ratio (default is 0.5)
#' @param hrf Hemodynamic response function to use (default is HRF_SPMG1)
#' @param seed Optional seed for reproducibility (default is NULL)
#'
#' @return A list containing:
#'   \itemize{
#'     \item clean: The simulated signals without noise (from sim_ts)
#'     \item noisy: The signals with added noise
#'     \item noise: The simulated noise component
#'     \item onsets: Trial onset times
#'     \item conditions: Condition labels for each trial
#'   }
#'
#' @examples
#' # Simulate a dataset with 3 conditions
#' data <- simulate_fmri_dataset(ncond = 3, TR = 2, snr = 0.5)
#' 
#' # Plot clean and noisy data
#' par(mfrow = c(2,1))
#' matplot(data$clean$mat[,1], data$clean$mat[,-1], type = "l",
#'         main = "Clean Signal", xlab = "Time (s)", ylab = "BOLD")
#' matplot(data$noisy[,1], data$noisy[,-1], type = "l",
#'         main = "Noisy Signal", xlab = "Time (s)", ylab = "BOLD")
#'
#' @export
simulate_fmri_dataset <- function(ncond, nreps = 12, TR = 1.5, snr = 0.5, 
                                 hrf = HRF_SPMG1, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Generate clean signals
  clean <- sim_ts(ncond = ncond, nreps = nreps, TR = TR, hrf = hrf)
  
  # Calculate noise level based on SNR
  signal_sd <- sd(as.vector(clean$mat[,-1]))
  noise_sd <- signal_sd / snr
  
  # Generate noise for each condition
  n_timepoints <- nrow(clean$mat)
  noise_mat <- replicate(ncond, 
                        simulate_fmri_noise(n = n_timepoints, TR = TR, 
                                         sd = noise_sd, seed = seed))
  
  # Combine signal and noise
  noisy_mat <- cbind(clean$mat[,1], clean$mat[,-1] + noise_mat)
  
  list(
    clean = clean,
    noisy = noisy_mat,
    noise = noise_mat,
    onsets = clean$onset,
    conditions = clean$condition
  )
}