#' Simulate fMRI Time Series
#'
#' This function simulates an fMRI time series for multiple experimental conditions with specified parameters.
#' It generates a realistic event-related design with randomized inter-stimulus intervals and condition orders.
#'
#' @param ncond The number of conditions to simulate.
#' @param hrf The hemodynamic response function to use (default is fmrihrf::HRF_SPMG1).
#' @param nreps The number of repetitions per condition (default is 12).
#' @param amps A vector of amplitudes for each condition (default is a vector of 1s with length ncond).
#' @param ampsd The standard deviation of the amplitudes (default is 0).
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
#' @importFrom assertthat assert_that
#' @examples
#' # Simulate 3 conditions with different amplitudes
#' sim <- simulate_bold_signal(ncond = 3, amps = c(1, 1.5, 2), TR = 2)
#' 
#' # Plot the simulated time series
#' matplot(sim$mat[,1], sim$mat[,-1], type = "l", 
#'         xlab = "Time (s)", ylab = "BOLD Response")
#' 
#' @export
simulate_bold_signal <- function(ncond, hrf=fmrihrf::HRF_SPMG1, nreps=12, amps=rep(1,ncond), isi=c(3,6), ampsd=0, TR=1.5) {
  assert_that(ncond > 0, msg = "ncond must be positive")
  assert_that(length(amps) == ncond,
              msg = "Length of 'amps' must equal 'ncond'")
  # Note: If ampsd > 0, amplitude variability is sampled *once per condition*,
  # meaning all trials of a given condition share the same sampled amplitude.
  assert_that(length(isi) == 2 && isi[2] > isi[1], 
              msg = "ISI must be a vector of length 2 with isi[2] > isi[1]")
  assert_that(TR > 0, msg = "TR must be positive")
  
  # Use robust condition naming
  cond <- paste0("Cond", 1:ncond) 
  trials <- sample(rep(cond, nreps))
  isis <- runif(length(trials), min = isi[1], max = isi[2])
  onset <- cumsum(isis)
  
  span <- attr(hrf, "span") %||% 12
  time <- seq(0, max(onset) + span, by = TR)
  ymat <- do.call(cbind, lapply(1:length(cond), function(i) {
    idx <- which(trials == cond[i])
    reg <- fmrihrf::regressor(onset[idx], hrf, amplitude=rnorm(1, mean=amps[i], sd=ampsd))
    fmrihrf::evaluate(reg, time)
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
#' noise <- simulate_noise_vector(n_timepoints, TR = 2)
#' plot(noise, type = "l", xlab = "Time Point", ylab = "Signal")
#' 
#' @export
simulate_noise_vector <- function(n, TR = 1.5, ar = c(0.3), ma = c(0.5), sd = 1, 
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
#' @param hrf Hemodynamic response function to use (default is fmrihrf::HRF_SPMG1)
#' @param seed Optional seed for reproducibility (default is NULL)
#'
#' @return A list containing:
#'   \itemize{
#'     \item clean: The simulated signals without noise (from simulate_bold_signal)
#'     \item noisy: The signals with added noise
#'     \item noise: The simulated noise component
#'     \item onsets: Trial onset times
#'     \item conditions: Condition labels for each trial
#'   }
#'
#' @examples
#' # Simulate a dataset with 3 conditions
#' data <- simulate_simple_dataset(ncond = 3, TR = 2, snr = 0.5)
#' 
#' # Plot clean and noisy data
#' par(mfrow = c(2,1))
#' matplot(data$clean$mat[,1], data$clean$mat[,-1], type = "l",
#'         main = "Clean Signal", xlab = "Time (s)", ylab = "BOLD")
#' matplot(data$noisy[,1], data$noisy[,-1], type = "l",
#'         main = "Noisy Signal", xlab = "Time (s)", ylab = "BOLD")
#'
#' @export
simulate_simple_dataset <- function(ncond, nreps = 12, TR = 1.5, snr = 0.5, 
                                 hrf = fmrihrf::HRF_SPMG1, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  # Generate clean signals
  clean <- simulate_bold_signal(ncond = ncond, nreps = nreps, TR = TR, hrf = hrf)
  
  # Calculate noise level based on SNR
  signal_sd <- sd(as.vector(clean$mat[,-1]))
  noise_sd <- signal_sd / snr
  
  # Generate noise for each condition
  n_timepoints <- nrow(clean$mat)
  noise_mat <- replicate(ncond, 
                        simulate_noise_vector(n = n_timepoints, TR = TR, 
                                         sd = noise_sd))
  
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

#' Simulate fMRI Time Courses, Return Shared Onsets + Column-Specific Amplitudes/Durations
#'
#' Generates \eqn{n} time-series (columns) with a single set of onsets, but
#' *resampled* amplitudes/durations for each column if \code{amplitude_sd>0}
#' or \code{duration_sd>0}. Each column also gets independent noise. The result
#' is a list containing:
#' \itemize{
#'   \item \code{time_series}: a \code{matrix_dataset} with \eqn{T \times n}.
#'         The \code{event_table} uses the first column's amplitude/duration draws.
#'   \item \code{ampmat}: an \eqn{n\_events \times n} matrix of per-column amplitudes.
#'   \item \code{durmat}: an \eqn{n\_events \times n} matrix of per-column durations.
#'   \item \code{hrf_info}: info about the HRF.
#'   \item \code{noise_params}: info about noise generation (type + AR coefficients + SD).
#' }
#'
#' @details
#' - If \code{noise_type="ar1"} and you do not provide \code{noise_ar}, we
#'   default to \code{c(0.3)}.
#' - If \code{noise_type="ar2"} and you do not provide a 2-element \code{noise_ar},
#'   we default to \code{c(0.3, 0.2)}.
#' - Onsets are either provided or generated once for all columns.
#' - **Amplitudes/durations** are re-sampled \emph{inside the loop} so each
#'   column can differ randomly. The final arrays \code{ampmat} and \code{durmat}
#'   each have one column per time-series.
#' - The \code{matrix_dataset}'s \code{event_table} records the first column's
#'   amplitudes/durations. If you need each column's, see \code{ampmat} and
#'   \code{durmat}.
#'
#' @param n Number of time-series (columns).
#' @param total_time Numeric. Total scan length (seconds).
#' @param TR Numeric. Repetition time (seconds).
#' @param hrf Hemodynamic response function, e.g. \code{fmrihrf::HRF_SPMG1}.
#' @param n_events Number of events (ignored if \code{onsets} is provided).
#' @param onsets Optional numeric vector of event onsets. If \code{NULL}, will be generated.
#' @param isi_dist One of \code{"even"}, \code{"uniform"}, or \code{"exponential"}.
#'   Default is \code{"even"} so events are evenly spaced from 0..total_time.
#' @param isi_min,isi_max For \code{isi_dist="uniform"}.
#' @param isi_rate For \code{isi_dist="exponential"}.
#' @param durations Numeric, scalar or length-\code{n_events}. If \code{duration_sd>0},
#'   random sampling is done per column.
#' @param duration_sd Numeric. If >0, random variation in durations.
#' @param duration_dist \code{"lognormal"} or \code{"gamma"} (strictly positive).
#' @param amplitudes Numeric, scalar or length-\code{n_events}. If \code{amplitude_sd>0},
#'   random sampling is done per column.
#' @param amplitude_sd Numeric. If >0, random variation in amplitudes.
#' @param amplitude_dist \code{"lognormal"}, \code{"gamma"}, or \code{"gaussian"} (can be negative).
#' @param single_trial If TRUE, each event is a separate single-trial regressor that gets summed.
#' @param noise_type \code{"none"}, \code{"white"}, \code{"ar1"}, or \code{"ar2"}.
#' @param noise_ar Numeric vector for AR(1) or AR(2). If missing or insufficient,
#'   defaults are used (0.3 for AR(1); c(0.3,0.2) for AR(2)).
#' @param noise_sd Std dev of the noise.
#' @param random_seed Optional integer for reproducibility.
#' @param verbose If TRUE, prints messages.
#'
#' @return A list containing:
#' \describe{
#'   \item{\code{time_series}}{A \code{matrix_dataset} with \eqn{T \times n} data
#'         and \code{event_table} for the *first* column's random draws.}
#'   \item{\code{ampmat}}{An \eqn{n\_events \times n} numeric matrix of amplitudes.}
#'   \item{\code{durmat}}{An \eqn{n\_events \times n} numeric matrix of durations.}
#'   \item{\code{hrf_info}}{A list with HRF metadata.}
#'   \item{\code{noise_params}}{A list describing noise generation.}
#' }
#'
#' @importFrom stats rnorm rexp runif rgamma rlnorm arima.sim
#' @importFrom assertthat assert_that
#' @export

# Internal helper for value resampling
#' @param base Base values to resample from
#' @param sd Standard deviation for resampling
#' @param dist Distribution to use for resampling
#' @param allow_negative Whether to allow negative values
#' @keywords internal
.resample_param <- function(base, sd, dist = c("lognormal", "gamma", "gaussian"),
                            allow_negative = FALSE) {
  dist <- match.arg(dist)
  out  <- base
  if (sd > 0) {
    out <- vapply(seq_along(base), function(i) {
      mu <- base[i]
      if (!allow_negative && mu <= 0 && dist != "gaussian") {
        stop("base values must be >0 for lognormal or gamma sampling")
      }
      switch(dist,
             lognormal = rlnorm(1, meanlog = log(mu), sdlog = sd),
             gamma = {
               shape_par <- (mu^2) / (sd^2)
               rate_par  <- mu / (sd^2)
               rgamma(1, shape = shape_par, rate = rate_par)
             },
             gaussian = rnorm(1, mean = mu, sd = sd))
    }, numeric(1))
  }
  out
}

simulate_fmri_matrix <- function(
    n                 = 1,
    total_time        = 240,
    TR                = 2,
    hrf               = fmrihrf::HRF_SPMG1,
    
    n_events          = 10,
    onsets            = NULL,
    isi_dist          = c("even", "uniform", "exponential"),
    isi_min           = 2,
    isi_max           = 6,
    isi_rate          = 0.25,
    
    durations         = 0,
    duration_sd       = 0,
    duration_dist     = c("lognormal", "gamma"),
    
    amplitudes        = 1,
    amplitude_sd      = 0,
    amplitude_dist    = c("lognormal", "gamma", "gaussian"),
    
    single_trial      = FALSE,
    noise_type        = c("none", "white", "ar1", "ar2"),
    noise_ar          = NULL,
    noise_sd          = 1.0,
    
    random_seed       = NULL,
    verbose           = FALSE,
    buffer            = 16
) {
  # ---------------------------
  # 0) Setup
  # ---------------------------
  # Parameter validation
  assert_that(n > 0, msg = "n must be positive")
  assert_that(total_time > 0, msg = "total_time must be positive")
  assert_that(TR > 0, msg = "TR must be positive")
  assert_that(n_events > 0, msg = "n_events must be positive")
  assert_that(isi_max > isi_min, msg = "isi_max must be greater than isi_min")
  
  if (!is.null(random_seed)) {
    set.seed(random_seed)
  }
  
  isi_dist        <- match.arg(isi_dist)
  noise_type      <- match.arg(noise_type)
  duration_dist   <- match.arg(duration_dist)
  amplitude_dist  <- match.arg(amplitude_dist)
  
  # Handle default AR params if user omitted or gave insufficient length
  if (noise_type == "ar1") {
    if (is.null(noise_ar) || length(noise_ar) < 1) {
      noise_ar <- 0.3
      if (verbose) message("Defaulting to noise_ar = 0.3 for AR(1).")
    }
  } else if (noise_type == "ar2") {
    if (is.null(noise_ar) || length(noise_ar) < 2) {
      noise_ar <- c(0.3, 0.2)
      if (verbose) message("Defaulting to noise_ar = c(0.3, 0.2) for AR(2).")
    }
  }
  
  # Original time grid without buffer (for event generation)
  effective_time <- total_time - buffer
  
  # Sample event onsets - ensure they only occur in the effective time, not the buffer
  if (isi_dist == "uniform") {
    isi_samples <- runif(n_events, min = isi_min, max = isi_max)
  } else if (isi_dist == "exponential") {
    isi_samples <- isi_min + rexp(n_events, rate = isi_rate)
  } else if (isi_dist == "even") {
    isi_samples <- rep(effective_time / n_events, n_events)
  }
  
  # Generate onsets cumulative ISIs, but ensure they fall within effective_time
  onsets <- cumsum(isi_samples)
  if (max(onsets) > effective_time) {
    # Keep only events that fit within effective time
    keepers <- which(onsets <= effective_time)
    onsets <- onsets[keepers]
    n_events <- length(onsets)
    message(sprintf("Reduced to %d events to fit within effective time", n_events))
  }
  
  # Final time grid with buffer
  n_time_points <- ceiling(total_time / TR)
  time_grid <- seq(0, by = TR, length.out = n_time_points)
  
  # ---------------------------
  # 1) Helper fns to sample durations/amplitudes for 1 column
  # ---------------------------
  do_sample_durations <- function() {
    if (length(durations) == 1L) {
      base_durs <- rep(durations, n_events)
    } else {
      if (length(durations) != n_events) {
        stop("durations must be length=1 or match n_events.")
      }
      base_durs <- durations
    }
    .resample_param(base_durs, duration_sd, duration_dist)
  }

  do_sample_amplitudes <- function() {
    if (length(amplitudes) == 1L) {
      base_amps <- rep(amplitudes, n_events)
    } else {
      if (length(amplitudes) != n_events) {
        stop("amplitudes must be length=1 or match n_events.")
      }
      base_amps <- amplitudes
    }
    .resample_param(base_amps, amplitude_sd, amplitude_dist,
                    allow_negative = (amplitude_dist == "gaussian"))
  }

  gen_noise <- function() {
    switch(noise_type,
           none  = rep(0, n_time_points),
           white = rnorm(n_time_points, 0, noise_sd),
           ar1   = arima.sim(model = list(ar = noise_ar),
                             n = n_time_points, sd = noise_sd),
           ar2   = arima.sim(model = list(ar = noise_ar),
                             n = n_time_points, sd = noise_sd))
  }
  
  # ---------------------------
  # 2) Build data matrix, col by col
  #    plus store amplitude/duration in big matrices
  # ---------------------------
  signal_list <- vector("list", n)
  ampmat <- matrix(NA, nrow=n_events, ncol=n)  # amplitude
  durmat <- matrix(NA, nrow=n_events, ncol=n)  # durations
  
  for (ii in seq_len(n)) {
    # sample for column ii
    this_dur <- do_sample_durations()
    this_amp <- do_sample_amplitudes()
    
    durmat[, ii] <- this_dur
    ampmat[, ii] <- this_amp
    
    # Build the BOLD signal using the full time grid (including buffer)
    if (!single_trial) {
      reg <- fmrihrf::regressor(
        onsets    = onsets,
        hrf       = hrf,
        duration  = this_dur,
        amplitude = this_amp
      )
      bold_signal <- fmrihrf::evaluate(reg, grid=time_grid)
    } else {
      bold_signal <- numeric(n_time_points)
      for (j in seq_along(onsets)) {
        sreg <- fmrihrf::single_trial_regressor(
          onsets    = onsets[j],
          hrf       = hrf,
          duration  = this_dur[j],
          amplitude = this_amp[j]
        )
        bold_signal <- bold_signal + fmrihrf::evaluate(sreg, time_grid)
      }
    }
    
    # Add noise to the entire signal including buffer
    eps <- gen_noise()
    noisy_tc <- bold_signal + eps
    signal_list[[ii]] <- noisy_tc
  }
  
  sim_matrix <- do.call(cbind, signal_list)
  
  # ---------------------------
  # 3) matrix_dataset
  #    event_table = first column's durations/amplitudes
  # ---------------------------
  event_tab <- data.frame(
    run       = 1,
    onset     = onsets,
    duration  = durmat[,1],
    amplitude = ampmat[,1]
  )
  
  ds <- fmridataset::matrix_dataset(
    datamat     = sim_matrix,
    TR          = TR,
    run_length  = n_time_points,
    event_table = event_tab
  )
  
  # ---------------------------
  # 4) Return
  # ---------------------------
  out <- list(
    time_series  = ds,         # matrix_dataset T x n
    ampmat       = ampmat,     # n_events x n
    durmat       = durmat,     # n_events x n
    hrf_info     = list(
      hrf_class = class(hrf),
      hrf_name  = attr(hrf, "name"),
      nbasis    = fmrihrf::nbasis(hrf),
      span      = attr(hrf, "span")
    ),
    noise_params = list(
      noise_type = noise_type,
      noise_ar   = noise_ar,
      noise_sd   = noise_sd
    )
  )
  
  return(out)
}
