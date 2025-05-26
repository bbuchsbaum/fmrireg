###############################################################################
# Regressor Utilities and Constructors for fMRI Modeling
###############################################################################

# Load required packages (make sure these are declared in DESCRIPTION)
library(assertthat)
library(RANN)
library(pracma)
library(memoise)

###############################################################################
# Helper Functions
###############################################################################

#' Ensure a vector is of length 1 or matches a supplied length.
#'
#' @param val The input vector.
#' @param len The target length.
#' @return A vector of length \code{len}.
#' @keywords internal
conform_len <- function(val, len) {
  name <- deparse(substitute(val))
  if (length(val) == 1) {
    rep(val, len)
  } else if (length(val) == len) {
    val
  } else {
    stop(paste(name, "must be of length 1 or the same length as onsets"))
  }
}

#' Capture additional arguments as a list.
#' @keywords internal
dots <- function(...) {
  eval(substitute(alist(...)))
}

###############################################################################
# RCPP Wrappers (Note: Rcpp functions are compiled separately)
###############################################################################

# neural_input_rcpp and evaluate_regressor_convolution are implemented in C++
# (see the provided Rcpp source code).

###############################################################################
# Null Regressor
###############################################################################

#' Create a null regressor.
#'
#' @param hrf A hemodynamic response function.
#' @param span The HRF span (in seconds).
#' @noRd
null_regressor <- function(hrf = HRF_SPMG1, span = 24) {
  ret <- list(onsets = NA, hrf = hrf, eval = hrf, duration = 0,
              amplitude = 0, span = span)
  class(ret) <- c("null_regressor", "regressor", "list")
  ret
}

###############################################################################
# Single Trial Regressor
###############################################################################

#' Create a Single Trial Regressor.
#'
#' Creates a regressor object for modeling a single trial event in an fMRI experiment.
#' This is useful for trial‐wise analyses. The regressor represents the predicted BOLD
#' response for a single event using a specified hemodynamic response function (HRF).
#'
#' @param onsets A numeric vector of event onsets in seconds (must be of length 1).
#' @param hrf A hemodynamic response function (e.g. \code{HRF_SPMG1}).
#' @param duration The duration of the event (default is 0).
#' @param amplitude The amplitude scaling (default is 1).
#' @param span The temporal window of the HRF impulse response (default is 24 seconds).
#'
#' @return A \code{single_trial_regressor} object (inherits from \code{regressor}).
#'
#' @seealso \code{\link{regressor}}, \code{\link{evaluate.single_trial_regressor}}, \code{\link{HRF_SPMG1}}
#'
#' @examples
#' # Create a single trial regressor with default HRF
#' reg <- single_trial_regressor(10, HRF_SPMG1)
#'
#' # Create a single trial regressor with a custom HRF and amplitude
#' reg <- single_trial_regressor(10, HRF_GAMMA, amplitude = 2)
#'
#' # Evaluate the regressor over a time grid
#' grid <- seq(0, 30, by = 0.5)
#' values <- evaluate(reg, grid)
#'
#' @export
single_trial_regressor <- function(onsets, hrf = HRF_SPMG1, duration = 0, amplitude = 1, span = 24) {
  assert_that(length(onsets) == 1, msg = "Length of 'onsets' must be 1 for single trial regressor")
  assert_that(length(duration) == 1, msg = "Length of 'duration' must be 1 for single trial regressor")
  assert_that(length(amplitude) == 1, msg = "Length of 'amplitude' must be 1 for single trial regressor")
  assert_that(is.function(hrf))
  
  ret <- list(onsets = onsets,
              hrf = hrf,
              eval = hrf,
              duration = duration,
              amplitude = amplitude,
              span = span)
  class(ret) <- c("single_trial_regressor", "regressor", "list")
  ret
}

###############################################################################
# Multiple Trial Regressor
###############################################################################

#' Construct a Regressor Object.
#'
#' Constructs a \code{regressor} object from a set of event onsets and a hemodynamic
#' response function. The resulting regressor can be evaluated over time to produce a
#' predicted BOLD response.
#'
#' @param onsets A numeric vector of event onsets (in seconds).
#' @param hrf A hemodynamic response function (e.g., \code{HRF_SPMG1}).
#' @param duration The duration of each event (default is 0). If a single value is provided,
#'   it is repeated to match the length of \code{onsets}.
#' @param amplitude The amplitude scaling vector (default is 1). If a single value is provided,
#'   it is repeated to match the length of \code{onsets}.
#' @param span The temporal window for the HRF impulse response (default is 40 seconds).
#' @param summate Logical; if \code{TRUE}, the HRF amplitudes for events are summed when the
#'   event duration is greater than 0.
#'
#' @return A \code{regressor} object.
#'
#' @export
#'
#' @examples
#' reg <- regressor(c(10, 12, 14, 16, 18, 40), HRF_SPMG1, duration = 3)
#' pred <- evaluate(reg, seq(0, 100, by = 2))
#' nbasis(reg) == 1
#'
#' reg2 <- regressor(c(10, 12, 14, 16, 18, 40), HRF_SPMG1, duration = 3, summate = FALSE)
#' pred2 <- evaluate(reg2, seq(0, 100, by = 2))
regressor <- function(onsets, hrf = HRF_SPMG1, duration = 0, amplitude = 1, span = 40, summate = TRUE) {
  # Ensure 'duration' and 'amplitude' are vectors of the correct length.
  duration <- conform_len(duration, length(onsets))
  amplitude <- conform_len(amplitude, length(onsets))
  
  # If any event has a duration longer than half the current span, update span.
  if (any(duration > span / 2)) {
    span <- max(duration) * 2
  }
  
  assert_that(is.function(hrf))
  assert_that(length(amplitude) == length(onsets))
  assert_that(length(duration) == length(onsets))
  
  # Keep only events with nonzero, non-NA amplitude.
  keep <- which(amplitude != 0 & !is.na(amplitude))
  if (length(keep) == 0) {
    warning("regressor: onsets vector has no non-NA, nonzero elements")
    ret <- list(onsets = NA, hrf = hrf, eval = hrf,
                duration = 0, amplitude = 0, span = span, summate = summate)
  } else {
    ret <- list(onsets = onsets[keep], hrf = hrf, eval = hrf,
                duration = duration[keep], amplitude = amplitude[keep],
                span = span, summate = summate)
  }
  
  class(ret) <- c("regressor", "list")
  ret
}

###############################################################################
# Evaluation Methods
###############################################################################

#' Evaluate a Single Trial Regressor.
#'
#' Evaluates a \code{single_trial_regressor} over a given time grid by convolving the
#' event (with its onset, duration, and amplitude) with its HRF.
#'
#' @param x A \code{single_trial_regressor} object.
#' @param grid Numeric vector specifying the time points (in seconds) for evaluation.
#' @param precision Numeric value specifying the sampling precision for HRF evaluation (default: 0.3).
#' @param ... Additional arguments passed to the HRF evaluation function.
#'
#' @return A numeric vector (or matrix if multiple basis functions) representing the predicted BOLD response.
#'
#' @export
#' @method evaluate single_trial_regressor
#'
#' @examples
#' reg <- single_trial_regressor(10, HRF_SPMG1)
#' grid <- seq(0, 30, by = 0.5)
#' response <- evaluate(reg, grid)
evaluate.single_trial_regressor <- function(x, grid, precision = 0.3, ...) {
  nb <- nbasis(x)
  dspan <- x$span / median(diff(grid))
  
  delta <- grid - x$onsets
  grid.idx <- which(delta >= 0 & delta <= x$span)
  relOns <- grid[grid.idx] - x$onsets
  
  # IMPORTANT FIX: Pass x$duration (not x$amplitude) as the duration.
  resp <- evaluate(x$hrf, relOns, amplitude = x$amplitude,
                   duration = x$duration, precision = precision, ...)
  
  outmat <- matrix(0, length(grid), nb)
  outmat[grid.idx, 1:nb] <- resp
  
  if (nb == 1) {
    outmat[, 1]
  } else {
    outmat
  }
}

#' Evaluate a Null Regressor.
#'
#' Returns an evaluation of a null regressor (i.e. zero signal) over the grid.
#'
#' @param x A \code{null_regressor} object.
#' @param grid Numeric vector specifying the time points.
#' @param precision Numeric precision (default: 0.3).
#' @param ... Additional arguments.
#'
#' @export
evaluate.null_regressor <- function(x, grid, precision = 0.3, ...) {
  nb <- nbasis(x)
  outmat <- matrix(0, length(grid), nb)
  if (nb == 1) {
    outmat[, 1]
  } else {
    outmat
  }
}

###############################################################################
# Fast Convolution Evaluations
###############################################################################

#' Fast evaluation of regressor using convolution (Rcpp fast path).
#'
#' This function uses pre–evaluation of the HRF on a fine grid, constructs a neural input
#' vector using the Rcpp implementation, and then convolves the neural input with the HRF
#' using the Rcpp routine \code{evaluate_regressor_convolution}.
#'
#' @param x A \code{regressor} object.
#' @param start The starting time of the grid.
#' @param end The ending time of the grid.
#' @param TR The repetition time (grid spacing).
#' @param precision Sampling precision (default: 0.3).
#'
#' @keywords internal
fastevalreg <- function(x, start, end, TR, precision = 0.3) {
  nb <- nbasis(x)
  grid <- seq(start, end, by = TR)
  finegrid <- seq(start, end, by = precision)
  
  valid <- x$onsets >= (grid[1] - 16) & x$onsets <= grid[length(grid)]
  valid.ons <- x$onsets[valid]
  valid.amp <- x$amplitude[valid]
  valid.durs <- x$duration[valid]
  
  time <- seq(0, attr(x$hrf, "span"), by = precision)
  samhrf <- evaluate(x$hrf, time)
  
  # Create the 'delta' signal on the fine grid.
  bins <- as.integer((valid.ons - finegrid[1]) / precision + 1)
  delta <- numeric(length(finegrid))
  delta[bins] <- valid.amp
  
  if (nb > 1) {
    lowres <- matrix(0, length(grid), nb)
    for (i in 1:nb) {
      highres <- pracma::conv(samhrf[, i], delta)
      lowres[, i] <- approx(finegrid, highres[1:length(finegrid)], xout = grid)$y
    }
    lowres
  } else {
    highres <- pracma::conv(samhrf, delta)
    approx(finegrid, highres[1:length(finegrid)], xout = grid, rule = 2)$y
  }
}

#' Evaluate a Regressor over a Sampling Grid.
#'
#' Evaluates a \code{regressor} object over a given time grid. When the fast convolution flag
#' is enabled and the grid is uniform with equal durations, the fast Rcpp convolution approach
#' is used.
#'
#' @param x A \code{regressor} object.
#' @param grid A numeric vector representing the sampling grid (in seconds).
#' @param precision Sampling precision for HRF evaluation (default: 0.33).
#' @param use_conv Logical; if \code{TRUE}, use the fast convolution approach (requires uniform grid and equal durations).
#' @param ... Additional arguments.
#'
#' @return A numeric vector (or matrix if multiple basis functions) representing the regressor signal.
#'
#' @import RANN
#' @export
#'
#' @examples
#' frame <- sampling_frame(blocklens = 100, TR = 2)
#' reg <- regressor(onsets = c(10, 12, 14, 16, 18, 20, 35, 47, 52, 68, 79, 86),
#'                  hrf = HRF_SPMG1, duration = runif(12) * 3, amplitude = runif(12))
#' e1 <- evaluate(reg, samples(frame))
#' e2 <- evaluate(reg, samples(frame), use_conv = TRUE)
evaluate.regressor <- function(x, grid, precision = 0.33, use_conv = FALSE, ...) {
  nb <- nbasis(x)
  
  # Determine valid events intersecting with the grid.
  valid <- x$onsets >= (grid[1] - 16) & x$onsets <= grid[length(grid)]
  valid.ons <- x$onsets[valid]
  valid.durs <- x$duration[valid]
  valid.amp <- x$amplitude[valid]
  
  # Precompute HRF fine grid evaluation.
  time <- seq(0, attr(x$hrf, "span"), by = precision)
  
  # Use fast convolution path if possible.
  if (use_conv &&
      length(grid) > 1 &&
      length(unique(valid.durs)) == 1 &&
      all(diff(grid) == grid[2] - grid[1])) {
    start <- grid[1]
    end <- grid[length(grid)]
    TR <- grid[2] - grid[1]
    return(fastevalreg(x, start, end, TR, precision))
  }
  
  # Otherwise, use a nearest–neighbor approach to place event responses.
  nidx <- if (length(grid) > 1) {
    apply(RANN::nn2(matrix(grid), matrix(x$onsets), k = 2)$nn.idx, 1, min)
  } else {
    apply(RANN::nn2(matrix(grid), matrix(x$onsets), k = 1)$nn.idx, 1, min)
  }
  
  # Only use indices corresponding to valid events.
  nidx <- nidx[valid]
  outmat <- matrix(0, length(grid), length(valid.ons) * nb)
  
  if (length(valid.ons) == 0) {
    warning("none of the regressor onsets intersect with sampling 'grid'; evaluating to zero at all times.")
    return(matrix(0, length(grid), nb))
  }
  
  for (i in seq_along(valid.ons)) {
    grid.idx <- seq(nidx[i], min(nidx[i] + (x$span / median(diff(grid))), length(grid)))
    relOns <- grid[grid.idx] - valid.ons[i]
    resp <- evaluate(x$hrf, relOns, amplitude = valid.amp[i],
                     duration = valid.durs[i], precision = precision,
                     summate = x$summate, ...)
    
    if (nb > 1) {
      start_col <- (i - 1) * nb + 1
      end_col <- i * nb
      outmat[grid.idx, start_col:end_col] <- resp
    } else {
      outmat[grid.idx, i] <- resp
    }
  }
  
  if (length(valid.ons) > 1) {
    if (nb == 1) {
      rowSums(outmat)
    } else {
      do.call(cbind, lapply(1:nb, function(i) {
        rowSums(outmat[, seq(i, by = nb, length.out = length(valid.ons))])
      }))
    }
  } else {
    if (nb == 1) {
      outmat[, 1]
    } else {
      outmat
    }
  }
}

#' Alternative evaluation using Rcpp convolution.
#'
#' This function uses the compiled C++ routine \code{evaluate_regressor_convolution}
#' to evaluate the regressor over the grid.
#'
#' @param x A \code{regressor} object.
#' @param grid A numeric vector for sampling.
#' @param precision Sampling precision (default: 0.3).
#'
#' @keywords internal
evaluate_regressor <- function(x, grid, precision = 0.3, ...) {
  hrf_span <- attr(x$hrf, "span")
  hrf_times <- seq(0, hrf_span, by = precision)
  samhrf <- evaluate(x$hrf, hrf_times)
  if (is.vector(samhrf)) {
    samhrf <- matrix(samhrf, ncol = 1)
  }
  
  start <- min(grid[1], x$onsets[1]) - hrf_span
  end <- max(grid[length(grid)], x$onsets[length(x$onsets)]) + hrf_span
  
  result <- evaluate_regressor_convolution(
    grid = grid,
    onsets = x$onsets,
    durations = x$duration,
    amplitudes = x$amplitude,
    hrf_values = samhrf,
    hrf_span = hrf_span,
    start = start,
    end = end,
    precision = precision
  )
  
  nb <- nbasis(x)
  if (nb == 1) as.vector(result) else result
}

###############################################################################
# S3 Methods for regressor objects
###############################################################################

#' @export
nbasis.regressor <- function(x) nbasis(x$hrf)

#' @export
nbasis.HRF <- function(x) attr(x, "nbasis")

#' @export
nbasis.AFNI_HRF <- function(x) attr(x, "nbasis")

#' @export
#' @rdname nbasis
nbasis.hrfspec <- function(x) nbasis(x$hrf)

#' @export
#' @rdname onsets
onsets.regressor <- function(x) x$onsets

#' @export
#' @rdname durations
durations.regressor <- function(x) x$duration

#' @export
#' @rdname amplitudes
amplitudes.regressor <- function(x) x$amplitude

#' Neural input for a regressor using the Rcpp routine.
#'
#' @param x A regressor object.
#' @param start The start time.
#' @param end The end time.
#' @param resolution The sampling resolution (default: 0.33).
#'
#' @export
neural_input.regressor <- function(x, start, end, resolution = 0.33, ...) {
  neural_input_rcpp(x, start, end, resolution)
}

#' (Legacy) Neural input function implemented in R.
#'
#' @param x A regressor object.
#' @param from The start time.
#' @param to The end time.
#' @param resolution The resolution.
#'
#' @keywords internal
neural_inputR <- function(x, from, to, resolution) {
  time <- seq(from + (resolution / 2), to - (resolution / 2), by = resolution)
  out <- numeric((to - from) / resolution)
  ons <- x$onsets
  dur <- x$duration
  amp <- x$amplitude
  
  for (i in seq_along(ons)) {
    on <- ons[i]
    d <- dur[i]
    startbin <- as.integer((on - from) / resolution) + 1
    endbin <- if (d > 0) {
      as.integer((on - from) / resolution + d / resolution) + 1
    } else {
      startbin
    }
    for (j in startbin:endbin) {
      out[j] <- out[j] + amp[i]
    }
  }
  
  ts(out, start = time[1], frequency = 1 / resolution)
}

###############################################################################
# Plot and Print Methods
###############################################################################

#' Plot a Regressor Object.
#'
#' Plots the regressor signal over time. Vertical lines mark event onsets.
#'
#' @param x A regressor object.
#' @param samples Numeric vector of times (in seconds) for plotting. If missing, an automatic range is computed.
#' @param add Logical; if \code{TRUE}, add to the current plot.
#' @param ... Additional graphical parameters passed to \code{plot}.
#'
#' @export
plot.regressor <- function(x, samples, add = FALSE, ...) {
  if (missing(samples)) {
    ons <- x$onsets
    samples <- seq(ons[1] - 5, ons[length(ons)] + 12, by = 1)
  }
  
  y <- evaluate(x, samples)
  if (add) {
    graphics::lines(samples, y)
  } else {
    plot(samples, y, type = "l", xlab = "Time", ylab = "Amplitude", ...)
  }
  srange <- range(samples)
  for (on in onsets(x)) {
    if (on >= srange[1] && on <= srange[2]) {
      graphics::abline(v = on, col = 2, lty = 2)
    }
  }
}

#' Print Method for Regressor Objects.
#'
#' Provides a summary of the regressor structure, temporal properties, and HRF details.
#'
#' @param x A regressor object.
#' @param ... Additional parameters (ignored).
#'
#' @export
print.regressor <- function(x, ...) {
  cat("\n═══ fMRI Regressor ═══\n")
  cat("\n📊 Structure:\n")
  cat("  • Type:", class(x)[1], "\n")
  cat("  • Number of basis functions:", nbasis(x), "\n")
  cat("  • Number of onsets:", length(onsets(x)), "\n")
  
  cat("\n⏱️  Temporal Properties:\n")
  cat("  • Onsets:", paste(head(round(onsets(x), 2), 3),
                           if (length(onsets(x)) > 3) "..." else "", collapse = ", "), "\n")
  if (!is.null(x$duration) && any(x$duration != 0)) {
    cat("  • Durations:", paste(head(round(x$duration, 2), 3),
                                if (length(x$duration) > 3) "..." else "", collapse = ", "), "\n")
  }
  if (!is.null(x$amplitude) && !all(x$amplitude == 1)) {
    cat("  • Amplitudes:", paste(head(round(x$amplitude, 2), 3),
                                 if (length(x$amplitude) > 3) "..." else "", collapse = ", "), "\n")
  }
  
  cat("\n🧠 HRF Properties:\n")
  cat("  • Type:", class(x$hrf)[1], "\n")
  cat("  • Span:", attr(x$hrf, "span"), "seconds\n")
  cat("\n")
}

###############################################################################
# Shift Method for Regressor Objects
###############################################################################

#' Shift a Regressor.
#'
#' Shifts the event onsets of a regressor by a specified amount.
#'
#' @param x A regressor object.
#' @param shift_amount Numeric value indicating the time (in seconds) to shift. Positive shifts right.
#' @param ... Additional arguments (ignored).
#'
#' @return A new regressor object with shifted onsets.
#'
#' @export
#' @family shift
shift.regressor <- function(x, shift_amount, ...) {
  if (!inherits(x, "regressor")) {
    stop("x must be an object of class 'regressor'")
  }
  if (!is.numeric(shift_amount)) {
    stop("shift_amount must be a numeric value")
  }
  
  shifted_onsets <- x$onsets + shift_amount
  shifted_regressor <- regressor(onsets = shifted_onsets,
                                 hrf = x$hrf,
                                 duration = x$duration,
                                 amplitude = x$amplitude,
                                 span = x$span,
                                 summate = x$summate)
  shifted_regressor
}

  
                