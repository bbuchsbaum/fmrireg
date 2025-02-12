#' null_regressor 
#' @param hrf an hrf function
#' @param span the hrf span
#' @noRd
null_regressor <- function(hrf=HRF_SPMG1, span=24) {
  ret <- list(onsets=NA,hrf=hrf, eval=hrf, duration=0,amplitude=0,span=span)
  class(ret) <- c("null_regressor", "regressor", "list")
  ret
}


#' Create a single trial regressor
#'
#' Creates a regressor object for modeling a single trial event in an fMRI experiment.
#' This is particularly useful for trial-wise analyses where each trial needs to be
#' modeled separately. The regressor represents the predicted BOLD response for a single
#' event using a specified hemodynamic response function (HRF).
#'
#' @param onsets the event onset in seconds, must be of length 1.
#' @param hrf a hemodynamic response function, e.g. \code{HRF_SPMG1}
#' @param duration duration of the event (default is 0)
#' @param amplitude scaling vector (default is 1)
#' @param span the temporal window of the impulse response function (default is 24)
#'
#' @return A single_trial_regressor object which inherits from regressor with the following components:
#'   \itemize{
#'     \item onsets: The onset time of the trial
#'     \item duration: The duration of the trial
#'     \item amplitude: The amplitude of the trial response
#'     \item hrf: The HRF object used for convolution
#'   }
#' @seealso 
#'   \code{\link{regressor}} for creating multiple trial regressors,
#'   \code{\link{evaluate.single_trial_regressor}} for evaluating the regressor,
#'   \code{\link{HRF_SPMG1}} for the default HRF
#'
#' @examples
#' # Create a single trial regressor with default HRF
#' reg <- single_trial_regressor(c(10), HRF_SPMG1)
#'
#' # Create a single trial regressor with custom HRF and amplitude
#' reg <- single_trial_regressor(c(10), HRF_GAMMA, amplitude = 2)
#'
#' # Evaluate the regressor over a time grid
#' grid <- seq(0, 30, by = 0.5)
#' values <- evaluate(reg, grid)
#'
#' @export
single_trial_regressor <- function(onsets, hrf=HRF_SPMG1, duration=0, amplitude=1, span=24) {
  assert_that(length(onsets) ==1, msg="length of 'onsets' must be 1 for single trial regressor")
  assert_that(length(duration) ==1, msg="length of 'duration' must be 1 for single trial regressor")
  assert_that(length(amplitude) ==1, msg="length of 'amplitude' must be 1 for single trial regressor")
  
  assertthat::assert_that(is.function(hrf))
  
  ret <- list(onsets=onsets,hrf=hrf, eval=hrf, duration=duration,amplitude=amplitude,span=span)  
  class(ret) <- c("single_trial_regressor", "regressor", "list")
  ret
}


#' construct a regressor object
#' 
#' construct a \code{regressor} object that can be used to generate regression variables
#' from a set of onset times and a hemodynamic response function. A \code{regressor} can be
#' evaluated at a set of times to generate a time-course appropriate for modeling an fMRI response.
#' 
#' @param onsets the event onsets in seconds
#' @param hrf a hemodynamic response function, e.g. \code{HRF_SPMG1} or costum \code{HRF}
#' @param duration duration of events (default is 0)
#' @param amplitude scaling vector (default is 1)
#' @param span the temporal window of the impulse response function (default is 24)
#' @param summate whether to summate hrf amplitude as a function of the duration of an event.
#' @return an S3 list of type \code{regressor}
#' @export
#' @examples 
#' 
#' reg <- regressor(c(10,12,14,16,18, 40), HRF_SPMG1, duration=3)
#' pred <- evaluate(reg, seq(0,100,by=2))
#' nbasis(reg) == 1
#' 
#' reg2 <- regressor(c(10,12,14,16,18, 40), HRF_SPMG1, duration=3, summate=FALSE)
#' pred2 <- evaluate(reg2, seq(0,100,by=2))
#' stopifnot(max(pred) > max(pred2))
regressor <- function(onsets, hrf=HRF_SPMG1, duration=0, amplitude=1, span=40, summate=TRUE) {
  if (length(duration) == 1) {
    duration = rep(duration, length(onsets))
  }
  
  if (length(amplitude) == 1) {
    amplitude = rep(as.vector(amplitude), length(onsets))
  }
  
  if (any(duration > span/2)) {
    span <- max(duration) * 2
  }
  
  assertthat::assert_that(is.function(hrf))
  assertthat::assert_that(length(amplitude) == length(onsets))
  assertthat::assert_that(length(duration) == length(onsets))
 
  
  keep <- which(amplitude != 0 & !is.na(amplitude))
  empty <- length(keep) == 0
  
  ret <- if (!empty) {
    list(onsets=onsets[keep],hrf=hrf, eval=hrf, duration=duration[keep],
         amplitude=amplitude[keep],span=span,summate=summate)  
  } else {
    warning("regressor: onsets vector has no non-NA elements")
    list(onsets=NA,hrf=hrf, eval=hrf, duration=0,amplitude=0,span=span,summate=summate)  
  }
  
  class(ret) <- c("regressor", "list")
  ret
}

#' @keywords internal
#' @noRd
dots <- function(...) {
  eval(substitute(alist(...)))
}



#' Evaluate a single trial regressor
#'
#' Method to evaluate a single trial regressor over a specified time grid.
#' This computes the predicted BOLD response for the trial at each time point
#' by convolving the trial event with its hemodynamic response function.
#'
#' @param x A single_trial_regressor object
#' @param grid Numeric vector specifying the time points at which to evaluate the regressor
#' @param precision Optional numeric value specifying the precision of the evaluation (default: 0.3)
#' @param ... Additional arguments passed to the HRF evaluation function
#'
#' @return A numeric vector containing the predicted BOLD response at each time point in the grid
#' @seealso 
#'   \code{\link{single_trial_regressor}} for creating single trial regressors,
#'   \code{\link{evaluate.HRF}} for HRF evaluation details
#'
#' @examples
#' # Create and evaluate a single trial regressor
#' reg <- single_trial_regressor(c(10), HRF_SPMG1)
#' grid <- seq(0, 30, by = 0.5)
#' response <- evaluate(reg, grid)
#'
#' # Plot the response
#' plot(grid, response, type = "l", 
#'      xlab = "Time (s)", ylab = "BOLD Response")
#'
#' @method evaluate single_trial_regressor
#' @export
evaluate.single_trial_regressor <- function(x, grid, precision=.3, ...) {
  nb <- nbasis(x)
  dspan <- x$span/median(diff(grid)) 
  
  delta <- grid - x$onsets 
  grid.idx <- which(delta >= 0 & delta <= x$span)
  relons <- grid[grid.idx] - x$onsets    
  resp <- evaluate(x$hrf, relons, amplitude=x$amplitude, duration=x$duration, precision=precision)   

  outmat <- matrix(0, length(grid), nb)
  outmat[grid.idx,1:nb] <- resp
  
  if (nb == 1) {
    outmat[,1]
  } else {
    outmat
  }
}


#' @export
evaluate.null_regressor <- function(x, grid, precision=.3, ...) {
  nb <- nbasis(x)
  dspan <- x$span/median(diff(grid)) 
  
  
  outmat <- matrix(0, length(grid), nb)
  
  if (nb == 1) {
    outmat[,1]
  } else {
    outmat
  }
}

#' @keywords internal 
#' @noRd
fastevalreg2 <- function(x, start, end, TR, precision=.3) {
  nb <- nbasis(x)
  #grid <- seq(start, end, by=TR)
  #finegrid <- seq(start, end, by=precision)
  #valid <- x$onsets >= (grid[1]-16) & x$onsets <= grid[length(grid)]
  
  #valid.ons <- x$onsets[valid]
  #valid.amp <- x$amplitude[valid]
  #valid.durs <- x$duration[valid]
  
  time <- seq(0,attr(x$hrf, "span"), by=precision)
  samhrf <- evaluate(x$hrf, time)
  
  ninput <- neural_input(x, start, end, resolution=precision)
  
  if (nb > 1) {
    lowres <- matrix(0, length(grid), nb)
    for (i in 1:nb) {
      highres <- pracma::conv(samhrf[,i],ninput$neural_input)
      lowres[,i] <- approx(ninput$time, highres[1:length(ninput$time)], xout=grid, rule=2)$y
    }
    lowres
  } else {
    highres <- pracma::conv(samhrf,ninput$neural_input)
    lowres <- approx(ninput$time, highres[1:length(ninput$time)], xout=grid, rule=2)$y
    lowres
  }
  
}

#' @importFrom pracma conv
#' @importFrom memoise memoise
#' @noRd
fastevalreg <- function(x, start, end, TR, precision=.3) {
  nb <- nbasis(x)
  grid <- seq(start, end, by=TR)
  finegrid <- seq(start, end, by=precision)
  valid <- x$onsets >= (grid[1]-16) & x$onsets <= grid[length(grid)]
  
  valid.ons <- x$onsets[valid]
  valid.amp <- x$amplitude[valid]
  valid.durs <- x$duration[valid]
  
  time <- seq(0,attr(x$hrf, "span"), by=precision)
  samhrf <- evaluate(x$hrf, time)
  
  bins <- as.integer((valid.ons - finegrid[1])/precision + 1)
  delta <- numeric(length(finegrid))
  delta[bins] <- x$amplitude
  
  
  #highres <- stats::convolve(samhrf, delta, type="open")
  
  if (nb > 1) {
    lowres <- matrix(0, length(grid), nb)
    for (i in 1:nb) {
      highres <- pracma::conv(samhrf[,i],delta)
      lowres[,i] <- approx(finegrid, highres[1:length(finegrid)], xout=grid)$y
    }
    lowres
  } else {

    highres <- pracma::conv(samhrf,delta)
    lowres <- approx(finegrid, highres[1:length(finegrid)], xout=grid, rule=2)$y
    lowres
  }
}


#' evaluate
#' 
#' evalute a regressor function over a sampling grid
#' 
#' @rdname evaluate
#' @param grid the sampling grid. A vector of real values in seconds.
#' @param precision the sampling precision for the hrf. This parameter is passed to \code{evaluate.HRF}
#' @param use_conv use fast convolution approach
#' @examples 
#' frame <- sampling_frame(blocklens=100, TR=2)
#' reg <- regressor(onsets=c(10,12, 14, 16, 18,20,35, 47,52, 68, 79,86), amp=runif(12), 
#' duration=runif(12)*3, hrf=HRF_SPMG1)
#' e1 = evaluate(reg, samples(frame))
#' e2 = evaluate(reg, samples(frame), use_conv=TRUE)
#' @import RANN
#' @export
evaluate.regressor <- function(x, grid, precision=.33, use_conv=FALSE, ...) {
 
  nb <- nbasis(x)
  dspan <- x$span/median(diff(grid)) 
  

  valid <- x$onsets >= (grid[1]-16) & x$onsets <= grid[length(grid)]
  valid.ons <- x$onsets[valid]
  valid.durs <- x$duration[valid]
  valid.amp <- x$amplitude[valid]
  
  time <- seq(0,attr(x$hrf, "span"), by=precision)
  #maxamp <- max(valid.amp)
  #samhrf <- evaluate(x$hrf, time, amplitude=maxamp)
  #maxheight <- max(samhrf)
  
  if (use_conv && all(valid.durs[1] == valid.durs) && all(diff(grid) == grid[2] - grid[1]) && length(grid) > 1) {
    start <- grid[1]
    end <- grid[length(grid)]
    TR <- grid[2] - grid[1]
    return(fastevalreg(x, start, end, TR, precision))
  } 

  nidx <- if (length(grid) > 1) {
    #apply(rflann::Neighbour(matrix(x$onsets), matrix(grid), k=2,build = "kdtree",cores=0, checks=1)$indices, 1, min)
    apply(RANN::nn2(matrix(grid), matrix(x$onsets), k=2)$nn.idx, 1, min)
  } else {
    #apply(rflann::Neighbour(matrix(x$onsets), matrix(grid), k=1,build = "kdtree",cores=0, checks=1)$indices, 1, min)
    apply(RANN::nn2(matrix(grid), matrix(x$onsets), k=1)$nn.idx, 1, min)
  }
  

  outmat <- matrix(0, length(grid), length(valid.ons) * nb)
  
  if (length(valid) == 0 || all(is.na(valid)) || all(!valid)) {
    warning("none of the regressor onsets intersect with sampling 'grid', evaluating to zero at all times.")
    return(matrix(0, length(grid), nb))
  }

  nidx <- nidx[valid]
  
  #browser()

  for (i in seq_along(valid.ons)) { 
    grid.idx <- seq(nidx[i], min(nidx[i] + dspan, length(grid)))             
    relOns <- grid[grid.idx] - valid.ons[i]    
    resp <- evaluate(x$hrf, relOns, amplitude=valid.amp[i], 
                     duration=valid.durs[i], 
                     precision=precision,
                     summate=x$summate)   
  
    if (nb > 1) {
      start <- (i-1) * nb + 1
      end <- i*nb 
      outmat[grid.idx,start:end] <- resp
    } else {
     
      outmat[grid.idx, i] <- resp
    }
  }
  
  
  if (length(valid.ons) > 1) {
    if (nb == 1) {
      rowSums(outmat)
    } else {
      do.call(cbind, lapply(1:nb, function(i) {
        rowSums(outmat[,seq(i, by=nb, length.out=length(valid.ons))])
      }))
    }
  } else {
    if (nb == 1) {
      outmat[,1]
    } else {
      outmat
    }
  }
}

#' @keywords internal
#' @noRd
evaluate_regressor <- function(x, grid, precision = 0.3, ...) {
  # Extract HRF span
  hrf_span <- attr(x$hrf, "span")
  
  # Generate HRF values over fine grid
  hrf_times <- seq(0, hrf_span, by = precision)
  samhrf <- evaluate(x$hrf, hrf_times)
  
  # Ensure samhrf is a matrix
  if (is.vector(samhrf)) {
    samhrf <- matrix(samhrf, ncol = 1)
  }
  
  # Determine start and end times
  start <- min(grid[1], x$onsets[1]) - hrf_span
  end <- max(grid[length(grid)], x$onsets[length(x$onsets)]) + hrf_span
  
  # Call the Rcpp function
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
  
  # Return result in the expected format
  nb <- nbasis(x)
  if (nb == 1) {
    return(as.vector(result))
  } else {
    return(result)
  }
}



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


#' @export
neural_input.regressor <- function(x, start, end, resolution=.33, ...) {
  neural_input_rcpp(x, start, end, resolution)
}

#' @keywords internal
#' @noRd
neural_inputR <- function(x, from, to, resolution) {
  time <- seq(from + (resolution/2), to - (resolution/2), by=resolution)
  out <- numeric( (to - from)/resolution)
  ons <- x$onsets
  dur <- x$duration
  amp <- x$amplitude
  
  for (i in seq_along(ons)) {
    on <- ons[i]
    d <- dur[i]
    startbin <- as.integer((on - from)/resolution) + 1

    if (d > 0) {
      endbin <- as.integer((on - from)/resolution + d/resolution) + 1
    } else {
      endbin <- startbin
    }
    
    for (j in startbin:endbin) {
      out[j] <- out[j] + amp[i]
    }
  }
  
  ts(out, start=time[1], frequency=1/resolution)
}

# #include <Rcpp.h>
# using namespace Rcpp;
# 
# // [[Rcpp::export]]
# List neural_input_rcpp(List x, double from, double to, double resolution) {
#   int n = (to - from) / resolution;
#   NumericVector time(n);
#   NumericVector out(n);
#   NumericVector ons = x["onsets"];
#   NumericVector dur = x["duration"];
#   NumericVector amp = x["amplitude"];
#   
#   for (int i = 0; i < ons.length(); i++) {
#     double on = ons[i];
#     double d = dur[i];
#     int startbin = (int) ((on - from) / resolution) + 1;
#     if (d > 0) {
#       int endbin = (int) ((on - from) / resolution + d / resolution) + 1;
#       for (int j = startbin; j <= endbin; j++) {
#         out[j-1] += amp[i];
#       }
#     } else {
#       out[startbin-1] += amp[i];
#     }
#   }
#   
#   for (int i = 0; i < n; i++) {
#     time[i] = from + (i + 0.5) * resolution;
#   }
#   
#   List result;
#   result["time"] = time;
#   result["neural_input"] = out;
#   return result;
# }


  


#' plot a regressor object
#' 
#' @export
#' @param x the object
#' @param samples the times in seconds along which to plot regressor function
#' @param add whether to add to existing plot
#' @param ... extra args to send to `plot` 
plot.regressor <- function(x, samples, add=FALSE, ...) {
  if (missing(samples)) {
    ons <- x$onsets
    samples <- seq(ons[1]-5, ons[length(ons)] + 12, by=1)
  }
  
  y <- evaluate(x, samples)
  if (add){
    graphics::lines(samples, y)
  } else {
    plot(samples, y, type='l', xlab="Time", ylab="Amplitude", ...)
  }
  srange <- range(samples)
  
  for (on in onsets(x)) {
    if (on >= srange[1] && on <= srange[2]) {
      graphics::abline(v=on, col=2, lty=2)
    }
  }
  
}

#' @export
print.regressor <- function(x, ...) {
  # Header with fancy box drawing
  cat("\n═══ fMRI Regressor ═══\n")
  
  # Basic information about the regressor
  cat("\n📊 Structure:\n")
  cat("  • Type:", class(x)[1], "\n")
  cat("  • Number of basis functions:", nbasis(x), "\n")
  cat("  • Number of onsets:", length(onsets(x)), "\n")
  
  # Temporal information
  cat("\n⏱️  Temporal Properties:\n")
  cat("  • Onsets:", paste(head(round(onsets(x), 2), 3), 
                          if(length(onsets(x)) > 3) "..." else "", collapse=", "), "\n")
  if (!is.null(x$duration) && any(x$duration != 0)) {
    cat("  • Durations:", paste(head(round(x$duration, 2), 3), 
                               if(length(x$duration) > 3) "..." else "", collapse=", "), "\n")
  }
  if (!is.null(x$amplitude) && !all(x$amplitude == 1)) {
    cat("  • Amplitudes:", paste(head(round(x$amplitude, 2), 3), 
                                if(length(x$amplitude) > 3) "..." else "", collapse=", "), "\n")
  }
  
  # HRF information
  cat("\n🧠 HRF Properties:\n")
  cat("  • Type:", class(x$hrf)[1], "\n")
  cat("  • Span:", attr(x$hrf, "span"), "seconds\n")
  
  # Add a nice footer
  cat("\n")
}

#' check vector is of length 1 or repeated for supplied length
#' @keywords internal
#' @noRd
conform_len <- function(val, len) {
  name <- deparse(substitute(val))
  if (length(val) == 1) {
    rep(val, len)
  } else if (length(val) == len) {
    val
  } else {
    stop(paste(name, "must be of length 1 or same length as onsets"))
  }
}


#' Shift method for the regressor class
#' 
#' @param x the regressor object
#' @param shift_amount A numeric value indicating the amount of time to shift the object by.
#'   Positive values will shift the object to the right, while negative values will shift it to the left.
#' @param ... extra args
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
  shifted_regressor <- regressor(onsets = shifted_onsets, hrf = x$hrf, duration = x$duration, amplitude = x$amplitude, span = x$span, summate = x$summate)
  
  return(shifted_regressor)
}


  
                