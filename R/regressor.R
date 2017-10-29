#' @importFrom RANN nn2
NULL


#' regressor 
#' 
#' construct a regressor function that can be used to generate a regression fucntion 
#' from a set of onset times and a hemodynamic response function
#' 
#' @param onset the event onsets in seconds
#' @param hrf a hemodynamic response function
#' @param duration duration of events (default is 0)
#' @param amplitude scaling vector (default is 1)
#' @param span the temporal window of the impulse response function (default is 24)
#' @return an S3 list of type \code{regressor}
#' @export
#' @examples 
#' 
#' reg <- regressor(c(10,12,14,16,18, 40), HRF_SPMG1)
#' pred <- evaluate(reg, seq(0,100,by=2))
#' nbasis(reg) == 1
regressor <- function(onsets, hrf, duration=0, amplitude=1, span=24) {
  if (length(duration) == 1) {
    duration = rep(duration, length(onsets))
  }
  
  if (length(amplitude) == 1) {
    amplitude = rep(as.vector(amplitude), length(onsets))
  }
  
  
  keep <- which(amplitude != 0)
  empty <- length(keep) == 0
  ret <- if (!empty) {
    list(onsets=onsets[keep],hrf=hrf, eval=hrf, duration=duration[keep],amplitude=amplitude[keep],span=span)  
  } else {
    list(onsets=NA,hrf=hrf, eval=hrf, duration=0,amplitude=0,span=span)  
  }
  class(ret) <- c("regressor", "list")
  ret
}


dots <- function(...) {
  eval(substitute(alist(...)))
}


#' evaluate
#' @rdname evaluate
#' @param grid the sampling grid. A vector of real values in seconds.
#' @param precision the sampling precision for the hrf. This parameter is passed to \code{evaluate.HRF}
#' @export
evaluate.regressor <- function(x, grid, precision=.1) {
  nb <- nbasis(x)
  dspan <- x$span/median(diff(grid)) 
  
  if (is.na(x$onsets) || length(x$onsets) == 0) {
    outmat <- matrix(0, length(grid), nb)
    return(outmat)
  }
    
  outmat <- matrix(0, length(grid), length(x$onsets) * nb)
  
  
  nidx <- if (length(grid) > 1) {
    apply(RANN::nn2(matrix(grid), matrix(x$onsets), k=2)$nn.idx, 1, min)
  } else {
    apply(RANN::nn2(matrix(grid), matrix(x$onsets), k=1)$nn.idx, 1, min)
  }
  
  valid <- x$onsets >= grid[1] & x$onsets < grid[length(grid)]
  
  if (all(!valid)) {
    warning("none of the regressor onsets intersect with sampling 'grid', evalauting to zero at all times.")
    return(outmat)
  }
  
  
  valid.ons <- x$onsets[valid]
  valid.durs <- x$duration[valid]
  valid.amp <- x$amplitude[valid]
  
  nidx <- nidx[valid]

  for (i in seq_along(valid.ons)) { 
    grid.idx <- seq(nidx[i], min(nidx[i] + dspan, length(grid)))             
    relOns <- grid[grid.idx] - valid.ons[i]    
    resp <- evaluate(x$hrf, relOns, amplitude=valid.amp[i], duration=valid.durs[i], precision=precision)   
  
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
    outmat
  }
}

#' @export
nbasis.regressor <- function(x) nbasis(x$hrf)

#' @export
nbasis.HRF <- function(x) attr(x, "nbasis")

#' @export
nbasis.hrfspec <- function(x) nbasis(x$hrf)

#' @export
onsets.regressor <- function(x) x$onsets

#' @export
durations.regressor <- function(x) x$duration

#' @export
amplitudes.regressor <- function(x) x$amplitude


#' @export
plot.regressor <- function(object, samples, add=FALSE, ...) {
  if (missing(samples)) {
    ons <- object$onsets
    samples <- seq(ons[1], ons[length(ons)], by=1)
  }
  
  y <- evaluate(object, samples)
  if (add){
    lines(samples, y)
  } else {
    plot(samples, y, type='l', xlab="Time", ylab="Amplitude", ...)
  }
  srange <- range(samples)
  
  for (on in onsets(object)) {
    if (on >= srange[1] && on <= srange[2]) {
      abline(v=on, col=2, lty=2)
    }
  }
  
}

#' @export
print.regressor <- function(object) {
  N <- min(c(6, length(onsets(object))))
  cat(paste("hemodynamic response function:", attr(object$hrf, "name")))
  cat("\n")
  cat(paste("onsets: ", paste(onsets(object)[1:N], collapse=" "), "..."))
  cat("\n")
  
  durs <- durations(object)
  amps <- amplitudes(object)
  
  if (all(durs == durs[1])) {
    cat(paste("durations: ", durs[1], "for all events"))
  } else {
    cat(paste("durations: ", paste(durs[1:N], collapse=" "), "..."))
  }
  
  cat("\n")
  
  if (!is.na(amps[1]) && all(amps == amps[1])) {
    cat(paste("amplitudes: ", amps[1], "for all events"))
  } else {     
    cat(paste("amplitudes: ", paste(amps[1:N], collapse=" "), "..."))
    
  }
  
  cat("\n")
}

#' check vector is of length 1 or repeated for supplied length
conformLen <- function(val, len) {
  name <- deparse(substitute(val))
  if (length(val) == 1) {
    rep(val, len)
  } else if (length(val) == len) {
    val
  } else {
    stop(paste(name, "must be of length 1 or same length as onsets"))
  }
}


  
                