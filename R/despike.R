

#' Despike Time Series Data
#' 
#' @description
#' Remove spikes from time series data using running median and a threshold.
#'
#' @param x A numeric vector containing the time series data.
#' @param k The window size for the running median (default: 7).
#' @param thresh The threshold for identifying spikes (default: 6).
#' 
#' @return A numeric vector with the spikes removed.
#' @export
despike <- function(x, k=7, thresh=6) {
  y <- runmed(x, k=k)
  delta <- x - y
  mad <- median(abs(delta))
  z <- abs(delta/mad)
  idx <- which(z > thresh)
  x[idx] <- y[idx]
  
  attr(x,  "idx") <- idx
  x
}



#' Compute Instantaneous Correlation between Two Signals
#'
#' Computes the instantaneous correlation between two signals using an exponential decay weighting function, allowing for efficient recursive computation.
#'
#' @param x Numeric vector. First time-series signal.
#' @param y Numeric vector. Second time-series signal.
#' @param eta Positive numeric value. Decay parameter for the exponential weighting function.
#' @param offset Integer. Temporal offset between the two signals. Defaults to 0.
#' @return Numeric vector of instantaneous correlation coefficients.
#' @examples
#' # Example signals
#' x <- sin(seq(0, 2 * pi, length.out = 100))
#' y <- cos(seq(0, 2 * pi, length.out = 100))
#' # Compute instantaneous correlation
#' rho <- instantaneous_correlation(x, y, eta = 0.1)
#' # Plot the result
#' plot(rho, type = 'l', ylab = 'Instantaneous Correlation', xlab = 'Time')
#' @export
instantaneous_correlation <- function(x, y, eta, offset = 0) {
  # Ensure inputs are numeric vectors of the same length
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric vectors.")
  }
  N <- length(x)
  if (length(y) != N) {
    stop("x and y must be of the same length.")
  }
  
  # Decay factor
  a <- exp(-eta)
  
  # Initialize variables
  x_mean <- numeric(N)
  y_mean <- numeric(N)
  s_xy <- numeric(N)
  s_xx <- numeric(N)
  s_yy <- numeric(N)
  rho <- numeric(N)
  
  # Initial values
  x_mean[1] <- x[1]
  y_index <- 1 - offset
  y_value <- ifelse(y_index < 1 || y_index > N, 0, y[y_index])
  y_mean[1] <- y_value
  s_xy[1] <- x[1] * y_value
  s_xx[1] <- x[1]^2
  s_yy[1] <- y_value^2
  rho[1] <- 0  # Correlation is undefined at first point
  
  # Recursive computation
  for (k in 2:N) {
    # Update indices and values with offset
    xk <- x[k]
    y_index <- k - offset
    yk <- ifelse(y_index < 1 || y_index > N, 0, y[y_index])
    
    # Update means
    x_mean[k] <- a * x_mean[k - 1] + (1 - a) * xk
    y_mean[k] <- a * y_mean[k - 1] + (1 - a) * yk
    
    # Update covariances
    s_xy[k] <- a * s_xy[k - 1] + (1 - a) * xk * yk
    s_xx[k] <- a * s_xx[k - 1] + (1 - a) * xk^2
    s_yy[k] <- a * s_yy[k - 1] + (1 - a) * yk^2
    
    # Compute covariance and variances
    cov_xy <- s_xy[k] - x_mean[k] * y_mean[k]
    var_xx <- s_xx[k] - x_mean[k]^2
    var_yy <- s_yy[k] - y_mean[k]^2
    
    # Compute instantaneous correlation
    denom <- sqrt(var_xx * var_yy)
    rho[k] <- ifelse(denom > 0, cov_xy / denom, 0)
  }
  
  return(rho)
}

#' Compute Instantaneous Correlation between Two Signals (Rcpp Version)
#'
#' Computes the instantaneous correlation between two signals using an exponential decay weighting function. This version uses Rcpp for improved performance.
#'
#' @param x Numeric vector. First time-series signal.
#' @param y Numeric vector. Second time-series signal.
#' @param eta Positive numeric value. Decay parameter for the exponential weighting function.
#' @param offset Integer. Temporal offset between the two signals. Defaults to 0.
#' @return Numeric vector of instantaneous correlation coefficients.
#' @examples
#' # Example signals
#' x <- sin(seq(0, 2 * pi, length.out = 100))
#' y <- cos(seq(0, 2 * pi, length.out = 100))
#' # Compute instantaneous correlation using Rcpp
#' rho <- instantaneous_correlation_rcpp(x, y, eta = 0.1)
#' # Plot the result
#' plot(rho, type = 'l', ylab = 'Instantaneous Correlation', xlab = 'Time')
#' @export
instantaneous_correlation_rcpp <- function(x, y, eta, offset = 0) {
  # Ensure inputs are numeric vectors of the same length
  if (!is.numeric(x) || !is.numeric(y)) {
    stop("x and y must be numeric vectors.")
  }
  N <- length(x)
  if (length(y) != N) {
    stop("x and y must be of the same length.")
  }
  #.Call('_instantaneous_correlation_rcpp', x, y, eta, offset)
}


### see roll_cor


## we were going to add dynamic psychophysiological interactions
## also, ppi between all brain regions

# dyncor <- function(x,y, window=7, symmetric=FALSE) {
#   assertthat::assert_that(length(x) == length(y))
#   halfwin <- as.integer(window/2)
#   furrr::future_map_dbl(1:length(x), function(i) {
#     start <- max(i-window, 1)
#     if (i > halfwin) {
#       cor(x[start:i], y[start:i])
#     } else {
#       NA
#     }
#     
#   })
# }
  
#}




