

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




