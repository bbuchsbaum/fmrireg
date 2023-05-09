#' @import splines
NULL

#' Soft-threshold function
#'
#' This function applies soft-thresholding to the input values, setting values below the threshold to zero
#' and shrinking the remaining values by the threshold amount.
#'
#' @param x A numeric vector of input values
#' @param threshold A non-negative threshold value for the soft-thresholding operation
#'
#' @return A numeric vector with the soft-thresholded values
#'
#' @examples
#' x <- c(-3, -2, -1, 0, 1, 2, 3)
#' threshold <- 1
#' soft_thresholded <- soft_threshold(x, threshold)
#' print(soft_thresholded)
#'
soft_threshold <- function(x, threshold) {
  if (threshold < 0) {
    stop("Threshold value should be non-negative.")
  }
  
  sign(x) * pmax(0, abs(x) - threshold)
}

#' Construct an HRF Instance
#' 
#' @description
#' `gen_hrf` takes a raw function `f(t)` and returns an HRF (Hemodynamic Response Function) instance.
#' 
#' @param hrf A function mapping from time to signal.
#' @param lag Optional lag in seconds.
#' @param width Optional block width in seconds.
#' @param precision Sampling precision in seconds.
#' @param summate Whether to allow each impulse response function to "add" up (default: TRUE).
#' @param normalize Rescale so that the peak of the output is 1 (default: FALSE).
#' @param name The assigned name of the generated HRF.
#' @param span The span of the HRF (maximum width in seconds after which function reverts to zero).
#' @param ... Extra parameters for the `hrf` function.
#' 
#' @return An instance of type `HRF` inheriting from `function`.
#' 
#' @examples 
#' 
#' ## Generate an HRF using SPMG1 canonical HRF, a lag of 3, and a width of 2.
#' grf <- gen_hrf(hrf_spmg1, lag=3, width=2)
#' grf(0:20)
#' 
#' hg <- purrr::partial(hrf_gaussian, mean=4, sd=1)
#' grf <- gen_hrf(hg, lag=1, width=2)
#' 
#' vals <- grf(0:20)
#' @export
gen_hrf <- function(hrf, lag=0, width=0, precision=.1, 
                    summate=TRUE, normalize=FALSE, name="gen_hrf", span=NULL, ...) {
  .orig <- list(...)
  
  if (width != 0) {
    assert_that(width > 0)
    #hrf <- gen_hrf_blocked(hrf, width=width, precision=precision, 
    #                       summate=summate, normalize=normalize, ...)
    hrf <- gen_hrf_blocked(hrf, width=width, precision=precision, summate=summate, normalize=normalize)
  }
  
  if (lag !=0) {
    hrf <- gen_hrf_lagged(hrf, lag=lag)
  }
  
  f <- if (length(.orig) > 0) {
    ret <- function(t) {
      do.call(hrf, c(list(t), .orig))
    }
    attr(ret, "params") <- .orig
    ret
  } else {
    hrf
  }
  
  vals <- f(0:2)
  
  nb <- if (is.vector(vals)) {
    1
  } else if (is.matrix(vals)) {
    ncol(vals)
  } else {
    stop("gen_hrf: constructed hrf is invalid")
  }
  
  if (is.null(span)) {
    span <- 16 + lag + (width*2)
  }
  HRF(f, name=name, nbasis=nb, span=span)
}


#' Generate an Empirical Hemodynamic Response Function
#' 
#' @description
#' `gen_empirical_hrf` generates an empirical hemodynamic response function (HRF) using provided time points and HRF values.
#' 
#' @param t Time points.
#' @param y Values of HRF at time `t[i]`.
#' @param name Name of the generated HRF.
#' 
#' @return An instance of type `HRF` inheriting from `function`.
#' 
#' @examples 
#' 
#' y <- -poly(0:24, 2)[,2]
#' emphrf <- gen_empirical_hrf(0:24, y)
#' ## plot(emphrf(seq(0,24,by=.5)), type='l')
#' @export
gen_empirical_hrf <- function(t, y, name="empirical_hrf") {
  f <- approxfun(t,y, yright=0, yleft=0)
  HRF(f, name=name, nbasis=1)
}



#' Generate an HRF Basis Set
#' 
#' @description
#' `gen_hrf_set` constructs an HRF basis set from one or more component functions.
#' This function is used to create arbitrary sets of basis functions for fMRI modeling.
#' 
#' @param ... A list of functions f(t) mapping from time to amplitude.
#' @param span The span in seconds of the HRF.
#' @param name The name of the HRF.
#' 
#' @return An instance of type `HRF` inheriting from `function`.
#' 
#' @family gen_hrf
#' 
#' @examples 
#' 
#' hrf1 <- hrf_spmg1 |> gen_hrf(lag=0)
#' hrf2 <- hrf_spmg1 |> gen_hrf(lag=3)
#' hrf3 <- hrf_spmg1 |> gen_hrf(lag=6)
#' 
#' hset <- gen_hrf_set(hrf1, hrf2, hrf3)
#' @export
gen_hrf_set <- function(..., span=32, name="hrf_set") {
  hrflist <- list(...)
  assertthat::assert_that(all(sapply(hrflist, is.function)))
  f <- function(t) {
    do.call("cbind", lapply(hrflist, function(fun) fun(t)))
  }
  
  HRF(f, name=name, span=span, nbasis=length(hrflist))
}



#' Generate an HRF (hemodynamic response function) library
#'
#' This internal function generates an HRF library by applying the provided HRF generating function (`fun`) to each row of the parameter grid (`pgrid`). Additional arguments can be passed to the generating function using `...`.
#'

#' @param fun A function that generates an HRF, given a set of parameters.
#' @param pgrid A data frame containing the parameter grid, where each row represents a set of parameters to be passed to the HRF generating function.
#' @param ... Additional arguments to be passed to the HRF generating function.
#' @family gen_hrf
#' @keywords internal
#' @return A list of HRFs generated by applying the HRF generating function to each row of the parameter grid.
gen_hrf_library <- function(fun, pgrid,...) {
  pnames <- names(pgrid)
  
  hrflist <- lapply(1:nrow(pgrid), function(i) {
    do.call(gen_hrf, c(fun, pgrid[i,],...))
  })
  
  do.call(gen_hrf_set, hrflist)

}


#' HRF Constructor Function
#'
#' @description
#' The `HRF` function creates an object representing a hemodynamic response function (HRF). It is a class constructor for HRFs.
#'
#' @param fun A function representing the hemodynamic response, mapping from time to BOLD response.
#' @param name A string specifying the name of the function.
#' @param nbasis An integer representing the number of basis functions, e.g., the columnar dimension of the HRF. Default is 1.
#' @param span A numeric value representing the span in seconds of the HRF. Default is 24.
#' @param param_names A character vector containing the names of the parameters for the HRF function.
#'
#' @return An HRF object with the specified properties.
#'
#' @examples
#' hrf <- HRF(hrf_gamma, "gamma", nbasis=1, param_names=c("shape", "rate"))
#' resp <- evaluate(hrf, seq(0, 24, by=1))
#'
#' @export
HRF <- function(fun, name, nbasis=1, span=24, param_names=NULL) {
  vals <- fun(seq(0,span))

  if (nbasis == 1) {
    peak <- max(vals, na.rm=TRUE)
  } else {
    peak <- max(apply(vals, 2, max, na.rm=TRUE))
  }
  
  scale_factor <- 1/peak
  
  
  structure(fun, name=name, 
            nbasis=as.integer(nbasis), 
            span=span,
            param_names=param_names, 
            scale_factor=scale_factor, 
            class=c("HRF", "function"))
  
}

#' AFNI HRF Constructor Function
#'
#' @description
#' The `AFNI_HRF` function creates an object representing an AFNI-specific hemodynamic response function (HRF). It is a class constructor for AFNI HRFs.
#'
#' @param name A string specifying the name of the AFNI HRF.
#' @param nbasis An integer representing the number of basis functions for the AFNI HRF.
#' @param params A list containing the parameter values for the AFNI HRF.
#'
#' @return An AFNI_HRF object with the specified properties.
#'
#' @seealso HRF
#'
#' @inheritParams HRF
#' @describeIn HRF-class AFNI hrf
#' @export
AFNI_HRF <- function(name, nbasis, params) {
  structure(name, 
            nbasis=as.integer(nbasis), 
            params=params, 
            class=c("AFNI_HRF"))
  
}


#' @export
as.character.AFNI_HRF <- function(x,...) {
  paste(x, "\\(", paste(attr(x, "params"), collapse=","), "\\)", sep="")
}



#' @importFrom numDeriv grad
#' @keywords internal
#' @noRd
makeDeriv <- function(HRF, n=1) {
  if (n == 1) {
    function(t) numDeriv::grad(HRF, t)
  } else {
    Recall(function(t) numDeriv::grad(HRF,t), n-1)
  }
}


#' Generate a Lagged HRF Function
#'
#' @description
#' The `gen_hrf_lagged` function takes an HRF function and applies a specified lag to it. This can be useful for modeling time-delayed hemodynamic responses.
#'
#' @param hrf A function representing the underlying HRF to be shifted.
#' @param lag A numeric value specifying the lag or delay in seconds to apply to the HRF. This can also be a vector of lags, in which case the function returns an HRF set.
#' @param normalize A logical value indicating whether to rescale the output so that the maximum absolute value is 1. Defaults to `FALSE`.
#' @param ... Extra arguments supplied to the `hrf` function.
#'
#' @return A function representing the lagged HRF. If `lag` is a vector of lags, the function returns an HRF set.
#' @family gen_hrf
#' @examples
#' hrf_lag5 <- gen_hrf_lagged(HRF_SPMG1, lag=5)
#' hrf_lag5(0:20)
#'
#' @export
gen_hrf_lagged <- function(hrf, lag=2, normalize=FALSE, ...) {
  force(hrf)
  # TODO deal with nbasis arg in ...
  if (length(lag)>1) {
    do.call(gen_hrf_set, lapply(lag, function(l) gen_hrf_lagged(hrf, l,...)))
  } else {
    function(t) {
      ret <- hrf(t-lag,...)
      if (normalize) {
        ret <- ret/max(abs(ret))
      } 
      
      ret
    }
  }
}

#' @export
#' @describeIn gen_hrf_lagged alias for gen_hrf_lagged
#' @family gen_hrf
hrf_lagged <- gen_hrf_lagged


#' Generate a Blocked HRF Function
#'
#' @description
#' The `gen_hrf_blocked` function creates a blocked HRF by convolving the input HRF with a boxcar function. This can be used to model block designs in fMRI analysis.
#'
#' @param hrf A function representing the hemodynamic response function. Default is `hrf_gaussian`.
#' @param width A numeric value specifying the width of the block in seconds. Default is 5.
#' @param precision A numeric value specifying the sampling resolution in seconds. Default is 0.1.
#' @param half_life A numeric value specifying the half-life of the exponential decay function, used to model response attenuation. Default is `Inf`, which means no decay.
#' @param summate A logical value indicating whether to allow each impulse response function to "add" up. Default is `TRUE`.
#' @param normalize A logical value indicating whether to rescale the output so that the peak of the output is 1. Default is `FALSE`.
#' @param ... Extra arguments passed to the HRF function.
#' @family gen_hrf
#'
#' @return A function representing the blocked HRF.
#'
#' @importFrom purrr partial
#' @export
gen_hrf_blocked <- function(hrf=hrf_gaussian, width=5, precision=.1, 
                            half_life=Inf, summate=TRUE, normalize=FALSE, ...) {
  force(hrf)
  purrr::partial(convolve_block, hrf=hrf, width=width, 
                 precision=precision, half_life=half_life, 
                 summate=summate, normalize=normalize, ...)
}

#' @export
#' @aliases gen_hrf_blocked
#' @describeIn gen_hrf_blocked alias for gen_hrf_blocked
hrf_blocked <- gen_hrf_blocked


#' Convolve hemodynamic response with a block duration
#'
#' This function convolves a hemodynamic response function (HRF) with a block duration, producing a time-varying response 
#' for the specified duration. The function supports various HRFs, block widths, precision levels, half-lives, and normalization.
#'
#' @param t A vector of time points in seconds at which the convolved response will be evaluated
#' @param hrf The hemodynamic response function to be convolved, provided as a function (default: hrf_gaussian)
#' @param width The fixed width of the response block in seconds (default: 5)
#' @param precision The sampling precision of the response in seconds (default: 0.2)
#' @param half_life The half-life of the exponential decay function in seconds, used to model attenuation (default: Inf)
#' @param summate A logical value indicating whether to sum the impulse response functions (default: TRUE)
#'                If FALSE, the function returns the maximum value at each time point.
#' @param normalize A logical value indicating whether to normalize the output so that its peak value is 1 (default: FALSE)
#' @param ... Additional arguments passed to the HRF function
#'
#' @return A vector of convolved responses at the specified time points
#'
#' @examples
#' time_points <- seq(0, 20, by = 0.5)
#' block_response <- convolve_block(time_points, hrf=hrf_spmg1, width=5, precision=0.2)
#' plot(time_points, block_response, type='l', main='Convolved Hemodynamic Response with Block Duration')
#'
#' @export
convolve_block <- function(t, hrf=hrf_gaussian, width=5, precision=.2, half_life=Inf, summate=TRUE, normalize=FALSE, ...) {
 
  hmat <- sapply(seq(0, width, by=precision), function(offset) {
    hrf(t-offset,...) * exp(-offset/half_life)
  })
  
  ret <- if (summate) {
    rowSums(hmat)
  } else {
    #r <- range(hmat[,1])
    #apply(hmat,1,function(vals) vals[which.max(abs(vals))])
    apply(hmat,1,function(vals) vals[which.max(vals)])
  }
  
  if (normalize) {
    ret <- ret/max(abs(ret))
  } 
  
  ret
}



# 
# find_best_cor <- function(y,ys) {
#   w =waveslim::dwt(y, wf="la16", n.levels=4)
#   g <- expand.grid(maxlevel=1:4, value=seq(.1,3,by=.2), hard=c(TRUE, FALSE))
#   ret <- lapply(1:nrow(g), function(i) {
#     ws=manual.thresh(w, max.level=g$maxlevel[i], value=g$value[i], hard=g$hard[i])
#     recon <- waveslim::idwt(ws)
#     cor(recon,ys)
#   })
#   g$cor = unlist(ret)
# 
# 
# }

# find_best_cor_dct <- function(y,ys) {
#   w =dct(y)
#   g <- expand.grid(value=seq(.01,.6,by=.01), hard=c(TRUE, FALSE))
#   ret <- lapply(1:nrow(g), function(i) {
#     wt <- w
#     if (g$hard[i]) {
#       wt[abs(wt) < g$value[i]] = 0
#     } else {
#       wt <- soft_threshold(wt, g$value[i])
#     }
#     recon <- idct(wt)
#     cor(recon,ys)
#   })
#   g$cor = unlist(ret)
# 
# }


# library(compiler)
# library(complex)
# 
# mdct4 <- function(x) {
#   N <- length(x)
#   if (N %% 4 != 0) {
#     stop("MDCT4 only defined for vectors of length multiple of four.")
#   }
#   M <- N %/% 2
#   N4 <- N %/% 4
#   
#   rot <- c(tail(x, N4), head(x, -N4))
#   rot[1:N4] <- -rot[1:N4]
#   t <- 0:(N4-1)
#   w <- exp(-1i * 2 * pi * (t + 1 / 8) / N)
#   c <- rot[2 * t + 1] - rot[N - 2 * t] - 1i * (rot[M + 2 * t + 1] - rot[M - 2 * t])
#   c <- (2 / sqrt(N)) * w * fft(0.5 * c * w, N4)
#   y <- numeric(M)
#   y[2 * t + 1] <- Re(c[t + 1])
#   y[M - 2 * t] <- -Im(c[t + 1])
#   return(y)
# }
# 
# imdct4 <- function(x) {
#   N <- length(x)
#   if (N %% 2 != 0) {
#     stop("iMDCT4 only defined for even-length vectors.")
#   }
#   M <- N %/% 2
#   N2 <- N * 2
#   
#   t <- 0:(M - 1)
#   w <- exp(-1i * 2 * pi * (t + 1 / 8) / N2)
#   c <- x[2 * t + 1] + 1i * x[N - 2 * t + 1]
#   c <- 0.5 * w * c
#   c <- fft(c, M)
#   c <- ((8 / sqrt(N2)) * w) * c
#   
#   rot <- numeric(N2)
#   rot[2 * t + 1] <- Re(c[t + 1])
#   rot[N + 2 * t + 1] <- Im(c[t + 1])
#   
#   t <- seq(1, N2, by = 2)
#   rot[t] <- -rot[N2 - t + 1]
#   
#   t <- 0:(3 * M - 1)
#   y <- numeric(N2)
#   y[t + 1] <- rot[t + M + 1]
#   t <- (3 * M):(N2 - 1)
#   y[t + 1] <- -rot[t - 3 * M + 1]
#   return(y)
# }
# 



#' HRF (hemodynamic response function) as a linear function of time
#'
#' The `hrf_time` function computes the value of an HRF, which is a simple linear function of time `t`, when `t` is greater than 0 and less than `maxt`.
#'
#' @param t A numeric value representing time in seconds.
#' @param maxt A numeric value representing the maximum time point in the domain. Default value is 22.
#' @return A numeric value representing the value of the HRF at the given time `t`.
#' @family hrf_functions
#' @export
#' @examples
#' # Compute the HRF value for t = 5 seconds with the default maximum time
#' hrf_val <- hrf_time(5)
#'
#' # Compute the HRF value for t = 5 seconds with a custom maximum time of 30 seconds
#' hrf_val_custom_maxt <- hrf_time(5, maxt = 30)
hrf_time <- function(t, maxt=22) {
  ifelse(t > 0 & t < maxt, t, 0)
}

# hrf_ident
# 
# @param t time in seconds
# @export
#hrf_ident <- function(t) {
#  ifelse( t == 0, 1, 0)
#}

#' B-spline HRF (hemodynamic response function)
#'
#' The `hrf_bspline` function computes the B-spline representation of an HRF (hemodynamic response function) at given time points `t`.
#'
#' @param t A vector of time points.
#' @param span A numeric value representing the temporal window over which the basis set spans. Default value is 20.
#' @param N An integer representing the number of basis functions. Default value is 5.
#' @param degree An integer representing the degree of the spline. Default value is 3.
#' @return A matrix representing the B-spline basis for the HRF at the given time points `t`.
#' @family hrf_functions
#' @examples
#' # Compute the B-spline HRF representation for time points from 0 to 20 with 0.5 increments
#' hrfb <- hrf_bspline(seq(0, 20, by = .5), N = 4, degree = 2)
#' @export
#' @importFrom splines bs
hrf_bspline <- function(t, span=20, N=5, degree=3) {
	
	ord <- 1 + degree
	nIknots <- N - ord + 1
	if (nIknots < 0) {
		nIknots <- 0
		#warning("'df' was too small; have used  ", ord - (1 - intercept))
	}
	
	knots <- if (nIknots > 0) {
				knots <- seq.int(from = 0, to = 1, length.out = nIknots + 2)[-c(1, nIknots + 2)]
				stats::quantile(seq(0,span), knots)
			} else {
				0
			}
	
	if (any(t < 0)) {
		t[t < 0] <- 0
	}
	
	if(any(t > span)) {
		t[t > span] <- 0
	}
	
	splines::bs(t, df=N, knots=knots, degree=degree, Boundary.knots=c(0,span))
}


#' Gamma HRF (hemodynamic response function)
#'
#' The `hrf_gamma` function computes the gamma density-based HRF (hemodynamic response function) at given time points `t`.
#'
#' @param t A vector of time points.
#' @param shape A numeric value representing the shape parameter for the gamma probability density function. Default value is 6.
#' @param rate A numeric value representing the rate parameter for the gamma probability density function. Default value is 1.
#' @return A numeric vector representing the gamma HRF at the given time points `t`.
#' @family hrf_functions
#' @examples
#' # Compute the gamma HRF representation for time points from 0 to 20 with 0.5 increments
#' hrf_gamma <- hrf_gamma(seq(0, 20, by = .5), shape = 6, rate = 1)
#' @export
hrf_gamma <- function(t, shape=6, rate=1) {
  stats::dgamma(t, shape=shape, rate=rate)
}

#' Gaussian HRF (hemodynamic response function)
#'
#' The `hrf_gaussian` function computes the Gaussian density-based HRF (hemodynamic response function) at given time points `t`.
#'
#' @param t A vector of time points.
#' @param mean A numeric value representing the mean of the Gaussian probability density function. Default value is 6.
#' @param sd A numeric value representing the standard deviation of the Gaussian probability density function. Default value is 2.
#' @return A numeric vector representing the Gaussian HRF at the given time points `t`.
#' @family hrf_functions
#' @examples
#' # Compute the Gaussian HRF representation for time points from 0 to 20 with 0.5 increments
#' hrf_gaussian <- hrf_gaussian(seq(0, 20, by = .5), mean = 6, sd = 2)
#' @export
hrf_gaussian <- function(t, mean=6, sd=2) {
	dnorm(t, mean=mean, sd=sd)
}



#' Mexican Hat HRF (hemodynamic response function)
#'
#' The `hrf_mexhat` function computes the Mexican hat wavelet-based HRF (hemodynamic response function) at given time points `t`.
#'
#' @param t A vector of time points.
#' @param mean A numeric value representing the mean of the Mexican hat wavelet. Default value is 6.
#' @param sd A numeric value representing the standard deviation of the Mexican hat wavelet. Default value is 2.
#' @return A numeric vector representing the Mexican hat wavelet-based HRF at the given time points `t`.
#' @family hrf_functions
#' @examples
#' # Compute the Mexican hat HRF representation for time points from 0 to 20 with 0.5 increments
#' hrf_mexhat <- hrf_mexhat(seq(0, 20, by = .5), mean = 6, sd = 2)
#' @export
hrf_mexhat <- function(t, mean = 6, sd = 2) {
  t0 <- t - mean
  a <- (1 - (t0 / sd)^2) * exp(-t0^2 / (2 * sd^2))
  scale <- sqrt(2 / (3 * sd * pi^(1/4)))
  return(scale * a)
}

#' hrf_spmg1
#'
#' A hemodynamic response function based on the SPM canonical double gamma parameterization.
#'
#' This function models the hemodynamic response using the canonical double gamma parameterization
#' in the SPM software. The HRF is defined by a linear combination of two gamma functions with
#' different exponents (P1 and P2) and amplitudes (A1 and A2). It is commonly used in fMRI data
#' analysis to estimate the BOLD (blood-oxygen-level-dependent) signal changes associated with
#' neural activity.
#'
#' @param t A vector of time points.
#' @param P1 The first exponent parameter (default: 5).
#' @param P2 The second exponent parameter (default: 15).
#' @param A1 Amplitude scaling factor for the positive gamma function component; normally fixed at .0833
#' @return A vector of HRF values at the given time points.
#' @family hrf_functions
#' @export
#' @examples
#' # Generate a time vector
#' time_points <- seq(0, 30, by=0.1)
#' # Compute the HRF values using the SPM canonical double gamma parameterization
#' hrf_values <- hrf_spmg1(time_points)
#' # Plot the HRF values
#' plot(time_points, hrf_values, type='l', main='SPM Canonical Double Gamma HRF')
hrf_spmg1 <- function(t, P1=5, P2=15,A1=.0833) {
 	ifelse(t < 0, 0, exp(-t)*(A1*t^P1 - 1.274527e-13*t^P2))
	
}

# objective_function <- function(params, desired_peak, desired_fwhm) {
#   P1 <- params[1]
#   P2 <- params[2]
#   hrf <- function(t) hrf_spmg1(t, P1, P2)
#   t_values <- seq(0, 30, by=0.1)
#   hrf_values <- sapply(t_values, hrf)
#   
#   # Calculate actual peak
#   actual_peak <- t_values[which.max(hrf_values)]
#   
#   # Calculate actual FWHM
#   half_max <- max(hrf_values) / 2
#   above_half_max <- t_values[hrf_values > half_max]
#   actual_fwhm <- max(above_half_max) - min(above_half_max)
#   
#   # Calculate the error based on both desired_peak and desired_fwhm
#   peak_error <- abs(actual_peak - desired_peak)
#   fwhm_error <- abs(actual_fwhm - desired_fwhm)
#   
#   return(peak_error + fwhm_error)
# }

#' @keywords internal
HRF_GAMMA <- HRF(hrf_gamma, "gamma", param_names=c("shape", "rate"))

#' @export
#' @describeIn HRF-class Gaussian hrf
HRF_GAUSSIAN <- HRF(hrf_gaussian, name="gaussian", param_names=c("mean", "sd"))

#' @export
#' @describeIn HRF-class B-spline hrf
HRF_BSPLINE <- HRF(gen_hrf(hrf_bspline), name="bspline", nbasis=5)

# @export
# @rdname HRF
# HRF_IDENT <- HRF(gen_hrf(hrf_ident), "ident", nbasis=1)


#' @keywords internal
#' @describeIn HRF-class SPMG1 hrf
#' @export 
HRF_SPMG1 <- HRF(hrf_spmg1, 
                 "SPMG1", param_names=c("A1", "A2"))


#' @describeIn HRF-class SPMG2 hrf
#' @export
HRF_SPMG2 <- HRF(gen_hrf_set(hrf_spmg1, makeDeriv(hrf_spmg1)), 
                 "SPMG2", nbasis=2, param_names=c("A1", "A2"))

#' @describeIn HRF-class SPMG3 hrf
#' @export
HRF_SPMG3 <- HRF(gen_hrf_set(hrf_spmg1, makeDeriv(hrf_spmg1), makeDeriv(makeDeriv(hrf_spmg1))), 
                 "SPMG3", nbasis=3, param_names=c("A1", "A2"))


#' evaluate.HRF
#'
#' This function evaluates a hemodynamic response function (HRF) for a given set of time points (grid) and other parameters.
#' It is used to generate the predicted BOLD (blood-oxygen-level-dependent) signal changes in fMRI data analysis.
#'
#' @param x The HRF function.
#' @param grid A vector of time points.
#' @param amplitude The scaling value for the event (default: 1).
#' @param duration The duration of the event (default: 0).
#' @param precision The temporal resolution used for computing summed responses when duration > 0 (default: 0.2).
#' @param summate Whether the HRF response increases its amplitude as a function of stimulus duration (default: TRUE).
#' @param normalize Scale output so that peak is 1 (default: FALSE).
#' @return A vector of HRF values at the specified time points.
#' @examples
#' hrf1 <- evaluate(HRF_SPMG1, grid=seq(0,20,by=1.5), duration=2, precision=.1)
#' hrf2 <- evaluate(HRF_SPMG1, grid=seq(0,20,by=1.5), duration=2, precision=.1, summate=FALSE)
#' @export
evaluate.HRF <- function(x, grid, amplitude=1, duration=0, precision=.2, summate=TRUE, normalize=FALSE, ...) {
  if (duration < precision) {
    if (normalize) {
      x(grid)*amplitude*attr(x, "scale_factor")   
    } else {
      x(grid)*amplitude
    }
  } else if (nbasis(x) == 1) {
    samples <- seq(0, duration, by=precision)
    #sfac <- attr(x, "scale_factor") 
    
    hmat <- sapply(samples, function(offset) {
      x(grid-offset)*amplitude
    })
    
    ret <- if (summate) {
      rowSums(hmat)
    } else {
      #rowMeans(hmat)
      apply(hmat,1,function(vals) vals[which.max(vals)])
    }
    
    if (normalize) {
      ret <- ret/max(abs(ret))
    }
    
    ret
  } else {
    mat <- Reduce("+", lapply(seq(0, duration, by=precision), function(offset) {
      x(grid-offset)*amplitude   
    }))
    
    if (normalize) {
      mat <- apply(mat, 2, function(vals) vals/max(abs(vals)))
    }
    
    mat
  }
}


#' evaluate.hrfspec
#'
#' This function evaluates a hemodynamic response function (HRF) specified by an hrfspec object for a given set of time points (grid) and other parameters.
#' It is a wrapper function that calls the evaluate.HRF function with the HRF function contained in the hrfspec object.
#'
#' @param x The hrfspec object containing the HRF function.
#' @param grid A vector of time points.
#' @param amplitude The scaling value for the event (default: 1).
#' @param duration The duration of the event (default: 0).
#' @param precision The temporal resolution used for computing summed responses when duration > 0 (default: 0.1).
#' @param ... Additional arguments to be passed to the evaluate.HRF function.
#' @return A vector of HRF values at the specified time points.
#' @examples
#' hrf_spec <- hrfspec(hrf = HRF_SPMG1)
#' hrf_values <- evaluate(hrfspec, grid=seq(0,20,by=1.5), duration=2, precision=.1)
#' @export
evaluate.hrfspec <- function(x, grid, amplitude=1, duration=0, precision=.1, ...) {
  evaluate(x$hrf, grid,amplitude, duration, precision)
}


#' @export
nbasis.HRF <- function(x) attr(x, "nbasis")

#' getHRF
#'
#' This function retrieves a specified hemodynamic response function (HRF) by name and creates an HRF object with specified properties.
#'
#' @param name The name of the HRF function. Available options: "gamma", "spmg1", "spmg2", "spmg3", "bspline", "gaussian", "tent", "bs".
#' @param nbasis The number of basis functions (if relevant, default: 5).
#' @param span The temporal window over which the basis sets span (default: 24).
#' @param lag The time lag parameter (default: 0).
#' @param width The width parameter (default: 0).
#' @param summate Whether the HRF response increases its amplitude as a function of stimulus duration (default: TRUE).
#' @param normalize Whether to scale output so that the peak is 1 (default: FALSE).
#' @param ... Additional arguments passed to the gen_hrf function.
#' @return An HRF object with the specified properties.
#' @keywords internal
#' @examples
#' hrf_obj <- getHRF(name = "gamma", nbasis = 5, span = 24)
getHRF <- function(name=c("gam", "gamma", "spmg1", "spmg2", 
                          "spmg3", "bspline", "gaussian", "tent", "bs"), 
                   nbasis=5, span=24,lag=0,width=0, summate=TRUE, normalize=FALSE, ...) {
  name <- match.arg(name)
	nb <- nbasis
	hrf <- switch(name,
			gamma=gen_hrf(hrf_gamma,lag=lag, span=span,width=width, summate=summate, normalize=normalize, name="gamma"),
			gam=gen_hrf(hrf_gamma,lag=lag,span=span,width=width, summate=summate, normalize=normalize, name="gamma"),
			gaussian=gen_hrf(HRF_GAUSSIAN,lag=lag,span=span, width=width, summate=summate, normalize=normalize, name="gaussian"),
			spmg1=gen_hrf(HRF_SPMG1,lag=lag,span=span, width=width, summate=summate, normalize=normalize, name="spmg1"),
			spmg2=gen_hrf(HRF_SPMG2,lag=lag,span=span, width=width, summate=summate, normalize=normalize, name="spmg2"),
			spmg3=gen_hrf(HRF_SPMG3,lag=lag,span=span, width=width, summate=summate, normalize=normalize, name="spmg3"),
			tent=gen_hrf(purrr::partial(hrf_bspline, degree=1, span=span, N=nb), lag=lag, width=width, summate=summate, normalize=normalize, name="tent", ...),
			bs=gen_hrf(purrr::partial(hrf_bspline,N=nb,degree=3,span=span), lag=lag, width=width, summate=summate, normalize=normalize, name="bspline", ...),
			bspline=gen_hrf(purrr::partial(hrf_bspline, N=nb, degree=3, span=span), lag=lag,width=width, summate=summate, normalize=normalize, name="bspline", ...)
	)
	
	if (is.null(hrf)) {
		stop("could not find create hrf named: ", name)
	}
	
	hrf
}

#' hrf_sine
#'
#' A hemodynamic response function using the Sine Basis Set.
#'
#' @param t A vector of times.
#' @param span The temporal window over which the basis sets span (default: 24).
#' @param N The number of basis functions (default: 5).
#' @return A matrix of sine basis functions.
#' @family hrf_functions
#' @export
#' @examples
#' hrf_sine_basis <- hrf_sine(seq(0, 20, by = 0.5), N = 4)
hrf_sine <- function(t, span = 24, N = 5) {
  sine_basis <- sapply(1:N, function(n) {
    sin(2 * pi * n * t / span)
  })
  return(sine_basis)
}

#' hrf_inv_logit
#'
#' A hemodynamic response function using the difference of two Inverse Logit functions.
#'
#' @param t A vector of times.
#' @param mu1 The time-to-peak for the rising phase (mean of the first logistic function).
#' @param s1 The width (slope) of the first logistic function.
#' @param mu2 The time-to-peak for the falling phase (mean of the second logistic function).
#' @param s2 The width (slope) of the second logistic function.
#' @param lag The time delay (default: 0).
#' @return A vector of the difference of two Inverse Logit HRF values.
#' @family hrf_functions
#' @export
#' @examples
#' hrf_inv_logit_basis <- hrf_diff_inv_logit(seq(0, 20, by = 0.5), mu1 = 6, s1 = 1, mu2 = 16, s2 = 1)
hrf_inv_logit <- function(t, mu1 = 6, s1 = 1, mu2 = 16, s2 = 1, lag = 0) {
  inv_logit1 <- 1 / (1 + exp(-(t - lag - mu1) / s1))
  inv_logit2 <- 1 / (1 + exp(-(t - lag - mu2) / s2))
  return(inv_logit1 - inv_logit2)
}


#' @keywords internal
make_hrf <- function(basis, lag, nbasis=1) {
  if (!is.numeric(lag) || length(lag) > 1) {
    stop("hrf: 'lag' must be a numeric scalar")
  }
  
  if (is.character(basis)) {
    getHRF(basis, nbasis=nbasis, lag=lag)
  } else if (inherits(basis, "HRF")) {
    if (lag > 0) {
      HRF(gen_hrf_lagged(basis, lag=lag), name=basis$name, nbasis=basis$nbasis)
    } else {
      basis
    }
    
  } else if (is.function(basis)) {
    test <- basis(1:10)
    nb <- if (is.vector(test)) {
      1
    } else {
      assert_that(is.matrix(test), "basis function must return vector or matrix")
      ncol(test) 
    }
    HRF(gen_hrf_lagged(basis,lag=lag), name="custom_hrf", nbasis=nb)
  } else {
    stop("invalid basis function: must be 1) character string indicating hrf type, e.g. 'gamma' 2) a function or 3) an object of class 'HRF': ", basis)
  }
  
}

#### TODO character variables need an "as.factor"



#' hemodynamic regressor specification function for model formulas.
#' 
#' This function is to be used in formulas for fitting functions, e.g. onsets ~ hrf(fac1,fac2) ...
#' 
#' 
#' @param ... the variable names, all of which must be present in the enclosing environment (e.g. an \code{event_model} object)
#' @param basis the impulse response function or the name of a pre-supplied function, one of: "gamma", "spmg1", "spmg2", "spmg3", "bspline", "gaussian".
#' @param onsets optional onsets override. If missing, onsets will be taken from the \code{event_model}
#' @param durations optional durations override. If missing, onsets will be taken from the \code{event_model}
#' @param prefix a character string that is prepended to the variable names and used to identify the term. 
#'               Can be used to disambiguate two \code{hrf} terms with the same variable(s) but different onsets or basis functions.
#' @param subset an expression indicating the subset of 'onsets' to keep
#' @param precision sampling precision in seconds
#' @param nbasis number of basis functions -- only used for hemodynamic response functions (e.g. bspline) that take a variable number of bases.
#' @param contrasts one or more \code{contrast_spec} objects created with the \code{contrast} function. 
#' If multiple contrasts are required, then these should be wrapped in a \code{list} or \code{contrast_set}.
#' @param id a  unique \code{character} identifier used to refer to term, otherwise will be determined from variable names.
#' @param lag a temporal offset in seconds which is added to onset before convolution
#' @param summate whether impulse amplitudes sum up when duration is greater than 0. 
#' @examples 
#' 
#' ## 'hrf' is typically used in the context of \code{formula}s.
#' 
#' form <- onsets ~ hrf(x) + hrf(y) + hrf(x,y)
#' 
#' @export
#' @importFrom rlang enquos enexpr
hrf <- function(..., basis="spmg1", onsets=NULL, durations=NULL, prefix=NULL, subset=NULL, precision=.3, 
                nbasis=1, contrasts=NULL, id=NULL, lag=0, summate=TRUE) {
  
  vars <- rlang::enquos(...)
  #vars <- as.list(substitute(...))[-1] 
  
  parsed <- parse_term(vars, "hrf")
  term <- parsed$term
  label <- parsed$label
  
  basis <- make_hrf(basis, lag, nbasis=nbasis)
  
  ret <- hrfspec(
    term,
    label=label,
    basis,
    basis=basis,         ## basis function of type "HRF"
    onsets=onsets,       ## optional onsets vector -- should this be lazy? onsets = ~ onsets
    durations=durations, ## optional durations vector
    prefix=prefix,       ## prefix
    subset=rlang::enexpr(subset), ## quoted subset expression
    precision=precision,
    contrasts=contrasts,
    summate=summate)
  
  class(ret) <- c("hrfspec", "list")
  ret
}


#' @keywords internal
hrfspec <- function(vars, label=NULL, basis=HRF_SPMG1, onsets=NULL, durations=NULL, prefix=NULL, 
                    subset=NULL, precision=.3, 
                    contrasts=NULL, id=NULL, summate=TRUE) {
  
  
  
  assert_that(inherits(basis, "HRF"))
  termname <- paste0(vars, collapse="::")
  
  varnames <- if (!is.null(prefix)) {
    paste0(prefix, "_", vars)
  } else {
    vars
  }
  
  if (is.null(id)) {
    id <- termname
  }
  
  if (is.null(label)) {
    label <- paste0("hrf(", paste0(varnames, collapse=","), ")")
  }
  
  cset <- if (inherits(contrasts, "contrast_spec")) {
    contrast_set(con1=contrasts)
  } else if (inherits(contrasts, "contrast_set")) {
    contrasts
  } 
  
  
  ret <- list(
    name=termname, 
    label=label,
    id=id, 
    vars=vars,
    varnames=varnames, 
    hrf=basis,
    onsets=onsets,
    durations=durations,
    prefix=prefix,
    subset=subset,
    precision=precision,
    contrasts=cset,
    summate=summate)
  
  class(ret) <- c("hrfspec", "list")
  ret
  
}


#' @keywords internal
construct_event_term <- function(x, model_spec, onsets) {
  
  ## TODO what if we are missing a block id?
  onsets <- if (!is.null(x$onsets)) x$onsets else model_spec$onsets
  durations <- if (!is.null(x$durations)) x$durations else model_spec$durations
  
  varlist <- lapply(seq_along(x$vars), function(i) {
    base::eval(parse(text=x$vars[[i]]), envir=model_spec$event_table, enclos=parent.frame())
  })
  
  names(varlist) <- x$varnames
  
  subs <- if (!is.null(x$subset)) {
    base::eval(x$subset, envir=model_spec$event_table, enclos=parent.frame()) 
  } else {
    rep(TRUE, length(onsets))
  }
  
  et <- event_term(varlist, onsets, model_spec$blockids, durations, subs)
}

#' @export
construct.hrfspec <- function(x, model_spec, ...) {
  ons <- if (!is.null(x$onsets)) x$onsets else model_spec$onsets
  et <- construct_event_term(x,model_spec, ons)
  
  ## could be lazy
  cterm <- convolve(et, x$hrf, model_spec$sampling_frame, summate=x$summate, precision=model_spec$precision)
  
  ret <- list(
    varname=et$varname,
    evterm=et,
    design_matrix=as.data.frame(cterm),
    sampling_frame=model_spec$sampling_frame,
    hrfspec=x,
    contrasts=x$contrasts,
    id=if(!is.null(x$id)) x$id else et$varname
  )
  
  class(ret) <- c("convolved_term", "fmri_term", "list") 
  ret
}


.hrf_parse <- function(..., prefix=NULL, basis=HRF_SPMG1, nbasis=1, lag=0, termsep=":") {
  vars <- as.list(substitute(list(...)))[-1] 
  #browser()
  if (length(vars) > 0) {
    parsed <- parse_term(vars, "hrf")
    term <- parsed$term
    label <- parsed$label
  } else {
    stop("hrf: must have at least one variable")
  }
  
  basis <- if (is.character(basis)) {
    getHRF(basis, nbasis=nbasis, lag=lag)
  } else if (inherits(basis, "HRF")) {
    if (lag > 0) {
      HRF(gen_hrf_lagged(basis, lag=lag), name=basis$name)
    } else {
      basis
    }
    
  } else if (is.function(basis)) {
    test <- basis(1:10)
    HRF(gen_hrf_lagged(basis,lag=lag,...), name="custom_hrf", nbasis=ncol(test))
  } else {
    stop("invalid basis function: must be 1) character string indicating hrf type, e.g. 'gamma' 2) a function or 3) an object of class 'HRF': ", basis)
  }
  
  
  varnames <- if (!is.null(prefix)) {
    paste0(prefix, "_", parsed$term)
  } else {
    parsed$term
    
  }
  
  termname <- paste0(varnames, collapse=termsep)
  
  list(vars=vars, parsed=parsed, term=term, label=label, basis=basis, varnames=varnames, termname=termname)
}


# construct_additive_event_term <- function(x, model_spec) {
#   ## TODO what if we are missing a block id?
#   onsets <- if (!is.null(x$onsets)) x$onsets else model_spec$onsets
#   durations <- if (!is.null(x$durations)) x$durations else model_spec$durations
#   
#   varlist <- lapply(seq_along(x$vars), function(i) {
#     base::eval(parse(text=x$vars[[i]]), envir=model_spec$event_table, enclos=parent.frame())
#   })
#   
#   names(varlist) <- x$varnames
#   
#   subs <- if (!is.null(x$subset)) {
#     base::eval(x$subset, envir=model_spec$event_table, enclos=parent.frame()) 
#   } else {
#     rep(TRUE, length(onsets))
#   }
#   
#   mat <- do.call(cbind, varlist)
#   vlist <- list(mat)
#   names(vlist) <- x$name
#   
#   et <- event_term(vlist, onsets, model_spec$blockids, durations, subs)
# }


# hrf_add <- function(..., basis=HRF_SPMG1, onsets=NULL, durations=NULL,
#                     prefix=NULL, subset=NULL, precision=.2, nbasis=1,contrasts=list(), id=NULL) {
#   parsed <- .hrf_parse(..., prefix=prefix, basis=basis, nbasis=nbasis, termsep="+")
# 
#   if (is.null(id)) {
#     id <- parsed$termname
#   }
# 
# 
# 
#   ret <- list(
#     name=parsed$termname,
#     varnames=parsed$varnames,
#     vars=parsed$term,
#     label=parsed$label,
#     hrf=parsed$basis,
#     onsets=onsets,
#     durations=durations,
#     prefix=prefix,
#     subset=substitute(subset),
#     precision=precision,
#     contrasts=contrasts)
# 
#   class(ret) <- c("hrf_add_spec", "hrfspec", "list")
#   ret
# }
# 
# # @export
# construct.hrf_add_spec <- function(x, model_spec) {
#   et <- construct_additive_event_term(x, model_spec)
#   cterm <- convolve(et, x$hrf, model_spec$sampling_frame, summate=x$summate)
# 
#   ret <- list(
#     varname=et$varname,
#     evterm=et,
#     design_matrix=cterm,
#     sampling_frame=model_spec$sampling_frame,
#     contrasts=x$contrasts,
#     hrfspec=x,
#     id=x$id
#   )
# 
#   class(ret) <- c("convolved_term", "fmri_term", "list")
#   ret
# }


#' trialwise
#' 
#' This function is to be used in formulas for fitting functions, e.g. onsets ~ trialwise() ...
#' 
#' @inheritParams hrf

#' 
#' @param label the label for the generated variable.
#' @param add_sum add the sum of all trialwise regressors to the set. 
#' This can be sued to model the average effect. 
#' 
#' @examples 
#' x <- trialwise(basis="gaussian", onsets=c(1,17,25), durations=c(1,2,3))
#' 
#' @export
trialwise <- function(label="trialwise", basis="spmg1", onsets=NULL, durations=NULL, 
                      prefix=NULL, subset=NULL, precision=.3, id=NULL, add_sum=FALSE) {
 
  termname = label
  
  if (is.null(id)) {
    id <- termname
  }  

  basis <- if (!inherits(basis, "HRF") && is.function(basis)) {
    gen_hrf(basis)
  } else if (is.character(basis)) {
    getHRF(basis)
  } else if (inherits(basis, "HRF")) {
    basis
  } else {
    stop(paste("illegal type for basis arg: ", class(basis)))
  }
  
  assert_that(inherits(basis, "HRF"))
  
  ret <- list(
    name=termname,
    varnames=list("trialwise"),
    label="trialwise",
    hrf=basis,
    onsets=onsets,
    durations=durations,
    prefix=prefix,
    subset=substitute(subset),
    precision=precision,
    contrasts=list(),
    add_sum=add_sum)
  
  class(ret) <- c("trialwisespec", "hrfspec", "list")
  ret
}

#' @export
construct.trialwisespec <- function(x, model_spec, ...) {
  
  ## compied almost verbatim from construct.hrfspec
  onsets <- if (!is.null(x$onsets)) x$onsets else model_spec$onsets
  durations <- if (!is.null(x$durations)) x$durations else model_spec$durations
  
 
  tind <- seq(1, length(onsets))
  trial_index <- formatC(seq(1, length(onsets)), width = nchar(as.character(max(tind))), format = "d", flag = "0")
  #trial_index <- factor(seq(1, length(onsets)))
  trial_index <- factor(trial_index)
  varlist <- list(trial_index)
  names(varlist) <- x$varname
  
  subs <- if (!is.null(x$subset)) base::eval(x$subset, envir=model_spec$event_table, enclos=parent.frame()) else rep(TRUE, length(onsets))
  et <- event_term(varlist, onsets, model_spec$blockids, durations, subs)
  cterm <- convolve(et, x$hrf, model_spec$sampling_frame)
  
  ret <- list(
    varname=et$varname,
    evterm=et,
    design_matrix=cterm,
    sampling_frame=model_spec$sampling_frame,
    contrasts=x$contrasts,
    hrfspec=x,
    id=x$id
  )
  
  class(ret) <- c("trialwise_convolved_term", "convolved_term", "fmri_term", "list") 
  ret
}


#' construct an native AFNI hrf specification for '3dDeconvolve' with the 'stim_times' argument.
#' 
#' @inheritParams hrf
#' @param start the start of the window for sin/poly/csplin models
#' @param stop the stop time for sin/poly/csplin models
#' @export
afni_hrf <- function(..., basis=c("spmg1", "block", "dmblock",           
                                  "tent",   "csplin", "poly",  "sin",        
                                  "gam", "spmg2", "spmg3", "wav"), 
                                  onsets=NULL, durations=NULL, prefix=NULL, subset=NULL, 
                                  nbasis=1, contrasts=NULL, id=NULL, 
                                  start=NULL, stop=NULL) {
  
  ## TODO cryptic error message when argument is mispelled and is then added to ...
  basis <- match.arg(basis)
  
  vars <- as.list(substitute(list(...)))[-1] 
  parsed <- parse_term(vars, "afni_hrf")
  term <- parsed$term
  label <- parsed$label
  
  hrf <- if (!is.null(durations)) {
    assert_that(length(durations) == 1, msg="afni_hrf does not currently accept variable durations")
    get_AFNI_HRF(basis, nbasis=nbasis, duration=durations[1], b=start, c=stop)
  } else {
    get_AFNI_HRF(basis, nbasis=nbasis, b=start, c=stop)
  }
  
  
  varnames <- if (!is.null(prefix)) {
    paste0(prefix, "_", term)
  } else {
    term
  }
  
  termname <- paste0(varnames, collapse="::")
  
  if (is.null(id)) {
    id <- termname
  }  
  
  cset <- if (inherits(contrasts, "contrast_spec")) {
    contrast_set(con1=contrasts)
  } else if (inherits(contrasts, "contrast_set")) {
    contrasts
  } 
  
  ret <- list(
    name=termname,
    id=id,
    varnames=varnames,
    vars=term,
    label=label,
    hrf=hrf,
    onsets=onsets,
    durations=durations,
    prefix=prefix,
    subset=substitute(subset),
    lag=lag,
    contrasts=cset)
  
  class(ret) <- c("afni_hrfspec", "hrfspec", "list")
  ret
  
}

#' construct an native AFNI hrf specification for '3dDeconvolve' and individually modulated events using the 'stim_times_IM' argument.
#' 
#' 
#' @param label name of regressor
#' @param start start of hrf (for multiple basis hrfs)
#' @param stop end of hrf (for multiple basis hrfs)
#' 
#' @inheritParams hrf
#' @examples 
#' 
#' 
#' tw <- afni_trialwise("trialwise", basis="gamma", onsets=seq(1,100,by=5))
#' 
#' @export
afni_trialwise <- function(label, basis=c("spmg1", "block", "dmblock", "gamma", "wav"), 
                     onsets=NULL, durations=0, subset=NULL, 
                      id=NULL, start=0, stop=22) {
  
  ## TODO cryptic error message when argument is mispelled and is then added to ...
  basis <- match.arg(basis)
  
  hrf <- if (!is.null(durations)) {
    assert_that(length(durations) == 1, msg="afni_trialwise does not currently accept variable durations")
    get_AFNI_HRF(basis, nbasis=1, duration=durations[1], b=start, c=stop)
  } else {
    get_AFNI_HRF(basis, nbasis=1, b=start, c=stop)
  }
  
  
  if (is.null(id)) {
    id <- label
  }  
  
  ret <- list(
    name=label,
    varname=label,
    id=id,
    hrf=hrf,
    onsets=onsets,
    durations=durations,
    subset=substitute(subset))
  
  class(ret) <- c("afni_trialwise_hrfspec", "hrfspec", "list")
  ret
  
}

#' @export
construct.afni_hrfspec <- function(x, model_spec, ...) {
  
  et <- construct_event_term(x, model_spec)
  
  ## do not convolve an afni term
  ##cterm <- convolve(et, x$hrf, model_spec$sampling_frame, summate=x$summate)
  
  ret <- list(
    varname=et$varname,
    evterm=et,
    sampling_frame=model_spec$sampling_frame,
    hrfspec=x,
    contrasts=x$contrasts,
    id=if(!is.null(x$id)) x$id else et$varname
  )
  
  class(ret) <- c("afni_hrf_convolved_term", "convolved_term", "fmri_term", "list") 
  ret
}


#' @export
construct.afni_trialwise_hrfspec <- function(x, model_spec, ...) {
  
  ## compied almost verbatim from construct.hrfspec
  onsets <- if (!is.null(x$onsets)) x$onsets else model_spec$onsets
  durations <- if (!is.null(x$durations)) x$durations else model_spec$durations
  
  trial_index <- factor(seq(1, length(onsets)))
  
  varlist <- list(trial_index)
  names(varlist) <- x$varname
  
  subs <- if (!is.null(x$subset)) base::eval(x$subset, envir=model_spec$event_table, enclos=parent.frame()) else rep(TRUE, length(onsets))
  
  et <- event_term(varlist, onsets, model_spec$blockids, durations, subs)
  #cterm <- convolve(et, x$hrf, model_spec$sampling_frame)
  
  ret <- list(
    varname=et$varname,
    evterm=et,
    sampling_frame=model_spec$sampling_frame,
    hrfspec=x,
    id=x$id
  )
  
  class(ret) <- c("afni_trialwise_convolved_term", "convolved_term", "fmri_term", "list") 
  ret
}


#' @keywords internal
AFNI_SPMG1 <- function(d=1) AFNI_HRF(name="SPMG1", nbasis=as.integer(1), params=list(d=d)) 

#' @keywords internal
AFNI_SPMG2 <- function(d=1) AFNI_HRF(name="SPMG2", nbasis=as.integer(2), params=list(d=d))

#' @keywords internal
AFNI_SPMG3 <- function(d=1) AFNI_HRF(name="SPMG3", nbasis=as.integer(3), params=list(d=d))

#' @keywords internal
AFNI_BLOCK <- function(d=1,p=1) AFNI_HRF(name="BLOCK", nbasis=as.integer(1), params=list(d=d,p=p))

#' @keywords internal
AFNI_dmBLOCK <- function(d=1,p=1) AFNI_HRF(name="dmBLOCK", nbasis=as.integer(1), params=list(d=d,p=p))

#' @keywords internal
AFNI_TENT <- function(b=0,c=18, n=10) AFNI_HRF(name="TENT", nbasis=as.integer(n), params=list(b=b,c=c,n=n))

#' @keywords internal
AFNI_CSPLIN <- function(b=0,c=18, n=6) AFNI_HRF(name="CSPLIN", nbasis=as.integer(n), params=list(b=b,c=c,n=n))

#' @keywords internal
AFNI_POLY <- function(b=0,c=18, n=10) AFNI_HRF(name="POLY", nbasis=as.integer(n), params=list(b=b,c=c,n=n))

#' @keywords internal
AFNI_SIN <- function(b=0,c=18, n=10) AFNI_HRF(name="SIN", nbasis=as.integer(n), params=list(b=b,c=c,n=n))

#' @keywords internal
AFNI_GAM <- function(p=8.6,q=.547) AFNI_HRF(name="GAM", nbasis=as.integer(1), params=list(p=p,q=q))

#' @keywords internal
AFNI_WAV <- function(d=1) AFNI_HRF(name="WAV", nbasis=as.integer(1), params=list(d=1))


#' @keywords internal
get_AFNI_HRF <- function(name, nbasis=1, duration=1, b=0, c=18) {
  hrf <- switch(name,
                gamma=AFNI_GAM(),
                spmg1=AFNI_SPMG1(d=duration),
                spmg2=AFNI_SPMG2(d=duration),
                spmg3=AFNI_SPMG3(d=duration),
                csplin=AFNI_CSPLIN(b=b,c=c, n=nbasis),
                poly=AFNI_POLY(b=b,c=c, n=nbasis),
                sine=AFNI_SIN(b=b,c=c,n=nbasis),
                wav=AFNI_WAV(),
                block=AFNI_BLOCK(d=duration),
                dmblock=AFNI_dmBLOCK())
  
  if (is.null(hrf)) {
    stop("could not find afni hrf named: ", name)
  }
  
  hrf
  
}


#' Hemodynamic Response Function with Half-Cosine Basis
#'
#' This function models a hemodynamic response function (HRF) using four half-period cosine basis functions.
#' The HRF consists of an initial dip, a rise to peak, a fall and undershoot, and a recovery to the baseline.
#'
#' @param t A vector of time values.
#' @param h1 Duration of the initial dip in seconds.
#' @param h2 Duration of the rise to peak in seconds.
#' @param h3 Duration of the fall and undershoot in seconds.
#' @param h4 Duration of the recovery to baseline in seconds.
#' @param f1 Height of the starting point.
#' @param f2 Height of the end point.
#' @return A vector of HRF values corresponding to the input time values.
#' @references Woolrich, M. W., Behrens, T. E., & Smith, S. M. (2004). Constrained linear basis sets for HRF modelling using Variational Bayes. NeuroImage, 21(4), 1748-1761.
#' @export
#' @export
hrf_half_cosine <- function(t, h1=1, h2=5, h3=7,h4=7, f1=0, f2=0) {
  rising_half_cosine <- function(t, f1, t0, w) {
    return(f1/2 * (1 - cos(pi * (t - t0) / w)))
  }
  
  falling_half_cosine <- function(t, f1, t0, w) {
    return(f1/2 * (1 + cos(pi * (t - t0) / w)))
  }
  
  ret = dplyr::case_when(
    t < 0 ~ 0,
    t <= h1 ~ falling_half_cosine(t, f1, 0, h1),
    (t > h1) & t <= (h1+h2) ~ rising_half_cosine(t, 1, h1, h2),
    (t > (h1+h2)) & t <= (h1+h2+h3) ~ falling_half_cosine(t,1,(h1+h2), h3),
    (t > (h1+h2+h3)) & t <= (h1+h2+h3+h4) ~ rising_half_cosine(t,f2,(h1+h2+h3), h4),
    (t > h1+h2+h3+h4) ~ f2,
  )
  return(ret)
}


#' HRF Toeplitz Matrix
#' 
#' @description
#' Create a Toeplitz matrix for hemodynamic response function (HRF) convolution.
#' 
#' @param hrf The hemodynamic response function.
#' @param time A numeric vector representing the time points.
#' @param len The length of the output Toeplitz matrix.
#' @param sparse Logical, if TRUE, the output Toeplitz matrix is returned as a sparse matrix (default: FALSE).
#' 
#' @return A Toeplitz matrix for HRF convolution.
hrf_toeplitz <- function(hrf, time, len, sparse=FALSE) {
  hreg <- hrf(time)
  padding <- len - length(hreg)
  H <- pracma::Toeplitz(c(hreg, rep(0, padding)), c(hreg[1], rep(0, len-1)))
  H <- Matrix(H, sparse=sparse)
  H
}

#neural_response_vector <- function(onsets, durations, amplitud)

# Z = as.matrix(H) %*% X
# H = Z_inv X


# construct an hrf that does not convolve it's argument with an response function
# 
# @inheritParams hrf
# @export
# hrf_identity <- function(x, subset=NULL, id=NULL, prefix=NULL) {
#  
#   vars <- substitute(x)
#   
#   term <- as.character(vars)
#   label <- term
#   
#   varnames <- if (!is.null(prefix)) {
#     paste0(prefix, "_", term)
#   } else {
#     term
#   }
#   
#   termname <- paste0(varnames, collapse="::")
#   
#   if (is.null(id)) {
#     id <- termname
#   }  
# 
#   ihrf <- HRF(identity, "ident", nbasis=1)
#   
#   ret <- list(
#     name=termname,
#     id=id,
#     varnames=varnames,
#     vars=term,
#     label=label,
#     hrf=ihrf,
#     prefix=prefix,
#     subset=substitute(subset)
#   )
#   
#   class(ret) <- c("identity_hrfspec", "hrfspec", "list")
#   ret
#   
# }

# @export
# construct.identity_hrfspec <- function(x, model_spec) {
#   
#   subs <- if (!is.null(x$subset)) {
#     base::eval(x$subset, envir=model_spec$event_table, enclos=parent.frame()) 
#   } else {
#     rep(TRUE, length(onsets))
#   }
#   
#   vals <- eval(x$name, envir=model_spec$event_table,enclos=parent.frame() )
#   matrix_term(x$name, vals)
#   
# }


# inv.logit <- plogis
# 
# hrf_logit <- function(t, a1=1, T1=3, T2=6, T3=9, D1=1, D2=-1, D3=1) {
#   a2 <- a1 * (((inv.logit(-T3)/D3) - (inv.logit(-T1)/D1)) / ((inv.logit(-T3)/D3) + (inv.logit(-T2/D2))))
#   print(a2)
#   a3 <- abs(a2) - abs(a1)
#   print(a3)
#   a1*inv.logit((t-T1)/D1) + a2*inv.logit((t-T2)/D2) + a3*inv.logit((t-T3)/D3)
# }
# 
# hrfg <- function(time, amp, mean, sd, c) {
#   amp*dnorm(time, mean=mean, sd=sd) + c
# }
# 
# getPred <- function(parS, xx) hrfg(xx, parS$amp, parS$mean, parS$sd, parS$c)
# residFun <- function(p, observed, xx) observed - getPred(p,xx)
# parStart <- list(amp=1, mean=6, sd=2,c=0)
# df1 <- data.frame(time=0:24, y=hrf_spmg1(0:24))
# nls.out <- nls.lm(par=parStart, fn = residFun, observed = df1$y,
#                   xx = df1$time, 
#                   lower=c(0, 3, .5, -200),
#                   upper=c(100, 12, 4, 200),
#                   control = nls.lm.control(nprint=1))
# 
# 
# 
# ret <- nls(y ~ hrfg(time,amp, lag, mean, sd), control=nls.control(maxiter=5000), data=df1,
#     start=list(amp=1, lag=0, mean=6, sd=2), 
#     lower=c(0,-2, 3, 0), upper=c(2,4,9,3), 
#     algorithm="port", trace=TRUE)

# hrf.logit <- function(t, a1=1, T1, T2, T3, D1, D2, D3) {
#    a2 <- a1 * (((inv.logit(-T3)/D3) - (inv.logit(-T1)/D1)) / ((inv.logit(-T3)/D3) + (inv.logit(-T2/D2))))
#    
#    a3 <- abs(a2) - abs(a1)
#    
#    a1*inv.logit((t-T1)/D1) + a2*inv.logit((t-T2)/D2) + a3*inv.logit((t-T3)/D3)
# }
# 
# hrf.logit2 <- function(t, a1, a2, a3, T1, T2, T3, D1, D2, D3) {
#   a1*inv.logit((t-T1)/D1) + a2*inv.logit((t-T2)/D2) + a3*inv.logit((t-T3)/D3)
# }
# 
# 
# getfun <- function(time) {
#   TIME <- time
#   ret <- function(par) {
#     hrf.logit2(TIME, par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9])
#   }
#   return(ret)
# }
# 
# minimize.fun <- function(yvals, fun, par) {
#   ypred <- fun(par)
#   return(sum((yvals-ypred)^2))
# }
# 
# get.minimizer <- function(yfun, yvals) {
#   YFUN <- yfun
#   YVALS <- yvals
#   ret <- function(par) {
#     ypred <- YFUN(par)
#     return(sum((YVALS-ypred)^2))
#   }
#   
#   return(ret)
# }
# 
# 
# shift.HRF <- function(HRF, shift) {
#   localShift <- shift
#   function(t) {
#     HRF(t+localShift)
#   }
# }
# 
# 
# 
# makeBlock <- function(HRF, duration) {
#   d1 <- duration
#   if (duration < 2) {
#     stop("duration must be greater than 1")
#   }
#   
#   funlist <- c(HRF, lapply(seq(-1, -(duration-1)), function(i) shift.HRF(HRF, i)))
#   function(t) {
#     ret <- numeric(length(t))
#     for (i in 1:length(t)) {
#       ret[i] <- sum(unlist(lapply(funlist, function(fun) fun(t[i]))))
#     }
#     ret
#     
#   }
# }
# 
# 
# makeBlockHRF <- function(eventOnset, duration, HRF) {
#   
#   localHRF <- HRF
#   localOnset <- eventOnset
#   onsets <- seq(eventOnset, eventOnset+duration, 1)
#   funlist <- lapply(onsets, function(onset) makeEventHRF(onset, localHRF))
#   
#   function(t) {
#     ret <- lapply(funlist, function(fun) fun(t)) 
#     Reduce("+", ret)
#   }
# }    
# 
# 
# makeEventHRF <- function(eventOnset, HRF, amp=1) {
#   
#   localHRF <- HRF
#   localOnset <- eventOnset
#   localAmp <- amp
#   function(t) {
#     
#     localHRF(t-localOnset)*amp
#     #for (i in 1:length(t)) {
#     #  ret[i] <- localHRF(t[i]-localOnset)*amp
#     #}
#     
#   }
# }
# 
# .makeEventHRF2 <- function(eventOnset, HRF, amp=1) {
#   
#   localHRF <- HRF
#   localOnset <- eventOnset
#   localAmp <- amp
#   function(t) {
#     ret <- numeric(length(t))
#     for (i in 1:length(t)) {
#       if (t[i] < localOnset) {
#         ret[i] <- 0
#       } else {
#         ret[i] <- localHRF(t[i]-localOnset)*amp
#       }
#     }
#     ret
#   }
# }
# 
# .makeEventHRF3 <- function(eventOnset, HRF, amp=1) {
#   
#   localHRF <- HRF
#   localOnset <- eventOnset
#   localAmp <- amp
#   function(t) {
#     if (t < localOnset) {
#       0
#     } else {
#       localHRF(t-localOnset)*amp
#     }
#     
#   }
# }


    
      
      
