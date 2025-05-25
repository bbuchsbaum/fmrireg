
#' @importFrom splines bs
#' @importFrom stats dgamma dnorm quantile
NULL

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
hrf_ident <- function(t) {
  ifelse( t == 0, 1, 0)
}

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
#' @param ... Additional arguments passed to `splines::bs`.
hrf_bspline <- function(t, span=24, N=5, degree=3, ...) {
	
	ord <- 1 + degree
	# Check if requested N is sufficient for the degree
	if (N < ord) {
	    warning(paste0("Requested N=", N, " basis functions is less than degree+1=", ord, ". ",
	                   "Using minimum required of ", ord, " basis functions."))
	    # We don't change N here, let splines::bs handle the df inconsistency if needed,
	    # but the warning informs the user.
	}
	
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
	
	splines::bs(t, df=N, knots=knots, degree=degree, Boundary.knots=c(0,span),...)
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
#' hrf_gamma_vals <- hrf_gamma(seq(0, 20, by = .5), shape = 6, rate = 1)
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
#' hrf_gaussian_vals <- hrf_gaussian(seq(0, 20, by = .5), mean = 6, sd = 2)
#' @export
hrf_gaussian <- function(t, mean=6, sd=2) {
	stats::dnorm(t, mean=mean, sd=sd)
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
#' hrf_mexhat_vals <- hrf_mexhat(seq(0, 20, by = .5), mean = 6, sd = 2)
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


# Fast analytic first derivative for hrf_spmg1
#' @keywords internal
#' @noRd
hrf_spmg1_deriv <- function(t, P1 = 5, P2 = 15, A1 = .0833) {
  C <- 1.274527e-13
  ret <- numeric(length(t))
  pos <- t >= 0
  if (any(pos)) {
    t_pos <- t[pos]
    ret[pos] <- exp(-t_pos) * (A1 * t_pos^(P1 - 1) * (P1 - t_pos) -
                                 C   * t_pos^(P2 - 1) * (P2 - t_pos))
  }
  ret
}

# Fast analytic second derivative for hrf_spmg1
#' @keywords internal
#' @noRd
hrf_spmg1_second_deriv <- function(t, P1 = 5, P2 = 15, A1 = .0833) {
  C <- 1.274527e-13
  ret <- numeric(length(t))
  pos <- t >= 0
  if (any(pos)) {
    t_pos <- t[pos]
    # Let D1 = A1 * t^(P1-1) * (P1 - t) and D2 = C * t^(P2-1) * (P2 - t)
    D1 <- A1 * t_pos^(P1 - 1) * (P1 - t_pos)
    D2 <- C   * t_pos^(P2 - 1) * (P2 - t_pos)
    # Their derivatives:
    D1_prime <- A1 * ((P1 - 1) * t_pos^(P1 - 2) * (P1 - t_pos) - t_pos^(P1 - 1))
    D2_prime <- C   * ((P2 - 1) * t_pos^(P2 - 2) * (P2 - t_pos) - t_pos^(P2 - 1))
    ret[pos] <- exp(-t_pos) * (D1_prime - D2_prime - (D1 - D2))
  }
  ret
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
#' hrf_inv_logit_basis <- hrf_inv_logit(seq(0, 20, by = 0.5), mu1 = 6, s1 = 1, mu2 = 16, s2 = 1)
hrf_inv_logit <- function(t, mu1 = 6, s1 = 1, mu2 = 16, s2 = 1, lag = 0) {
  inv_logit1 <- 1 / (1 + exp(-(t - lag - mu1) / s1))
  inv_logit2 <- 1 / (1 + exp(-(t - lag - mu2) / s2))
  return(inv_logit1 - inv_logit2)
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

#' Fourier basis for HRF modeling
#'
#' Generates a set of Fourier basis functions (sine and cosine pairs) over a given span.
#'
#' @param t A vector of time points.
#' @param span The temporal window over which the basis functions span (default: 24).
#' @param nbasis The number of basis functions (default: 5). Should be even for full sine-cosine pairs.
#' @return A matrix of Fourier basis functions with nbasis columns.
#' @export
hrf_fourier <- function(t, span = 24, nbasis = 5) {
  freqs <- ceiling(seq_len(nbasis) / 2)
  basis <- sapply(seq_len(nbasis), function(k) {
    n <- freqs[k]
    if (k %% 2 == 1) {
      sin(2 * pi * n * t / span)
    } else {
      cos(2 * pi * n * t / span)
    }
  })
  return(basis)
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
#' @export
hrf_toeplitz <- function(hrf, time, len, sparse=FALSE) {
  hreg <- hrf(time)
  padding <- len - length(hreg)
  H <- pracma::Toeplitz(c(hreg, rep(0, padding)), c(hreg[1], rep(0, len-1)))
  H <- Matrix::Matrix(H, sparse=sparse)
  H
}


#' Generate Daguerre spherical basis functions
#' 
#' @description
#' Creates a set of Daguerre spherical basis functions. These are orthogonal 
#' polynomials on [0,∞) with respect to the weight function w(x) = x^2 * exp(-x).
#' They are particularly useful for modeling hemodynamic responses as they naturally
#' decay to zero and can capture various response shapes.
#'
#' @param t Time points at which to evaluate the basis functions
#' @param n_basis Number of basis functions to generate (default: 3)
#' @param scale Scale parameter for the time axis (default: 1)
#' @return A matrix with columns containing the basis functions
#' @keywords internal
#' @noRd
daguerre_basis <- function(t, n_basis = 3, scale = 1) {
  # Scale time
  x <- t/scale
  
  # Initialize basis matrix
  basis <- matrix(0, length(x), n_basis)
  
  # First basis function (n=0)
  basis[,1] <- exp(-x/2)
  
  if(n_basis > 1) {
    # Second basis function (n=1)
    basis[,2] <- (1 - x) * exp(-x/2)
  }
  
  if(n_basis > 2) {
    # Higher order basis functions using recurrence relation
    for(n in 3:n_basis) {
      k <- n - 1
      basis[,n] <- ((2*k - 1 - x) * basis[,n-1] - (k - 1) * basis[,n-2]) / k
    }
  }
  
  # Normalize basis functions
  for(i in 1:n_basis) {
    # Avoid division by zero if a basis function is all zero
    max_abs_val <- max(abs(basis[,i]))
    if (max_abs_val > 1e-10) {
      basis[,i] <- basis[,i] / max_abs_val
    }
  }
  
  basis
}

