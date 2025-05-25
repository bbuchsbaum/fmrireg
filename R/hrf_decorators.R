#' Lag an HRF Object
#'
#' Creates a new HRF object by applying a temporal lag to an existing HRF object.
#'
#' @param hrf The HRF object (of class `HRF`) to lag.
#' @param lag The time lag in seconds to apply. Positive values shift the response later in time.
#'
#' @return A new HRF object representing the lagged function.
#'
#' @family HRF_decorator_functions
#' @export
#' @examples
#' lagged_spmg1 <- lag_hrf(HRF_SPMG1, 5)
#' # Evaluate at time 10; equivalent to HRF_SPMG1(10 - 5)
#' lagged_spmg1(10)
#' HRF_SPMG1(5)
lag_hrf <- function(hrf, lag) {
  assertthat::assert_that(inherits(hrf, "HRF"), msg = "Input 'hrf' must be an HRF object.")
  assertthat::assert_that(is.numeric(lag) && length(lag) == 1, msg = "'lag' must be a single numeric value.")

  # Original attributes
  orig_name <- attr(hrf, "name")
  orig_span <- attr(hrf, "span")
  orig_nbasis <- nbasis(hrf)
  orig_params <- attr(hrf, "params")

  # Create the lagged function
  lagged_func <- function(t) {
    hrf(t - lag)
  }

  # Create new HRF object using as_hrf
  as_hrf(
    f = lagged_func,
    name = paste0(orig_name, "_lag(", lag, ")"),
    nbasis = orig_nbasis,
    span = orig_span + max(0, lag), # Increase span if lag is positive
    params = c(orig_params, list(.lag = lag)) # Add lag to params for bookkeeping
  )
}


#' Create a Blocked HRF Object
#'
#' Creates a new HRF object representing a response to a sustained (blocked)
#' stimulus by convolving the input HRF with a boxcar function of a given width.
#'
#' @param hrf The HRF object (of class `HRF`) to block.
#' @param width The width of the block in seconds.
#' @param precision The sampling precision in seconds used for the internal convolution (default: 0.1).
#' @param half_life The half-life of an optional exponential decay applied during the block (default: Inf, meaning no decay).
#' @param summate Logical; if TRUE (default), the responses from each time point within the block are summed. If FALSE, the maximum response at each time point is taken.
#' @param normalize Logical; if TRUE, the resulting blocked HRF is scaled so that its peak value is 1 (default: FALSE).
#'
#' @return A new HRF object representing the blocked function.
#'
#' @family HRF_decorator_functions
#' @export
#' @examples
#' blocked_spmg1 <- block_hrf(HRF_SPMG1, width = 5)
#' t_vals <- seq(0, 30, by = 0.5)
#' plot(t_vals, HRF_SPMG1(t_vals), type = 'l', col = "blue", ylab = "Response", xlab = "Time")
#' lines(t_vals, blocked_spmg1(t_vals), col = "red")
#' legend("topright", legend = c("Original", "Blocked (width=5)"), col = c("blue", "red"), lty = 1)
block_hrf <- function(hrf, width, precision = 0.1, half_life = Inf, summate = TRUE, normalize = FALSE) {
  assertthat::assert_that(inherits(hrf, "HRF"), msg = "Input 'hrf' must be an HRF object.")
  assertthat::assert_that(is.numeric(width) && length(width) == 1 && width >= 0, msg = "'width' must be a single non-negative numeric value.")
  assertthat::assert_that(is.numeric(precision) && length(precision) == 1 && precision > 0, msg = "'precision' must be a single positive numeric value.")
  assertthat::assert_that(is.numeric(half_life) && length(half_life) == 1 && half_life > 0, msg = "'half_life' must be a single positive numeric value (use Inf for no decay).")
  assertthat::assert_that(is.logical(summate) && length(summate) == 1, msg = "'summate' must be a single logical value.")
  assertthat::assert_that(is.logical(normalize) && length(normalize) == 1, msg = "'normalize' must be a single logical value.")

  # Original attributes
  orig_name <- attr(hrf, "name")
  orig_span <- attr(hrf, "span")
  orig_nbasis <- nbasis(hrf)
  orig_params <- attr(hrf, "params")

  # Create the blocked function
  blocked_func <- function(t) {
    if (width < precision) {
      # If width is negligible, just return the original hrf value
      res <- hrf(t)
    } else {
      samples <- seq(0, width, by = precision)
      hmat_list <- lapply(samples, function(offset) {
        decay_factor <- if (is.infinite(half_life)) 1 else exp(-log(2) * offset / half_life)
        hrf(t - offset) * decay_factor
      })
      
      # Combine results for each time point t across offsets
      if (orig_nbasis == 1) {
        hmat <- do.call(cbind, hmat_list) # Matrix with rows=time, cols=offsets
        res <- if (summate) {
          rowSums(hmat)
        } else {
          apply(hmat, 1, function(vals) vals[which.max(vals)])
        }
      } else {
         # For multi-basis, sum matrices if summate=TRUE (difficult to define max sensibly)
         if (!summate) {
           warning("'summate = FALSE' is not fully supported for multi-basis HRFs in block_hrf; returning sum.")
         }
         res <- Reduce("+", hmat_list)
      }
    }
    
    # Apply normalization if requested
    if (normalize) {
      if (orig_nbasis == 1) {
        peak_val <- max(abs(res), na.rm = TRUE)
        if (!is.na(peak_val) && peak_val != 0) {
          res <- res / peak_val
        }
      } else {
        # Normalize each basis column independently
        res <- apply(res, 2, function(basis_col) {
           peak_val <- max(abs(basis_col), na.rm = TRUE)
           if (!is.na(peak_val) && peak_val != 0) {
             basis_col / peak_val
           } else {
             basis_col
           }
        })
      }
    }
    return(res)
  }

  # Store parameters used for blocking
  block_params <- list(
      .width = width,
      .precision = precision,
      .half_life = half_life,
      .summate = summate,
      .normalize = normalize
  )
  
  # Create new HRF object using as_hrf
  as_hrf(
    f = blocked_func,
    name = paste0(orig_name, "_block(w=", width, ")"),
    nbasis = orig_nbasis,
    span = orig_span + width, # Span increases by the block width
    params = c(orig_params, block_params) # Add block params for bookkeeping
  )
}


#' Normalise an HRF Object
#'
#' Creates a new HRF object whose output is scaled such that the maximum absolute
#' value of the response is 1.
#'
#' @param hrf The HRF object (of class `HRF`) to normalise.
#'
#' @return A new HRF object representing the normalised function.
#' @details For multi-basis HRFs, each basis function (column) is normalised independently.
#'
#' @family HRF_decorator_functions
#' @export
#' @examples
#' # Create a gaussian HRF with a peak value != 1
#' gauss_unnorm <- as_hrf(function(t) 5 * dnorm(t, 6, 2), name="unnorm_gauss")
#' # Normalise it
#' gauss_norm <- normalise_hrf(gauss_unnorm)
#' t_vals <- seq(0, 20, by = 0.1)
#' max(gauss_unnorm(t_vals)) # Peak is > 1
#' max(gauss_norm(t_vals))   # Peak is 1
normalise_hrf <- function(hrf) {
  assertthat::assert_that(inherits(hrf, "HRF"), msg = "Input 'hrf' must be an HRF object.")

  # Original attributes
  orig_name <- attr(hrf, "name")
  orig_span <- attr(hrf, "span")
  orig_nbasis <- nbasis(hrf)
  orig_params <- attr(hrf, "params")

  # Create the normalised function
  normalised_func <- function(t) {
    res <- hrf(t)
    if (orig_nbasis == 1) {
      peak_val <- max(abs(res), na.rm = TRUE)
      if (!is.na(peak_val) && peak_val != 0) {
        res <- res / peak_val
      }
    } else if (is.matrix(res)) {
      # Normalise each basis column independently
      res <- apply(res, 2, function(basis_col) {
        peak_val <- max(abs(basis_col), na.rm = TRUE)
        if (!is.na(peak_val) && peak_val != 0) {
          basis_col / peak_val
        } else {
          basis_col
        }
      })
    }
    # If it's not numeric or matrix (e.g., NULL or error result), return as is
    return(res)
  }

  # Create new HRF object using as_hrf
  as_hrf(
    f = normalised_func,
    name = paste0(orig_name, "_norm"),
    nbasis = orig_nbasis,
    span = orig_span,
    params = c(orig_params, list(.normalised = TRUE)) # Add flag for bookkeeping
  )
} 