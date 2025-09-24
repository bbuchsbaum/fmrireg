#' Legacy AR Function Compatibility Layer
#'
#' These functions provide backward compatibility for code that still
#' references the old AR implementation. They now delegate to fmriAR.
#'
#' @keywords internal
#' @noRd
NULL

#' Estimate AR parameters (legacy compatibility)
#'
#' @param residuals_vec Numeric vector of residuals
#' @param p_order Integer AR order
#' @return Numeric vector of AR coefficients
#' @keywords internal
#' @noRd
estimate_ar_parameters <- function(residuals_vec, p_order) {
  stopifnot(is.numeric(residuals_vec))
  stopifnot(length(p_order) == 1L, p_order >= 1L)

  if (length(residuals_vec) == 0L) {
    stop("Residual vector must contain at least one value", call. = FALSE)
  }

  if (anyNA(residuals_vec)) {
    stop("NA values detected in 'residuals_vec' for estimate_ar_parameters")
  }

  # Try fmriAR first
  residuals_mat <- matrix(residuals_vec, ncol = 1)
  plan <- tryCatch(
    fmriAR::fit_noise(
      resid = residuals_mat,
      method = "ar",
      p = p_order,
      p_max = p_order
    ),
    error = function(e) NULL
  )

  if (!is.null(plan) && !is.null(plan$phi) && length(plan$phi) > 0L) {
    phi <- plan$phi[[1]]
    if (length(phi) > 0L) {
      return(as.numeric(phi))
    }
  }

  # Fall back to Yule-Walker implementation
  fallback_phi <- .estimate_ar(residuals_vec, p_order)
  if (length(fallback_phi) == 0L) {
    return(numeric(0))
  }

  as.numeric(fallback_phi)
}

#' Apply AR whitening transform (legacy compatibility)
#'
#' @param X Design matrix
#' @param Y Data matrix
#' @param phi AR coefficients
#' @param exact_first Logical for exact first sample scaling
#' @return List with whitened X and Y
#' @keywords internal
#' @noRd
ar_whiten_transform <- function(X, Y, phi, exact_first = FALSE, ...) {
  if (anyNA(X) || anyNA(Y)) {
    stop("NA values detected in 'X' or 'Y' for ar_whiten_transform")
  }

  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)

  # Debug
  if (getOption("fmrireg.debug.ar", FALSE)) {
    message("ar_whiten_transform (via fmriAR) called with phi=", paste(phi, collapse=","),
            " exact_first=", exact_first)
    message("  X dims: ", nrow(X), "x", ncol(X))
    message("  Y dims: ", nrow(Y), "x", ncol(Y))
  }

  # Normalize phi/theta for fmriAR plan creation
  phi_list <- if (is.list(phi)) phi else list(phi)
  theta_list <- replicate(length(phi_list), numeric(0), simplify = FALSE)
  pooling_mode <- if (length(phi_list) > 1L) "run" else "global"
  phi_arg <- if (identical(pooling_mode, "global")) phi_list[[1]] else phi_list

  # Use fmriAR (ensure theta provided to avoid empty list issues)
  plan <- fmriAR::compat$plan_from_phi(
    phi = phi_arg,
    theta = theta_list,
    pooling = pooling_mode,
    exact_first = exact_first
  )

  if (is.null(plan$theta) || length(plan$theta) == 0L) {
    plan$theta <- theta_list
  }

  result <- fmriAR::whiten_apply(
    plan = plan,
    X = X,
    Y = Y,
    parallel = TRUE
  )

  list(X = result$X, Y = result$Y)
}

#' Estimate AR via Yule-Walker (legacy compatibility)
#'
#' @param residuals_vec Numeric vector
#' @param p_order AR order
#' @return Numeric vector of AR coefficients
#' @keywords internal
#' @noRd
.estimate_ar <- function(residuals_vec, p_order) {
  stopifnot(is.numeric(residuals_vec))
  stopifnot(length(p_order) == 1L, p_order >= 1L)

  # Center residuals
  res <- residuals_vec - mean(residuals_vec, na.rm = TRUE)

  if (stats::sd(res, na.rm = TRUE) < sqrt(.Machine$double.eps)) {
    return(numeric(0))
  }

  # Use stats::ar for Yule-Walker estimation (matching old behavior)
  ar_fit <- stats::ar(res, order.max = p_order, method = "yule-walker",
                       demean = FALSE, aic = FALSE)

  # Return coefficients
  if (length(ar_fit$ar) > 0) {
    as.numeric(ar_fit$ar)
  } else {
    numeric(0)
  }
}

#' AR whitening in-place (legacy compatibility)
#'
#' Compatibility wrapper for code expecting ar_whiten_inplace function
#'
#' @param Y Data matrix
#' @param X Design matrix
#' @param phi AR coefficients
#' @param exact_first_ar1 Logical for exact first scaling
#' @return List with whitened Y and X
#' @keywords internal
#' @noRd
ar_whiten_inplace <- function(Y, X, phi = NULL, exact_first_ar1 = FALSE, ...) {
  dots <- list(...)
  if (is.null(phi) && !is.null(dots$phi_coeffs)) {
    phi <- dots$phi_coeffs
    dots$phi_coeffs <- NULL
  }
  if (is.null(phi)) {
    phi <- numeric(0)
  }

  # Note: fmriAR doesn't modify in place, but we maintain the interface
  phi_list <- if (is.list(phi)) phi else list(phi)
  theta_list <- replicate(length(phi_list), numeric(0), simplify = FALSE)
  pooling_mode <- if (length(phi_list) > 1L) "run" else "global"
  phi_arg <- if (identical(pooling_mode, "global")) phi_list[[1]] else phi_list

  plan <- fmriAR::compat$plan_from_phi(
    phi = phi_arg,
    theta = theta_list,
    pooling = pooling_mode,
    exact_first = exact_first_ar1
  )

  if (is.null(plan$theta) || length(plan$theta) == 0L) {
    plan$theta <- theta_list
  }

  result <- fmriAR::whiten_apply(
    plan = plan,
    X = X,
    Y = Y,
    parallel = TRUE
  )

  # Return in expected order (Y, X)
  list(Y = result$Y, X = result$X)
}
