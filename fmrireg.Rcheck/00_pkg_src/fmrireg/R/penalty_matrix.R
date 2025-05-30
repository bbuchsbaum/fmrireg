#' Penalty Matrix Methods for HRF Objects
#'
#' @description
#' This file contains S3 methods for generating penalty matrices for different
#' types of HRF basis functions. Penalty matrices encode shape priors that
#' improve regularized HRF estimation by discouraging implausible shapes.
#'
#' @name penalty_matrix_methods
NULL

# Helper functions for difference matrices ----

#' Create first-order difference matrix
#' @param d Dimension (number of basis functions)
#' @return First-order difference matrix
#' @keywords internal
#' @noRd
first_difference <- function(d) {
  if (d <= 1) {
    return(matrix(0, nrow = 0, ncol = d))
  }
  diff(diag(d), differences = 1)
}

#' Create second-order difference matrix  
#' @param d Dimension (number of basis functions)
#' @return Second-order difference matrix
#' @keywords internal
#' @noRd
second_difference <- function(d) {
  if (d <= 2) {
    return(matrix(0, nrow = 0, ncol = d))
  }
  diff(diag(d), differences = 2)
}

# Default method ----

#' Default penalty matrix method
#' 
#' @description
#' Returns an identity matrix (ridge penalty) when no specific penalty structure
#' is known for the HRF type. This method also handles name-based dispatch for
#' HRF objects since they all share the same class structure.
#' 
#' @param x An HRF object
#' @param ... Additional arguments passed to specific penalty methods
#' @return Identity matrix of dimension nbasis(x) × nbasis(x), or appropriate penalty matrix based on HRF name
#' @export
penalty_matrix.HRF <- function(x, ...) {
  # First try name-based dispatch
  hrf_name <- attr(x, "name")
  
  if (!is.null(hrf_name)) {
    # Normalize name for comparison
    name_lower <- tolower(hrf_name)
    
    if (name_lower == "fir") {
      return(penalty_matrix_fir(x, ...))
    } else if (name_lower %in% c("bspline", "bs")) {
      return(penalty_matrix_bspline(x, ...))
    } else if (grepl("spmg", name_lower)) {
      return(penalty_matrix_spmg(x, ...))
    } else if (name_lower %in% c("cosine", "fourier")) {
      return(penalty_matrix_fourier(x, ...))
    }
  }
  
  # Default: identity matrix (ridge penalty)
  d <- nbasis(x)
  diag(d)
}

# Specific penalty matrix implementations ----

#' FIR penalty matrix
#' @param hrf FIR HRF object
#' @param order Derivative order for roughness (1 or 2)
#' @param scale Logical; if TRUE, scales R so that diag(R) ≈ 1
#' @return Penalty matrix for FIR basis
#' @keywords internal
#' @noRd
penalty_matrix_fir <- function(hrf, order = 2, scale = TRUE) {
  d <- nbasis(hrf)
  
  if (d <= order) {
    # Not enough basis functions for the requested order
    return(diag(d))
  }
  
  D <- switch(as.character(order),
              "1" = first_difference(d),
              "2" = second_difference(d),
              stop("order must be 1 or 2"))
  
  R <- crossprod(D)  # D^T D
  
  if (scale && any(diag(R) > 0)) {
    R <- R / mean(diag(R)[diag(R) > 0])  # Scale so diagonal ≈ 1
  }
  
  R
}

#' B-spline penalty matrix
#' @param hrf B-spline HRF object
#' @param order Derivative order for roughness (1 or 2)
#' @param scale Logical; if TRUE, scales R so that diag(R) ≈ 1
#' @return Penalty matrix for B-spline basis
#' @keywords internal
#' @noRd
penalty_matrix_bspline <- function(hrf, order = 2, scale = TRUE) {
  # B-splines use the same difference-based penalty as FIR
  penalty_matrix_fir(hrf, order = order, scale = scale)
}

#' SPM canonical + derivatives penalty matrix
#' @param hrf SPM HRF object (canonical + derivatives)
#' @param shrink_deriv Shrinkage factor for derivative columns
#' @param shrink_disp Shrinkage factor for dispersion derivative (if present)
#' @param ... Additional arguments
#' @return Penalty matrix for SPM basis
#' @keywords internal
#' @noRd
penalty_matrix_spmg <- function(hrf, shrink_deriv = 4, shrink_disp = NULL, ...) {
  d <- nbasis(hrf)
  
  if (d == 1) {
    # Single canonical HRF - no penalty needed
    return(diag(1))
  } else if (d == 2) {
    # Canonical + temporal derivative
    return(diag(c(1, shrink_deriv)))
  } else if (d == 3) {
    # Canonical + temporal + dispersion derivatives
    if (is.null(shrink_disp)) {
      shrink_disp <- shrink_deriv
    }
    return(diag(c(1, shrink_deriv, shrink_disp)))
  } else {
    # Unexpected number of basis functions - fall back to identity
    return(diag(d))
  }
}

#' Fourier/Cosine penalty matrix
#' @param hrf Fourier/Cosine HRF object
#' @param order Derivative order for roughness (1 or 2)
#' @param scale Logical; if TRUE, scales R so that diag(R) ≈ 1
#' @return Penalty matrix for Fourier basis
#' @keywords internal
#' @noRd
penalty_matrix_fourier <- function(hrf, order = 2, scale = TRUE) {
  # Fourier bases use the same difference-based penalty as FIR
  # This penalizes high-frequency coefficients
  penalty_matrix_fir(hrf, order = order, scale = scale)
}

# Additional convenience functions ----

#' Create penalty matrix for any HRF-like object
#' 
#' @description
#' Convenience wrapper that can handle different types of HRF specifications
#' including hrfspec objects from model formulas.
#' 
#' @param x An HRF object, hrfspec, or other HRF-like object
#' @param ... Additional arguments passed to penalty_matrix methods
#' @return Appropriate penalty matrix
#' @export
penalty_matrix.hrfspec <- function(x, ...) {
  penalty_matrix(x$hrf, ...)
} 