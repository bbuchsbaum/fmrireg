#' CFALS Design Utilities
#'
#' Helper functions for interfacing the CF-ALS engine with the
#' fmrireg HRF basis system.
#'
#' @name cfals_design_utils
NULL

#' Reconstruction matrix for an HRF basis
#'
#' Returns a matrix \eqn{\Phi} that converts basis coefficients into a
#' sampled HRF shape.
#'
#' @param hrf An object of class `HRF`.
#' @param sframe A `sampling_frame` object or numeric vector of times.
#' @return A numeric matrix with one column per basis function.
#' @export
reconstruction_matrix <- function(hrf, sframe) {
  UseMethod("reconstruction_matrix")
}

#' @export
reconstruction_matrix.HRF <- function(hrf, sframe) {
  grid <- if (inherits(sframe, "sampling_frame")) {
    seq(0, attr(hrf, "span"), by = sframe$TR[1])
  } else {
    as.numeric(sframe)
  }
  vals <- evaluate(hrf, grid)
  if (is.vector(vals)) matrix(vals, ncol = 1L) else as.matrix(vals)
}

#' Penalty matrix for an HRF basis
#'
#' Provides a ridge or smoothness penalty matrix for the basis
#' coefficients. The default method returns an identity matrix.
#'
#' @param hrf An object of class `HRF`.
#' @return A square numeric matrix.
#' @export
penalty_matrix <- function(hrf) {
  UseMethod("penalty_matrix")
}

#' @export
penalty_matrix.HRF <- function(hrf) {
  diag(nbasis(hrf))
}

#' Project design and data matrices to the null space of confounds
#'
#' Projects both the data matrix `Y` and each design matrix in
#' `X_list` using QR decomposition of the confound matrix.  The
#' projection can optionally use LAPACK's QR implementation for
#' improved numerical stability.
#'
#' @param Y Numeric matrix of BOLD data (time points \eqn{\times}
#'   voxels).
#' @param X_list A list of design matrices with the same number of
#'   rows as `Y`.
#' @param confounds Optional confound matrix with matching rows.
#' @param lapack_qr Logical; passed to `qr()` as the `LAPACK`
#'   argument.
#' @return A list with projected `X_list` and `Y` matrices.
#' @export
project_confounds <- function(Y, X_list, confounds = NULL, lapack_qr = TRUE) {
  if (is.null(confounds)) {
    return(list(X_list = X_list, Y = Y))
  }
  qrZ <- qr(confounds, LAPACK = lapack_qr)
  Xp <- lapply(X_list, function(X) qr.resid(qrZ, X))
  Yp <- qr.resid(qrZ, Y)
  list(X_list = Xp, Y = Yp)
}

#' Create design matrices for CFALS estimation
#'
#' Convenience helper that constructs the list of design matrices for
#' a given `event_model` and HRF basis.  It also returns useful
#' metadata such as the number of basis functions and conditions as
#' well as a reconstruction matrix for converting HRF coefficients to
#' sampled shapes and a normalised reference HRF vector for sign
#' alignment.
#'
#' @param event_model An object of class `event_model`.
#' @param hrf_basis An `HRF` basis object.
#' @return A list with elements `X_list`, `d`, `k`, `Phi`, and
#'   `h_ref_shape_norm`.
#' @export
create_fmri_design <- function(event_model, hrf_basis) {
  if (!inherits(event_model, "event_model")) {
    stop("'event_model' must be an 'event_model' object")
  }
  if (!inherits(hrf_basis, "HRF")) {
    stop("'hrf_basis' must be an object of class 'HRF'")
  }

  sframe <- event_model$sampling_frame
  sample_times <- samples(sframe, global = TRUE)

  reg_lists <- lapply(event_model$terms, regressors.event_term,
                      hrf = hrf_basis,
                      sampling_frame = sframe,
                      summate = FALSE,
                      drop.empty = TRUE)
  regs <- unlist(reg_lists, recursive = FALSE)
  cond_names <- names(regs)
  X_list <- lapply(regs, function(r)
    evaluate(r, sample_times, precision = sframe$precision))
  names(X_list) <- cond_names

  Phi <- reconstruction_matrix(hrf_basis, sframe)
  h_ref <- drop(Phi[, 1])
  h_ref <- h_ref / max(abs(h_ref))

  list(X_list = X_list,
       d = nbasis(hrf_basis),
       k = length(X_list),
       Phi = Phi,
       h_ref_shape_norm = h_ref)
}
