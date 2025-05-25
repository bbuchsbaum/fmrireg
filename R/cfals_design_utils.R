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
#' @param X_list A list of design matrices.
#' @param Y Data matrix with matching rows.
#' @param confounds Optional confound matrix with matching rows.
#' @return A list with projected `X_list` and `Y` matrices.
#' @export
project_confounds <- function(X_list, Y, confounds = NULL) {
  if (is.null(confounds)) return(list(X_list = X_list, Y = Y))
  qrZ <- qr(confounds)
  Xp <- lapply(X_list, function(X) qr.resid(qrZ, X))
  Yp <- qr.resid(qrZ, Y)
  list(X_list = Xp, Y = Yp)
}
