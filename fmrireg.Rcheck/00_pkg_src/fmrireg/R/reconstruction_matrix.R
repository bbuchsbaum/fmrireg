#' Reconstruction matrix for an HRF basis
#'
#' Returns a matrix \eqn{\Phi} that converts basis coefficients into a
#' sampled HRF shape.
#'
#' @param hrf An object of class `HRF`.
#' @param sframe A `sampling_frame` object or numeric vector of times.
#' @param precision Optional sampling interval in seconds when `sframe`
#'   is a `sampling_frame`. Defaults to the TR of `sframe`.
#' @return A numeric matrix with one column per basis function.
#' @export
reconstruction_matrix <- function(hrf, sframe, ...) {
  UseMethod("reconstruction_matrix")
}

#' @export
#' @rdname reconstruction_matrix
reconstruction_matrix.HRF <- function(hrf, sframe, precision = NULL, ...) {
  grid <- if (inherits(sframe, "sampling_frame")) {
    if (is.null(precision)) precision <- sframe$TR[1]
    seq(0, attr(hrf, "span"), by = precision)
  } else {
    as.numeric(sframe)
  }
  vals <- evaluate(hrf, grid)
  if (is.vector(vals)) matrix(vals, ncol = 1L) else as.matrix(vals)
}
