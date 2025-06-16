#' Design Matrix Method for fmri_lm Objects
#'
#' Extract the design matrix from an fmri_lm object by delegating to its model component.
#'
#' @param x An fmri_lm object
#' @param ... Additional arguments passed to the design_matrix method for the model
#' @return The design matrix from the fmri_lm object's model
#' @export
#' @method design_matrix fmri_lm
design_matrix.fmri_lm <- function(x, ...) {
  design_matrix(x$model, ...)
}