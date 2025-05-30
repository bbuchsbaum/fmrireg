#' Extract Estimated AR Parameters from fmri_lm Fit
#'
#' Retrieves the estimated autoregressive parameters from a fitted
#' fMRI linear model that used AR error modeling.
#'
#' @param object An object of class \code{fmri_lm}
#' @param ... Additional arguments (currently unused)
#'
#' @return A numeric vector of AR coefficients if AR modeling was used,
#'   or NULL if no AR modeling was performed (cor_struct = "iid").
#'   For runwise estimation without global AR, returns the average of
#'   run-specific estimates.
#'
#' @export
#' @examples
#' \dontrun{
#' # Fit model with AR(1) errors
#' fit <- fmri_lm(onset ~ hrf(cond), dataset = dset, cor_struct = "ar1")
#' ar_parameters(fit)  # Extract estimated AR(1) coefficient
#' }
ar_parameters <- function(object, ...) {
  UseMethod("ar_parameters")
}

#' @rdname ar_parameters
#' @export
ar_parameters.fmri_lm <- function(object, ...) {
  # Check if AR was used
  cfg <- attr(object, "config")
  if (is.null(cfg) || cfg$ar$struct == "iid") {
    return(NULL)
  }
  
  # TODO: This requires storing AR parameters during fitting
  # For now, return NULL with a message
  message("AR parameter extraction not yet implemented. This feature will be added in a future version.")
  return(NULL)
}