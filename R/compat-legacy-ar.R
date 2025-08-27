#' Legacy Compatibility Functions for AR Modeling
#'
#' @description
#' Provides backward compatibility for functions that have been renamed
#' or refactored. These aliases ensure existing code continues to work.
#'
#' @keywords internal
#' @noRd
NULL

#' Alias for internal preproject function
#'
#' @description
#' Provides a consistent interface for the fast preprojection function.
#' Both .fast_preproject and fast_preproject point to the same implementation.
#'
#' @param ... All arguments passed to the internal function
#' @return Same as .fast_preproject
#' @keywords internal
#' @noRd
fast_preproject <- function(X) {
  # Call the actual implementation in fmri_lm_internal.R
  .fast_preproject(X)
}