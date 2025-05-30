#' Default value operator
#'
#' @description
#' Returns left side if not NULL, otherwise right side
#'
#' @param lhs Left hand side
#' @param rhs Right hand side (default value)
#'
#' @return lhs if not NULL, otherwise rhs
#' @keywords internal
#' @noRd
`%||%` <- function(lhs, rhs) {
  if (!is.null(lhs)) lhs else rhs
}