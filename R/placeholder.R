
#' Deprecated helper: estimate()
#'
#' This function is deprecated. Use \code{estimate_betas()} instead.
#'
#' @param ... Ignored.
#' @return No return value; always errors with a deprecation message.
#' @export
estimate <- function(...) {
  stop("estimate() is deprecated. Please use estimate_betas() instead.")
}
