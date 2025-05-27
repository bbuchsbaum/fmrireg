#' GLM Context Object
#'
#' Internal S3 class for storing matrices and intermediate values used in
#' linear model fitting. The context bundles the design matrix, data matrix
#' and related metadata so helper functions can operate with a single object.
#'
#' @param X Design matrix (time points \eqn{\times} predictors).
#' @param Y Data matrix (time points \eqn{\times} observations).
#' @param proj Preprojection object created by \code{.fast_preproject(X)}.
#' @param phi_hat Optional autoregressive parameter estimates.
#' @param sigma_robust_scale Optional robust scale estimate.
#' @param robust_weights Optional vector of robust weights.
#'
#' @return A list with class \code{glm_context}.
#' @keywords internal
#' @noRd
glm_context <- function(X = NULL, Y = NULL,
                        proj = NULL,
                        phi_hat = NULL,
                        sigma_robust_scale = NULL,
                        robust_weights = NULL) {
  if (!is.null(X) && !is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (!is.null(Y) && !is.matrix(Y)) {
    Y <- as.matrix(Y)
  }

  structure(list(
    X = X,
    Y = Y,
    proj = proj,
    phi_hat = phi_hat,
    sigma_robust_scale = sigma_robust_scale,
    robust_weights = robust_weights
  ), class = c("glm_context", "list"))
}

#' Test if object is a glm_context
#'
#' @param x Object to test
#' @return Logical TRUE if object inherits from \code{glm_context}
#' @keywords internal
#' @noRd
is.glm_context <- function(x) {
  inherits(x, "glm_context")
}
