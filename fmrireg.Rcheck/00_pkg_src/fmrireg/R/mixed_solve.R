#' Mixed Model Solver using Rcpp and roptim
#'
#' This function solves a mixed model using Rcpp and roptim for optimization.
#' It estimates variance components in a mixed model, potentially speeding up computations compared to the pure R implementation.
#'
#' @param y Response vector.
#' @param Z Design matrix for random effects (default: identity matrix of size n).
#' @param K Kinship matrix (default: NULL).
#' @param X Design matrix for fixed effects (default: vector of ones).
#' @param method Optimization method, either "REML" or "ML" (default: "REML").
#' @param bounds Bounds for the optimizer (default: c(1e-9, 1e9)).
#' @param SE Logical, whether to return standard errors (default: FALSE).
#' @param return_Hinv Logical, whether to return the inverse of H (default: FALSE).
#' @return A list containing:
#'   \item{Vu}{Estimated variance component for random effects.}
#'   \item{Ve}{Estimated variance component for residuals.}
#'   \item{beta}{Estimated fixed effects coefficients.}
#'   \item{u}{Estimated random effects coefficients.}
#'   \item{LL}{Log-likelihood of the model.}
#'   \item{beta.SE}{Standard errors of fixed effects coefficients (if SE = TRUE).}
#'   \item{u.SE}{Standard errors of random effects coefficients (if SE = TRUE).}
#'   \item{Hinv}{Inverse of H (if return_Hinv = TRUE).}
#' @examples
#' \dontrun{
#' # Example usage with random data
#' set.seed(123)
#' n <- 100
#' y <- rnorm(n)
#' Z <- matrix(rnorm(n * 5), n, 5)
#' K <- diag(5)
#' X <- matrix(1, n, 1)
#' result <- mixed_solve_cpp(y, Z, K, X)
#' }
#' @export
mixed_solve_cpp <- function(y,
                            Z = NULL,
                            K = NULL,
                            X = NULL,
                            method = "REML",
                            bounds = c(1e-9, 1e9),
                            SE = FALSE,
                            return_Hinv = FALSE) {

  result <- mixed_solve_internal(y, Z, K, X, method, bounds, SE, return_Hinv)
  return(result)
}