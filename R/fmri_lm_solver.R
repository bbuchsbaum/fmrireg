#' Core GLM solver
#'
#' Internal low-level solver used by the linear modelling code. It operates
#' on a `glm_context` containing design matrix `X`, data matrix `Y` and the
#' preprojection object `proj` from `.fast_preproject()`. If
#' `glm_ctx$robust_weights` is `NULL` an ordinary (or prewhitened) least
#' squares fit is performed. When `robust_weights` is supplied the matrices in
#' the context are assumed to already have the weights applied and weighted
#' least squares is performed.
#'
#' @param glm_ctx A `glm_context` object.
#' @param return_fitted Logical; return fitted values.
#'
#' @return List with components `betas`, `rss`, `sigma2` and optionally
#'   `fitted` when `return_fitted` is `TRUE`.
#' @keywords internal
#' @noRd
solve_glm_core <- function(glm_ctx, return_fitted = FALSE) {
  stopifnot(is.glm_context(glm_ctx))
  X <- glm_ctx$X
  Y <- glm_ctx$Y
  proj <- glm_ctx$proj

  if (is.null(X) || is.null(Y) || is.null(proj)) {
    stop("solve_glm_core: context must contain 'X', 'Y' and 'proj'")
  }

  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)

  if (nrow(X) != nrow(Y)) {
    stop("solve_glm_core: X and Y dimensions do not match")
  }

  if (isTRUE(proj$is_full_rank) && !is.null(proj$qr)) {
    # Full-rank path: solve from QR directly for numerical stability.
    betas <- tryCatch(
      qr.coef(proj$qr, Y),
      error = function(e) NULL
    )
    if (is.null(betas)) {
      if (is.null(proj$Pinv) || ncol(X) != nrow(proj$Pinv)) {
        stop("solve_glm_core: QR solve failed and projection matrix dimensions do not match")
      }
      betas <- proj$Pinv %*% Y
    }
  } else {
    if (is.null(proj$Pinv) || ncol(X) != nrow(proj$Pinv)) {
      stop("solve_glm_core: X and projection matrix dimensions do not match")
    }
    betas <- proj$Pinv %*% Y
  }

  if (return_fitted) {
    fitted <- X %*% betas
    residuals <- Y - fitted
    rss <- colSums(residuals^2)
  } else {
    # Memory-lean path: avoid allocating n x V fitted/residual matrices when
    # only RSS/sigma2 are needed.
    XtY <- crossprod(X, Y)
    yTy <- colSums(Y * Y)
    rss <- yTy - colSums(betas * XtY)
    rss <- pmax(rss, 0)
    fitted <- NULL
  }

  sigma2 <- rss / proj$dfres
  
  # Add rank information to output if available
  result <- list(
    betas = betas,
    rss = rss,
    sigma2 = sigma2,
    fitted = fitted,
    dfres = proj$dfres  # Always include dfres
  )
  
  # Add rank info if available in projection
  if (!is.null(proj$rank)) {
    result$rank <- proj$rank
    result$is_full_rank <- proj$is_full_rank
    
    # If rank deficient, mark coefficients as potentially unreliable
    if (!proj$is_full_rank) {
      attr(result$betas, "rank_deficient") <- TRUE
    }
  }
  
  result
}
