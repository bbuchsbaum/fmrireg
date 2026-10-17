#' Core GLM solver
#'
#' Internal low-level solver used by the linear modelling code. It operates
#' on a `glm_context` containing design matrix `X`, data matrix `Y` and the
#' preprojection object `proj` from `.fast_preproject()`, and performs least
#' squares on them as given: callers whiten (AR) or weight (robust, volume
#' weights) `X` and `Y` beforehand.
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

  # Aliased (rank-deficient) coefficient positions: used both to zero only those
  # coefficients below and to mark the result for downstream NA handling.
  aliased <- as.integer(proj$aliased %||% attr(proj$XtXinv, "aliased") %||% integer(0))

  # Exact RSS from the QR, when available (set below).
  rss_qr <- NULL
  if (!is.null(proj$qr)) {
    # Solve from the same tolerance-aware QR used for rank diagnostics, in one
    # pass: Q'Y gives the coefficients (R^-1 of its first `rank` rows) and the
    # exact residual sum of squares (squared norm of its remaining rows). This
    # replaces qr.coef() plus the normal-equation identity yTy - b'X'y, which
    # loses precision on near-collinear designs with a large response offset.
    # Aliased (rank-deficient) coefficients are zero, so X %*% betas stays
    # finite; they are marked via attributes and NA'd downstream. Non-finite
    # betas arising from NA in Y (missing-data voxels) are left as NA so such
    # voxels surface as NA rather than as a 0-coefficient fit.
    qr_fit <- tryCatch({
      qrx <- proj$qr
      rnk <- as.integer(proj$rank %||% qrx$rank)
      n_obs <- nrow(qrx$qr)
      qty <- qr.qty(qrx, Y)
      if (is.null(dim(qty))) qty <- matrix(qty, ncol = ncol(Y))
      b <- matrix(0, nrow = ncol(qrx$qr), ncol = ncol(Y))
      if (rnk > 0L) {
        keep <- seq_len(rnk)
        R <- qr.R(qrx)[keep, keep, drop = FALSE]
        b[qrx$pivot[keep], ] <- backsolve(R, qty[keep, , drop = FALSE])
      }
      r_ss <- if (n_obs > rnk) {
        colSums(qty[seq.int(rnk + 1L, n_obs), , drop = FALSE]^2)
      } else {
        rep(0, ncol(Y))
      }
      list(betas = b, rss = r_ss)
    }, error = function(e) NULL)
    if (!is.null(qr_fit)) {
      betas <- qr_fit$betas
      rss_qr <- qr_fit$rss
    } else {
      if (is.null(proj$Pinv) || ncol(X) != nrow(proj$Pinv)) {
        stop("solve_glm_core: QR solve failed and projection matrix dimensions do not match")
      }
      betas <- proj$Pinv %*% Y
    }
    if (is.null(dim(betas))) {
      betas <- matrix(betas, ncol = 1L)
    }
    if (length(aliased) > 0L) {
      betas[aliased, ] <- 0
    }
  } else {
    if (is.null(proj$Pinv) || ncol(X) != nrow(proj$Pinv)) {
      stop("solve_glm_core: X and projection matrix dimensions do not match")
    }
    betas <- proj$Pinv %*% Y
  }

  if (!is.null(colnames(X)) && nrow(betas) == length(colnames(X))) {
    rownames(betas) <- colnames(X)
  }
  # Mark rank-deficient fits for downstream NA handling.
  if (length(aliased) > 0L || (!is.null(proj$rank) && !proj$is_full_rank)) {
    attr(betas, "rank_deficient") <- TRUE
    attr(betas, "aliased") <- aliased
    attr(betas, "estimable") <- as.integer(proj$estimable %||% integer(0))
  }

  if (return_fitted) {
    fitted <- X %*% betas
    rss <- rss_qr %||% colSums((Y - fitted)^2)
  } else if (!is.null(rss_qr)) {
    rss <- rss_qr
    fitted <- NULL
  } else {
    # Memory-lean path: avoid allocating n x V fitted/residual matrices when
    # only RSS/sigma2 are needed.
    XtY <- crossprod(X, Y)
    yTy <- colSums(Y * Y)
    rss <- yTy - colSums(betas * XtY)
    # yTy - b'X'y loses all precision when the fit is nearly exact: the RSS is
    # then below the rounding error of yTy and can come out as 0 or negative,
    # giving se = 0. Recompute those voxels from explicit residuals.
    inexact <- which(is.finite(rss) & rss <= yTy * 1e-8)
    if (length(inexact) > 0L) {
      resid <- Y[, inexact, drop = FALSE] - X %*% betas[, inexact, drop = FALSE]
      rss[inexact] <- colSums(resid^2)
    }
    rss <- pmax(rss, 0)
    fitted <- NULL
  }

  has_residual_df <- length(proj$dfres) == 1L &&
    is.finite(proj$dfres) && proj$dfres > 0
  sigma2 <- if (has_residual_df) {
    rss / proj$dfres
  } else {
    # Residual variance is undefined for saturated designs. Returning NA
    # explicitly avoids platform-dependent NaN/Inf results from 0/0 or x/0.
    rep(NA_real_, ncol(Y))
  }
  
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
    result$aliased <- aliased
    result$estimable <- as.integer(proj$estimable %||% integer(0))
  }
  
  result
}
