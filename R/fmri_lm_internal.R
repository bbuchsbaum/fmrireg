# Internal Utilities for fMRI Linear Models
# Low-level utilities used throughout the fmri_lm implementation


#' Check if object is a formula
#' @keywords internal
#' @noRd
is.formula <- function(x) {
  inherits(x, "formula")
}

.rank_deficiency_message <- function(rank, p, aliased, varnames = NULL, context = "Design matrix") {
  aliased_names <- if (!is.null(varnames) && length(varnames) >= max(aliased, 0L)) {
    varnames[aliased]
  } else {
    paste0("#", aliased)
  }

  if (length(aliased_names) > 8L) {
    aliased_names <- c(aliased_names[seq_len(8L)], "...")
  }

  sprintf(
    "%s is rank deficient: rank %d < %d columns. Aliased columns: %s",
    context,
    rank,
    p,
    paste(aliased_names, collapse = ", ")
  )
}

.design_rank_info <- function(X, tol = 1e-7) {
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  p <- ncol(X)
  n <- nrow(X)
  if (p == 0L || n == 0L) {
    stop(".fast_preproject: design matrix must have at least one row and one column", call. = FALSE)
  }

  qrx <- qr(X, tol = tol, LAPACK = FALSE)
  rank <- qrx$rank
  pivot <- qrx$pivot %||% seq_len(p)
  estimable <- if (rank > 0L) pivot[seq_len(rank)] else integer(0)
  aliased <- if (rank < p) pivot[seq.int(rank + 1L, p)] else integer(0)

  list(
    qr = qrx,
    rank = rank,
    pivot = pivot,
    estimable = estimable,
    aliased = aliased,
    is_full_rank = (rank == p),
    tol = tol
  )
}

.attach_rank_attrs <- function(x, rank_info) {
  attr(x, "rank") <- rank_info$rank
  attr(x, "is_full_rank") <- rank_info$is_full_rank
  attr(x, "estimable") <- rank_info$estimable
  attr(x, "aliased") <- rank_info$aliased
  attr(x, "rank_tolerance") <- rank_info$tol
  x
}

#' Fast Pre-projection of Design Matrix
#'
#' @description
#' Computes the projection components needed for fast least squares: the
#' tolerance-aware QR decomposition, (X'X)^-1 and the residual df.
#'
#' The p x n pseudo-inverse is *not* formed here: building it costs an n x n
#' solve, and the solver never needs it (betas and RSS come from one
#' `qr.qty()` pass in `solve_glm_core()`). Callers that do need it read
#' `proj$Pinv`, which the `$` method for `fmri_preproject` computes on demand
#' from the stored QR (`.preproject_pinv()`). Only `$` is lazy:
#' `proj[["Pinv"]]` is `NULL` unless the SVD fallback stored it eagerly.
#'
#' @param X Design matrix
#' @param tol Rank tolerance for the LINPACK QR.
#' @return An `fmri_preproject` list containing:
#'   - qr: QR decomposition of X
#'   - XtXinv: (X'X)^-1 matrix (zero rows/columns for aliased coefficients)
#'   - dfres: Residual degrees of freedom (n - rank)
#'   - Pinv: pseudo-inverse, computed on `$` access
#' @keywords internal
#' @noRd
.fast_preproject <- function(X, tol = 1e-7) {
  # Ensure X is a numeric matrix
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  if (anyNA(X) || any(!is.finite(X))) {
    stop(".fast_preproject: design matrix contains NA/Inf values", call. = FALSE)
  }
  if (!is.numeric(X)) {
    orig_colnames <- colnames(X)
    X <- matrix(as.numeric(X), nrow = nrow(X), ncol = ncol(X))
    if (!is.null(orig_colnames)) {
      colnames(X) <- orig_colnames
    }
  }
  
  # QR decomposition for rank detection and stable computation. Use the
  # LINPACK path here because its rank is tolerance-aware and matches lm.fit.
  rank_info <- .design_rank_info(X, tol = tol)
  qr_decomp <- rank_info$qr
  rank <- rank_info$rank
  p <- ncol(X)
  n <- nrow(X)
  pivot <- rank_info$pivot
  estimable <- rank_info$estimable
  aliased <- rank_info$aliased
  
  # Check for rank deficiency
  if (rank < p) {
    warning(.rank_deficiency_message(rank, p, aliased, colnames(X)), call. = FALSE)
  }
  
  # Pinv stays NULL (computed lazily from the QR) unless the SVD fallback
  # below has to produce it.
  Pinv <- NULL
  if (rank > 0L) {
    XtXinv <- tryCatch({
      R <- qr.R(qr_decomp)[seq_len(rank), seq_len(rank), drop = FALSE]
      Rinv <- backsolve(R, diag(rank))
      out <- matrix(0, nrow = p, ncol = p)
      out[estimable, estimable] <- tcrossprod(Rinv)
      out
    }, error = function(e) NULL)

    if (is.null(XtXinv)) {
      X_est <- X[, estimable, drop = FALSE]
      svd_result <- svd(X_est)
      d <- svd_result$d
      sv_tol <- tol * max(d)
      pos <- d > sv_tol
      U <- svd_result$u[, pos, drop = FALSE]
      V <- svd_result$v[, pos, drop = FALSE]
      D_inv <- diag(1/d[pos], nrow = sum(pos))
      Pinv_est <- V %*% D_inv %*% t(U)
      XtXinv_est <- V %*% D_inv^2 %*% t(V)
      Pinv <- matrix(0, nrow = p, ncol = n)
      XtXinv <- matrix(0, nrow = p, ncol = p)
      Pinv[estimable, ] <- Pinv_est
      XtXinv[estimable, estimable] <- XtXinv_est
    }
  } else {
    Pinv <- matrix(0, nrow = p, ncol = n)
    XtXinv <- matrix(0, nrow = p, ncol = p)
  }

  if (!is.null(colnames(X))) {
    dimnames(XtXinv) <- list(colnames(X), colnames(X))
    if (!is.null(Pinv)) rownames(Pinv) <- colnames(X)
  }

  # Return everything needed for fast operations
  structure(
    list(
      qr = qr_decomp,
      Pinv = if (is.null(Pinv)) NULL else .attach_rank_attrs(Pinv, rank_info),
      XtXinv = .attach_rank_attrs(XtXinv, rank_info),
      dfres = n - rank,
      rank = rank,
      is_full_rank = (rank == p),
      estimable = estimable,
      aliased = aliased,
      rank_tolerance = tol,
      varnames = colnames(X)
    ),
    class = c("fmri_preproject", "list")
  )
}

#' Pseudo-inverse of a pre-projected design
#'
#' Computes the p x n matrix P with P %*% Y equal to the least-squares
#' coefficients, as R^-1 Q1' from the stored QR (rows of aliased
#' coefficients are zero). This avoids the n x n identity solve that
#' `qr.coef(qr, diag(n))` performs.
#' @keywords internal
#' @noRd
.preproject_pinv <- function(proj) {
  eager <- .subset2(proj, "Pinv")
  if (!is.null(eager)) return(eager)
  qrx <- .subset2(proj, "qr")
  rank <- .subset2(proj, "rank")
  estimable <- .subset2(proj, "estimable")
  n <- nrow(qrx$qr)
  p <- ncol(qrx$qr)
  Pinv <- matrix(0, nrow = p, ncol = n)
  if (rank > 0L) {
    keep <- seq_len(rank)
    Q1 <- qr.Q(qrx)[, keep, drop = FALSE]
    R <- qr.R(qrx)[keep, keep, drop = FALSE]
    Pinv[estimable, ] <- backsolve(R, t(Q1))
  }
  varnames <- .subset2(proj, "varnames")
  if (!is.null(varnames)) rownames(Pinv) <- varnames
  attr(Pinv, "rank") <- rank
  attr(Pinv, "is_full_rank") <- .subset2(proj, "is_full_rank")
  attr(Pinv, "estimable") <- estimable
  attr(Pinv, "aliased") <- .subset2(proj, "aliased")
  attr(Pinv, "rank_tolerance") <- .subset2(proj, "rank_tolerance")
  Pinv
}

#' Lazy field access for pre-projected designs
#'
#' `proj$Pinv` computes the pseudo-inverse on demand; every other field is
#' returned as stored (exact matching).
#' @param x An `fmri_preproject` object.
#' @param name Field name.
#' @keywords internal
#' @noRd
#' @exportS3Method base::`$`
`$.fmri_preproject` <- function(x, name) {
  if (identical(name, "Pinv")) {
    return(.preproject_pinv(x))
  }
  .subset2(x, name)
}

#' Fast Matrix-based Linear Model
#'
#' @description
#' Low-level function for fast least squares computation.
#' This is being phased out in favor of solve_glm_core.
#'
#' @param X Design matrix
#' @param Y Response matrix (n x v)
#' @param proj Pre-computed projection from .fast_preproject
#' @param return_fitted Whether to return fitted values
#' @return List with regression results
#' @keywords internal
#' @noRd
.fast_lm_matrix <- function(X, Y, proj, return_fitted = FALSE) {
  ctx <- glm_context(X = X, Y = Y, proj = proj)
  res <- solve_glm_core(ctx, return_fitted = return_fitted)
  res$sigma <- sqrt(res$sigma2)
  res
}

#' Classify a response vector before voxelwise fitting
#'
#' Degenerate and non-finite responses are retained in their original voxel
#' positions, but they do not enter covariance or test-statistic calculations.
#' Keeping this classification in one helper makes the sequential and future
#' execution paths obey the same contract.
#'
#' @keywords internal
#' @noRd
.fmri_lm_voxel_status <- function(y) {
  if (any(!is.finite(y))) {
    return("nonfinite")
  }
  if (all(y == 0)) {
    return("all_zero")
  }
  if (length(y) < 2L || stats::var(y) <= .Machine$double.eps) {
    return("constant")
  }
  "ok"
}

#' Inference df for one voxel
#'
#' Without robust weights this is the residual df (whitening already accounts
#' for AR structure, see `calculate_effective_df()`). With robust weights it is
#' the single robust formula `.robust_effective_df()`: sum(w) - rank, with the
#' design rank (not its width) so aliased columns are not subtracted.
#' @keywords internal
#' @noRd
.voxel_inference_df <- function(dfres, rank, robust_weights = NULL) {
  if (is.null(robust_weights)) return(dfres)
  .robust_effective_df(robust_weights, rank)
}

#' Stack per-voxel covariance matrices
#'
#' Returns the columns `vec(XtXinv_v)` (p^2 x V; NA for voxels that are
#' skipped) and the mask of usable voxels: a p x p matrix and a finite scale.
#' @keywords internal
#' @noRd
.voxelwise_cov_stack <- function(XtXinv_list, scale, p, check_dim = TRUE) {
  V <- length(XtXinv_list)
  ok <- vapply(seq_len(V), function(v) {
    M <- XtXinv_list[[v]]
    !is.null(M) && is.finite(scale[v]) &&
      (!check_dim || (is.matrix(M) && all(dim(M) == c(p, p))))
  }, logical(1))
  S <- matrix(NA_real_, p * p, V)
  if (any(ok)) {
    S[, ok] <- vapply(XtXinv_list[ok], function(M) as.numeric(M), numeric(p * p))
  }
  list(S = S, ok = ok)
}

#' Per-voxel quadratic forms l' M_v k for stacked M_v
#' @keywords internal
#' @noRd
.voxelwise_bilinear <- function(S, l, k = l) {
  as.numeric(crossprod(as.numeric(outer(l, k)), S))
}

#' Beta statistics when each voxel has its own projection matrix
#'
#' Helper for voxelwise AR fitting where every voxel yields a distinct
#' whitened design matrix.  Computes standard errors and t statistics using
#' the per-voxel \code{XtXinv} matrices, for all voxels at once.
#'
#' @keywords internal
#' @noRd
beta_stats_matrix_voxelwise <- function(Betas, XtXinv_list, sigma, dfres,
                                        varnames,
                                        robust_weights_list = NULL,
                                        ar_order = 0) {
  V <- ncol(Betas)
  p <- nrow(Betas)

  est_mat  <- matrix(NA_real_, V, p)
  se_mat   <- matrix(NA_real_, V, p)
  t_mat    <- matrix(NA_real_, V, p)
  prob_mat <- matrix(NA_real_, V, p)

  stack <- .voxelwise_cov_stack(XtXinv_list, sigma, p)
  ok <- stack$ok
  if (any(ok)) {
    diag_idx <- (seq_len(p) - 1L) * p + seq_len(p)
    se_scal <- t(sqrt(stack$S[diag_idx, ok, drop = FALSE]))   # n_ok x p
    df_eff <- if (is.null(robust_weights_list)) {
      rep(dfres, sum(ok))
    } else {
      vapply(which(ok), function(v) {
        .voxel_inference_df(dfres, attr(XtXinv_list[[v]], "rank") %||% p,
                            robust_weights_list[[v]])
      }, numeric(1))
    }
    se_ok <- se_scal * sigma[ok]
    est_ok <- t(Betas[, ok, drop = FALSE])
    t_ok <- ifelse(abs(se_ok) < .Machine$double.eps^0.5, 0, est_ok / se_ok)
    est_mat[ok, ] <- est_ok
    se_mat[ok, ] <- se_ok
    t_mat[ok, ] <- t_ok
    prob_mat[ok, ] <- 2 * pt(-abs(t_ok), df_eff)
  }

  colnames(est_mat)  <- varnames
  colnames(se_mat)   <- varnames
  colnames(t_mat)    <- varnames
  colnames(prob_mat) <- varnames

  tibble::tibble(
    type = "beta",
    name = "parameter_estimates",
    stat_type = "tstat",
    df.residual = dfres,
    conmat = list(NULL),
    colind = list(NULL),
    data = list(tibble::tibble(
      estimate = list(est_mat),
      se = list(se_mat),
      stat = list(t_mat),
      prob = list(prob_mat),
      sigma = list(sigma)
    ))
  )
}

#' Per-voxel F quadratic forms u_v' (L M_v L')^-1 u_v
#'
#' Batched Cholesky across voxels (loops run over the r x r entries, not over
#' voxels). Voxels whose factorisation is not clearly well conditioned fall
#' back to the per-voxel `solve()` used by `.fast_F_contrast()`, including its
#' warning for a singular hypothesis matrix.
#' @keywords internal
#' @noRd
.voxelwise_F_quadform <- function(S, L, U) {
  r <- nrow(L)
  nv <- ncol(U)
  # Mv[[i]][[j]][v] = (L M_v L')_{ij}
  Mv <- lapply(seq_len(r), function(i) {
    lapply(seq_len(r), function(j) .voxelwise_bilinear(S, L[i, ], L[j, ]))
  })
  if (r == 1L) {
    m <- Mv[[1L]][[1L]]
    qf <- (U[1L, ] * (1 / m)) * U[1L, ]
    bad <- !is.finite(m) | m == 0
  } else {
    Lc <- vector("list", r)
    for (i in seq_len(r)) Lc[[i]] <- vector("list", r)
    scale <- do.call(pmax, lapply(seq_len(r), function(i) abs(Mv[[i]][[i]])))
    bad <- !is.finite(scale) | scale <= 0
    z <- vector("list", r)
    qf <- numeric(nv)
    for (j in seq_len(r)) {
      d <- Mv[[j]][[j]]
      if (j > 1L) for (k in seq_len(j - 1L)) d <- d - Lc[[j]][[k]]^2
      bad <- bad | !is.finite(d) | d <= scale * 1e-6
      Lc[[j]][[j]] <- sqrt(pmax(d, 0))
      if (j < r) {
        for (i in (j + 1L):r) {
          a <- Mv[[i]][[j]]
          if (j > 1L) for (k in seq_len(j - 1L)) a <- a - Lc[[i]][[k]] * Lc[[j]][[k]]
          Lc[[i]][[j]] <- a / Lc[[j]][[j]]
        }
      }
      zj <- U[j, ]
      if (j > 1L) for (k in seq_len(j - 1L)) zj <- zj - Lc[[j]][[k]] * z[[k]]
      z[[j]] <- zj / Lc[[j]][[j]]
      qf <- qf + z[[j]]^2
    }
  }
  for (v in which(bad)) {
    M <- matrix(vapply(seq_len(r * r), function(ij) {
      Mv[[(ij - 1L) %% r + 1L]][[(ij - 1L) %/% r + 1L]][v]
    }, numeric(1)), r, r)
    Cinv <- tryCatch(solve(M), error = function(e) {
      warning(paste("Singular matrix in F-contrast computation (L(X'X)^-1 L'). Details:", e$message))
      matrix(NaN, nrow = r, ncol = r)
    })
    qf[v] <- if (any(is.nan(Cinv))) NaN else {
      tmp <- crossprod(U[, v], Cinv)
      sum(tmp * U[, v])
    }
  }
  qf
}

#' Contrast statistics with voxelwise projection matrices
#'
#' Computes t and F contrasts when each voxel has a distinct \code{XtXinv},
#' for all voxels at once: the per-voxel variances \eqn{l' M_v l} (and the
#' r x r matrices \eqn{L M_v L'}) are read off the stacked \eqn{M_v} with one
#' matrix product per contrast. Values equal those of calling
#' `.fast_t_contrast()` / `.fast_F_contrast()` voxel by voxel.
#'
#' @keywords internal
#' @noRd
fit_lm_contrasts_voxelwise <- function(Betas, sigma2, XtXinv_list,
                                       conlist, fconlist, dfres,
                                       robust_weights_list = NULL,
                                       ar_order = 0) {
  p <- nrow(Betas)
  V <- ncol(Betas)

  results <- list()
  stack <- .voxelwise_cov_stack(XtXinv_list, sigma2, p, check_dim = FALSE)
  ok <- stack$ok
  S_ok <- stack$S[, ok, drop = FALSE]
  B_ok <- Betas[, ok, drop = FALSE]
  # Inference df exactly as .fast_t_contrast()/.fast_F_contrast() derive it.
  df_ok <- if (is.null(robust_weights_list)) {
    rep(.effective_df_from_rank(dfres, p, integer(0), NULL, ar_order), sum(ok))
  } else {
    vapply(which(ok), function(v) {
      .effective_df_from_rank(dfres, p, integer(0), robust_weights_list[[v]],
                              ar_order)
    }, numeric(1))
  }

  for (nm in names(conlist)) {
    l <- conlist[[nm]]
    colind <- attr(l, "colind")
    full_l <- matrix(0, nrow = 1, ncol = p)
    full_l[, colind] <- l

    est <- se <- stat <- prob <- sigma_out <- rep(NA_real_, V)
    if (any(ok)) {
      e <- drop(full_l %*% B_ok)
      s2 <- .voxelwise_bilinear(S_ok, full_l[1L, ])
      se_ok <- sqrt(pmax(0, s2 * sigma2[ok]))
      tval <- .safe_tstat(e, se_ok)
      est[ok] <- e
      se[ok] <- se_ok
      stat[ok] <- tval
      prob[ok] <- 2 * pt(-abs(tval), df_ok)
      sigma_out[ok] <- sqrt(pmax(0, sigma2[ok]))
    }

    results[[nm]] <- tibble::tibble(
      type = "contrast",
      name = nm,
      stat_type = "tstat",
      df.residual = dfres,
      conmat = list(l),
      colind = list(colind),
      data = list(tibble::tibble(
        estimate = est,
        se = se,
        stat = stat,
        prob = prob,
        sigma = sigma_out
      ))
    )
  }

  for (nm in names(fconlist)) {
    L <- fconlist[[nm]]
    colind <- attr(L, "colind")
    # Normalise to hypotheses-in-rows so cell-oriented (e.g. oneway, k x (k-1))
    # weights build the correct r x p matrix; see .orient_fcontrast().
    Cw <- .orient_fcontrast(L, colind)
    full_L <- matrix(0, nrow = nrow(Cw), ncol = p)
    full_L[, colind] <- Cw
    r <- nrow(full_L)

    est <- se <- stat <- prob <- rep(NA_real_, V)
    if (any(ok)) {
      U <- full_L %*% B_ok
      qf <- .voxelwise_F_quadform(S_ok, full_L, U)
      estimate <- qf / r
      s2 <- sigma2[ok]
      Fval <- ifelse(!(s2 > 0) | is.nan(qf), NaN, estimate / s2)
      est[ok] <- estimate
      se[ok] <- s2
      stat[ok] <- Fval
      prob[ok] <- pf(Fval, r, df_ok, lower.tail = FALSE)
    }

    results[[nm]] <- tibble::tibble(
      type = "Fcontrast",
      name = nm,
      stat_type = "Fstat",
      df.residual = dfres,
      conmat = list(L),
      colind = list(colind),
      data = list(tibble::tibble(
        estimate = est,
        se = se,
        stat = stat,
        prob = prob
      ))
    )
  }

  results
}
