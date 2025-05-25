#' Confound-Free ALS HRF Estimation Engine
#'
#' Internal helper implementing the CF-ALS algorithm described in
#' `data-raw/CFALS_proposal.md`. Inputs should already be projected
#' to the confound null space.
#'
#' @param X_list_proj list of k design matrices (n x d each)
#' @param Y_proj numeric matrix of projected BOLD data (n x v)
#' @param lambda_b ridge penalty for beta-update
#' @param lambda_h ridge penalty for h-update
#' @param fullXtX_flag logical; if TRUE use cross-condition terms in h-update
#' @param h_ref_shape_norm optional reference HRF shape for sign alignment
#' @param max_alt number of alternating updates after initialization
#' @param epsilon_svd tolerance for singular value screening
#' @param epsilon_scale tolerance for scale in identifiability step
#' @return list with matrices `h` (d x v) and `beta` (k x v). The
#'   matrix `h` has an attribute `"iterations"` recording the number
#'   of alternating updates performed.
#' @keywords internal
#' @noRd
cf_als_engine <- function(X_list_proj, Y_proj,
                          lambda_b = 10,
                          lambda_h = 1,
                          fullXtX_flag = FALSE,
                          h_ref_shape_norm = NULL,
                          max_alt = 1,
                          epsilon_svd = 1e-8,
                          epsilon_scale = 1e-8) {
  stopifnot(is.list(X_list_proj), length(X_list_proj) >= 1)
  n <- nrow(Y_proj)
  v <- ncol(Y_proj)
  d <- ncol(X_list_proj[[1]])
  k <- length(X_list_proj)
  for (X in X_list_proj) {
    if (nrow(X) != n) stop("Design matrices must have same rows as Y_proj")
    if (ncol(X) != d) stop("All design matrices must have the same column count")
  }
  if (lambda_b < 0 || lambda_h < 0) {
    stop("lambda_b and lambda_h must be non-negative")
  }

  cholSolve <- function(M, b, eps = max(epsilon_svd, epsilon_scale)) {
    L <- tryCatch(chol(M),
                  error = function(e) chol(M + eps * diag(nrow(M))))
    backsolve(L, forwardsolve(t(L), b))
  }

  init <- ls_svd_engine(X_list_proj, Y_proj,
                        lambda_init = 0,
                        h_ref_shape_norm = h_ref_shape_norm,
                        epsilon_svd = epsilon_svd,
                        epsilon_scale = epsilon_scale)
  h_current <- init$h
  b_current <- init$beta

  XtX_list <- lapply(X_list_proj, crossprod)
  XtY_list <- lapply(X_list_proj, function(X) crossprod(X, Y_proj))

  if (fullXtX_flag) {
    XtX_full_list <- matrix(vector("list", k * k), k, k)
    for (l in seq_len(k)) {
      for (m in seq_len(k)) {
        XtX_full_list[[l, m]] <- crossprod(X_list_proj[[l]], X_list_proj[[m]])
      }
    }
  } else {
    XtX_full_list <- NULL
  }

  for (iter in seq_len(max_alt)) {
    for (vx in seq_len(v)) {
      h_vx <- h_current[, vx]
      DhTy_vx <- vapply(seq_len(k), function(c)
        crossprod(h_vx, XtY_list[[c]][, vx]), numeric(1))
      G_vx <- matrix(0.0, k, k)
      for (l in seq_len(k)) {
        if (fullXtX_flag) {
          for (m in seq_len(k)) {
            term <- XtX_full_list[[l, m]]
            G_vx[l, m] <- crossprod(h_vx, term %*% h_vx)
          }
        } else {
          G_vx[l, l] <- crossprod(h_vx, XtX_list[[l]] %*% h_vx)
        }
      }
      b_current[, vx] <- cholSolve(G_vx + lambda_b * diag(k), DhTy_vx)
    }

    for (vx in seq_len(v)) {
      b_vx <- b_current[, vx]
      lhs <- lambda_h * diag(d)
      rhs <- numeric(d)
      for (l in seq_len(k)) {
        rhs <- rhs + b_vx[l] * XtY_list[[l]][, vx]
        if (fullXtX_flag) {
          for (m in seq_len(k)) {
            lhs <- lhs + b_vx[l] * b_vx[m] * XtX_full_list[[l, m]]
          }
        } else {
          lhs <- lhs + b_vx[l]^2 * XtX_list[[l]]
        }
      }
      h_current[, vx] <- cholSolve(lhs, rhs)
    }
  }

  scl <- apply(abs(h_current), 2, max)
  flip <- rep(1.0, v)
  if (!is.null(h_ref_shape_norm)) {
    align <- colSums(h_current * h_ref_shape_norm)
    flip[align < 0 & scl > epsilon_scale] <- -1.0
  }
  eff_scl <- pmax(scl, epsilon_scale)
  h_final <- sweep(h_current, 2, flip / eff_scl, "*")
  b_final <- sweep(b_current, 2, flip * eff_scl, "*")
  zero_idx <- scl <= epsilon_scale
  if (any(zero_idx)) {
    h_final[, zero_idx] <- 0
    b_final[, zero_idx] <- 0
  }

  attr(h_final, "iterations") <- max_alt
  list(h = h_final, beta = b_final)
}

