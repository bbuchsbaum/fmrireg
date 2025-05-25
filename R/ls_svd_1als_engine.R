#' LS+SVD+1ALS Rank-1 HRF Estimation Engine
#'
#' Internal helper implementing the LS+SVD+1ALS algorithm described in
#' `data-raw/LSS+SVD_proposal.md`. Inputs should be projected to the
#' confound null space.
#'
#' @param X_list_proj list of k design matrices (n x d each)
#' @param Y_proj numeric matrix of projected BOLD data (n x v)
#' @param lambda_init ridge penalty for initial GLM solve
#' @param lambda_b ridge penalty for \eqn{\beta}-update
#' @param lambda_h ridge penalty for \eqn{h}-update
#' @param fullXtX_flag logical; if TRUE use full cross-terms in h-update
#' @param h_ref_shape_norm optional reference HRF shape for sign alignment
#' @param svd_backend backend for SVD in the initialization step
#' @param epsilon_svd tolerance for singular value screening
#' @param epsilon_scale tolerance for scale in identifiability step
#' @return list with matrices `h` (d x v), `beta` (k x v) and the
#'         initial estimates `h_ls_svd`, `beta_ls_svd`
#' @keywords internal
#' @noRd
ls_svd_1als_engine <- function(X_list_proj, Y_proj,
                               lambda_init = 1,
                               lambda_b = 10,
                               lambda_h = 1,
                               fullXtX_flag = FALSE,
                               h_ref_shape_norm = NULL,
                               svd_backend = c("base_R"),
                               epsilon_svd = 1e-8,
                               epsilon_scale = 1e-8) {

  if (lambda_init < 0 || lambda_b < 0 || lambda_h < 0)
    stop("Lambdas must be non-negative")
  stopifnot(is.list(X_list_proj), length(X_list_proj) >= 1)
  d <- ncol(X_list_proj[[1]])
  if (!is.null(h_ref_shape_norm) && length(h_ref_shape_norm) != d)
    stop("`h_ref_shape_norm` must be length d")

  cholSolve <- function(M, B) {
    R <- tryCatch(chol(M),
                  error = function(e) chol(M + 1e-6 * diag(nrow(M))))
    backsolve(R, forwardsolve(t(R), B))
  }
  init <- ls_svd_engine(X_list_proj, Y_proj,
                        lambda_init = lambda_init,
                        h_ref_shape_norm = h_ref_shape_norm,
                        svd_backend = svd_backend,
                        epsilon_svd = epsilon_svd,
                        epsilon_scale = epsilon_scale)

  h_current <- init$h
  b_current <- init$beta
  k <- length(X_list_proj)
  d <- ncol(X_list_proj[[1]])
  v <- ncol(Y_proj)

  XtX_list <- lapply(X_list_proj, crossprod)
  XtY_list <- lapply(X_list_proj, function(X) crossprod(X, Y_proj))

  XtX_full_list <- NULL
  if (fullXtX_flag) {
    XtX_full_list <- matrix(vector("list", k * k), k, k)
    for (l in seq_len(k)) {
      for (m in seq_len(k)) {
        XtX_full_list[[l, m]] <- crossprod(X_list_proj[[l]], X_list_proj[[m]])
      }
    }
  }

  H_als <- matrix(0.0, d, v)
  B_als <- matrix(0.0, k, v)

  for (vx in seq_len(v)) {
    h_vx <- h_current[, vx]
    DhTy_vx <- vapply(seq_len(k), function(c)
      crossprod(h_vx, XtY_list[[c]][, vx]), numeric(1))
    if (fullXtX_flag) {
      G_vx <- matrix(0.0, k, k)
      for (l in seq_len(k)) {
        for (m in seq_len(k)) {
          G_vx[l, m] <- crossprod(h_vx, XtX_full_list[[l, m]] %*% h_vx)
        }
      }
    } else {
      diag_vals <- vapply(seq_len(k), function(c)
        crossprod(h_vx, XtX_list[[c]] %*% h_vx), numeric(1))
      G_vx <- diag(diag_vals, k)
    }
    B_als[, vx] <- cholSolve(G_vx + lambda_b * diag(k), DhTy_vx)
  }

  b_current <- B_als

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
    H_als[, vx] <- cholSolve(lhs, rhs)
  }

  scl <- apply(abs(H_als), 2, max)
  flip <- rep(1.0, v)
  if (!is.null(h_ref_shape_norm)) {
    align <- as.numeric(crossprod(h_ref_shape_norm, H_als))
    flip[align < 0 & scl > epsilon_scale] <- -1.0
  }
  eff_scl <- pmax(scl, epsilon_scale)
  H_final <- sweep(H_als, 2, flip / eff_scl, "*")
  B_final <- sweep(B_als, 2, flip * eff_scl, "*")
  zero_idx <- scl <= epsilon_scale
  if (any(zero_idx)) {
    H_final[, zero_idx] <- 0
    B_final[, zero_idx] <- 0
  }

  dimnames(H_final) <- list(NULL, colnames(Y_proj))
  dimnames(B_final) <- list(names(X_list_proj), colnames(Y_proj))

  list(h = H_final, beta = B_final,
       h_ls_svd = init$h, beta_ls_svd = init$beta)
}

