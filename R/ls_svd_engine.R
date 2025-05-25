#' LS+SVD Rank-1 HRF Estimation Engine
#'
#' Internal helper implementing the LS+SVD algorithm described in
#' `data-raw/LSS+SVD_proposal.md`. The design matrices should already
#' be projected to the confound null space.
#'
#' @param X_list_proj list of k design matrices (n x d each)
#' @param Y_proj      numeric matrix of projected BOLD data (n x v)
#' @param lambda_init ridge penalty for initial GLM solve
#' @param h_ref_shape_norm optional reference HRF shape for sign alignment
#' @param svd_backend currently ignored, placeholder for future backends
#' @param epsilon_svd tolerance for singular value screening
#' @param epsilon_scale tolerance for scale in identifiability step
#' @return list with matrices `h` (d x v), `beta` (k x v) and
#'         `Gamma_hat` (d*k x v)
#' @keywords internal
#' @noRd
ls_svd_engine <- function(X_list_proj, Y_proj, lambda_init = 1, 
                          h_ref_shape_norm = NULL,
                          svd_backend = c("base_R"),
                          epsilon_svd = 1e-8,
                          epsilon_scale = 1e-8) {
  svd_backend <- match.arg(svd_backend)
  stopifnot(is.list(X_list_proj), length(X_list_proj) >= 1)
  n <- nrow(Y_proj)
  v <- ncol(Y_proj)
  d <- ncol(X_list_proj[[1]])
  k <- length(X_list_proj)
  for (X in X_list_proj) {
    if (nrow(X) != n) stop("Design matrices must have same rows as Y_proj")
    if (ncol(X) != d) stop("All design matrices must have the same column count")
  }
  if (!is.null(h_ref_shape_norm) && length(h_ref_shape_norm) != d)
    stop("`h_ref_shape_norm` must have length d")

  cholSolve <- function(M, B, eps = max(epsilon_svd, epsilon_scale)) {
    L <- tryCatch(chol(M),
                  error = function(e) chol(M + eps * diag(nrow(M))))
    backsolve(L, forwardsolve(t(L), B))
  }

  Xbig <- do.call(cbind, X_list_proj)
  XtX  <- crossprod(Xbig)
  Xty  <- crossprod(Xbig, Y_proj)
  XtX_ridge <- XtX + lambda_init * diag(ncol(Xbig))
  Gamma_hat <- cholSolve(XtX_ridge, Xty)

  H_out <- matrix(0.0, d, v)
  B_out <- matrix(0.0, k, v)

  for (vx in seq_len(v)) {
    G_vx <- matrix(Gamma_hat[, vx], nrow = d, ncol = k)
    sv <- svd(G_vx, nu = 1, nv = 1)
    if (length(sv$d) && sv$d[1] > epsilon_svd) {
      s1 <- sqrt(sv$d[1])
      H_out[, vx] <- sv$u[, 1] * s1
      B_out[, vx] <- sv$v[, 1] * s1
    }
  }

  scl <- apply(abs(H_out), 2, max)
  flip <- rep(1.0, v)
  if (!is.null(h_ref_shape_norm)) {
    align <- as.numeric(crossprod(h_ref_shape_norm, H_out))
    flip[align < 0 & scl > epsilon_scale] <- -1.0
  }
  eff_scl <- pmax(scl, epsilon_scale)
  H_final <- sweep(H_out, 2, flip / eff_scl, "*")
  B_final <- sweep(B_out, 2, flip * eff_scl, "*")
  zero_idx <- scl <= epsilon_scale
  if (any(zero_idx)) {
    H_final[, zero_idx] <- 0
    B_final[, zero_idx] <- 0
  }

  dimnames(H_final) <- list(NULL, colnames(Y_proj))
  dimnames(B_final) <- list(names(X_list_proj), colnames(Y_proj))

  list(h = H_final, beta = B_final, Gamma_hat = Gamma_hat)
}

