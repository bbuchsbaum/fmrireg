#' Build a time sketch matrix S (m x T)
#'
#' All sketches are normalised so that \eqn{E\|S r\|^2 = \|r\|^2}.
#'
#' @param Tlen Integer time length
#' @param ctrl List(method = "gaussian"|"countsketch"|"srht"|"ihs", m, iters)
#' @return A dense or sparse sketch matrix for Gaussian/CountSketch methods;
#'   `NULL` for `"srht"` and `"ihs"` because those methods are applied via plans/operators.
#' @keywords internal
make_time_sketch <- function(Tlen, ctrl) {
  stopifnot(is.list(ctrl))
  method <- match.arg(ctrl$method %||% "gaussian", c("gaussian", "countsketch", "srht", "ihs"))
  m <- ctrl$m
  stopifnot(!is.null(m), m > 0L, m <= Tlen)
  if (method == "gaussian") {
    S <- matrix(stats::rnorm(m * Tlen), nrow = m, ncol = Tlen) / sqrt(m)
  } else if (method == "countsketch") {
    rows <- sample.int(m, Tlen, replace = TRUE)
    signs <- sample(c(-1, 1), Tlen, replace = TRUE)
    S <- Matrix::sparseMatrix(i = rows, j = seq_len(Tlen), x = signs, dims = c(m, Tlen))
  } else {
    # SRHT / IHS do not materialize S; callers should use srht_apply()/ihs_latent_solve()
    S <- NULL
  }
  S
}

# --- SRHT helpers ---

#' Build SRHT plan
#'
#' The plan's `scale` is `1 / sqrt(m)`, which makes the sketch normalised
#' (\eqn{E\|S r\|^2 = \|r\|^2}) because the Walsh-Hadamard transform is
#' applied unnormalised.
#' @param Tlen Integer series length.
#' @param m Integer number of sketch rows.
#' @return A list with the sign, permutation, row-selection and scale
#'   components used by [srht_apply()].
#' @keywords internal
make_srht_plan <- function(Tlen, m) {
  stopifnot(m > 0L, m <= Tlen)
  list(T = Tlen, m = m,
       signs = sample(c(-1, 1), Tlen, replace = TRUE),
       perm  = sample.int(Tlen) - 1L,   # zero-based for C++
       rows  = sort(sample.int(Tlen, m, replace = FALSE)) - 1L,
       scale = 1 / sqrt(m))
}

#' Apply SRHT to a T x k matrix using a plan
#' @param M Numeric matrix with `plan$T` rows.
#' @param plan A plan from [make_srht_plan()].
#' @return The `plan$m` x `ncol(M)` sketched matrix.
#' @keywords internal
srht_apply <- function(M, plan) {
  .Call('_fmrireg_cpp_srht_apply', PACKAGE = 'fmrireg', M,
        as.integer(plan$rows),
        as.numeric(plan$signs),
        as.integer(plan$perm),
        as.numeric(plan$scale))
}

#' Apply the adjoint of an SRHT plan (S' B) to an m x k matrix
#' @keywords internal
#' @noRd
srht_adjoint <- function(B, plan) {
  .Call('_fmrireg_cpp_srht_adjoint', PACKAGE = 'fmrireg', as.matrix(B),
        as.integer(plan$rows),
        as.numeric(plan$signs),
        as.integer(plan$perm),
        as.numeric(plan$scale))
}

#' Iterative Hessian Sketch (multi-RHS)
#'
#' Returns the coefficients `M`, the exact \eqn{(X'X)^{-1}} as `Ginv`, the
#' exact full-data `residuals`, the number of iterations used and whether the
#' tolerance was met. With `tol = 0` exactly `iters` iterations are run;
#' otherwise iteration stops once every coefficient moves by less than `tol`
#' OLS standard errors.
#' @param X Numeric design matrix (T x p).
#' @param Z Numeric response matrix (T x k).
#' @param m Integer number of sketch rows per Hessian sketch.
#' @param iters Maximum number of iterations (at least 1).
#' @param tol Non-negative stopping tolerance in OLS standard errors.
#' @return A list with `M`, `Ginv`, `residuals`, `iters` and `converged`.
#' @keywords internal
ihs_latent_solve <- function(X, Z, m, iters = 3L, tol = 0) {
  .Call('_fmrireg_cpp_ihs_latent', PACKAGE = 'fmrireg', as.matrix(X), as.matrix(Z),
        as.integer(m), as.integer(iters), as.numeric(tol))
}

# --- Sketch operators and solvers used by the latent_sketch engine ---

#' Internal: build a time-sketch operator
#'
#' Returns `apply(A)` (S A), `adjoint(B)` (S' B) and `gram()` (S S', m x m)
#' for the sketch-and-solve methods; `NULL` for `"ihs"`, which draws its own
#' sketches.
#' @keywords internal
#' @noRd
.lowrank_sketch_operator <- function(Tlen, sk) {
  method <- sk$method
  if (identical(method, "ihs")) return(NULL)
  if (identical(method, "srht")) {
    plan <- make_srht_plan(Tlen, sk$m)
    apply_fn <- function(A) srht_apply(as.matrix(A), plan)
    adjoint_fn <- function(B) srht_adjoint(B, plan)
    gram_fn <- function() {
      K <- apply_fn(adjoint_fn(diag(plan$m)))
      (K + t(K)) / 2
    }
  } else {
    S <- make_time_sketch(Tlen, sk)
    apply_fn <- function(A) as.matrix(S %*% A)
    adjoint_fn <- function(B) as.matrix(Matrix::crossprod(S, B))
    gram_fn <- function() as.matrix(Matrix::tcrossprod(S))
  }
  list(method = method, m = sk$m, apply = apply_fn, adjoint = adjoint_fn,
       gram = gram_fn)
}

#' Internal: inverse of a small SPD matrix with a ridge fallback
#' @keywords internal
#' @noRd
.lowrank_spd_inverse <- function(G) {
  tryCatch(chol2inv(chol(G)), error = function(e) {
    ridge <- 1e-6 * sum(diag(G)) / max(1L, ncol(G))
    chol2inv(chol(G + diag(ridge, ncol(G))))
  })
}

#' Internal: sketch-and-solve with its conditional-on-S inference
#'
#' Solves \eqn{b_s = (X'S'SX)^{-1} X'S'S z} and returns what honest inference
#' for that estimator needs. With \eqn{z = Xb + e}, \eqn{e \sim (0, \sigma^2 I)}
#' and \eqn{X_s = SX}, \eqn{G = X_s'X_s}, \eqn{K = SS'}:
#'
#' * \eqn{b_s - b = G^{-1} X_s' S e}, so
#'   \eqn{Var(b_s | S) = \sigma^2 G^{-1} X_s' K X_s G^{-1}} (returned as
#'   `cov_unscaled`). For an SRHT with \eqn{T} a power of two,
#'   \eqn{K = (T/m) I}, and this reduces to \eqn{(T/m)\sigma^2 G^{-1}}: the
#'   sketched fit discards information, so its standard errors are about
#'   \eqn{\sqrt{T/m}} times the full-data OLS ones.
#' * The sketched residual is \eqn{r_s = P S e} with
#'   \eqn{P = I - X_s G^{-1} X_s'}, so \eqn{E\|r_s\|^2 = \sigma^2 tr(PK)}:
#'   \eqn{\hat\sigma^2 = \|r_s\|^2 / tr(PK)} is unbiased (`kappa = tr(PK)`).
#' * \eqn{\|r_s\|^2} is a quadratic form in \eqn{e}; its Satterthwaite degrees
#'   of freedom are \eqn{tr(PK)^2 / tr(PKPK)} (`df`), which equals \eqn{m - p}
#'   when \eqn{K} is a multiple of the identity.
#'
#' All traces use only m x p and p x p products besides \eqn{K} itself.
#' @keywords internal
#' @noRd
.lowrank_sketch_solve <- function(X, Z, op) {
  Xs <- op$apply(X)
  Zs <- op$apply(Z)
  G <- crossprod(Xs)
  Ginv <- .lowrank_spd_inverse(G)
  M <- Ginv %*% crossprod(Xs, Zs)
  resid <- Zs - Xs %*% M

  K <- op$gram()
  U <- K %*% Xs                        # m x p
  W <- crossprod(Xs, U)                # X_s' K X_s
  GW <- Ginv %*% W
  cov_unscaled <- GW %*% Ginv
  cov_unscaled <- (cov_unscaled + t(cov_unscaled)) / 2

  # Residual traces through an orthonormal basis Q of col(X_s), so H = QQ'
  # respects the design rank even when G needed a ridge:
  # tr(PK) = tr(K) - tr(Q'KQ), tr(PKPK) = tr(KK) - 2||KQ||^2 + ||Q'KQ||^2.
  qx <- qr(Xs)
  Q <- qr.Q(qx)[, seq_len(qx$rank), drop = FALSE]
  KQ <- K %*% Q
  QKQ <- crossprod(Q, KQ)
  kappa <- sum(diag(K)) - sum(diag(QKQ))                    # tr(PK)
  tr_PKPK <- sum(K * K) - 2 * sum(KQ * KQ) + sum(QKQ * QKQ)
  df <- if (kappa > 0 && tr_PKPK > 0) kappa^2 / tr_PKPK else 1
  list(M = M, cov_unscaled = cov_unscaled, residuals = resid,
       kappa = kappa, df = max(df, 1), iters = NA_integer_, converged = NA)
}

#' Internal: IHS solve with exact OLS inference quantities
#' @keywords internal
#' @noRd
.lowrank_ihs_solve <- function(X, Z, sk) {
  iters <- as.integer(sk$iters %||% .lowrank_ihs_default_iters)
  tol <- as.numeric(sk$tol %||% .lowrank_ihs_default_tol)
  sol <- ihs_latent_solve(X, Z, m = sk$m, iters = iters, tol = tol)
  rank_x <- as.integer(Matrix::rankMatrix(X)[1])
  df <- max(1L, nrow(X) - rank_x)
  list(M = sol$M, cov_unscaled = sol$Ginv, residuals = sol$residuals,
       kappa = df, df = df, iters = as.integer(sol$iters),
       converged = if (tol > 0) isTRUE(sol$converged) else NA)
}

.lowrank_ihs_default_iters <- 100L
.lowrank_ihs_default_tol <- 1e-3

#' Internal: dispatch a time-sketched solve
#' @keywords internal
#' @noRd
.lowrank_time_solve <- function(X, Z, sk, op) {
  if (identical(sk$method, "ihs")) {
    .lowrank_ihs_solve(X, Z, sk)
  } else {
    .lowrank_sketch_solve(X, Z, op)
  }
}
