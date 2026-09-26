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

#' Internal: estimability of a (whitened) design for the sketch engine
#'
#' Judges estimability on the full design with the same tolerance-aware
#' pivoted QR (`.design_rank_info()`, default tolerance) that the exact OLS
#' fast path uses in `.fast_preproject()`, so sketched and exact fits alias
#' the same columns. A sketch is never used to decide estimability: with
#' `m` close to `p` it can create spurious rank loss, and it cannot remove
#' a genuine one.
#' @keywords internal
#' @noRd
.lowrank_rank_info <- function(X, warn = FALSE) {
  info <- .design_rank_info(X)
  if (info$rank == 0L) {
    stop("latent_sketch engine: the design matrix has no estimable column.",
         call. = FALSE)
  }
  if (warn && length(info$aliased) > 0L) {
    warning(.rank_deficiency_message(info$rank, ncol(X), info$aliased,
                                     colnames(X)), call. = FALSE)
  }
  info
}

#' Internal: embed an estimable-columns solve in the full coefficient space
#'
#' Keep-but-aliased: every declared column keeps its position. Aliased
#' coefficients and their covariance rows/columns are held at zero internally
#' and marked through the rank attributes on `cov_unscaled`, exactly as
#' `.fast_preproject()` marks `XtXinv` for exact fits; the statistics and
#' contrast code report them as NA (contrasts with a named warning).
#' @keywords internal
#' @noRd
.lowrank_expand_solution <- function(sol, info, est, varnames, p) {
  M <- matrix(0, p, ncol(sol$M))
  M[est, ] <- sol$M
  cov_unscaled <- matrix(0, p, p)
  cov_unscaled[est, est] <- sol$cov_unscaled
  if (!is.null(varnames)) {
    rownames(M) <- varnames
    dimnames(cov_unscaled) <- list(varnames, varnames)
  }
  sol$M <- M
  sol$cov_unscaled <- .attach_rank_attrs(cov_unscaled, info)
  sol$rank <- info$rank
  sol$aliased <- info$aliased
  sol
}

#' Internal: sketch-and-solve with its conditional-on-S inference
#'
#' Solves \eqn{b_s = (X'S'SX)^{-1} X'S'S z} and returns what honest inference
#' for that estimator needs. `X` must have full column rank: the caller,
#' `.lowrank_time_solve()`, passes only the estimable columns. With
#' \eqn{z = Xb + e}, \eqn{e \sim (0, \sigma^2 I)}
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
#' \eqn{G} is inverted through the pivoted QR of \eqn{X_s}, never through a
#' Cholesky factor of \eqn{G}, which can succeed on a numerically singular
#' \eqn{G} and return garbage. If \eqn{X_s} is numerically rank deficient
#' under the QR tolerance used for the full design, the sketch has too few
#' rows for this design and the solve stops. There is no ridge or
#' pseudo-inverse fallback: either would silently change the estimator that
#' the reported inference describes.
#'
#' All traces use only m x p and p x p products besides \eqn{K} itself.
#' @keywords internal
#' @noRd
.lowrank_sketch_solve <- function(X, Z, op, tol = 1e-7) {
  Xs <- op$apply(X)
  Zs <- op$apply(Z)
  k <- ncol(Xs)
  qx <- qr(Xs, tol = tol, LAPACK = FALSE)
  if (qx$rank < k) {
    stop(sprintf(paste0(
      "time_sketch$m = %d is too small for this design: the sketched design ",
      "(%d rows x %d estimable columns) has numerical rank %d. Increase ",
      "`time_sketch$m`%s."), as.integer(op$m), nrow(Xs), k, qx$rank,
      if (identical(op$method, "countsketch")) {
        paste0(", or use method \"srht\" or \"gaussian\": CountSketch maps ",
               "each scan to one row, so sparse regressors such as spike ",
               "(scrubbing) columns collide")
      } else ""),
      call. = FALSE)
  }
  Rinv <- backsolve(qr.R(qx), diag(k))
  Ginv <- matrix(0, k, k)
  Ginv[qx$pivot, qx$pivot] <- tcrossprod(Rinv)
  M <- qr.coef(qx, Zs)
  if (is.null(dim(M))) M <- matrix(M, nrow = k)
  resid <- qr.resid(qx, Zs)
  if (is.null(dim(resid))) resid <- matrix(resid, ncol = ncol(M))

  K <- op$gram()
  U <- K %*% Xs                        # m x p
  W <- crossprod(Xs, U)                # X_s' K X_s
  GW <- Ginv %*% W
  cov_unscaled <- GW %*% Ginv
  cov_unscaled <- (cov_unscaled + t(cov_unscaled)) / 2

  # Residual traces through the orthonormal basis Q of col(X_s) (full column
  # rank, checked above), so H = QQ' has the estimable rank:
  # tr(PK) = tr(K) - tr(Q'KQ), tr(PKPK) = tr(KK) - 2||KQ||^2 + ||Q'KQ||^2.
  Q <- qr.Q(qx)
  KQ <- K %*% Q
  QKQ <- crossprod(Q, KQ)
  kappa <- sum(diag(K)) - sum(diag(QKQ))                    # tr(PK)
  tr_PKPK <- sum(K * K) - 2 * sum(KQ * KQ) + sum(QKQ * QKQ)
  df <- if (kappa > 0 && tr_PKPK > 0) kappa^2 / tr_PKPK else 1
  list(M = M, cov_unscaled = cov_unscaled, residuals = resid,
       kappa = kappa, df = max(df, 1), iters = NA_integer_, converged = NA)
}

#' Internal: IHS solve with exact OLS inference quantities
#'
#' `X` must have full column rank (see `.lowrank_time_solve()`), so the
#' residual df is `nrow(X) - ncol(X)`. Each IHS step inverts an m x p
#' sketched Hessian, which is singular when `m < p`; the kernel no longer
#' falls back to a pseudo-inverse there (that returned coefficients hundreds
#' of OLS standard errors off with no warning), so both cases stop with an
#' error naming `m`. The reported covariance is `(X'X)^{-1}` from the QR of
#' `X`, the same factorisation exact fits use.
#' @keywords internal
#' @noRd
.lowrank_ihs_solve <- function(X, Z, sk, tol_qr = 1e-7) {
  k <- ncol(X)
  too_small <- function(detail) {
    stop(sprintf(paste0(
      "time_sketch$m = %d is too small for this design: %s. Increase ",
      "`time_sketch$m`."), as.integer(sk$m), detail), call. = FALSE)
  }
  if (sk$m < k) {
    too_small(sprintf("IHS needs at least as many sketch rows as the %d estimable columns", k))
  }
  iters <- as.integer(sk$iters %||% .lowrank_ihs_default_iters)
  tol <- as.numeric(sk$tol %||% .lowrank_ihs_default_tol)
  sol <- tryCatch(
    ihs_latent_solve(X, Z, m = sk$m, iters = iters, tol = tol),
    error = function(e) {
      if (grepl("sketched Gram", conditionMessage(e), fixed = TRUE)) {
        too_small("a sketched Hessian of the estimable design is singular")
      }
      stop(e)
    }
  )
  qx <- qr(X, tol = tol_qr, LAPACK = FALSE)
  Rinv <- backsolve(qr.R(qx), diag(k))
  XtXinv <- matrix(0, k, k)
  XtXinv[qx$pivot, qx$pivot] <- tcrossprod(Rinv)
  df <- max(1L, nrow(X) - k)
  list(M = sol$M, cov_unscaled = XtXinv, residuals = sol$residuals,
       kappa = df, df = df, iters = as.integer(sol$iters),
       converged = if (tol > 0) isTRUE(sol$converged) else NA)
}

.lowrank_ihs_default_iters <- 100L
.lowrank_ihs_default_tol <- 1e-3

#' Internal: dispatch a time-sketched solve on the estimable columns
#'
#' Estimability is judged on the full (whitened) `X`; the solver sees only
#' the estimable columns, and its solution is embedded back into the declared
#' coefficient space by `.lowrank_expand_solution()`.
#' @keywords internal
#' @noRd
.lowrank_time_solve <- function(X, Z, sk, op, warn = FALSE) {
  X <- as.matrix(X)
  info <- .lowrank_rank_info(X, warn = warn)
  est <- sort(info$estimable)
  Xe <- X[, est, drop = FALSE]
  sol <- if (identical(sk$method, "ihs")) {
    .lowrank_ihs_solve(Xe, Z, sk, tol_qr = info$tol)
  } else {
    .lowrank_sketch_solve(Xe, Z, op, tol = info$tol)
  }
  .lowrank_expand_solution(sol, info, est, colnames(X), ncol(X))
}
