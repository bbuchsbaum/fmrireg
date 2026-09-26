#' Build a time sketch matrix S (m x T)
#'
#' All sketches are normalised so that \eqn{E\|S r\|^2 = \|r\|^2}.
#'
#' @param Tlen Integer time length
#' @param ctrl List(method = "gaussian"|"countsketch"|"srht", m)
#' @return A dense or sparse sketch matrix for Gaussian/CountSketch methods;
#'   `NULL` for `"srht"`, which is applied through a plan ([srht_apply()]).
#' @keywords internal
make_time_sketch <- function(Tlen, ctrl) {
  stopifnot(is.list(ctrl))
  method <- match.arg(ctrl$method %||% "gaussian", c("gaussian", "countsketch", "srht"))
  m <- ctrl$m
  stopifnot(!is.null(m), m > 0L, m <= Tlen)
  if (method == "gaussian") {
    S <- matrix(stats::rnorm(m * Tlen), nrow = m, ncol = Tlen) / sqrt(m)
  } else if (method == "countsketch") {
    rows <- sample.int(m, Tlen, replace = TRUE)
    signs <- sample(c(-1, 1), Tlen, replace = TRUE)
    S <- Matrix::sparseMatrix(i = rows, j = seq_len(Tlen), x = signs, dims = c(m, Tlen))
  } else {
    # SRHT does not materialise S; callers use srht_apply()/srht_adjoint().
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
  .Call('_fmrireg_cpp_srht_apply', PACKAGE = 'fmrireg', as.matrix(M),
        as.integer(plan$rows),
        as.numeric(plan$signs),
        as.integer(plan$perm),
        as.numeric(plan$scale),
        .sketch_num_threads())
}

#' Apply the adjoint of an SRHT plan (S' B) to an m x k matrix
#' @keywords internal
#' @noRd
srht_adjoint <- function(B, plan) {
  .Call('_fmrireg_cpp_srht_adjoint', PACKAGE = 'fmrireg', as.matrix(B),
        as.integer(plan$rows),
        as.numeric(plan$signs),
        as.integer(plan$perm),
        as.numeric(plan$scale),
        .sketch_num_threads())
}

#' Internal: thread count for the column-parallel SRHT kernels
#'
#' Follows the other compiled kernels: `getOption("fmrireg.num_threads", 0)`,
#' where 0 (or an unusable value) means the OpenMP default.
#' @keywords internal
#' @noRd
.sketch_num_threads <- function() {
  n <- suppressWarnings(as.integer(getOption("fmrireg.num_threads", 0L))[1L])
  if (length(n) != 1L || is.na(n) || n < 0L) 0L else n
}

# --- Sketch operators and solvers used by the latent_sketch engine ---

#' Internal: build a time-sketch operator
#'
#' Returns `apply(A)` (S A), `adjoint(B)` (S' B) and `gram()` (S S', m x m)
#' for the sketch-and-solve methods; `NULL` for the exact solve (`"ihs"`,
#' deprecated). `gram()` is computed once and memoised: the sketch is fixed
#' for the operator's lifetime, and by_cluster fits call it once per cluster.
#' @keywords internal
#' @noRd
.lowrank_sketch_operator <- function(Tlen, sk) {
  method <- sk$method
  if (identical(method, "ihs")) return(NULL)
  if (identical(method, "srht")) {
    plan <- make_srht_plan(Tlen, sk$m)
    apply_fn <- function(A) srht_apply(as.matrix(A), plan)
    adjoint_fn <- function(B) srht_adjoint(B, plan)
    compute_gram <- function() {
      K <- apply_fn(adjoint_fn(diag(plan$m)))
      (K + t(K)) / 2
    }
  } else {
    S <- make_time_sketch(Tlen, sk)
    apply_fn <- function(A) as.matrix(S %*% A)
    adjoint_fn <- function(B) as.matrix(Matrix::crossprod(S, B))
    compute_gram <- function() as.matrix(Matrix::tcrossprod(S))
  }
  K_cache <- NULL
  gram_fn <- function() {
    if (is.null(K_cache)) K_cache <<- compute_gram()
    K_cache
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
  # Coefficients and residuals through the explicit orthonormal factor Q
  # (m x k): two BLAS products over the response columns. qr.coef() and
  # qr.resid() apply the Householder reflections one response column at a
  # time in LINPACK and dominated the solve for r ~ 1e4-1e5 columns.
  Q <- qr.Q(qx)
  QtZ <- crossprod(Q, Zs)              # k x r
  M <- matrix(0, k, ncol(Zs))
  M[qx$pivot, ] <- Rinv %*% QtZ
  resid <- Zs - Q %*% QtZ

  K <- op$gram()
  U <- K %*% Xs                        # m x p
  W <- crossprod(Xs, U)                # X_s' K X_s
  GW <- Ginv %*% W
  cov_unscaled <- GW %*% Ginv
  cov_unscaled <- (cov_unscaled + t(cov_unscaled)) / 2

  # Residual traces through the orthonormal basis Q of col(X_s) (full column
  # rank, checked above), so H = QQ' has the estimable rank:
  # tr(PK) = tr(K) - tr(Q'KQ), tr(PKPK) = tr(KK) - 2||KQ||^2 + ||Q'KQ||^2.
  KQ <- K %*% Q
  QKQ <- crossprod(Q, KQ)
  kappa <- sum(diag(K)) - sum(diag(QKQ))                    # tr(PK)
  tr_PKPK <- sum(K * K) - 2 * sum(KQ * KQ) + sum(QKQ * QKQ)
  df <- if (kappa > 0 && tr_PKPK > 0) kappa^2 / tr_PKPK else 1
  list(M = M, cov_unscaled = cov_unscaled, residuals = resid,
       kappa = kappa, df = max(df, 1))
}

#' Internal: exact OLS solve for the latent_sketch engine
#'
#' Used for the deprecated `method = "ihs"`. In this multi-response regime
#' (p small, many response columns) the cost of least squares is the single
#' pass X'Z, and an iterative Hessian sketch has to repeat that pass for
#' every exact gradient, so it can never beat one-step OLS. The solve goes
#' through `.fast_preproject()`, the projection exact `fmri_lm()` fits use,
#' so coefficients, the keep-but-aliased rank handling, `(X'X)^{-1}` and the
#' `T - rank` residual degrees of freedom are those of exact OLS. The rank
#' warning is raised (or not) by the caller through `.lowrank_rank_info()`,
#' so the identical warning from `.fast_preproject()` is muffled here.
#'
#' Coefficients and residual sums of squares come from `solve_glm_core()`,
#' the solver of the exact joint (chunkwise) OLS path, so voxel-space fits
#' reproduce that path exactly. `rss` is what `.lowrank_sigma2()` uses for
#' voxel data; the residuals are kept for latent datasets, whose voxel
#' residual variances need the latent residual cross-products.
#' @keywords internal
#' @noRd
.lowrank_exact_solve <- function(X, Z, info) {
  rank_msg <- if (length(info$aliased)) {
    .rank_deficiency_message(info$rank, ncol(X), info$aliased, colnames(X))
  }
  proj <- withCallingHandlers(
    .fast_preproject(X, tol = info$tol),
    warning = function(w) {
      if (identical(conditionMessage(w), rank_msg)) invokeRestart("muffleWarning")
    }
  )
  res <- solve_glm_core(glm_context(X = X, Y = Z, proj = proj))
  M <- res$betas
  attributes(M) <- list(dim = dim(M), dimnames = dimnames(M))
  df <- as.numeric(proj$dfres)
  # A saturated design has no residual df: the variance is NA, as in the
  # exact path, rather than rss / 0.
  list(M = M, cov_unscaled = proj$XtXinv, residuals = Z - X %*% M,
       rss = res$rss, kappa = if (df > 0) df else NA_real_, df = df)
}

#' Internal: dispatch a time-sketched solve on the estimable columns
#'
#' Estimability is judged on the full (whitened) `X`; the solver sees only
#' the estimable columns, and its solution is embedded back into the declared
#' coefficient space by `.lowrank_expand_solution()`. The exact solve
#' (`"ihs"`) works on the full design directly, like exact fits.
#' @keywords internal
#' @noRd
.lowrank_time_solve <- function(X, Z, sk, op, warn = FALSE) {
  X <- as.matrix(X)
  info <- .lowrank_rank_info(X, warn = warn)
  if (identical(sk$method, "ihs")) {
    return(.lowrank_exact_solve(X, as.matrix(Z), info))
  }
  est <- sort(info$estimable)
  Xe <- X[, est, drop = FALSE]
  sol <- .lowrank_sketch_solve(Xe, Z, op, tol = info$tol)
  .lowrank_expand_solution(sol, info, est, colnames(X), ncol(X))
}
