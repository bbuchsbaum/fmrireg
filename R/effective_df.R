#' Calculate Effective Degrees of Freedom
#'
#' Residual degrees of freedom for a fit, adjusted for robust down-weighting and
#' -- when the information is supplied -- for residual correlation that survives
#' prewhitening.
#'
#' @param n Integer. Number of observations.
#' @param p Integer. Number of parameters.
#' @param robust_weights Numeric vector of robust weights (NULL if not robust).
#' @param ar_order Integer. Order of the AR model, or 0. Retained for backward
#'   compatibility and for reporting; see Details for why it does not by itself
#'   reduce the degrees of freedom.
#' @param method Character. `"residual"` (default) uses `n - p`.
#'   `"satterthwaite"` computes the Satterthwaite approximation and requires
#'   both `X` and `resid_cov`. `"simple"` is accepted as an alias for
#'   `"residual"`.
#' @param X Numeric matrix. The design actually fitted (already whitened, if
#'   whitening was applied). Required for `method = "satterthwaite"`.
#' @param resid_cov Numeric matrix. Correlation (or covariance) of the errors of
#'   the model as fitted -- that is, *after* any whitening. Required for
#'   `method = "satterthwaite"`. Supply the identity, or omit, when whitening is
#'   believed exact.
#' @param XtXinv Optional precomputed `(X'X)^-1` for `method = "satterthwaite"`
#'   (lets callers hoist it out of a voxel loop).
#'
#' @return Numeric. Effective residual degrees of freedom, at least 1.
#'
#' @details
#' Write \eqn{R = I - X(X'X)^{-1}X'} for the residual-forming matrix of the model
#' as fitted, and \eqn{V} for the correlation of that model's errors. Then
#' \eqn{e'e} is a quadratic form whose Satterthwaite degrees of freedom are
#' \deqn{\nu = \mathrm{tr}(RV)^2 / \mathrm{tr}(RVRV).}
#'
#' When whitening succeeds, \eqn{V = I} and this reduces exactly to
#' \eqn{\mathrm{rank}(R) = n - p}. That is why `ar_order` alone does not reduce
#' the degrees of freedom: fmrireg's AR path *whitens* the data and then fits by
#' ordinary least squares, so if the filter is adequate the whitened residuals
#' already carry \eqn{n - p} degrees of freedom, and reducing them again would
#' double-count a correction that has already been applied.
#'
#' Any shortfall comes from the filter being imperfect, not from the AR order,
#' and it is captured by passing the post-whitening `resid_cov` and using
#' `method = "satterthwaite"`.
#'
#' Robust down-weighting is handled by [.robust_effective_df()], the single
#' implementation used by every fitting path: an observation with weight
#' \eqn{w_i} counts as \eqn{w_i} of an observation, so the residual df is
#' \eqn{\sum_i w_i - p}. Zero weights contribute nothing, and 0/1 weights give
#' exactly the residual df of `stats::lm(weights = w)`. With
#' `method = "residual"` the weights alone fix the effective sample (`n` is
#' not used); with `method = "satterthwaite"` the Satterthwaite df is scaled by
#' \eqn{(\sum_i w_i - p) / (n_w - p)}, with \eqn{n_w} the number of weights.
#'
#' @section History:
#' Before 2026-08-07 this function multiplied the residual degrees of freedom by
#' \eqn{1 / (1 + 2\sum_k (1 - k/n))} whenever `ar_order > 0`. That factor is the
#' variance inflation one would obtain if *every* autocorrelation equalled 1, it
#' never depended on the fitted AR coefficients, and it was applied on top of
#' whitening. At n = 200, p = 12, AR(1) it returned 62.9 regardless of the true
#' autocorrelation, where a coefficient-aware calculation spans roughly 170 to 10
#' as phi runs from 0.05 to 0.9. The result was conservative -- it cost power
#' rather than inflating false positives -- but it had no basis, and
#' `method = "satterthwaite"` was a verbatim duplicate of `method = "simple"`.
#'
#' @keywords internal
#' @noRd
calculate_effective_df <- function(n, p, robust_weights = NULL,
                                   ar_order = 0,
                                   method = c("residual", "satterthwaite", "simple"),
                                   X = NULL, resid_cov = NULL, XtXinv = NULL) {
  method <- match.arg(method)
  if (identical(method, "simple")) method <- "residual"

  n <- as.numeric(n)
  p <- as.numeric(p)
  if (!is.finite(n) || !is.finite(p) || n <= p) {
    return(1)
  }

  if (identical(method, "satterthwaite")) {
    if (is.null(X) || is.null(resid_cov)) {
      stop("method = 'satterthwaite' requires both `X` and `resid_cov`.",
           call. = FALSE)
    }
    df_base <- .satterthwaite_df(X, resid_cov, XtXinv = XtXinv)
    if (is.null(robust_weights)) return(max(df_base, 1))
    return(.robust_effective_df(robust_weights, p, df_base = df_base))
  }

  if (is.null(robust_weights)) return(max(n - p, 1))
  .robust_effective_df(robust_weights, p)
}

#' Effective residual df under robust (observation) weights
#'
#' The one robust effective-df formula used by every fitting path:
#' \deqn{\nu = \sum_i w_i - r,}
#' where \eqn{r} is the rank of the fitted design. A weight-\eqn{w} observation
#' counts as \eqn{w} of an observation; zero weights contribute nothing, so 0/1
#' weights reproduce the residual df of `stats::lm(weights = w)`.
#'
#' When another df (`df_base`, e.g. a Satterthwaite or HAC df computed for the
#' unweighted residual count) is to be adjusted, it is scaled by the same
#' ratio, \eqn{\nu_{base} (\sum_i w_i - r) / (n_w - r)}, which reduces to the
#' formula above when \eqn{\nu_{base} = n_w - r}.
#'
#' @param weights Observation weights (non-finite entries are ignored).
#' @param rank Rank of the fitted design.
#' @param df_base Optional df to scale; `NULL` returns \eqn{\sum w - r}.
#' @return Effective df, at least 1.
#' @keywords internal
#' @noRd
.robust_effective_df <- function(weights, rank, df_base = NULL) {
  w <- as.numeric(weights)
  w <- w[is.finite(w)]
  rank <- as.numeric(rank)
  if (!length(w) || !is.finite(rank)) return(1)
  df_w <- sum(w) - rank
  if (is.null(df_base)) return(max(df_w, 1))
  denom <- length(w) - rank
  if (denom <= 0) return(1)
  max(df_base * df_w / denom, 1)
}

#' Satterthwaite degrees of freedom for a residual sum of squares
#'
#' `tr(RV)^2 / tr(RVRV)` with `R = I - X(X'X)^-1 X'`. Computed without forming
#' `R` explicitly: with `H = X(X'X)^-1X'`, `tr(RV) = tr(V) - tr(HV)` and
#' `tr(RVRV) = tr(VV) - 2 tr(HVV) + tr(HVHV)`, each of which needs only `n x p`
#' and `p x p` products.
#'
#' @keywords internal
#' @noRd
.satterthwaite_xtxinv <- function(X) {
  XtX <- crossprod(X)
  tryCatch(solve(XtX), error = function(e) {
    sv <- svd(XtX)
    keep <- sv$d > max(sv$d) * .Machine$double.eps^0.5
    sv$v[, keep, drop = FALSE] %*% ((1 / sv$d[keep]) * t(sv$u[, keep, drop = FALSE]))
  })
}

.satterthwaite_df <- function(X, V, XtXinv = NULL) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (is.null(dim(V))) V <- diag(as.numeric(V), n)
  V <- as.matrix(V)
  stopifnot(nrow(V) == n, ncol(V) == n)

  if (is.null(XtXinv)) XtXinv <- .satterthwaite_xtxinv(X)

  A  <- X %*% XtXinv          # n x p; H = A X'
  VX <- V %*% X               # n x p
  B  <- crossprod(X, VX)      # p x p; X' V X

  # tr(RV)   = tr(V) - tr(HV)
  # tr(RVRV) = tr(VV) - 2 tr(HVV) + tr(HVHV),  using tr(VHV) = tr(HVV)
  tr_RV   <- sum(diag(V)) - sum(A * VX)
  tr_HVV  <- sum(A * (V %*% VX))
  BG      <- B %*% XtXinv
  tr_RVRV <- sum(V * t(V)) - 2 * tr_HVV + sum(BG * t(BG))

  if (!is.finite(tr_RV) || !is.finite(tr_RVRV) || tr_RVRV <= 0) return(1)
  max(tr_RV^2 / tr_RVRV, 1)
}

