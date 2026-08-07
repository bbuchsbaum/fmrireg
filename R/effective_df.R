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
#' Robust down-weighting is handled separately and multiplicatively, by
#' \eqn{\sum_i w_i / n}.
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
                                   X = NULL, resid_cov = NULL) {
  method <- match.arg(method)
  if (identical(method, "simple")) method <- "residual"

  n <- as.numeric(n)
  p <- as.numeric(p)
  if (!is.finite(n) || !is.finite(p) || n <= p) {
    return(1)
  }

  df_base <- if (identical(method, "satterthwaite")) {
    if (is.null(X) || is.null(resid_cov)) {
      stop("method = 'satterthwaite' requires both `X` and `resid_cov`.",
           call. = FALSE)
    }
    .satterthwaite_df(X, resid_cov)
  } else {
    n - p
  }

  df_robust_adjust <- 1
  if (!is.null(robust_weights)) {
    w <- robust_weights[is.finite(robust_weights)]
    if (length(w)) {
      # Down-weighting every observation to zero leaves no effective data, so the
      # adjustment is 0 and the floor below takes over.
      df_robust_adjust <- sum(w) / length(w)
    }
  }

  max(df_base * df_robust_adjust, 1)
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
.satterthwaite_df <- function(X, V) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (is.null(dim(V))) V <- diag(as.numeric(V), n)
  V <- as.matrix(V)
  stopifnot(nrow(V) == n, ncol(V) == n)

  XtX <- crossprod(X)
  XtXinv <- tryCatch(solve(XtX), error = function(e) {
    sv <- svd(XtX)
    keep <- sv$d > max(sv$d) * .Machine$double.eps^0.5
    sv$v[, keep, drop = FALSE] %*% ((1 / sv$d[keep]) * t(sv$u[, keep, drop = FALSE]))
  })

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

calculate_sandwich_variance <- function(X, residuals, XtXinv, weights = NULL) {
  n <- nrow(X)
  p <- ncol(X)
  
  # Squared residuals (meat of sandwich)
  resid2 <- residuals^2
  
  if (!is.null(weights)) {
    # Weight the squared residuals
    resid2 <- resid2 * weights
  }
  
  # Calculate meat: X' diag(resid2) X
  # Correct computation without sqrt (resid2 is already squared)
  if (is.matrix(residuals) && ncol(residuals) > 1) {
    # Multiple response case
    # Average over voxels for now (could be improved)
    resid2_avg <- rowMeans(resid2)
    meat <- crossprod(X, X * resid2_avg)  # X' diag(resid2) X
  } else {
    meat <- crossprod(X, X * as.vector(resid2))  # X' diag(resid2) X
  }
  
  # Sandwich: (X'X)^-1 meat (X'X)^-1
  sandwich <- XtXinv %*% meat %*% XtXinv
  
  # Small sample correction
  correction <- n / (n - p)
  
  sandwich * correction
}

#' Adjust Statistics for Effective Degrees of Freedom
#'
#' Updates t-statistics and p-values using effective degrees of freedom
#' rather than nominal degrees of freedom.
#'
#' @param t_stats Vector of t-statistics.
#' @param df_nominal Nominal degrees of freedom.
#' @param df_effective Effective degrees of freedom.
#'
#' @return List with adjusted t-statistics and p-values.
#'
#' @keywords internal
#' @noRd
adjust_stats_for_effective_df <- function(t_stats, df_nominal, df_effective) {
  # Adjustment factor for t-statistics
  # Based on the relationship between t-distributions with different df
  adjust_factor <- sqrt(df_effective / df_nominal)
  
  # Adjusted t-statistics
  t_adj <- t_stats * adjust_factor
  
  # P-values using effective df
  p_values <- 2 * pt(abs(t_adj), df = df_effective, lower.tail = FALSE)
  
  list(
    t_stats = t_adj,
    p_values = p_values,
    df = df_effective
  )
}
