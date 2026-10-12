#' Soft Subspace Projection for Nuisance Removal
#'
#' Ridge-regularized projection to remove nuisance variance without
#' explicitly choosing components or adding confound regressors.
#' This is "CompCor without choosing components."
#'
#' @section Conceptual Overview:
#' Traditional CompCor extracts K principal components from white matter/CSF
#' and regresses them out. This requires choosing K, which is arbitrary.
#'
#' Soft subspace projection instead treats the entire WM/CSF timeseries as a
#' nuisance basis and removes it with shrinkage. Each nuisance direction is
#' partially removed proportional to its variance, avoiding the hard keep/delete
#' decision of PCA truncation.
#'
#' @section The Projection Operator:
#' Given nuisance matrix N (time x nuisance_voxels), the soft projection is:
#' \deqn{P_\lambda = I - N (N^T N + \lambda I)^{-1} N^T}
#'
#' Applied to data Y and design X:
#' \itemize{
#'   \item \code{Y_clean = P_lambda \%*\% Y}
#'
#'   \item \code{X_clean = P_lambda \%*\% X} (important to avoid bias)
#' }
#'
#' @section Lambda Selection:
#' The shrinkage parameter lambda controls aggressiveness:
#' \itemize{
#'   \item Small lambda: aggressive removal (risk removing signal)
#'   \item Large lambda: gentle removal (risk leaving nuisance)
#' }
#'
#' Three selection methods are available:
#' \describe{
#'   \item{"auto"}{Singular value heuristic: \code{lambda = median(d^2)} where d
#'     are singular values of N. Fast, stable, no tuning. Components with variance
#'     below median are heavily shrunk; above median lightly shrunk.}
#'   \item{"gcv"}{Generalized cross-validation minimizes \code{RSS/(1-df/n)^2}
#'     for ridge regression of Y on N. Finds lambda giving best leave-one-out
#'     prediction without actually doing LOO. Requires Y.
#'     Computational cost is O(k) per evaluation where k = min(n, p).}
#'   \item{numeric}{User-specified lambda value.}
#' }
#'
#' @section Typical Usage:
#' The soft subspace projection workflow has three steps:
#' \enumerate{
#'   \item Extract nuisance timeseries from WM/CSF mask (or provide pre-computed)
#'   \item Create the soft projection operator
#'   \item Apply to both data Y and design matrix X before GLM fitting
#' }
#'
#' For use within \code{fmri_lm()}, see the convenience parameter
#' \code{nuisance_projection} or the more detailed \code{soft_subspace_options}.
#'
#' @section When to Use:
#' Soft subspace projection is most beneficial when:
#' \itemize{
#'   \item You have physiological noise from WM/CSF that isn't captured by motion parameters

#'   \item Traditional CompCor requires arbitrary component selection
#'   \item You want to avoid adding many confound regressors to the design
#' }
#'
#' Consider alternatives when:
#' \itemize{
#'   \item Motion parameters alone sufficiently control artifacts
#'   \item You prefer explicit confound regressors for interpretability
#'   \item Data has very few timepoints relative to nuisance dimensions
#' }
#'
#' @name soft_subspace_projection
NULL

#' Create Soft Projection Operator
#'
#' Computes the soft (ridge-regularized) projection matrix that removes
#' variance explainable by the nuisance subspace while avoiding overfitting.
#'
#' @param N Nuisance matrix (time x nuisance_voxels). Typically extracted
#'   from white matter and CSF voxels. Can have many columns (thousands);
#'   computational cost depends on min(nrow, ncol) due to SVD.
#' @param lambda Ridge penalty controlling shrinkage strength:
#'   \describe{
#'     \item{"auto"}{(Default) Uses median of squared singular values.
#'       Fast, stable, requires no tuning.
#'       Recommended for most use cases.}
#'     \item{"gcv"}{Generalized cross-validation. Optimizes prediction of Y
#'       from N with ridge penalty. Requires Y parameter. More principled but
#'       slower.}
#'     \item{numeric}{User-specified value. Larger = less aggressive removal.}
#'   }
#' @param Y Optional data matrix for GCV-based lambda selection. Required if
#'   \code{lambda = "gcv"}, ignored otherwise.
#' @return A list with class \code{"soft_projection"} containing:
#'   \item{P_lambda}{Function that applies the projection to a matrix}
#'   \item{lambda}{The selected/specified lambda value}
#'   \item{method}{Method used: "singular_value_heuristic", "gcv", or "user_specified"}
#'   \item{effective_df}{Effective degrees of freedom removed (sum of shrinkage factors)}
#'   \item{n_nuisance}{Number of nuisance columns in N}
#'   \item{n_timepoints}{Number of timepoints}
#' @export
#' @examples
#' # Create nuisance matrix (e.g., from WM/CSF voxels)
#' set.seed(123)
#' N <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
#' proj <- soft_projection(N, lambda = "auto")
#'
#' # Apply to data and design
#' Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
#' Y_clean <- proj$P_lambda(Y)
#'
#' # Full workflow: project both data and design
#' X <- cbind(1, rnorm(100), rnorm(100))  # intercept + 2 predictors
#' cleaned <- apply_soft_projection(proj, Y, X)
#' # Now fit GLM with cleaned$Y and cleaned$X
#'
#' # Using GCV for lambda selection (data-driven)
#' proj_gcv <- soft_projection(N, lambda = "gcv", Y = Y)
#' print(proj_gcv)  # Shows selected lambda and effective df
soft_projection <- function(N, lambda = "auto", Y = NULL) {
  if (!is.matrix(N)) {
    N <- as.matrix(N)
  }

  n <- nrow(N)
  p <- ncol(N)

  # Compute SVD of N for efficient operations (only U and d are used)
  svd_N <- svd(N, nu = min(n, p), nv = 0)
  d <- svd_N$d
  U <- svd_N$u

  # Select lambda
  if (is.character(lambda)) {
    lambda <- match.arg(lambda, c("auto", "gcv"))

    if (lambda == "auto") {
      # Singular value heuristic: median of squared singular values
      lambda_val <- median(d^2)
      method <- "singular_value_heuristic"
    } else if (lambda == "gcv") {
      if (is.null(Y)) {
        warning("GCV requires Y; falling back to auto")
        lambda_val <- median(d^2)
        method <- "singular_value_heuristic"
      } else {
        lambda_val <- .select_lambda_gcv(N, Y, svd_N)
        method <- "gcv"
      }
    }
  } else {
    lambda_val <- as.numeric(lambda)
    method <- "user_specified"
    if (lambda_val < 0) {
      stop("lambda must be non-negative")
    }
  }

  # Compute effective degrees of freedom removed
  # df = sum(d^2 / (d^2 + lambda))
  effective_df <- sum(d^2 / (d^2 + lambda_val))

  # Create projection function using SVD for efficiency
  # P_lambda = I - U * diag(d^2/(d^2+lambda)) * U^T
  shrink_factors <- d^2 / (d^2 + lambda_val)

  P_lambda <- function(X) {
    if (!is.matrix(X)) {
      X <- as.matrix(X)
    }
    # X - U %*% diag(shrink_factors) %*% t(U) %*% X
    UX <- crossprod(U, X)  # k x cols
    X - U %*% (shrink_factors * UX)
  }

  structure(
    list(
      P_lambda = P_lambda,
      lambda = lambda_val,
      method = method,
      effective_df = effective_df,
      svd = svd_N,
      n_nuisance = p,
      n_timepoints = n
    ),
    class = "soft_projection"
  )
}

#' @export
print.soft_projection <- function(x, ...) {
  cat("<soft_projection>\n")
  cat("  Lambda:", format(x$lambda, digits = 4),
      "(", x$method, ")\n", sep = "")
  cat("  Nuisance dimensions:", x$n_nuisance, "\n")
  cat("  Effective df removed:", format(x$effective_df, digits = 2), "\n")
  invisible(x)
}

#' Select Lambda via Generalized Cross-Validation
#'
#' Uses GCV to select the ridge penalty that minimizes prediction error
#' of Y from N.
#'
#' @param N Nuisance matrix.
#' @param Y Data matrix.
#' @param svd_N Pre-computed SVD of N.
#' @return Optimal lambda value.
#' @keywords internal
#' @noRd
.select_lambda_gcv <- function(N, Y, svd_N) {
  d <- svd_N$d
  U <- svd_N$u
  n <- nrow(N)

  # GCV(lambda) = ||Y - Y_hat||^2 / (1 - df/n)^2, df = sum(d^2/(d^2+lambda)),
  # Y_hat = U diag(s) U'Y. Y enters only through ||Y||^2 and ||u_k'Y||^2,
  # computed once here; each evaluation is then O(k).
  summ <- .gcv_summaries(U, Y)

  # Search over reasonable range of lambda
  # From very small (near full projection) to very large (minimal projection)
  log_lambda_range <- log(c(min(d^2) / 100, max(d^2) * 100))

  opt <- stats::optimize(.gcv_score_closed_form, interval = log_lambda_range,
                         d = d, summ = summ, n = n)
  exp(opt$minimum)
}

#' GCV sufficient statistics of Y for the soft projection
#'
#' @param U Orthonormal left singular vectors of N (n x k).
#' @param Y Data matrix (n x V).
#' @return list(yy = ||Y||_F^2, uy2 = ||u_k'Y||^2 for each k).
#' @keywords internal
#' @noRd
.gcv_summaries <- function(U, Y) {
  Y <- as.matrix(Y)
  UY <- crossprod(U, Y)  # k x nvox
  list(yy = sum(Y^2), uy2 = rowSums(UY^2))
}

#' Closed-form GCV score at log(lambda)
#'
#' RSS(lambda) = ||Y||^2 - sum_k (2 s_k - s_k^2) ||u_k'Y||^2 with
#' s_k = d_k^2 / (d_k^2 + lambda), since U has orthonormal columns.
#' @keywords internal
#' @noRd
.gcv_score_closed_form <- function(log_lambda, d, summ, n) {
  lam <- exp(log_lambda)
  shrink <- d^2 / (d^2 + lam)
  df <- sum(shrink)
  rss <- summ$yy - sum((2 * shrink - shrink^2) * summ$uy2)

  denom <- (1 - df / n)^2
  if (denom < 1e-10) return(Inf)

  rss / denom
}

#' Apply Soft Projection to Data and Design
#'
#' Applies the soft projection to both the data matrix Y and design matrix X.
#' Both must be projected to avoid bias in coefficient estimates.
#'
#' @param proj A soft_projection object from \code{soft_projection()}.
#' @param Y Data matrix (time x voxels).
#' @param X Design matrix (time x predictors).
#' @return List with projected Y and X.
#' @export
#' @examples
#' set.seed(123)
#' N <- matrix(rnorm(100 * 20), nrow = 100, ncol = 20)
#' Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
#' X <- cbind(1, rnorm(100))
#'
#' proj <- soft_projection(N, lambda = "auto")
#' cleaned <- apply_soft_projection(proj, Y, X)
#' # Use cleaned$Y and cleaned$X for GLM fitting
apply_soft_projection <- function(proj, Y, X) {
  if (!inherits(proj, "soft_projection")) {
    stop("proj must be a soft_projection object")
  }

  list(
    Y = proj$P_lambda(Y),
    X = proj$P_lambda(X)
  )
}

#' Residual degrees of freedom of an OLS fit after soft projection
#'
#' After \code{Y <- P Y}, \code{X <- P X} with
#' \eqn{P = I - U S U'} (\eqn{S = diag(d^2/(d^2+\lambda))}), the OLS residual
#' is \eqn{(I - H) P \epsilon}, where \eqn{H} is the hat matrix of
#' \eqn{PX}. With iid noise \eqn{E[RSS] = \sigma^2 tr((I - H) P^2)}, so the
#' unbiased residual df is
#' \deqn{n - \sum_j g_j - rank(PX) + \sum_j g_j \|Q' u_j\|^2,\quad g_j = 2 s_j - s_j^2,}
#' with \eqn{Q} an orthonormal basis of \eqn{PX}. It reduces to
#' \eqn{n - p} when \eqn{\lambda \to \infty} and to \eqn{n - k - p} for a hard
#' projection onto the complement of a \eqn{k}-dimensional nuisance space
#' that \eqn{X} avoids. The usual \code{n - rank(PX)} ignores the projection
#' and overstates the df, so \eqn{\hat\sigma^2} is biased low.
#'
#' With volume weights applied before the projection (\code{Y <- P D Y},
#' \eqn{D = diag(\sqrt{w})}), the WLS noise model gives
#' \eqn{Var(D\epsilon) = \sigma^2 A}, \eqn{A = diag(w > 0)}, so
#' \eqn{E[RSS] = \sigma^2 tr((I - H) P A P)}: the iid value minus
#' \eqn{\sum_{i: w_i = 0} \|(I - H) P e_i\|^2}. With every weight positive
#' this is the iid value; with no projection it is \code{sum(w > 0) - rank},
#' the volume-weighted df used without projection.
#'
#' @param proj A \code{soft_projection} object.
#' @param X_proj The projected design \code{proj$P_lambda(X)}.
#' @param tol Rank tolerance for \code{qr()}.
#' @param active Optional logical vector, \code{FALSE} for zero-weight volumes.
#' @return Scalar residual degrees of freedom (not necessarily an integer).
#' @keywords internal
#' @noRd
.soft_projection_residual_df <- function(proj, X_proj, tol = 1e-7,
                                         active = NULL) {
  if (!inherits(proj, "soft_projection")) {
    stop("proj must be a soft_projection object", call. = FALSE)
  }
  X_proj <- as.matrix(X_proj)
  n <- nrow(X_proj)
  d <- proj$svd$d
  U <- proj$svd$u[, seq_along(d), drop = FALSE]
  s <- d^2 / (d^2 + proj$lambda)
  g <- 2 * s - s^2

  qx <- qr(X_proj, tol = tol)
  r <- qx$rank
  Q <- if (r > 0L) qr.Q(qx)[, seq_len(r), drop = FALSE] else matrix(0, n, 0L)
  leverage_u <- if (r > 0L) {
    colSums(crossprod(Q, U)^2)                  # ||Q' u_j||^2
  } else {
    numeric(length(g))
  }
  df <- n - sum(g) - r + sum(g * leverage_u)

  zero <- if (is.null(active)) integer(0) else which(!as.logical(active))
  if (length(zero)) {
    # Columns P e_i of the zero-weight volumes, then (I - H) P e_i.
    E0 <- matrix(0, n, length(zero))
    E0[cbind(zero, seq_along(zero))] <- 1
    PE <- proj$P_lambda(E0)
    if (r > 0L) PE <- PE - Q %*% crossprod(Q, PE)
    df <- df - sum(PE * PE)
  }
  df
}

#' Residual degrees of freedom after soft projection and AR whitening
#'
#' The runwise AR path projects (\eqn{PY}, \eqn{PX}), estimates the AR model
#' from the projected OLS residuals and fits by OLS on the whitened data
#' \eqn{WPY}, \eqn{WPX}. With noise \eqn{\epsilon \sim (0, \sigma^2\Sigma)}
#' and \eqn{\Sigma = (W'W)^{-1}} -- the covariance the whitening filter
#' assumes, so that \eqn{W\epsilon} is white -- the whitened residual is
#' \eqn{R W P \epsilon} with \eqn{R = I - H_w} the residual projector of
#' \eqn{WPX}, hence
#' \deqn{E[RSS] = \sigma^2 tr(R M M'),\quad M = W P W^{-1} = I - A S B',}
#' with \eqn{A = WU}, \eqn{B = W^{-T}U} (\eqn{P = I - USU'}). Expanding,
#' \deqn{tr(RMM') = (n - r) - 2 tr(S B'RA) + tr(S B'B S A'RA),}
#' which needs only \eqn{n \times k} and \eqn{k \times k} products besides
#' one triangular solve with \eqn{W}. Without projection (\eqn{S = 0}) it is
#' \eqn{n - r}, the df the AR path uses; with \eqn{W = I} it is
#' \code{.soft_projection_residual_df()}.
#'
#' @param proj A \code{soft_projection} object.
#' @param X_white The whitened projected design \eqn{WPX}.
#' @param W The \eqn{n \times n} whitening matrix (the whitening transform
#'   applied to the identity).
#' @param tol Rank tolerance for \code{qr()}.
#' @return Scalar residual degrees of freedom.
#' @keywords internal
#' @noRd
.soft_projection_ar_residual_df <- function(proj, X_white, W, tol = 1e-7) {
  if (!inherits(proj, "soft_projection")) {
    stop("proj must be a soft_projection object", call. = FALSE)
  }
  X_white <- as.matrix(X_white)
  W <- as.matrix(W)
  n <- nrow(X_white)
  if (!identical(dim(W), c(n, n))) {
    stop("W must be an n x n whitening matrix", call. = FALSE)
  }
  d <- proj$svd$d
  U <- proj$svd$u[, seq_along(d), drop = FALSE]
  s <- d^2 / (d^2 + proj$lambda)

  A <- W %*% U
  B <- if (all(W[upper.tri(W)] == 0)) {
    forwardsolve(W, U, transpose = TRUE)        # W' B = U, W lower triangular
  } else {
    solve(t(W), U)
  }

  qx <- qr(X_white, tol = tol)
  r <- qx$rank
  if (r > 0L) {
    Q <- qr.Q(qx)[, seq_len(r), drop = FALSE]
    QA <- crossprod(Q, A)
    QB <- crossprod(Q, B)
    BRA <- crossprod(B, A) - crossprod(QB, QA)
    ARA <- crossprod(A) - crossprod(QA)
  } else {
    BRA <- crossprod(B, A)
    ARA <- crossprod(A)
  }
  G <- crossprod(B)
  SGS <- (s * G) * rep(s, each = length(s))     # S G S
  (n - r) - 2 * sum(s * diag(BRA)) + sum(SGS * ARA)
}

#' Residual-ACVF bias correction for AR estimation after soft projection
#'
#' AR coefficients are estimated from the OLS residuals of the projected run,
#' \eqn{r = G\epsilon} with \eqn{G = (I - H_{PX})P}. Their expected
#' autocovariance at lag \eqn{h} is linear in the true autocovariances,
#' \eqn{E\hat\gamma_r(h) = \sum_k A_{hk}\gamma(k)}; fmriAR's design-based
#' correction builds \eqn{A} for \eqn{G = I - H_X} only, so a projected fit
#' would be corrected for the wrong operator and the AR coefficients biased
#' (towards zero for AR(1) noise). This builds \eqn{A} for the projected
#' operator with the same convention as fmriAR (residuals demeaned over the
#' kept volumes, lag pairs within censoring segments), for
#' `fmriAR::fit_noise(acvf_correction = )`.
#'
#' @param proj A \code{soft_projection} object.
#' @param X_proj The projected run design.
#' @param max_lag Correction lag budget.
#' @param censor Optional censored volume indices (or logical mask).
#' @param tol Rank tolerance for \code{qr()}.
#' @return A \code{(max_lag + 1)} square matrix.
#' @keywords internal
#' @noRd
.soft_projection_acvf_correction <- function(proj, X_proj, max_lag,
                                             censor = NULL, tol = 1e-7) {
  X_proj <- as.matrix(X_proj)
  n <- nrow(X_proj)
  max_lag <- max(0L, as.integer(max_lag))
  G <- proj$P_lambda(diag(n))
  qx <- qr(X_proj, tol = tol)
  if (qx$rank > 0L) {
    Q <- qr.Q(qx)[, seq_len(qx$rank), drop = FALSE]
    G <- G - Q %*% crossprod(Q, G)
  }
  if (is.logical(censor)) censor <- which(censor)
  keep <- setdiff(seq_len(n), as.integer(censor %||% integer(0)))
  nv <- length(keep)
  if (nv < 2L) return(diag(max_lag + 1L))
  seg_id <- cumsum(c(1L, as.integer(diff(keep) > 1L)))
  R <- G[keep, , drop = FALSE]
  R <- R - rep(colMeans(R), each = nv)

  # R T_k with T_k the symmetric lag-k shift (T_0 = I).
  lag_apply <- function(B, k) {
    if (k == 0L) return(B)
    out <- matrix(0, nrow(B), n)
    if (k >= n) return(out)
    hi <- seq.int(k + 1L, n)
    lo <- seq.int(1L, n - k)
    out[, hi] <- out[, hi] + B[, lo, drop = FALSE]
    out[, lo] <- out[, lo] + B[, hi, drop = FALSE]
    out
  }
  pairs <- lapply(0:max_lag, function(h) {
    if (h == 0L) return(list(a = seq_len(nv), b = seq_len(nv)))
    if (nv <= h) return(NULL)
    hi <- seq.int(h + 1L, nv)
    lo <- seq.int(1L, nv - h)
    ok <- seg_id[hi] == seg_id[lo]
    list(a = hi[ok], b = lo[ok])
  })
  A <- diag(max_lag + 1L)
  for (k in 0:max_lag) {
    SkR <- lag_apply(R, k)
    for (h in 0:max_lag) {
      pr <- pairs[[h + 1L]]
      if (is.null(pr) || !length(pr$a)) next
      A[h + 1L, k + 1L] <- sum(R[pr$a, , drop = FALSE] *
                                 SkR[pr$b, , drop = FALSE]) / length(pr$a)
    }
  }
  A
}

#' Unscaled coefficient covariance of a fit after soft projection
#'
#' The projected (and possibly whitened) noise is not white, so the model
#' covariance \eqn{(X_f'X_f)^{-1}} of the fitted design \eqn{X_f} is not the
#' covariance of the estimate. With fitted noise \eqn{T\epsilon},
#' \eqn{Var(\epsilon) = \sigma^2\Sigma}, the estimate has
#' \deqn{Var(\hat\beta) = \sigma^2 C X_f' T\Sigma T' X_f C,\quad C = (X_f'X_f)^{-1}.}
#' `Z` is a factor with \eqn{Z'Z = X_f'T\Sigma T'X_f}: \eqn{Z = A^{1/2}PX_f}
#' for iid noise (\eqn{T = PD}, \eqn{A} the non-zero volume weights) and
#' \eqn{Z = W^{-T}PW'X_f} for AR noise (\eqn{T = WP}, \eqn{\Sigma = (W'W)^{-1}}).
#' Together with the expectation-matched residual df this makes
#' \eqn{\hat\sigma^2 c'C_s c} unbiased for the contrast variance. The rank
#' attributes of `XtXinv` (aliased columns) are kept.
#'
#' @param XtXinv \eqn{(X_f'X_f)^{-1}} from `.fast_preproject()`.
#' @param Z Factor described above (\eqn{n \times p}).
#' @return Unscaled covariance with the attributes of `XtXinv`.
#' @keywords internal
#' @noRd
.soft_projection_coef_cov <- function(XtXinv, Z) {
  C <- as.matrix(XtXinv)
  Z <- as.matrix(Z)
  Z[!is.finite(Z)] <- 0
  C0 <- C
  C0[!is.finite(C0)] <- 0
  out <- C0 %*% crossprod(Z) %*% C0
  out <- (out + t(out)) / 2
  out[!is.finite(C)] <- C[!is.finite(C)]
  attributes(out) <- attributes(XtXinv)
  out
}

#' ACVF correction for one projected run, or NULL without projection
#'
#' The lag budget is fmriAR's design-based budget, capped so that about three
#' residual df (after projection) back each corrected lag.
#' @keywords internal
#' @noRd
.soft_run_acvf_correction <- function(proj, X_proj, ar_order, censor = NULL) {
  if (is.null(proj)) return(NULL)
  budget <- .ar_correction_lag_budget(X_proj, target_order = ar_order,
                                      censor = censor)
  df_soft <- .soft_projection_residual_df(proj, X_proj)
  budget <- as.integer(max(ar_order, min(budget, floor(df_soft / 3))))
  .soft_projection_acvf_correction(proj, X_proj, budget, censor = censor)
}

#' Extract Nuisance Timeseries from Mask
#'
#' Extracts voxel timeseries from regions defined by a mask (e.g., WM/CSF).
#' This is the typical input for soft subspace projection.
#'
#' @param dataset An fmri_frame object.
#' @param mask A binary mask (logical vector or 3D array) indicating nuisance voxels,
#'   or a file path to a NIfTI mask.
#' @param run Optional run index to extract data from a specific run.
#' @return Matrix of nuisance timeseries (time x nuisance_voxels).
#' @export
extract_nuisance_timeseries <- function(dataset, mask, run = NULL) {
  # Handle mask input
  if (is.character(mask)) {
    # Load from file
    if (!requireNamespace("RNifti", quietly = TRUE)) {
      stop("RNifti package required to load NIfTI masks")
    }
    mask_img <- RNifti::readNifti(mask)
    mask <- as.logical(as.vector(mask_img))
  }

  if (is.array(mask) && length(dim(mask)) == 3) {
    mask <- as.logical(as.vector(mask))
  }

  # Get data
  Y <- .dset_data_matrix(dataset)

  if (!is.null(run)) {
    # Extract run-specific data
    sframe <- .dset_sampling_frame(dataset)
    blocklens <- fmrihrf::blocklens(sframe)
    run_start <- if (run == 1) 1 else sum(blocklens[1:(run - 1)]) + 1
    run_end <- sum(blocklens[1:run])
    Y <- Y[run_start:run_end, , drop = FALSE]
  }

  # Subset to nuisance voxels
  if (length(mask) != ncol(Y)) {
    stop("Mask length (", length(mask), ") must match number of voxels (",
         ncol(Y), ")")
  }

  Y[, mask, drop = FALSE]
}

#' Soft Subspace Control Options
#'
#' Creates a configuration object for soft subspace projection.
#'
#' @param enabled Logical. Whether to apply soft subspace projection.
#' @param nuisance_mask Path to NIfTI mask or logical vector indicating nuisance voxels.
#' @param nuisance_matrix Pre-computed nuisance timeseries matrix (alternative to mask).
#' @param lambda Ridge penalty: numeric, "auto", or "gcv".
#' @param warn_redundant Logical. Warn if baseline model contains nuisance terms.
#' @return A list of class "soft_subspace_options".
#' @export
#' @examples
#' # Using a mask file
#' opts <- soft_subspace_options(
#'   enabled = TRUE,
#'   nuisance_mask = "path/to/wm_csf_mask.nii.gz",
#'   lambda = "auto"
#' )
#'
#' # Using pre-computed nuisance matrix
#' N <- matrix(rnorm(100 * 20), 100, 20)
#' opts <- soft_subspace_options(
#'   enabled = TRUE,
#'   nuisance_matrix = N,
#'   lambda = 0.5
#' )
soft_subspace_options <- function(enabled = FALSE,
                                   nuisance_mask = NULL,
                                   nuisance_matrix = NULL,
                                   lambda = "auto",
                                   warn_redundant = TRUE) {
  if (enabled) {
    if (is.null(nuisance_mask) && is.null(nuisance_matrix)) {
      stop("Either nuisance_mask or nuisance_matrix must be provided when enabled=TRUE")
    }
    if (!is.null(nuisance_mask) && !is.null(nuisance_matrix)) {
      warning("Both nuisance_mask and nuisance_matrix provided; using nuisance_matrix")
    }
  }

  structure(
    list(
      enabled = enabled,
      nuisance_mask = nuisance_mask,
      nuisance_matrix = nuisance_matrix,
      lambda = lambda,
      warn_redundant = warn_redundant
    ),
    class = "soft_subspace_options"
  )
}

#' Check for Redundant Nuisance Regressors
#'
#' Checks if the baseline model contains nuisance regressors when
#' soft subspace projection is enabled, and issues a warning.
#'
#' @param baseline_model A baseline_model object.
#' @param soft_opts Soft subspace options.
#' @return Invisible NULL. Issues warning if redundancy detected.
#' @keywords internal
#' @noRd
.check_redundant_nuisance <- function(baseline_model, soft_opts) {
  if (!soft_opts$enabled || !soft_opts$warn_redundant) {
    return(invisible(NULL))
  }

  # Check if baseline model has nuisance terms
  has_nuisance <- !is.null(baseline_model$nuisance_term)

  if (has_nuisance) {
    warning(
      "Soft subspace projection is enabled, but baseline model contains nuisance regressors.\n",
      "Soft projection removes variance explainable by noise voxels.\n",
      "Additional nuisance regressors may be redundant.\n",
      "If intentional (e.g., motion parameters), set warn_redundant=FALSE to suppress.",
      call. = FALSE
    )
  }

  invisible(NULL)
}

#' Apply Soft Projection in Pipeline
#'
#' Internal function called by fmri_lm to apply soft subspace projection
#' to data and design matrix.
#'
#' @param Y Data matrix.
#' @param X Design matrix.
#' @param dataset fmri_frame for extracting nuisance (if using mask).
#' @param soft_opts Soft subspace options.
#' @param run Optional run index.
#' @return List with projected Y, X, and projection details.
#' @keywords internal
#' @noRd
.apply_soft_projection_pipeline <- function(Y, X, dataset, soft_opts, run = NULL) {
  if (!soft_opts$enabled) {
    return(list(Y = Y, X = X, projection = NULL))
  }

  # Get nuisance matrix
  if (!is.null(soft_opts$nuisance_matrix)) {
    N <- soft_opts$nuisance_matrix
    if (!is.null(run)) {
      # Need to subset to run
      sframe <- .dset_sampling_frame(dataset)
      blocklens <- fmrihrf::blocklens(sframe)
      run_start <- if (run == 1) 1 else sum(blocklens[1:(run - 1)]) + 1
      run_end <- sum(blocklens[1:run])
      N <- N[run_start:run_end, , drop = FALSE]
    }
  } else {
    N <- extract_nuisance_timeseries(dataset, soft_opts$nuisance_mask, run = run)
  }

  # Check dimensions
  if (nrow(N) != nrow(Y)) {
    stop("Nuisance matrix rows (", nrow(N), ") must match data rows (", nrow(Y), ")")
  }

  # Create projection
  proj <- soft_projection(N, lambda = soft_opts$lambda, Y = Y)

  # Apply to data and design
  cleaned <- apply_soft_projection(proj, Y, X)

  list(Y = cleaned$Y, X = cleaned$X, projection = proj)
}
