# Thin R Wrappers for C++ Meta-Analysis Functions

# Internal: "reml" has never had its own estimator; it is accepted for
# backward compatibility and runs the Paule-Mandel estimator. Say so and
# record the estimator that actually ran.
.meta_resolve_method <- function(method) {
  if (identical(method, "reml")) {
    message("method = 'reml' is not implemented; using the Paule-Mandel ",
            "estimator (method = 'pm') instead.")
    method <- "pm"
  }
  method
}

# Internal: relative weights (S x P) for weights = "equal" / "custom".
.meta_relative_weights <- function(weights, weights_custom, S, P) {
  if (identical(weights, "equal")) {
    return(matrix(1, nrow = S, ncol = P))
  }
  w <- .fmri_ttest_validate_weights_custom(weights_custom, S, P)
  if (is.matrix(w)) w else matrix(w, nrow = S, ncol = P)
}

# Internal: the single meta-regression fitter used by every fmri_meta /
# fmri_ttest path that has subject-level effects.
#
# weights = "ivw": the given sampling variances V are absolute; `method`
#   chooses the between-subject variance estimator (fe, dl, pm).
# weights = "equal" / "custom": the weights are *relative*. The model is the
#   weighted least-squares model Var(y_i) = sigma^2 / w_i with the residual
#   scale sigma^2 estimated per feature (as lm(weights = w) does), so `method`
#   and the supplied variances are not used and tau2 / Q / I2 are NA.
#
# Always returns the packed covariance of the coefficients (cov_tri).
.meta_fit_core <- function(Y, V, X, method = "pm", robust = "none",
                           weights = "ivw", weights_custom = NULL,
                           huber_c = 1.345, robust_iter = 2,
                           n_threads = getOption("fmrireg.num_threads", 0)) {
  Y <- as.matrix(Y)
  X <- as.matrix(X)
  storage.mode(X) <- "double"
  S <- nrow(Y)
  P <- ncol(Y)
  if (identical(weights, "ivw")) {
    V <- as.matrix(V)
    res <- meta_fit_cov_cpp(Y, V, X, method, robust, huber_c, robust_iter, n_threads)
    res$ok <- as.logical(res$ok)
    return(res)
  }

  W <- .meta_relative_weights(weights, weights_custom, S, P)
  V_rel <- 1 / W
  # Weighted least squares with unit scale: Q_fe is the weighted residual
  # sum of squares, so sigma^2 = Q_fe / df.
  res <- meta_fit_cov_cpp(Y, V_rel, X, "fe", "none", huber_c, robust_iter, n_threads)
  sigma2 <- as.numeric(res$Q_fe) / as.numeric(res$df)
  if (identical(robust, "huber")) {
    V_s <- sweep(V_rel, 2L, sigma2, "*")
    res <- meta_fit_cov_cpp(Y, V_s, X, "fe", robust, huber_c, robust_iter, n_threads)
  } else {
    scale <- sqrt(sigma2)
    res$se <- sweep(res$se, 2L, scale, "*")
    res$z <- sweep(res$z, 2L, scale, "/")
    res$cov_tri <- sweep(res$cov_tri, 2L, sigma2, "*")
  }
  res$sigma2 <- sigma2
  res$tau2[] <- NA_real_
  res$Q_fe[] <- NA_real_
  res$I2_fe[] <- NA_real_
  res$ok <- as.logical(res$ok)
  res
}

# Internal: contrasts c' beta with SE sqrt(c' Var(beta) c) from packed
# covariance triangles. C is K x J, beta is K x P, cov_tri is tsize x P.
.meta_contrasts_from_tri <- function(beta, cov_tri, C) {
  C <- as.matrix(C)
  K <- nrow(C)
  idx <- which(upper.tri(diag(K), diag = TRUE), arr.ind = TRUE)
  # packed order is row-major over the upper triangle (a <= b)
  idx <- idx[order(idx[, "row"], idx[, "col"]), , drop = FALSE]
  mult <- ifelse(idx[, "row"] == idx[, "col"], 1, 2)
  est <- crossprod(C, beta)                                   # J x P
  coef_w <- (C[idx[, "row"], , drop = FALSE] * C[idx[, "col"], , drop = FALSE]) * mult
  vv <- crossprod(coef_w, cov_tri)                             # J x P
  se <- sqrt(ifelse(is.finite(vv) & vv >= 0, vv, NA_real_))
  list(estimate = est, se = se, z = est / se)
}

#' Fit Meta-Analysis Models
#'
#' Low-level function that calls the C++ meta-analysis implementation.
#' This is typically called internally by higher-level functions like fmri_meta().
#'
#' @param Y Numeric matrix of effect sizes (subjects x features)
#' @param V Numeric matrix of variances (subjects x features)
#' @param X Numeric matrix; design matrix (subjects x predictors), including intercept
#' @param method Character scalar; meta-analysis method: "pm" (Paule-Mandel), "dl" (DerSimonian-Laird),
#'   "fe" (fixed-effects). "reml" is not implemented: it is accepted for
#'   backward compatibility, runs Paule-Mandel with a message, and the result
#'   records \code{method = "pm"}.
#' @param robust Character scalar; robust estimation method: "none" or "huber"
#' @param huber_c Numeric scalar; tuning constant for Huber M-estimator (default: 1.345).
#'   Smaller values provide more robust estimates but may reduce efficiency.
#' @param robust_iter Integer scalar; number of IRLS iterations for robust estimation (default: 2)
#' @param n_threads Integer scalar; number of OpenMP threads (0 = use all available)
#'
#' @return List with components:
#'   \item{beta}{Numeric matrix of coefficients (predictors x features)}
#'   \item{se}{Numeric matrix of standard errors (predictors x features)}
#'   \item{z}{Numeric matrix of z-scores (predictors x features)}
#'   \item{tau2}{Numeric vector of between-study variance estimates}
#'   \item{Q_fe}{Numeric vector of Q statistics from fixed-effects model}
#'   \item{I2_fe}{Numeric vector of I² statistics from fixed-effects model}
#'   \item{df}{Numeric vector of degrees of freedom}
#'   \item{ok}{Logical vector indicating successful fits}
#'
#' @seealso \code{\link{fmri_meta}}
#' @export
fmri_meta_fit <- function(Y, V, X, 
                         method = c("pm", "dl", "fe", "reml"),
                         robust = c("none", "huber"),
                         huber_c = 1.345, 
                         robust_iter = 2,
                         n_threads = getOption("fmrireg.num_threads", 0)) {
  
  method <- .meta_resolve_method(match.arg(method))
  robust <- match.arg(robust)
  
  # Validate inputs
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.matrix(V)) V <- as.matrix(V)
  if (!is.matrix(X)) X <- as.matrix(X)
  
  if (nrow(Y) != nrow(V)) {
    stop("Y and V must have the same number of rows (subjects)", call. = FALSE)
  }
  
  if (ncol(Y) != ncol(V)) {
    stop("Y and V must have the same number of columns (features)", call. = FALSE)
  }
  
  if (nrow(X) != nrow(Y)) {
    stop("X must have the same number of rows as Y (subjects)", call. = FALSE)
  }
  
  # Check for valid variances
  if (any(V <= 0, na.rm = TRUE)) {
    warning("Non-positive variances detected; these will be treated as missing", call. = FALSE)
    V[V <= 0] <- NA
  }
  
  # Call C++ implementation
  result <- .Call("_fmrireg_meta_fit_cpp", 
                  Y, V, X, method, robust, huber_c, robust_iter, n_threads,
                  PACKAGE = "fmrireg")
  
  # Add method and robust info to result
  result$method <- method
  result$robust <- robust
  
  # Convert OK vector to logical
  result$ok <- as.logical(result$ok)
  
  return(result)
}

#' Fit Meta-Analysis Models with Exact Contrasts
#'
#' Low-level function that calls the C++ meta-analysis implementation and returns
#' exact contrast statistics c' (X' W X)^(-1) c for provided contrasts.
#'
#' @param Y,V,X See fmri_meta_fit
#' @param Cmat Numeric matrix K x J (columns are contrasts over predictors)
#' @inheritParams fmri_meta_fit
#' @return List with base meta outputs plus c_beta, c_se, c_z
#' @export
fmri_meta_fit_contrasts <- function(Y, V, X, Cmat,
                                   method = c("pm", "dl", "fe", "reml"),
                                   robust = c("none", "huber"),
                                   huber_c = 1.345,
                                   robust_iter = 2,
                                   n_threads = getOption("fmrireg.num_threads", 0)) {
  method <- .meta_resolve_method(match.arg(method))
  robust <- match.arg(robust)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.matrix(V)) V <- as.matrix(V)
  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Cmat)) Cmat <- as.matrix(Cmat)
  if (nrow(Cmat) != ncol(X)) stop("nrow(Cmat) must equal ncol(X)", call. = FALSE)
  result <- meta_fit_contrasts_cpp(Y, V, X, Cmat, method, robust, huber_c, robust_iter, n_threads)
  result$method <- method
  result$robust <- robust
  result$ok <- as.logical(result$ok)
  result
}

#' Fit Meta-Analysis and return packed covariance per voxel
#'
#' @inheritParams fmri_meta_fit
#' @return List with base outputs and cov_tri (tsize x P) where tsize = K*(K+1)/2
#' @export
fmri_meta_fit_cov <- function(Y, V, X,
                             method = c("pm", "dl", "fe", "reml"),
                             robust = c("none", "huber"),
                             huber_c = 1.345,
                             robust_iter = 2,
                             n_threads = getOption("fmrireg.num_threads", 0)) {
  method <- .meta_resolve_method(match.arg(method))
  robust <- match.arg(robust)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.matrix(V)) V <- as.matrix(V)
  if (!is.matrix(X)) X <- as.matrix(X)
  result <- meta_fit_cov_cpp(Y, V, X, method, robust, huber_c, robust_iter, n_threads)
  result$method <- method
  result$robust <- robust
  result$ok <- as.logical(result$ok)
  result
}

#' Compute Effective Sample Size for Meta-Analysis
#'
#' Computes the effective sample size based on the heterogeneity estimate.
#' This is useful for understanding the impact of between-study heterogeneity.
#'
#' @param v Numeric vector of within-study variances
#' @param tau2 Numeric scalar; between-study variance (tau-squared)
#'
#' @return Numeric scalar; effective sample size
#'
#' @seealso \code{\link{fmri_meta}}
#' @examples
#' meta_effective_n(v = rep(0.05, 3), tau2 = 0.01)
#' @export
meta_effective_n <- function(v, tau2) {
  w <- 1 / (v + tau2)
  n_eff <- sum(w)^2 / sum(w^2)
  return(n_eff)
}

#' Convert T-statistics to Effect Sizes and Variances
#'
#' Converts t-statistics and degrees of freedom to standardized mean differences
#' (Cohen's d) and their sampling variances for meta-analysis.
#'
#' @param t Numeric vector or matrix of t-statistics
#' @param df Numeric scalar or vector; degrees of freedom (matching t)
#' @param n Numeric scalar; sample size per group (for two-sample t-tests)
#'
#' @return List with components:
#'   \item{d}{Numeric vector or matrix; standardized mean differences}
#'   \item{v}{Numeric vector or matrix; sampling variances}
#'
#' @seealso \code{\link{fmri_meta}}
#' @examples
#' t_to_d(t = 2, df = 18)
#' t_to_d(t = 2, df = 18, n = 20)
#' @export
t_to_d <- function(t, df, n = NULL) {
  if (is.null(n)) {
    # One-sample or paired t-test
    n <- df + 1
    d <- t * sqrt(1/n)
    v <- 1/n + d^2/(2*n)
  } else {
    # Two-sample t-test
    d <- 2 * t / sqrt(df)
    v <- 4/n + d^2/(2*df)
  }
  
  return(list(d = d, v = v))
}

#' Convert Correlation to Fisher's Z
#'
#' Transforms correlations to Fisher's Z scale for meta-analysis.
#'
#' @param r Numeric vector or matrix of correlations
#' @param n Integer scalar; sample size
#'
#' @return List with components:
#'   \item{z}{Numeric vector or matrix; Fisher's Z transformed correlations}
#'   \item{v}{Numeric vector or matrix; sampling variances}
#'
#' @seealso \code{\link{fmri_meta}}
#' @examples
#' r_to_z(r = 0.4, n = 30)
#' @export
r_to_z <- function(r, n) {
  z <- 0.5 * log((1 + r) / (1 - r))
  v <- 1 / (n - 3)
  return(list(z = z, v = v))
}

#' Back-transform Fisher's Z to Correlation
#'
#' @param z Numeric vector or matrix; Fisher's Z values
#' @return Numeric vector or matrix; correlations
#'
#' @seealso \code{\link{fmri_meta}}
#' @examples
#' z_to_r(0.2)
#' @export
z_to_r <- function(z) {
  (exp(2 * z) - 1) / (exp(2 * z) + 1)
}
