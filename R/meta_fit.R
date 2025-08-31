# Thin R Wrappers for C++ Meta-Analysis Functions

#' Fit Meta-Analysis Models
#'
#' Low-level function that calls the C++ meta-analysis implementation.
#' This is typically called internally by higher-level functions like fmri_meta().
#'
#' @param Y Numeric matrix of effect sizes (subjects x features)
#' @param V Numeric matrix of variances (subjects x features)
#' @param X Numeric matrix; design matrix (subjects x predictors), including intercept
#' @param method Character scalar; meta-analysis method: "pm" (Paule-Mandel), "dl" (DerSimonian-Laird),
#'   "fe" (fixed-effects), or "reml" (REML, uses PM solver)
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
  
  method <- match.arg(method)
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
#' @export
z_to_r <- function(z) {
  (exp(2 * z) - 1) / (exp(2 * z) + 1)
}