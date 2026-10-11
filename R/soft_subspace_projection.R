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
#' @param proj A \code{soft_projection} object.
#' @param X_proj The projected design \code{proj$P_lambda(X)}.
#' @param tol Rank tolerance for \code{qr()}.
#' @return Scalar residual degrees of freedom (not necessarily an integer).
#' @keywords internal
#' @noRd
.soft_projection_residual_df <- function(proj, X_proj, tol = 1e-7) {
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
  if (r > 0L) {
    Q <- qr.Q(qx)[, seq_len(r), drop = FALSE]
    leverage_u <- colSums(crossprod(Q, U)^2)   # ||Q' u_j||^2
  } else {
    leverage_u <- numeric(length(g))
  }
  n - sum(g) - r + sum(g * leverage_u)
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
