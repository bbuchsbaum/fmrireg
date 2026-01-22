#' Volume Quality Metrics and Temporal Weighting
#'
#' Functions for computing volume-level quality metrics (DVARS) and
#' converting them to weights for weighted least squares fitting.
#' This implements "soft scrubbing" - downweighting bad volumes rather
#' than hard censoring them.
#'
#' @section Conceptual Overview:
#' Traditional fMRI artifact removal uses hard censoring ("scrubbing"),
#' where volumes exceeding a threshold are completely removed. This
#' loses temporal information and can create discontinuities.
#'
#' Soft scrubbing instead assigns each volume a weight between 0 and 1
#' based on its quality. High-quality volumes receive full weight;
#' artifacts receive reduced weight. This preserves temporal continuity
#' while still downweighting problematic data.
#'
#' @section DVARS as Quality Metric:
#' DVARS (Derivative of VARiance across voxels) measures the root mean square
#' of the temporal derivative across all voxels:
#' \deqn{DVARS_t = \sqrt{\frac{1}{V} \sum_v (Y_{t,v} - Y_{t-1,v})^2}}
#'
#' Interpretation:
#' \itemize{
#'   \item Low DVARS: signal changed smoothly from previous volume (good)
#'   \item High DVARS: signal changed rapidly (possible artifact or motion)
#'   \item Normalized DVARS: values ~1 are typical; >1.5 suggests artifacts
#' }
#'
#' @section Weighting Methods:
#' Three methods convert DVARS to weights, offering different trade-offs:
#'
#' \describe{
#'   \item{inverse_squared}{
#'     Formula: \code{w = 1 / (1 + dvars^2)}
#'
#'     Properties: Smooth, continuous decay. A volume with DVARS=1 gets
#'     weight 0.5. Most conservative choice - provides gentle downweighting
#'     even for moderately elevated DVARS.
#'
#'     Use when: You want smooth, gradual artifact handling without
#'     sharp transitions.
#'   }
#'
#'   \item{soft_threshold}{
#'     Formula: Sigmoid decay above threshold
#'
#'     Properties: Volumes below threshold get full weight; above
#'     threshold, weights decay smoothly. Steepness parameter controls
#'     how rapidly weights drop.
#'
#'     Use when: You have a clear idea of what "acceptable" DVARS looks
#'     like (e.g., 1.5x median) and want to preserve good volumes fully.
#'   }
#'
#'   \item{tukey}{
#'     Formula: Tukey bisquare \code{(1 - u^2)^2} for \code{|u| <= 1}
#'
#'     Properties: Complete downweighting for extreme values. Volumes
#'     beyond 2x threshold get zero weight. Most aggressive choice.
#'
#'     Use when: You have clear artifacts that should be fully excluded,
#'     similar to hard scrubbing but with smooth transitions.
#'   }
#' }
#'
#' @section When to Use Volume Weighting vs Robust Fitting:
#' Both approaches handle artifacts, but through different mechanisms:
#'
#' \strong{Volume weighting} downweights entire timepoints uniformly across
#' all voxels. Best when artifacts affect the whole brain (head motion,
#' scanner spikes).
#'
#' \strong{Robust fitting} (Huber/Tukey bisquare) downweights outlier
#' residuals voxel-by-voxel. Best when artifacts are spatially localized
#' or when temporal structure matters.
#'
#' \strong{Combined approach}: Use volume weighting for global artifacts
#' plus robust fitting for residual voxel-level outliers. Set
#' \code{robust = TRUE, volume_weights = TRUE} in \code{fmri_lm()}.
#'
#' @section Typical Usage:
#' For use within \code{fmri_lm()}, see the convenience parameter
#' \code{volume_weights} or the more detailed \code{volume_weights_options}.
#'
#' @name volume_quality
NULL

#' Compute DVARS (Temporal Derivative of Timecourses)
#'
#' DVARS measures the root mean square of the temporal derivative across
#' all voxels, providing a single quality metric per volume. High DVARS
#' indicates rapid signal changes, often associated with motion artifacts.
#'
#' @param Y Numeric matrix of fMRI data (time x voxels).
#' @param normalize Logical. If TRUE, normalize by median DVARS. Default TRUE.
#' @return Numeric vector of length nrow(Y) with DVARS values.
#'   The first timepoint is set to the median of subsequent values.
#' @export
#' @examples
#' # Simulate fMRI data with a motion spike
#' set.seed(123)
#' Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
#' Y[50, ] <- Y[50, ] + 5  # Add spike at volume 50
#' dvars <- compute_dvars(Y)
#' # dvars[50] will be elevated compared to other volumes
#'
#' # Identify suspect volumes (normalized DVARS > 1.5)
#' suspect <- which(dvars > 1.5)
#' cat("Suspect volumes:", suspect, "\n")
compute_dvars <- function(Y, normalize = TRUE) {
  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }

  n <- nrow(Y)
  if (n < 2) {
    stop("DVARS requires at least 2 timepoints")
  }


  # Temporal derivative: Y[t] - Y[t-1]
  dY <- diff(Y)  # (n-1) x voxels

  # RMS across voxels for each timepoint
  dvars_raw <- sqrt(rowMeans(dY^2))

  # First timepoint has no derivative; use median of subsequent

  dvars <- c(median(dvars_raw), dvars_raw)

  if (normalize) {
    med <- median(dvars)
    if (med > 0) {
      dvars <- dvars / med
    }
  }

  dvars
}

#' Convert DVARS to Volume Weights
#'
#' Transforms DVARS quality metrics into weights for weighted least squares
#' fitting. Volumes with high DVARS receive lower weights, implementing
#' "soft scrubbing" without hard censoring thresholds.
#'
#' @param dvars Numeric vector of DVARS values (from \code{compute_dvars}).
#' @param method Character. Weighting method: "inverse_squared" (default),
#'   "soft_threshold", or "tukey".
#' @param threshold Numeric. For "soft_threshold", the DVARS value above which
#'   weights decay. Default 1.5 (1.5x median if normalized).
#' @param steepness Numeric. For "soft_threshold", controls decay rate. Default 2.
#' @return Numeric vector of weights in the interval from 0 to 1, same length as dvars.
#' @export
#' @examples
#' set.seed(123)
#' Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
#' Y[50, ] <- Y[50, ] + 5
#' dvars <- compute_dvars(Y)
#'
#' # Compare different weighting methods
#' w_inv <- dvars_to_weights(dvars, method = "inverse_squared")
#' w_soft <- dvars_to_weights(dvars, method = "soft_threshold")
#' w_tukey <- dvars_to_weights(dvars, method = "tukey")
#'
#' # Check weight at the artifact volume
#' cat("Weights at volume 50:\n")
#' cat("  inverse_squared:", round(w_inv[50], 3), "\n")
#' cat("  soft_threshold:", round(w_soft[50], 3), "\n")
#' cat("  tukey:", round(w_tukey[50], 3), "\n")
dvars_to_weights <- function(dvars,
                              method = c("inverse_squared", "soft_threshold", "tukey"),
                              threshold = 1.5,
                              steepness = 2) {
  method <- match.arg(method)

  if (any(dvars < 0)) {
    stop("DVARS values must be non-negative")
  }

  weights <- switch(method,
    "inverse_squared" = {
      # w = 1 / (1 + dvars^2)
      # Simple and smooth; dvars=1 gives w=0.5
      1 / (1 + dvars^2)
    },
    "soft_threshold" = {
      # Sigmoid-like decay above threshold
      # w = 1 for dvars < threshold, decays smoothly above
      1 / (1 + ((pmax(dvars, threshold) - threshold) / threshold)^steepness)
    },
    "tukey" = {
      # Tukey bisquare: complete downweighting for extreme values
      c_tukey <- threshold * 2  # Cutoff at 2x threshold
      u <- dvars / c_tukey
      ifelse(abs(u) <= 1, (1 - u^2)^2, 0)
    }
  )

  # Ensure weights are in [0, 1]
  weights <- pmax(pmin(weights, 1), 0)

  # Normalize so mean weight is approximately 1 (preserves scale)
  weights <- weights / mean(weights)

  weights
}

#' Compute Volume Quality Weights from Data
#'
#' Convenience function that computes DVARS and converts to weights in one step.
#' This is the main user-facing function for volume quality weighting.
#'
#' @param Y Numeric matrix of fMRI data (time x voxels).
#' @param method Weighting method passed to \code{dvars_to_weights}.
#' @param threshold Threshold passed to \code{dvars_to_weights}.
#' @param return_dvars Logical. If TRUE, return both weights and DVARS.
#' @return If return_dvars is FALSE, a numeric vector of weights.
#'   If TRUE, a list with components "weights" and "dvars".
#' @export
#' @examples
#' set.seed(123)
#' Y <- matrix(rnorm(100 * 50), nrow = 100, ncol = 50)
#' Y[50, ] <- Y[50, ] + 5  # Add artifact
#'
#' # One-step computation of weights
#' result <- volume_weights(Y, return_dvars = TRUE)
#' cat("DVARS at artifact:", round(result$dvars[50], 2), "\n")
#' cat("Weight at artifact:", round(result$weights[50], 3), "\n")
#'
#' # With Tukey method for more aggressive downweighting
#' w_tukey <- volume_weights(Y, method = "tukey")
volume_weights <- function(Y,
                           method = "inverse_squared",
                           threshold = 1.5,
                           return_dvars = FALSE) {
  dvars <- compute_dvars(Y, normalize = TRUE)
  weights <- dvars_to_weights(dvars, method = method, threshold = threshold)

  if (return_dvars) {
    list(weights = weights, dvars = dvars)
  } else {
    weights
  }
}

#' Apply Volume Weights to GLM
#'
#' Internal function to apply volume-level weights to the design matrix
#' and data for weighted least squares fitting.
#'
#' @param X Design matrix (time x predictors).
#' @param Y Data matrix (time x voxels).
#' @param weights Numeric vector of volume weights.
#' @return List with weighted X and Y matrices.
#' @keywords internal
#' @noRd
apply_volume_weights <- function(X, Y, weights) {
  if (length(weights) != nrow(X)) {
    stop("Length of weights must equal number of timepoints")
  }

  # Square root of weights for WLS transformation
  sqrt_w <- sqrt(weights)

  # Weight both design and data
  X_weighted <- X * sqrt_w
  Y_weighted <- sweep(Y, 1, sqrt_w, `*`)

  list(X = X_weighted, Y = Y_weighted, sqrt_weights = sqrt_w)
}
