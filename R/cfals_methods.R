#' Methods for fmrireg_cfals_fit Objects
#'
#' Basic utilities for inspecting the results
#' of `fmrireg_hrf_cfals`.
#'
#' @param x,object An `fmrireg_cfals_fit` object.
#' @param vox Index of the voxel to plot.
#' @param ... Additional arguments passed to underlying functions.
#' @export
print.fmrireg_cfals_fit <- function(x, ...) {
  cat("\nfmrireg CF-ALS Fit\n")
  cat("==================\n")
  info <- x$design_info
  cat(sprintf("Voxels: %d\n", info$v))
  cat(sprintf("Time points: %d\n", info$n))
  cat(sprintf("Conditions: %d\n", info$k))
  cat(sprintf("Basis functions: %d\n", info$d))
  invisible(x)
}

#' @export
summary.fmrireg_cfals_fit <- function(object, ...) {
  res <- list(r2 = object$gof_per_voxel,
              design = object$design_info,
              lambdas = object$lambda_used)
  class(res) <- "summary.fmrireg_cfals_fit"
  res
}

#' @export
residuals.fmrireg_cfals_fit <- function(object, ...) {
  object$residuals
}

#' @export
plot.fmrireg_cfals_fit <- function(x, vox = 1, ...) {
  if (vox < 1 || vox > ncol(x$reconstructed_hrfs))
    stop("'vox' out of range")
  hrf <- x$reconstructed_hrfs[, vox]
  plot(hrf, type = "l", xlab = "Time index", ylab = "Amplitude",
       main = paste("Reconstructed HRF - voxel", vox), ...)
}
