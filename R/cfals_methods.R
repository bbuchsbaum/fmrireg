#' Construct an \code{fmrireg_cfals_fit} object
#'
#' Simple constructor used by [fmrireg_cfals()] to package the
#' results returned by the various CFALS engines.
#'
#' @param h_coeffs Matrix of HRF basis coefficients (d \eqn{\times} v).
#' @param beta_amps Matrix of condition amplitudes (k \eqn{\times} v).
#' @param method Character string indicating the estimation method.
#' @param lambdas Numeric vector of regularisation parameters.
#' @param call The matched call to the wrapper function.
#' @param hrf_basis HRF basis object used for the estimation.
#' @param design_info List with design metadata (d, k, n, v, fullXtX).
#' @param residuals Residual matrix from the projected data fit.
#' @param recon_hrf Matrix of reconstructed HRF shapes.
#' @param gof Numeric vector of goodness-of-fit statistics per voxel.
#' @return An object of class \code{fmrireg_cfals_fit}.
#' @keywords internal
fmrireg_cfals_fit <- function(h_coeffs, beta_amps, method, lambdas, call,
                              hrf_basis, design_info, residuals,
                              recon_hrf = NULL, gof = NULL) {
  out <- list(h_coeffs = h_coeffs,
              beta_amps = beta_amps,
              method_used = method,
              lambdas = lambdas,
              call = call,
              hrf_basis_used = hrf_basis,
              design_info = design_info,
              residuals = residuals,
              reconstructed_hrfs = recon_hrf,
              gof_per_voxel = gof)
  class(out) <- c("fmrireg_cfals_fit", "list")
  out
}

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
              lambdas = object$lambdas)
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
