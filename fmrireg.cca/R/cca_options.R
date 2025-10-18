#' Options for the CCA engine
#'
#' @param mode Character. "scale-only" (default) or "full" for orientation steering.
#' @param spatial Character. "3d" (default) or "2d".
#' @param fwhm Numeric. Base Gaussian FWHM in mm (default 5).
#' @param alpha Numeric or NULL. Alpha for temporal alpha-mix; NULL = auto.
#' @param scales Numeric vector of scale multipliers (default c(1.0)).
#' @param ar1_prewhite Logical. Apply AR(1) prewhitening before CCA (default TRUE).
#' @param shrink Numeric. Diagonal loading for covariance shrinkage (default 0.02).
#' @param constrain Logical. Enforce non-negativity constraints (default TRUE).
#' @return A named list of options suitable to pass as `cca = cca_options(...)`.
#' @export
cca_options <- function(mode = c("scale-only", "full"), spatial = c("3d", "2d"),
                        fwhm = 5, alpha = NULL, scales = c(1.0),
                        ar1_prewhite = TRUE, shrink = 0.02, constrain = TRUE) {
  list(
    mode = match.arg(mode),
    spatial = match.arg(spatial),
    fwhm = fwhm,
    alpha = alpha,
    scales = scales,
    ar1_prewhite = isTRUE(ar1_prewhite),
    shrink = shrink,
    constrain = isTRUE(constrain)
  )
}
