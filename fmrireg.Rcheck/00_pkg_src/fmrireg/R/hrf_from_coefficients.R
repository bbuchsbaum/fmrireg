#' Combine HRF Basis with Coefficients
#'
#' Create a new HRF by linearly weighting the basis functions of an existing HRF.
#' Useful when coefficients have been estimated for an FIR/bspline/SPMG3 basis
#' and one wants a single functional HRF.
#'
#' @param hrf  An object of class `HRF`.
#' @param h    Numeric vector of length `nbasis(hrf)` giving the weights.
#' @param name Optional name for the resulting HRF.
#' @param ...  Reserved for future extensions.
#'
#' @return A new `HRF` object with `nbasis = 1`.
#' @export
hrf_from_coefficients <- function(hrf, h, ...) UseMethod("hrf_from_coefficients")

#' @export
#' @rdname hrf_from_coefficients
hrf_from_coefficients.HRF <- function(hrf, h, name = NULL, ...) {
  nb <- nbasis(hrf)
  if (length(h) != nb) {
    stop("length(h) must equal nbasis(hrf)")
  }
  weighted_fun <- function(t) {
    vals <- hrf(t)
    if (is.matrix(vals)) {
      drop(vals %*% as.numeric(h))
    } else {
      vals * h[1L]
    }
  }
  if (is.null(name)) {
    name <- paste0(attr(hrf, "name"), "_from_coef")
  }
  as_hrf(
    f      = weighted_fun,
    name   = name,
    nbasis = 1L,
    span   = attr(hrf, "span"),
    params = c(attr(hrf, "params"), list(coefficients = h))
  )
}
