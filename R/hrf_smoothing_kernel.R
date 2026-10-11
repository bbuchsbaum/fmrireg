
#' Compute an HRF smoothing kernel
#'
#' This function computes a temporal similarity matrix from a series of hemodynamic response functions.
#'
#' @param len The number of scans.
#' @param TR The repetition time (default is 2 seconds).
#' @param form the `trialwise` formula expression, see examples.
#' @param buffer_scans The number of scans to buffer before and after the event.
#' @param normalise Whether to normalise the kernel.
#' @param method "gram" (inner products of the design rows) or "cosine"
#'   (the same, scaled to unit diagonal). With `normalise = TRUE` both give
#'   `K[i, j] / sqrt(K[i, i] * K[j, j])`, a symmetric matrix with unit diagonal.
#' @export
#' @examples
#' form <- onset ~ trialwise(basis="gaussian")
#' sk <- hrf_smoothing_kernel(100, TR=1.5, form)
#' @return a smoothing matrix
hrf_smoothing_kernel <- function(len, TR = 2,
                                 form  = onset ~ trialwise(),
                                 buffer_scans = 3L,
                                 normalise = TRUE,
                                 method = c("gram", "cosine")) {

  method <- match.arg(method)
  n_buf  <- as.integer(buffer_scans)

  sf   <- sampling_frame(len + 2 * n_buf, TR)
  dfx  <- data.frame(onset = fmrihrf::samples(sf), block = 1L)
  em   <- event_model(form, data = dfx, block = ~ block, sampling_frame = sf)
  X    <- as.matrix(design_matrix(em))

  # Similarity between scans (rows of X). "cosine" is the scan-by-scan cosine,
  # i.e. the Gram matrix scaled symmetrically to unit diagonal.
  K <- tcrossprod(X)
  if (method == "cosine" || normalise) {
    K <- .unit_diagonal(K)
  }

  keep <- (n_buf + 1):(n_buf + len)
  K[keep, keep, drop = FALSE]
}

#' Scale a Gram matrix to unit diagonal, symmetrically
#'
#' `K[i, j] / sqrt(K[i, i] K[j, j])`. A scan with no signal (zero row of the
#' design, so `K[i, i] == 0`) is similar only to itself: its row and column are
#' 0 with a 1 on the diagonal, rather than NaN.
#' @keywords internal
#' @noRd
.unit_diagonal <- function(K) {
  d <- diag(K)
  s <- ifelse(d > 0, 1 / sqrt(d), 0)
  K <- K * outer(s, s)
  K <- (K + t(K)) / 2
  diag(K) <- 1
  K
}
