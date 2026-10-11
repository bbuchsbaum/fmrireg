
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

#' Design kernel for a given design matrix
#'
#' This internal function calculates a design kernel for the given design matrix (`dmat`). It can also compute the derivative of the design kernel if the `deriv` parameter is set to TRUE.
#' @noRd
#' @param dmat A design matrix.
#' @param deriv A logical value indicating whether to compute the derivative (default is FALSE).
#' @keywords internal
#' @examples
#' onsets <- sort(runif(25, 0, 200))
#' fac <- factor(sample(letters[1:4], length(onsets), replace=TRUE))
#' sframe <- sampling_frame(200, TR=1)
#' des <- data.frame(onsets=onsets, fac=fac, constant = factor(rep(1, length(fac))), block=rep(1, length(fac)))
#' emod <- event_model(onsets ~ hrf(constant) + hrf(fac), data=des, sampling_frame=sframe, block = ~ block)
#' dmat <- design_matrix(emod)
design_kernel <- function(dmat, deriv=FALSE) {
  dd <- rbind(rep(0, ncol(dmat)), apply(dmat, 2, diff))
}