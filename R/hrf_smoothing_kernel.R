
#' Compute an HRF smoothing kernel
#'
#' This function computes a temporal similarity matrix from a series of hemodynamic response functions.
#'
#' @param len The number of scans.
#' @param TR The repetition time (default is 2 seconds).
#' @param form the `trialwise` formula expression, see examples.
#' @importFrom proxy simil
#' @export
#' @examples
#' form <- onsets ~ trialwise(basis="gaussian")
#' sk <- hrf_smoothing_kernel(100, TR=1.5, form)
hrf_smoothing_kernel <- function(len, TR=2, form) {
  buffer <- 6
  sframe <- sampling_frame(len+buffer,TR)
  onsets <- samples(sframe)
  
  dfx <- data.frame(onsets=samples(sframe), block=rep(1,len+buffer))
  
  #em <- event_model(onsets ~ trialwise(basis={{basis}}), block=~ block, 
  #                  data=dfx, sampling_frame=sframe)
  
  em <- event_model(form, block=~ block, 
                                      data=dfx, sampling_frame=sframe)
  
  dmat <- as.matrix(design_matrix(em))
  #m <- as.matrix(proxy::simil(t(dmat), "cosine"))
  m <- dmat %*% t(dmat)
  ret <- m[(buffer/2):((len-1)+buffer/2),(buffer/2):((len-1)+buffer/2)]
  ret
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