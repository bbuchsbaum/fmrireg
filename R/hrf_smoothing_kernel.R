
#' compute an hrf smoothing kernel
#' 
#' A function to compute a temporal similarity matrix from a series of hemodynamic response functions
#' 
#' @param len the number of scans
#' @param TR the repetition time
#' @param basis the hemodynamic response function
#  @importFrom proxy simil
#' @export
hrf_smoothing_kernel <- function(len, TR=2, basis="spmg1") {
  buffer <- 6
  sframe <- sampling_frame(len+buffer,TR)
  onsets <- samples(sframe)
  dfx <- data.frame(onsets=samples(sframe), block=rep(1,len+buffer))
  em <- event_model(onsets ~ trialwise(basis="spmg1"), block=~ block, 
                    data=dfx, sampling_frame=sframe)
  
  dmat <- as.matrix(design_matrix(em))
  #m <- as.matrix(proxy::simil(t(dmat), "cosine"))
  m <- dmat %*% t(dmat)
  ret <- m[(buffer/2):(len+buffer/2),(buffer/2):(len+buffer/2)]
  ret
}


#'  
#' 
#' onsets <- sort(runif(25, 0, 200))
#' fac <- factor(sample(letters[1:4], length(onsets), replace=TRUE))
#' sframe <- sampling_frame(200, TR=1)
#' des <- data.frame(onsets=onsets, fac=fac, constant = factor(rep(1, length(fac))), block=rep(1, length(fac)))
#' emod <- event_model(onsets ~ hrf(constant) + hrf(fac), data=des, sampling_frame=sframe, block = ~ block)
#' dmat <- design_matrix(emod)
#' @noRd
#' @param dmat a design matrix
#' @param deriv compute derivative
#' @keywords internal
design_kernel <- function(dmat, deriv=FALSE) {
  dd <- rbind(rep(0, ncol(dmat)), apply(dmat, 2, diff))
}