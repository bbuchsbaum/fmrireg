hrf_smoothing_kernel <- function(len, TR=2) {
  sframe <- sampling_frame(len,TR)
  onsets <- samples(sframe)
  dfx <- data.frame(onsets=samples(sframe), block=rep(1,len))
  em <- event_model(onsets ~ trialwise(basis="spmg1"), block=~ block, 
                    data=dfx, sampling_frame=sframe)
  
  dmat <- as.matrix(design_matrix(em))
  cov(dmat)
}