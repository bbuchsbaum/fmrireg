
#' @export
fmri_latent_lm <- function(formula, block, baseline_model=NULL, dataset, 
                    durations, drop_empty=TRUE, contrasts=NULL, robust=FALSE,
                    strategy="chunkwise", nchunks=1, ...) {
  
  
  assert_that(inherits(dataset, "latent_dataset"))
  result <- fmri_lm(formula, block, baseline_model=baseline_model, dataset, 
          durations, drop_empty=drop_empty, contrasts=contrasts, robust=robust, 
          strategy=strategy, nchunks=nchunks,...)
  result$dataset <- dataset
  class(result) <- c("fmri_latent_lm", class(result))
  result
}

#' @export
coef.fmri_latent_lm <- function(x, type=c("estimates", "contrasts"), recon=FALSE) {
  bvals <- coef.fmri_lm(x, type=type)
  
  
  if (recon) {
    lds <- x$dataset$lvec@loadings
    out <- t(as.matrix(bvals)) %*% t(lds)
    sp <- space(x$dataset$lvec@mask)
    SparseNeuroVec(as.matrix(out), neuroim2::add_dim(sp, nrow(out)), mask=x$dataset$lvec@mask)
  } else {
    bvals
  }
}

# stats.fmri_latent_lm <-
#   function(x,
#            type = c("estimates", "contrasts", "F"),
#            recon = TRUE) {
#     bvals <- stats.fmri_lm(x, type = type)
#     errs <- standard_error.fmri_lm(x, type = "estimates")
#     lds <- x$dataset$lvec@loadings
#     if (recon) {
#       if (type == "F") {
#         stop("cannot reconstruct an F-statistic")
#       }
#       out <- t(as.matrix(bvals)) %*% t(lds)
#       sp <- space(x$dataset$lvec@mask)
#       SparseNeuroVec(as.matrix(out),
#                      neuroim2::add_dim(sp, nrow(out)),
#                      mask = x$dataset$lvec@mask)
#     } else {
#       bvals
#     }
#   }
