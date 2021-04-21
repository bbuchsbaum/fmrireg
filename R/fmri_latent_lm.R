
#' @export
#' @inheritParams fmri_lm
fmri_latent_lm <- function(formula, block, baseline_model=NULL, dataset, 
                    durations, drop_empty=TRUE, robust=FALSE, 
                    autocor=c("auto", "ar1", "ar2", "arma", "none"), 
                    ...) {
  
  autocor <- match.arg(autocor)
  ### might be easier here to rewrite specifically for latent_lm, handling bootstrap and potential prewhitening.
  assert_that(inherits(dataset, "latent_dataset"))
  
  model <- create_fmri_model(formula, block, baseline_model,dataset, durations, drop_empty)
  ret <- fmri_lm_fit(model, dataset, strategy="chunkwise", robust=robust, 
                     contrasts=contrasts, nchunks=1, autocor=autocor)
  
  #result <- fmri_lm(formula, block, baseline_model=baseline_model, dataset, 
  #        durations, drop_empty=drop_empty, robust=robust, 
  #        strategy="chunkwise", nchunks=1, autocor=autocor, ...)
  #result$dataset <- dataset
  ret$dataset <- dataset
  class(ret) <- c("fmri_latent_lm", class(ret))
  ret
}

#' @keywords internal
chunkwise_lm.latent_dataset <- function(dset, model, conlist, nchunks, robust=FALSE, verbose=FALSE, 
                                        autocor=c("auto", "ar1", "ar2", "arma", "none")) {
  autocor <- match.arg(autocor)
  form <- get_formula(model)
  tmats <- term_matrices(model)
  data_env <- list2env(tmats)
  data_env[[".y"]] <- rep(0, nrow(tmats[[1]]))
  modmat <- model.matrix(as.formula(form), data_env)
  
  basismat <- get_data(dset)
  
  wmat <- if (autocor != "none") {
    message("whitening components")
    auto_whiten(basismat, modmat, autocor)
  } else {
    basismat
  }
  
  lmfun <- if (robust) multiresponse_rlm else multiresponse_lm
  data_env[[".y"]] <- as.matrix(wmat)
  
  ret <- lmfun(form, data_env, conlist, attr(tmats,"varnames"), fcon=NULL, modmat=modmat)
  event_indices=attr(tmats, "event_term_indices")
  baseline_indices=attr(tmats, "baseline_term_indices")
  wrap_chunked_lm_results(list(ret), event_indices)
}

#' @export
coef.fmri_latent_lm <- function(x, type=c("estimates", "contrasts"), recon=FALSE) {
  bvals <- coef.fmri_lm(x, type=type)
  
  if (recon) {
    lds <- x$dataset$lvec@loadings
    out <- t(as.matrix(bvals)) %*% t(lds)
    out <- as.matrix(t(out))
    as_tibble(out)
    #sp <- space(get_mask(x$dataset))
    #SparseNeuroVec(as.matrix(out), neuroim2::add_dim(sp, nrow(out)), mask=x$dataset$lvec@mask)
  } else {
    bvals
  }
}

stats.fmri_latent_lm <-
  function(x,
           type = c("estimates", "contrasts", "F"),
           recon = TRUE) {
    bvals <- stats.fmri_lm(x, type = type)
    #errs <- standard_error.fmri_lm(x, type = "estimates")
    lds <- x$dataset$lvec@loadings
    if (recon) {
      if (type == "F") {
        stop("cannot reconstruct an F-statistic")
      }
      out <- as.matrix(t(t(as.matrix(bvals)) %*% t(lds)))
      #sp <- space(x$dataset$lvec@mask)
      #SparseNeuroVec(as.matrix(out),
      #               neuroim2::add_dim(sp, nrow(out)),
      #               mask = x$dataset$lvec@mask)
      as_tibble(out)
    } else {
      bvals
    }
  }
