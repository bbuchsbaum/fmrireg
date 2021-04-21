
#' @export
#' @inheritParams fmri_lm
fmri_latent_lm <- function(formula, block, baseline_model=NULL, dataset, 
                    durations, drop_empty=TRUE, robust=FALSE, 
                    autocor=c("auto", "ar1", "ar2", "arma", "none"), 
                    ...) {
  
  autocor <- match.arg(autocor)
  assert_that(inherits(dataset, "latent_dataset"))
  
  model <- create_fmri_model(formula, block, baseline_model,dataset, durations, drop_empty)
  ret <- fmri_lm_fit(model, dataset, strategy="chunkwise", robust=robust, 
                     contrasts=contrasts, nchunks=1, autocor=autocor)

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
  
  rss <- colSums(as.matrix(ret$fit$residuals^2))
  rdf <- ret$fit$df.residual
  resvar <- rss/rdf
  sigma <- sqrt(resvar)
  
  wrap_chunked_lm_results(list(ret), event_indices) %>% purrr::list_modify(event_indices=event_indices,
                                                                           baseline_indices=baseline_indices,
                                                                           sigma=sigma,
                                                                           residuals=resid(ret$fit),
                                                                           df.residual=ret$fit$df.residual,
                                                                           qr=stats:::qr.lm(ret$fit))
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

#' @export
standard_error.fmri_latent_lm <- function(x, type=c("estimates", "contrasts"), recon=FALSE) {
  type <- match.arg(type)
  if (!recon) {
    standard_error.fmri_lm(x,type)
  } else {
    R <- x$result$residuals
    CR <- cov(R) * (nrow(R) - 1)/x$result$df.residual

    Qr <- x$result$qr
    cov.unscaled <- chol2inv(Qr$qr)
    lds <- x$dataset$lvec@loadings
  
    if (type == "estimates") {
      ret <- do.call(cbind, lapply(x$result$event_indices, function(i) {
        sqrt(rowSums((lds %*% (CR * cov.unscaled[i,i])) * lds))
      }))
      colnames(ret) <- conditions(x$model$event_model)
      as_tibble(ret,.name_repair="check_unique")
    } else {
      ret <- do.call(cbind, lapply(x$result$conmats, function(cmat) {
        csc <- t(cmat) %*% cov.unscaled %*% cmat
        sqrt(rowSums((lds %*% (CR * csc[1,1])) * lds))
      }))
      colnames(ret) <- names(x$bcons)
      as_tibble(ret,.name_repair="check_unique")
    }
  }
  
  ret
  
}

#' @export
stats.fmri_latent_lm <- function(x,type = c("estimates", "contrasts"),recon = FALSE) {
    if (!recon) {
      stats.fmri_lm(x,type)
    } else {
      if (type == "estimates") {
        bvals <- coef(x, type = "estimates", recon=TRUE)
        errs <- standard_error(x, type = "estimates", recon=TRUE)
        as_tibble(bvals/errs)
      } else {
        cvals <- coef(x, type="contrasts", recon=TRUE)
        errs <- standard_error(x, type = "contrasts", recon=TRUE)
        as_tibble(cvals/errs)
      }
    }
}
