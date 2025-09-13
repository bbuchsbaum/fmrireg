#' Fast fMRI Regression Model Estimation from a Latent Component Dataset
#'
#' This function estimates a regression model for fMRI data using a latent component dataset.
#' The dataset must be of type `latent_dataset`, which itself requires a `LatentNeuroVec` input.
#'
#' @param formula A formula specifying the regression model.
#' @param block A factor indicating the block structure of the data.
#' @param baseline_model An optional baseline model.
#' @param dataset A dataset of class 'latent_dataset'.
#' @param durations The duration of events in the dataset.
#' @param drop_empty Whether to drop empty events from the model. Default is TRUE.
#' @param robust Whether to use robust regression methods. Default is FALSE.
#' @param autocor The autocorrelation correction method to use on components.
#'   One of 'none', 'auto', 'ar1', 'ar2', or 'arma'. Default is 'none'.
#' @param bootstrap Whether to compute bootstrapped parameter estimates. Default is FALSE.
#' @param nboot The number of bootstrap iterations. Default is 1000.
#' @param ... Additional arguments.
#'
#' @return An object of class 'fmri_latent_lm' containing the regression model and dataset.
#'
#' @note This method is currently experimental.
#'
#' @export
#'
#' @examples
#' 
#'
#' # Estimate the fMRI regression model using the latent dataset
#' #result <- fmri_latent_lm(formula = formula, block = block, dataset = dset,
#' #                          durations = NULL, drop_empty = TRUE, robust = FALSE)
#'
#' # Print the result
#' #print(result)
fmri_latent_lm <- function(formula, block, baseline_model=NULL, dataset, 
                    durations, drop_empty=TRUE, robust=FALSE, 
                    autocor=c("none", "auto", "ar1", "ar2", "arma"), 
                    bootstrap=FALSE, nboot=1000,
                    ...) {
  
  autocor <- match.arg(autocor)
  assert_that(inherits(dataset, "latent_dataset"))
  
  model <- create_fmri_model(formula, block, baseline_model,dataset, drop_empty=drop_empty)
  # The new fmri_lm function handles all these options
  ret <- fmri_lm(formula, block, baseline_model=baseline_model, dataset=dataset,
                 durations=durations, drop_empty=drop_empty, 
                 strategy="chunkwise", nchunks=1,
                 robust=robust)

  ret$dataset <- dataset
  class(ret) <- c("fmri_latent_lm", class(ret))
  ret
}


#' Internal function for performing chunkwise linear regression on latent datasets
#'
#' @keywords internal
#' @details This function is an internal helper function for fmri_latent_lm, which performs
#'          chunkwise linear regression on latent datasets. The function handles different
#'          autocorrelation options, as well as robust regression and bootstrapping.
#'
#' @param dset A latent dataset object.
#' @param model The fmri model object.
#' @param contrast_objects A list of contrast objects.
#' @param nchunks The number of chunks to use for the regression.
#' @param cfg Configuration object containing robust and AR options.
#' @param progress Logical, if TRUE, show progress bar.
#' @param phi_fixed Fixed AR coefficients (optional).
#' @param sigma_fixed Fixed sigma values (optional).
#' @param verbose Logical, if TRUE, additional output will be printed.
#' @param ... Additional arguments (for compatibility).
#' @return A list containing the results of the chunkwise linear regression.
#' @seealso fmri_latent_lm
#' @noRd
chunkwise_lm.latent_dataset <- function(dset, model, contrast_objects, nchunks, cfg,
                                        verbose = FALSE, use_fast_path = FALSE, progress = FALSE,
                                        phi_fixed = NULL,
                                        sigma_fixed = NULL) {
  # Extract options from cfg if provided
  robust <- if (!is.null(cfg) && !is.null(cfg$robust)) {
    r <- as.logical(cfg$robust)[1]
    if (is.na(r)) FALSE else r
  } else {
    FALSE
  }
  autocor <- if (!is.null(cfg) && !is.null(cfg$ar$cor_struct)) {
    cfg$ar$cor_struct
  } else {
    "none"
  }
  bootstrap <- FALSE  # Not currently supported through cfg
  nboot <- 1000
  boot_rows <- FALSE
  
  # Map conlist to maintain backward compatibility
  conlist <- contrast_objects
  form <- get_formula(model)
  tmats <- term_matrices(model)
  data_env <- list2env(tmats)
  data_env[[".y"]] <- rep(0, nrow(tmats[[1]]))
  modmat <- model.matrix(as.formula(form), data_env)
  
  basismat <- fmridataset::get_data(dset)

  #wmat <- if (autocor != "none") {
  #  message("whitening components")
  #  auto_whiten(basismat, modmat, autocor)
  #} else {
  #  basismat
  #}
  
  #data_env[[".y"]] <- as.matrix(wmat)
  data_env[[".y"]] <- as.matrix(basismat)
  
  event_indices=attr(tmats, "event_term_indices")
  baseline_indices=attr(tmats, "baseline_term_indices")
  
  if (bootstrap) {
    lmfun <- multiresponse_bootstrap_lm
    ret <- lmfun(form, data_env, conlist, attr(tmats,"varnames"), fcon=NULL, 
                 modmat=modmat, 
                 nboot=nboot, 
                 boot_rows=boot_rows,
                 event_indices=event_indices)
  } else if (autocor != "none") {
    
    ## need to split by run
    lmfun <- multiresponse_arma
   
    ret <- lmfun(form, data_env, conlist, attr(tmats,"varnames"), fcon=NULL, 
                 modmat=modmat, blockids=model$event_model$sampling_frame$blockids, autocor=autocor)
    unpack_chunkwise(ret, event_indices, baseline_indices) %>% purrr::list_modify(event_indices=event_indices,
                                                                                  baseline_indices=baseline_indices)
    
  } else {
    lmfun <- if (robust) multiresponse_rlm else multiresponse_lm
    ret <- lmfun(form, data_env, conlist, attr(tmats,"varnames"), fcon=NULL, 
                 modmat=modmat)
    
    rss <- colSums(as.matrix(ret$fit$residuals^2))
    rdf <- ret$fit$df.residual
    resvar <- rss/rdf
    sigma <- sqrt(resvar)
    
    unpack_chunkwise(list(ret), event_indices, baseline_indices) %>% purrr::list_modify(event_indices=event_indices,
                                                                                        baseline_indices=baseline_indices,
                                                                                        sigma=sigma,
                                                                                        residuals=resid(ret$fit),
                                                                                        df.residual=ret$fit$df.residual,
                                                                                        qr=qr.lm(ret$fit))
    
  }

}

  
#' @keywords internal
#' @noRd
tibble_to_neurovec <- function(dset, tab, mask) {
  sp <- neuroim2::space(fmridataset::get_mask(dset))
  neuroim2::SparseNeuroVec(as.matrix(tab), neuroim2::add_dim(sp, nrow(tab)), mask=mask)
}

#' @export
coef.fmri_latent_lm <- function(object, type=c("estimates", "contrasts"), recon=FALSE, comp=0, ...) {
  bvals <- coef.fmri_lm(object, type=type)
  if (recon) {
    # Extract LatentNeuroVec from dataset
    lvec <- if (!is.null(object$dataset$lvec)) {
      object$dataset$lvec
    } else if (!is.null(object$dataset$backend) && !is.null(object$dataset$backend$data)) {
      object$dataset$backend$data[[1]]
    } else {
      stop("Cannot find LatentNeuroVec in latent_dataset")
    }
    lds <- lvec@loadings
    comp <- if (length(comp) == 1 && comp == 0) {
      seq(1, ncol(lds)) 
    } else {
      assertthat::assert_that(all(comp > 0) && all(comp <= ncol(lds)))
      comp
    }
    
    lds <- lds[,comp,drop=FALSE]
    
    out <- t(as.matrix(bvals)[comp,]) %*% t(lds)
    out <- as.matrix(t(out))
    tibble::as_tibble(out)
  } else {
    bvals
  }
}

#' @method standard_error fmri_latent_lm
#' @rdname standard_error
#' @export
#' @importFrom Matrix rowSums t
standard_error.fmri_latent_lm <- function(x, type=c("estimates", "contrasts"), recon=FALSE,...) {
  type <- match.arg(type)
 
  if (!recon) {
    standard_error.fmri_lm(x, type)
  } else {
    R <- x$result$residuals
    CR <- cov(R) * (nrow(R) - 1)/x$result$df.residual
    
    Qr <- x$result$qr
    cov.unscaled <- chol2inv(Qr$qr)
    # Extract LatentNeuroVec from dataset
    lvec <- if (!is.null(x$dataset$lvec)) {
      x$dataset$lvec
    } else if (!is.null(x$dataset$backend) && !is.null(x$dataset$backend$data)) {
      x$dataset$backend$data[[1]]
    } else {
      stop("Cannot find LatentNeuroVec in latent_dataset")
    }
    lds <- lvec@loadings
  
    if (type == "estimates") {
      ret <- do.call(cbind, lapply(x$result$event_indices, function(i) {
        sqrt(rowSums((lds %*% (CR * cov.unscaled[i,i])) * lds))
      }))
      colnames(ret) <- conditions(x$model$event_model)
      tibble::as_tibble(ret, .name_repair="check_unique")
    } else {
      # Handle contrasts case
      contrasts_df <- x$result$contrasts %>% dplyr::filter(type == "contrast")
      
      if (nrow(contrasts_df) == 0) {
        return(tibble())  # Return empty tibble if no contrasts
      }
      
      ret <- do.call(cbind, lapply(seq_len(nrow(contrasts_df)), function(i) {
        cmat <- as.matrix(contrasts_df$conmat[[i]])  # Get contrast matrix directly from tibble
        cind <- contrasts_df$colind[[i]]  # Get column indices from tibble
        
        csc <- t(cmat) %*% cov.unscaled[cind, cind, drop=FALSE] %*% cmat
        sqrt(rowSums(as.matrix((lds %*% (CR * csc[1,1])) * lds)))
      }))
      
      colnames(ret) <- contrasts_df$name  # Use names from the contrasts tibble
      tibble::as_tibble(ret, .name_repair="check_unique")
    }
  }
}

#' @export
stats.fmri_latent_lm <- function(x,type = c("estimates", "contrasts"), recon = FALSE, ...) {
    if (!recon) {
      stats.fmri_lm(x,type)
    } else {
      if (type == "estimates") {
        bvals <- coef(x, type = "estimates", recon=TRUE)
        errs <- standard_error(x, type = "estimates", recon=TRUE)
        tibble::as_tibble(bvals/errs)
      } else {
        cvals <- coef(x, type="contrasts", recon=TRUE)
        errs <- standard_error(x, type = "contrasts", recon=TRUE)
        tibble::as_tibble(cvals/errs)
      }
    }
}
