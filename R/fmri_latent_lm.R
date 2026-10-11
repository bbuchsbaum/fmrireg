#' Fast fMRI Regression Model Estimation from a Latent Component Dataset
#'
#' This function estimates a regression model for fMRI data using a latent component dataset.
#' The dataset must be an `fmri_frame` whose feature space is a
#' `fmridataset::basis_space` (see [latent_frame()]), so that the assay holds
#' component scores and the space carries the spatial loadings.
#'
#' @param formula A formula specifying the regression model.
#' @param block A factor indicating the block structure of the data.
#' @param baseline_model An optional baseline model.
#' @param dataset A latent `fmri_frame` (see [latent_frame()]).
#' @param durations The duration of events in the dataset.
#' @param drop_empty Whether to drop empty events from the model. Default is TRUE.
#' @param robust Must be FALSE. Robust weights are voxel-specific and cannot be
#'   computed on component scores, so any robust request is an error.
#' @param autocor The autocorrelation model. 'none' (OLS), 'ar1' or 'ar2'. AR
#'   models use one set of AR coefficients per run shared across voxels,
#'   estimated so that the fit equals a dense fit of the reconstructed voxel
#'   data. 'auto' and 'arma' are not implemented and are an error.
#' @param bootstrap Must be FALSE; bootstrapping is not implemented for latent
#'   datasets and TRUE is an error.
#' @param nboot Unused (see `bootstrap`).
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
  assert_that(.dset_is_latent(dataset),
              msg = "'dataset' must be an fmri_frame over a basis_space (see latent_frame())")
  if (autocor %in% c("auto", "arma")) {
    stop("`fmri_latent_lm()` does not implement autocor = '", autocor,
         "'; choose 'none', 'ar1', or 'ar2'.", call. = FALSE)
  }
  if (!identical(robust, FALSE) && !identical(robust, "none")) {
    stop("Robust fitting is not supported for latent datasets: robust weights ",
         "are voxel-specific and cannot be computed on component scores. ",
         "Fit the reconstructed voxel data instead.", call. = FALSE)
  }
  if (!isFALSE(bootstrap)) {
    stop("`bootstrap` is not supported for latent datasets.", call. = FALSE)
  }
  control <- fmri_lm_control(
    estimation = estimation_spec("joint"),
    noise = noise_spec(if (autocor == "none") "iid" else autocor)
  )
  ret <- fmri_lm(formula, block, baseline_model=baseline_model, dataset=dataset,
                 durations=durations, drop_empty=drop_empty,
                 control = control, compute = compute_spec(voxel_chunks = 1L))

  ret$dataset <- dataset
  class(ret) <- c("fmri_latent_lm", class(ret))
  ret
}


#' Internal function for performing linear regression on latent datasets
#'
#' @keywords internal
#' @details This function is an internal helper function for fmri_latent_lm,
#'   which fits the model to the component scores of a latent dataset.
#'
#'   Because voxel data are \eqn{Y = Z L'} (scores \eqn{Z}, loadings \eqn{L}),
#'   any fit that is linear in \eqn{Y} can be done on \eqn{Z} and mapped back
#'   by \eqn{L}. That covers OLS and GLS with a temporal covariance shared
#'   across voxels (\code{noise_spec(pooling = "run")} or \code{"global"}). The
#'   shared AR parameters are estimated on \eqn{Z Q \Lambda^{1/2}}, where
#'   \eqn{L'L = Q \Lambda Q'}: those columns have the same time-by-time
#'   cross-product as the voxel residuals, so the pooled autocovariance -- and
#'   therefore phi -- equals the one a dense fit of the reconstructed data
#'   would estimate. Robust (voxel-specific, nonlinear) weights, voxelwise or
#'   parcel AR, and MA terms have no such reduction and are rejected.
#'
#' @param x An `fmri_frame` whose feature space is a `basis_space`.
#' @param model The fmri model object.
#' @param contrast_objects A list of contrast objects.
#' @param nchunks Ignored (components are fitted in one block).
#' @param cfg An `fmri_lm_control` object.
#' @param phi_fixed Ignored: a shared phi passed down by `fmri_lm_fit()` is
#'   estimated from component scores, not voxel residuals, so it is re-estimated
#'   here on the voxel-equivalent residuals.
#' @param sigma_fixed Ignored (robust fitting is unsupported).
#' @param ... Additional arguments (for compatibility).
#' @return A list containing the unpacked regression results.
#' @seealso fmri_latent_lm
#' @noRd
chunkwise_lm_latent <- function(x, model, contrast_objects, nchunks, cfg,
                                verbose = FALSE, use_fast_path = FALSE, progress = FALSE,
                                phi_fixed = NULL,
                                sigma_fixed = NULL, ...) {
  dset <- x
  .latent_check_supported_config(cfg)

  form <- get_formula(model)
  tmats <- term_matrices(model)
  data_env <- list2env(tmats)
  data_env[[".y"]] <- rep(0, nrow(tmats[[1]]))
  modmat <- model.matrix(as.formula(form), data_env)

  basismat <- as.matrix(.dset_data_matrix(dset))
  data_env[[".y"]] <- basismat

  event_indices <- attr(tmats, "event_term_indices")
  baseline_indices <- attr(tmats, "baseline_term_indices")

  if (!is.null(cfg) && .temporal_noise_enabled(cfg$ar)) {
    return(.latent_ar_fit(
      dset = dset, modmat = modmat, basismat = basismat, cfg = cfg,
      contrast_objects = contrast_objects,
      event_indices = event_indices, baseline_indices = baseline_indices
    ))
  }

  ret <- multiresponse_lm(form, data_env, contrast_objects, attr(tmats, "varnames"),
                          fcon = NULL, modmat = modmat)

  rss <- colSums(as.matrix(ret$fit$residuals^2))
  rdf <- ret$fit$df.residual
  resvar <- rss / rdf
  sigma <- sqrt(resvar)

  unpack_chunkwise(list(ret), event_indices, baseline_indices) %>%
    purrr::list_modify(event_indices = event_indices,
                       baseline_indices = baseline_indices,
                       sigma = sigma,
                       residuals = resid(ret$fit),
                       df.residual = ret$fit$df.residual,
                       qr = qr.lm(ret$fit))
}


#' Internal: reject fitting options with no exact latent-space reduction
#' @keywords internal
#' @noRd
.latent_check_supported_config <- function(cfg) {
  if (is.null(cfg)) return(invisible(TRUE))
  if (.fmri_lm_robust_enabled(cfg$robust)) {
    stop("Robust fitting is not supported for latent datasets: robust weights ",
         "are voxel-specific and cannot be computed on component scores. ",
         "Fit the reconstructed voxel data instead.", call. = FALSE)
  }
  ar <- cfg$ar
  if (!is.null(ar) && .temporal_noise_enabled(ar)) {
    if (as.integer(ar$q %||% 0L) > 0L) {
      stop("MA (ARMA) noise models are not supported for latent datasets; ",
           "use an AR structure ('ar1', 'ar2', 'arp').", call. = FALSE)
    }
    if (isTRUE(ar$voxelwise) || isTRUE(ar$by_cluster) ||
        identical(ar$pooling, "parcel")) {
      stop("Voxelwise or parcel-wise AR models are not supported for latent ",
           "datasets; use noise_spec(pooling = 'run') or 'global' (one AR ",
           "model shared across voxels).", call. = FALSE)
    }
  }
  if (isTRUE(cfg$volume_weights$enabled)) {
    stop("Volume weighting is not supported for latent datasets.", call. = FALSE)
  }
  invisible(TRUE)
}


#' Internal: shared-AR GLS fit on latent component scores
#' @keywords internal
#' @noRd
.latent_ar_fit <- function(dset, modmat, basismat, cfg, contrast_objects,
                           event_indices, baseline_indices) {
  run_indices <- lapply(.dset_run_chunks(dset), `[[`, "row_ind")

  # Voxel-equivalent columns for AR estimation (see chunkwise_lm_latent).
  L <- as.matrix(.dset_loadings(dset))
  eg <- eigen(crossprod(L), symmetric = TRUE)
  root <- eg$vectors %*% diag(sqrt(pmax(eg$values, 0)), nrow = length(eg$values))
  Y_equiv <- basismat %*% root

  cfg_est <- cfg
  cfg_est$ar$phi <- NULL
  cfg_est$ar$phi_fixed <- NULL
  est <- solve_integrated_glm(X = modmat, Y = Y_equiv, config = cfg_est,
                              run_indices = run_indices)
  phi <- est$ar_coef %||% est$phi_hat
  if (is.null(phi) || !length(unlist(phi))) {
    stop("Failed to estimate AR parameters for the latent dataset.", call. = FALSE)
  }

  # GLS is linear in the response once phi is fixed: fit the scores directly.
  cfg_fit <- cfg
  cfg_fit$ar$phi_fixed <- phi
  result <- solve_integrated_glm(X = modmat, Y = basismat, config = cfg_fit,
                                 run_indices = run_indices)

  betas <- as.matrix(result$betas %||% result$coefficients)
  XtXinv <- result$XtXinv
  # [[ ]] not $: `result$sigma` partial-matches `sigma2` (a variance).
  sigma_vec <- as.numeric(result[["sigma"]] %||% sqrt(result[["sigma2"]]))
  if (length(sigma_vec) == 1L) sigma_vec <- rep(sigma_vec, ncol(betas))
  dfres <- result[["dfres"]] %||% result[["df_residual"]]
  ar_order <- result$ar_order %||% .target_ar_order(.normalize_ar_options(cfg$ar)) %||% 0L

  simple_con <- Filter(function(x) inherits(x, "contrast") && !inherits(x, "Fcontrast"),
                       contrast_objects)
  f_con <- Filter(function(x) inherits(x, "Fcontrast"), contrast_objects)
  weights_of <- function(cons) {
    w <- lapply(cons, function(x) {
      ww <- x$weights
      attr(ww, "colind") <- attr(x, "colind")
      ww
    })
    if (length(cons)) names(w) <- vapply(cons, `[[`, character(1), "name")
    w
  }

  bstats <- beta_stats_matrix(betas, XtXinv, sigma_vec, dfres, colnames(modmat),
                              ar_order = ar_order)
  conres <- fit_lm_contrasts_fast(betas, sigma_vec^2, XtXinv,
                                  weights_of(simple_con), weights_of(f_con),
                                  dfres, ar_order = ar_order)

  out <- unpack_chunkwise(
    list(list(bstats = bstats, contrasts = conres,
              event_indices = event_indices, baseline_indices = baseline_indices)),
    event_indices, baseline_indices
  )
  out$event_indices <- event_indices
  out$baseline_indices <- baseline_indices
  out$sigma <- sigma_vec
  out$df.residual <- dfres
  out$cov.unscaled <- XtXinv
  out$ar_coef <- phi
  out
}


#' @keywords internal
#' @noRd
tibble_to_neurovec <- function(dset, tab, mask) {
  sp <- .fmri_dataset_mask_space(dset, "latent coefficient reconstruction")$space
  neuroim2::SparseNeuroVec(as.matrix(tab), neuroim2::add_dim(sp, nrow(tab)), mask=mask)
}

#' @export
coef.fmri_latent_lm <- function(object, type=c("estimates", "contrasts", "betas"), recon=FALSE, comp=0, ...) {
  type <- match.arg(type)
  # coef.fmri_lm uses "betas"; accept the latent-friendly alias "estimates".
  coef_type <- if (identical(type, "estimates")) "betas" else type
  bvals <- coef.fmri_lm(object, type = coef_type)
  if (recon) {
    lds <- .dset_loadings(object$dataset)
    comp <- if (length(comp) == 1 && comp == 0) {
      seq(1, ncol(lds)) 
    } else {
      assertthat::assert_that(all(comp > 0) && all(comp <= ncol(lds)))
      comp
    }
    
    lds <- lds[,comp,drop=FALSE]
    bmat <- as.matrix(bvals)

    # coef.fmri_lm() is components x terms for every `type` (#227), so the
    # reconstruction is loadings (voxels x comp) %*% coefficients (comp x terms).
    if (nrow(bmat) < max(comp)) {
      stop("Cannot align latent coefficients with component loadings for reconstruction",
           call. = FALSE)
    }
    out <- as.matrix(lds %*% bmat[comp, , drop = FALSE])
    if (!is.null(colnames(bmat))) {
      colnames(out) <- colnames(bmat)
    }
    tibble::as_tibble(out)
  } else {
    bvals
  }
}

#' @method standard_error fmri_latent_lm
#' @rdname standard_error
#' @export
#' @importFrom Matrix rowSums t
standard_error.fmri_latent_lm <- function(x, type = c("betas", "contrasts"), recon = FALSE, ...) {
  type <- .accessor_family(type[1], c("betas", "contrasts"), "standard_error")

  if (!recon) {
    standard_error.fmri_lm(x, type)
  } else {
    R <- x$result$residuals
    CR <- cov(R) * (nrow(R) - 1)/x$result$df.residual
    
    Qr <- x$result$qr
    cov.unscaled <- chol2inv(Qr$qr)
    lds <- .dset_loadings(x$dataset)
  
    if (type == "betas") {
      ret <- do.call(cbind, lapply(x$result$event_indices, function(i) {
        sqrt(rowSums((lds %*% (CR * cov.unscaled[i,i])) * lds))
      }))
      # one name per design column (conditions() is per condition, not per basis)
      colnames(ret) <- colnames(design_matrix(x$model))[unlist(x$result$event_indices)]
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
stats.fmri_latent_lm <- function(x, type = c("betas", "contrasts"), recon = FALSE, ...) {
    type <- .accessor_family(type[1], c("betas", "contrasts"), "stats",
                             deprecate_estimates = TRUE)
    if (!recon) {
      stats.fmri_lm(x, type)
    } else {
      if (type == "betas") {
        bvals <- coef(x, type = "betas", recon = TRUE)
        errs <- standard_error(x, type = "betas", recon = TRUE)
        tibble::as_tibble(bvals/errs)
      } else {
        cvals <- coef(x, type="contrasts", recon=TRUE)
        errs <- standard_error(x, type = "contrasts", recon=TRUE)
        tibble::as_tibble(cvals/errs)
      }
    }
}
