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
#' @param durations The duration of events in the dataset. Default is 0
#'   (impulse events), as in [fmri_lm()].
#' @param drop_empty Whether to drop empty events from the model. Default is TRUE.
#' @param robust Must be FALSE. Robust weights are voxel-specific and cannot be
#'   computed on component scores, so any robust request is an error.
#' @param autocor The autocorrelation model used when `control` is `NULL`.
#'   'none' (OLS), 'ar1' or 'ar2'. AR models use one set of AR coefficients per
#'   run shared across all components, estimated so that the fit equals a dense
#'   fit of the reconstructed voxel data. 'auto' and 'arma' are not implemented
#'   and are an error. When `control` is supplied, `autocor` may be omitted; if
#'   it is also given it must name the same structure as `control`.
#' @param bootstrap Must be FALSE; bootstrapping is not implemented for latent
#'   datasets and TRUE is an error.
#' @param nboot Unused (see `bootstrap`).
#' @param control Optional [fmri_lm_control()] object. `NULL` (the default)
#'   fits `estimation_spec("joint")` with the noise model named by `autocor`,
#'   shared across components. Supply a control to choose the estimation scope
#'   and noise model explicitly, e.g.
#'   `fmri_lm_control(estimation = estimation_spec("runwise_meta"),
#'   noise = noise_spec("ar1", voxelwise = TRUE))` for a separate AR model per
#'   component (see Details).
#' @param ... Not used. Any argument passed here is an error, so a misspelled
#'   or unsupported argument is never silently dropped.
#'
#' @details
#' **Default (shared noise) model.** With `control = NULL` the model is fitted
#' jointly over runs and, for AR noise, uses one AR model per run shared by
#' all components. Because voxel data are `Y = Z %*% t(L)` (scores `Z`, loadings
#' `L`) this is exactly the fit a dense GLS of the reconstructed voxel data
#' would give, so voxel-level standard errors and statistics are available
#' with `recon = TRUE`.
#'
#' **Per-component noise.** Latent components can differ strongly in temporal
#' autocorrelation. A `control` with `estimation_spec("runwise_meta")` (and,
#' optionally, `noise_spec(..., voxelwise = TRUE)`) treats each component as
#' a separate response with its own noise model and run weights, i.e. the
#' same fit `fmri_lm()` gives for the score matrix as a `matrix_frame()`.
#' That model has no exact voxel-space equivalent: `coef(recon = TRUE)`
#' returns the loadings applied to the component estimates, but
#' `standard_error()` and `stats()` with `recon = TRUE` are an error because
#' the component noise models do not combine into one voxel covariance.
#' Under joint estimation only shared (run- or global-pooled) AR models are
#' supported.
#'
#' @return An object of class 'fmri_latent_lm' containing the regression model and dataset.
#'
#' @note This method is currently experimental.
#'
#' @seealso [latent_frame()], [fmri_lm()], [fmri_lm_control()]
#' @export
#'
#' @examples
#' 
#'
#' # Estimate the fMRI regression model using the latent dataset
#' #result <- fmri_latent_lm(formula = formula, block = block, dataset = dset,
#' #                          durations = NULL, drop_empty = TRUE, robust = FALSE)
#'
#' # Per-component AR(1) noise, run-wise estimation
#' #ctrl <- fmri_lm_control(estimation = estimation_spec("runwise_meta"),
#' #                        noise = noise_spec("ar1", voxelwise = TRUE))
#' #result <- fmri_latent_lm(formula, block = block, dataset = dset,
#' #                         control = ctrl)
#'
#' # Print the result
#' #print(result)
fmri_latent_lm <- function(formula, block, baseline_model=NULL, dataset, 
                    durations = 0, drop_empty=TRUE, robust=FALSE, 
                    autocor=c("none", "auto", "ar1", "ar2", "arma"), 
                    bootstrap=FALSE, nboot=1000, control = NULL,
                    ...) {
  dots <- list(...)
  if (length(dots)) {
    nms <- names(dots) %||% rep("", length(dots))
    nms[!nzchar(nms)] <- "<unnamed>"
    stop("Unused argument(s) to `fmri_latent_lm()`: ",
         paste(nms, collapse = ", "),
         ". Statistical options belong in `control = fmri_lm_control(...)`.",
         call. = FALSE)
  }

  autocor_missing <- missing(autocor)
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

  if (is.null(control)) {
    control <- fmri_lm_control(
      estimation = estimation_spec("joint"),
      noise = noise_spec(if (autocor == "none") "iid" else autocor)
    )
  } else {
    if (!inherits(control, "fmri_lm_control")) {
      stop("`control` must be an object created by `fmri_lm_control()`.",
           call. = FALSE)
    }
    ctrl_struct <- (control$noise %||% control$ar)$struct %||% "iid"
    want_struct <- if (autocor == "none") "iid" else autocor
    if (!autocor_missing && !identical(ctrl_struct, want_struct)) {
      stop("`autocor = '", autocor, "'` conflicts with the noise model in ",
           "`control` ('", ctrl_struct, "'). Specify the noise model in ",
           "`control` only.", call. = FALSE)
    }
    if (.fmri_lm_robust_enabled(control$robust)) {
      stop("Robust fitting is not supported for latent datasets: robust weights ",
           "are voxel-specific and cannot be computed on component scores. ",
           "Fit the reconstructed voxel data instead.", call. = FALSE)
    }
  }

  ret <- fmri_lm(formula, block, baseline_model=baseline_model, dataset=dataset,
                 durations=durations, drop_empty=drop_empty,
                 control = control, compute = compute_spec(voxel_chunks = 1L))

  ret$dataset <- dataset
  class(ret) <- c("fmri_latent_lm", class(ret))
  ret
}


#' Internal: refuse voxel-space inference that the latent fit cannot support
#'
#' Voxel standard errors from component results need one temporal covariance
#' shared by all components (joint fit, shared AR). Run-wise pooling and
#' per-component AR weight each component differently, so no voxel covariance
#' follows from the component results.
#' @keywords internal
#' @noRd
.latent_check_recon_inference <- function(x) {
  strategy <- x[["strategy"]] %||% attr(x, "strategy")
  if (!identical(strategy, "chunkwise") || is.null(x$result[["resid_cov"]])) {
    stop("Voxel-space standard errors and statistics (`recon = TRUE`) need a ",
         "noise model shared across components (the default ",
         "`estimation_spec(\"joint\")` fit). This fit used run-wise ",
         "estimation or per-component noise, which has no voxel-space ",
         "covariance. Use `recon = FALSE` for component-level inference.",
         call. = FALSE)
  }
  invisible(TRUE)
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
      contrast_prep = prepare_fmri_lm_contrasts(model),
      event_indices = event_indices, baseline_indices = baseline_indices
    ))
  }

  ret <- multiresponse_lm(form, data_env, contrast_objects, attr(tmats, "varnames"),
                          fcon = NULL, modmat = modmat)

  R <- as.matrix(ret$fit$residuals)
  rss <- colSums(R^2)
  rdf <- ret$fit$df.residual
  resvar <- rss / rdf
  sigma <- sqrt(resvar)

  unpack_chunkwise(list(ret), event_indices, baseline_indices) %>%
    purrr::list_modify(event_indices = event_indices,
                       baseline_indices = baseline_indices,
                       sigma = sigma,
                       residuals = resid(ret$fit),
                       df.residual = ret$fit$df.residual,
                       qr = qr.lm(ret$fit),
                       resid_cov = crossprod(R) / rdf)
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
.latent_ar_fit <- function(dset, modmat, basismat, cfg, contrast_prep,
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

  # The built-in fitters name contrasts by their own `name`.
  stats <- .external_fmri_lm_stats(
    betas, XtXinv, sigma_vec, sigma_vec^2, dfres, colnames(modmat),
    .prepared_contrast_weights(contrast_prep, names_from = "name"),
    ar_order = ar_order
  )
  bstats <- stats$betas
  conres <- stats$contrasts

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
  # Component residual covariance in the whitened domain (see
  # standard_error.fmri_latent_lm()).
  Rw <- result$Y_white - result$X_white %*% betas
  out$resid_cov <- crossprod(Rw) / dfres
  out$ar_coef <- phi
  out
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
    .latent_check_recon_inference(x)
    # Voxel coefficients are L B, so Var(L b_i) = C[i, i] * L Sigma L' with
    # C the unscaled coefficient covariance and Sigma the component residual
    # covariance (both in the fitted, whitened domain for AR fits).
    CR <- x$result[["resid_cov"]]
    cov.unscaled <- x$result[["cov.unscaled"]] %||% chol2inv(x$result[["qr"]]$qr)
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
