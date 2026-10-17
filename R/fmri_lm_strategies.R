# Shared Strategy Components for fMRI Linear Models
# Common functions used by both runwise and chunkwise strategies

#' Extract Censor Vector from Dataset
#'
#' @description
#' Extracts and formats censor information from a dataset object.
#' Handles both binary (0/1) vectors and integer index vectors.
#'
#' @param dataset fmri_frame object (may carry a censor column)
#' @param run_num Optional run number for run-specific extraction
#' @param n_time Optional total number of timepoints (for validation)
#' @return Integer vector of 1-based timepoint indices to censor, or NULL if no censoring
#' @keywords internal
#' @noRd
extract_censor_from_dataset <- function(dataset, run_num = NULL, n_time = NULL) {
  # Only a frame can carry a censor column. Anything else has no censoring to
  # report rather than being an error: this is consulted on every fit now that
  # an unset censor resolves against the dataset, so a non-frame must not blow
  # up inside fmridataset's accessors.
  if (is.null(dataset) || !.is_fmri_frame(dataset)) {
    return(NULL)
  }

  censor <- .dset_censor(dataset)
  if (is.null(censor)) {
    return(NULL)
  }

  if (!is.null(run_num)) {
    censor <- censor[.dset_run_rows(dataset)[[run_num]]]
  }
  .censor_indices(censor)
}

#' Censor specification as 1-based indices
#'
#' A logical or 0/1 mask becomes the indices of its censored (TRUE/1)
#' entries; anything else is taken to be indices already.
#' @return Integer indices, or NULL when nothing is censored.
#' @keywords internal
#' @noRd
.censor_indices <- function(censor) {
  if (is.logical(censor)) {
    censor <- which(censor)
  } else if (is.numeric(censor) && all(censor %in% c(0, 1))) {
    censor <- which(censor == 1)
  }
  if (length(censor) == 0) {
    return(NULL)
  }
  as.integer(censor)
}

#' Resolve Censor Configuration
#'
#' @description
#' Resolves censor specification from config and/or dataset.
#'
#' An unset `censor` is not a decision to ignore censoring. When the dataset
#' carries a censor column the frame is the single source of truth, so an
#' unset spec resolves against it. Requiring `ar_options$censor = "auto"` to
#' activate the frame's own censoring meant a `matrix_frame(censor = )` was
#' silently discarded by every fit. Pass `censor = "none"` to opt out.
#'
#' @param cfg fmri_lm_control object
#' @param dataset fmri_frame object
#' @param run_num Optional run number for run-specific extraction
#' @param n_time Number of timepoints for this run/chunk
#' @return Integer vector of 1-based timepoint indices to censor, or NULL
#' @keywords internal
#' @noRd
resolve_censor <- function(cfg, dataset = NULL, run_num = NULL, n_time = NULL) {
  censor_spec <- cfg$ar$censor

  # Explicit opt-out: ignore any censor column the dataset carries.
  if (is.character(censor_spec) && identical(censor_spec, "none")) {
    return(NULL)
  }

  # Unset falls back to the dataset, which is what "auto" asks for explicitly.
  if (is.null(censor_spec) || identical(censor_spec, "auto")) {
    return(extract_censor_from_dataset(dataset, run_num, n_time))
  }

  # Explicit censor vector provided in config
  censor <- censor_spec

  # Subset to run if needed and we have global censor vector
  if (!is.null(run_num) && !is.null(dataset) && .is_fmri_frame(dataset)) {
    run_rows <- .dset_run_rows(dataset)
    rows <- run_rows[[run_num]]
    total_time <- sum(lengths(run_rows))
    is_mask <- is.logical(censor) ||
      (is.numeric(censor) && all(censor %in% c(0, 1)))

    if (is_mask && length(censor) == total_time) {
      # Global mask covering all runs
      censor <- censor[rows]
    } else if (is.numeric(censor) && !is_mask && max(censor) > length(rows)) {
      # Indices beyond the run length are global: keep this run's, made local.
      run_start <- rows[1L]
      censor <- censor[censor >= run_start & censor <= rows[length(rows)]] - run_start + 1
    }
  }

  .censor_indices(censor)
}

#' Split contrast objects into t- and F-contrast weight lists
#'
#' F-contrasts also inherit "contrast", so they are excluded from the simple
#' (t) list; routing one through the t-contrast engine fails on its multi-row
#' weight matrix. Each weight carries its contrast's `colind` attribute.
#' @return List with `simple` and `f`, each a named list of weights.
#' @keywords internal
#' @noRd
.contrast_weight_lists <- function(contrast_objects) {
  weights_of <- function(cons) {
    w <- lapply(cons, function(x) {
      w <- x$weights
      attr(w, "colind") <- attr(x, "colind")
      w
    })
    names(w) <- sapply(cons, `[[`, "name")
    w
  }
  list(
    simple = weights_of(Filter(function(x) inherits(x, "contrast") &&
                                 !inherits(x, "Fcontrast"), contrast_objects)),
    f = weights_of(Filter(function(x) inherits(x, "Fcontrast"), contrast_objects))
  )
}

#' Preprocess Run Data
#'
#' @description
#' Applies preprocessing steps (volume weighting and/or soft subspace projection)
#' to the data and design matrices before GLM fitting.
#'
#' @param X Design matrix (time x predictors)
#' @param Y Data matrix (time x voxels)
#' @param cfg fmri_lm_control object
#' @param dataset fmri_frame object (for extracting nuisance if using mask)
#' @param run_num Run number (for run-specific nuisance extraction)
#' @return List with preprocessed X, Y, and preprocessing metadata
#' @details A volume-weight vector covering all runs is sliced to this run's
#'   rows, taken from the dataset's temporal schema.
#' @keywords internal
#' @noRd
preprocess_run_data <- function(X, Y, cfg, dataset = NULL, run_num = NULL) {
  # Missing-data policy (shared chokepoint for all runwise/chunkwise per-run
  # fitting). Under the default "error" policy, fail fast on non-finite response
  # values rather than silently fitting on complete cases (lm na.omit) or
  # producing invalid estimates. "propagate" passes through; the fast solver then
  # surfaces affected voxels as NA coefficients.
  if (identical(cfg$na_action %||% "error", "error") && any(!is.finite(Y))) {
    stop(
      sprintf(
        "Non-finite values (NA/NaN/Inf) in the response data for run %s. Set %s in fmri_lm_control() to fit voxels with missing data as NA, or clean/impute the data.",
        if (is.null(run_num)) "?" else run_num,
        'na_action = "propagate"'
      ),
      call. = FALSE
    )
  }

  preprocess_info <- list(
    volume_weights = NULL,
    dvars = NULL,
    soft_projection = NULL
  )

  # Step 1: Compute and apply volume weights if enabled
  if (cfg$volume_weights$enabled) {
    if (!is.null(cfg$volume_weights$weights)) {
      # Use pre-computed weights
      vol_weights <- cfg$volume_weights$weights
      if (!is.null(run_num) && is.matrix(vol_weights)) {
        # If weights is a matrix with run-specific weights, extract this run
        vol_weights <- vol_weights[, run_num]
      } else if (length(vol_weights) != nrow(X) && !is.null(run_num)) {
        # A full-length vector covers every run: take this run's rows.
        row_ind <- if (!is.null(dataset) && .is_fmri_frame(dataset)) {
          .dset_run_rows(dataset)[[run_num]]
        }
        if (!is.null(row_ind) && length(row_ind) == nrow(X) &&
            max(row_ind) <= length(vol_weights)) {
          vol_weights <- vol_weights[row_ind]
        }
      }
    } else {
      # Compute DVARS-based weights
      result <- volume_weights(Y,
                               method = cfg$volume_weights$method,
                               threshold = cfg$volume_weights$threshold,
                               return_dvars = TRUE)
      vol_weights <- result$weights
      preprocess_info$dvars <- result$dvars
    }

    # Apply volume weights
    weighted <- apply_volume_weights(X, Y, vol_weights)
    X <- weighted$X
    Y <- weighted$Y
    preprocess_info$volume_weights <- vol_weights
  }

  # Step 2: Apply soft subspace projection if enabled
  if (cfg$soft_subspace$enabled) {
    # Get nuisance matrix
    if (!is.null(cfg$soft_subspace$nuisance_matrix)) {
      N <- cfg$soft_subspace$nuisance_matrix
      # Subset to run if needed
      if (!is.null(run_num) && !is.null(dataset)) {
        N <- N[.dset_run_rows(dataset)[[run_num]], , drop = FALSE]
      }
    } else if (!is.null(cfg$soft_subspace$nuisance_mask) && !is.null(dataset)) {
      N <- extract_nuisance_timeseries(dataset, cfg$soft_subspace$nuisance_mask, run = run_num)
    } else {
      stop("soft_subspace enabled but no nuisance_matrix or nuisance_mask provided")
    }

    # Create and apply projection
    proj <- soft_projection(N, lambda = cfg$soft_subspace$lambda, Y = Y)
    cleaned <- apply_soft_projection(proj, Y, X)
    X <- cleaned$X
    Y <- cleaned$Y
    preprocess_info$soft_projection <- proj
  }

  list(X = X, Y = Y, preprocess_info = preprocess_info)
}

#' Design matrix of one run, or of the whole model when `run` is NULL
#'
#' Evaluates the model formula on `term_matrices(model, run)`; pass `tmats`
#' when the caller already holds them.
#' @keywords internal
#' @noRd
.run_design <- function(model, run = NULL, tmats = term_matrices(model, run)) {
  data_env <- list2env(tmats)
  data_env[[".y"]] <- rep(0, nrow(tmats[[1]]))
  model.matrix(get_formula(model), data_env)
}

#' Run-level design and (preprocessed) data
#'
#' Reads the run's data once, builds the run design from
#' `term_matrices(model, run)`, and applies volume weighting / soft
#' projection. `active` marks the observations that take part in the fit
#' (volume weight > 0).
#' @keywords internal
#' @noRd
.prepare_run_inputs <- function(run_chunk, model, cfg, dataset = NULL) {
  Y_run <- as.matrix(.chunk_data(run_chunk))
  X_run <- .run_design(model, run_chunk$chunk_num)
  n_time_run <- nrow(X_run)

  preproc <- preprocess_run_data(X_run, Y_run, cfg, dataset, run_chunk$chunk_num)
  vw <- preproc$preprocess_info$volume_weights
  list(
    X = preproc$X,
    Y = preproc$Y,
    n_time = n_time_run,
    volume_weights = vw,
    active = if (is.null(vw)) NULL else as.numeric(vw) > 0,
    soft_projection = preproc$preprocess_info$soft_projection
  )
}

#' Process a Single Run with Standard OLS/GLS
#'
#' @description
#' Common logic for processing a single run with standard (non-robust) fitting.
#' Handles both IID and AR error structures.
#'
#' @param run_chunk Run data from exec_strategy
#' @param model The fmri_model object
#' @param cfg fmri_lm_control object
#' @param phi_fixed Optional fixed AR parameters
#' @param contrast_weights List of contrast weight vectors/matrices
#' @param dataset Optional fmri_frame for nuisance extraction
#' @return List with fitting results
#' @keywords internal
#' @noRd
process_run_standard <- function(run_chunk, model, cfg, phi_fixed = NULL,
                                 simple_conlist_weights, fconlist_weights,
                                 dataset = NULL) {
  inputs <- .prepare_run_inputs(run_chunk, model, cfg, dataset)
  X_run <- inputs$X
  Y_run <- inputs$Y
  n_time_run <- inputs$n_time

  proj_run <- .fast_preproject(X_run)

  # Resolve censor for this run (for AR estimation)
  censor_run <- resolve_censor(cfg, dataset, run_chunk$chunk_num, n_time_run)

  # AR modeling if needed
  ar_order <- .target_ar_order(cfg$ar) %||% 0L

  if (.temporal_noise_enabled(cfg$ar)) {
    if (as.integer(cfg$ar$q %||% 0L) > 0L) {
      if (!is.null(phi_fixed)) {
        stop("Fixed global AR parameters are not supported for MA models.",
             call. = FALSE)
      }
      temporal_fit <- .iterative_ar_gls_via_fmriAR(
        X = X_run,
        Y = Y_run,
        cfg = cfg$ar,
        max_iter = cfg$ar$iter_gls,
        censor = censor_run
      )
      X_w <- temporal_fit$X_white
      Y_w <- temporal_fit$Y_white
      proj_w <- .fast_preproject(X_w)
      gls <- solve_glm_core(
        glm_context(X = X_w, Y = Y_w, proj = proj_w),
        return_fitted = TRUE
      )

      return(list(
        betas = gls$betas,
        sigma2 = gls$sigma2,
        XtXinv = proj_w$XtXinv,
        dfres = proj_w$dfres,
        rss = gls$rss,
        phi_hat = temporal_fit$ar_coef,
        ma_hat = temporal_fit$ma_coef,
        ar_plan = temporal_fit$plan,
        temporal_diagnostics = list(
          iterations = temporal_fit$iterations,
          converged = temporal_fit$converged
        ),
        X_final = X_w,
        Y_final = Y_w,
        proj_final = proj_w,
        ar_order = ar_order
      ))
    }

    # Estimate AR parameters (excluding censored timepoints)
    phi_hat_run <- if (!is.null(phi_fixed)) {
      phi_fixed
    } else {
      glm_ctx_orig <- glm_context(X = X_run, Y = Y_run, proj = proj_run)
      initial_fit <- solve_glm_core(glm_ctx_orig, return_fitted = TRUE)
      resid_ols <- Y_run - initial_fit$fitted
      .estimate_shared_ar_parameters(
        resid_ols, ar_order, cfg$ar, censor = censor_run,
        design = X_run
      )
    }

    # Any design correction is derived for the initial OLS residual operator.
    # Hold that estimate fixed for the final GLS solve; applying it to later
    # GLS residuals would undo the wrong projection.
    tmp <- ar_whiten_transform(X_run, Y_run, phi_hat_run, cfg$ar$exact_first,
                               censor = censor_run)
    X_w <- tmp$X
    Y_w <- tmp$Y
    proj_w <- .fast_preproject(X_w)

    glm_ctx_whitened <- glm_context(
      X = X_w, Y = Y_w, proj = proj_w, phi_hat = phi_hat_run
    )
    gls <- solve_glm_core(glm_ctx_whitened)
    
    # Return whitened results
    list(
      betas = gls$betas,
      sigma2 = gls$sigma2,
      XtXinv = proj_w$XtXinv,
      dfres = proj_w$dfres,
      rss = gls$rss,
      phi_hat = phi_hat_run,
      ma_hat = NULL,
      temporal_diagnostics = NULL,
      X_final = X_w,
      Y_final = Y_w,
      proj_final = proj_w,
      ar_order = ar_order
    )
  } else {
    # IID case. Zero-weight volumes are omitted by WLS and carry no residual
    # degrees of freedom (as in stats::lm.wfit()).
    if (!is.null(inputs$active)) {
      proj_run$dfres <- as.integer(sum(inputs$active) - proj_run$rank)
    }
    # Soft projection shrinks the nuisance subspace of both X and Y, so the
    # residual df is tr((I - H) P^2) rather than n - rank.
    if (!is.null(inputs$soft_projection)) {
      proj_run$dfres <- .soft_projection_residual_df(inputs$soft_projection, X_run)
    }
    glm_ctx_orig <- glm_context(X = X_run, Y = Y_run, proj = proj_run)
    ols <- solve_glm_core(glm_ctx_orig)
    
    list(
      betas = ols$betas,
      sigma2 = ols$sigma2,
      XtXinv = proj_run$XtXinv,
      dfres = proj_run$dfres,
      rss = ols$rss,
      phi_hat = NULL,
      ma_hat = NULL,
      temporal_diagnostics = NULL,
      X_final = X_run,
      Y_final = Y_run,
      proj_final = proj_run,
      ar_order = 0L
    )
  }
}

#' Process a Single Run with Robust Fitting
#'
#' @description
#' Common logic for processing a single run with robust fitting (no AR).
#'
#' @param run_chunk Run data from exec_strategy
#' @param model The fmri_model object
#' @param cfg fmri_lm_control object
#' @param sigma_fixed Optional fixed robust scale estimate
#' @param dataset Optional fmri_frame for nuisance extraction
#' @return List with robust fitting results
#' @keywords internal
#' @noRd
process_run_robust <- function(run_chunk, model, cfg, sigma_fixed = NULL, dataset = NULL) {
  inputs <- .prepare_run_inputs(run_chunk, model, cfg, dataset)
  X_run <- inputs$X
  Y_run <- inputs$Y

  proj_run <- .fast_preproject(X_run)
  glm_ctx_orig <- glm_context(X = X_run, Y = Y_run, proj = proj_run)

  # Determine sigma_fixed based on scope
  sigma_fixed_for_run <- if (cfg$robust$scale_scope == "global" && !is.null(sigma_fixed)) {
    sigma_fixed
  } else {
    NULL
  }

  robust_fit <- robust_iterative_fitter(
    initial_glm_ctx = glm_ctx_orig,
    cfg_robust_options = cfg$robust,
    X_orig_for_resid = X_run,
    sigma_fixed = sigma_fixed_for_run,
    obs_weights = inputs$volume_weights
  )

  .robust_run_result(robust_fit, X_run, Y_run, proj_run, ncol(Y_run),
                     ar_order = 0L)
}

#' Package a robust run fit
#'
#' Every robust path reports the WLS scale of the final weighted solve
#' (`robust_fit$sigma2`, residual df = non-zero weights minus rank), so
#' runwise and chunkwise fits agree and the fit is reproducible with
#' `stats::lm(weights = robust_weights)` on the (whitened) data. The MAD scale
#' that produced the weights is kept as `sigma_robust` for diagnostics.
#' @keywords internal
#' @noRd
.robust_run_result <- function(robust_fit, X_final, Y_final, proj_final, nvox,
                               ar_order = 0L, phi_hat = NULL) {
  sigma_components <- robust_fit$sigma_robust_scale_components
  if (is.null(sigma_components)) {
    sigma_components <- rep(robust_fit$sigma_robust_scale_final, nvox)
  } else {
    sigma_components <- as.numeric(sigma_components)
    if (length(sigma_components) == 1L) {
      sigma_components <- rep(sigma_components, nvox)
    }
  }
  out <- list(
    betas = robust_fit$betas_robust,
    sigma2 = robust_fit$sigma2,
    XtXinv = robust_fit$XtWXi_final,
    dfres = robust_fit$dfres,
    rank = robust_fit$rank,
    rss = robust_fit$rss,
    robust_weights = robust_fit$robust_weights_final,
    sigma_robust = robust_fit$sigma_robust_scale_final,
    sigma_robust_components = sigma_components,
    X_final = X_final,
    Y_final = Y_final,
    proj_final = proj_final,
    ar_order = ar_order
  )
  if (!is.null(phi_hat)) out$phi_hat <- phi_hat
  out
}

#' Process a Single Run with Combined AR + Robust
#'
#' @description
#' Common logic for processing a single run with both AR and robust fitting.
#' Implements the "whiten then robustly weight" pipeline.
#'
#' With `cfg$robust$reestimate_phi = TRUE` the AR parameters are re-estimated
#' from the *unwhitened* residuals of the robust fit (`Y - X beta_robust`):
#' the whitened residuals carry (approximately) no autocorrelation, so
#' estimating from them would return phi ~ 0 and undo the whitening. The data
#' are then re-whitened, the robust weights re-derived on the re-whitened data
#' and the final betas, weights and scale all come from that last fit.
#'
#' @param run_chunk Run data from exec_strategy
#' @param model The fmri_model object
#' @param cfg fmri_lm_control object
#' @param phi_fixed Optional fixed AR parameters
#' @param sigma_fixed Optional fixed robust scale estimate
#' @param dataset Optional fmri_frame for nuisance extraction
#' @return List with combined fitting results
#' @keywords internal
#' @noRd
process_run_ar_robust <- function(run_chunk, model, cfg, phi_fixed = NULL, sigma_fixed = NULL,
                                  dataset = NULL) {
  pad_phi <- function(phi, order) {
    if (order <= 0L) {
      return(numeric(0))
    }
    if (is.null(phi)) {
      return(rep(0, order))
    }
    phi_vec <- as.numeric(phi)
    if (!length(phi_vec)) {
      return(rep(0, order))
    }
    if (length(phi_vec) < order) {
      phi_vec <- c(phi_vec, rep(0, order - length(phi_vec)))
    }
    phi_vec[seq_len(order)]
  }
  inputs <- .prepare_run_inputs(run_chunk, model, cfg, dataset)
  X_run <- inputs$X
  Y_run <- inputs$Y
  n_time_run <- inputs$n_time

  proj_run <- .fast_preproject(X_run)

  # Resolve censor for this run (for AR estimation)
  censor_run <- resolve_censor(cfg, dataset, run_chunk$chunk_num, n_time_run)

  # AR order
  ar_order <- .target_ar_order(cfg$ar) %||% 0L

  # Step 1: Initial OLS for AR parameter estimation
  glm_ctx_orig <- glm_context(X = X_run, Y = Y_run, proj = proj_run)

  phi_hat_run <- if (!is.null(phi_fixed)) {
    pad_phi(phi_fixed, ar_order)
  } else {
    initial_fit <- solve_glm_core(glm_ctx_orig, return_fitted = TRUE)
    resid_ols <- Y_run - initial_fit$fitted
    pad_phi(.estimate_shared_ar_parameters(
      resid_ols, ar_order, cfg$ar, censor = censor_run, design = X_run
    ), ar_order)
  }

  sigma_fixed_for_run <- if (cfg$robust$scale_scope == "global" && !is.null(sigma_fixed)) {
    sigma_fixed
  } else {
    NULL
  }

  # Steps 2-4: whiten, then robust IRLS on the whitened data.
  whiten_and_fit <- function(phi) {
    tmp <- ar_whiten_transform(X_run, Y_run, phi, cfg$ar$exact_first,
                               censor = censor_run)
    proj_w <- .fast_preproject(tmp$X)
    fit <- robust_iterative_fitter(
      initial_glm_ctx = glm_context(X = tmp$X, Y = tmp$Y, proj = proj_w,
                                    phi_hat = phi),
      cfg_robust_options = cfg$robust,
      X_orig_for_resid = tmp$X,
      sigma_fixed = sigma_fixed_for_run,
      obs_weights = inputs$volume_weights
    )
    list(X = tmp$X, Y = tmp$Y, proj = proj_w, fit = fit)
  }

  stage <- whiten_and_fit(phi_hat_run)

  # Step 5: Optional re-estimation of AR parameters from the unwhitened
  # residuals of the robust fit, then a full refit on the re-whitened data.
  if (isTRUE(cfg$robust$reestimate_phi) && is.null(phi_fixed) && ar_order > 0L) {
    resid_robust <- Y_run - X_run %*% stage$fit$betas_robust
    phi_hat_run <- pad_phi(.estimate_shared_ar_parameters(
      resid_robust, ar_order, cfg$ar, censor = censor_run
    ), ar_order)
    stage <- whiten_and_fit(phi_hat_run)
  }

  .robust_run_result(stage$fit, stage$X, stage$Y, stage$proj, ncol(Y_run),
                     ar_order = ar_order, phi_hat = phi_hat_run)
}

#' Prepare Chunkwise Pre-computation
#'
#' @description
#' Pre-computes transformations needed for chunkwise processing.
#' This includes AR whitening and/or robust weighting of the global design matrix.
#'
#' @param model The fmri_model object
#' @param dataset The fmri_frame object
#' @param cfg fmri_lm_control object
#' @param phi_fixed Optional fixed AR parameters
#' @param sigma_fixed Optional fixed robust scale
#' @param modmat Optional global design matrix (`.run_design(model)`)
#' @return List with pre-computed matrices and metadata
#' @keywords internal
#' @noRd
prepare_chunkwise_matrices <- function(model, dataset, cfg, phi_fixed = NULL, sigma_fixed = NULL,
                                       modmat = NULL) {
  modmat_orig <- modmat %||% .run_design(model)

  # Run structure. Chunks are lazy: each run's data is read once, inside the
  # per-run fit below, and never retained.
  run_chunks <- .dset_run_chunks(dataset)
  run_row_inds <- lapply(run_chunks, `[[`, "row_ind")

  # AR order
  ar_order <- .target_ar_order(cfg$ar) %||% 0L

  ar_modeling <- cfg$ar$struct != "iid"
  robust_modeling <- .fmri_lm_robust_enabled(cfg$robust)

  # Pre-compute per-run information
  run_info <- vector("list", length(run_chunks))

  for (ri in seq_along(run_chunks)) {
    rch <- run_chunks[[ri]]

    # Fit noise locally, but transform rows of the joint design so nuisance
    # coefficients remain run-specific, including unequal nuisance widths.
    X_run_joint <- modmat_orig[run_row_inds[[ri]], , drop = FALSE]
    n_time_run <- length(run_row_inds[[ri]])

    # Resolve censor for this run (it feeds AR estimation and whitening)
    censor_run <- if (ar_modeling) {
      resolve_censor(cfg, dataset, rch$chunk_num, n_time_run)
    } else {
      NULL
    }

    info <- list(
      phi_hat = NULL,
      weights = NULL,
      sqrtw = NULL,
      sigma = NULL,
      X_orig = X_run_joint,
      dummyX0 = matrix(0, n_time_run, 0),
      row_indices = run_row_inds[[ri]],
      censor = censor_run
    )

    if (ar_modeling && robust_modeling) {
      # Combined AR + Robust (censor handled inside process_run_ar_robust)
      res <- process_run_ar_robust(rch, model, cfg, phi_fixed, sigma_fixed, dataset)
      info$phi_hat <- res$phi_hat
      info$weights <- res$robust_weights
      info$sigma <- res$sigma_robust
    } else if (ar_modeling) {
      info$phi_hat <- phi_fixed %||% .run_ar_phi(model, rch, cfg, censor_run)
    } else if (robust_modeling) {
      res <- process_run_robust(rch, model, cfg, sigma_fixed, dataset)
      info$weights <- res$robust_weights
      info$sigma <- res$sigma_robust
    }
    if (!is.null(info$weights)) info$sqrtw <- sqrt(info$weights)
    run_info[[ri]] <- info
  }

  # Build global transformed matrices: whiten (AR) then weight (robust).
  X_global_final <- if (ar_modeling || robust_modeling) {
    do.call(rbind, lapply(run_info, function(info) {
      X_t <- info$X_orig
      if (ar_modeling) {
        X_t <- ar_whiten_transform(X_t, matrix(0, nrow(X_t), 0), info$phi_hat,
                                   cfg$ar$exact_first, censor = info$censor)$X
      }
      if (robust_modeling) X_t <- X_t * info$sqrtw
      X_t
    }))
  } else {
    modmat_orig
  }

  # Pre-compute projection
  proj_global_final <- .fast_preproject(X_global_final)
  if (robust_modeling) {
    # Same residual df as the per-run robust fits: zero-weight observations
    # carry no information (stats::lm.wfit()).
    w_all <- unlist(lapply(run_info, `[[`, "weights"), use.names = FALSE)
    proj_global_final$dfres <- as.integer(sum(w_all > 0) - proj_global_final$rank)
  }

  list(
    X_global = X_global_final,
    proj_global = proj_global_final,
    run_info = run_info,
    run_row_inds = run_row_inds,
    ar_modeling = ar_modeling,
    ar_order = ar_order
  )
}

#' Shared AR coefficients of one run, from all of its voxels
#'
#' Fits the run's own design by OLS and estimates the design-corrected shared
#' AR coefficients from the residuals. Keep the estimate fixed for the GLS
#' solve: a GLS residual update has a different residual-forming operator and
#' cannot reuse the design correction.
#' @keywords internal
#' @noRd
.run_ar_phi <- function(model, run_chunk, cfg, censor = NULL) {
  X_run <- .run_design(model, run_chunk$chunk_num)
  Y_run <- as.matrix(.chunk_data(run_chunk))
  ols <- solve_glm_core(
    glm_context(X = X_run, Y = Y_run, proj = .fast_preproject(X_run)),
    return_fitted = TRUE
  )
  .estimate_shared_ar_parameters(
    Y_run - ols$fitted, .target_ar_order(cfg$ar) %||% 0L, cfg$ar,
    censor = censor, design = X_run
  )
}

#' Process a Chunk of Data
#'
#' @description
#' Processes a single chunk of voxel data using pre-computed matrices.
#'
#' @param chunk_data Matrix of voxel data for this chunk
#' @param precomp Pre-computed matrices from prepare_chunkwise_matrices
#' @param cfg fmri_lm_control object
#' @return List with chunk results
#' @keywords internal
#' @noRd
process_chunk <- function(chunk_data, precomp, cfg) {
  Ymat <- as.matrix(chunk_data)
  
  # Transform Y data to match X transformations
  if (precomp$ar_modeling || .fmri_lm_robust_enabled(cfg$robust)) {
    Ymat_final <- matrix(0, nrow = nrow(Ymat), ncol = ncol(Ymat))

    for (ri in seq_along(precomp$run_info)) {
      rows <- precomp$run_row_inds[[ri]]
      Y_chunk_segment <- Ymat[rows, , drop = FALSE]
      censor_ri <- precomp$run_info[[ri]]$censor

      if (precomp$ar_modeling && .fmri_lm_robust_enabled(cfg$robust)) {
        # AR + Robust: whiten then weight
        phi_hat <- precomp$run_info[[ri]]$phi_hat
        Y_whitened <- ar_whiten_transform(precomp$run_info[[ri]]$dummyX0, Y_chunk_segment, phi_hat, cfg$ar$exact_first, censor = censor_ri)$Y
        Ymat_final[rows, ] <- Y_whitened * precomp$run_info[[ri]]$sqrtw
      } else if (precomp$ar_modeling) {
        # AR only: whiten
        phi_hat <- precomp$run_info[[ri]]$phi_hat
        Ymat_final[rows, ] <- ar_whiten_transform(precomp$run_info[[ri]]$dummyX0, Y_chunk_segment, phi_hat, cfg$ar$exact_first, censor = censor_ri)$Y
      } else if (.fmri_lm_robust_enabled(cfg$robust)) {
        # Robust only: weight
        Ymat_final[rows, ] <- Y_chunk_segment * precomp$run_info[[ri]]$sqrtw
      }
    }
  } else {
    Ymat_final <- Ymat
  }
  
  # Create GLM context and solve
  glm_ctx_chunk <- glm_context(X = precomp$X_global, Y = Ymat_final, proj = precomp$proj_global)
  res <- solve_glm_core(glm_ctx_chunk)
  
  list(
    betas = res$betas,
    sigma2 = res$sigma2,
    rss = res$rss,
    residuals = if (!identical(cfg$variance$method, "model") ||
                    identical(cfg$variance$df, "satterthwaite")) {
      Ymat_final - precomp$X_global %*% res$betas
    } else {
      NULL
    }
  )
}
