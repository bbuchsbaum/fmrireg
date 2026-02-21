# Shared Strategy Components for fMRI Linear Models
# Common functions used by both runwise and chunkwise strategies

#' Extract Censor Vector from Dataset
#'
#' @description
#' Extracts and formats censor information from a dataset object.
#' Handles both binary (0/1) vectors and integer index vectors.
#'
#' @param dataset fmri_dataset object (may contain $censor field)
#' @param run_num Optional run number for run-specific extraction
#' @param n_time Optional total number of timepoints (for validation)
#' @return Integer vector of 1-based timepoint indices to censor, or NULL if no censoring
#' @keywords internal
#' @noRd
extract_censor_from_dataset <- function(dataset, run_num = NULL, n_time = NULL) {
  if (is.null(dataset) || is.null(dataset$censor)) {
    return(NULL)
  }

  censor <- dataset$censor

  # Subset to run if needed

  if (!is.null(run_num) && !is.null(dataset$sampling_frame)) {
    blocklens <- fmrihrf::blocklens(dataset$sampling_frame)
    run_start <- if (run_num == 1) 1 else sum(blocklens[1:(run_num - 1)]) + 1
    run_end <- sum(blocklens[1:run_num])
    censor <- censor[run_start:run_end]
  }

  # Convert to indices: if binary (0/1), convert to indices where 1 = censored
  if (is.logical(censor)) {
    censor <- which(censor)
  } else if (is.numeric(censor) && all(censor %in% c(0, 1))) {
    censor <- which(censor == 1)
  }
  # Otherwise assume censor is already integer indices

  if (length(censor) == 0) {
    return(NULL)
  }

  as.integer(censor)
}

#' Resolve Censor Configuration
#'
#' @description
#' Resolves censor specification from config and/or dataset.
#' Handles "auto" mode which extracts from dataset$censor.
#'
#' @param cfg fmri_lm_config object
#' @param dataset fmri_dataset object
#' @param run_num Optional run number for run-specific extraction
#' @param n_time Number of timepoints for this run/chunk
#' @return Integer vector of 1-based timepoint indices to censor, or NULL
#' @keywords internal
#' @noRd
resolve_censor <- function(cfg, dataset = NULL, run_num = NULL, n_time = NULL) {
  censor_spec <- cfg$ar$censor

  if (is.null(censor_spec)) {
    return(NULL)
  }

  # "auto" means extract from dataset
  if (is.character(censor_spec) && identical(censor_spec, "auto")) {
    return(extract_censor_from_dataset(dataset, run_num, n_time))
  }

  # Explicit censor vector provided in config
  censor <- censor_spec

  # Subset to run if needed and we have global censor vector
  if (!is.null(run_num) && !is.null(dataset) && !is.null(dataset$sampling_frame)) {
    blocklens <- fmrihrf::blocklens(dataset$sampling_frame)
    total_time <- sum(blocklens)

    # Check if censor is global (covers all runs)
    if (is.logical(censor) && length(censor) == total_time) {
      run_start <- if (run_num == 1) 1 else sum(blocklens[1:(run_num - 1)]) + 1
      run_end <- sum(blocklens[1:run_num])
      censor <- censor[run_start:run_end]
    } else if (is.numeric(censor) && all(censor %in% c(0, 1)) && length(censor) == total_time) {
      run_start <- if (run_num == 1) 1 else sum(blocklens[1:(run_num - 1)]) + 1
      run_end <- sum(blocklens[1:run_num])
      censor <- censor[run_start:run_end]
    }
    # If censor is indices, they may be global or local - assume global if max > run length
    if (is.numeric(censor) && !all(censor %in% c(0, 1))) {
      run_start <- if (run_num == 1) 1 else sum(blocklens[1:(run_num - 1)]) + 1
      run_end <- sum(blocklens[1:run_num])
      run_len <- run_end - run_start + 1
      if (max(censor) > run_len) {
        # Global indices - filter to this run and adjust
        censor <- censor[censor >= run_start & censor <= run_end] - run_start + 1
      }
    }
  }

  # Convert to indices
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

#' Preprocess Run Data
#'
#' @description
#' Applies preprocessing steps (volume weighting and/or soft subspace projection)
#' to the data and design matrices before GLM fitting.
#'
#' @param X Design matrix (time x predictors)
#' @param Y Data matrix (time x voxels)
#' @param cfg fmri_lm_config object
#' @param dataset fmri_dataset object (for extracting nuisance if using mask)
#' @param run_num Run number (for run-specific nuisance extraction)
#' @return List with preprocessed X, Y, and preprocessing metadata
#' @keywords internal
#' @noRd
preprocess_run_data <- function(X, Y, cfg, dataset = NULL, run_num = NULL) {
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
        sframe <- dataset$sampling_frame
        blocklens <- fmrihrf::blocklens(sframe)
        run_start <- if (run_num == 1) 1 else sum(blocklens[1:(run_num - 1)]) + 1
        run_end <- sum(blocklens[1:run_num])
        N <- N[run_start:run_end, , drop = FALSE]
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

#' Process a Single Run with Standard OLS/GLS
#'
#' @description
#' Common logic for processing a single run with standard (non-robust) fitting.
#' Handles both IID and AR error structures.
#'
#' @param run_chunk Run data from exec_strategy
#' @param model The fmri_model object
#' @param cfg fmri_lm_config object
#' @param phi_fixed Optional fixed AR parameters
#' @param contrast_weights List of contrast weight vectors/matrices
#' @param dataset Optional fmri_dataset for nuisance extraction
#' @return List with fitting results
#' @keywords internal
#' @noRd
process_run_standard <- function(run_chunk, model, cfg, phi_fixed = NULL,
                                 simple_conlist_weights, fconlist_weights,
                                 dataset = NULL) {
  # Extract data
  Y_run <- as.matrix(run_chunk$data)
  tmats_run <- term_matrices(model, run_chunk$chunk_num)

  # Build design matrix
  data_env_run <- list2env(tmats_run)
  n_time_run <- nrow(tmats_run[[1]])
  data_env_run[[".y"]] <- rep(0, n_time_run)
  form <- get_formula(model)
  X_run <- model.matrix(form, data_env_run)

  # Apply preprocessing (volume weights and/or soft subspace projection)
  preproc <- preprocess_run_data(X_run, Y_run, cfg, dataset, run_chunk$chunk_num)
  X_run <- preproc$X
  Y_run <- preproc$Y

  proj_run <- .fast_preproject(X_run)

  # Resolve censor for this run (for AR estimation)
  censor_run <- resolve_censor(cfg, dataset, run_chunk$chunk_num, n_time_run)

  # AR modeling if needed
  ar_order <- switch(cfg$ar$struct,
                     ar1 = 1L,
                     ar2 = 2L,
                     arp = cfg$ar$p,
                     iid = 0L)

  if (cfg$ar$struct != "iid") {
    # Estimate AR parameters (excluding censored timepoints)
    phi_hat_run <- if (!is.null(phi_fixed)) {
      phi_fixed
    } else {
      glm_ctx_orig <- glm_context(X = X_run, Y = Y_run, proj = proj_run)
      initial_fit <- solve_glm_core(glm_ctx_orig, return_fitted = TRUE)
      resid_ols <- Y_run - initial_fit$fitted
      .estimate_ar_parameters_routed(rowMeans(resid_ols), ar_order, censor = censor_run)
    }

    # Iterative GLS
    for (iter in seq_len(cfg$ar$iter_gls)) {
      tmp <- ar_whiten_transform(X_run, Y_run, phi_hat_run, cfg$ar$exact_first, censor = censor_run)
      X_w <- tmp$X
      Y_w <- tmp$Y
      proj_w <- .fast_preproject(X_w)

      glm_ctx_whitened <- glm_context(X = X_w, Y = Y_w, proj = proj_w, phi_hat = phi_hat_run)
      gls <- solve_glm_core(glm_ctx_whitened)

      # Update phi if needed for next iteration
      if (is.null(phi_fixed) && iter < cfg$ar$iter_gls) {
        resid_gls <- Y_run - X_run %*% gls$betas
        phi_hat_run <- .estimate_ar_parameters_routed(rowMeans(resid_gls), ar_order, censor = censor_run)
      }
    }
    
    # Return whitened results
    list(
      betas = gls$betas,
      sigma2 = gls$sigma2,
      XtXinv = proj_w$XtXinv,
      dfres = proj_w$dfres,
      rss = gls$rss,
      phi_hat = phi_hat_run,
      X_final = X_w,
      Y_final = Y_w,
      proj_final = proj_w,
      ar_order = ar_order
    )
  } else {
    # IID case
    glm_ctx_orig <- glm_context(X = X_run, Y = Y_run, proj = proj_run)
    ols <- solve_glm_core(glm_ctx_orig)
    
    list(
      betas = ols$betas,
      sigma2 = ols$sigma2,
      XtXinv = proj_run$XtXinv,
      dfres = proj_run$dfres,
      rss = ols$rss,
      phi_hat = NULL,
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
#' @param cfg fmri_lm_config object
#' @param sigma_fixed Optional fixed robust scale estimate
#' @param dataset Optional fmri_dataset for nuisance extraction
#' @return List with robust fitting results
#' @keywords internal
#' @noRd
process_run_robust <- function(run_chunk, model, cfg, sigma_fixed = NULL, dataset = NULL) {
  # Extract data
  Y_run <- as.matrix(run_chunk$data)
  tmats_run <- term_matrices(model, run_chunk$chunk_num)

  # Build design matrix
  data_env_run <- list2env(tmats_run)
  n_time_run <- nrow(tmats_run[[1]])
  data_env_run[[".y"]] <- rep(0, n_time_run)
  form <- get_formula(model)
  X_run <- model.matrix(form, data_env_run)

  # Apply preprocessing (volume weights and/or soft subspace projection)
  preproc <- preprocess_run_data(X_run, Y_run, cfg, dataset, run_chunk$chunk_num)
  X_run <- preproc$X
  Y_run <- preproc$Y

  proj_run <- .fast_preproject(X_run)
  
  # Create GLM context
  glm_ctx_orig <- glm_context(X = X_run, Y = Y_run, proj = proj_run)
  
  # Determine sigma_fixed based on scope
  sigma_fixed_for_run <- if (cfg$robust$scale_scope == "global" && !is.null(sigma_fixed)) {
    sigma_fixed
  } else {
    NULL
  }
  
  # Robust fitting
  robust_fit <- robust_iterative_fitter(
    initial_glm_ctx = glm_ctx_orig,
    cfg_robust_options = cfg$robust,
    X_orig_for_resid = X_run,
    sigma_fixed = sigma_fixed_for_run
  )

  sigma_components <- robust_fit$sigma_robust_scale_components
  if (is.null(sigma_components)) {
    sigma_components <- rep(robust_fit$sigma_robust_scale_final, ncol(Y_run))
  } else {
    sigma_components <- as.numeric(sigma_components)
    if (length(sigma_components) == 1L) {
      sigma_components <- rep(sigma_components, ncol(Y_run))
    }
  }
  
  list(
    betas = robust_fit$betas_robust,
    sigma2 = sigma_components^2,
    XtXinv = robust_fit$XtWXi_final,
    dfres = robust_fit$dfres,
    rss = colSums((Y_run - X_run %*% robust_fit$betas_robust)^2),
    robust_weights = robust_fit$robust_weights_final,
    sigma_robust = robust_fit$sigma_robust_scale_final,
    sigma_robust_components = sigma_components,
    X_final = X_run,
    Y_final = Y_run,
    proj_final = proj_run,
    ar_order = 0L
  )
}

#' Process a Single Run with Combined AR + Robust
#'
#' @description
#' Common logic for processing a single run with both AR and robust fitting.
#' Implements the "whiten then robustly weight" pipeline.
#'
#' @param run_chunk Run data from exec_strategy
#' @param model The fmri_model object
#' @param cfg fmri_lm_config object
#' @param phi_fixed Optional fixed AR parameters
#' @param sigma_fixed Optional fixed robust scale estimate
#' @param dataset Optional fmri_dataset for nuisance extraction
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
  # Extract data
  Y_run <- as.matrix(run_chunk$data)
  tmats_run <- term_matrices(model, run_chunk$chunk_num)

  # Build design matrix
  data_env_run <- list2env(tmats_run)
  n_time_run <- nrow(tmats_run[[1]])
  data_env_run[[".y"]] <- rep(0, n_time_run)
  form <- get_formula(model)
  X_run <- model.matrix(form, data_env_run)

  # Apply preprocessing (volume weights and/or soft subspace projection)
  preproc <- preprocess_run_data(X_run, Y_run, cfg, dataset, run_chunk$chunk_num)
  X_run <- preproc$X
  Y_run <- preproc$Y

  proj_run <- .fast_preproject(X_run)

  # Resolve censor for this run (for AR estimation)
  censor_run <- resolve_censor(cfg, dataset, run_chunk$chunk_num, n_time_run)

  # AR order
  ar_order <- switch(cfg$ar$struct,
                     ar1 = 1L,
                     ar2 = 2L,
                     arp = cfg$ar$p,
                     iid = 0L)

  # Step 1: Initial OLS for AR parameter estimation
  glm_ctx_orig <- glm_context(X = X_run, Y = Y_run, proj = proj_run)

  phi_hat_run <- if (!is.null(phi_fixed)) {
    pad_phi(phi_fixed, ar_order)
  } else {
    initial_fit <- solve_glm_core(glm_ctx_orig, return_fitted = TRUE)
    resid_ols <- Y_run - initial_fit$fitted
    pad_phi(.estimate_ar_parameters_routed(rowMeans(resid_ols), ar_order, censor = censor_run), ar_order)
  }

  # Step 2: AR whitening
  tmp <- ar_whiten_transform(X_run, Y_run, phi_hat_run, cfg$ar$exact_first, censor = censor_run)
  X_run_w <- tmp$X
  Y_run_w <- tmp$Y
  proj_run_w <- .fast_preproject(X_run_w)
  
  # Step 3: Create whitened GLM context
  glm_ctx_whitened <- glm_context(X = X_run_w, Y = Y_run_w, proj = proj_run_w, phi_hat = phi_hat_run)
  
  # Determine sigma_fixed
  sigma_fixed_for_run <- if (cfg$robust$scale_scope == "global" && !is.null(sigma_fixed)) {
    sigma_fixed
  } else {
    NULL
  }
  
  # Step 4: Robust fitting on whitened data
  robust_fit_run <- robust_iterative_fitter(
    initial_glm_ctx = glm_ctx_whitened,
    cfg_robust_options = cfg$robust,
    X_orig_for_resid = X_run_w,
    sigma_fixed = sigma_fixed_for_run
  )

  sigma_components <- robust_fit_run$sigma_robust_scale_components
  if (is.null(sigma_components)) {
    sigma_components <- rep(robust_fit_run$sigma_robust_scale_final, ncol(Y_run))
  } else {
    sigma_components <- as.numeric(sigma_components)
    if (length(sigma_components) == 1L) {
      sigma_components <- rep(sigma_components, ncol(Y_run))
    }
  }
  
  # Step 5: Optional re-estimation of AR parameters
  if (!is.null(cfg$robust$reestimate_phi) && cfg$robust$reestimate_phi && is.null(phi_fixed)) {
    # Calculate robust residuals on whitened data
    resid_robust_w <- Y_run_w - X_run_w %*% robust_fit_run$betas_robust

    # De-whiten residuals (approximate)
    phi_hat_run_updated <- pad_phi(.estimate_ar_parameters_routed(rowMeans(resid_robust_w), ar_order, censor = censor_run), ar_order)

    # Re-whiten with updated phi
    tmp2 <- ar_whiten_transform(X_run, Y_run, phi_hat_run_updated, cfg$ar$exact_first, censor = censor_run)
    X_run_w2 <- tmp2$X
    Y_run_w2 <- tmp2$Y
    
    # Apply robust weights to re-whitened data
    sqrtw <- sqrt(robust_fit_run$robust_weights_final)
    X_run_wr <- X_run_w2 * sqrtw
    Y_run_wr <- sweep(Y_run_w2, 1, sqrtw, `*`)
    proj_run_wr <- .fast_preproject(X_run_wr)
    
    # Final WLS fit
    glm_ctx_final_wls <- glm_context(X = X_run_wr, Y = Y_run_wr, proj = proj_run_wr)
    final_fit <- solve_glm_core(glm_ctx_final_wls)
    
    list(
      betas = final_fit$betas,
      sigma2 = sigma_components^2,
      XtXinv = proj_run_wr$XtXinv,
      dfres = proj_run_wr$dfres,
      rss = final_fit$rss,
      robust_weights = robust_fit_run$robust_weights_final,
      sigma_robust = robust_fit_run$sigma_robust_scale_final,
      sigma_robust_components = sigma_components,
      phi_hat = phi_hat_run_updated,
      X_final = X_run_w2,
      Y_final = Y_run_w2,
      proj_final = .fast_preproject(X_run_w2),
      ar_order = ar_order
    )
  } else {
    # Use robust fit results directly
    list(
      betas = robust_fit_run$betas_robust,
      sigma2 = sigma_components^2,
      XtXinv = robust_fit_run$XtWXi_final,
      dfres = robust_fit_run$dfres,
      rss = colSums((Y_run_w - X_run_w %*% robust_fit_run$betas_robust)^2),
      robust_weights = robust_fit_run$robust_weights_final,
      sigma_robust = robust_fit_run$sigma_robust_scale_final,
      sigma_robust_components = sigma_components,
      phi_hat = phi_hat_run,
      X_final = X_run_w,
      Y_final = Y_run_w,
      proj_final = proj_run_w,
      ar_order = ar_order
    )
  }
}

#' Pool Results Across Runs
#'
#' @description
#' Combines results from multiple runs into overall statistics.
#'
#' @param run_results List of results from process_run_* functions
#' @return List with pooled statistics
#' @keywords internal
#' @noRd
pool_run_results <- function(run_results) {
  # Extract components
  rdf_vals <- sapply(run_results, `[[`, "dfres")
  sigma_mats <- lapply(seq_along(run_results), function(i) {
    as.matrix(run_results[[i]]$sigma2) * rdf_vals[i]
  })
  
  # Pool sigma
  sigma_mat <- do.call(rbind, sigma_mats)
  sigma <- sqrt(colSums(sigma_mat) / sum(rdf_vals))
  
  # Pool RSS
  rss <- colSums(do.call(rbind, lapply(run_results, function(x) as.matrix(x$rss))))
  
  # Overall degrees of freedom
  rdf <- sum(rdf_vals)
  
  # Overall residual variance
  resvar <- rss / rdf
  
  list(
    sigma = sigma,
    rss = rss,
    rdf = rdf,
    resvar = resvar
  )
}

#' Prepare Chunkwise Pre-computation
#'
#' @description
#' Pre-computes transformations needed for chunkwise processing.
#' This includes AR whitening and/or robust weighting of the global design matrix.
#'
#' @param model The fmri_model object
#' @param dataset The fmri_dataset object
#' @param cfg fmri_lm_config object
#' @param phi_fixed Optional fixed AR parameters
#' @param sigma_fixed Optional fixed robust scale
#' @return List with pre-computed matrices and metadata
#' @keywords internal
#' @noRd
prepare_chunkwise_matrices <- function(model, dataset, cfg, phi_fixed = NULL, sigma_fixed = NULL) {
  # Get global design matrix
  form <- get_formula(model)
  tmats <- term_matrices(model)
  data_env <- list2env(tmats)
  data_env[[".y"]] <- rep(0, nrow(tmats[[1]]))
  modmat_orig <- model.matrix(as.formula(form), data_env)
  
  # Get run structure
  chunk_iter <- exec_strategy("runwise")(dataset)
  run_chunks <- collect_chunks(chunk_iter)
  
  run_row_inds <- lapply(run_chunks, `[[`, "row_ind")
  
  # AR order
  ar_order <- switch(cfg$ar$struct,
                     ar1 = 1L,
                     ar2 = 2L,
                     arp = cfg$ar$p,
                     iid = 0L)
  
  ar_modeling <- cfg$ar$struct != "iid"
  
  # Pre-compute per-run information
  run_info <- vector("list", length(run_chunks))
  
  for (ri in seq_along(run_chunks)) {
    rch <- run_chunks[[ri]]
    
    # Get run-specific matrices
    tmats_run <- term_matrices(model, rch$chunk_num)
    data_env_run <- list2env(tmats_run)
    n_time_run <- nrow(tmats_run[[1]])
    data_env_run[[".y"]] <- rep(0, n_time_run)
    X_run_orig <- model.matrix(form, data_env_run)
    Y_run_full <- as.matrix(rch$data)
    proj_run_orig <- .fast_preproject(X_run_orig)
    
    # Resolve censor for this run
    censor_run <- resolve_censor(cfg, dataset, rch$chunk_num, n_time_run)

    if (ar_modeling && cfg$robust$type != FALSE) {
      # Combined AR + Robust (censor handled inside process_run_ar_robust)
      res <- process_run_ar_robust(rch, model, cfg, phi_fixed, sigma_fixed, dataset)
      run_info[[ri]] <- list(
        phi_hat = res$phi_hat,
        weights = res$robust_weights,
        sqrtw = sqrt(res$robust_weights),
        sigma = res$sigma_robust,
        X_orig = X_run_orig,
        Y_orig = Y_run_full,
        dummyX0 = matrix(0, nrow(X_run_orig), 0),
        row_indices = run_row_inds[[ri]],
        censor = censor_run
      )
    } else if (ar_modeling) {
      # AR only
      phi_hat_run <- if (!is.null(phi_fixed)) {
        phi_fixed
      } else {
        glm_ctx <- glm_context(X = X_run_orig, Y = Y_run_full, proj = proj_run_orig)
        ols <- solve_glm_core(glm_ctx, return_fitted = TRUE)
        resid_ols <- Y_run_full - ols$fitted
        .estimate_ar_parameters_routed(rowMeans(resid_ols), ar_order, censor = censor_run)
      }

      run_info[[ri]] <- list(
        phi_hat = phi_hat_run,
        weights = NULL,
        sqrtw = NULL,
        sigma = NULL,
        X_orig = X_run_orig,
        Y_orig = Y_run_full,
        dummyX0 = matrix(0, nrow(X_run_orig), 0),
        row_indices = run_row_inds[[ri]],
        censor = censor_run
      )
    } else if (cfg$robust$type != FALSE) {
      # Robust only
      res <- process_run_robust(rch, model, cfg, sigma_fixed, dataset)
      run_info[[ri]] <- list(
        phi_hat = NULL,
        weights = res$robust_weights,
        sqrtw = sqrt(res$robust_weights),
        sigma = res$sigma_robust,
        X_orig = X_run_orig,
        Y_orig = Y_run_full,
        dummyX0 = matrix(0, nrow(X_run_orig), 0),
        row_indices = run_row_inds[[ri]]
      )
    }
  }
  
  # Build global transformed matrices
  X_global_final <- if (ar_modeling && cfg$robust$type != FALSE) {
    # AR + Robust: whiten then weight
    X_global_list <- vector("list", length(run_chunks))
    for (ri in seq_along(run_chunks)) {
      X_orig <- run_info[[ri]]$X_orig
      phi_hat <- run_info[[ri]]$phi_hat
      censor_ri <- run_info[[ri]]$censor

      # Whiten
      dummyY <- matrix(0, nrow(X_orig), 0)
      X_whitened <- ar_whiten_transform(X_orig, dummyY, phi_hat, cfg$ar$exact_first, censor = censor_ri)$X

      # Weight
      X_whitened_weighted <- X_whitened * run_info[[ri]]$sqrtw
      X_global_list[[ri]] <- X_whitened_weighted
    }
    do.call(rbind, X_global_list)
  } else if (ar_modeling) {
    # AR only: just whiten
    X_w_list <- vector("list", length(run_chunks))
    for (ri in seq_along(run_chunks)) {
      X_orig <- run_info[[ri]]$X_orig
      phi_hat <- run_info[[ri]]$phi_hat
      censor_ri <- run_info[[ri]]$censor
      dummyY <- matrix(0, nrow(X_orig), 0)
      X_w_list[[ri]] <- ar_whiten_transform(X_orig, dummyY, phi_hat, cfg$ar$exact_first, censor = censor_ri)$X
    }
    do.call(rbind, X_w_list)
  } else if (cfg$robust$type != FALSE) {
    # Robust only: just weight
    X_weighted_list <- vector("list", length(run_chunks))
    for (ri in seq_along(run_chunks)) {
      X_weighted_list[[ri]] <- run_info[[ri]]$X_orig * run_info[[ri]]$sqrtw
    }
    do.call(rbind, X_weighted_list)
  } else {
    # No transformation
    modmat_orig
  }
  
  # Pre-compute projection
  proj_global_final <- .fast_preproject(X_global_final)
  
  list(
    X_global = X_global_final,
    proj_global = proj_global_final,
    run_info = run_info,
    run_row_inds = run_row_inds,
    ar_modeling = ar_modeling,
    ar_order = ar_order
  )
}

#' Process a Chunk of Data
#'
#' @description
#' Processes a single chunk of voxel data using pre-computed matrices.
#'
#' @param chunk_data Matrix of voxel data for this chunk
#' @param precomp Pre-computed matrices from prepare_chunkwise_matrices
#' @param cfg fmri_lm_config object
#' @return List with chunk results
#' @keywords internal
#' @noRd
process_chunk <- function(chunk_data, precomp, cfg) {
  Ymat <- as.matrix(chunk_data)
  
  # Transform Y data to match X transformations
  if (precomp$ar_modeling || cfg$robust$type != FALSE) {
    Ymat_final <- matrix(0, nrow = nrow(Ymat), ncol = ncol(Ymat))

    for (ri in seq_along(precomp$run_info)) {
      rows <- precomp$run_row_inds[[ri]]
      Y_chunk_segment <- Ymat[rows, , drop = FALSE]
      censor_ri <- precomp$run_info[[ri]]$censor

      if (precomp$ar_modeling && cfg$robust$type != FALSE) {
        # AR + Robust: whiten then weight
        phi_hat <- precomp$run_info[[ri]]$phi_hat
        Y_whitened <- ar_whiten_transform(precomp$run_info[[ri]]$dummyX0, Y_chunk_segment, phi_hat, cfg$ar$exact_first, censor = censor_ri)$Y
        Ymat_final[rows, ] <- Y_whitened * precomp$run_info[[ri]]$sqrtw
      } else if (precomp$ar_modeling) {
        # AR only: whiten
        phi_hat <- precomp$run_info[[ri]]$phi_hat
        Ymat_final[rows, ] <- ar_whiten_transform(precomp$run_info[[ri]]$dummyX0, Y_chunk_segment, phi_hat, cfg$ar$exact_first, censor = censor_ri)$Y
      } else if (cfg$robust$type != FALSE) {
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
    rss = res$rss
  )
}
