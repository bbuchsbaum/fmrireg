# Chunkwise Linear Model Strategy
# Implementation of the chunkwise fitting strategy using modular components


#' Perform Chunkwise Linear Modeling on fMRI Dataset
#'
#' This function performs a chunkwise linear model analysis on an fMRI dataset,
#' splitting the dataset into chunks and running the linear model on each chunk.
#'
#' @param x An \code{fmri_dataset} object.
#' @param model The \code{fmri_model} used for the analysis.
#' @param contrast_objects The list of full contrast objects.
#' @param nchunks The number of chunks to divide the dataset into.
#' @param cfg An \code{fmri_lm_config} object containing all fitting options.
#' @param verbose Logical. Whether to display progress messages (default is \code{FALSE}).
#' @param use_fast_path Logical. If \code{TRUE}, use matrix-based computation for speed. Default is \code{FALSE}.
#' @param progress Logical. Display a progress bar for chunk processing. Default is \code{FALSE}.
#' @param phi_fixed Optional fixed AR parameters.
#' @param sigma_fixed Optional fixed robust scale estimate.
#' @return A list containing the unpacked chunkwise results.
#' @keywords internal
chunkwise_lm.fmri_dataset <- function(x, model, contrast_objects, nchunks, cfg,
                                      verbose = FALSE, use_fast_path = FALSE, progress = FALSE,
                                      phi_fixed = NULL,
                                      sigma_fixed = NULL, ...) {
  dset <- x
  
  # Validate config
  assert_that(inherits(cfg, "fmri_lm_config"), msg = "'cfg' must be an 'fmri_lm_config' object")
  
  # Get chunks
  chunk_iter <- exec_strategy("chunkwise", nchunks = nchunks)(dset)
  chunks <- collect_chunks(chunk_iter)
  
  # Progress bar setup
  if (progress) {
    pb <- cli::cli_progress_bar("Fitting chunks", total = length(chunks), clear = FALSE)
    on.exit(cli::cli_progress_done(id = pb), add = TRUE)
  }
  
  # Extract model components
  form <- get_formula(model)
  tmats <- term_matrices(model)
  vnames <- attr(tmats, "varnames")
  event_indices <- attr(tmats, "event_term_indices")
  baseline_indices <- attr(tmats, "baseline_term_indices")
  
  # Global design matrix for Vu calculation
  data_env <- list2env(tmats)
  data_env[[".y"]] <- rep(0, nrow(tmats[[1]]))
  modmat <- model.matrix(as.formula(form), data_env)
  proj_global <- .fast_preproject(modmat)
  Vu <- proj_global$XtXinv
  
  # Build run indices for AR-aware slow path
  run_indices <- NULL
  if (cfg$ar$struct != "iid") {
    run_indices <- lapply(
      collect_chunks(exec_strategy("runwise")(dset)),
      `[[`,
      "row_ind"
    )
  }

  # Process chunks
  if (use_fast_path) {
    # Fast path implementation
    result <- chunkwise_lm_fast(
      dset = dset,
      chunks = chunks,
      model = model,
      cfg = cfg,
      contrast_objects = contrast_objects,
      event_indices = event_indices,
      baseline_indices = baseline_indices,
      Vu = Vu,
      phi_fixed = phi_fixed,
      sigma_fixed = sigma_fixed,
      verbose = verbose,
      progress = progress
    )
  } else {
    # Slow path implementation
    result <- chunkwise_lm_slow(
      chunks = chunks,
      model = model,
      cfg = cfg,
      contrast_objects = contrast_objects,
      tmats = tmats,
      vnames = vnames,
      event_indices = event_indices,
      baseline_indices = baseline_indices,
      Vu = Vu,
      modmat = modmat,
      run_indices = run_indices,
      verbose = verbose,
      progress = progress
    )
  }
  
  result
}

#' Chunkwise LM Fast Path
#' @keywords internal
#' @noRd
chunkwise_lm_fast <- function(dset, chunks, model, cfg, contrast_objects,
                              event_indices, baseline_indices, Vu,
                              phi_fixed = NULL, sigma_fixed = NULL,
                              verbose = FALSE, progress = FALSE) {
  
  # Separate contrast types
  simple_conlist <- Filter(function(x) inherits(x, "contrast"), contrast_objects)
  fconlist <- Filter(function(x) inherits(x, "Fcontrast"), contrast_objects)
  
  # Extract weights with colind attributes
  simple_conlist_weights <- lapply(simple_conlist, function(x) {
    w <- x$weights
    attr(w, "colind") <- attr(x, "colind")
    w
  })
  names(simple_conlist_weights) <- sapply(simple_conlist, `[[`, "name")
  
  fconlist_weights <- lapply(fconlist, function(x) {
    w <- x$weights
    attr(w, "colind") <- attr(x, "colind")
    w
  })
  names(fconlist_weights) <- sapply(fconlist, `[[`, "name")
  
  # Check if we need special handling
  ar_modeling <- cfg$ar$struct != "iid"
  robust_modeling <- cfg$robust$type != FALSE
  
  if (ar_modeling || robust_modeling) {
    # Pre-computation phase for AR and/or robust
    if (verbose) message("Pre-computing transformed matrices...")
    
    precomp <- prepare_chunkwise_matrices(model, dset, cfg, phi_fixed, sigma_fixed)
    
    # Process chunks with pre-computed matrices
    cres <- vector("list", length(chunks))
    
    for (i in seq_along(chunks)) {
      ym <- chunks[[i]]
      if (verbose) message("Processing chunk (fast path) ", ym$chunk_num)
      
      # Process chunk using pre-computed matrices
      chunk_res <- process_chunk(ym$data, precomp, cfg)
      
      # Calculate statistics
      actual_vnames <- colnames(precomp$X_global)
      sigma_vec <- sqrt(chunk_res$sigma2)
      
      # Build global robust weights aligned with precomp$X_global row order.
      robust_weights_for_stats <- if (robust_modeling) {
        run_weights <- lapply(precomp$run_info, `[[`, "weights")
        if (!any(vapply(run_weights, is.null, logical(1)))) {
          unlist(run_weights, use.names = FALSE)
        } else {
          NULL
        }
      } else {
        NULL
      }
      
      # Beta statistics
      bstats <- beta_stats_matrix(
        chunk_res$betas,
        precomp$proj_global$XtXinv,
        sigma_vec,
        precomp$proj_global$dfres,
        actual_vnames,
        robust_weights = robust_weights_for_stats,
        ar_order = precomp$ar_order
      )
      
      # Contrast statistics
      conres <- fit_lm_contrasts_fast(
        chunk_res$betas,
        chunk_res$sigma2,
        precomp$proj_global$XtXinv,
        simple_conlist_weights,
        fconlist_weights,
        precomp$proj_global$dfres,
        robust_weights = robust_weights_for_stats,
        ar_order = precomp$ar_order
      )
      
      cres[[i]] <- list(
        bstats = bstats,
        contrasts = conres,
        event_indices = event_indices,
        baseline_indices = baseline_indices
      )
      
      if (progress) cli::cli_progress_update()
    }
    
  } else {
    # Simple OLS case - no pre-computation needed
    ar_order <- 0L
    form <- get_formula(model)
    tmats <- term_matrices(model)
    data_env <- list2env(tmats)
    data_env[[".y"]] <- rep(0, nrow(tmats[[1]]))
    modmat <- model.matrix(as.formula(form), data_env)

    # Check if preprocessing is requested but strategy may not fully support it
    if (cfg$volume_weights$enabled || cfg$soft_subspace$enabled) {
      warning("Preprocessing (volume_weights/soft_subspace) with chunkwise OLS ",
              "requires use of runwise strategy or AR/robust options for full support.",
              " Consider using strategy='runwise' with use_fast_path=TRUE.", call. = FALSE)
    }

    proj <- .fast_preproject(modmat)

    cres <- vector("list", length(chunks))

    for (i in seq_along(chunks)) {
      ym <- chunks[[i]]
      if (verbose) message("Processing chunk (fast path) ", ym$chunk_num)

      Ymat <- as.matrix(ym$data)
      
      # Create GLM context and solve
      glm_ctx_chunk <- glm_context(X = modmat, Y = Ymat, proj = proj)
      res <- solve_glm_core(glm_ctx_chunk)
      
      # Calculate statistics
      actual_vnames <- colnames(modmat)
      sigma_vec <- sqrt(res$sigma2)
      
      # Beta statistics
      bstats <- beta_stats_matrix(
        res$betas,
        proj$XtXinv,
        sigma_vec,
        proj$dfres,
        actual_vnames,
        robust_weights = NULL,
        ar_order = ar_order
      )
      
      # Contrast statistics
      conres <- fit_lm_contrasts_fast(
        res$betas,
        res$sigma2,
        proj$XtXinv,
        simple_conlist_weights,
        fconlist_weights,
        proj$dfres,
        robust_weights = NULL,
        ar_order = ar_order
      )
      
      cres[[i]] <- list(
        bstats = bstats,
        contrasts = conres,
        event_indices = event_indices,
        baseline_indices = baseline_indices
      )
      
      if (progress) cli::cli_progress_update()
    }
  }
  
  # Unpack results
  out <- unpack_chunkwise(cres, event_indices, baseline_indices)
  out$cov.unscaled <- Vu
  ar_coef_list <- lapply(cres, function(x) x$ar_coef)
  ar_coef_list <- Filter(function(x) !is.null(x) && length(unlist(x)) > 0, ar_coef_list)
  out$ar_coef <- if (length(ar_coef_list) > 0) ar_coef_list else NULL
  out
}

#' Chunkwise LM Slow Path
#' @keywords internal
#' @noRd
chunkwise_lm_slow <- function(chunks, model, cfg, contrast_objects,
                              tmats, vnames, event_indices, baseline_indices,
                              Vu, modmat, run_indices = NULL,
                              verbose = FALSE, progress = FALSE) {
  
  # Determine fitting function
  lmfun <- if (cfg$robust$type != FALSE) multiresponse_rlm else multiresponse_lm
  
  # Setup data environment
  data_env <- list2env(tmats)
  form <- get_formula(model)
  proj_modmat <- .fast_preproject(modmat)

  # Prepare contrast weights for integrated solver path
  simple_conlist <- Filter(function(x) inherits(x, "contrast"), contrast_objects)
  simple_conlist_weights <- lapply(simple_conlist, function(x) {
    w <- x$weights
    attr(w, "colind") <- attr(x, "colind")
    w
  })
  if (length(simple_conlist) > 0) {
    names(simple_conlist_weights) <- sapply(simple_conlist, `[[`, "name")
  }

  fconlist <- Filter(function(x) inherits(x, "Fcontrast"), contrast_objects)
  fconlist_weights <- lapply(fconlist, function(x) {
    w <- x$weights
    attr(w, "colind") <- attr(x, "colind")
    w
  })
  if (length(fconlist) > 0) {
    names(fconlist_weights) <- sapply(fconlist, `[[`, "name")
  }

  needs_integrated <- !(cfg$ar$struct %in% c("iid", "none")) || cfg$robust$type != FALSE
  if (needs_integrated && is.null(run_indices)) {
    run_indices <- list(seq_len(nrow(modmat)))
  }

  cres <- vector("list", length(chunks))
  
  for (i in seq_along(chunks)) {
    ym <- chunks[[i]]
    if (verbose) message("Processing chunk ", ym$chunk_num)

    Y_chunk <- as.matrix(ym$data)

    if (needs_integrated) {
      result <- solve_integrated_glm(
        X = modmat,
        Y = Y_chunk,
        config = cfg,
        run_indices = run_indices
      )

      betas <- result$betas %||% result$coefficients
      if (!is.matrix(betas)) {
        betas <- as.matrix(betas)
      }

      XtXinv <- result$XtXinv %||% proj_modmat$XtXinv
      sigma_vec <- result$sigma %||% sqrt(result$sigma2)
      sigma_vec <- as.numeric(sigma_vec)
      if (length(sigma_vec) == 1L) {
        sigma_vec <- rep(sigma_vec, ncol(betas))
      }
      dfres <- result$dfres %||% result$df_residual %||% proj_modmat$dfres
      ar_order <- result$ar_order %||% switch(cfg$ar$struct,
        ar1 = 1L,
        ar2 = 2L,
        arp = cfg$ar$p,
        0L
      )

      bstats <- beta_stats_matrix(
        betas,
        XtXinv,
        sigma_vec,
        dfres,
        colnames(modmat),
        robust_weights = result$robust_weights,
        ar_order = ar_order
      )

      conres <- fit_lm_contrasts_fast(
        betas,
        sigma_vec^2,
        XtXinv,
        simple_conlist_weights,
        fconlist_weights,
        dfres,
        robust_weights = result$robust_weights,
        ar_order = ar_order
      )

      cres[[i]] <- list(
        bstats = bstats,
        contrasts = conres,
        event_indices = event_indices,
        baseline_indices = baseline_indices,
        ar_coef = result$ar_coef %||% result$phi_hat
      )
    } else {
      data_env[[".y"]] <- Y_chunk

      ret <- lmfun(form, data_env, contrast_objects, vnames, fcon = NULL, modmat = modmat)

      cres[[i]] <- list(
        bstats = ret$bstats,
        contrasts = ret$contrasts,
        event_indices = event_indices,
        baseline_indices = baseline_indices,
        ar_coef = NULL
      )
    }

    if (progress) cli::cli_progress_update()
  }
  
  # Unpack results
  out <- unpack_chunkwise(cres, event_indices, baseline_indices)
  out$cov.unscaled <- Vu
  ar_coef_list <- lapply(cres, function(x) x$ar_coef)
  ar_coef_list <- Filter(function(x) !is.null(x) && length(unlist(x)) > 0, ar_coef_list)
  out$ar_coef <- if (length(ar_coef_list) > 0) ar_coef_list else NULL
  out
}
