# Chunkwise Linear Model Strategy
# Implementation of the chunkwise fitting strategy using modular components


#' Perform Chunkwise Linear Modeling on fMRI Dataset
#'
#' This function performs a chunkwise linear model analysis on an fMRI dataset,
#' splitting the dataset into chunks and running the linear model on each chunk.
#'
#' @param x An \code{fmri_frame} object.
#' @param model The \code{fmri_model} used for the analysis.
#' @param contrast_objects The list of full contrast objects.
#' @param nchunks The number of chunks to divide the dataset into.
#' @param cfg An \code{fmri_lm_control} object containing all fitting options.
#' @param verbose Logical. Whether to display progress messages (default is \code{FALSE}).
#' @param use_fast_path Logical. If \code{TRUE}, use matrix-based computation for speed. Default is \code{TRUE}.
#' @param progress Logical. Display a progress bar for chunk processing. Default is \code{FALSE}.
#' @param parallel_chunks Logical. If \code{TRUE}, process chunks with
#'   \code{future.apply::future_lapply()} using the active future plan.
#' @param phi_fixed Optional fixed AR parameters.
#' @param sigma_fixed Optional fixed robust scale estimate.
#' @return A list containing the unpacked chunkwise results.
#' @keywords internal
#' @export
chunkwise_lm.fmri_frame <- function(x, model, contrast_objects, nchunks, cfg,
                                      verbose = FALSE, use_fast_path = TRUE, progress = FALSE,
                                      parallel_chunks = FALSE,
                                      phi_fixed = NULL,
                                      sigma_fixed = NULL, ...) {
  dset <- x
  
  # Validate config
  assert_that(inherits(cfg, "fmri_lm_control"), msg = "'cfg' must be an 'fmri_lm_control' object")
  assert_that(
    is.logical(parallel_chunks) && length(parallel_chunks) == 1L && !is.na(parallel_chunks),
    msg = "'parallel_chunks' must be TRUE or FALSE"
  )
  
  # Get chunks
  chunks <- .dset_chunks(dset, nchunks)
  
  # Progress bar setup
  if (progress && !parallel_chunks) {
    pb <- cli::cli_progress_bar("Fitting chunks", total = length(chunks), clear = FALSE)
    on.exit(cli::cli_progress_done(id = pb), add = TRUE)
  }
  
  # Extract model components
  tmats <- term_matrices(model)
  event_indices <- attr(tmats, "event_term_indices")
  baseline_indices <- attr(tmats, "baseline_term_indices")
  modmat <- .run_design(model, tmats = tmats)

  if (use_fast_path) {
    chunkwise_lm_fast(
      dset = dset,
      chunks = chunks,
      model = model,
      cfg = cfg,
      contrast_objects = contrast_objects,
      event_indices = event_indices,
      baseline_indices = baseline_indices,
      modmat = modmat,
      phi_fixed = phi_fixed,
      sigma_fixed = sigma_fixed,
      verbose = verbose,
      progress = progress,
      parallel_chunks = parallel_chunks
    )
  } else {
    chunkwise_lm_slow(
      dset = dset,
      chunks = chunks,
      model = model,
      cfg = cfg,
      contrast_objects = contrast_objects,
      tmats = tmats,
      vnames = attr(tmats, "varnames"),
      event_indices = event_indices,
      baseline_indices = baseline_indices,
      modmat = modmat,
      # Run rows come from the temporal schema; no data is read.
      run_indices = if (cfg$ar$struct != "iid") .dset_run_rows(dset),
      phi_fixed = phi_fixed,
      verbose = verbose,
      progress = progress,
      parallel_chunks = parallel_chunks
    )
  }
}

#' Apply a function to chunk data, optionally through the active future plan
#' @keywords internal
#' @noRd
.chunkwise_apply <- function(chunks, FUN, parallel_chunks = FALSE, progress = FALSE) {
  if (!isTRUE(parallel_chunks)) {
    return(lapply(seq_along(chunks), function(i) {
      out <- FUN(i, chunks[[i]])
      if (progress) cli::cli_progress_update()
      out
    }))
  }

  if (!requireNamespace("future.apply", quietly = TRUE)) {
    stop("parallel_chunks requires the future.apply package.", call. = FALSE)
  }

  if (progress && requireNamespace("progressr", quietly = TRUE)) {
    return(progressr::with_progress({
      p <- progressr::progressor(steps = length(chunks))
      future.apply::future_lapply(seq_along(chunks), function(i) {
        out <- FUN(i, chunks[[i]])
        p()
        out
      }, future.seed = TRUE)
    }))
  }

  future.apply::future_lapply(seq_along(chunks), function(i) {
    FUN(i, chunks[[i]])
  }, future.seed = TRUE)
}

#' Chunkwise LM Fast Path
#' @keywords internal
#' @noRd
chunkwise_lm_fast <- function(dset, chunks, model, cfg, contrast_objects,
                              event_indices, baseline_indices, modmat,
                              phi_fixed = NULL, sigma_fixed = NULL,
                              verbose = FALSE, progress = FALSE,
                              parallel_chunks = FALSE) {
  
  weights <- .contrast_weight_lists(contrast_objects)
  simple_conlist_weights <- weights$simple
  fconlist_weights <- weights$f

  # Check if we need special handling
  ar_modeling <- cfg$ar$struct != "iid"
  robust_modeling <- .fmri_lm_robust_enabled(cfg$robust)

  if (isTRUE(cfg$volume_weights$enabled) &&
      (ar_modeling || robust_modeling)) {
    stop(
      "Volume weighting is not supported by joint AR or robust fitting; ",
      "use `estimation_spec('runwise_meta')` with the matrix backend.",
      call. = FALSE
    )
  }
  
  precomp <- NULL
  if (ar_modeling || robust_modeling) {
    # Pre-computation phase for AR and/or robust
    if (verbose) message("Pre-computing transformed matrices...")
    
    precomp <- prepare_chunkwise_matrices(model, dset, cfg, phi_fixed, sigma_fixed,
                                          modmat = modmat)
    
    # Process chunks with pre-computed matrices
    cres <- .chunkwise_apply(chunks, function(i, ym) {
      if (verbose) message("Processing chunk (fast path) ", ym$chunk_num)
      
      # Process chunk using pre-computed matrices
      chunk_res <- process_chunk(.chunk_data(ym), precomp, cfg)
      
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
      
      list(
        bstats = bstats,
        contrasts = conres,
        event_indices = event_indices,
        baseline_indices = baseline_indices,
        inference_residuals = if (!identical(cfg$variance$method, "model") ||
                                 identical(cfg$variance$df, "satterthwaite")) {
          chunk_res$residuals
        } else {
          NULL
        }
      )
    }, parallel_chunks = parallel_chunks, progress = progress)
    
  } else {
    # Simple IID OLS. Explicit volume weights are shared across voxel chunks,
    # so transform the design once and each response chunk with the same
    # square-root weights.
    ar_order <- 0L
    explicit_weights <- cfg$volume_weights$weights
    has_explicit_weights <- isTRUE(cfg$volume_weights$enabled) &&
      !is.null(explicit_weights)

    if (isTRUE(cfg$volume_weights$enabled) && is.null(explicit_weights)) {
      stop(
        "Computed volume weights are not supported by joint IID OLS. ",
        "Supply explicit weights with `weights_spec(values = ...)` or use ",
        "`estimation_spec('runwise_meta')`.",
        call. = FALSE
      )
    }
    if (has_explicit_weights && length(explicit_weights) != nrow(modmat)) {
      stop("Length of explicit volume weights must equal the number of timepoints.",
           call. = FALSE)
    }
    if (isTRUE(cfg$soft_subspace$enabled)) {
      warning("Soft-subspace preprocessing with chunkwise OLS requires the ",
              "runwise matrix backend for full support.", call. = FALSE)
    }

    sqrt_volume_weights <- if (has_explicit_weights) sqrt(explicit_weights) else NULL
    fit_modmat <- if (has_explicit_weights) {
      modmat * sqrt_volume_weights
    } else {
      modmat
    }
    proj <- .fast_preproject(fit_modmat)
    if (has_explicit_weights) {
      # Zero-weight observations are omitted by WLS and therefore do not
      # contribute residual degrees of freedom (matching stats::lm.wfit()).
      proj$dfres <- as.integer(sum(explicit_weights > 0) - proj$rank)
    }

    cres <- .chunkwise_apply(chunks, function(i, ym) {
      if (verbose) message("Processing chunk (fast path) ", ym$chunk_num)

      Ymat <- as.matrix(.chunk_data(ym))
      fit_Ymat <- if (has_explicit_weights) {
        sweep(Ymat, 1L, sqrt_volume_weights, `*`)
      } else {
        Ymat
      }
      
      # Create GLM context and solve
      glm_ctx_chunk <- glm_context(X = fit_modmat, Y = fit_Ymat, proj = proj)
      res <- solve_glm_core(glm_ctx_chunk)
      
      # Calculate statistics
      actual_vnames <- colnames(fit_modmat)
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
      
      list(
        bstats = bstats,
        contrasts = conres,
        event_indices = event_indices,
        baseline_indices = baseline_indices,
        inference_residuals = if (!identical(cfg$variance$method, "model") ||
                                 identical(cfg$variance$df, "satterthwaite")) {
          fit_Ymat - fit_modmat %*% res$betas
        } else {
          NULL
        }
      )
    }, parallel_chunks = parallel_chunks, progress = progress)
  }
  
  # Unpack results
  out <- unpack_chunkwise(cres, event_indices, baseline_indices)
  proj_fit <- if (is.null(precomp)) proj else precomp$proj_global
  out$cov.unscaled <- proj_fit$XtXinv
  out$rank <- proj_fit$rank
  robust_weights <- if (!is.null(precomp) && robust_modeling) {
    unlist(lapply(precomp$run_info, `[[`, "weights"), use.names = FALSE)
  }
  out$robust_weights <- robust_weights
  inference_residuals <- lapply(cres, `[[`, "inference_residuals")
  if (length(inference_residuals) &&
      all(vapply(inference_residuals, Negate(is.null), logical(1)))) {
    censor_global <- rep(FALSE, nrow(modmat))
    if (!is.null(precomp)) {
      for (ri in seq_along(precomp$run_info)) {
        local <- precomp$run_info[[ri]]$censor
        rows <- precomp$run_row_inds[[ri]]
        if (is.null(local)) next
        if (is.logical(local)) {
          censor_global[rows] <- local
        } else {
          local <- as.integer(local)
          local <- local[local >= 1L & local <= length(rows)]
          censor_global[rows[local]] <- TRUE
        }
      }
    }
    out$inference_context <- list(
      X = if (is.null(precomp)) fit_modmat else precomp$X_global,
      residuals = do.call(cbind, inference_residuals),
      runs = .dset_run_rows(dset),
      censor = censor_global,
      robust_weights = robust_weights
    )
  }
  out$ar_coef <- if (!is.null(precomp) && isTRUE(precomp$ar_modeling)) {
    lapply(precomp$run_info, `[[`, "phi_hat")
  }
  out
}

#' Chunkwise LM Slow Path
#' @keywords internal
#' @noRd
chunkwise_lm_slow <- function(dset, chunks, model, cfg, contrast_objects,
                              tmats, vnames, event_indices, baseline_indices,
                              modmat, run_indices = NULL, phi_fixed = NULL,
                              verbose = FALSE, progress = FALSE,
                              parallel_chunks = FALSE) {

  if (isTRUE(cfg$volume_weights$enabled)) {
    stop(
      "Volume weighting is not supported by the reference fitting backend; ",
      "use `compute_spec(backend = 'matrix')`.",
      call. = FALSE
    )
  }
  
  form <- get_formula(model)
  proj_modmat <- .fast_preproject(modmat)

  weights <- .contrast_weight_lists(contrast_objects)
  simple_conlist_weights <- weights$simple
  fconlist_weights <- weights$f

  needs_integrated <- !identical(cfg$ar$struct, "iid") || .fmri_lm_robust_enabled(cfg$robust)
  if (needs_integrated && is.null(run_indices)) {
    run_indices <- list(seq_len(nrow(modmat)))
  }
  if (!identical(cfg$ar$struct, "iid")) {
    # Estimate phi once per run from every voxel, as the matrix backend does,
    # and hold it fixed in each chunk's solve. Estimating inside each chunk
    # (from at most 64 of its voxels) made the fit depend on the chunking.
    cfg$ar$phi_fixed <- phi_fixed %||% lapply(.dset_run_chunks(dset), function(rch) {
      .run_ar_phi(model, rch, cfg,
                  censor = resolve_censor(cfg, dset, rch$chunk_num, length(rch$row_ind)))
    })
  }

  cres <- .chunkwise_apply(chunks, function(i, ym) {
    if (verbose) message("Processing chunk ", ym$chunk_num)

    Y_chunk <- as.matrix(.chunk_data(ym))

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
      sigma_vec <- result[["sigma"]] %||% sqrt(result[["sigma2"]])
      sigma_vec <- as.numeric(sigma_vec)
      if (length(sigma_vec) == 1L) {
        sigma_vec <- rep(sigma_vec, ncol(betas))
      }
      dfres <- result$dfres %||% result$df_residual %||% proj_modmat$dfres
      ar_order <- result$ar_order %||% .target_ar_order(cfg$ar) %||% 0L

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

      list(
        bstats = bstats,
        contrasts = conres,
        event_indices = event_indices,
        baseline_indices = baseline_indices,
        ar_coef = result$ar_coef %||% result$phi_hat,
        XtXinv = XtXinv,
        nvox = ncol(betas)
      )
    } else {
      chunk_data_env <- list2env(tmats)
      chunk_data_env[[".y"]] <- Y_chunk

      ret <- multiresponse_lm(form, chunk_data_env, contrast_objects, vnames, fcon = NULL, modmat = modmat)

      list(
        bstats = ret$bstats,
        contrasts = ret$contrasts,
        event_indices = event_indices,
        baseline_indices = baseline_indices,
        ar_coef = NULL
      )
    }
  }, parallel_chunks = parallel_chunks, progress = progress)
  
  # Unpack results
  out <- unpack_chunkwise(cres, event_indices, baseline_indices)
  out$cov.unscaled <- proj_modmat$XtXinv
  if (needs_integrated) {
    # Each chunk was solved on its own whitened and/or robustly weighted
    # design, and its scale (sigma) is on that scale. The covariance stored
    # for post-hoc inference must be that chunk's (X'X)^-1, not the
    # unwhitened design's.
    chunk_cov <- lapply(cres, `[[`, "XtXinv")
    shared <- all(vapply(chunk_cov, function(m) {
      isTRUE(all.equal(m, chunk_cov[[1L]], tolerance = 1e-12,
                       check.attributes = FALSE))
    }, logical(1)))
    if (shared) {
      out$cov.unscaled <- chunk_cov[[1L]]
    } else {
      # One basis per chunk: record it per voxel (list elements share the
      # chunk's matrix, so this costs one pointer per voxel).
      nvox_chunk <- vapply(cres, `[[`, integer(1), "nvox")
      out$covariance_by_voxel <- rep(chunk_cov, times = nvox_chunk)
      out$cov.unscaled <- chunk_cov[[1L]]
    }
  }
  ar_coef_list <- lapply(cres, function(x) x$ar_coef)
  ar_coef_list <- Filter(function(x) !is.null(x) && length(unlist(x)) > 0, ar_coef_list)
  out$ar_coef <- if (length(ar_coef_list) > 0) ar_coef_list else NULL
  out
}
