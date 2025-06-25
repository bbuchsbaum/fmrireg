# Runwise Linear Model Strategy
# Implementation of the runwise fitting strategy using modular components
#' @importFrom future.apply future_lapply


#' Perform Runwise Linear Modeling on fMRI Dataset
#'
#' This function performs a runwise linear model analysis on an fMRI dataset,
#' running the linear model on each run separately and then pooling results.
#'
#' @param dset An \code{fmri_dataset} object.
#' @param model The \code{fmri_model} used for the analysis.
#' @param contrast_objects The list of full contrast objects.
#' @param cfg An \code{fmri_lm_config} object containing all fitting options.
#' @param verbose Logical. Whether to display progress messages (default is \code{FALSE}).
#' @param use_fast_path Logical. If \code{TRUE}, use matrix-based computation for speed. Default is \code{FALSE}.
#' @param progress Logical. Display a progress bar for run processing. Default is \code{FALSE}.
#' @param phi_fixed Optional fixed AR parameters.
#' @param sigma_fixed Optional fixed robust scale estimate.
#' @param parallel_voxels Logical. If TRUE, process voxels in parallel using future.apply.
#' @return A list containing the combined results from runwise linear model analysis.
#' @keywords internal
runwise_lm <- function(dset, model, contrast_objects, cfg, verbose = FALSE,
                       use_fast_path = FALSE, progress = FALSE,
                       phi_fixed = NULL,
                       sigma_fixed = NULL,
                       parallel_voxels = FALSE) {
  
  # Validate config
  assert_that(inherits(cfg, "fmri_lm_config"), msg = "'cfg' must be an 'fmri_lm_config' object")
  
  # Get run chunks
  chunk_iter <- exec_strategy("runwise")(dset)
  chunks <- collect_chunks(chunk_iter)
  
  # Progress bar setup
  if (progress) {
    pb <- cli::cli_progress_bar("Fitting runs", total = length(chunks), clear = FALSE)
    on.exit(cli::cli_progress_done(id = pb), add = TRUE)
  }
  
  # Extract model components
  form <- get_formula(model)
  tmats <- term_matrices(model)
  event_indices <- attr(tmats, "event_term_indices")
  baseline_indices <- attr(tmats, "baseline_term_indices")
  
  # Global design matrix for Vu calculation
  modmat_global <- design_matrix(model)
  proj_global <- .fast_preproject(modmat_global)
  Vu <- proj_global$XtXinv
  
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
  
  # Determine if we're using voxelwise AR
  if (!use_fast_path && cfg$ar$voxelwise && cfg$ar$struct != "iid") {
    # Voxelwise AR path
    result <- runwise_lm_voxelwise(
      chunks = chunks,
      model = model,
      cfg = cfg,
      simple_conlist_weights = simple_conlist_weights,
      fconlist_weights = fconlist_weights,
      event_indices = event_indices,
      baseline_indices = baseline_indices,
      Vu = Vu,
      verbose = verbose,
      progress = progress,
      parallel_voxels = parallel_voxels
    )
    return(result)
  }
  
  # Process runs
  cres <- if (use_fast_path) {
    runwise_lm_fast(
      chunks = chunks,
      model = model,
      cfg = cfg,
      simple_conlist_weights = simple_conlist_weights,
      fconlist_weights = fconlist_weights,
      event_indices = event_indices,
      baseline_indices = baseline_indices,
      phi_fixed = phi_fixed,
      sigma_fixed = sigma_fixed,
      verbose = verbose,
      progress = progress
    )
  } else {
    runwise_lm_slow(
      chunks = chunks,
      model = model,
      cfg = cfg,
      contrast_objects = contrast_objects,
      event_indices = event_indices,
      baseline_indices = baseline_indices,
      verbose = verbose,
      progress = progress
    )
  }
  
  # Pool results across runs
  pool_runwise_results(cres, event_indices, baseline_indices, Vu)
}

#' Runwise LM Fast Path
#' @keywords internal
#' @noRd
runwise_lm_fast <- function(chunks, model, cfg, simple_conlist_weights, fconlist_weights,
                           event_indices, baseline_indices, phi_fixed = NULL, 
                           sigma_fixed = NULL, verbose = FALSE, progress = FALSE) {
  
  cres <- vector("list", length(chunks))
  ar_order <- switch(cfg$ar$struct,
                     ar1 = 1L,
                     ar2 = 2L,
                     arp = cfg$ar$p,
                     iid = 0L)
  
  for (i in seq_along(chunks)) {
    ym <- chunks[[i]]
    if (verbose) message("Processing run (fast path) ", ym$chunk_num)
    
    # Skip empty runs
    if (ncol(ym$data) == 0 || nrow(ym$data) == 0) {
      warning(paste("Skipping empty run", ym$chunk_num))
      next
    }
    
    # Determine which processing function to use
    if (cfg$robust$type != FALSE && cfg$ar$struct != "iid") {
      # Combined AR + Robust
      res <- process_run_ar_robust(ym, model, cfg, phi_fixed, sigma_fixed)
    } else if (cfg$robust$type != FALSE) {
      # Robust only
      res <- process_run_robust(ym, model, cfg, sigma_fixed)
    } else {
      # Standard (with or without AR)
      res <- process_run_standard(ym, model, cfg, phi_fixed, 
                                 simple_conlist_weights, fconlist_weights)
    }
    
    # Calculate statistics
    actual_vnames <- colnames(res$X_final)
    sigma_vec <- sqrt(res$sigma2)
    
    # Extract robust weights if available
    robust_weights_for_stats <- if (cfg$robust$type != FALSE && !is.null(res$robust_weights)) {
      res$robust_weights
    } else {
      NULL
    }
    
    # Beta statistics
    bstats <- beta_stats_matrix(
      res$betas, 
      res$XtXinv, 
      sigma_vec,
      res$dfres, 
      actual_vnames,
      robust_weights = robust_weights_for_stats,
      ar_order = res$ar_order
    )
    
    # Contrast statistics
    conres <- fit_lm_contrasts_fast(
      res$betas, 
      res$sigma2, 
      res$XtXinv,
      simple_conlist_weights, 
      fconlist_weights,
      res$dfres,
      robust_weights = robust_weights_for_stats,
      ar_order = res$ar_order
    )
    
    cres[[i]] <- list(
      conres = conres,
      bstats = bstats,
      event_indices = event_indices,
      baseline_indices = baseline_indices,
      rss = res$rss,
      rdf = res$dfres,
      resvar = res$sigma2,
      sigma = sigma_vec
    )
    
    if (progress) cli::cli_progress_update()
  }
  
  # Filter out NULL results from skipped empty runs
  Filter(Negate(is.null), cres)
}

#' Runwise LM Slow Path
#' @keywords internal
#' @noRd
runwise_lm_slow <- function(chunks, model, cfg, contrast_objects, 
                           event_indices, baseline_indices,
                           verbose = FALSE, progress = FALSE) {
  
  # Determine fitting function
  lmfun <- if (cfg$robust$type != FALSE) multiresponse_rlm else multiresponse_lm
  
  cres <- vector("list", length(chunks))
  form <- get_formula(model)
  
  for (i in seq_along(chunks)) {
    ym <- chunks[[i]]
    if (verbose) message("Processing run (slow path) ", ym$chunk_num)
    
    tmats <- term_matrices(model, ym$chunk_num)
    vnames <- attr(tmats, "varnames")
    
    data_env <- list2env(tmats)
    data_env$.y <- as.matrix(ym$data)
    ret <- lmfun(form, data_env, contrast_objects, vnames, fcon = NULL)
    
    rss <- colSums(as.matrix(ret$fit$residuals^2))
    rdf <- ret$fit$df.residual
    resvar <- rss / rdf
    sigma <- sqrt(resvar)
    
    cres[[i]] <- list(
      conres = ret$contrasts,
      bstats = ret$bstats,
      event_indices = event_indices,
      baseline_indices = baseline_indices,
      rss = rss,
      rdf = rdf,
      resvar = resvar,
      sigma = sigma
    )
    
    if (progress) cli::cli_progress_update()
  }
  
  cres
}

#' Runwise LM with Voxelwise AR
#' @keywords internal
#' @noRd
runwise_lm_voxelwise <- function(chunks, model, cfg, simple_conlist_weights, fconlist_weights,
                                event_indices, baseline_indices, Vu,
                                verbose = FALSE, progress = FALSE, 
                                parallel_voxels = FALSE) {
  
  ar_order <- switch(cfg$ar$struct,
                     ar1 = 1L,
                     ar2 = 2L,
                     arp = cfg$ar$p,
                     iid = 0L)
  
  if (verbose) message("Using voxelwise AR(", ar_order, ") modeling...")
  
  cres <- vector("list", length(chunks))
  
  for (i in seq_along(chunks)) {
    ym <- chunks[[i]]
    if (verbose) message("Processing run ", ym$chunk_num, " with voxelwise AR")
    
    tmats <- term_matrices(model, ym$chunk_num)
    vnames <- attr(tmats, "varnames")
    
    data_env <- list2env(tmats)
    Y_run <- as.matrix(ym$data)
    data_env$.y <- Y_run
    form <- get_formula(model)
    X_run <- model.matrix(form, data_env)
    
    n_voxels <- ncol(Y_run)
    n_timepoints <- nrow(Y_run)
    p <- ncol(X_run)
    
    # Storage for voxel-wise results
    betas_voxelwise <- matrix(NA_real_, p, n_voxels)
    sigma_voxelwise <- numeric(n_voxels)
    rss_voxelwise <- numeric(n_voxels)
    phi_voxelwise <- matrix(NA_real_, ar_order, n_voxels)
    XtXinv_list <- vector("list", n_voxels)
    
    # Setup for parallel processing if requested
    if (parallel_voxels && requireNamespace("future.apply", quietly = TRUE)) {
      if (verbose) message("Processing voxels in parallel...")
      
      # Use progressr for parallel progress
      if (progress && requireNamespace("progressr", quietly = TRUE)) {
        p_prog <- progressr::progressor(steps = n_voxels)
      }
      
      # Process voxels in parallel
      voxel_results <- future.apply::future_lapply(seq_len(n_voxels), function(v) {
        y_voxel <- Y_run[, v]
        
        # Skip if all zeros or constant
        if (var(y_voxel) < .Machine$double.eps) {
          return(NULL)
        }
        
        # Initial OLS fit
        Y_voxel_mat <- matrix(y_voxel, ncol = 1)
        proj_voxel <- .fast_preproject(X_run)
        glm_ctx_voxel <- glm_context(X = X_run, Y = Y_voxel_mat, proj = proj_voxel)
        initial_fit <- solve_glm_core(glm_ctx_voxel)
        
        # Estimate voxel-specific AR parameters
        resid_ols <- Y_voxel_mat - X_run %*% initial_fit$betas
        phi_voxel <- estimate_ar_parameters(resid_ols[, 1], ar_order)
        
        # AR whitening for this voxel
        tmp <- ar_whiten_transform(X_run, Y_voxel_mat, phi_voxel, cfg$ar$exact_first)
        X_voxel_w <- tmp$X
        Y_voxel_w <- tmp$Y
        
        # Fit on whitened data
        if (cfg$robust$type != FALSE) {
          # Robust fit on whitened data
          proj_voxel_w <- .fast_preproject(X_voxel_w)
          glm_ctx_voxel_w <- glm_context(X = X_voxel_w, Y = Y_voxel_w, 
                                         proj = proj_voxel_w, phi_hat = phi_voxel)
          robust_fit <- robust_iterative_fitter(
            initial_glm_ctx = glm_ctx_voxel_w,
            cfg_robust_options = cfg$robust,
            X_orig_for_resid = X_voxel_w
          )
          
          list(
            betas = robust_fit$betas_robust[, 1],
            rss = sum((Y_voxel_w - X_voxel_w %*% robust_fit$betas_robust)^2),
            sigma = robust_fit$sigma_robust_scale_final,
            phi = phi_voxel,
            XtXinv = robust_fit$XtWXi_final,
            robust_weights = robust_fit$robust_weights_final
          )
        } else {
          # OLS on whitened data
          proj_voxel_w <- .fast_preproject(X_voxel_w)
          glm_ctx_voxel_w <- glm_context(X = X_voxel_w, Y = Y_voxel_w, 
                                         proj = proj_voxel_w, phi_hat = phi_voxel)
          fit_w <- solve_glm_core(glm_ctx_voxel_w)
          
          list(
            betas = fit_w$betas[, 1],
            rss = fit_w$rss,
            sigma = sqrt(fit_w$sigma2),
            phi = phi_voxel,
            XtXinv = proj_voxel_w$XtXinv,
            robust_weights = NULL
          )
        }
      }, future.seed = TRUE)
      
      # Collect results
      for (v in seq_len(n_voxels)) {
        if (!is.null(voxel_results[[v]])) {
          betas_voxelwise[, v] <- voxel_results[[v]]$betas
          rss_voxelwise[v] <- voxel_results[[v]]$rss
          sigma_voxelwise[v] <- voxel_results[[v]]$sigma
          phi_voxelwise[, v] <- voxel_results[[v]]$phi
          XtXinv_list[[v]] <- voxel_results[[v]]$XtXinv
        }
      }
      
    } else {
      # Sequential processing
      for (v in seq_len(n_voxels)) {
        y_voxel <- Y_run[, v]
        
        # Skip if all zeros or constant
        if (var(y_voxel) < .Machine$double.eps) {
          next
        }
        
        # Initial OLS fit
        Y_voxel_mat <- matrix(y_voxel, ncol = 1)
        proj_voxel <- .fast_preproject(X_run)
        glm_ctx_voxel <- glm_context(X = X_run, Y = Y_voxel_mat, proj = proj_voxel)
        initial_fit <- solve_glm_core(glm_ctx_voxel)
        
        # Estimate voxel-specific AR parameters
        resid_ols <- Y_voxel_mat - X_run %*% initial_fit$betas
        phi_voxel <- estimate_ar_parameters(resid_ols[, 1], ar_order)
        phi_voxelwise[, v] <- phi_voxel
        
        # AR whitening for this voxel
        tmp <- ar_whiten_transform(X_run, Y_voxel_mat, phi_voxel, cfg$ar$exact_first)
        X_voxel_w <- tmp$X
        Y_voxel_w <- tmp$Y
        
        # Fit on whitened data
        if (cfg$robust$type != FALSE) {
          # Robust fit on whitened data
          proj_voxel_w <- .fast_preproject(X_voxel_w)
          glm_ctx_voxel_w <- glm_context(X = X_voxel_w, Y = Y_voxel_w, 
                                         proj = proj_voxel_w, phi_hat = phi_voxel)
          robust_fit <- robust_iterative_fitter(
            initial_glm_ctx = glm_ctx_voxel_w,
            cfg_robust_options = cfg$robust,
            X_orig_for_resid = X_voxel_w
          )
          
          betas_voxelwise[, v] <- robust_fit$betas_robust[, 1]
          rss_voxelwise[v] <- sum((Y_voxel_w - X_voxel_w %*% robust_fit$betas_robust)^2)
          sigma_voxelwise[v] <- robust_fit$sigma_robust_scale_final
          XtXinv_list[[v]] <- robust_fit$XtWXi_final
        } else {
          # OLS on whitened data
          proj_voxel_w <- .fast_preproject(X_voxel_w)
          glm_ctx_voxel_w <- glm_context(X = X_voxel_w, Y = Y_voxel_w, 
                                         proj = proj_voxel_w, phi_hat = phi_voxel)
          fit_w <- solve_glm_core(glm_ctx_voxel_w)
          
          betas_voxelwise[, v] <- fit_w$betas[, 1]
          rss_voxelwise[v] <- fit_w$rss
          sigma_voxelwise[v] <- sqrt(fit_w$sigma2)
          XtXinv_list[[v]] <- proj_voxel_w$XtXinv
        }
        
        if (progress && v %% 100 == 0) {
          message("  Processed ", v, " / ", n_voxels, " voxels")
        }
      }
    }
    
    # Calculate statistics and contrasts
    rdf <- n_timepoints - p
    
    # Beta statistics with voxelwise XtXinv
    bstats <- beta_stats_matrix_voxelwise(
      betas_voxelwise,
      XtXinv_list,
      sigma_voxelwise,
      rdf,
      vnames,
      robust_weights_list = if (cfg$robust$type != FALSE) {
        lapply(XtXinv_list, function(x) attr(x, "robust_weights"))
      } else NULL,
      ar_order = ar_order
    )
    
    # Contrasts with voxelwise XtXinv
    conres <- fit_lm_contrasts_voxelwise(
      betas_voxelwise,
      sigma_voxelwise^2,
      XtXinv_list,
      simple_conlist_weights,
      fconlist_weights,
      rdf,
      robust_weights_list = if (cfg$robust$type != FALSE) {
        lapply(XtXinv_list, function(x) attr(x, "robust_weights"))
      } else NULL,
      ar_order = ar_order
    )
    
    cres[[i]] <- list(
      conres = conres,
      bstats = bstats,
      event_indices = event_indices,
      baseline_indices = baseline_indices,
      rss = rss_voxelwise,
      rdf = rdf,
      resvar = rss_voxelwise / rdf,
      sigma = sigma_voxelwise
    )
    
    if (progress) cli::cli_progress_update()
  }
  
  # Pool results
  pool_runwise_results(cres, event_indices, baseline_indices, Vu)
}

#' Pool Runwise Results
#' @keywords internal
#' @noRd
pool_runwise_results <- function(cres, event_indices, baseline_indices, Vu) {
  # Extract components for pooling
  bstats_list <- lapply(cres, `[[`, "bstats")
  conres_list <- lapply(cres, `[[`, "conres")
  
  # Compute overall statistics
  rdf_vals <- unlist(lapply(cres, `[[`, "rdf"))
  sigma_mat <- do.call(rbind, lapply(seq_along(cres), function(i) {
    as.matrix(cres[[i]]$sigma^2) * rdf_vals[i]
  }))
  sigma <- sqrt(colSums(sigma_mat) / sum(rdf_vals))
  rss <- colSums(do.call(rbind, lapply(cres, function(x) as.matrix(x$rss))))
  rdf <- sum(rdf_vals)
  resvar <- rss / rdf
  
  # Pool over runs
  if (length(cres) > 1) {
    meta_con <- meta_contrasts(conres_list)
    # Include all beta indices (event + baseline)
    all_indices <- c(event_indices, baseline_indices)
    meta_beta <- meta_betas(bstats_list, all_indices)
    
    list(
      contrasts = meta_con,
      betas = meta_beta,
      event_indices = event_indices,
      baseline_indices = baseline_indices,
      cov.unscaled = Vu,
      sigma = sigma,
      rss = rss,
      rdf = rdf,
      resvar = resvar
    )
  } else {
    # Single run - need to combine contrasts into single tibble format
    # The conres_list[[1]] is a list of contrast tibbles, but we need a single tibble
    single_contrasts <- conres_list[[1]]
    
    if (length(single_contrasts) > 0) {
      # Combine the list of contrast tibbles into a single tibble
      combined_contrasts <- dplyr::bind_rows(single_contrasts)
    } else {
      # Empty contrasts
      combined_contrasts <- tibble::tibble()
    }
    
    list(
      contrasts = combined_contrasts,  # Single tibble with all contrasts
      betas = bstats_list[[1]],
      event_indices = event_indices,
      baseline_indices = baseline_indices,
      cov.unscaled = Vu,
      sigma = cres[[1]]$sigma,
      rss = cres[[1]]$rss,
      rdf = cres[[1]]$rdf,
      resvar = cres[[1]]$resvar
    )
  }
}
