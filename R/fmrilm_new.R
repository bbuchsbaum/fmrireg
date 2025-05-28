#' @title Core fMRI Linear Model Functions
#' @description Main API functions for fitting linear models to fMRI data
#' @keywords internal
#' @importFrom assertthat assert_that

#' Fit a Linear Model to fMRI Data
#'
#' @description
#' This function fits a linear model to fMRI data using various strategies and options.
#' It supports both standard and robust regression, with optional AR error modeling.
#'
#' @param formula A formula specifying the linear model to fit to the fMRI data
#' @param block A factor or character vector indicating the block structure
#' @param baseline_model An optional baseline model object
#' @param dataset An fmri_dataset object containing the fMRI data
#' @param strategy The fitting strategy: "runwise" or "chunkwise"
#' @param nchunks Number of data chunks for "chunkwise" strategy
#' @param robust Robust regression method: FALSE (none), "huber", or "bisquare"
#' @param robust_options List of robust fitting options or NULL to use defaults
#' @param ar_options List of AR modeling options or NULL to use defaults
#' @param extra_nuisance Additional nuisance regressors (NULL, matrix, or formula)
#' @param keep_extra_nuisance_in_model Whether to include extra nuisance in the model
#' @param use_fast_path Use optimized matrix computations
#' @param durations Event durations (default 0 for impulse events)
#' @param drop_empty Drop empty factor levels from block
#' @param verbose Print progress messages
#' @param progress Show progress bar
#' @return An fmri_lm object containing the fitted model and results
#' @export
fmri_lm <- function(formula, block, baseline_model = NULL, dataset, durations = 0, drop_empty = TRUE,
                    strategy = c("runwise", "chunkwise"), nchunks = 10,
                    robust = c(FALSE, "huber", "bisquare"),
                    robust_options = NULL,
                    ar_options = NULL,
                    extra_nuisance = NULL,
                    keep_extra_nuisance_in_model = FALSE,
                    use_fast_path = FALSE,
                    verbose = FALSE,
                    progress = FALSE) {
  
  # Validate inputs
  strategy <- match.arg(strategy)
  robust <- if (is.logical(robust) && !robust) FALSE else match.arg(robust)
  
  # Create configuration object
  cfg <- fmri_lm_control(
    robust_options = robust_options,
    ar_options = ar_options
  )
  
  # Override robust type if specified at top level
  if (!is.null(robust)) {
    cfg$robust$type <- robust
  }
  
  # Create the fMRI model
  fmrimod <- create_fmri_model(
    formula = formula,
    block = block,
    baseline_model = baseline_model,
    dataset = dataset,
    drop_empty = drop_empty,
    durations = durations
  )
  
  # Fit the model
  result <- fmri_lm_fit(
    fmrimod = fmrimod,
    dataset = dataset,
    strategy = strategy,
    nchunks = nchunks,
    cfg = cfg,
    extra_nuisance = extra_nuisance,
    keep_extra_nuisance_in_model = keep_extra_nuisance_in_model,
    use_fast_path = use_fast_path,
    verbose = verbose,
    progress = progress
  )
  
  # Return fmri_lm object
  structure(
    list(
      model = fmrimod,
      dataset = dataset,
      result = result
    ),
    class = "fmri_lm",
    strategy = strategy,
    config = cfg
  )
}

#' Fit an fMRI Linear Model
#'
#' @description
#' Internal function that performs the actual model fitting.
#'
#' @param fmrimod An fmri_model object
#' @param dataset An fmri_dataset object
#' @param strategy Fitting strategy: "runwise" or "chunkwise"
#' @param nchunks Number of chunks for chunkwise strategy
#' @param cfg An fmri_lm_config object with all fitting options
#' @param extra_nuisance Additional nuisance regressors
#' @param keep_extra_nuisance_in_model Whether to keep extra nuisance
#' @param use_fast_path Use optimized computations
#' @param verbose Print progress messages
#' @param progress Show progress bar
#' @return List with fitting results
#' @keywords internal
fmri_lm_fit <- function(fmrimod, dataset, strategy = c("runwise", "chunkwise"),
                        nchunks = 10, cfg,
                        extra_nuisance = NULL,
                        keep_extra_nuisance_in_model = FALSE,
                        use_fast_path = FALSE, verbose = FALSE, progress = FALSE) {
  
  strategy <- match.arg(strategy)
  
  # Validate config
  assert_that(inherits(cfg, "fmri_lm_config"), msg = "'cfg' must be an 'fmri_lm_config' object")
  
  # Error checking
  assert_that(inherits(fmrimod, "fmri_model"), msg = "'fmrimod' must be an 'fmri_model' object")
  assert_that(inherits(dataset, "fmri_dataset"), msg = "'dataset' must be an 'fmri_dataset' object")
  assert_that(is.logical(use_fast_path), msg = "'use_fast_path' must be logical")
  if (strategy == "chunkwise") {
    assert_that(is.numeric(nchunks) && nchunks > 0, msg = "'nchunks' must be a positive number")
  }
  
  # Get contrast info grouped by term
  contrast_info_by_term <- contrast_weights(fmrimod$event_model)
  full_design_colnames <- colnames(design_matrix(fmrimod))
  processed_conlist <- list()
  
  # Process contrasts term by term to correctly assign column indices
  for (term_name in names(contrast_info_by_term)) {
    term_contrasts <- contrast_info_by_term[[term_name]]
    
    # Get col_indices from design matrix instead of term_indices from event model
    col_indices <- attr(fmrimod$event_model$design_matrix, "col_indices")
    
    if (length(term_contrasts) > 0 && !is.null(col_indices) && !is.null(col_indices[[term_name]])) {
      # Get the column indices directly from col_indices
      colind <- col_indices[[term_name]]
      
      if (length(colind) == 0) {
        warning(paste("No column indices found for term:", term_name))
        next
      }
      
      # Apply colind attribute to each contrast spec within this term
      processed_term_contrasts <- lapply(term_contrasts, function(con_spec) {
        if (inherits(con_spec, "contrast") || inherits(con_spec, "Fcontrast")) {
          # Set the colind attribute on the contrast weights for the slow path
          attr(con_spec$weights, "colind") <- colind
          # Also set it directly on the contrast object
          attr(con_spec, "colind") <- colind
        } else {
          warning(paste("Item in contrast list for term", term_name, "is not a contrast or Fcontrast object."))
        }
        con_spec
      })
      processed_conlist <- c(processed_conlist, processed_term_contrasts)
    } else if (length(term_contrasts) > 0 && (is.null(col_indices) || is.null(col_indices[[term_name]]))) {
      warning(paste("Contrasts found for term '", term_name, "' but col_indices are missing in the event model design matrix."))
    }
  }
  
  # Handle global AR and robust scale estimation
  phi_global <- NULL
  sigma_global <- NULL
  
  if (cfg$ar$global && cfg$ar$struct != "iid") {
    ar_order <- switch(cfg$ar$struct,
                       ar1 = 1L,
                       ar2 = 2L,
                       arp = cfg$ar$p)
    
    run_chunks <- exec_strategy("runwise")(dataset)
    form <- get_formula(fmrimod)
    resid_vec <- numeric(0)
    
    for (rch in run_chunks) {
      tmats_run <- term_matrices(fmrimod, rch$chunk_num)
      data_env_run <- list2env(tmats_run)
      n_time_run <- nrow(tmats_run[[1]])
      data_env_run[[".y"]] <- rep(0, n_time_run)
      X_run <- model.matrix(form, data_env_run)
      proj_run <- .fast_preproject(X_run)
      Y_run <- as.matrix(rch$data)
      
      glm_ctx <- glm_context(X = X_run, Y = Y_run, proj = proj_run)
      ols <- solve_glm_core(glm_ctx, return_fitted = TRUE)
      resid_vec <- c(resid_vec, rowMeans(Y_run - ols$fitted))
    }
    
    phi_global <- estimate_ar_parameters(resid_vec, ar_order)
    cfg$ar$iter_gls <- 1L
  }
  
  if (cfg$robust$type != FALSE && cfg$robust$scale_scope == "global") {
    run_chunks <- exec_strategy("runwise")(dataset)
    form <- get_formula(fmrimod)
    row_med_all <- numeric(0)
    
    for (rch in run_chunks) {
      tmats_run <- term_matrices(fmrimod, rch$chunk_num)
      data_env_run <- list2env(tmats_run)
      n_time_run <- nrow(tmats_run[[1]])
      data_env_run[[".y"]] <- rep(0, n_time_run)
      X_run <- model.matrix(form, data_env_run)
      proj_run <- .fast_preproject(X_run)
      Y_run <- as.matrix(rch$data)
      
      glm_ctx <- glm_context(X = X_run, Y = Y_run, proj = proj_run)
      initial_fit <- solve_glm_core(glm_ctx)
      resid <- Y_run - X_run %*% initial_fit$betas
      row_med_all <- c(row_med_all, apply(abs(resid), 1, median))
    }
    
    sigma_global <- median(row_med_all) / 0.6745
  }
  
  # Dispatch to appropriate strategy
  if (strategy == "chunkwise") {
    chunkwise_lm.fmri_dataset(
      dset = dataset,
      model = fmrimod,
      contrast_objects = processed_conlist,
      nchunks = nchunks,
      cfg = cfg,
      verbose = verbose,
      use_fast_path = use_fast_path,
      progress = progress,
      phi_fixed = phi_global,
      sigma_fixed = sigma_global,
      extra_nuisance = extra_nuisance,
      keep_extra_nuisance_in_model = keep_extra_nuisance_in_model
    )
  } else {
    runwise_lm(
      dset = dataset,
      model = fmrimod,
      contrast_objects = processed_conlist,
      cfg = cfg,
      verbose = verbose,
      use_fast_path = use_fast_path,
      progress = progress,
      phi_fixed = phi_global,
      sigma_fixed = sigma_global,
      extra_nuisance = extra_nuisance,
      keep_extra_nuisance_in_model = keep_extra_nuisance_in_model,
      parallel_voxels = cfg$ar$voxelwise
    )
  }
}

#' Fit Contrasts for Linear Model
#'
#' @description
#' This function calculates contrasts for a fitted linear model.
#'
#' @param fit The fitted linear model object.
#' @param conlist The list of contrasts.
#' @param fcon The F-contrasts (currently not used).
#' @param vnames The names of the variables in the linear model.
#' @param se Whether to compute standard errors (default: TRUE).
#' @return A list containing contrasts, beta statistics, and the fitted model.
#' @keywords internal
fit_lm_contrasts <- function(fit, conlist, fcon, vnames, se = TRUE) {
  conres <- if (!is.null(conlist)) {
    ret <- lapply(conlist, function(con) {
      # Extract colind from the contrast object's attributes
      colind <- attr(con, "colind")
      if (is.null(colind)) {
        warning(paste("Missing colind attribute for contrast:", con$name %||% "unnamed"))
        return(NULL) # Skip this contrast
      }
      estimate_contrast(con, fit, colind)
    })
    # Filter out NULL results
    ret <- ret[!sapply(ret, is.null)]
    names(ret) <- sapply(conlist[!sapply(ret, is.null)], function(x) x$name %||% "unnamed")
    ret
  } else {
    list()
  }
  
  bstats <- beta_stats(fit, vnames, se = se)
  list(contrasts = conres, bstats = bstats, fit = fit)
}

#' Fit Multiresponse Linear Model
#'
#' This function fits a linear model to multiple responses in an fMRI dataset.
#'
#' @param form The formula used to define the linear model.
#' @param data_env The environment containing the data to be used in the linear model.
#' @param conlist The list of contrasts used in the analysis.
#' @param vnames The names of the variables used in the linear model.
#' @param fcon The F-contrasts used in the analysis.
#' @param modmat The model matrix (default is \code{NULL}, which will calculate the model matrix using the formula).
#' @return A list containing the results from the multiresponse linear model analysis.
#' @keywords internal
multiresponse_lm <- function(form, data_env, conlist, vnames, fcon, modmat = NULL) {
  lm_fit <- if (is.null(modmat)) {
    lm(as.formula(form), data = data_env)
  } else {
    lm.fit(modmat, data_env$.y)
  }
  
  # Use the actual column names from the model matrix instead of vnames
  actual_vnames <- if (is.null(modmat)) {
    names(coef(lm_fit))
  } else {
    colnames(modmat)
  }
  
  fit_lm_contrasts(lm_fit, conlist, fcon, actual_vnames)
}

#' Unpack Chunkwise Results
#' @keywords internal
#' @noRd
unpack_chunkwise <- function(cres, event_indices, baseline_indices) {
  # Beta Processing
  cbetas <- lapply(cres, function(x) x$bstats) %>% dplyr::bind_rows()
  dat_beta <- cbetas$data %>% dplyr::bind_rows()
  
  # Check validity
  valid_estimates_idx <- sapply(dat_beta$estimate, function(x) !is.null(x) && is.matrix(x) && nrow(x) > 0)
  if (!any(valid_estimates_idx)) {
      stop("No valid beta estimates found across chunks in unpack_chunkwise.")
  }
  dat_beta_valid <- dat_beta[valid_estimates_idx, , drop = FALSE]

  # Concatenate beta results across chunks
  estimate_beta <- do.call(rbind, dat_beta_valid$estimate)
  se_beta <- do.call(rbind, dat_beta_valid$se)
  stat_beta <- do.call(rbind, dat_beta_valid$stat)
  prob_beta <- do.call(rbind, dat_beta_valid$prob)
  sigma_beta <- do.call(c, dat_beta_valid$sigma)
  
  # Re-package combined beta results
  cbetas_out <- dplyr::tibble(
    type = cbetas$type[1],
    stat_type = cbetas$stat_type[1],
    df.residual = cbetas$df.residual[1],
    conmat = list(NULL),
    colind = list(NULL),
    data = list(
      dplyr::tibble(
        estimate = list(estimate_beta),   
        se = list(se_beta),               
        stat = list(stat_beta),            
        prob = list(prob_beta),            
        sigma = list(sigma_beta)          
      )
    )
  )

  # Contrast Processing
  ncon <- if (length(cres) > 0 && !is.null(cres[[1]]$contrasts) && length(cres[[1]]$contrasts) > 0) {
      length(cres[[1]]$contrasts)
  } else { 0 }

  if (ncon > 0) {
    contab <- lapply(cres, function(x) { 
        cons <- x$contrasts 
        if (!is.list(cons)) cons <- list(cons)
        if (length(cons) > 0 && !is.null(names(cons))) {
             dplyr::bind_rows(cons, .id = "contrast_internal_name")
        } else if (length(cons) > 0) {
             warning("Contrast list per chunk lacks names, attempting bind_rows without .id")
             dplyr::bind_rows(cons)
        } else {
             dplyr::tibble()
        }
    }) %>% dplyr::bind_rows()

    if (nrow(contab) == 0) {
        con <- dplyr::tibble()
    } else {
        grouping_vars <- intersect(c("name", "type"), names(contab))
        if (length(grouping_vars) == 0) stop("Cannot group contrasts: 'name' or 'type' column missing.")
        
        gsplit <- contab %>% dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>% dplyr::group_split()

        con <- lapply(gsplit, function(g) {
            dat <- g$data %>% dplyr::bind_rows() 
            
            estimate_full <- dat$estimate
            se_full <- dat$se
            stat_full <- dat$stat
            prob_full <- dat$prob
            sigma_full <- if ("sigma" %in% names(dat)) dat$sigma else NULL
            
            combined_data_tibble <- dplyr::tibble(
                estimate = list(estimate_full), 
                se = list(se_full),             
                stat = list(stat_full),          
                prob = list(prob_full)           
            )
            if (!is.null(sigma_full)) {
                combined_data_tibble$sigma = list(sigma_full)
            }

            g %>% dplyr::select(-data) %>% dplyr::slice_head() %>% 
                dplyr::mutate(data = list(combined_data_tibble))
                
        }) %>% dplyr::bind_rows() 
    }
  } else {
    con <- dplyr::tibble()
  }

  list(
    betas = cbetas_out,
    contrasts = con,
    event_indices = event_indices,
    baseline_indices = baseline_indices
  )
}