#' fmriAR Integration Adapter
#'
#' Internal functions for integrating fmriAR package functionality into fmrireg
#' GLM pipeline. This module handles configuration translation, AR estimation,
#' and whitening operations via fmriAR.
#'
#' @keywords internal
#' @noRd
NULL

#' Normalize AR options
#'
#' The AR options reaching the adapter are a `noise_spec()` (legacy user
#' spellings are mapped once, at the boundary, by
#' `.fmri_lm_noise_from_legacy()`). This fills the defaults the adapter relies
#' on and validates them.
#'
#' @param ar_opts AR options list or `noise_spec()`, or NULL for iid.
#' @return List with `struct` and `shared_estimator` always set.
#' @keywords internal
#' @noRd
.normalize_ar_options <- function(ar_opts) {
  ar_opts <- as.list(ar_opts %||% list())
  ar_opts$struct <- ar_opts$struct %||% "iid"
  if (identical(ar_opts$struct, "arp") && is.null(ar_opts$p)) {
    stop("p must be specified for struct='arp'")
  }
  ar_opts$shared_estimator <- ar_opts$shared_estimator %||% "pooled_acvf"
  if (!ar_opts$shared_estimator %in% c("pooled_acvf", "mean_series")) {
    stop("shared_estimator must be 'pooled_acvf' or 'mean_series'", call. = FALSE)
  }
  ar_opts
}

#' Select the residual representation for a shared temporal covariance
#' @keywords internal
#' @noRd
.shared_ar_residual_input <- function(residuals, ar_opts) {
  resid_mat <- as.matrix(residuals)
  estimator <- .normalize_ar_options(ar_opts)$shared_estimator
  if (identical(estimator, "mean_series") && ncol(resid_mat) > 1L) {
    return(matrix(rowMeans(resid_mat), ncol = 1L))
  }
  resid_mat
}

#' Select whether the OLS design correction is valid for a shared estimator
#'
#' Always the design itself (see .estimate_shared_ar_parameters()); kept for
#' the rrr engine, its remaining caller.
#' @keywords internal
#' @noRd
.shared_ar_correction_design <- function(residuals, ar_opts, design) {
  design
}

#' Estimate shared AR parameters under the configured spatial estimand
#' @keywords internal
#' @noRd
.estimate_shared_ar_parameters <- function(residuals, ar_order, ar_opts,
                                           run_indices = NULL, censor = NULL,
                                           design = NULL) {
  # OLS residualization is linear across response columns:
  # rowMeans((I - P_X)Y) = (I - P_X)rowMeans(Y), so the matching design stays
  # the exact residual-forming operator for the mean_series estimator too.
  .estimate_ar_parameters_routed(
    .shared_ar_residual_input(residuals, ar_opts),
    ar_order = ar_order,
    run_indices = run_indices,
    censor = censor,
    design = design,
    ar_options = ar_opts
  )
}

#' Estimate fixed-order AR coefficients via fmriAR
#' @keywords internal
#' @noRd
.estimate_ar_parameters_routed <- function(residuals_vec, ar_order,
                                           run_indices = NULL, censor = NULL,
                                           design = NULL, ar_options = NULL) {
  if (is.null(ar_order) || ar_order <= 0L) {
    return(numeric(0))
  }

  resid_mat <- as.matrix(residuals_vec)
  if (ncol(resid_mat) == 0L || nrow(resid_mat) <= ar_order) {
    return(rep(0, ar_order))
  }

  ar_cfg <- switch(as.integer(ar_order),
    `1` = list(struct = "ar1"),
    `2` = list(struct = "ar2"),
    list(struct = "arp", p = as.integer(ar_order))
  )
  if (!is.null(ar_options)) {
    normalized_options <- .normalize_ar_options(ar_options)
    ar_cfg$global <- isTRUE(normalized_options$global)
    ar_cfg$exact_first <- isTRUE(normalized_options$exact_first)
  }

  plan <- tryCatch(
    .estimate_ar_via_fmriAR(
      residuals = resid_mat,
      cfg = ar_cfg,
      run_indices = run_indices,
      censor = censor,
      design = design
    ),
    error = function(e) {
      if (inherits(e, "fmriAR_design_residual_mismatch")) stop(e)
      NULL
    }
  )

  phi <- NULL
  if (!is.null(plan)) {
    if (!is.null(plan$phi) && length(plan$phi) > 0L) {
      phi <- if (length(plan$phi) == 1L) {
        as.numeric(plan$phi[[1L]])
      } else {
        lapply(plan$phi, as.numeric)
      }
    } else if (!is.null(plan$phi_by_parcel) && length(plan$phi_by_parcel) > 0L) {
      phi <- if (length(plan$phi_by_parcel) == 1L) {
        as.numeric(plan$phi_by_parcel[[1L]])
      } else {
        lapply(plan$phi_by_parcel, as.numeric)
      }
    }
  }

  # No estimate (fmriAR failed or returned no coefficients): no whitening.
  if (is.null(phi) || !length(phi)) {
    phi <- rep(0, ar_order)
  }

  normalize_phi <- function(value) {
    value <- as.numeric(value)
    if (length(value) < ar_order) {
      value <- c(value, rep(0, ar_order - length(value)))
    } else if (length(value) > ar_order) {
      value <- value[seq_len(ar_order)]
    }
    value[!is.finite(value)] <- 0
    value
  }

  if (is.list(phi)) {
    return(lapply(phi, normalize_phi))
  }
  normalize_phi(phi)
}


#' Determine whether temporal whitening is requested
#' @keywords internal
#' @noRd
.temporal_noise_enabled <- function(ar_opts) {
  ar_opts <- .normalize_ar_options(ar_opts)
  !identical(ar_opts$struct, "iid") || as.integer(ar_opts$q %||% 0L) > 0L
}

#' Determine AR order implied by configuration
#'
#' @param ar_opts Normalized AR options list
#' @return Integer AR order (>=0) or NULL if unspecified
#' @keywords internal
#' @noRd
.target_ar_order <- function(ar_opts) {
  ord <- switch(ar_opts$struct,
    "iid" = 0L,
    "ar1" = 1L,
    "ar2" = 2L,
    "arp" = as.integer(ar_opts$p %||% NA_integer_),
    NULL
  )
  if (!is.null(ord) && !is.na(ord)) {
    return(as.integer(ord))
  }
  NULL
}

#' Build per-timepoint run labels from run index list
#'
#' @param n Number of time points
#' @param run_indices Optional list of integer indices per run
#' @return Integer vector of length n with run labels, or NULL for single run
#' @keywords internal
#' @noRd
.build_run_labels <- function(n, run_indices) {
  if (is.null(run_indices) || !length(run_indices)) {
    return(NULL)
  }

  flat <- unlist(run_indices, use.names = FALSE)
  if (!is.numeric(flat) || anyNA(flat) || any(!is.finite(flat)) ||
      any(flat != as.integer(flat)) || length(flat) != n ||
      any(flat < 1L | flat > n) || anyDuplicated(flat) ||
      !setequal(as.integer(flat), seq_len(n)) ||
      any(vapply(run_indices, function(idx) {
        idx <- as.integer(idx)
        length(idx) > 1L && any(diff(idx) != 1L)
      }, logical(1L)))) {
    stop("'run_indices' must partition every timepoint exactly once into contiguous runs",
         call. = FALSE)
  }

  runs <- integer(n)
  for (i in seq_along(run_indices)) {
    idx <- as.integer(run_indices[[i]])
    runs[idx] <- i
  }

  as.integer(runs)
}

#' Derive run indices from a model sampling frame
#' @keywords internal
#' @noRd
.model_run_indices <- function(model, n) {
  sframe <- tryCatch(model$event_model$sampling_frame, error = function(e) NULL)
  if (is.null(sframe)) {
    return(NULL)
  }

  runs <- tryCatch(fmrihrf::blockids(sframe), error = function(e) NULL)
  if (!is.numeric(runs) || length(runs) != n || anyNA(runs)) {
    return(NULL)
  }

  split_idx <- split(seq_len(n), runs)
  lapply(split_idx, as.integer)
}

#' Construct the residual-forming design for separately fitted runs
#' @keywords internal
#' @noRd
.block_diagonal_design <- function(designs) {
  if (!is.list(designs) || !length(designs)) {
    stop("'designs' must be a non-empty list of matrices", call. = FALSE)
  }
  designs <- lapply(designs, function(x) {
    x <- as.matrix(x)
    storage.mode(x) <- "double"
    x
  })
  nr <- vapply(designs, nrow, integer(1L))
  nc <- vapply(designs, ncol, integer(1L))
  out <- matrix(0, sum(nr), sum(nc))
  row_offset <- 0L
  col_offset <- 0L
  out_names <- character(0L)

  for (i in seq_along(designs)) {
    rows <- row_offset + seq_len(nr[[i]])
    cols <- col_offset + seq_len(nc[[i]])
    out[rows, cols] <- designs[[i]]
    block_names <- colnames(designs[[i]])
    if (is.null(block_names)) block_names <- paste0("x", seq_len(nc[[i]]))
    out_names <- c(out_names, paste0("run", i, "::", block_names))
    row_offset <- row_offset + nr[[i]]
    col_offset <- col_offset + nc[[i]]
  }
  colnames(out) <- out_names
  out
}

#' Split residuals according to pooling strategy
#'
#' @param n Number of time points
#' @param pooling Pooling mode ("global" or "run")
#' @param run_indices Optional run index list
#' @return List of integer index vectors
#' @keywords internal
#' @noRd
.split_run_slices <- function(n, pooling, run_indices = NULL) {
  if (identical(pooling, "run") && !is.null(run_indices) && length(run_indices)) {
    lapply(run_indices, as.integer)
  } else {
    list(seq_len(n))
  }
}

#' Estimate AR coefficients of fixed order from residual matrices
#'
#' @param residuals Residual matrix (time x voxels)
#' @param order Target AR order
#' @param pooling Pooling strategy used for estimation
#' @param run_indices Optional run index list (used when pooling == "run")
#' @return List of numeric vectors containing AR coefficients per pool
#' @keywords internal
#' @noRd
.estimate_phi_fixed_order <- function(residuals, order, pooling, run_indices = NULL) {
  if (order <= 0L) {
    return(if (identical(pooling, "run") && !is.null(run_indices) && length(run_indices)) {
      replicate(length(run_indices), numeric(0), simplify = FALSE)
    } else {
      list(numeric(0))
    })
  }

  n <- nrow(residuals)
  slices <- .split_run_slices(n, pooling, run_indices)

  lapply(slices, function(idx) {
    block <- residuals[idx, , drop = FALSE]
    n_eff <- nrow(block)
    if (n_eff <= 1L) {
      return(rep(0, order))
    }
    max_lag <- min(order, n_eff - 1L)
    if (max_lag <= 0L) {
      return(rep(0, order))
    }

    # One pooled fmriAR pass per block avoids a per-voxel loop. Explicit numeric
    # p is a fixed-order AR fit; importantly, this path pools autocovariances
    # across columns instead of fitting the spatial mean series used by ARMA.
    plan <- tryCatch(
      fmriAR::fit_noise(
        resid = block,
        method = "ar",
        p = order,
        p_max = order,
        pooling = "global",
        exact_first = "ar1"
      ),
      error = function(e) NULL
    )

    if (!is.null(plan) && !is.null(plan$phi) && length(plan$phi) > 0L) {
      phi <- as.numeric(plan$phi[[1]])
      if (length(phi) < order) {
        phi <- c(phi, rep(0, order - length(phi)))
      } else if (length(phi) > order) {
        phi <- phi[seq_len(order)]
      }
      return(phi)
    }

    # No explicit phi from fmriAR: zero coefficients (no whitening).
    rep(0, order)
  })
}

#' Choose a stable covariance-tail budget for OLS residual correction
#' @keywords internal
#' @noRd
.ar_correction_lag_budget <- function(design, target_order,
                                      run_indices = NULL, censor = NULL,
                                      max_lag = 25L) {
  design <- as.matrix(design)
  n <- nrow(design)
  order <- suppressWarnings(as.integer(target_order %||% 1L))
  if (!is.finite(order) || order < 1L) order <- 1L
  max_lag <- suppressWarnings(as.integer(max_lag))
  if (!is.finite(max_lag) || max_lag < 1L) max_lag <- 25L

  slices <- if (is.null(run_indices) || !length(run_indices)) {
    list(seq_len(n))
  } else {
    lapply(run_indices, as.integer)
  }
  if (is.logical(censor)) censor <- which(censor)
  censor <- as.integer(censor %||% integer(0L))

  rdf <- vapply(slices, function(idx) {
    keep <- setdiff(idx, censor)
    if (!length(keep)) return(0L)
    as.integer(length(keep) - qr(design[keep, , drop = FALSE])$rank)
  }, integer(1L))

  # Projection mixes the entire covariance tail into lags 0:p. Five lags is a
  # useful floor for low-order AR, while 2p+1 grows with the target. Reserve
  # roughly three residual df per corrected lag and cap the quadratic solve.
  desired <- max(5L, 2L * order + 1L)
  supported <- max(0L, floor(min(rdf) / 3L))
  as.integer(max(order, min(max_lag, desired, supported)))
}

#' Convert fmrireg AR configuration to fmriAR parameters
#'
#' @param cfg An fmri_lm_control object or AR options list
#' @param n_runs Number of runs in the data
#' @return List of fmriAR-compatible parameters
#' @keywords internal
#' @noRd
.fmrireg_to_fmriAR_config <- function(cfg, n_runs = NULL) {
  # Extract AR options from either config object or direct list
  ar_opts <- if (inherits(cfg, "fmri_lm_control")) {
    cfg$ar
  } else if (is.list(cfg) && !is.null(cfg$ar)) {
    cfg$ar
  } else if (is.list(cfg)) {
    cfg
  } else {
    list(struct = "iid")
  }

  ar_opts <- .normalize_ar_options(ar_opts)

  if (is.null(n_runs) || !is.finite(n_runs) || n_runs <= 0L) {
    n_runs <- 1L
  }

  # Map correlation structure to fmriAR parameters. A positive q selects the
  # ARMA estimator; q = 0 retains the segmented/design-corrected AR path.
  q <- as.integer(ar_opts$q %||% 0L)
  if (!is.finite(q) || q < 0L) {
    stop("q must be a nonnegative integer", call. = FALSE)
  }
  method <- if (q > 0L) "arma" else "ar"
  p <- .target_ar_order(ar_opts) %||% 0L
  if (!is.finite(p) || p < 0L) {
    p <- 0L
  }

  # Determine pooling mode
  pooling <- if (isTRUE(ar_opts$global) || n_runs <= 1) {
    "global"
  } else {
    "run"
  }

  # Build fmriAR parameter list
  list(
    method = method,
    p = p,
    p_max = if (p == 0L) 0L else max(6L, p),  # Allow auto-selection up to 6
    q = q,
    pooling = pooling,
    exact_first = if (isFALSE(ar_opts$exact_first)) "none" else "ar1",
    # Number of iterations handled externally
    iter = as.integer(ar_opts$iter_gls %||% 1L)
  )
}

#' Estimate AR parameters using fmriAR
#'
#' @param residuals Residual matrix (time x voxels)
#' @param cfg fmri_lm_control object or AR options list
#' @param run_indices List of indices for each run
#' @param censor Integer vector of 1-based timepoint indices to exclude from
#'   AR estimation, or NULL for no censoring. When provided, these timepoints
#'   are excluded from autocorrelation computation and the time series is
#'
#'   segmented at censor points.
#' @return fmriAR_plan object containing AR parameters
#' @keywords internal
#' @noRd
.estimate_ar_via_fmriAR <- function(residuals, cfg, run_indices = NULL,
                                    censor = NULL, design = NULL) {
  n_runs <- if (is.null(run_indices)) 1L else length(run_indices)

  ar_opts <- if (inherits(cfg, "fmri_lm_control")) cfg$ar else cfg
  ar_opts <- .normalize_ar_options(ar_opts)

  # Convert configuration for fmriAR operations
  ar_params <- .fmrireg_to_fmriAR_config(ar_opts, n_runs)
  runs <- .build_run_labels(nrow(residuals), run_indices)
  target_order <- .target_ar_order(ar_opts)
  correction_max_lag <- 25L
  if (!is.null(design)) {
    design <- as.matrix(design)
    correction_max_lag <- .ar_correction_lag_budget(
      design = design,
      target_order = target_order %||% ar_params$p %||% 1L,
      run_indices = run_indices,
      censor = censor
    )
  }

  # If explicit AR order requested (including 0), estimate via Yule-Walker
  # Note: fixed-order estimation doesn't currently support censor; use fmriAR path
  if (ar_params$q == 0L && !is.null(target_order) &&
      is.null(censor) && is.null(design) &&
      (is.null(run_indices) || length(run_indices) <= 1L)) {
    phi_list <- .estimate_phi_fixed_order(residuals, target_order, ar_params$pooling, run_indices)
    theta_list <- replicate(length(phi_list), numeric(0), simplify = FALSE)
    phi_for_plan <- phi_list
    theta_for_plan <- theta_list

    return(fmriAR::compat$plan_from_phi(
      phi = phi_for_plan,
      theta = theta_for_plan,
      runs = runs,
      pooling = ar_params$pooling,
      exact_first = ar_params$exact_first == "ar1"
    ))
  }

  # Use fmriAR for AR estimation (supports censor)
  fmriAR::fit_noise(
    resid = residuals,
    runs = runs,
    method = ar_params$method,
    p = if (!is.null(target_order)) target_order else ar_params$p,
    q = ar_params$q,
    p_max = if (!is.null(target_order)) target_order else ar_params$p_max,
    exact_first = ar_params$exact_first,
    pooling = ar_params$pooling,
    censor = censor,
    # fmriAR's residual-bias correction is currently AR-only. MA support is
    # therefore restricted at the public boundary to uncensored run pooling,
    # and the ARMA estimator receives the OLS/GLS residuals directly.
    design = if (identical(ar_params$method, "ar")) design else NULL,
    correction_max_lag = correction_max_lag
  )
}

#' Apply AR whitening using fmriAR
#'
#' @param X Design matrix (time x predictors)
#' @param Y Data matrix (time x voxels)
#' @param plan fmriAR_plan object from fit_noise
#' @param run_indices List of indices for each run
#' @param censor Integer vector of 1-based timepoint indices that were censored
#'   during AR estimation. Passed to whiten_apply for consistent handling.
#' @return List with whitened X and Y matrices
#' @keywords internal
#' @noRd
.apply_ar_whitening_via_fmriAR <- function(X, Y, plan, run_indices = NULL, censor = NULL) {
  # Create run labels if needed
  runs <- .build_run_labels(nrow(X), run_indices)

  # Ensure theta list matches phi structure for pure AR models
  if (is.null(plan$theta) || length(plan$theta) == 0L) {
    n_blocks <- if (is.list(plan$phi)) length(plan$phi) else 1L
    plan$theta <- replicate(n_blocks, numeric(0), simplify = FALSE)
  }

  parallel_whiten <- getOption("fmrireg.ar.parallel_whiten", TRUE)
  if (!is.logical(parallel_whiten) || length(parallel_whiten) != 1L || is.na(parallel_whiten)) {
    parallel_whiten <- TRUE
  }

  # Apply whitening (censor is passed for consistent segment handling)
  result <- fmriAR::whiten_apply(
    plan = plan,
    X = X,
    Y = Y,
    runs = runs,
    parallel = isTRUE(parallel_whiten),
    censor = censor
  )

  # Return in expected format
  list(X = result$X, Y = result$Y)
}

#' Whiten X and Y with fixed AR coefficients
#'
#' @param X Design matrix
#' @param Y Data matrix
#' @param phi AR coefficients: one vector (shared), or a list with one vector
#'   per run.
#' @param exact_first Logical for exact first sample scaling
#' @param censor Integer vector of 1-based indices that were censored during
#'   AR estimation. Passed to whiten_apply for consistent segment handling.
#' @param run_indices Optional list of timepoint indices for each run. Whitening
#'   is reset at every run boundary.
#' @param ... Additional arguments (ignored)
#' @return List with whitened X and Y
#' @keywords internal
#' @noRd
ar_whiten_transform <- function(X, Y, phi, exact_first = FALSE, censor = NULL,
                                run_indices = NULL, ...) {
  if (anyNA(X) || anyNA(Y)) {
    stop("NA values detected in 'X' or 'Y' for ar_whiten_transform")
  }

  if (!is.matrix(X)) X <- as.matrix(X)
  if (!is.matrix(Y)) Y <- as.matrix(Y)

  phi_list <- if (is.list(phi)) phi else list(phi)
  pooling_mode <- if (length(phi_list) > 1L) "run" else "global"
  plan <- fmriAR::compat$plan_from_phi(
    phi = if (identical(pooling_mode, "global")) phi_list[[1]] else phi_list,
    theta = replicate(length(phi_list), numeric(0), simplify = FALSE),
    runs = .build_run_labels(nrow(X), run_indices),
    pooling = pooling_mode,
    exact_first = exact_first
  )

  .apply_ar_whitening_via_fmriAR(X, Y, plan, run_indices, censor = censor)
}

#' Estimate AR parameters of one residual series
#'
#' Used by the voxelwise AR path. fmriAR estimates the coefficients; the
#' Yule-Walker fallback (`.estimate_ar()`) covers series fmriAR cannot fit.
#'
#' @param residuals_vec Numeric vector of residuals
#' @param p_order Integer AR order
#' @param censor Integer vector of 1-based indices to exclude from AR estimation,
#'   or NULL for no censoring. Censored timepoints are excluded when computing
#'   autocorrelations.
#' @return Numeric vector of AR coefficients
#' @keywords internal
#' @noRd
estimate_ar_parameters <- function(residuals_vec, p_order, censor = NULL) {
  stopifnot(is.numeric(residuals_vec))
  stopifnot(length(p_order) == 1L, p_order >= 1L)

  if (length(residuals_vec) == 0L) {
    stop("Residual vector must contain at least one value", call. = FALSE)
  }

  if (anyNA(residuals_vec)) {
    stop("NA values detected in 'residuals_vec' for estimate_ar_parameters")
  }

  # A constant (zero-variance) residual series has no estimable AR structure.
  # Short-circuit before consulting fmriAR so the result is deterministic across
  # fmriAR versions (some return phi = 0 here rather than an empty vector).
  res_sd <- stats::sd(residuals_vec)
  if (!is.na(res_sd) && res_sd < sqrt(.Machine$double.eps)) {
    return(numeric(0))
  }

  # Try fmriAR first (supports censor)
  residuals_mat <- matrix(residuals_vec, ncol = 1)
  plan <- tryCatch(
    fmriAR::fit_noise(
      resid = residuals_mat,
      method = "ar",
      p = p_order,
      p_max = p_order,
      censor = censor
    ),
    error = function(e) NULL
  )

  if (!is.null(plan) && !is.null(plan$phi) && length(plan$phi) > 0L) {
    phi <- plan$phi[[1]]
    if (length(phi) > 0L) {
      return(as.numeric(phi))
    }
  }

  # Fall back to Yule-Walker implementation (excludes censored manually if needed)
  if (!is.null(censor) && length(censor) > 0) {
    # For fallback, simply exclude censored points from estimation
    valid_idx <- setdiff(seq_along(residuals_vec), censor)
    if (length(valid_idx) <= p_order) {
      return(rep(0, p_order))
    }
    residuals_vec <- residuals_vec[valid_idx]
  }
  fallback_phi <- .estimate_ar(residuals_vec, p_order)
  if (length(fallback_phi) == 0L) {
    return(numeric(0))
  }

  as.numeric(fallback_phi)
}

#' Estimate AR via Yule-Walker (fallback for estimate_ar_parameters())
#'
#' @param residuals_vec Numeric vector
#' @param p_order AR order
#' @return Numeric vector of AR coefficients
#' @keywords internal
#' @noRd
.estimate_ar <- function(residuals_vec, p_order) {
  stopifnot(is.numeric(residuals_vec))
  stopifnot(length(p_order) == 1L, p_order >= 1L)

  # Center residuals
  res <- residuals_vec - mean(residuals_vec, na.rm = TRUE)

  if (stats::sd(res, na.rm = TRUE) < sqrt(.Machine$double.eps)) {
    return(numeric(0))
  }

  # Use stats::ar for Yule-Walker estimation (matching old behavior)
  ar_fit <- stats::ar(res, order.max = p_order, method = "yule-walker",
                       demean = FALSE, aic = FALSE)

  # Return coefficients
  if (length(ar_fit$ar) > 0) {
    as.numeric(ar_fit$ar)
  } else {
    numeric(0)
  }
}

#' Iterative AR-GLS estimation using fmriAR
#'
#' @param X Design matrix
#' @param Y Data matrix
#' @param cfg fmri_lm_control object
#' @param run_indices List of run indices
#' @param max_iter Maximum iterations (overrides config)
#' @param censor Integer vector of 1-based timepoint indices to exclude from
#'   AR estimation. These points are excluded when computing autocorrelations
#'   but are still whitened and included in the final matrices.
#' @return List with AR plan and whitened matrices
#' @keywords internal
#' @noRd
.iterative_ar_gls_via_fmriAR <- function(X, Y, cfg, run_indices = NULL, max_iter = NULL,
                                         censor = NULL, tol = 5e-3) {
  ar_params <- .fmrireg_to_fmriAR_config(cfg, length(run_indices))
  n_iter <- max_iter %||% ar_params$iter %||% 1L
  # AR residual-bias correction is tied to the initial OLS residual operator.
  # ARMA has no such correction and may continue iterative refinement.
  if (identical(ar_params$method, "ar") && n_iter > 0L) n_iter <- 1L
  min_iter <- as.integer(cfg$min_iter %||% 1L)
  if (!is.finite(min_iter) || min_iter < 1L) min_iter <- 1L

  n_vox <- ncol(Y)
  est_max_vox <- cfg$estimation_max_vox %||% getOption("fmrireg.ar.estimation_max_vox", 64L)
  est_max_vox <- suppressWarnings(as.integer(est_max_vox))
  if (!is.finite(est_max_vox) || est_max_vox <= 0L) est_max_vox <- n_vox
  if (is.na(est_max_vox)) est_max_vox <- n_vox
  est_n <- min(n_vox, est_max_vox)
  est_cols <- if (est_n >= n_vox) {
    seq_len(n_vox)
  } else {
    unique(as.integer(round(seq(1, n_vox, length.out = est_n))))
  }
  Y_est <- Y[, est_cols, drop = FALSE]

  # Handle no AR case
  if ((ar_params$p == 0L && ar_params$q == 0L) || n_iter == 0L) {
    return(list(
      plan = NULL,
      X_white = X,
      Y_white = Y,
      ar_coef = NULL,
      ma_coef = NULL,
      iterations = 0L,
      converged = NA,
      censor = censor
    ))
  }

  # Initial residuals
  coef <- tryCatch(base::qr.solve(X, Y_est), error = function(e) {
    .fast_preproject(X)$Pinv %*% Y_est
  })
  residuals <- Y_est - X %*% coef

  plan <- NULL
  prev_parameters <- NULL
  converged <- FALSE

  # Iterative refinement
  for (iter in seq_len(n_iter)) {
    # Estimate AR from current residuals, excluding censored timepoints
    plan <- .estimate_ar_via_fmriAR(
      .shared_ar_residual_input(residuals, cfg),
      cfg, run_indices, censor = censor, design = X
    )
    phi_now <- if (!is.null(plan) && !is.null(plan$phi) && length(plan$phi)) {
      unlist(lapply(plan$phi, as.numeric), use.names = FALSE)
    } else {
      numeric(0)
    }
    theta_now <- if (!is.null(plan) && !is.null(plan$theta) && length(plan$theta)) {
      unlist(lapply(plan$theta, as.numeric), use.names = FALSE)
    } else {
      numeric(0)
    }
    parameters_now <- c(phi_now, theta_now)

    # For intermediate iterations, whiten only a representative voxel subset to
    # update AR coefficients. Apply full whitening once at the end.
    if (iter < n_iter) {
      whitened_est <- .apply_ar_whitening_via_fmriAR(X, Y_est, plan, run_indices, censor = censor)
    }

    # Stop early when AR coefficients stabilize.
    if (iter >= min_iter &&
        !is.null(prev_parameters) &&
        length(parameters_now) > 0L &&
        length(prev_parameters) == length(parameters_now) &&
        is.finite(tol) &&
        max(abs(parameters_now - prev_parameters)) <= tol) {
      converged <- TRUE
      break
    }
    prev_parameters <- parameters_now

    # Update residuals only for uncorrected ARMA refinement. Corrected AR has
    # n_iter = 1 because later GLS residuals use a different operator.
    if (iter < n_iter) {
      coef <- tryCatch(
        base::qr.solve(whitened_est$X, whitened_est$Y),
        error = function(e) .fast_preproject(whitened_est$X)$Pinv %*% whitened_est$Y
      )
      residuals <- Y_est - X %*% coef
    }
  }

  # Final whitening on full data using the last AR plan.
  # If no plan was estimated (edge case), passthrough original matrices.
  if (is.null(plan)) {
    X_white <- X
    Y_white <- Y
  } else {
    whitened_full <- .apply_ar_whitening_via_fmriAR(X, Y, plan, run_indices, censor = censor)
    X_white <- whitened_full$X
    Y_white <- whitened_full$Y
  }

  # Extract AR coefficients for compatibility
  ar_coef <- if (!is.null(plan)) {
    # Extract phi from plan
    if (!is.null(plan$phi)) {
      plan$phi  # List of AR coefficients
    } else if (!is.null(plan$phi_by_parcel)) {
      # Average across parcels if parcel-based
      phi_list <- plan$phi_by_parcel
      list(rowMeans(sapply(phi_list, function(x) c(x, rep(0, max(lengths(phi_list)) - length(x))))))
    } else {
      NULL
    }
  } else {
    NULL
  }
  ma_coef <- if (!is.null(plan)) {
    if (!is.null(plan$theta)) {
      plan$theta
    } else if (!is.null(plan$theta_by_parcel)) {
      theta_list <- plan$theta_by_parcel
      list(rowMeans(sapply(
        theta_list,
        function(x) c(x, rep(0, max(lengths(theta_list)) - length(x)))
      )))
    } else {
      NULL
    }
  } else {
    NULL
  }

  list(
    plan = plan,
    X_white = X_white,
    Y_white = Y_white,
    ar_coef = ar_coef,
    ma_coef = ma_coef,
    iterations = as.integer(iter),
    converged = if (n_iter <= 1L) NA else converged,
    censor = censor
  )
}
