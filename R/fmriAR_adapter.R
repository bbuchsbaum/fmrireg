#' fmriAR Integration Adapter
#'
#' Internal functions for integrating fmriAR package functionality into fmrireg
#' GLM pipeline. This module handles configuration translation, AR estimation,
#' and whitening operations via fmriAR.
#'
#' @keywords internal
#' @noRd
NULL

# Internal operator
`%||%` <- function(x, y) if (!is.null(x)) x else y

#' Normalize AR options to ensure struct/cor_struct coherence
#'
#' @param ar_opts List with AR configuration
#' @return List with normalized fields
#' @keywords internal
#' @noRd
.normalize_ar_options <- function(ar_opts) {
  if (is.null(ar_opts)) {
    return(list(struct = "iid", cor_struct = "none"))
  }

  # Work on a copy to avoid side effects
  ar_opts <- as.list(ar_opts)

  # Prefer explicit struct, otherwise fall back to cor_struct
  struct <- ar_opts$struct
  cor_struct <- ar_opts$cor_struct

  if (is.null(struct) && !is.null(cor_struct)) {
    struct <- cor_struct
  }

  if (is.null(struct)) {
    struct <- "iid"
  }

  # Treat "none" as iid for compatibility
  if (identical(struct, "none")) {
    struct <- "iid"
  }

  # Handle convenience shorthand like "ar5" by mapping to arp
  if (!identical(struct, "arp") && startsWith(as.character(struct), "ar")) {
    p_val <- suppressWarnings(as.integer(sub("ar", "", as.character(struct))))
    if (!is.na(p_val) && p_val >= 0L) {
      if (p_val == 0L) {
        struct <- "iid"
      } else if (p_val > 4L) {
        struct <- "arp"
        ar_opts$p <- ar_opts$p %||% p_val
      }
    }
  }

  if (identical(struct, "arp")) {
    if (is.null(ar_opts$p)) {
      # Try to derive p from cor_struct like "ar3"
      if (!is.null(cor_struct) && startsWith(cor_struct, "ar")) {
        ar_opts$p <- suppressWarnings(as.integer(sub("ar", "", cor_struct)))
      }
    }
    if (is.null(ar_opts$p)) {
      stop("p must be specified for struct='arp'")
    }
  }

  # Ensure cor_struct is populated for downstream consumers
  if (is.null(cor_struct)) {
    cor_struct <- switch(struct,
      "iid" = "none",
      "arp" = paste0("ar", ar_opts$p),
      struct
    )
  }

  # Map iter (legacy name) to iter_gls when needed
  if (is.null(ar_opts$iter_gls) && !is.null(ar_opts$iter)) {
    ar_opts$iter_gls <- ar_opts$iter
  }

  ar_opts$struct <- struct
  ar_opts$cor_struct <- cor_struct

  ar_opts
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
    "ar3" = 3L,
    "ar4" = 4L,
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

  runs <- integer(n)
  for (i in seq_along(run_indices)) {
    idx <- as.integer(run_indices[[i]])
    runs[idx] <- i
  }

  if (any(runs == 0L)) {
    runs[runs == 0L] <- 1L
  }

  as.integer(runs)
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

    accum <- rep(0, order)
    used  <- 0L

    for (j in seq_len(ncol(block))) {
      resid_vec <- block[, j]
      resid_vec <- resid_vec[is.finite(resid_vec)]
      if (length(resid_vec) <= order) {
        next
      }

      phi <- estimate_ar_parameters(resid_vec, order)
      if (length(phi) < order) {
        phi <- c(phi, rep(0, order - length(phi)))
      } else if (length(phi) > order) {
        phi <- phi[seq_len(order)]
      }

      accum <- accum + phi
      used <- used + 1L
    }

    if (used == 0L) {
      return(rep(0, order))
    }

    accum / used
  })
}

#' Convert fmrireg AR configuration to fmriAR parameters
#'
#' @param cfg An fmri_lm_config object or ar_options list
#' @param n_runs Number of runs in the data
#' @return List of fmriAR-compatible parameters
#' @keywords internal
#' @noRd
.fmrireg_to_fmriAR_config <- function(cfg, n_runs = NULL) {
  # Extract AR options from either config object or direct list
  ar_opts <- if (inherits(cfg, "fmri_lm_config")) {
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

  # Map correlation structure to fmriAR parameters
  method <- "ar"  # Default to AR (not ARMA)
  p <- switch(as.character(ar_opts$struct),
    "iid" = 0L,
    "ar1" = 1L,
    "ar2" = 2L,
    "ar3" = 3L,
    "ar4" = 4L,
    "arp" = as.integer(ar_opts$p %||% stop("p must be specified for struct='arp'")),
    0L  # Default to no AR
  )
  p <- as.integer(p)
  if (!is.finite(p) || p < 0L) {
    p <- 0L
  }

  # Determine pooling mode
  pooling <- if (isTRUE(ar_opts$global) || n_runs <= 1) {
    "global"
  } else {
    "run"
  }

  # Handle voxelwise option (future: will use parcel-based approach)
  if (isTRUE(ar_opts$voxelwise)) {
    message("Note: voxelwise AR will be implemented via parcel pooling in fmriAR")
    # For now, fall back to run-based pooling
  }

  # Build fmriAR parameter list
  list(
    method = method,
    p = p,
    p_max = if (p == 0L) 0L else max(6L, p),  # Allow auto-selection up to 6
    q = 0L,  # No MA terms for now
    pooling = pooling,
    exact_first = if (isFALSE(ar_opts$exact_first)) "none" else "ar1",
    # Number of iterations handled externally
    iter = as.integer(ar_opts$iter_gls %||% ar_opts$iter %||% 1L)
  )
}

#' Estimate AR parameters using fmriAR
#'
#' @param residuals Residual matrix (time x voxels)
#' @param cfg fmri_lm_config object or AR options list
#' @param run_indices List of indices for each run
#' @return fmriAR_plan object containing AR parameters
#' @keywords internal
#' @noRd
.estimate_ar_via_fmriAR <- function(residuals, cfg, run_indices = NULL) {
  n_runs <- if (is.null(run_indices)) 1L else length(run_indices)

  ar_opts <- if (inherits(cfg, "fmri_lm_config")) cfg$ar else cfg
  ar_opts <- .normalize_ar_options(ar_opts)

  # Convert configuration for fmriAR operations
  ar_params <- .fmrireg_to_fmriAR_config(ar_opts, n_runs)
  runs <- .build_run_labels(nrow(residuals), run_indices)
  target_order <- .target_ar_order(ar_opts)

  # If explicit AR order requested (including 0), estimate via Yule-Walker
  if (!is.null(target_order)) {
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

  # Fallback to fmriAR auto-selection when order not specified
  fmriAR::fit_noise(
    resid = residuals,
    runs = runs,
    method = ar_params$method,
    p = ar_params$p,
    q = ar_params$q,
    p_max = ar_params$p_max,
    exact_first = ar_params$exact_first,
    pooling = ar_params$pooling
  )
}

#' Apply AR whitening using fmriAR
#'
#' @param X Design matrix (time x predictors)
#' @param Y Data matrix (time x voxels)
#' @param plan fmriAR_plan object from fit_noise
#' @param run_indices List of indices for each run
#' @return List with whitened X and Y matrices
#' @keywords internal
#' @noRd
.apply_ar_whitening_via_fmriAR <- function(X, Y, plan, run_indices = NULL) {
  # Create run labels if needed
  runs <- .build_run_labels(nrow(X), run_indices)

  # Ensure theta list matches phi structure for pure AR models
  if (is.null(plan$theta) || length(plan$theta) == 0L) {
    n_blocks <- if (is.list(plan$phi)) length(plan$phi) else 1L
    plan$theta <- replicate(n_blocks, numeric(0), simplify = FALSE)
  }

  # Apply whitening
  result <- fmriAR::whiten_apply(
    plan = plan,
    X = X,
    Y = Y,
    runs = runs,
    parallel = FALSE
  )

  # Return in expected format
  list(X = result$X, Y = result$Y)
}

#' Iterative AR-GLS estimation using fmriAR
#'
#' @param X Design matrix
#' @param Y Data matrix
#' @param cfg fmri_lm_config object
#' @param run_indices List of run indices
#' @param max_iter Maximum iterations (overrides config)
#' @return List with AR plan and whitened matrices
#' @keywords internal
#' @noRd
.iterative_ar_gls_via_fmriAR <- function(X, Y, cfg, run_indices = NULL, max_iter = NULL) {
  ar_params <- .fmrireg_to_fmriAR_config(cfg, length(run_indices))
  n_iter <- max_iter %||% ar_params$iter %||% 1L

  # Handle no AR case
  if (ar_params$p == 0 || ar_params$p_max == 0 || n_iter == 0) {
    return(list(
      plan = NULL,
      X_white = X,
      Y_white = Y,
      ar_coef = NULL
    ))
  }

  # Initial residuals
  coef <- tryCatch(base::qr.solve(X, Y), error = function(e) {
    qr.coef(qr(X, LAPACK = TRUE), Y)
  })
  residuals <- Y - X %*% coef

  plan <- NULL
  X_white <- X
  Y_white <- Y

  # Iterative refinement
  for (iter in seq_len(n_iter)) {
    # Estimate AR from current residuals
    plan <- .estimate_ar_via_fmriAR(residuals, cfg, run_indices)

    # Apply whitening
    whitened <- .apply_ar_whitening_via_fmriAR(X, Y, plan, run_indices)
    X_white <- whitened$X
    Y_white <- whitened$Y

    # Update residuals if more iterations (use original scale residuals)
    if (iter < n_iter) {
      coef <- tryCatch(
        base::qr.solve(X_white, Y_white),
        error = function(e) qr.coef(qr(X_white, LAPACK = TRUE), Y_white)
      )
      residuals <- Y - X %*% coef
    }
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

  list(
    plan = plan,
    X_white = X_white,
    Y_white = Y_white,
    ar_coef = ar_coef
  )
}

#' Get AR order from plan or config
#'
#' @param plan fmriAR_plan object or NULL
#' @param cfg fmri_lm_config object
#' @return Integer AR order
#' @keywords internal
#' @noRd
.get_ar_order <- function(plan = NULL, cfg = NULL) {
  if (!is.null(plan) && !is.null(plan$order)) {
    return(plan$order[["p"]] %||% 0L)
  }

  if (!is.null(cfg)) {
    ar_opts <- if (inherits(cfg, "fmri_lm_config")) cfg$ar else cfg
    ar_opts <- .normalize_ar_options(ar_opts)
    return(switch(ar_opts$struct,
      "ar1" = 1L,
      "ar2" = 2L,
      "ar3" = 3L,
      "ar4" = 4L,
      "arp" = as.integer(ar_opts$p %||% 0L),
      0L
    ))
  }

  0L
}

#' Compute effective degrees of freedom for AR models
#'
#' Maintains compatibility with existing fmrireg approach
#'
#' @param n Sample size
#' @param p Number of parameters
#' @param plan fmriAR_plan object or NULL
#' @param ar_coef Alternative AR coefficients if plan not available
#' @return Effective degrees of freedom
#' @keywords internal
#' @noRd
.compute_ar_effective_df_compat <- function(n, p, plan = NULL, ar_coef = NULL) {
  # Extract phi coefficients
  phi <- NULL

  if (!is.null(plan)) {
    if (!is.null(plan$phi) && length(plan$phi) > 0) {
      # Average across runs if multiple
      phi_list <- plan$phi
      if (length(phi_list) == 1) {
        phi <- phi_list[[1]]
      } else {
        # Pool across runs
        max_len <- max(lengths(phi_list))
        phi_mat <- sapply(phi_list, function(x) c(x, rep(0, max_len - length(x))))
        phi <- rowMeans(phi_mat)
      }
    }
  } else if (!is.null(ar_coef)) {
    if (is.list(ar_coef)) {
      max_len <- max(lengths(ar_coef))
      phi_mat <- sapply(ar_coef, function(x) c(x, rep(0, max_len - length(x))))
      phi <- rowMeans(phi_mat)
    } else {
      phi <- ar_coef
    }
  }

  if (is.null(phi) || length(phi) == 0) {
    return(n - p)
  }

  # Effective sample size via lag-correlation inflation:
  # n_eff = n / (1 + 2 * sum_{k>=1} (1-k/n) rho_k)
  rho <- tryCatch(
    {
      if (length(phi) == 1L) {
        phi[1]^(seq_len(n - 1L))
      } else {
        stats::ARMAacf(ar = phi, ma = numeric(0), lag.max = n - 1L)[-1]
      }
    },
    error = function(e) NULL
  )
  if (is.null(rho) || length(rho) == 0) {
    return(n - p)
  }
  k <- seq_along(rho)
  denom <- 1 + 2 * sum((1 - k / n) * rho)
  if (!is.finite(denom) || denom <= 0) {
    return(n - p)
  }
  effective_n <- n / denom
  effective_n <- min(max(effective_n, p + 1), n)

  max(effective_n - p, 1)
}
