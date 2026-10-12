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
    return(list(struct = "iid", cor_struct = "none",
                shared_estimator = "pooled_acvf"))
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

  ar_opts$shared_estimator <- ar_opts$shared_estimator %||% "pooled_acvf"
  if (!ar_opts$shared_estimator %in% c("pooled_acvf", "mean_series")) {
    stop("shared_estimator must be 'pooled_acvf' or 'mean_series'", call. = FALSE)
  }

  ar_opts$struct <- struct
  ar_opts$cor_struct <- cor_struct

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
#' @keywords internal
#' @noRd
.shared_ar_correction_design <- function(residuals, ar_opts, design) {
  if (is.null(design)) return(NULL)
  # OLS residualization is linear across response columns:
  # rowMeans((I - P_X)Y) = (I - P_X)rowMeans(Y).  Consequently the matching
  # design remains the exact residual-forming operator for mean_series too.
  design
}

#' Estimate shared AR parameters under the configured spatial estimand
#' @keywords internal
#' @noRd
.estimate_shared_ar_parameters <- function(residuals, ar_order, ar_opts,
                                           run_indices = NULL, censor = NULL,
                                           design = NULL) {
  correction_design <- .shared_ar_correction_design(residuals, ar_opts, design)
  .estimate_ar_parameters_routed(
    .shared_ar_residual_input(residuals, ar_opts),
    ar_order = ar_order,
    run_indices = run_indices,
    censor = censor,
    design = correction_design,
    ar_options = ar_opts
  )
}

#' Determine whether temporal whitening is requested
#' @keywords internal
#' @noRd
.temporal_noise_enabled <- function(ar_opts) {
  ar_opts <- .normalize_ar_options(ar_opts)
  !identical(ar_opts$struct, "iid") || as.integer(ar_opts$q %||% 0L) > 0L
}

#' AR order of a noise configuration
#'
#' The one mapping from an AR structure to its order, used by every fitting
#' path. Accepts an `fmri_lm_control` object or an AR options list (with
#' `struct` or legacy `cor_struct`). An unknown structure is an error: the
#' former per-file `switch()` copies silently returned order 0 (no whitening)
#' or `NULL` for it.
#'
#' @param ar_opts `fmri_lm_control` object or AR options list.
#' @return Integer AR order (0 for iid).
#' @keywords internal
#' @noRd
.ar_order <- function(ar_opts) {
  if (inherits(ar_opts, "fmri_lm_control")) ar_opts <- ar_opts$ar
  ar_opts <- .normalize_ar_options(ar_opts)
  struct <- as.character(ar_opts$struct)
  ord <- switch(struct,
    "iid" = 0L,
    "ar1" = 1L,
    "ar2" = 2L,
    "ar3" = 3L,
    "ar4" = 4L,
    "arp" = suppressWarnings(as.integer(ar_opts$p)),
    NULL
  )
  if (is.null(ord)) {
    stop(sprintf(paste0("Unknown AR structure '%s'; expected one of 'iid', ",
                        "'ar1', 'ar2', 'ar3', 'ar4' or 'arp'."), struct),
         call. = FALSE)
  }
  if (length(ord) != 1L || is.na(ord) || ord < 0L ||
      (identical(struct, "arp") && ord < 1L)) {
    stop("struct = 'arp' requires `p` to be a positive integer.", call. = FALSE)
  }
  ord
}

#' Determine AR order implied by configuration
#'
#' @param ar_opts Normalized AR options list
#' @return Integer AR order (>= 0); see `.ar_order()`.
#' @keywords internal
#' @noRd
.target_ar_order <- function(ar_opts) {
  .ar_order(ar_opts)
}

#' Whether the first observation gets exact AR scaling
#'
#' The documented default (noise_spec(), fmri_lm_control()) is `FALSE`; a
#' missing value means `FALSE` on every path.
#' @keywords internal
#' @noRd
.ar_exact_first <- function(ar_opts) {
  isTRUE(ar_opts$exact_first)
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

    # If fmriAR plan generation did not yield explicit phi coefficients, default
    # to zero coefficients to keep the fast routed path and avoid per-voxel
    # fallback estimation overhead. Legacy fallback can be re-enabled for
    # diagnostics via option fmrireg.ar.fixed_order_legacy_fallback = TRUE.
    if (!isTRUE(getOption("fmrireg.ar.fixed_order_legacy_fallback", FALSE))) {
      return(rep(0, order))
    }

    # Legacy fallback: aggregate per-voxel fixed-order estimates.
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
  p <- .ar_order(ar_opts)

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
    q = q,
    pooling = pooling,
    exact_first = if (.ar_exact_first(ar_opts)) "ar1" else "none",
    # Number of iterations handled externally
    iter = as.integer(ar_opts$iter_gls %||% ar_opts$iter %||% 1L)
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

  # Apply whitening (censor is passed for consistent segment handling)
  result <- fmriAR::whiten_apply(
    plan = plan,
    X = X,
    Y = Y,
    runs = runs,
    parallel = .ar_parallel_whiten(Y),
    censor = censor
  )

  # Return in expected format
  list(X = result$X, Y = result$Y)
}

#' Columns used to estimate a shared AR model
#'
#' At most `estimation_max_vox` (option `fmrireg.ar.estimation_max_vox`,
#' default 64) evenly spaced response columns.
#' @param n_vox Number of response columns.
#' @param cfg AR options.
#' @keywords internal
#' @noRd
.ar_estimation_columns <- function(n_vox, cfg) {
  est_max_vox <- cfg$estimation_max_vox %||% getOption("fmrireg.ar.estimation_max_vox", 64L)
  est_max_vox <- suppressWarnings(as.integer(est_max_vox))
  if (!is.finite(est_max_vox) || est_max_vox <= 0L) est_max_vox <- n_vox
  if (is.na(est_max_vox)) est_max_vox <- n_vox
  est_n <- min(n_vox, est_max_vox)
  if (est_n >= n_vox) {
    seq_len(n_vox)
  } else {
    unique(as.integer(round(seq(1, n_vox, length.out = est_n))))
  }
}

#' Iterative AR-GLS estimation of a shared AR plan
#'
#' The estimation half of `.iterative_ar_gls_via_fmriAR()`: OLS residuals of
#' the estimation columns, AR estimation, and (ARMA only) GLS refinement.
#' @param X Design matrix.
#' @param Y_est The estimation columns of the data
#'   (`Y[, .ar_estimation_columns(ncol(Y), cfg)]`).
#' @return List with `plan` (`NULL` when no AR model is fitted),
#'   `iterations` and `converged`.
#' @keywords internal
#' @noRd
.estimate_ar_gls_plan <- function(X, Y_est, cfg, run_indices = NULL,
                                  max_iter = NULL, censor = NULL, tol = 5e-3) {
  ar_params <- .fmrireg_to_fmriAR_config(cfg, length(run_indices))
  n_iter <- max_iter %||% ar_params$iter %||% 1L
  # AR residual-bias correction is tied to the initial OLS residual operator.
  # ARMA has no such correction and may continue iterative refinement.
  if (identical(ar_params$method, "ar") && n_iter > 0L) n_iter <- 1L
  min_iter <- as.integer(cfg$min_iter %||% 1L)
  if (!is.finite(min_iter) || min_iter < 1L) min_iter <- 1L

  # Handle no AR case
  if ((ar_params$p == 0L && ar_params$q == 0L) || n_iter == 0L) {
    return(list(plan = NULL, iterations = 0L, converged = NA))
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
    correction_design <- .shared_ar_correction_design(residuals, cfg, X)
    plan <- .estimate_ar_via_fmriAR(
      .shared_ar_residual_input(residuals, cfg),
      cfg, run_indices, censor = censor, design = correction_design
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

  list(plan = plan, iterations = as.integer(iter),
       converged = if (n_iter <= 1L) NA else converged)
}

#' AR and MA coefficient lists of an fmriAR plan
#' @keywords internal
#' @noRd
.ar_plan_coefficients <- function(plan) {
  if (is.null(plan)) return(list(ar_coef = NULL, ma_coef = NULL))
  ar_coef <- if (!is.null(plan$phi)) {
    plan$phi  # List of AR coefficients
  } else if (!is.null(plan$phi_by_parcel)) {
    # Average across parcels if parcel-based
    phi_list <- plan$phi_by_parcel
    list(rowMeans(sapply(phi_list, function(x) c(x, rep(0, max(lengths(phi_list)) - length(x))))))
  } else {
    NULL
  }
  ma_coef <- if (!is.null(plan$theta)) {
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
  list(ar_coef = ar_coef, ma_coef = ma_coef)
}

#' Whether fmriAR::whiten_apply() should use threads
#'
#' `getOption("fmrireg.ar.parallel_whiten", TRUE)` (an invalid value counts
#' as `TRUE`), and never for a response with at most one column: spinning up
#' threads for one series (e.g. per-voxel whitening) costs more than it saves.
#' @param Y The response being whitened.
#' @keywords internal
#' @noRd
.ar_parallel_whiten <- function(Y) {
  parallel_whiten <- getOption("fmrireg.ar.parallel_whiten", TRUE)
  if (!is.logical(parallel_whiten) || length(parallel_whiten) != 1L ||
      is.na(parallel_whiten)) {
    parallel_whiten <- TRUE
  }
  n_cols <- if (is.null(dim(Y))) 1L else ncol(Y)
  isTRUE(parallel_whiten) && n_cols > 1L
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
#' @param plan Optional AR plan estimated beforehand (e.g. once over all
#'   voxels by `.estimate_ar_gls_plan()`); when supplied, no estimation is done
#'   and the data are whitened with it.
#' @return List with AR plan and whitened matrices
#' @keywords internal
#' @noRd
.iterative_ar_gls_via_fmriAR <- function(X, Y, cfg, run_indices = NULL, max_iter = NULL,
                                         censor = NULL, tol = 5e-3, plan = NULL) {
  if (is.null(plan)) {
    est <- .estimate_ar_gls_plan(
      X, Y[, .ar_estimation_columns(ncol(Y), cfg), drop = FALSE], cfg,
      run_indices = run_indices, max_iter = max_iter, censor = censor, tol = tol
    )
    plan <- est$plan
  } else {
    est <- list(iterations = 0L, converged = NA)
  }

  # Final whitening on full data using the last AR plan.
  # If no plan was estimated (no AR model), passthrough original matrices.
  if (is.null(plan)) {
    X_white <- X
    Y_white <- Y
  } else {
    whitened_full <- .apply_ar_whitening_via_fmriAR(X, Y, plan, run_indices, censor = censor)
    X_white <- whitened_full$X
    Y_white <- whitened_full$Y
  }
  coefs <- .ar_plan_coefficients(plan)

  list(
    plan = plan,
    X_white = X_white,
    Y_white = Y_white,
    ar_coef = coefs$ar_coef,
    ma_coef = coefs$ma_coef,
    iterations = est$iterations,
    converged = est$converged,
    censor = censor
  )
}

#' Get AR order from plan or config
#'
#' @param plan fmriAR_plan object or NULL
#' @param cfg fmri_lm_control object
#' @return Integer AR order
#' @keywords internal
#' @noRd
.get_ar_order <- function(plan = NULL, cfg = NULL) {
  if (!is.null(plan) && !is.null(plan$order)) {
    return(plan$order[["p"]] %||% 0L)
  }

  if (!is.null(cfg)) {
    return(.ar_order(cfg))
  }

  0L
}
