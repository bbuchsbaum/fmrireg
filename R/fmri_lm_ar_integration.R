#' AR Integration for GLM
#'
#' Functions to integrate AR modeling with the GLM solver pipeline
#' Now delegates to fmriAR package for all AR operations
#'
#' @keywords internal
#' @noRd
NULL

#' Apply AR whitening to GLM context
#'
#' @description
#' Transforms the design matrix X and response Y in a GLM context
#' using estimated or provided AR coefficients. Delegates to fmriAR
#' for all AR operations.
#'
#' @param glm_ctx A glm_context object
#' @param ar_options List with AR configuration:
#'   - cor_struct: "ar1", "ar2", etc. or "none"
#'   - iter: Number of AR estimation iterations
#'   - phi: Optional pre-specified AR coefficients
#' @param run_indices List of indices for each run
#'
#' @return Updated glm_context with whitened matrices and AR info
#' @keywords internal
#' @noRd
whiten_glm_context <- function(glm_ctx, ar_options, run_indices = NULL) {
  stopifnot(is.glm_context(glm_ctx))

  ar_opts <- if (is.null(ar_options)) NULL else .normalize_ar_options(ar_options)
  ar_struct <- if (is.null(ar_opts)) "iid" else ar_opts$struct

  if (is.null(ar_opts) || identical(ar_struct, "iid")) {
    return(glm_ctx)
  }

  X <- glm_ctx$X
  Y <- glm_ctx$Y

  if (anyNA(X) || anyNA(Y)) {
    stop("NA values detected in X or Y", call. = FALSE)
  }

  # Get residuals for AR estimation
  if (is.null(glm_ctx$residuals)) {
    if (!is.null(glm_ctx$proj) && !is.null(glm_ctx$proj$Pinv)) {
      coef <- glm_ctx$proj$Pinv %*% Y
    } else {
      coef <- tryCatch(
        qr.coef(qr(X, LAPACK = TRUE), Y),
        error = function(e) base::qr.solve(X, Y)
      )
    }
    residuals <- Y - X %*% coef
  } else {
    residuals <- glm_ctx$residuals
  }

  # If phi provided, use it directly
  if (!is.null(ar_opts$phi)) {
    phi_input <- ar_opts$phi
    if (!is.list(phi_input)) {
      phi_input <- list(phi_input)
    }
    run_count <- if (is.null(run_indices)) 1L else length(run_indices)
    if (run_count > 1L && length(phi_input) == 1L && !isTRUE(ar_opts$global)) {
      phi_input <- rep(phi_input, run_count)
    }
    theta_input <- if (!is.null(ar_opts$theta)) {
      if (is.list(ar_opts$theta)) {
        ar_opts$theta
      } else {
        list(ar_opts$theta)
      }
    } else {
      replicate(length(phi_input), numeric(0), simplify = FALSE)
    }
    if (run_count > 1L && length(theta_input) == 1L && !isTRUE(ar_opts$global)) {
      theta_input <- rep(theta_input, run_count)
    }
    # Create plan from provided coefficients
    runs <- NULL
    if (!is.null(run_indices)) {
      n <- nrow(X)
      runs <- integer(n)
      for (i in seq_along(run_indices)) {
        runs[run_indices[[i]]] <- i
      }
    }

    pooling_mode <- if (isTRUE(ar_opts$global) || run_count <= 1L) "global" else "run"
    plan <- fmriAR::compat$plan_from_phi(
      phi = phi_input,
      theta = theta_input,
      runs = runs,
      pooling = pooling_mode,
      exact_first = ar_opts$exact_first %||% TRUE
    )
  } else {
    # Estimate AR using fmriAR
    plan <- .estimate_ar_via_fmriAR(residuals, ar_opts, run_indices)
  }

  # Apply whitening
  whitened <- .apply_ar_whitening_via_fmriAR(X, Y, plan, run_indices)

  # Extract phi for compatibility
  phi_list <- if (!is.null(plan$phi)) {
    plan$phi
  } else if (!is.null(plan$phi_by_parcel)) {
    # Average across parcels
    list(rowMeans(sapply(plan$phi_by_parcel,
                        function(x) c(x, rep(0, max(lengths(plan$phi_by_parcel)) - length(x))))))
  } else {
    list(numeric(0))
  }

  # Update projection for whitened X
  proj_new <- .fast_preproject(whitened$X)

  # Create new context with whitened data
  glm_context(
    X = whitened$X,
    Y = whitened$Y,
    proj = proj_new,
    phi_hat = phi_list,
    sigma_robust_scale = glm_ctx$sigma_robust_scale,
    robust_weights = glm_ctx$robust_weights
  )
}

#' Iterative AR estimation and whitening
#'
#' @description
#' Performs iterative AR parameter estimation and whitening using fmriAR.
#' Alternates between estimating AR parameters from residuals
#' and re-fitting the whitened model.
#'
#' @param glm_ctx Initial GLM context
#' @param ar_options AR configuration options
#' @param run_indices Run structure
#' @param max_iter Maximum iterations
#' @param tol Convergence tolerance for AR parameters
#'
#' @return List with final results including AR coefficients
#' @keywords internal
#' @noRd
iterative_ar_solve <- function(glm_ctx, ar_options, run_indices = NULL,
                               max_iter = NULL, tol = 1e-4) {

  ar_opts <- if (is.null(ar_options)) NULL else .normalize_ar_options(ar_options)

  if (is.null(max_iter)) {
    max_iter <- if (is.null(ar_opts)) 1 else ar_opts$iter_gls %||% ar_opts$iter %||% 1
  }

  ar_struct <- if (is.null(ar_opts)) "iid" else ar_opts$struct

  # Initial OLS fit
  result <- solve_glm_core(glm_ctx, return_fitted = TRUE)

  if (is.null(ar_opts) || identical(ar_struct, "iid") || max_iter == 0) {
    return(result)
  }

  # Use fmriAR for iterative AR-GLS
  ar_result <- .iterative_ar_gls_via_fmriAR(
    X = glm_ctx$X,
    Y = glm_ctx$Y,
    cfg = ar_opts,
    run_indices = run_indices,
    max_iter = max_iter
  )

  # Solve final whitened system
  proj_white <- .fast_preproject(ar_result$X_white)
  glm_ctx_white <- glm_context(
    X = ar_result$X_white,
    Y = ar_result$Y_white,
    proj = proj_white
  )

  result <- solve_glm_core(glm_ctx_white, return_fitted = TRUE)

  result$XtXinv <- proj_white$XtXinv

  # Add AR info
  result$ar_coef <- ar_result$ar_coef
  result$phi_hat <- ar_result$ar_coef
  result$ar_order <- .get_ar_order(ar_result$plan, ar_opts)
  result$ar_plan <- ar_result$plan  # Store for downstream use

  result
}

#' Compute effective degrees of freedom for AR models
#'
#' @description
#' Adjusts degrees of freedom to account for autocorrelation
#' in the residuals. Uses Satterthwaite-type approximation.
#' Now delegates to fmriAR adapter for consistency.
#'
#' @param n Sample size
#' @param p Number of parameters
#' @param phi AR coefficients (can be list or vector)
#' @param plan Optional fmriAR_plan object
#' @param n_runs Number of runs (for compatibility)
#' @param penalize_ar Whether to penalize for AR estimation
#'
#' @return Effective degrees of freedom
#' @keywords internal
#' @noRd
compute_ar_effective_df <- function(n, p, phi = NULL, plan = NULL,
                                   n_runs = 1, penalize_ar = FALSE) {
  # Delegate to adapter for consistency
  df <- .compute_ar_effective_df_compat(n, p, plan, phi)

  # Optional penalty for AR estimation (not standard)
  if (penalize_ar && (!is.null(phi) || !is.null(plan))) {
    ar_order <- if (!is.null(plan)) {
      .get_ar_order(plan)
    } else if (is.list(phi)) {
      length(phi[[1]])
    } else {
      length(phi)
    }
    df <- df - ar_order
  }

  max(df, 1)
}
