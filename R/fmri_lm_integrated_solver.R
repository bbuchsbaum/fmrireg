#' Integrated GLM Solver with AR and Robust Options
#'
#' Main entry point for solving GLM with optional AR whitening and robust fitting
#'
#' @keywords internal
#' @noRd
NULL

#' Solve GLM with integrated AR and robust options
#'
#' @description
#' Coordinates the complete GLM solving pipeline including:
#' - Initial OLS fit
#' - Optional AR whitening (iterative)
#' - Optional robust fitting (IRLS)
#' - Combined AR+Robust ("whiten then robustly weight")
#'
#' @param X Design matrix
#' @param Y Response matrix
#' @param config fmri_lm_control object with all options
#' @param run_indices List of indices for each run
#' @param fixed Optional noise model estimated beforehand, as returned by
#'   `.integrated_noise_estimate()`: `ar_plan` (an fmriAR plan) and/or
#'   `robust_weights`. When supplied, the AR model is not re-estimated and the
#'   robust weights are not re-iterated; `Y` is fitted with them, so a voxel's
#'   result does not depend on which other voxels are in `Y`.
#'
#' @return List with full results including coefficients, residuals,
#'   standard errors, AR parameters, robust weights, etc.
#' @keywords internal
#' @noRd
solve_integrated_glm <- function(X, Y, config, run_indices = NULL,
                                 fixed = NULL) {
  stopifnot(inherits(config, "fmri_lm_control"))

  # Initial projection and context
  proj <- .fast_preproject(X)
  glm_ctx <- glm_context(X = X, Y = Y, proj = proj)

  # Determine solving strategy and extract ar_options in the format
  # iterative_ar_solve() expects.
  method <- .integrated_method(config)
  ar_options <- .integrated_ar_options(config)

  result <- switch(method,
    "ols" = solve_ols_pipeline(glm_ctx),
    "ar" = solve_ar_pipeline(glm_ctx, ar_options, run_indices,
                             ar_plan = fixed$ar_plan),
    "robust" = solve_robust_pipeline(glm_ctx, config$robust,
                                     weights_fixed = fixed$robust_weights),
    "ar_robust" = solve_ar_robust_pipeline(glm_ctx, config, run_indices,
                                           ar_plan = fixed$ar_plan,
                                           weights_fixed = fixed$robust_weights),
    stop("Unknown method: ", method)
  )

  # Add method info
  result$method <- method
  result$config <- config

  # Compute standard errors
  result <- compute_standard_errors(result, glm_ctx)

  result
}

#' Solving strategy of the integrated solver
#' @keywords internal
#' @noRd
.integrated_method <- function(config) {
  method <- config$method %||% "auto"
  if (method == "auto") {
    # Determine method based on options
    # Handle both config structures (ar_options vs ar)
    ar_struct <- config$ar$struct %||%
                 (config$ar_options$cor_struct %||% "iid")
    robust_type <- config$robust$type %||% config$robust %||% FALSE
    if (identical(robust_type, "none")) robust_type <- FALSE

    temporal_noise <- (ar_struct != "iid" && ar_struct != "none") ||
      as.integer(config$ar$q %||% 0L) > 0L
    if (!isFALSE(robust_type) && temporal_noise) {
      method <- "ar_robust"
    } else if (!isFALSE(robust_type)) {
      method <- "robust"
    } else if (temporal_noise) {
      method <- "ar"
    } else {
      method <- "ols"
    }
  }
  method
}

#' AR options of the integrated solver, in iterative_ar_solve() format
#' @keywords internal
#' @noRd
.integrated_ar_options <- function(config) {
  if (!is.null(config$ar)) {
    list(
      cor_struct = config$ar$struct,
      iter = config$ar$iter_gls,
      exact_first = config$ar$exact_first,
      p = config$ar$p,
      q = config$ar$q %||% 0L,
      global = config$ar$global %||% FALSE,
      shared_estimator = config$ar$shared_estimator %||% "pooled_acvf",
      phi = config$ar$phi %||% config$ar$phi_fixed,
      theta = config$ar$theta %||% config$ar$theta_fixed
    )
  } else {
    config$ar_options
  }
}

#' Robust options of the integrated solver
#' @keywords internal
#' @noRd
.integrated_robust_options <- function(robust_options) {
  if (is.character(robust_options)) {
    list(
      type = robust_options,
      k_huber = 1.345,
      c_tukey = 4.685,
      max_iter = 10,
      scale_scope = "voxel"
    )
  } else {
    robust_options
  }
}

#' Estimate the integrated solver's noise model once over all voxels
#'
#' The AR plan and the robust weights of `solve_integrated_glm()` are shared
#' by all response columns, and estimating them pools over the columns of
#' `Y`: the AR plan uses (at most `fmrireg.ar.estimation_max_vox`, default 64)
#' evenly spaced columns, the robust weights use row medians over all
#' columns. Fitting voxel chunks separately therefore made the reference
#' backend's results depend on the chunking. This function estimates both
#' exactly as one `solve_integrated_glm()` call over the whole dataset does,
#' for `solve_integrated_glm(fixed = )`.
#'
#' Only the AR estimation columns are read for the AR plan. The robust
#' weights depend on every voxel's residuals, so a robust model reads all
#' columns once.
#'
#' @param X Design matrix.
#' @param config fmri_lm_control object.
#' @param run_indices List of run row indices.
#' @param n_vox Total number of response columns.
#' @param read_columns Function of a column index vector (`NULL` for all)
#'   returning those columns of the data as a matrix.
#' @return `NULL` for OLS, else a list with `ar_plan` and `robust_weights`.
#' @keywords internal
#' @noRd
.integrated_noise_estimate <- function(X, config, run_indices, n_vox,
                                       read_columns) {
  method <- .integrated_method(config)
  if (identical(method, "ols")) return(NULL)

  plan <- NULL
  if (method %in% c("ar", "ar_robust")) {
    # Mirrors the estimation branch of iterative_ar_solve().
    ar_opts <- .normalize_ar_options(.integrated_ar_options(config))
    max_iter <- ar_opts$iter_gls %||% ar_opts$iter %||% 1
    if (.temporal_noise_enabled(ar_opts) && max_iter != 0 &&
        is.null(ar_opts$phi)) {
      cols <- .ar_estimation_columns(n_vox, ar_opts)
      plan <- .estimate_ar_gls_plan(
        X, as.matrix(read_columns(cols)), ar_opts,
        run_indices = run_indices, max_iter = max_iter,
        tol = ar_opts$tol %||% 5e-3
      )$plan
    }
  }

  weights <- NULL
  if (method %in% c("robust", "ar_robust")) {
    full <- solve_integrated_glm(X, as.matrix(read_columns(NULL)), config,
                                 run_indices, fixed = list(ar_plan = plan))
    weights <- full$robust_weights
  }
  list(ar_plan = plan, robust_weights = weights)
}

#' OLS pipeline
#' @keywords internal
#' @noRd
solve_ols_pipeline <- function(glm_ctx) {
  result <- solve_glm_core(glm_ctx, return_fitted = TRUE)
  result$residuals <- glm_ctx$Y - result$fitted
  # Add XtXinv for contrast computation
  result$XtXinv <- glm_ctx$proj$XtXinv
  result
}

#' AR pipeline
#' @keywords internal
#' @noRd
solve_ar_pipeline <- function(glm_ctx, ar_options, run_indices,
                              ar_plan = NULL) {
  result <- iterative_ar_solve(glm_ctx, ar_options, run_indices,
                               ar_plan = ar_plan)

  # Effective df: the data were whitened, so the whitened residuals carry
  # n - rank df (calculate_effective_df(); reducing them again for the AR
  # coefficients would count the correction twice).
  if (!is.null(result$ar_plan) || !is.null(result$ar_coef) || !is.null(result$ma_coef)) {
    rank_w <- if (!is.null(result$X_white)) {
      .fast_preproject(result$X_white)$rank
    } else {
      glm_ctx$proj$rank
    }
    result$effective_df <- calculate_effective_df(
      nrow(glm_ctx$X), rank_w, method = "residual"
    )
  }

  # Add XtXinv if not present
  if (is.null(result$XtXinv)) {
    if (!is.null(result$X_white)) {
      result$XtXinv <- .fast_preproject(result$X_white)$XtXinv
    } else {
      warning("AR solver did not return whitened XtXinv; using unwhitened fallback.")
      result$XtXinv <- glm_ctx$proj$XtXinv
    }
  }

  result
}

#' Robust pipeline
#' @keywords internal
#' @noRd
solve_robust_pipeline <- function(glm_ctx, robust_options,
                                  weights_fixed = NULL) {
  robust_config <- .integrated_robust_options(robust_options)

  # Run robust fitting
  robust_result <- robust_iterative_fitter(
    initial_glm_ctx = glm_ctx,
    cfg_robust_options = robust_config,
    X_orig_for_resid = glm_ctx$X,
    sigma_fixed = robust_config$sigma_fixed %||% NULL,
    weights_fixed = weights_fixed
  )

  # WLS scale of the final weighted solve (the one robust scale estimator;
  # see robust_iterative_fitter()).
  sigma2 <- as.numeric(robust_result$sigma2)

  # Format output
  list(
    betas = robust_result$betas_robust,
    sigma2 = sigma2,
    robust_weights = robust_result$robust_weights_final,
    XtXinv = robust_result$XtWXi_final,
    dfres = robust_result$dfres,
    fitted = glm_ctx$X %*% robust_result$betas_robust,
    residuals = glm_ctx$Y - glm_ctx$X %*% robust_result$betas_robust
  )
}

#' AR + Robust pipeline ("whiten then robustly weight")
#' @keywords internal
#' @noRd
solve_ar_robust_pipeline <- function(glm_ctx, config, run_indices,
                                     ar_plan = NULL, weights_fixed = NULL) {
  ar_options <- .integrated_ar_options(config)

  # Step 1: Iterative AR whitening
  ar_result <- iterative_ar_solve(
    glm_ctx,
    ar_options,
    run_indices,
    return_solution = FALSE,
    ar_plan = ar_plan
  )

  # Step 2: Reuse final AR-whitened matrices from iterative_ar_solve.
  # Avoid a redundant second whitening pass.
  proj_white <- .fast_preproject(ar_result$X_white)
  glm_ctx_white <- glm_context(
    X = ar_result$X_white,
    Y = ar_result$Y_white,
    proj = proj_white
  )

  # Step 3: Robust fitting on whitened data
  robust_config <- .integrated_robust_options(config$robust)

  robust_result <- robust_iterative_fitter(
    initial_glm_ctx = glm_ctx_white,
    cfg_robust_options = robust_config,
    X_orig_for_resid = glm_ctx_white$X,
    sigma_fixed = robust_config$sigma_fixed %||% NULL,
    weights_fixed = weights_fixed
  )

  # WLS scale of the final weighted solve (the one robust scale estimator;
  # see robust_iterative_fitter()).
  sigma2 <- as.numeric(robust_result$sigma2)

  # Combine results
  list(
    betas = robust_result$betas_robust,
    sigma2 = sigma2,
    robust_weights = robust_result$robust_weights_final,
    ar_coef = ar_result$ar_coef,
    ma_coef = ar_result$ma_coef,
    ar_plan = ar_result$ar_plan,
    XtXinv = robust_result$XtWXi_final,
    dfres = robust_result$dfres,
    fitted = glm_ctx_white$X %*% robust_result$betas_robust,
    residuals = glm_ctx_white$Y - glm_ctx_white$X %*% robust_result$betas_robust,
    # Whitened, then robustly weighted: the one robust effective df,
    # sum(w) - rank (calculate_effective_df() -> .robust_effective_df()).
    effective_df = calculate_effective_df(
      nrow(glm_ctx$X), robust_result$rank,
      robust_weights = robust_result$robust_weights_final,
      method = "residual"
    )
  )
}

#' Compute standard errors
#' @keywords internal
#' @noRd
compute_standard_errors <- function(result, glm_ctx) {
  # Extract variance-covariance matrix
  if (!is.null(result$XtXinv)) {
    vcov <- result$XtXinv
  } else if (!is.null(glm_ctx$proj$XtXinv)) {
    vcov <- glm_ctx$proj$XtXinv
  } else {
    stop("Cannot compute standard errors: no variance-covariance matrix")
  }

  # Use appropriate sigma2
  if (!is.null(result$sigma2)) {
    if (length(result$sigma2) == 1) {
      # Single sigma2 for all voxels
      result$standard_errors <- sqrt(diag(vcov) * result$sigma2)
    } else {
      # Vector of sigma2 (one per voxel)
      result$standard_errors <- sqrt(diag(vcov) %o% result$sigma2)
    }
  }

  result
}
