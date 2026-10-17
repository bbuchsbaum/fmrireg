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
#' @param config fmri_lm_control object with all options. `config$ar$phi_fixed`
#'   holds the AR coefficients fixed instead of estimating them.
#' @param run_indices List of indices for each run
#'
#' @return List with `betas`, `sigma2`, `XtXinv` (in the fitted, whitened
#'   and/or weighted domain), `dfres` and, as applicable, `rss`, `ar_coef`,
#'   `ma_coef`, `phi_hat`, `ar_order`, `ar_plan` and `robust_weights`.
#' @keywords internal
#' @noRd
solve_integrated_glm <- function(X, Y, config, run_indices = NULL) {
  stopifnot(inherits(config, "fmri_lm_control"))

  glm_ctx <- glm_context(X = X, Y = Y, proj = .fast_preproject(X))
  robust <- .fmri_lm_robust_enabled(config$robust)
  temporal <- .temporal_noise_enabled(config$ar)

  if (robust && temporal) {
    solve_ar_robust_pipeline(glm_ctx, config, run_indices)
  } else if (robust) {
    solve_robust_pipeline(glm_ctx, config$robust)
  } else if (temporal) {
    solve_ar_pipeline(glm_ctx, config$ar, run_indices)
  } else {
    solve_ols_pipeline(glm_ctx)
  }
}

#' OLS pipeline
#' @keywords internal
#' @noRd
solve_ols_pipeline <- function(glm_ctx) {
  result <- solve_glm_core(glm_ctx, return_fitted = TRUE)
  result$fitted <- NULL
  result$XtXinv <- glm_ctx$proj$XtXinv
  result
}

#' AR pipeline
#' @keywords internal
#' @noRd
solve_ar_pipeline <- function(glm_ctx, ar_options, run_indices) {
  result <- iterative_ar_solve(glm_ctx, ar_options, run_indices)

  if (is.null(result$XtXinv)) {
    warning("AR solver did not return whitened XtXinv; using unwhitened fallback.")
    result$XtXinv <- glm_ctx$proj$XtXinv
  }
  # X_white/Y_white stay: the latent AR fit forms its whitened residual
  # covariance from them (standard_error(<fmri_latent_lm>, recon = TRUE)).
  result$fitted <- NULL
  result
}

#' Robust pipeline
#' @keywords internal
#' @noRd
solve_robust_pipeline <- function(glm_ctx, robust_options) {
  robust_result <- robust_iterative_fitter(
    initial_glm_ctx = glm_ctx,
    cfg_robust_options = robust_options,
    X_orig_for_resid = glm_ctx$X
  )

  # WLS scale of the final weighted solve (the one robust scale estimator;
  # see robust_iterative_fitter()).
  list(
    betas = robust_result$betas_robust,
    sigma2 = as.numeric(robust_result$sigma2),
    robust_weights = robust_result$robust_weights_final,
    XtXinv = robust_result$XtWXi_final,
    dfres = robust_result$dfres
  )
}

#' AR + Robust pipeline ("whiten then robustly weight")
#' @keywords internal
#' @noRd
solve_ar_robust_pipeline <- function(glm_ctx, config, run_indices) {
  # Step 1: AR estimation and whitening; reuse the whitened matrices.
  ar_result <- iterative_ar_solve(
    glm_ctx,
    config$ar,
    run_indices,
    return_solution = FALSE
  )

  glm_ctx_white <- glm_context(
    X = ar_result$X_white,
    Y = ar_result$Y_white,
    proj = .fast_preproject(ar_result$X_white)
  )

  # Step 2: Robust fitting on whitened data
  robust_result <- robust_iterative_fitter(
    initial_glm_ctx = glm_ctx_white,
    cfg_robust_options = config$robust,
    X_orig_for_resid = glm_ctx_white$X
  )

  # WLS scale of the final weighted solve (the one robust scale estimator;
  # see robust_iterative_fitter()).
  list(
    betas = robust_result$betas_robust,
    sigma2 = as.numeric(robust_result$sigma2),
    robust_weights = robust_result$robust_weights_final,
    ar_coef = ar_result$ar_coef,
    ma_coef = ar_result$ma_coef,
    ar_plan = ar_result$ar_plan,
    XtXinv = robust_result$XtWXi_final,
    dfres = robust_result$dfres
  )
}
