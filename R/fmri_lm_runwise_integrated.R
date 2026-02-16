#' Process a single run using the integrated solver
#'
#' This function replaces the complex manual AR+Robust logic in runwise_lm
#' with a call to the integrated solver.
#'
#' @keywords internal
#' @noRd
process_run_integrated <- function(X_run, Y_run, cfg, phi_fixed = NULL, 
                                   sigma_fixed = NULL, conlist_weights = NULL,
                                   fconlist_weights = NULL, vnames = NULL) {
  
  # Debug message
  if (getOption("fmrireg.debug", FALSE)) {
    message("INTEGRATED SOLVER CALLED: AR=", cfg$ar$struct, ", Robust=", cfg$robust$type)
    message("  X_run dim: ", paste(dim(X_run), collapse="x"))
    message("  Y_run dim: ", paste(dim(Y_run), collapse="x"))
  }
  
  # Create a temporary config for this run
  run_cfg <- cfg
  
  # If we have fixed parameters, inject them
  if (!is.null(phi_fixed)) {
    # Store fixed phi in canonical field consumed by integrated AR pipeline.
    run_cfg$ar$phi <- phi_fixed
    # Keep legacy alias for compatibility with older callers.
    run_cfg$ar$phi_fixed <- phi_fixed
  }
  
  if (!is.null(sigma_fixed) && run_cfg$robust$scale_scope == "global") {
    run_cfg$robust$sigma_fixed <- sigma_fixed
  }
  
  # Call the integrated solver
  # Run indices not needed for single run
  result <- solve_integrated_glm(X = X_run, Y = Y_run, config = run_cfg, 
                                 run_indices = NULL)
  
  # Extract what we need for runwise_lm compatibility
  betas <- result$coefficients %||% result$betas
  
  # Get XtXinv - it should be in the result
  XtXinv <- result$XtXinv
  
  if (getOption("fmrireg.debug", FALSE)) {
    message("  Result names: ", paste(names(result), collapse=", "))
    message("  betas class: ", class(betas), " dim: ", paste(dim(betas), collapse="x"))
  }
  
  # Ensure betas is a matrix
  if (!is.matrix(betas)) {
    betas <- as.matrix(betas)
  }

  dfres <- result$dfres %||% result$df_residual %||% (nrow(X_run) - ncol(X_run))
  
  # Calculate RSS from residuals
  if (!is.null(result$residuals)) {
    resid_mat <- as.matrix(result$residuals)
    rss <- colSums(resid_mat^2)
  } else {
    # If residuals not provided, calculate them
    fitted <- X_run %*% betas
    resid_mat <- Y_run - fitted
    rss <- colSums(resid_mat^2)
  }
  
  # Calculate sigma - use robust scale if available
  sigma2_vec <- result$sigma2
  if (is.null(sigma2_vec) && !is.null(result$sigma_robust)) {
    sigma2_vec <- rep(result$sigma_robust^2, ncol(Y_run))
  }

  if (!is.null(sigma2_vec)) {
    sigma2_vec <- as.numeric(sigma2_vec)
    if (length(sigma2_vec) == 1L) {
      sigma2_vec <- rep(sigma2_vec, ncol(Y_run))
    }
    sigma_vec <- sqrt(pmax(0, sigma2_vec))
  } else {
    sigma_vec <- sqrt(pmax(0, rss / dfres))
    sigma2_vec <- sigma_vec^2
  }
  
  # Compute contrasts if weights provided
  contrasts_result <- NULL
  if (!is.null(conlist_weights) && length(conlist_weights) > 0) {
    contrasts_result <- fit_lm_contrasts_fast(
      B = betas,
      sigma2 = sigma2_vec,
      XtXinv = XtXinv,
      conlist = conlist_weights,
      fconlist = fconlist_weights,
      df = dfres,
      robust_weights = result$robust_weights,
      ar_order = if (!is.null(result$ar_order)) {
        as.integer(result$ar_order)
      } else if (!is.null(result$phi_hat) && length(result$phi_hat) > 0L) {
        if (is.list(result$phi_hat)) length(result$phi_hat[[1]]) else length(result$phi_hat)
      } else {
        0L
      }
    )
  }
  
  # Build return structure compatible with runwise_lm
  list(
    betas = betas,
    sigma = sigma_vec,
    rss = rss,  # Already calculated above
    dfres = dfres,
    XtXinv = XtXinv,
    contrasts = contrasts_result,
    phi_hat = result$phi_hat %||% result$ar_coef,  # AR parameters if estimated/fixed
    robust_weights = result$robust_weights  # Robust weights if used
  )
}
