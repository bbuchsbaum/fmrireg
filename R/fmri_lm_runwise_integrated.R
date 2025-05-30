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
    # Store the fixed phi for the solver to use
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
  if (!is.null(result$sigma_robust)) {
    sigma_vec <- rep(result$sigma_robust, ncol(Y_run))
  } else {
    # Standard sigma calculation
    dfres <- nrow(X_run) - ncol(X_run)
    if (!is.null(result$effective_df)) {
      dfres <- result$effective_df
    }
    sigma_vec <- sqrt(rss / dfres)
  }
  
  # Compute contrasts if weights provided
  contrasts_result <- NULL
  if (!is.null(conlist_weights) && length(conlist_weights) > 0) {
    # Use the existing contrast computation function
    if (exists("fit_lm_contrasts_fast", mode = "function")) {
      contrasts_result <- fit_lm_contrasts_fast(
        betas = betas,
        sigma2 = mean(sigma_vec^2),
        XtXinv = XtXinv,
        conlist = conlist_weights,
        fcon = fconlist_weights,
        dfres = result$df_residual %||% (nrow(X_run) - ncol(X_run)),
        robust_weights = result$robust_weights,
        ar_order = if (!is.null(result$phi_hat)) length(result$phi_hat) else 0
      )
    }
  }
  
  # Build return structure compatible with runwise_lm
  list(
    betas = betas,
    sigma = sigma_vec,
    rss = rss,  # Already calculated above
    dfres = result$df_residual %||% (nrow(X_run) - ncol(X_run)),
    XtXinv = XtXinv,
    contrasts = contrasts_result,
    phi_hat = result$phi_hat,  # AR parameters if estimated
    robust_weights = result$robust_weights  # Robust weights if used
  )
}