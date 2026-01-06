#' Robust Iterative Fitting Engine
#'
#' Performs Iteratively Reweighted Least Squares (IRLS) for robust regression
#' using either Huber or Tukey bisquare weights. This function takes an initial
#' GLM context (typically from OLS on original or whitened data) and iteratively
#' refines the fit using robust weights.
#'
#' @param initial_glm_ctx A \code{glm_context} object containing initial fit results
#' @param cfg_robust_options List of robust fitting options from \code{fmri_lm_config$robust}
#' @param X_orig_for_resid The design matrix corresponding to Y in initial_glm_ctx
#'   before any robust weighting (needed for residual calculation)
#' @param sigma_fixed Optional fixed sigma value (for global scale estimation)
#'
#' @return A list containing:
#' \describe{
#'   \item{betas_robust}{Final robust coefficient estimates}
#'   \item{XtWXi_final}{(X'WX)^-1 from the final weighted iteration}
#'   \item{sigma_robust_scale_final}{Final robust scale estimate}
#'   \item{robust_weights_final}{Final robust weights}
#'   \item{dfres}{Residual degrees of freedom}
#' }
#'
#' @details
#' The IRLS algorithm:
#' 1. Calculate residuals using X_orig_for_resid and current betas
#' 2. Estimate robust scale (MAD) unless sigma_fixed is provided
#' 3. Calculate robust weights based on scaled residuals
#' 4. Create weighted GLM context and solve
#' 5. Iterate until max_iter reached
#'
#' @keywords internal
#' @importFrom matrixStats rowMedians
#' @noRd
robust_iterative_fitter <- function(initial_glm_ctx, 
                                    cfg_robust_options,
                                    X_orig_for_resid,
                                    sigma_fixed = NULL) {
  
  # Validate inputs
  stopifnot(is.glm_context(initial_glm_ctx))
  stopifnot(is.list(cfg_robust_options))
  stopifnot(!is.null(X_orig_for_resid))
  
  # Extract options
  psi_type <- cfg_robust_options$type
  if (isFALSE(psi_type)) {
    stop("robust_iterative_fitter called with robust type = FALSE")
  }
  
  k_huber <- cfg_robust_options$k_huber
  c_tukey <- cfg_robust_options$c_tukey
  max_iter <- cfg_robust_options$max_iter
  scale_scope <- cfg_robust_options$scale_scope
  
  # Ensure matrices
  if (!is.matrix(X_orig_for_resid)) X_orig_for_resid <- as.matrix(X_orig_for_resid)
  Y <- initial_glm_ctx$Y
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  
  # Check for NAs
  if (anyNA(Y) || anyNA(X_orig_for_resid)) {
    stop("NA values detected in 'X_orig_for_resid' or 'Y' for robust_iterative_fitter")
  }
  
  # Initial fit from context (could be OLS or whitened OLS)
  initial_fit <- solve_glm_core(initial_glm_ctx, return_fitted = TRUE)
  betas <- initial_fit$betas
  
  # Initialize final values
  XtWXi_final <- initial_glm_ctx$proj$XtXinv
  dfres <- initial_glm_ctx$proj$dfres
  
  # IRLS iterations
  for (it in seq_len(max_iter)) {
    # Calculate residuals using original (unweighted) X
    resid <- Y - X_orig_for_resid %*% betas
    
    # Calculate row-wise median absolute residuals
    row_med <- matrixStats::rowMedians(abs(resid))
    
    # Estimate or use fixed sigma
    if (is.null(sigma_fixed)) {
      sigma_hat <- 1.4826 * median(row_med)  # MAD estimate
      if (sigma_hat <= .Machine$double.eps) sigma_hat <- .Machine$double.eps
    } else {
      sigma_hat <- sigma_fixed
    }
    
    # Scaled residuals
    u <- row_med / sigma_hat
    
    # Calculate weights based on psi function
    # Protect against division by zero when residuals are exactly 0
    u_safe <- pmax(abs(u), .Machine$double.eps)
    w <- switch(psi_type,
                huber = pmin(1, k_huber / u_safe),
                bisquare = ifelse(abs(u) <= c_tukey,
                                  (1 - (u / c_tukey)^2)^2,
                                  0))
    
    # Apply sqrt of weights for weighted least squares
    sqrtw <- sqrt(w)
    Xw <- X_orig_for_resid * sqrtw
    Yw <- sweep(Y, 1, sqrtw, `*`)
    
    # Create weighted projection
    proj_w <- .fast_preproject(Xw)
    
    # Create weighted GLM context
    glm_ctx_weighted <- glm_context(
      X = Xw,
      Y = Yw,
      proj = proj_w,
      robust_weights = w  # Store weights in context
    )
    
    # Solve weighted least squares
    fit_w <- solve_glm_core(glm_ctx_weighted, return_fitted = (it < max_iter))
    
    # Update betas and XtWXi
    betas <- fit_w$betas
    XtWXi_final <- proj_w$XtXinv
  }
  
  # Final scale estimation
  resid_final <- Y - X_orig_for_resid %*% betas
  row_med_final <- matrixStats::rowMedians(abs(resid_final))
  sigma_robust <- 1.4826 * median(row_med_final)
  if (sigma_robust <= .Machine$double.eps) sigma_robust <- .Machine$double.eps
  
  # Final weights for output
  u_final <- row_med_final / sigma_robust
  # Protect against division by zero when residuals are exactly 0
  u_final_safe <- pmax(abs(u_final), .Machine$double.eps)
  w_final <- switch(psi_type,
                    huber = pmin(1, k_huber / u_final_safe),
                    bisquare = ifelse(abs(u_final) <= c_tukey,
                                      (1 - (u_final / c_tukey)^2)^2,
                                      0))
  
  # Return results
  list(
    betas_robust = betas,
    XtWXi_final = XtWXi_final,
    sigma_robust_scale_final = sigma_robust,
    robust_weights_final = w_final,
    dfres = dfres
  )
}