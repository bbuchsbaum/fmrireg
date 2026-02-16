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
  tol <- as.numeric(cfg_robust_options$tol %||% 1e-4)
  min_iter <- as.integer(cfg_robust_options$min_iter %||% 1L)
  if (!is.finite(tol) || tol < 0) tol <- 1e-4
  if (!is.finite(min_iter) || min_iter < 1L) min_iter <- 1L
  scope_opt <- cfg_robust_options$scale_scope
  if (is.null(scope_opt)) scope_opt <- "run"
  scale_scope <- tolower(as.character(scope_opt))
  if (identical(scale_scope, "local")) scale_scope <- "voxel"
  if (!(scale_scope %in% c("run", "global", "voxel"))) {
    scale_scope <- "run"
  }
  
  # Ensure matrices
  if (!is.matrix(X_orig_for_resid)) X_orig_for_resid <- as.matrix(X_orig_for_resid)
  Y <- initial_glm_ctx$Y
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  
  # Check for NAs
  if (anyNA(Y) || anyNA(X_orig_for_resid)) {
    stop("NA values detected in 'X_orig_for_resid' or 'Y' for robust_iterative_fitter")
  }
  
  # Initial fit from context (could be OLS or whitened OLS)
  initial_fit <- solve_glm_core(initial_glm_ctx, return_fitted = FALSE)
  betas <- initial_fit$betas
  
  # Initialize final values
  XtWXi_final <- initial_glm_ctx$proj$XtXinv
  dfres <- initial_glm_ctx$proj$dfres

  compute_scale <- function(resid_mat, scope, sigma_fixed_value = NULL) {
    eps <- .Machine$double.eps
    abs_resid <- abs(resid_mat)

    if (!is.null(sigma_fixed_value)) {
      sigma_vec <- as.numeric(sigma_fixed_value)
      if (length(sigma_vec) == 1L) {
        sigma_vec <- rep(sigma_vec, ncol(abs_resid))
      } else if (length(sigma_vec) != ncol(abs_resid)) {
        stop("sigma_fixed must have length 1 or match number of response columns")
      }
      sigma_vec[!is.finite(sigma_vec) | sigma_vec <= eps] <- eps
    } else if (identical(scope, "voxel")) {
      sigma_vec <- 1.4826 * matrixStats::colMedians(abs_resid)
      sigma_vec[!is.finite(sigma_vec) | sigma_vec <= eps] <- eps
    } else {
      row_med <- matrixStats::rowMedians(abs_resid)
      sigma_scalar <- 1.4826 * stats::median(row_med)
      if (!is.finite(sigma_scalar) || sigma_scalar <= eps) {
        sigma_scalar <- eps
      }
      sigma_vec <- rep(sigma_scalar, ncol(abs_resid))
    }

    u <- matrixStats::rowMedians(sweep(abs_resid, 2, sigma_vec, "/"))
    sigma_scalar <- stats::median(sigma_vec)
    if (!is.finite(sigma_scalar) || sigma_scalar <= eps) {
      sigma_scalar <- eps
    }

    list(u = u, sigma_scalar = sigma_scalar, sigma_components = sigma_vec)
  }
  
  # IRLS iterations
  for (it in seq_len(max_iter)) {
    betas_prev <- betas

    # Calculate residuals using original (unweighted) X
    resid <- Y - X_orig_for_resid %*% betas

    scale_info <- compute_scale(resid, scale_scope, sigma_fixed)
    u <- scale_info$u
    
    # Calculate weights based on psi function
    w <- switch(psi_type,
                huber = pmin(1, k_huber / abs(u)),
                bisquare = ifelse(abs(u) <= c_tukey,
                                  (1 - (u / c_tukey)^2)^2,
                                  0))
    
    # Apply sqrt of weights for weighted least squares
    sqrtw <- sqrt(w)
    Xw <- X_orig_for_resid * sqrtw
    Yw <- Y * sqrtw
    
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
    fit_w <- solve_glm_core(glm_ctx_weighted, return_fitted = FALSE)

    # Update betas and XtWXi
    betas <- fit_w$betas
    XtWXi_final <- proj_w$XtXinv

    # Convergence check on coefficient updates.
    if (it >= min_iter) {
      scale <- max(1, max(abs(betas_prev)))
      delta_beta <- max(abs(betas - betas_prev)) / scale
      if (is.finite(delta_beta) && delta_beta <= tol) {
        break
      }
    }
  }
  
  # Final scale estimation/reporting under configured scope.
  resid_final <- Y - X_orig_for_resid %*% betas
  scale_final <- compute_scale(resid_final, scale_scope, sigma_fixed)
  sigma_robust <- scale_final$sigma_scalar
  
  # Final weights for output
  u_final <- scale_final$u
  w_final <- switch(psi_type,
                    huber = pmin(1, k_huber / abs(u_final)),
                    bisquare = ifelse(abs(u_final) <= c_tukey,
                                      (1 - (u_final / c_tukey)^2)^2,
                                      0))
  
  # Return results
  list(
    betas_robust = betas,
    XtWXi_final = XtWXi_final,
    sigma_robust_scale_final = sigma_robust,
    sigma_robust_scale_components = scale_final$sigma_components,
    robust_weights_final = w_final,
    dfres = dfres,
    scale_scope = scale_scope
  )
}
