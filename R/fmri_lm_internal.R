#' @title Internal Utilities for fMRI Linear Models
#' @description Low-level utilities used throughout the fmri_lm implementation
#' @keywords internal

#' Check if object is a formula
#' @keywords internal
#' @noRd
is.formula <- function(x) {
  inherits(x, "formula")
}

#' Fast Pre-projection of Design Matrix
#'
#' @description
#' Computes the projection matrix components needed for fast least squares.
#' This includes the QR decomposition and (X'X)^-1.
#'
#' @param X Design matrix
#' @return List containing:
#'   - qr: QR decomposition of X
#'   - XtXinv: (X'X)^-1 matrix
#'   - dfres: Residual degrees of freedom
#' @keywords internal
#' @noRd
.fast_preproject <- function(X) {
  # Ensure X is a matrix
  if (!is.matrix(X)) {
    X <- as.matrix(X)
  }
  XtX   <- crossprod(X)        # p × p
  # Add small ridge for stability if needed, but try direct first
  # Rchol <- tryCatch(chol(XtX), error = function(e) chol(XtX + diag(ncol(XtX)) * 1e-10))
  Rchol <- chol(XtX)           # p × p  upper‑triangular
  Pinv  <- backsolve(Rchol, t(X), transpose = TRUE)  # (Rchol^-1)' Xᵀ = (Rchol'^-1) Xᵀ = (XtX)^-1 Xᵀ -> p x n
  
  # (XᵀX)⁻¹ via Cholesky
  XtXinv <- chol2inv(Rchol)
  
  # QR decomposition
  qr_decomp <- qr(X)
  
  # Return everything needed for fast operations
  list(
    qr = qr_decomp,
    Pinv = Pinv,
    XtXinv = XtXinv,
    dfres = nrow(X) - qr_decomp$rank
  )
}

#' Fast Matrix-based Linear Model
#'
#' @description
#' Low-level function for fast least squares computation.
#' This is being phased out in favor of solve_glm_core.
#'
#' @param X Design matrix
#' @param Y Response matrix (n x v)
#' @param proj Pre-computed projection from .fast_preproject
#' @param return_fitted Whether to return fitted values
#' @return List with regression results
#' @keywords internal
#' @noRd
.fast_lm_matrix <- function(X, Y, proj, return_fitted = FALSE) {
  # Use pre-computed components
  B <- proj$Pinv %*% Y  # p × V
  
  if (return_fitted) {
    fitted <- X %*% B     # n × V
    resid  <- Y - fitted  # n × V
  } else {
    resid <- Y - X %*% B  # n × V
  }
  
  rss    <- colSums(resid^2)  # V-vector
  sigma2 <- rss / proj$dfres  # V-vector
  
  ret <- list(
    betas = B,
    rss = rss,
    sigma2 = sigma2,
    dfres = proj$dfres
  )
  
  if (return_fitted) {
    ret$fitted <- fitted
  }
  
  ret
}

#' Meta-analysis of Beta Statistics Across Runs
#'
#' @description
#' Combines beta statistics from multiple runs using fixed-effects meta-analysis.
#'
#' @param bstats_list List of beta statistics tibbles from each run
#' @param event_indices Indices of event-related parameters
#' @return Combined beta statistics tibble
#' @keywords internal
#' @noRd
meta_betas <- function(bstats_list, event_indices) {
  # Extract data from each run
  estimates <- lapply(bstats_list, function(x) x$data[[1]]$estimate[[1]])
  ses <- lapply(bstats_list, function(x) x$data[[1]]$se[[1]])
  
  # Number of runs and parameters
  n_runs <- length(estimates)
  n_params <- ncol(estimates[[1]])
  n_voxels <- nrow(estimates[[1]])
  
  # Pre-allocate results
  meta_estimate <- matrix(0, n_voxels, n_params)
  meta_se <- matrix(0, n_voxels, n_params)
  
  # Fixed-effects meta-analysis for each parameter
  for (p in 1:n_params) {
    # Extract parameter p across runs
    est_p <- sapply(estimates, function(x) x[, p])
    se_p <- sapply(ses, function(x) x[, p])
    
    # Inverse variance weights
    w_p <- 1 / se_p^2
    
    # Weighted average
    meta_estimate[, p] <- rowSums(est_p * w_p) / rowSums(w_p)
    meta_se[, p] <- 1 / sqrt(rowSums(w_p))
  }
  
  # Compute meta t-statistics
  meta_tstat <- meta_estimate / meta_se
  
  # Degrees of freedom (sum across runs minus parameters)
  df_total <- sum(sapply(bstats_list, function(x) x$df.residual))
  
  # P-values
  meta_prob <- 2 * pt(-abs(meta_tstat), df_total)
  
  # Return in same format as single-run results
  tibble::tibble(
    type = "beta",
    name = "parameter_estimates",
    stat_type = "tstat",
    df.residual = df_total,
    conmat = list(NULL),
    colind = list(NULL),
    data = list(tibble::tibble(
      estimate = list(meta_estimate),
      se = list(meta_se),
      stat = list(meta_tstat),
      prob = list(meta_prob),
      sigma = list(sqrt(rowMeans(sapply(bstats_list, function(x) x$data[[1]]$sigma[[1]]^2))))
    ))
  )
}

#' Meta-analysis of Contrasts Across Runs
#'
#' @description
#' Combines contrast results from multiple runs using fixed-effects meta-analysis.
#'
#' @param conres_list List of contrast result lists from each run
#' @return Combined contrast results
#' @keywords internal
#' @noRd
meta_contrasts <- function(conres_list) {
  # Get all unique contrast names
  all_names <- unique(unlist(lapply(conres_list, names)))
  
  # Process each contrast
  meta_results <- lapply(all_names, function(con_name) {
    # Extract this contrast from each run
    con_runs <- lapply(conres_list, function(x) x[[con_name]])
    con_runs <- Filter(Negate(is.null), con_runs)
    
    if (length(con_runs) == 0) return(NULL)
    
    # Get contrast type
    con_type <- con_runs[[1]]$type[1]
    stat_type <- con_runs[[1]]$stat_type[1]
    
    # Extract estimates and SEs
    estimates <- lapply(con_runs, function(x) x$data[[1]]$estimate)
    ses <- lapply(con_runs, function(x) x$data[[1]]$se)
    
    if (con_type == "contrast") {
      # Simple contrast - use inverse variance weighting
      est_mat <- do.call(cbind, estimates)
      se_mat <- do.call(cbind, ses)
      
      w_mat <- 1 / se_mat^2
      meta_est <- rowSums(est_mat * w_mat) / rowSums(w_mat)
      meta_se <- 1 / sqrt(rowSums(w_mat))
      
      # T-statistics and p-values
      df_total <- sum(sapply(con_runs, function(x) x$df.residual[1]))
      meta_t <- meta_est / meta_se
      meta_p <- 2 * pt(-abs(meta_t), df_total)
      
      # Return tibble
      tibble::tibble(
        type = con_type,
        name = con_name,
        stat_type = stat_type,
        df.residual = df_total,
        conmat = con_runs[[1]]$conmat,
        colind = con_runs[[1]]$colind,
        data = list(tibble::tibble(
          estimate = meta_est,
          se = meta_se,
          stat = meta_t,
          prob = meta_p
        ))
      )
    } else {
      # F-contrast - combine F-statistics
      f_stats <- sapply(con_runs, function(x) x$data[[1]]$stat)
      df1 <- nrow(con_runs[[1]]$conmat[[1]])
      df2_total <- sum(sapply(con_runs, function(x) x$df.residual[1]))
      
      # Average F-statistics (approximation)
      meta_f <- rowMeans(f_stats)
      meta_p <- pf(meta_f, df1, df2_total, lower.tail = FALSE)
      
      tibble::tibble(
        type = con_type,
        name = con_name,
        stat_type = stat_type,
        df.residual = df2_total,
        conmat = con_runs[[1]]$conmat,
        colind = con_runs[[1]]$colind,
        data = list(tibble::tibble(
          estimate = rowMeans(do.call(cbind, estimates)),
          se = rowMeans(do.call(cbind, ses)),
          stat = meta_f,
          prob = meta_p
        ))
      )
    }
  })
  
  # Remove NULLs and return
  Filter(Negate(is.null), meta_results)
}