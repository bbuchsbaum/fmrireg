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

#' Beta statistics when each voxel has its own projection matrix
#'
#' Helper for voxelwise AR fitting where every voxel yields a distinct
#' whitened design matrix.  Computes standard errors and t statistics using
#' the per-voxel \code{XtXinv} matrices.
#'
#' @keywords internal
#' @noRd
beta_stats_matrix_voxelwise <- function(Betas, XtXinv_list, sigma, dfres,
                                        varnames,
                                        robust_weights_list = NULL,
                                        ar_order = 0) {
  V <- ncol(Betas)
  p <- nrow(Betas)

  est_mat  <- matrix(NA_real_, V, p)
  se_mat   <- matrix(NA_real_, V, p)
  t_mat    <- matrix(NA_real_, V, p)
  prob_mat <- matrix(NA_real_, V, p)

  for (v in seq_len(V)) {
    XtXinv <- XtXinv_list[[v]]
    se_scal <- sqrt(diag(XtXinv))

    rw <- if (!is.null(robust_weights_list)) robust_weights_list[[v]] else NULL

    df_eff <- if (!is.null(rw) || ar_order > 0) {
      n <- dfres + p
      calculate_effective_df(n, p, rw, ar_order, method = "simple")
    } else {
      dfres
    }

    se_vec <- se_scal * sigma[v]
    est_vec <- Betas[, v]
    t_vec <- ifelse(abs(se_vec) < .Machine$double.eps^0.5, 0, est_vec / se_vec)
    p_vec <- 2 * pt(-abs(t_vec), df_eff)

    est_mat[v, ]  <- est_vec
    se_mat[v, ]   <- se_vec
    t_mat[v, ]    <- t_vec
    prob_mat[v, ] <- p_vec
  }

  colnames(est_mat)  <- varnames
  colnames(se_mat)   <- varnames
  colnames(t_mat)    <- varnames
  colnames(prob_mat) <- varnames

  tibble::tibble(
    type = "beta",
    name = "parameter_estimates",
    stat_type = "tstat",
    df.residual = dfres,
    conmat = list(NULL),
    colind = list(NULL),
    data = list(tibble::tibble(
      estimate = list(est_mat),
      se = list(se_mat),
      stat = list(t_mat),
      prob = list(prob_mat),
      sigma = list(sigma)
    ))
  )
}

#' Contrast statistics with voxelwise projection matrices
#'
#' Computes t and F contrasts when each voxel has a distinct \code{XtXinv}.
#'
#' @keywords internal
#' @noRd
fit_lm_contrasts_voxelwise <- function(Betas, sigma2, XtXinv_list,
                                       conlist, fconlist, dfres,
                                       robust_weights_list = NULL,
                                       ar_order = 0) {
  p <- nrow(Betas)
  V <- ncol(Betas)

  results <- list()

  for (nm in names(conlist)) {
    l <- conlist[[nm]]
    colind <- attr(l, "colind")
    full_l <- matrix(0, nrow = 1, ncol = p)
    full_l[, colind] <- l

    est <- se <- stat <- prob <- sigma_out <- numeric(V)

    for (v in seq_len(V)) {
      rw <- if (!is.null(robust_weights_list)) robust_weights_list[[v]] else NULL
      res <- .fast_t_contrast(Betas[, v, drop = FALSE], sigma2[v],
                              XtXinv_list[[v]], full_l, dfres, rw, ar_order)
      est[v] <- res$estimate
      se[v]  <- res$se
      stat[v] <- res$stat
      prob[v] <- res$prob
      sigma_out[v] <- res$sigma
    }

    results[[nm]] <- tibble::tibble(
      type = "contrast",
      name = nm,
      stat_type = "tstat",
      df.residual = dfres,
      conmat = list(l),
      colind = list(colind),
      data = list(tibble::tibble(
        estimate = est,
        se = se,
        stat = stat,
        prob = prob,
        sigma = sigma_out
      ))
    )
  }

  for (nm in names(fconlist)) {
    L <- fconlist[[nm]]
    colind <- attr(L, "colind")
    full_L <- matrix(0, nrow = nrow(L), ncol = p)
    full_L[, colind] <- L

    est <- se <- stat <- prob <- numeric(V)

    for (v in seq_len(V)) {
      rw <- if (!is.null(robust_weights_list)) robust_weights_list[[v]] else NULL
      res <- .fast_F_contrast(Betas[, v, drop = FALSE], sigma2[v],
                              XtXinv_list[[v]], full_L, dfres, rw, ar_order)
      est[v]  <- res$estimate
      se[v]   <- res$se
      stat[v] <- res$stat
      prob[v] <- res$prob
    }

    results[[nm]] <- tibble::tibble(
      type = "Fcontrast",
      name = nm,
      stat_type = "Fstat",
      df.residual = dfres,
      conmat = list(L),
      colind = list(colind),
      data = list(tibble::tibble(
        estimate = est,
        se = se,
        stat = stat,
        prob = prob
      ))
    )
  }

  results
}