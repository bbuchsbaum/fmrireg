# Internal Utilities for fMRI Linear Models
# Low-level utilities used throughout the fmri_lm implementation


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
  
  # QR decomposition for rank detection and stable computation
  qr_decomp <- qr(X, LAPACK = TRUE)
  rank <- qr_decomp$rank
  p <- ncol(X)
  n <- nrow(X)
  
  # Check for rank deficiency
  if (rank < p) {
    warning(sprintf("Design matrix is rank deficient: rank = %d, ncol = %d", rank, p))
  }
  
  # Use QR for stable computation
  if (rank == p) {
    # Full rank: use Cholesky for efficiency
    XtX <- crossprod(X)
    Rchol <- tryCatch(chol(XtX), error = function(e) {
      # Fallback to SVD if Cholesky fails
      svd_result <- svd(X)
      d <- svd_result$d
      tol <- max(dim(X)) * .Machine$double.eps * max(d)
      pos <- d > tol
      V <- svd_result$v[, pos, drop = FALSE]
      D_inv <- diag(1/d[pos], nrow = sum(pos))
      XtXinv <- V %*% D_inv^2 %*% t(V)
      return(list(XtXinv = XtXinv, method = "svd"))
    })
    
    if (is.list(Rchol)) {
      # SVD was used
      XtXinv <- Rchol$XtXinv
    } else {
      # Cholesky succeeded
      XtXinv <- chol2inv(Rchol)
    }
    Pinv <- XtXinv %*% t(X)
  } else {
    # Rank deficient: use SVD-based pseudoinverse
    svd_result <- svd(X)
    d <- svd_result$d
    tol <- max(dim(X)) * .Machine$double.eps * max(d)
    pos <- d > tol
    
    # Moore-Penrose pseudoinverse
    U <- svd_result$u[, pos, drop = FALSE]
    V <- svd_result$v[, pos, drop = FALSE]
    D_inv <- diag(1/d[pos], nrow = sum(pos))
    Pinv <- V %*% D_inv %*% t(U)
    XtXinv <- V %*% D_inv^2 %*% t(V)
  }
  
  # Return everything needed for fast operations
  list(
    qr = qr_decomp,
    Pinv = Pinv,
    XtXinv = XtXinv,
    dfres = n - rank,
    rank = rank,
    is_full_rank = (rank == p)
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

#' Contrast statistics using stored QR decompositions
#'
#' Computes voxelwise t and F contrasts when each voxel stores only the
#' QR decomposition of its design matrix. This avoids keeping full
#' \code{XtXinv} matrices in memory.
#'
#' @param Betas p x V matrix of regression coefficients.
#' @param qr_list List of QR objects, one per voxel.
#' @param sigma Numeric vector of residual standard deviations per voxel.
#' @param conlist Named list of t contrast vectors (with \code{colind} attribute).
#' @param fconlist Named list of F contrast matrices (with \code{colind} attribute).
#' @param dfres Residual degrees of freedom.
#' @param robust_weights_list Optional list of robust weights per voxel.
#' @param ar_order AR model order used for effective df calculation.
#' @return List of tibbles, one per contrast.
#' @keywords internal
#' @noRd
fit_lm_contrasts_voxelwise_qr <- function(Betas, qr_list, sigma,
                                          conlist, fconlist, dfres,
                                          robust_weights_list = NULL,
                                          ar_order = 0) {
  p <- nrow(Betas)
  V <- ncol(Betas)

  results <- list()

  for (nm in names(conlist)) {
    l <- conlist[[nm]]
    colind <- attr(l, "colind")
    full_l <- numeric(p)
    full_l[colind] <- l

    est <- se <- stat <- prob <- sigma_out <- numeric(V)

    for (v in seq_len(V)) {
      qr_v <- qr_list[[v]]
      rw <- if (!is.null(robust_weights_list)) robust_weights_list[[v]] else NULL

      Qtl <- qr.qty(qr_v, full_l)
      var_con <- sum(Qtl[1:qr_v$rank]^2) * sigma[v]^2

      est[v] <- sum(full_l * Betas[, v])
      se[v]  <- sqrt(var_con)

      df_eff <- if (!is.null(rw) || ar_order > 0) {
        n <- dfres + p
        calculate_effective_df(n, p, rw, ar_order, method = "simple")
      } else {
        dfres
      }

      stat[v] <- ifelse(se[v] < .Machine$double.eps^0.5, 0, est[v] / se[v])
      prob[v] <- 2 * pt(-abs(stat[v]), df_eff)
      sigma_out[v] <- sigma[v]
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
      qr_v <- qr_list[[v]]
      rw <- if (!is.null(robust_weights_list)) robust_weights_list[[v]] else NULL

      RinvLt <- backsolve(qr.R(qr_v), t(full_L), upper.tri = TRUE)
      M <- crossprod(RinvLt)
      Cinv <- tryCatch(solve(M), error = function(e) {
        warning("Singular matrix in F contrast computation")
        matrix(NaN, nrow(M), ncol(M))
      })

      LB <- full_L %*% Betas[, v]

      qf <- if (any(is.nan(Cinv))) {
        NaN
      } else {
        drop(t(LB) %*% Cinv %*% LB)
      }

      est[v] <- qf / nrow(full_L)
      se[v]  <- sigma[v]^2

      df_eff <- if (!is.null(rw) || ar_order > 0) {
        n <- dfres + p
        calculate_effective_df(n, p, rw, ar_order, method = "simple")
      } else {
        dfres
      }

      Fval <- ifelse(abs(se[v]) < .Machine$double.eps^0.5 || is.nan(qf),
                      NaN, est[v] / se[v])
      stat[v] <- Fval
      prob[v] <- pf(Fval, nrow(full_L), df_eff, lower.tail = FALSE)
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

#' Initialize storage for voxelwise contrast results
#'
#' Creates pre-allocated numeric vectors for each contrast so that results can
#' be filled in incrementally while processing voxels in chunks.
#'
#' @param conlist Named list of t-contrast vectors.
#' @param fconlist Named list of F-contrast matrices.
#' @param n_voxels Number of voxels to allocate storage for.
#' @keywords internal
#' @noRd
initialize_contrast_storage <- function(conlist, fconlist, n_voxels) {
  storage <- list()

  for (nm in names(conlist)) {
    l <- conlist[[nm]]
    storage[[nm]] <- list(
      type = "contrast",
      stat_type = "tstat",
      conmat = l,
      colind = attr(l, "colind"),
      estimate = numeric(n_voxels),
      se = numeric(n_voxels),
      stat = numeric(n_voxels),
      prob = numeric(n_voxels),
      sigma = numeric(n_voxels)
    )
  }

  for (nm in names(fconlist)) {
    L <- fconlist[[nm]]
    storage[[nm]] <- list(
      type = "Fcontrast",
      stat_type = "Fstat",
      conmat = L,
      colind = attr(L, "colind"),
      estimate = numeric(n_voxels),
      se = numeric(n_voxels),
      stat = numeric(n_voxels),
      prob = numeric(n_voxels)
    )
  }

  storage
}

#' Store contrast results for a single voxel
#'
#' Helper used by `fit_lm_contrasts_voxelwise_chunked()` to compute and store
#' contrast statistics for one voxel using its QR decomposition.
#'
#' @param storage Contrast storage list from `initialize_contrast_storage()`.
#' @param voxel_index Integer voxel index being processed.
#' @param qr_or_XtXinv QR decomposition or precomputed XtX inverse for this voxel.
#' @param beta_w Regression coefficients for this voxel.
#' @param sigma_w Residual standard deviation for this voxel.
#' @param conlist Named list of t-contrasts.
#' @param fconlist Named list of F-contrasts.
#' @param dfres Residual degrees of freedom.
#' @param ar_order AR model order.
#' @keywords internal
#' @noRd
store_voxel_contrasts <- function(storage, voxel_index, qr_or_XtXinv, beta_w,
                                  sigma_w, conlist, fconlist, dfres,
                                  ar_order = 0, robust_weights = NULL) {
  XtXinv <- if (is.matrix(qr_or_XtXinv)) {
    qr_or_XtXinv
  } else {
    chol2inv(qr.R(qr_or_XtXinv))
  }
  Bv <- matrix(beta_w, ncol = 1)

  for (nm in names(conlist)) {
    l <- conlist[[nm]]
    colind <- attr(l, "colind")
    full_l <- matrix(0, nrow = 1, ncol = nrow(XtXinv))
    full_l[, colind] <- l

    res <- .fast_t_contrast(Bv, sigma_w^2, XtXinv, full_l, dfres,
                            robust_weights = robust_weights,
                            ar_order = ar_order)

    storage[[nm]]$estimate[voxel_index] <- res$estimate
    storage[[nm]]$se[voxel_index] <- res$se
    storage[[nm]]$stat[voxel_index] <- res$stat
    storage[[nm]]$prob[voxel_index] <- res$prob
    storage[[nm]]$sigma[voxel_index] <- res$sigma
  }

  for (nm in names(fconlist)) {
    L <- fconlist[[nm]]
    colind <- attr(L, "colind")
    full_L <- matrix(0, nrow = nrow(L), ncol = nrow(XtXinv))
    full_L[, colind] <- L

    res <- .fast_F_contrast(Bv, sigma_w^2, XtXinv, full_L, dfres,
                            robust_weights = robust_weights,
                            ar_order = ar_order)

    storage[[nm]]$estimate[voxel_index] <- res$estimate
    storage[[nm]]$se[voxel_index] <- res$se
    storage[[nm]]$stat[voxel_index] <- res$stat
    storage[[nm]]$prob[voxel_index] <- res$prob
  }

  storage
}

#' Format contrast storage into result tibbles
#'
#' Converts the contrast storage structure returned by
#' `initialize_contrast_storage()` into the list-of-tibbles format used by other
#' contrast functions.
#'
#' @param storage List returned by `initialize_contrast_storage()`.
#' @param dfres Residual degrees of freedom.
#' @keywords internal
#' @noRd
format_contrast_results <- function(storage, dfres) {
  results <- vector("list", length(storage))
  names(results) <- names(storage)

  for (nm in names(storage)) {
    s <- storage[[nm]]
    dat <- tibble::tibble(
      estimate = s$estimate,
      se = s$se,
      stat = s$stat,
      prob = s$prob,
      sigma = s$sigma %||% NULL
    )

    if (all(is.na(s$sigma))) {
      dat$sigma <- NULL
    }

    results[[nm]] <- tibble::tibble(
      type = s$type,
      name = nm,
      stat_type = s$stat_type,
      df.residual = dfres,
      conmat = list(s$conmat),
      colind = list(s$colind),
      data = list(dat)
    )
  }

  results
}

#' Chunked voxelwise contrast computation
#'
#' Implements a memory-efficient voxelwise contrast engine by processing voxels
#' in small chunks. Each voxel is whitened and analysed independently, and only
#' temporary matrices for the current chunk are kept in memory.
#'
#' @param X_run Design matrix for the current run.
#' @param Y_run Data matrix (time points \eqn{\times} voxels) for the run.
#' @param phi_matrix AR coefficients matrix (order \eqn{\times} voxels).
#' @param conlist Named list of t-contrast vectors.
#' @param fconlist Named list of F-contrast matrices.
#' @param chunk_size Number of voxels to process per chunk.
#' @keywords internal
#' @noRd
fit_lm_contrasts_voxelwise_chunked <- function(X_run, Y_run, phi_matrix,
                                               conlist, fconlist,
                                               robust_options = NULL,
                                               chunk_size = 100) {
  if (!is.matrix(X_run)) X_run <- as.matrix(X_run)
  if (!is.matrix(Y_run)) Y_run <- as.matrix(Y_run)

  n_voxels <- ncol(Y_run)
  if (is.null(dim(phi_matrix))) {
    phi_matrix <- matrix(phi_matrix, ncol = n_voxels)
  }
  if (ncol(phi_matrix) != n_voxels) {
    stop("phi_matrix columns must equal number of voxels")
  }

  ar_order <- nrow(phi_matrix)
  dfres <- nrow(X_run) - qr(X_run)$rank
  n_chunks <- ceiling(n_voxels / chunk_size)

  storage <- initialize_contrast_storage(conlist, fconlist, n_voxels)

  for (chunk_idx in seq_len(n_chunks)) {
    idx <- ((chunk_idx - 1) * chunk_size + 1):min(chunk_idx * chunk_size,
                                                n_voxels)

    for (v in idx) {
      phi_v <- phi_matrix[, v]
      tmp <- ar_whiten_transform(X_run, Y_run[, v, drop = FALSE],
                                 phi_v, exact_first = TRUE)
      X_w <- tmp$X
      Y_w <- tmp$Y

      if (!is.null(robust_options) && robust_options$type != FALSE) {
        proj_w <- .fast_preproject(X_w)
        ctx_w <- glm_context(X = X_w, Y = Y_w, proj = proj_w, phi_hat = phi_v)
        rfit <- robust_iterative_fitter(ctx_w, robust_options, X_w)
        beta_w <- rfit$betas_robust
        XtXinv <- rfit$XtWXi_final
        sigma_w <- rfit$sigma_robust_scale_final
        rw <- rfit$robust_weights_final
      } else {
        qr_w <- qr(X_w)
        beta_w <- qr.coef(qr_w, Y_w)
        XtXinv <- chol2inv(qr.R(qr_w))
        sigma_w <- sqrt(sum(qr.resid(qr_w, Y_w)^2) /
                        max(1, nrow(X_w) - qr_w$rank))
        rw <- NULL
      }

      storage <- store_voxel_contrasts(storage, v, XtXinv, beta_w, sigma_w,
                                       conlist, fconlist, dfres, ar_order, rw)
    }
  }

  format_contrast_results(storage, dfres)
}
