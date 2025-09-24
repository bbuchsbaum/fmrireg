#' Helper Functions for fmri_ttest
#'
#' Support functions for paired differences, sign flipping, and wrapper functions
#' for OLS and meta-analysis with voxelwise covariates.

#' Compute Paired Within-Subject Differences
#'
#' Creates within-subject differences (A - B) for paired t-tests from two blocks
#' with identical subjects and features.
#'
#' @param blkA First group_data block
#' @param blkB Second group_data block
#' @param rho Optional within-subject correlation between A and B. Can be:
#'   \itemize{
#'     \item Scalar: Same correlation for all subjects and features
#'     \item Vector of length S: Per-subject correlations
#'     \item Vector of length P: Per-feature correlations
#'     \item Matrix (S x P): Subject-feature specific correlations
#'   }
#'   Default is 0 (independence).
#'
#' @return A new block with Y = Y_A - Y_B and propagated variance if available
#' @export
paired_diff_block <- function(blkA, blkB, rho = 0) {
  if (!is.matrix(blkA$Y) || !is.matrix(blkB$Y)) {
    stop("Both blocks must contain Y matrices", call. = FALSE)
  }
  
  if (!identical(dim(blkA$Y), dim(blkB$Y))) {
    stop("Blocks must have identical dimensions", call. = FALSE)
  }
  
  # Check subjects match if available
  if (!is.null(blkA$meta$subjects) && !is.null(blkB$meta$subjects)) {
    if (!identical(blkA$meta$subjects, blkB$meta$subjects)) {
      stop("paired_diff_block: subjects must match between blocks", call. = FALSE)
    }
  }
  
  # Compute difference
  Y <- blkA$Y - blkB$Y
  
  # Propagate variance if available
  V <- NULL
  if (!is.null(blkA$V) && !is.null(blkB$V)) {
    VA <- blkA$V
    VB <- blkB$V
    
    # Broadcast rho to S x P
    S <- nrow(Y)
    P <- ncol(Y)
    
    rho_mat <- if (length(rho) == 1L) {
      matrix(rho, nrow = S, ncol = P)
    } else if (is.vector(rho) && length(rho) == S) {
      matrix(rho, nrow = S, ncol = P)
    } else if (is.vector(rho) && length(rho) == P) {
      matrix(rho, nrow = S, ncol = P, byrow = TRUE)
    } else {
      as.matrix(rho)
    }
    
    # Var(A - B) = Var(A) + Var(B) - 2*rho*sqrt(Var(A)*Var(B))
    V <- VA + VB - 2 * rho_mat * sqrt(pmax(0, VA) * pmax(0, VB))
    V[!is.finite(V)] <- NA_real_
  }
  
  # Create output block
  result <- list(
    Y = Y,
    V = V,
    T = NULL,
    DF = NULL,
    index = blkA$index,
    meta = blkA$meta,
    covars = blkA$covars,
    feature = blkA$feature
  )
  
  # Update contrast name if available
  if (!is.null(result$meta$contrast)) {
    contrast_a <- result$meta$contrast
    contrast_b <- if (!is.null(blkB$meta$contrast)) blkB$meta$contrast else "B"
    result$meta$contrast <- paste0(contrast_a, "_minus_", contrast_b)
  }
  
  result
}

#' Flip Sign of Coefficients
#'
#' Reverses the sign of coefficient-like outputs in a fit object.
#' Useful for switching between A-B and B-A conventions.
#'
#' @param fit An fmri_ttest_fit or similar object
#' @param coef Character vector of coefficient names to flip (default: all)
#'
#' @return Modified fit object with flipped signs
#'
#' @export
flip_sign <- function(fit, coef = NULL) {
  elements <- c("beta", "t", "z", "z_contrast", "p_contrast")
  
  for (nm in elements) {
    if (!is.null(fit[[nm]])) {
      if (is.null(coef)) {
        # Flip all
        fit[[nm]] <- -fit[[nm]]
      } else {
        # Flip specific coefficients
        if (is.matrix(fit[[nm]]) && !is.null(rownames(fit[[nm]]))) {
          idx <- which(rownames(fit[[nm]]) %in% coef)
          if (length(idx) > 0) {
            fit[[nm]][idx, ] <- -fit[[nm]][idx, ]
          }
        } else if (is.vector(fit[[nm]]) && nm %in% c("z_contrast", "p_contrast")) {
          fit[[nm]] <- -fit[[nm]]
        }
      }
    }
  }
  
  fit
}

#' Extended Meta-Analysis Fit with Voxelwise Covariate
#'
#' Wrapper for meta-analysis that supports an optional voxelwise covariate.
#' This extends the basic fmri_meta_fit to handle per-voxel covariates.
#'
#' @param Y Matrix of effect sizes (S x P)
#' @param V Matrix of variances (S x P)
#' @param X Design matrix (S x K)
#' @param method Meta-analysis method
#' @param robust Robust estimation method
#' @param huber_c Huber tuning constant
#' @param robust_iter Number of IRLS iterations
#' @param voxelwise Optional voxelwise covariate matrix (S x P)
#' @param center_voxelwise Logical; center voxelwise covariate per feature
#' @param voxel_name Name for voxelwise coefficient
#' @param n_threads Number of threads
#'
#' @return List with meta-analysis results
#'
#' @export
fmri_meta_fit_extended <- function(Y, V, X,
                                  method = c("pm", "dl", "fe", "reml"),
                                  robust = c("none", "huber"),
                                  huber_c = 1.345,
                                  robust_iter = 2,
                                  voxelwise = NULL,
                                  center_voxelwise = TRUE,
                                  voxel_name = "voxel_cov",
                                  n_threads = getOption("fmrireg.num_threads", 0)) {
  
  method <- match.arg(method)
  robust <- match.arg(robust)
  
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.matrix(V)) V <- as.matrix(V)
  if (!is.matrix(X)) X <- as.matrix(X)
  
  if (nrow(Y) != nrow(V) || ncol(Y) != ncol(V)) {
    stop("Y and V must have the same dimensions", call. = FALSE)
  }
  
  if (nrow(X) != nrow(Y)) {
    stop("X must have same number of rows as Y", call. = FALSE)
  }
  
  if (is.null(voxelwise)) {
    # Standard meta-analysis
    out <- fmri_meta_fit(
      Y = Y, V = V, X = X,
      method = method, robust = robust,
      huber_c = huber_c, robust_iter = robust_iter,
      n_threads = n_threads
    )
  } else {
    # Meta-analysis with voxelwise covariate
    C <- as.matrix(voxelwise)
    if (nrow(C) != nrow(Y) || ncol(C) != ncol(Y)) {
      stop("Voxelwise covariate must match Y dimensions", call. = FALSE)
    }
    
    if (isTRUE(center_voxelwise)) {
      C <- sweep(C, 2L, colMeans(C, na.rm = TRUE), FUN = "-")
    }
    
    out <- meta_fit_vcov_cpp(
      Y = Y, V = V, X = X, C = C,
      method = method, robust = robust,
      huber_c = huber_c, robust_iter = robust_iter,
      n_threads = as.integer(n_threads)
    )
    
    # Add row names
    rn <- c(colnames(X), voxel_name)
    rownames(out$beta) <- rownames(out$se) <- rownames(out$z) <- rn
  }
  
  # Add method and robust info
  out$method <- method
  out$robust <- robust
  
  # Convert OK to logical
  out$ok <- as.logical(out$ok)
  
  out
}

#' OLS Fit with Optional Voxelwise Covariate
#'
#' Wrapper for OLS t-tests that supports an optional voxelwise covariate.
#'
#' @param Y Outcome matrix (S x P)
#' @param X Design matrix (S x K)
#' @param voxelwise Optional voxelwise covariate matrix (S x P)
#' @param center_voxelwise Logical; center voxelwise covariate per feature
#' @param voxel_name Name for voxelwise coefficient
#'
#' @return List with OLS results
#'
#' @export
fmri_ols_fit <- function(Y, X,
                        voxelwise = NULL,
                        center_voxelwise = TRUE,
                        voxel_name = "voxel_cov") {
  
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  if (!is.matrix(X)) X <- as.matrix(X)
  
  if (nrow(Y) != nrow(X)) {
    stop("Y and X must have same number of rows", call. = FALSE)
  }
  
  if (is.null(voxelwise)) {
    # Standard OLS
    ols <- ols_t_cpp(Y, X)
    rn <- colnames(X)
    if (!is.null(rn)) {
      rownames(ols$beta) <- rownames(ols$se) <- rownames(ols$t) <- rn
    }
  } else {
    # OLS with voxelwise covariate
    C <- as.matrix(voxelwise)
    if (nrow(C) != nrow(Y) || ncol(C) != ncol(Y)) {
      stop("Voxelwise covariate must match Y dimensions", call. = FALSE)
    }
    
    if (isTRUE(center_voxelwise)) {
      C <- sweep(C, 2L, colMeans(C, na.rm = TRUE), FUN = "-")
    }
    
    ols <- ols_t_vcov_cpp(Y, X, C)
    rn <- c(colnames(X), voxel_name)
    if (!is.null(rn)) {
      rownames(ols$beta) <- rownames(ols$se) <- rownames(ols$t) <- rn
    }
  }
  
  ols
}

#' Extract Coefficient Image from fmri_ttest_fit
#'
#' Creates a NeuroVol image from coefficients in an fmri_ttest_fit object.
#'
#' @param object An fmri_ttest_fit object
#' @param coef Character or integer; coefficient to extract
#' @param statistic Character string; type of statistic to extract:
#'   \itemize{
#'     \item "estimate": Coefficient estimates (beta values)
#'     \item "se": Standard errors (if available)
#'     \item "z": Z-scores
#'     \item "p": P-values
#'   }
#' @param ... Additional arguments (e.g., mask to apply)
#'
#' @return NeuroVol object or numeric vector
#'
#' @export
coef_image.fmri_ttest_fit <- function(object, coef = 1, 
                                      statistic = c("estimate", "se", "z", "p"),
                                      ...) {
  statistic <- match.arg(statistic)
  
  # Extract additional arguments
  dots <- list(...)
  mask <- dots$mask
  
  # Map statistic names to internal names used in fmri_ttest_fit
  type_map <- c(
    "estimate" = "beta",
    "se" = "se",
    "z" = "z", 
    "p" = "p"
  )
  
  type <- type_map[statistic]
  
  # Get the requested data
  data_mat <- object[[type]]
  if (is.null(data_mat)) {
    # If z-scores not precomputed, compute from t-statistics
    if (statistic == "z" && !is.null(object$t)) {
      data_mat <- object$t
      # Convert t to z (will be done below)
    } else if (statistic == "se" && is.null(object$se)) {
      stop("Standard errors not available in object", call. = FALSE)
    } else {
      stop(paste0("Statistic '", statistic, "' not available in object"), call. = FALSE)
    }
  }
  
  # Extract coefficient
  if (is.character(coef)) {
    if (!coef %in% rownames(data_mat)) {
      stop(paste0("Coefficient '", coef, "' not found"), call. = FALSE)
    }
    values <- data_mat[coef, ]
  } else {
    if (coef > nrow(data_mat)) {
      stop("Coefficient index out of bounds", call. = FALSE)
    }
    values <- data_mat[coef, ]
  }
  
  # Convert t to z if needed
  if (statistic == "z" && is.null(object$z) && !is.null(object$t)) {
    # Convert t-statistics to z-scores
    df <- object$df
    if (!is.null(df)) {
      values <- qnorm(pt(abs(values), df = df, lower.tail = FALSE)) * sign(values)
    } else {
      warning("No df available, treating t as z", call. = FALSE)
    }
  }
  
  # Apply mask if provided
  if (!is.null(mask)) {
    if (inherits(mask, "NeuroVol")) {
      vol_data <- array(NA_real_, dim = dim(mask))
      mask_indices <- which(mask > 0)
      vol_data[mask_indices] <- values
      return(neuroim2::NeuroVol(vol_data, space(mask)))
    }
  }
  
  values
}
