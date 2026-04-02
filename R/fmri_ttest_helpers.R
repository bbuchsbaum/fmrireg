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

.fmri_ttest_materialize_effects <- function(gd, context = "analysis") {
  if (inherits(gd, "group_data_nifti")) {
    dat_full <- read_nifti_full(gd)
    if (is.null(dat_full$beta)) {
      stop(
        sprintf("group_data_nifti missing beta paths for %s", context),
        call. = FALSE
      )
    }
    return(dat_full$beta)
  }

  if (inherits(gd, "group_data_h5")) {
    arr <- read_h5_full(gd, stat = "beta")
    if (is.null(arr)) {
      stop(
        sprintf("group_data_h5 missing beta statistic for %s", context),
        call. = FALSE
      )
    }
    return(t(arr[, , which(dimnames(arr)$stat == "beta")]))
  }

  stop("This data source cannot be materialized into betas for fmri_ttest()", call. = FALSE)
}

.fmri_ttest_validate_weights_custom <- function(weights_custom, n_subjects, n_features) {
  if (is.null(weights_custom)) {
    stop("weights='custom' requires 'weights_custom'", call. = FALSE)
  }

  if (is.vector(weights_custom)) {
    if (length(weights_custom) != n_subjects) {
      stop("weights_custom must be length S or SxP to match data", call. = FALSE)
    }
    weights_custom <- as.numeric(weights_custom)
    if (any(!is.finite(weights_custom)) || any(weights_custom <= 0)) {
      stop("weights_custom must contain only finite positive values", call. = FALSE)
    }
    return(weights_custom)
  }

  if (is.matrix(weights_custom)) {
    if (!all(dim(weights_custom) == c(n_subjects, n_features))) {
      stop("weights_custom must be length S or SxP to match data", call. = FALSE)
    }
    if (any(!is.finite(weights_custom)) || any(weights_custom <= 0)) {
      stop("weights_custom must contain only finite positive values", call. = FALSE)
    }
    return(weights_custom)
  }

  stop("weights_custom must be length S or SxP to match data", call. = FALSE)
}

.fmri_ttest_group_term <- function(X, covars) {
  has_group <- "group" %in% colnames(covars) &&
    is.factor(covars$group) &&
    nlevels(covars$group) == 2

  if (!has_group) {
    return(NULL)
  }

  candidates <- grep("^group[^:]*$", colnames(X), value = TRUE)
  if (length(candidates) != 1L) {
    return(NULL)
  }

  list(
    raw_name = candidates[[1]],
    canonical_name = "group",
    levels = levels(covars$group),
    has_group = TRUE
  )
}

.fmri_ttest_rename_rows <- function(x, old, new) {
  if (is.null(x) || is.null(old) || is.null(rownames(x))) {
    return(x)
  }

  idx <- which(rownames(x) == old)
  if (length(idx) == 1L) {
    rownames(x)[idx] <- new
  }
  x
}

.fmri_ttest_normalize_group_rows <- function(res, group_info) {
  if (is.null(group_info)) {
    return(res)
  }

  for (nm in c("beta", "se", "t", "z", "p", "q", "df")) {
    res[[nm]] <- .fmri_ttest_rename_rows(
      res[[nm]],
      old = group_info$raw_name,
      new = group_info$canonical_name
    )
  }

  res
}

.fmri_ttest_canonical_coef_names <- function(coef_names, group_info) {
  if (is.null(group_info)) {
    return(coef_names)
  }

  coef_names[coef_names == group_info$raw_name] <- group_info$canonical_name
  coef_names
}

.fmri_ttest_apply_group_sign <- function(res, group_info, target_sign, source_sign) {
  if (is.null(group_info) || identical(target_sign, source_sign)) {
    return(res)
  }

  idx <- which(rownames(res$beta) == group_info$canonical_name)
  if (length(idx) != 1L) {
    return(res)
  }

  for (nm in c("beta", "t", "z")) {
    if (!is.null(res[[nm]])) {
      res[[nm]][idx, ] <- -res[[nm]][idx, ]
    }
  }

  if (!is.null(res$z_contrast)) {
    res$z_contrast <- -res$z_contrast
  }

  if (!is.null(res$t_contrast)) {
    res$t_contrast <- -res$t_contrast
  }

  res
}

.fmri_ttest_resolve_contrast <- function(contrast, coef_names, group_info = NULL) {
  if (is.null(contrast)) {
    return(NULL)
  }

  weights <- stats::setNames(rep(0, length(coef_names)), coef_names)

  if (is.null(names(contrast))) {
    if (length(contrast) != length(coef_names)) {
      stop("Unnamed contrast must have length equal to number of coefficients", call. = FALSE)
    }
    weights[] <- as.numeric(contrast)
    return(weights)
  }

  name_map <- stats::setNames(coef_names, coef_names)
  if (!is.null(group_info)) {
    name_map[[group_info$raw_name]] <- group_info$canonical_name
    name_map[[group_info$canonical_name]] <- group_info$canonical_name
  }

  unknown <- setdiff(names(contrast), names(name_map))
  if (length(unknown) > 0) {
    stop(
      sprintf(
        "Unknown contrast coefficient(s): %s",
        paste(unknown, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  for (nm in names(contrast)) {
    weights[[name_map[[nm]]]] <- weights[[name_map[[nm]]]] + as.numeric(contrast[[nm]])
  }

  weights
}

.fmri_ttest_raw_contrast_weights <- function(canonical_weights, coef_names, group_info, target_sign, source_sign) {
  raw_weights <- canonical_weights
  names(raw_weights) <- coef_names

  if (!is.null(group_info)) {
    raw_weights[[group_info$raw_name]] <- canonical_weights[[group_info$canonical_name]]
    raw_weights <- raw_weights[coef_names]
    if (!identical(target_sign, source_sign)) {
      raw_weights[[group_info$raw_name]] <- -raw_weights[[group_info$raw_name]]
    }
  }

  raw_weights
}

.fmri_ttest_single_coef_contrast <- function(res, weights) {
  nz <- which(abs(weights) > 0)
  if (length(nz) != 1L) {
    return(NULL)
  }

  idx <- nz[[1]]
  scale <- unname(weights[[idx]])
  coef_name <- names(weights)[[idx]]

  estimate <- scale * res$beta[coef_name, ]
  se <- if (!is.null(res$se)) abs(scale) * res$se[coef_name, ] else NULL
  t_stat <- if (!is.null(res$t)) sign(scale) * res$t[coef_name, ] else NULL
  z_stat <- if (!is.null(res$z)) sign(scale) * res$z[coef_name, ] else NULL
  p_val <- if (!is.null(res$p)) res$p[coef_name, ] else NULL
  df <- if (!is.null(res$df)) res$df[coef_name, ] else NULL

  list(
    estimate = estimate,
    se = se,
    t = t_stat,
    z = z_stat,
    p = p_val,
    df = df
  )
}

.fmri_ttest_exact_ols_contrast <- function(ols, X, raw_weights) {
  XtX_inv <- solve(crossprod(X))
  scale <- as.numeric(crossprod(raw_weights, XtX_inv %*% raw_weights))
  diag_inv <- diag(XtX_inv)
  ref_idx <- which(is.finite(diag_inv) & diag_inv > 0)[1]
  sigma2 <- (ols$se[ref_idx, ]^2) / diag_inv[[ref_idx]]
  estimate <- colSums(ols$beta * raw_weights)
  se <- sqrt(scale * sigma2)
  t_stat <- estimate / se
  df <- rep(ols$df, length(estimate))
  p_val <- 2 * stats::pt(abs(t_stat), df = df, lower.tail = FALSE)
  z_stat <- stats::qnorm(pmax(1e-300, 1 - p_val / 2)) * sign(t_stat)

  list(
    estimate = estimate,
    se = se,
    t = t_stat,
    z = z_stat,
    p = p_val,
    df = df
  )
}

.fmri_ttest_store_contrast <- function(res, contrast_stats) {
  if (is.null(contrast_stats)) {
    return(res)
  }

  res$beta_contrast <- contrast_stats$estimate
  res$se_contrast <- contrast_stats$se
  res$t_contrast <- contrast_stats$t
  res$z_contrast <- contrast_stats$z
  res$p_contrast <- contrast_stats$p
  res$df_contrast <- contrast_stats$df
  res
}

.fmri_ttest_feature_group <- function(gd, block = NULL) {
  if (!is.null(block) && !is.null(block$feature) && !is.null(block$feature$group)) {
    return(block$feature$group)
  }

  attr_group <- attr(gd, "fmrireg_feature_group", exact = TRUE)
  if (!is.null(attr_group)) {
    return(attr_group)
  }

  if (inherits(gd, c("gds", "gds_plan", "group_data_gds"))) {
    gd_resolved <- if (inherits(gd, "gds_plan")) {
      tryCatch(fmrigds::compute(gd), error = function(e) gd)
    } else {
      gd
    }
    rd <- tryCatch(fmrigds::row_data(gd_resolved), error = function(e) NULL)
    if (!is.null(rd)) {
      if ("feature_group" %in% names(rd)) {
        return(rd$feature_group)
      }
      if ("group" %in% names(rd)) {
        return(rd$group)
      }
    }
  }

  NULL
}

.fmri_ttest_sample_labels <- function(gd, block = NULL) {
  if (!is.null(block) && !is.null(block$feature) && !is.null(block$feature$label)) {
    return(block$feature$label)
  }

  attr_labels <- attr(gd, "fmrireg_sample_labels", exact = TRUE)
  if (!is.null(attr_labels)) {
    return(attr_labels)
  }

  if (inherits(gd, c("gds", "gds_plan", "group_data_gds"))) {
    gd_resolved <- if (inherits(gd, "gds_plan")) {
      tryCatch(fmrigds::compute(gd), error = function(e) gd)
    } else {
      gd
    }
    rd <- tryCatch(fmrigds::row_data(gd_resolved), error = function(e) NULL)
    if (!is.null(rd)) {
      if ("sample" %in% names(rd)) {
        return(as.character(rd$sample))
      }
      if (!is.null(rownames(rd))) {
        return(rownames(rd))
      }
    }
  }

  NULL
}

.fmri_ttest_apply_mc <- function(res, mc, alpha, feature_group) {
  if (is.null(mc)) {
    return(res)
  }

  if (identical(mc, "spatial_fdr") && is.null(feature_group)) {
    stop("mc='spatial_fdr' requires feature grouping metadata", call. = FALSE)
  }

  adjust_method <- switch(mc, bh = "BH", by = "BY", NULL)

  if (!is.null(res$z_contrast)) {
    if (!is.null(adjust_method)) {
      res$q_contrast <- stats::p.adjust(res$p_contrast, method = adjust_method)
    } else {
      out <- spatial_fdr(
        z = as.numeric(res$z_contrast),
        group = feature_group,
        alpha = alpha
      )
      res$q_contrast <- out$q
      res$reject_contrast <- out$reject
    }
    return(res)
  }

  corrs <- vector("list", nrow(res$z))
  rejects <- vector("list", nrow(res$z))
  for (i in seq_len(nrow(res$z))) {
    if (!is.null(adjust_method)) {
      corrs[[i]] <- stats::p.adjust(res$p[i, ], method = adjust_method)
    } else {
      out <- spatial_fdr(
        z = as.numeric(res$z[i, ]),
        group = feature_group,
        alpha = alpha
      )
      corrs[[i]] <- out$q
      rejects[[i]] <- out$reject
    }
  }

  res$q <- do.call(rbind, corrs)
  rownames(res$q) <- rownames(res$z)
  if (identical(mc, "spatial_fdr")) {
    res$reject <- do.call(rbind, rejects)
    rownames(res$reject) <- rownames(res$z)
  }
  res
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
  elements <- c("beta", "t", "z", "z_contrast", "t_contrast")
  
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
        } else if (is.vector(fit[[nm]]) && nm %in% c("z_contrast", "t_contrast")) {
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
      if (is.matrix(df)) {
        if (is.character(coef)) {
          df <- df[coef, ]
        } else {
          df <- df[coef, ]
        }
      }
      p_two_sided <- 2 * stats::pt(abs(values), df = df, lower.tail = FALSE)
      values <- stats::qnorm(pmax(1e-300, 1 - p_two_sided / 2)) * sign(values)
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
