# Core Meta-Analysis Functions for Group fMRI Analysis
# Implements voxelwise and ROI-based meta-analysis with various methods

#' Fit Group-Level Meta-Analysis
#'
#' Performs voxelwise or ROI-based meta-analysis on group fMRI data using
#' fixed-effects, random-effects, or robust methods. Supports meta-regression
#' with covariates for group comparisons and other moderator analyses.
#'
#' @param data A group_data object created by \code{\link{group_data}}
#' @param formula Formula specifying the meta-regression model. Default is ~ 1 (intercept only).
#'   Use ~ 1 + group for group comparisons, or include continuous covariates.
#' @param method Character string specifying the meta-analysis method:
#'   \itemize{
#'     \item "fe": Fixed-effects (inverse variance weighted)
#'     \item "pm": Paule-Mandel random-effects (default, good for whole-brain)
#'     \item "dl": DerSimonian-Laird random-effects
#'     \item "reml": Restricted maximum likelihood random-effects
#'   }
#' @param robust Character string specifying robust estimation:
#'   \itemize{
#'     \item "none": No robust estimation (default)
#'     \item "huber": Huber M-estimator with IRLS
#'     \item "t": Student-t mixture model (for heavy-tailed distributions)
#'   }
#' @param weights Character string specifying weighting scheme:
#'   \itemize{
#'     \item "ivw": Inverse variance weighting (default)
#'     \item "equal": Equal weights for all subjects
#'     \item "custom": User-provided weights (must supply weights argument)
#'   }
#' @param weights_custom Numeric vector or matrix of custom weights (required if
#'   `weights = "custom"`). If a vector, length must equal the number of subjects.
#'   If a matrix, it must be subjects x features.
#'
#' @param combine For t-statistic-only data, combination method ("stouffer",
#'   "fisher", or "lancaster"). Stouffer combines z-scores and supports
#'   equal, inverse-variance, or custom weighting (via `weights`). Fisher uses
#'   equal weights. Lancaster implements a weighted Fisher method by mapping
#'   weights to per-subject degrees of freedom.
#' @param contrasts Optional numeric vector or matrix specifying fit-time exact contrasts. If a
#'   vector is provided, its names must match the column names of the design matrix X. A matrix
#'   should have columns corresponding to predictors and rows corresponding to contrasts.
#' @param return_cov Optional. If set to "tri", returns the packed upper-triangular Var(beta)
#'   per feature under `$cov` to enable exact post-hoc contrasts via `contrast()`.
#' @param chunk_size Number of voxels to process at once (default: 10000)
#' @param n_threads Number of parallel threads to use. Defaults to fmrireg.num_threads option.
#' @param verbose Logical. Print progress messages (default: TRUE)
#'
#' @return An fmri_meta object containing:
#'   \itemize{
#'     \item coefficients: Meta-regression coefficients
#'     \item se: Standard errors
#'     \item tau2: Between-study variance (for random-effects)
#'     \item I2: I-squared heterogeneity statistic
#'     \item Q: Cochran's Q statistic
#'     \item model: Model specification
#'     \item data: Input group_data object
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Simple fixed-effects meta-analysis
#' fit <- fmri_meta(gd, method = "fe")
#'
#' # Random-effects with group comparison
#' fit <- fmri_meta(gd, formula = ~ 1 + group, method = "pm")
#'
#' # Robust meta-regression with continuous covariate
#' fit <- fmri_meta(gd, formula = ~ 1 + age + sex, method = "reml", robust = "huber")
#'
#' # Stouffer's Z for t-statistics only
#' fit <- fmri_meta(gd_tstat, combine = "stouffer")
#'
#' # Exact post-hoc contrasts by storing covariance
#' fit_cov <- fmri_meta(gd, formula = ~ 1 + group, method = "pm", return_cov = "tri")
#' con <- contrast(fit_cov, c("(Intercept)" = 0, group = 1))
#'
#' # Exact fit-time contrast without storing covariance
#' fit_con <- fmri_meta(gd, formula = ~ 1 + group, method = "pm",
#'                      contrasts = c("(Intercept)" = 0, group = 1))
#' }
fmri_meta <- function(data,
                      formula = ~ 1,
                      method = c("pm", "fe", "dl", "reml"),
                      robust = c("none", "huber", "t"),
                      weights = c("ivw", "equal", "custom"),
                      weights_custom = NULL,
                      combine = NULL,
                      contrasts = NULL,
                      return_cov = NULL,
                      chunk_size = 10000,
                      n_threads = getOption("fmrireg.num_threads", 0),
                      verbose = TRUE) {
  
  # Validate input
  validate_group_data(data)
  
  # Match arguments
  method <- match.arg(method)
  robust <- match.arg(robust)
  weights <- match.arg(weights)
  if (!is.null(return_cov)) return_cov <- match.arg(return_cov, c("tri"))
  
  # Check for custom weights
  if (weights == "custom" && is.null(weights_custom)) {
    stop("Must provide 'weights_custom' when weights = 'custom'", call. = FALSE)
  }
  
  # Set up parallel processing
  if (n_threads == 0) {
    n_threads <- parallel::detectCores() - 1
  }
  n_threads <- max(1, n_threads)
  
  # Dispatch based on data type
  UseMethod("fmri_meta")
}

#' @export
fmri_meta.group_data_h5 <- function(data,
                                    formula = ~ 1,
                                    method = c("pm", "fe", "dl", "reml"),
                                    robust = c("none", "huber", "t"),
                                    weights = c("ivw", "equal", "custom"),
                                    weights_custom = NULL,
                                    combine = NULL,
                                    contrasts = NULL,
                                    return_cov = NULL,
                                    chunk_size = 10000,
                                    n_threads = getOption("fmrireg.num_threads", 0),
                                    verbose = TRUE) {
  
  method <- match.arg(method)
  robust <- match.arg(robust)
  weights <- match.arg(weights)
  if (!is.null(return_cov)) return_cov <- match.arg(return_cov, c("tri"))
  
  # Parse formula and create design matrix
  model_info <- parse_meta_formula(formula, data)
  
  # Build contrast matrix if provided
  Cmat <- NULL; contrast_names <- NULL
  if (!is.null(contrasts)) {
    if (is.matrix(contrasts)) {
      # Expect rows = contrasts, cols = predictors
      if (ncol(contrasts) != ncol(model_info$X)) stop("contrasts columns must match predictors", call. = FALSE)
      Cmat <- t(contrasts)  # K x J
      contrast_names <- rownames(contrasts)
    } else if (is.numeric(contrasts)) {
      # Named numeric vector
      w <- rep(0, ncol(model_info$X)); names(w) <- colnames(model_info$X)
      if (!is.null(names(contrasts))) {
        w[names(contrasts)] <- contrasts
      } else if (length(contrasts) == length(w)) {
        w <- contrasts
      } else {
        stop("Unnamed contrast vector must have length equal to number of predictors", call. = FALSE)
      }
      Cmat <- matrix(w, ncol = 1)
      contrast_names <- "contrast1"
    } else {
      stop("Unsupported 'contrasts' type; provide matrix or numeric vector", call. = FALSE)
    }
  }

  # Initialize results storage
  n_voxels <- data$n_voxels
  n_coef <- ncol(model_info$X)
  
  results <- list(
    coefficients = matrix(NA_real_, n_voxels, n_coef),
    se = matrix(NA_real_, n_voxels, n_coef),
    tau2 = rep(NA_real_, n_voxels),
    I2 = rep(NA_real_, n_voxels),
    Q = rep(NA_real_, n_voxels),
    Q_df = rep(NA_real_, n_voxels)
  )
  if (!is.null(Cmat)) {
    n_con <- ncol(Cmat)
    results$con_est <- matrix(NA_real_, n_voxels, n_con)
    results$con_se  <- matrix(NA_real_, n_voxels, n_con)
    results$con_z   <- matrix(NA_real_, n_voxels, n_con)
  }
  
  # Set column names
  colnames(results$coefficients) <- colnames(model_info$X)
  colnames(results$se) <- colnames(model_info$X)
  
  # Process in chunks
  n_chunks <- ceiling(n_voxels / chunk_size)
  
  if (verbose) {
    cat("Processing", n_voxels, "voxels in", n_chunks, "chunks using", n_threads, "threads\n")
    pb <- txtProgressBar(min = 0, max = n_chunks, style = 3)
  }
  
  for (chunk_idx in seq_len(n_chunks)) {
    # Get voxel indices for this chunk
    start_idx <- (chunk_idx - 1) * chunk_size + 1
    end_idx <- min(chunk_idx * chunk_size, n_voxels)
    voxel_indices <- start_idx:end_idx
    
    # Read data for this chunk
    chunk_data <- read_h5_chunk(data, voxel_indices)
    
    # Fit meta-analysis for each voxel in chunk
    chunk_results <- fit_meta_chunk(
      chunk_data,
      model_info$X,
      method = if (!is.null(combine)) paste0("combine:", combine) else method,
      robust = robust,
      weights = weights,
      weights_custom = weights_custom,
      n_threads = n_threads,
      combine = combine,
      contrasts = Cmat,
      return_cov = return_cov
    )
    
    # Store results
    results$coefficients[voxel_indices, ] <- chunk_results$coefficients
    results$se[voxel_indices, ] <- chunk_results$se
    results$tau2[voxel_indices] <- chunk_results$tau2
    results$I2[voxel_indices] <- chunk_results$I2
    results$Q[voxel_indices] <- chunk_results$Q
    results$Q_df[voxel_indices] <- chunk_results$Q_df
    if (!is.null(Cmat) && !is.null(chunk_results$con_est)) {
      results$con_est[voxel_indices, ] <- t(chunk_results$con_est)
      results$con_se[voxel_indices, ]  <- t(chunk_results$con_se)
      results$con_z[voxel_indices, ]   <- t(chunk_results$con_z)
    }
    if (identical(return_cov, "tri") && !is.null(chunk_results$cov_tri)) {
      if (is.null(results$cov_tri)) {
        results$cov_tri <- matrix(NA_real_, nrow = nrow(chunk_results$cov_tri), ncol = n_voxels)
      }
      results$cov_tri[, voxel_indices] <- chunk_results$cov_tri
    }
    
    if (verbose) {
      setTxtProgressBar(pb, chunk_idx)
    }
  }
  
  if (verbose) {
    close(pb)
    cat("\nMeta-analysis complete\n")
  }
  
  # Create fmri_meta object
  structure(
    list(
      coefficients = results$coefficients,
      se = results$se,
      tau2 = results$tau2,
      I2 = results$I2,
      Q = results$Q,
      Q_df = results$Q_df,
      model = model_info,
      method = if (!is.null(combine)) paste0("combine:", combine) else method,
      robust = robust,
      weights = weights,
      data = data,
      formula = formula,
      n_voxels = n_voxels,
      n_subjects = n_subjects(data),
      contrasts = if (!is.null(Cmat)) list(
        names = if (!is.null(contrast_names)) contrast_names else paste0("c", seq_len(ncol(Cmat))),
        estimate = results$con_est,
        se = results$con_se,
        z = results$con_z
      ) else NULL,
      cov = if (identical(return_cov, "tri")) list(
        type = "tri",
        tri = results$cov_tri,
        coef_names = colnames(model_info$X)
      ) else NULL
    ),
    class = "fmri_meta"
  )
}

#' @export
fmri_meta.group_data_nifti <- function(data,
                                       formula = ~ 1,
                                       method = c("pm", "fe", "dl", "reml"),
                                       robust = c("none", "huber", "t"),
                                       weights = c("ivw", "equal", "custom"),
                                       weights_custom = NULL,
                                       combine = NULL,
                                       contrasts = NULL,
                                       return_cov = NULL,
                                       chunk_size = 10000,
                                       n_threads = getOption("fmrireg.num_threads", 0),
                                       verbose = TRUE) {
  
  # Similar to group_data_h5 but using read_nifti_chunk
  method <- match.arg(method)
  robust <- match.arg(robust)
  weights <- match.arg(weights)
  if (!is.null(return_cov)) return_cov <- match.arg(return_cov, c("tri"))
  
  # Parse formula and create design matrix
  model_info <- parse_meta_formula(formula, data)
  
  # Determine voxel indices to process
  if (!is.null(data$mask_data)) {
    voxel_indices <- which(as.vector(data$mask_data) > 0)
    n_voxels <- length(voxel_indices)
  } else {
    n_voxels <- data$n_voxels
    voxel_indices <- seq_len(n_voxels)
  }
  
  # Build contrast matrix if provided
  Cmat <- NULL; contrast_names <- NULL
  if (!is.null(contrasts)) {
    if (is.matrix(contrasts)) {
      if (ncol(contrasts) != ncol(model_info$X)) stop("contrasts columns must match predictors", call. = FALSE)
      Cmat <- t(contrasts)
      contrast_names <- rownames(contrasts)
    } else if (is.numeric(contrasts)) {
      w <- rep(0, ncol(model_info$X)); names(w) <- colnames(model_info$X)
      if (!is.null(names(contrasts))) {
        w[names(contrasts)] <- contrasts
      } else if (length(contrasts) == length(w)) {
        w <- contrasts
      } else {
        stop("Unnamed contrast vector must have length equal to number of predictors", call. = FALSE)
      }
      Cmat <- matrix(w, ncol = 1)
      contrast_names <- "contrast1"
    } else {
      stop("Unsupported 'contrasts' type; provide matrix or numeric vector", call. = FALSE)
    }
  }

  # Initialize results
  n_coef <- ncol(model_info$X)
  results <- list(
    coefficients = matrix(NA_real_, n_voxels, n_coef),
    se = matrix(NA_real_, n_voxels, n_coef),
    tau2 = rep(NA_real_, n_voxels),
    I2 = rep(NA_real_, n_voxels),
    Q = rep(NA_real_, n_voxels),
    Q_df = rep(NA_real_, n_voxels)
  )
  if (!is.null(Cmat)) {
    n_con <- ncol(Cmat)
    results$con_est <- matrix(NA_real_, n_voxels, n_con)
    results$con_se  <- matrix(NA_real_, n_voxels, n_con)
    results$con_z   <- matrix(NA_real_, n_voxels, n_con)
  }
  
  colnames(results$coefficients) <- colnames(model_info$X)
  colnames(results$se) <- colnames(model_info$X)
  
  # Process in chunks
  n_chunks <- ceiling(n_voxels / chunk_size)
  
  if (verbose) {
    cat("Processing", n_voxels, "voxels in", n_chunks, "chunks using", n_threads, "threads\n")
    pb <- txtProgressBar(min = 0, max = n_chunks, style = 3)
  }
  
  for (chunk_idx in seq_len(n_chunks)) {
    # Get voxel indices for this chunk
    start_idx <- (chunk_idx - 1) * chunk_size + 1
    end_idx <- min(chunk_idx * chunk_size, n_voxels)
    chunk_voxel_indices <- voxel_indices[start_idx:end_idx]
    
    # Read data for this chunk
    chunk_data <- read_nifti_chunk(data, chunk_voxel_indices)
    
    # Fit meta-analysis
    chunk_results <- fit_meta_chunk(
      chunk_data,
      model_info$X,
      method = if (!is.null(combine)) paste0("combine:", combine) else method,
      robust = robust,
      weights = weights,
      weights_custom = weights_custom,
      n_threads = n_threads,
      combine = combine,
      contrasts = Cmat,
      return_cov = return_cov
    )
    
    # Store results
    idx_range <- start_idx:end_idx
    results$coefficients[idx_range, ] <- chunk_results$coefficients
    results$se[idx_range, ] <- chunk_results$se
    results$tau2[idx_range] <- chunk_results$tau2
    results$I2[idx_range] <- chunk_results$I2
    results$Q[idx_range] <- chunk_results$Q
    results$Q_df[idx_range] <- chunk_results$Q_df
    if (identical(return_cov, "tri") && !is.null(chunk_results$cov_tri)) {
      if (is.null(results$cov_tri)) {
        results$cov_tri <- matrix(NA_real_, nrow = nrow(chunk_results$cov_tri), ncol = n_voxels)
      }
      results$cov_tri[, idx_range] <- chunk_results$cov_tri
    }
    if (!is.null(Cmat) && !is.null(chunk_results$con_est)) {
      results$con_est[idx_range, ] <- t(chunk_results$con_est)
      results$con_se[idx_range, ]  <- t(chunk_results$con_se)
      results$con_z[idx_range, ]   <- t(chunk_results$con_z)
    }
    
    if (verbose) {
      setTxtProgressBar(pb, chunk_idx)
    }
  }
  
  if (verbose) {
    close(pb)
    cat("\nMeta-analysis complete\n")
  }
  
  # Create fmri_meta object
  structure(
    list(
      coefficients = results$coefficients,
      se = results$se,
      tau2 = results$tau2,
      I2 = results$I2,
      Q = results$Q,
      Q_df = results$Q_df,
      model = model_info,
      method = method,
      robust = robust,
      weights = weights,
      data = data,
      formula = formula,
      n_voxels = n_voxels,
      n_subjects = n_subjects(data),
      voxel_indices = voxel_indices,  # Store which voxels were analyzed
      contrasts = if (!is.null(Cmat)) list(
        names = if (!is.null(contrast_names)) contrast_names else paste0("c", seq_len(ncol(Cmat))),
        estimate = results$con_est,
        se = results$con_se,
        z = results$con_z
      ) else NULL,
      cov = if (identical(return_cov, "tri")) list(
        type = "tri",
        tri = results$cov_tri,
        coef_names = colnames(model_info$X)
      ) else NULL
    ),
    class = "fmri_meta"
  )
}

#' @export
fmri_meta.group_data_csv <- function(data,
                                     formula = ~ 1,
                                     method = c("pm", "fe", "dl", "reml"),
                                     robust = c("none", "huber", "t"),
                                     weights = c("ivw", "equal", "custom"),
                                     weights_custom = NULL,
                                     combine = NULL,
                                     contrasts = NULL,
                                     return_cov = NULL,   # For consistency with other methods
                                     chunk_size = NULL,  # Not used for CSV
                                     n_threads = 1,       # Single-threaded for ROIs
                                     verbose = TRUE) {
  
  method <- match.arg(method)
  robust <- match.arg(robust)
  weights <- match.arg(weights)
  
  # Parse formula and create design matrix
  model_info <- parse_meta_formula(formula, data)
  
  # Get ROIs to analyze
  rois <- if (!is.null(data$rois)) data$rois else "overall"
  n_rois <- length(rois)
  
  # Initialize results
  n_coef <- ncol(model_info$X)
  results <- list(
    coefficients = matrix(NA_real_, n_rois, n_coef),
    se = matrix(NA_real_, n_rois, n_coef),
    tau2 = rep(NA_real_, n_rois),
    I2 = rep(NA_real_, n_rois),
    Q = rep(NA_real_, n_rois),
    Q_df = rep(NA_real_, n_rois),
    roi_names = rois
  )
  
  colnames(results$coefficients) <- colnames(model_info$X)
  colnames(results$se) <- colnames(model_info$X)
  rownames(results$coefficients) <- rois
  rownames(results$se) <- rois
  
  if (verbose) {
    cat("Analyzing", n_rois, "ROI(s)\n")
  }
  
  # Analyze each ROI
  for (i in seq_len(n_rois)) {
    roi <- rois[i]
    
    # Extract data for this ROI
    roi_data <- extract_csv_data(data, roi = if (roi != "overall") roi else NULL)
    
    # Fit meta-analysis
    if (!is.null(roi_data$beta) && !is.null(roi_data$se)) {
      roi_result <- fit_meta_single(
        y = roi_data$beta,
        se = roi_data$se,
        X = model_info$X,
        method = method,
        robust = robust,
        weights = weights,
        weights_custom = weights_custom
      )
      
      # Store results
      results$coefficients[i, ] <- roi_result$coefficients
      results$se[i, ] <- roi_result$se
      results$tau2[i] <- roi_result$tau2
      results$I2[i] <- roi_result$I2
      results$Q[i] <- roi_result$Q
      results$Q_df[i] <- roi_result$Q_df
    } else if (!is.null(combine) && !is.null(roi_data$t)) {
      # Combine t-statistics (t-only ROI path). Require intercept-only model.
      if (ncol(model_info$X) != 1) {
        stop("combine methods with t-statistics require intercept-only model (~ 1)", call. = FALSE)
      }
      tt <- as.numeric(roi_data$t)
      dfi <- roi_data$df
      if (is.null(dfi)) {
        stop("t-statistics provided without 'df'; cannot combine", call. = FALSE)
      }
      # Use same combine helper with optional weighting
      cm <- match.arg(tolower(combine), c("stouffer", "fisher", "lancaster"))
      # Expand to Sx1 matrices for helper
      tmat <- matrix(tt, ncol = 1)
      se_mat <- NULL  # No SE for ROI t-only
      # Build weights for ROI if requested
      zc_vec <- combine_t_statistics(tmat, df = dfi, method = cm,
                                     weights = weights, weights_custom = weights_custom,
                                     se_mat = se_mat)
      zc <- as.numeric(zc_vec)
      # Store as z with se=1
      results$coefficients[i, ] <- zc
      results$se[i, ] <- 1
      results$tau2[i] <- NA_real_
      results$I2[i] <- NA_real_
      results$Q[i] <- NA_real_
      results$Q_df[i] <- NA_real_
    }
    
    if (verbose && i %% 10 == 0) {
      cat(".")
    }
  }
  
  if (verbose) {
    cat("\nMeta-analysis complete\n")
  }
  
  # Create fmri_meta object
  structure(
    list(
      coefficients = results$coefficients,
      se = results$se,
      tau2 = results$tau2,
      I2 = results$I2,
      Q = results$Q,
      Q_df = results$Q_df,
      model = model_info,
      method = method,
      robust = robust,
      weights = weights,
      data = data,
      formula = formula,
      n_rois = n_rois,
      roi_names = rois,
      n_subjects = n_subjects(data)
    ),
    class = c("fmri_meta_roi", "fmri_meta")
  )
}

#' Parse Meta-Analysis Formula
#'
#' @param formula Formula object
#' @param data Group data object
#' @return List with design matrix X and other model info
#' @keywords internal
parse_meta_formula <- function(formula, data) {
  # Get covariates
  covariates <- get_covariates(data)
  
  if (is.null(covariates)) {
    # No covariates, just intercept
    n_subjects <- n_subjects(data)
    X <- matrix(1, nrow = n_subjects, ncol = 1)
    colnames(X) <- "(Intercept)"
  } else {
    # Create model matrix from formula and covariates
    X <- model.matrix(formula, data = covariates)
  }
  
  list(
    X = X,
    formula = formula,
    terms = terms(formula)
  )
}

#' Fit Meta-Analysis for a Chunk of Voxels
#'
#' @param chunk_data List with beta and se matrices
#' @param X Design matrix
#' @param method Meta-analysis method
#' @param robust Robust estimation method
#' @param weights Weighting scheme
#' @param weights_custom Custom weights
#' @param n_threads Number of threads
#' @return List with meta-analysis results
#' @keywords internal
fit_meta_chunk <- function(chunk_data, X, method, robust, weights, weights_custom, n_threads, combine = NULL, contrasts = NULL, return_cov = NULL) {
  if (!is.null(return_cov)) return_cov <- match.arg(return_cov, c("tri"))

  # 1) If we have beta/SE, use standard meta-regression
  if (!is.null(chunk_data$beta) && !is.null(chunk_data$se)) {
    Y <- chunk_data$beta  # subjects x voxels
    V <- chunk_data$se^2  # Convert SE to variance

    # Handle custom weights if provided
    if (weights == "equal") {
      V[] <- 1
    } else if (weights == "custom" && !is.null(weights_custom)) {
      # w = 1/v => v = 1/w
      V <- matrix(1/weights_custom, nrow = nrow(V), ncol = ncol(V))
    }

    # Set robust parameters
    huber_c <- if (robust == "huber") 1.345 else 1.345
    robust_iter <- if (robust == "huber") 2 else 0

    # Call C++ implementation (with or without contrasts)
    if (!is.null(contrasts) || identical(return_cov, "tri")) {
      # Compute and return covariance triangles; derive contrasts in R if requested
      result <- fmri_meta_fit_cov(
        Y = Y,
        V = V,
        X = X,
        method = method,
        robust = robust,
        huber_c = huber_c,
        robust_iter = robust_iter,
        n_threads = n_threads
      )
      out <- list(
        coefficients = t(result$beta),
        se = t(result$se),
        tau2 = result$tau2,
        I2 = result$I2_fe,
        Q = result$Q_fe,
        Q_df = result$df,
        cov_tri = result$cov_tri
      )
      if (!is.null(contrasts)) {
        # Compute contrasts: est = C' beta, se = sqrt(diag(C' V C)) using packed cov
        K <- ncol(X)
        tsize <- nrow(result$cov_tri)
        # reconstruct per voxel cov on the fly for c'Var c
        C <- contrasts
        J <- ncol(C)
        P <- ncol(Y)
        est <- matrix(NA_real_, J, P)
        se  <- matrix(NA_real_, J, P)
        for (p in seq_len(P)) {
          # unpack upper triangle to matrix
          Vtri <- result$cov_tri[, p]
          Vbeta <- matrix(0, ncol(X), ncol(X))
          idx <- 1L
          for (a in seq_len(ncol(X))) {
            for (b in a:ncol(X)) {
              Vbeta[a, b] <- Vtri[idx]
              Vbeta[b, a] <- Vtri[idx]
              idx <- idx + 1L
            }
          }
          bp <- out$coefficients[p, ]  # 1 x K
          for (j in seq_len(J)) {
            cw <- C[, j]
            est[j, p] <- sum(cw * bp)
            vv <- as.numeric(t(cw) %*% Vbeta %*% cw)
            se[j, p] <- if (is.finite(vv) && vv >= 0) sqrt(vv) else NA_real_
          }
        }
        out$con_est <- est
        out$con_se  <- se
        out$con_z   <- est / se
      }
      return(out)
    } else {
      result <- fmri_meta_fit(
        Y = Y,
        V = V,
        X = X,
        method = method,
        robust = robust,
        huber_c = huber_c,
        robust_iter = robust_iter,
        n_threads = n_threads
      )
      return(list(
        coefficients = t(result$beta),
        se = t(result$se),
        tau2 = result$tau2,
        I2 = result$I2_fe,
        Q = result$Q_fe,
        Q_df = result$df
      ))
    }
  }

  # 2) If t-statistics are available and no beta/SE, allow combine methods
  if (!is.null(chunk_data$t)) {
    # Only intercept-only model is supported for t-only combine
    if (ncol(X) != 1) {
      stop("combine methods with t-statistics require intercept-only model (~ 1)", call. = FALSE)
    }
    tmat <- chunk_data$t           # subjects x voxels
    df <- if (!is.null(chunk_data$df)) chunk_data$df else attr(chunk_data$t, "df")
    if (is.null(df)) {
      stop("t-statistics provided without 'df'; cannot combine", call. = FALSE)
    }
    if (is.null(combine)) {
      stop("Only t-statistics available; please set combine = 'stouffer' or 'fisher'", call. = FALSE)
    }
    method_combine <- match.arg(tolower(combine), c("stouffer", "fisher", "lancaster"))
    se_mat <- if (!is.null(chunk_data$se)) chunk_data$se else NULL
    zc <- combine_t_statistics(tmat, df, method = method_combine,
                               weights = weights, weights_custom = weights_custom,
                               se_mat = se_mat)
    # Return as a pseudo meta result: beta = z, se = 1 so z = beta/se
    return(list(
      coefficients = matrix(zc, ncol = 1),
      se = matrix(1, nrow = length(zc), ncol = 1),
      tau2 = rep(NA_real_, length(zc)),
      I2 = rep(NA_real_, length(zc)),
      Q = rep(NA_real_, length(zc)),
      Q_df = rep(NA_real_, length(zc))
    ))
  }

  stop("fit_meta_chunk: No usable data found (need beta+se or t+df)", call. = FALSE)
}

#' Fit Meta-Analysis for a Single Voxel/ROI
#'
#' @param y Effect sizes
#' @param se Standard errors
#' @param X Design matrix
#' @param method Meta-analysis method
#' @param robust Robust estimation method
#' @param weights Weighting scheme
#' @param weights_custom Custom weights
#' @return List with meta-analysis results
#' @keywords internal
fit_meta_single <- function(y, se, X, method, robust, weights, weights_custom) {
  # Prepare data as matrices (single column for one voxel/ROI)
  Y <- matrix(y, ncol = 1)
  V <- matrix(se^2, ncol = 1)  # Convert SE to variance
  
  # Handle custom weights if provided
  if (weights == "equal") {
    # Equal weights - modify variances to be equal
    V[] <- 1
  } else if (weights == "custom" && !is.null(weights_custom)) {
    # Custom weights - convert to effective variances
    V <- matrix(1/weights_custom, ncol = 1)
  }
  
  # Set robust parameters
  huber_c <- if (robust == "huber") 1.345 else 1.345
  robust_iter <- if (robust == "huber") 2 else 0
  
  # Call C++ implementation
  result <- fmri_meta_fit(
    Y = Y,
    V = V,
    X = X,
    method = method,
    robust = robust,
    huber_c = huber_c,
    robust_iter = robust_iter,
    n_threads = 1  # Single threaded for single voxel
  )
  
  # Extract results for single voxel
  list(
    coefficients = as.vector(result$beta),
    se = as.vector(result$se),
    tau2 = result$tau2[1],
    I2 = result$I2_fe[1],
    Q = result$Q_fe[1],
    Q_df = result$df[1]
  )
}

#' Estimate Between-Study Variance (tau2)
#'
#' @param y Effect sizes
#' @param se Standard errors
#' @param X Design matrix
#' @param method Estimation method
#' @return tau2 estimate
#' @keywords internal
estimate_tau2 <- function(y, se, X, method) {
  n <- length(y)
  p <- ncol(X)
  w <- 1 / se^2
  
  # Fit fixed-effects model
  W <- diag(w)
  XtWX <- t(X) %*% W %*% X
  XtWy <- t(X) %*% W %*% y
  
  if (det(XtWX) > 1e-10) {
    coef <- solve(XtWX) %*% XtWy
    fitted <- X %*% coef
    residuals <- y - fitted
    Q <- sum(w * residuals^2)
  } else {
    return(0)
  }
  
  if (method == "dl") {
    # DerSimonian-Laird
    df <- n - p
    if (Q > df) {
      C <- sum(w) - sum(w^2) / sum(w)
      tau2 <- (Q - df) / C
    } else {
      tau2 <- 0
    }
  } else if (method == "pm") {
    # Paule-Mandel (iterative)
    tau2 <- 0
    for (iter in 1:50) {
      w_new <- 1 / (se^2 + tau2)
      W_new <- diag(w_new)
      XtWX_new <- t(X) %*% W_new %*% X
      XtWy_new <- t(X) %*% W_new %*% y
      
      if (det(XtWX_new) > 1e-10) {
        coef_new <- solve(XtWX_new) %*% XtWy_new
        fitted_new <- X %*% coef_new
        residuals_new <- y - fitted_new
        Q_new <- sum(w_new * residuals_new^2)
        
        # Update tau2
        tau2_new <- tau2 * Q_new / (n - p)
        
        if (abs(tau2_new - tau2) < 1e-6) {
          break
        }
        tau2 <- tau2_new
      } else {
        break
      }
    }
  } else if (method == "reml") {
    # Simplified REML (would use metafor for full implementation)
    tau2 <- max(0, var(y) - mean(se^2))
  }
  
  return(max(0, tau2))
}

#' Apply Huber Robust Weights
#'
#' @param y Effect sizes
#' @param se Standard errors
#' @param X Design matrix
#' @param w Initial weights
#' @return Updated weights
#' @keywords internal
apply_huber_weights <- function(y, se, X, w) {
  # Fit initial model
  W <- diag(w)
  XtWX <- t(X) %*% W %*% X
  XtWy <- t(X) %*% W %*% y
  
  if (det(XtWX) > 1e-10) {
    coef <- solve(XtWX) %*% XtWy
    fitted <- X %*% coef
    residuals <- y - fitted
    
    # Standardize residuals
    sigma <- sqrt(sum(w * residuals^2) / (length(y) - ncol(X)))
    std_residuals <- residuals / (sigma * sqrt(1/w))
    
    # Huber weights (k = 1.345 for 95% efficiency)
    k <- 1.345
    huber_w <- ifelse(abs(std_residuals) <= k, 1, k / abs(std_residuals))
    
    # Update weights
    w <- w * huber_w
  }
  
  return(w)
}
