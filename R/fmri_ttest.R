#' fmrireg t-tests for Group Analysis
#'
#' Performs group-level t-tests on fMRI data with support for one-sample, two-sample,
#' paired, and ANCOVA designs. Provides both meta-analysis and classical t-test engines,
#' mirroring AFNI 3dttest++ functionality within the fmrireg framework.
#'
#' @param gd A group_data object or data frame with subject data
#' @param formula R formula for between-subjects design. Examples:
#'   \itemize{
#'     \item ~ 1: One-sample t-test
#'     \item ~ 1 + group: Two-sample t-test
#'     \item ~ 1 + group + age: ANCOVA with covariate
#'   }
#' @param engine Character string; analysis engine to use:
#'   \itemize{
#'     \item "auto": Automatically select based on data (default)
#'     \item "meta": Use inverse-variance meta-analysis (requires SE)
#'     \item "classic": Use OLS/Student t-test
#'     \item "welch": Use Welch t-test for unequal variances
#'   }
#' @param paired Logical; if TRUE, perform paired t-test on differences (default: FALSE)
#' @param mu0 Numeric scalar; constant to test against for one-sample tests (default: 0)
#' @param contrast Optional named numeric vector for linear combination of coefficients
#' @param mc Character string or NULL; multiple comparisons correction:
#'   \itemize{
#'     \item NULL: No correction (default)
#'     \item "bh": Benjamini-Hochberg FDR
#'     \item "by": Benjamini-Yekutieli FDR
#'     \item "spatial_fdr": Spatially-aware FDR
#'   }
#' @param alpha Numeric scalar; FDR level if mc is not NULL (default: 0.05)
#' @param sign Character string; sign convention for group differences:
#'   \itemize{
#'     \item "AminusB": A - B (default)
#'     \item "BminusA": B - A
#'   }
#' @param mask Optional mask object or path
#' @param voxelwise_cov Optional S x P matrix of voxelwise covariates
#' @param center_voxelwise Logical; center voxelwise covariate per feature (default: TRUE)
#' @param voxel_name Character string; name for voxelwise coefficient (default: "voxel_cov")
#' @param weights Character string for meta-engine weighting: "ivw" (default,
#'   uses provided SE), "equal" (equal weights), or "custom" (supply
#'   `weights_custom`). Ignored for classic/welch engines.
#' @param weights_custom Numeric vector (length S) or matrix (S x P) of custom
#'   weights when `weights = "custom"`.
#' @param combine Optional. When using the meta engine with t-only inputs
#'   (i.e., per-subject t-statistics and df), specify the t-combination
#'   method: "stouffer", "fisher", or "lancaster". Passed through to
#'   `fmri_meta()` when delegating to the meta engine on group_data_*.
#'
#' @return An fmri_ttest_fit object containing:
#'   \itemize{
#'     \item beta: Matrix of coefficients
#'     \item se: Matrix of standard errors (if available)
#'     \item t: Matrix of t-statistics (classic engine)
#'     \item z: Matrix of z-scores
#'     \item p: Matrix of p-values
#'     \item df: Matrix of degrees of freedom
#'     \item q: Matrix of FDR-adjusted p-values (if mc is used)
#'     \item z_contrast: Vector of contrast z-scores (if contrast is used)
#'     \item p_contrast: Vector of contrast p-values (if contrast is used)
#'   }
#'
#' @examples
#' \dontrun{
#' # One-sample t-test
#' fit <- fmri_ttest(gd, formula = ~ 1)
#' 
#' # Two-sample t-test
#' fit <- fmri_ttest(gd, formula = ~ 1 + group)
#' 
#' # Paired t-test
#' fit <- fmri_ttest(gd_diff, formula = ~ 1, paired = TRUE)
#' 
#' # ANCOVA with age covariate
#' fit <- fmri_ttest(gd, formula = ~ 1 + group + age)
#' 
#' # With spatial FDR correction
#' fit <- fmri_ttest(gd, formula = ~ 1 + group, mc = "spatial_fdr", alpha = 0.05)
#' }
#'
#' @export
fmri_ttest <- function(gd,
                      formula = ~ 1,
                      engine = c("auto", "meta", "classic", "welch"),
                      paired = FALSE,
                      mu0 = 0,
                      contrast = NULL,
                      mc = NULL,
                      alpha = 0.05,
                      sign = c("AminusB", "BminusA"),
                      mask = NULL,
                      voxelwise_cov = NULL,
                      center_voxelwise = TRUE,
                      voxel_name = "voxel_cov",
                      weights = c("ivw", "equal", "custom"),
                      weights_custom = NULL,
                      combine = NULL) {
  # Fast-path for gds-backed objects: delegate to fmrigds reducers
  if (inherits(gd, "gds") || inherits(gd, "group_data_gds")) {
    gds_fn <- get0("fmri_ttest.group_data_gds", envir = asNamespace("fmrireg"),
                   ifnotfound = NULL)
    if (is.null(gds_fn) && requireNamespace("fmrigds", quietly = TRUE)) {
      gds_fn <- get0("fmri_ttest.group_data_gds",
                     envir = asNamespace("fmrigds"), ifnotfound = NULL)
    }
    if (!is.null(gds_fn)) {
      return(gds_fn(gd,
                    formula = formula,
                    engine = match.arg(engine),
                    paired = paired,
                    mu0 = mu0,
                    contrast = contrast,
                    mc = mc,
                    alpha = alpha,
                    sign = match.arg(sign),
                    mask = mask,
                    voxelwise_cov = voxelwise_cov,
                    center_voxelwise = center_voxelwise,
                    voxel_name = voxel_name,
                    weights = match.arg(weights),
                    weights_custom = weights_custom,
                    combine = combine))
    }
    stop("GDS-backed group_data requires the 'fmrigds' package.", call. = FALSE)
  }
  
  engine <- match.arg(engine)
  sign <- match.arg(sign)
  weights <- match.arg(weights)
  
  if (!is.null(mc)) {
    mc <- match.arg(mc, c("bh", "by", "spatial_fdr"))
  }
  
  # 1) Normalize input to matrices Y (S x P), optional V (S x P), covars
  if (inherits(gd, "group_data")) {
    if (!is.null(gd$blocks)) {
      if (length(gd$blocks) != 1) {
        stop("fmri_ttest expects a single contrast block", call. = FALSE)
      }
      B <- gd$blocks[[1]]
      Y <- B$Y  # S x P effects (betas or derived)
      V <- if (!is.null(B$V)) B$V else NULL  # S x P sampling variances if available
      covars <- B$covars
      feature_group <- .fmri_ttest_feature_group(gd, B)
    } else {
      # Modern group_data_* path: leave Y/V unset and use stored covariates
      B <- NULL
      Y <- NULL
      V <- NULL
      covars <- gd$covariates
      feature_group <- .fmri_ttest_feature_group(gd, NULL)
    }
  } else {
    # Convert data frame to group_data
    if (!inherits(gd, "data.frame")) {
      stop("gd must be a group_data object or data frame", call. = FALSE)
    }
    gd <- group_data(gd, mask = mask)
    B <- gd$blocks[[1]]
    Y <- B$Y
    V <- if (!is.null(B$V)) B$V else NULL
    covars <- B$covars
    feature_group <- .fmri_ttest_feature_group(gd, B)
  }
  
  # 2) Subtract mu0 from Y for testing against constant
  if (!is.null(Y)) Y <- Y - mu0
  
  # 3) Build design matrix X from formula
  if (is.null(covars)) {
    covars <- data.frame(intercept = rep(1, if (!is.null(Y)) nrow(Y) else n_subjects(gd)))
  }
  
  X <- stats::model.matrix(stats::as.formula(formula), data = covars)
  storage.mode(X) <- "double"
  group_info <- .fmri_ttest_group_term(X, covars)
  K <- ncol(X)
  S <- nrow(X)
  P <- if (!is.null(Y)) ncol(Y) else NA_integer_
  
  if (!is.null(Y) && S != nrow(Y)) {
    stop("Design matrix rows must match data rows", call. = FALSE)
  }
  
  # 4) Choose engine
  use_meta <- if (engine == "auto") {
    !is.null(V) || (inherits(gd, "group_data") && is.null(Y))
  } else {
    engine == "meta"
  }
  
  res <- list()
  raw_contrast_weights <- NULL
  canonical_contrast_weights <- NULL
  exact_contrast <- NULL
  
  if (paired) {
    # Paired: assume Y contains within-subject differences
    message("paired=TRUE: performing one-sample test on differences vs mu0")
    if (is.null(Y)) {
      Y <- .fmri_ttest_materialize_effects(gd, "paired analyses")
      P <- ncol(Y)
      Y <- Y - mu0
    }
    
    # One-sample test (intercept only)
    Xp <- matrix(1, nrow = nrow(Y), ncol = 1, dimnames = list(NULL, "(Intercept)"))
    
    if (!is.null(voxelwise_cov)) {
      ols <- fmri_ols_fit(Y, Xp, voxelwise = voxelwise_cov,
                          center_voxelwise = center_voxelwise,
                          voxel_name = voxel_name)
    } else {
      ols <- ols_t_cpp(Y, Xp)
      rownames(ols$beta) <- "(Intercept)"
      rownames(ols$se) <- "(Intercept)"
      rownames(ols$t) <- "(Intercept)"
    }
    
    T <- ols$t[1, , drop = FALSE]
    df <- rep(ols$df, P)
    Pval <- 2 * stats::pt(abs(as.numeric(T)), df = df, lower.tail = FALSE)
    Z <- stats::qnorm(pmax(1e-300, 1 - Pval/2)) * sign(as.numeric(T))
    
    res$beta <- ols$beta
    res$se <- ols$se
    res$t <- T
    res$df <- matrix(df, nrow = 1, ncol = P)
    res$z <- matrix(Z, nrow = 1, ncol = P)
    res$p <- matrix(Pval, nrow = 1, ncol = P)
    
    rownames(res$beta) <- "(Intercept)"
    rownames(res$se) <- "(Intercept)"
    rownames(res$t) <- "(Intercept)"
    rownames(res$df) <- "(Intercept)"
    rownames(res$z) <- "(Intercept)"
    rownames(res$p) <- "(Intercept)"
    
  } else if (use_meta) {
    # Meta-regression using existing infrastructure
    method <- getOption("fmrireg.meta.method", "pm")
    robust <- getOption("fmrireg.meta.robust", "none")
    
    if (inherits(gd, "group_data") && is.null(Y)) {
      # Delegate to fmri_meta for group_data_* inputs
      if (!is.null(voxelwise_cov) && !is.null(contrast)) {
        stop("contrast is not supported with voxelwise_cov in meta analyses", call. = FALSE)
      }
      if (!is.null(contrast)) {
        canonical_contrast_weights <- .fmri_ttest_resolve_contrast(
          contrast,
          coef_names = .fmri_ttest_canonical_coef_names(colnames(X), group_info),
          group_info = group_info
        )
        raw_contrast_weights <- .fmri_ttest_raw_contrast_weights(
          canonical_contrast_weights,
          coef_names = colnames(X),
          group_info = group_info,
          target_sign = sign,
          source_sign = "BminusA"
        )
      }
      meta_fit <- fmri_meta(
        gd, formula = stats::as.formula(formula),
        method = method, robust = robust,
        weights = weights, weights_custom = weights_custom,
        contrasts = if (!is.null(raw_contrast_weights)) matrix(raw_contrast_weights, ncol = 1) else NULL,
        combine = combine,
        verbose = FALSE
      )
      res$beta <- t(meta_fit$coefficients)
      res$se   <- t(meta_fit$se)
      res$z    <- res$beta / res$se
      res$p    <- 2 * stats::pnorm(abs(res$z), lower.tail = FALSE)
      res$df   <- matrix(Inf, nrow = nrow(res$z), ncol = ncol(res$z))
      rownames(res$beta) <- colnames(X)
      rownames(res$se) <- colnames(X)
      rownames(res$z) <- colnames(X)
      rownames(res$df) <- colnames(X)
      if (!is.null(meta_fit$contrasts)) {
        exact_contrast <- list(
          estimate = as.numeric(meta_fit$contrasts$estimate[, 1]),
          se = as.numeric(meta_fit$contrasts$se[, 1]),
          z = as.numeric(meta_fit$contrasts$z[, 1]),
          p = 2 * stats::pnorm(abs(meta_fit$contrasts$z[, 1]), lower.tail = FALSE),
          df = rep(Inf, nrow(meta_fit$contrasts$estimate))
        )
      }
    } else if (!is.null(voxelwise_cov)) {
      if (!is.null(contrast)) {
        stop("contrast is not supported with voxelwise_cov in meta analyses", call. = FALSE)
      }
      # Update meta_fit wrapper to handle voxelwise
      out <- fmri_meta_fit_extended(
        Y = Y, V = V, X = X,
        method = method, robust = robust,
        voxelwise = voxelwise_cov,
        center_voxelwise = center_voxelwise,
        voxel_name = voxel_name,
        n_threads = getOption("fmrireg.num_threads", 0)
      )
    } else {
      # Apply weighting variants for meta
      V_eff <- V
      if (weights == "equal") {
        V_eff[] <- 1
      } else if (weights == "custom") {
        weights_custom <- .fmri_ttest_validate_weights_custom(weights_custom, nrow(V_eff), ncol(V_eff))
        if (is.vector(weights_custom)) {
          V_eff <- matrix(1 / weights_custom, nrow = nrow(V_eff), ncol = ncol(V_eff))
        } else {
          V_eff <- 1 / weights_custom
        }
      }
      if (!is.null(contrast)) {
        canonical_contrast_weights <- .fmri_ttest_resolve_contrast(
          contrast,
          coef_names = .fmri_ttest_canonical_coef_names(colnames(X), group_info),
          group_info = group_info
        )
        raw_contrast_weights <- .fmri_ttest_raw_contrast_weights(
          canonical_contrast_weights,
          coef_names = colnames(X),
          group_info = group_info,
          target_sign = sign,
          source_sign = "BminusA"
        )
        out <- fmri_meta_fit_contrasts(
          Y = Y, V = V_eff, X = X, Cmat = matrix(raw_contrast_weights, ncol = 1),
          method = method, robust = robust,
          n_threads = getOption("fmrireg.num_threads", 0)
        )
        exact_contrast <- list(
          estimate = as.numeric(out$c_beta[1, ]),
          se = as.numeric(out$c_se[1, ]),
          z = as.numeric(out$c_z[1, ]),
          p = 2 * stats::pnorm(abs(out$c_z[1, ]), lower.tail = FALSE),
          df = rep(Inf, ncol(out$c_beta))
        )
      } else {
        out <- fmri_meta_fit(
          Y = Y, V = V_eff, X = X,
          method = method, robust = robust,
          n_threads = getOption("fmrireg.num_threads", 0)
        )
      }
    }
    
    if (exists("out")) {
      res$beta <- out$beta
      res$se <- out$se
      res$z <- out$z
      res$p <- 2 * stats::pnorm(abs(out$z), lower.tail = FALSE)
      res$df <- matrix(Inf, nrow = nrow(out$z), ncol = ncol(out$z))
      if (is.null(rownames(res$beta))) rownames(res$beta) <- rownames(out$z) %||% colnames(X)
      if (is.null(rownames(res$se))) rownames(res$se) <- rownames(res$beta)
      if (is.null(rownames(res$z))) rownames(res$z) <- rownames(res$beta)
      rownames(res$df) <- rownames(res$beta)
    }
    
  } else {
    # Classic Student t or Welch
    # If Y is not yet materialized (e.g., group_data_nifti/h5), load betas
    if (is.null(Y)) {
      Y <- .fmri_ttest_materialize_effects(gd, "classic/welch engines")
      P <- ncol(Y)
      Y <- Y - mu0
    }
    has_group <- "group" %in% colnames(covars) && 
                 is.factor(covars$group) && 
                 nlevels(covars$group) == 2
    
    if (engine == "welch" && has_group && K <= 2) {
      # Welch t-test
      g <- as.integer(covars$group)
      w <- welch_t_cpp(Y, g)
      
      T <- as.numeric(w$t)
      df <- as.numeric(w$df)
      Pval <- 2 * stats::pt(abs(T), df = df, lower.tail = FALSE)
      Z <- stats::qnorm(pmax(1e-300, 1 - Pval/2)) * sign(T)
      
      res$beta <- rbind(
        "(Intercept)" = (w$muA + w$muB) / 2,
        group = w$muA - w$muB
      )
      rownames(res$beta) <- c("(Intercept)", "group")
      res$se <- matrix(NA_real_, nrow = 2, ncol = P,
                       dimnames = list(c("(Intercept)", "group"), NULL))
      res$t <- rbind("(Intercept)" = NA_real_, group = T)
      res$df <- rbind("(Intercept)" = NA_real_, group = df)
      res$z <- rbind("(Intercept)" = NA_real_, group = Z)
      res$p <- rbind("(Intercept)" = NA_real_, group = Pval)
      if (!is.null(contrast)) {
        canonical_contrast_weights <- .fmri_ttest_resolve_contrast(
          contrast,
          coef_names = rownames(res$beta),
          group_info = list(raw_name = "group", canonical_name = "group")
        )
        nz <- which(abs(canonical_contrast_weights) > 0)
        if (length(nz) > 1L || !identical(names(canonical_contrast_weights)[nz], "group")) {
          stop("Welch contrasts currently support only the group coefficient", call. = FALSE)
        }
        exact_contrast <- .fmri_ttest_single_coef_contrast(res, canonical_contrast_weights)
      }
      
    } else {
      # General OLS t-test/ANCOVA
      if (!is.null(voxelwise_cov)) {
        ols <- fmri_ols_fit(Y, X, voxelwise = voxelwise_cov,
                           center_voxelwise = center_voxelwise,
                           voxel_name = voxel_name)
      } else {
        ols <- ols_t_cpp(Y, X)
        rownames(ols$beta) <- colnames(X)
        rownames(ols$se) <- colnames(X)
        rownames(ols$t) <- colnames(X)
      }
      
      # Convert t to z and p
      df <- rep(ols$df, P)
      Pval <- apply(ols$t, 1, function(tt) {
        2 * stats::pt(abs(tt), df = df, lower.tail = FALSE)
      })
      Pval <- t(Pval)
      rownames(Pval) <- rownames(ols$t)
      
      Z <- matrix(NA_real_, nrow = nrow(ols$t), ncol = P,
                  dimnames = dimnames(ols$t))
      for (i in seq_len(nrow(Z))) {
        Zi <- stats::qnorm(pmax(1e-300, 1 - Pval[i,]/2)) * sign(ols$t[i,])
        Z[i,] <- Zi
      }
      
      res$beta <- ols$beta
      res$se <- ols$se
      res$t <- ols$t
      res$df <- matrix(ols$df, nrow = nrow(ols$t), ncol = ncol(ols$t),
                       byrow = TRUE, dimnames = dimnames(ols$t))
      res$z <- Z
      res$p <- Pval
      if (!is.null(contrast)) {
        canonical_coef_names <- rownames(res$beta)
        if (!is.null(group_info)) {
          canonical_coef_names[canonical_coef_names == group_info$raw_name] <- group_info$canonical_name
        }
        canonical_contrast_weights <- .fmri_ttest_resolve_contrast(
          contrast,
          coef_names = canonical_coef_names,
          group_info = group_info
        )
        if (!is.null(voxelwise_cov) && sum(abs(canonical_contrast_weights) > 0) > 1L) {
          stop("contrast with multiple coefficients is not supported when voxelwise_cov is used", call. = FALSE)
        }
        if (is.null(voxelwise_cov)) {
          raw_contrast_weights <- .fmri_ttest_raw_contrast_weights(
            canonical_contrast_weights,
            coef_names = rownames(ols$beta),
            group_info = group_info,
            target_sign = sign,
            source_sign = "BminusA"
          )
          exact_contrast <- .fmri_ttest_exact_ols_contrast(ols, X, raw_contrast_weights)
        } else {
          exact_contrast <- .fmri_ttest_single_coef_contrast(
            res,
            stats::setNames(canonical_contrast_weights, canonical_coef_names)
          )
        }
      }
    }
  }
  
  res <- .fmri_ttest_normalize_group_rows(res, group_info)
  if (!is.null(exact_contrast)) {
    res <- .fmri_ttest_store_contrast(res, exact_contrast)
  }
  source_sign <- if (use_meta || (!paired && engine != "welch")) "BminusA" else "AminusB"
  res <- .fmri_ttest_apply_group_sign(res, group_info, target_sign = sign, source_sign = source_sign)

  res <- .fmri_ttest_apply_mc(res, mc, alpha, feature_group)

  # Store metadata
  res$call <- match.call()
  res$formula <- formula
  res$engine <- if (use_meta) "meta" else engine
  res$n_subjects <- S
  res$n_features <- P
  res$mc <- mc
  res$alpha <- alpha
  if (!is.null(group_info)) {
    res$group_levels <- group_info$levels
    res$sign <- sign
  }
  
  class(res) <- c("fmri_ttest_fit", "list")
  res
}

#' Print method for fmri_ttest_fit
#' 
#' @param x An fmri_ttest_fit object
#' @param ... Additional print arguments
#' @return Invisibly returns the input object x
#' @export
print.fmri_ttest_fit <- function(x, ...) {
  cat("fmri_ttest Results\n")
  cat("==================\n")
  cat("Formula:", deparse(x$formula), "\n")
  cat("Engine:", x$engine, "\n")
  cat("Subjects:", x$n_subjects, "\n")
  cat("Features:", x$n_features, "\n")
  
  if (!is.null(x$beta)) {
    cat("\nCoefficients:\n")
    print(summary(x$beta))
  }
  
  invisible(x)
}

#' Summary method for fmri_ttest_fit
#' 
#' @param object An fmri_ttest_fit object
#' @param ... Additional summary arguments
#' @return Invisibly returns the input object
#' @export
summary.fmri_ttest_fit <- function(object, ...) {
  print(object)
  
  if (!is.null(object$q)) {
    cat("\nMultiple comparisons correction applied\n")
    cat("Method:", object$mc %||% "unknown", "\n")
    cat("Significant features (FDR <", object$alpha %||% 0.05, "):\n")
    for (i in seq_len(nrow(object$q))) {
      n_sig <- sum(object$q[i,] < (object$alpha %||% 0.05), na.rm = TRUE)
      cat("  ", rownames(object$q)[i], ":", n_sig, "\n")
    }
  }
  
  invisible(object)
}
