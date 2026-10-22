#' fmrireg t-tests for Group Analysis
#'
#' Performs group-level t-tests on fMRI data with support for one-sample, two-sample,
#' paired, and ANCOVA designs. Provides both meta-analysis and classical t-test engines,
#' mirroring AFNI 3dttest++ functionality within the fmrireg framework. The
#' classic and Welch engines are fmrigds's `"ols:voxelwise"` and
#' `"ttest:welch"` reducers; the meta engine is fmrireg's meta-regression
#' (see [fmri_meta()]).
#'
#' @param gd A gds-backed `group_data()` object (or an fmrigds `gds`), a legacy
#'   group_data object, or a data frame accepted by [group_data()]. gds data
#'   must hold one subject-level contrast; every engine and option works on it,
#'   except that the meta engine on gds data supports neither `mu0` nor
#'   `voxelwise_cov`.
#' @param formula R formula for between-subjects design. Examples:
#'   \itemize{
#'     \item ~ 1: One-sample t-test
#'     \item ~ 1 + group: Two-sample t-test
#'     \item ~ 1 + group + age: ANCOVA with covariate
#'   }
#' @param engine Character string; analysis engine to use:
#'   \itemize{
#'     \item "auto": "meta" when sampling variances are available (always for
#'       gds data), otherwise "classic" (default). The fit records the engine
#'       that ran.
#'     \item "meta": Use inverse-variance meta-analysis (requires SE)
#'     \item "classic": Use OLS/Student t-test. A subject with a non-finite
#'       value for a feature is dropped from that feature only, so degrees of
#'       freedom can differ between features.
#'     \item "welch": Welch t-test for unequal variances; needs
#'       `formula = ~ 1 + group` with a two-level factor `group`. Non-finite
#'       subjects are dropped per feature.
#'   }
#' @param paired Logical; if TRUE, the effects are within-subject differences and
#'   are tested against `mu0` with the classic one-sample engine (default: FALSE)
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
#'   `weights_custom`). "equal"/"custom" are relative weights with an estimated
#'   residual scale (weighted least squares; see \code{\link{fmri_meta}}).
#'   Ignored for classic/welch engines.
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
  engine <- match.arg(engine)
  sign <- match.arg(sign)
  weights <- match.arg(weights)
  if (!is.null(mc)) {
    mc <- match.arg(mc, c("bh", "by", "spatial_fdr"))
  }
  formula <- stats::as.formula(formula)

  inp <- .fmri_ttest_input(gd, formula, mask)
  gd <- inp$gd
  Y <- inp$Y
  V <- inp$V
  X <- inp$X
  S <- nrow(X)
  if (!is.null(Y) && nrow(Y) != S) {
    stop("Design matrix rows must match data rows", call. = FALSE)
  }
  if (qr(X)$rank < ncol(X)) {
    stop("The design matrix for `formula` is rank deficient (", qr(X)$rank,
         " of ", ncol(X), " columns are estimable).", call. = FALSE)
  }
  group_info <- .fmri_ttest_group_term(X, inp$covars)

  # File-backed sources (gds, group_data_*) are read on demand; the meta engine
  # hands them to fmri_meta(). Every other engine works on the subjects x
  # features effect matrix.
  file_backed <- is.null(Y)
  use_meta <- !isTRUE(paired) &&
    (engine == "meta" || (engine == "auto" && (file_backed || !is.null(V))))
  if (!use_meta && engine == "meta") engine <- "classic"
  if (!use_meta && engine == "auto") engine <- "classic"
  if (isTRUE(paired)) {
    if (engine == "welch") {
      stop("paired = TRUE tests within-subject differences; it cannot use engine = 'welch'.",
           call. = FALSE)
    }
    message("paired=TRUE: performing one-sample test on differences vs mu0")
    X <- matrix(1, nrow = S, ncol = 1L, dimnames = list(NULL, "(Intercept)"))
    group_info <- NULL
  }
  if (!use_meta && file_backed) {
    Y <- .fmri_ttest_materialize_effects(gd, paste(engine, "engine"))
    if (nrow(Y) != S) stop("Design matrix rows must match data rows", call. = FALSE)
  }
  if (!is.null(Y)) Y <- Y - mu0
  P <- if (!is.null(Y)) ncol(Y) else NA_integer_

  # Every engine estimates the group coefficient as B - A (second level minus
  # first); results are reported in `sign`. Contrast weights are mapped onto
  # the raw coefficients and into that convention once, here.
  source_sign <- "BminusA"
  raw_weights <- NULL
  if (!is.null(contrast)) {
    canonical <- .fmri_ttest_resolve_contrast(
      contrast,
      coef_names = .fmri_ttest_canonical_coef_names(colnames(X), group_info),
      group_info = group_info
    )
    raw_weights <- .fmri_ttest_raw_contrast_weights(
      canonical, coef_names = colnames(X), group_info = group_info,
      target_sign = sign, source_sign = source_sign
    )
  }

  exact_contrast <- NULL
  analysis_method <- engine
  if (use_meta) {
    method <- .meta_resolve_method(getOption("fmrireg.meta.method", "pm"))
    robust <- getOption("fmrireg.meta.robust", "none")
    analysis_method <- if (!is.null(combine)) paste0("combine:", combine) else method
    if (!is.null(voxelwise_cov) && !is.null(contrast)) {
      stop("contrast is not supported with voxelwise_cov in meta analyses", call. = FALSE)
    }

    if (file_backed) {
      if (!isTRUE(all.equal(mu0, 0))) {
        stop("mu0 is not supported by the meta engine for file-backed group data; ",
             "use engine = 'classic'.", call. = FALSE)
      }
      if (!is.null(voxelwise_cov)) {
        stop("voxelwise_cov is not supported by the meta engine for file-backed group data; ",
             "use engine = 'classic'.", call. = FALSE)
      }
      meta_fit <- fmri_meta(
        gd, formula = formula,
        method = method, robust = robust,
        weights = weights, weights_custom = weights_custom,
        combine = combine,
        contrasts = if (!is.null(raw_weights)) {
          matrix(raw_weights, nrow = 1L, dimnames = list("contrast1", colnames(X)))
        },
        verbose = FALSE
      )
      if (inherits(meta_fit, "fmri_meta_list")) {
        stop("fmri_ttest(engine = 'meta') tests one subject-level contrast; these data have ",
             length(meta_fit), " (", paste(names(meta_fit), collapse = ", "),
             "). Use fmri_meta(), which fits each of them.", call. = FALSE)
      }
      fit <- list(beta = t(meta_fit$coefficients), se = t(meta_fit$se))
      if (!is.null(meta_fit$contrasts)) {
        exact_contrast <- .fmri_ttest_z_contrast(
          meta_fit$contrasts$estimate[, 1], meta_fit$contrasts$se[, 1]
        )
      }
    } else if (!is.null(voxelwise_cov)) {
      if (!identical(weights, "ivw")) {
        stop("weights = '", weights, "' is not supported with voxelwise_cov in meta analyses",
             call. = FALSE)
      }
      if (is.null(V)) {
        stop("engine = 'meta' with voxelwise_cov requires sampling variances", call. = FALSE)
      }
      fit <- fmri_meta_fit_extended(
        Y = Y, V = V, X = X,
        method = method, robust = robust,
        voxelwise = voxelwise_cov,
        center_voxelwise = center_voxelwise,
        voxel_name = voxel_name,
        n_threads = .fmrireg_num_threads()
      )
    } else {
      # One fitter for every weighting scheme: "ivw" uses V as absolute
      # sampling variances; "equal"/"custom" are relative weights with an
      # estimated residual scale (weighted least squares).
      if (is.null(V) && identical(weights, "ivw")) {
        stop("engine = 'meta' with weights = 'ivw' requires sampling variances", call. = FALSE)
      }
      fit <- .meta_fit_core(
        Y = Y, V = if (is.null(V)) matrix(1, nrow(Y), ncol(Y)) else V, X = X,
        method = method, robust = robust,
        weights = weights, weights_custom = weights_custom,
        n_threads = .fmrireg_num_threads()
      )
      if (!is.null(raw_weights)) {
        cs <- .meta_contrasts_from_tri(fit$beta, fit$cov_tri, matrix(raw_weights, ncol = 1))
        exact_contrast <- .fmri_ttest_z_contrast(cs$estimate[1, ], cs$se[1, ])
      }
    }
    rn <- rownames(fit$beta) %||% colnames(X)
    res <- list(beta = fit$beta, se = fit$se)
    res$z <- res$beta / res$se
    res$p <- 2 * stats::pnorm(abs(res$z), lower.tail = FALSE)
    res$df <- matrix(Inf, nrow = nrow(res$z), ncol = ncol(res$z))
    for (nm in c("beta", "se", "z", "p", "df")) rownames(res[[nm]]) <- rn
    P <- ncol(res$beta)

  } else if (engine == "welch") {
    if (is.null(group_info) || ncol(X) != 2L) {
      stop("engine = 'welch' requires a two-level factor named 'group' and no ",
           "other terms (formula = ~ 1 + group).", call. = FALSE)
    }
    if (!is.null(voxelwise_cov)) {
      stop("voxelwise_cov is not supported by engine = 'welch'", call. = FALSE)
    }
    res <- .fmri_ttest_welch(Y, inp$covars$group, coef_names = colnames(X))
    if (!is.null(raw_weights)) {
      if (raw_weights[[1]] != 0) {
        stop("Welch contrasts currently support only the group coefficient", call. = FALSE)
      }
      exact_contrast <- .fmri_ttest_single_coef_contrast(res, raw_weights)
    }

  } else {
    # Student t / ANCOVA: fmrigds's voxelwise OLS, or fmrireg's fitter when a
    # voxelwise covariate adds a per-feature design column.
    if (!is.null(voxelwise_cov)) {
      ols <- fmri_ols_fit(Y, X, voxelwise = voxelwise_cov,
                          center_voxelwise = center_voxelwise,
                          voxel_name = voxel_name)
    } else {
      ols <- .fmri_ttest_ols(Y, X)
    }
    res <- .fmri_ttest_ols_stats(ols, P)
    if (!is.null(raw_weights)) {
      if (is.null(voxelwise_cov)) {
        cs <- .meta_contrasts_from_tri(ols$beta, ols$cov_tri, matrix(raw_weights, ncol = 1))
        exact_contrast <- .fmri_ttest_t_contrast(cs$estimate[1, ], cs$se[1, ], ols$df)
      } else {
        if (sum(raw_weights != 0) > 1L) {
          stop("contrast with multiple coefficients is not supported when voxelwise_cov is used",
               call. = FALSE)
        }
        exact_contrast <- .fmri_ttest_single_coef_contrast(res, raw_weights)
      }
    }
  }

  if (!is.null(inp$sample_labels) && length(inp$sample_labels) == P) {
    for (nm in c("beta", "se", "t", "z", "p", "df")) {
      if (!is.null(res[[nm]])) colnames(res[[nm]]) <- inp$sample_labels
    }
  }
  res <- .fmri_ttest_normalize_group_rows(res, group_info)
  res <- .fmri_ttest_store_contrast(res, exact_contrast)
  res <- .fmri_ttest_apply_group_sign(res, group_info, target_sign = sign, source_sign = source_sign)
  res <- .fmri_ttest_apply_mc(res, mc, alpha, inp$feature_group)

  res$call <- match.call()
  res$formula <- formula
  res$engine <- if (use_meta) "meta" else engine
  res$method <- analysis_method
  res$combine <- if (use_meta) combine else NULL
  res$weights <- if (use_meta) weights else NULL
  res$n_subjects <- S
  res$n_features <- P
  res$mc <- mc
  res$alpha <- alpha
  res$sign <- sign
  res$roi_names <- inp$sample_labels
  if (!is.null(group_info)) res$group_levels <- group_info$levels

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
