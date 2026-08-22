# Typed configuration specifications for fmri_lm().

.fmri_lm_scalar <- function(x, name, type = NULL, finite = FALSE) {
  if (length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be a single non-missing value.", name), call. = FALSE)
  }
  if (!is.null(type) && !isTRUE(do.call(type, list(x)))) {
    stop(sprintf("`%s` has the wrong type.", name), call. = FALSE)
  }
  if (finite && !is.finite(x)) {
    stop(sprintf("`%s` must be finite.", name), call. = FALSE)
  }
  x
}

.fmri_lm_flag <- function(x, name) {
  .fmri_lm_scalar(x, name, "is.logical")
  x
}

.fmri_lm_number <- function(x, name, lower = -Inf, strictly = FALSE,
                            integer = FALSE) {
  .fmri_lm_scalar(x, name, "is.numeric", finite = TRUE)
  if ((strictly && x <= lower) || (!strictly && x < lower)) {
    relation <- if (strictly) ">" else ">="
    stop(sprintf("`%s` must be %s %s.", name, relation, lower), call. = FALSE)
  }
  if (integer && x != as.integer(x)) {
    stop(sprintf("`%s` must be an integer.", name), call. = FALSE)
  }
  if (integer) as.integer(x) else as.numeric(x)
}

.fmri_lm_choice <- function(x, choices, name) {
  .fmri_lm_scalar(x, name, "is.character")
  if (!x %in% choices) {
    stop(sprintf("`%s` must be one of: %s.", name,
                 paste0("'", choices, "'", collapse = ", ")), call. = FALSE)
  }
  x
}

.fmri_lm_check_keys <- function(x, allowed, name) {
  if (is.null(x)) return(invisible(TRUE))
  if (!is.list(x)) {
    stop(sprintf("`%s` must be a list.", name), call. = FALSE)
  }
  if (length(x) && (is.null(names(x)) || any(!nzchar(names(x))))) {
    stop(sprintf("Every element of `%s` must be named.", name), call. = FALSE)
  }
  unknown <- setdiff(names(x), allowed)
  if (length(unknown)) {
    stop(sprintf("Unknown %s %s: %s.", name,
                 if (length(unknown) == 1L) "key" else "keys",
                 paste0("`", unknown, "`", collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

.new_fmri_lm_spec <- function(x, class) {
  structure(x, class = c(class, "fmri_lm_spec"))
}

#' Statistical fitting scope
#'
#' @param scope Fit a single joint model or fit runs separately and combine
#'   them by fixed-effect meta-analysis.
#' @param run_pooling Runwise combination rule. Only fixed-effect pooling is
#'   currently supported.
#' @export
estimation_spec <- function(scope = c("joint", "runwise_meta"),
                            run_pooling = "fixed_effect") {
  scope <- match.arg(scope)
  run_pooling <- .fmri_lm_choice(run_pooling, "fixed_effect", "run_pooling")
  .new_fmri_lm_spec(list(scope = scope, run_pooling = run_pooling),
                    "fmri_lm_estimation_spec")
}

#' Temporal-noise specification
#'
#' @param struct One of iid, ar1, ar2, or arp.
#' @param p AR order for `struct = "arp"`.
#' @param q Nonnegative moving-average order. Positive values request an
#'   ARMA(p, q) model, with the AR order determined by `struct` and `p`.
#'   The current MA-capable path supports runwise meta-estimation with run
#'   pooling, without censoring, parcel pooling, or voxelwise covariance
#'   estimation.
#' @param iter_gls Maximum number of GLS refinement iterations. Standard
#'   design-corrected AR fitting uses the initial OLS residuals once and holds
#'   that estimate fixed for the GLS solve. Iterative refinement remains
#'   available to ARMA models (`q > 0`).
#' @param pooling Temporal covariance pooling scope. With built-in shared AR
#'   estimation, `"run"` estimates one coefficient vector per run and
#'   `"global"` pools across runs.
#' @param shared_estimator Spatial estimator for a shared temporal covariance.
#'   `"pooled_acvf"` (the default) pools residual autocovariances across voxels
#'   and therefore targets a typical voxel covariance. `"mean_series"` first
#'   averages residual values across voxels and targets the coherent spatial
#'   component; it is experimental, can be much more autocorrelated than an
#'   individual voxel. Because OLS projection is linear across response
#'   columns, the matching design correction remains valid after averaging.
#' @param parcels Optional parcel labels for parcel pooling.
#' @param voxelwise Estimate temporal covariance separately by voxel. In the
#'   built-in fitter this currently requires AR-only runwise meta-estimation,
#'   `pooling = "run"`, `iter_gls = 1`, no censoring, and no volume weighting
#'   or soft-subspace projection. Robust fitting is supported, but robust AR
#'   re-estimation is not. Registered engines may define broader capabilities.
#' @param exact_first Use exact first-observation AR scaling.
#' @param censor Optional censor indices, logical mask, or `"auto"`.
#' @param shrink_c0 Parcel shrinkage constant used by supporting engines.
#' @export
noise_spec <- function(struct = c("iid", "ar1", "ar2", "arp"),
                       p = NULL, q = 0L, iter_gls = 1L,
                       pooling = c("run", "global", "parcel"),
                       shared_estimator = c("pooled_acvf", "mean_series"),
                       parcels = NULL, voxelwise = FALSE,
                       exact_first = FALSE, censor = NULL,
                       shrink_c0 = 100L) {
  struct <- match.arg(struct)
  pooling <- match.arg(pooling)
  shared_estimator <- match.arg(shared_estimator)
  q <- .fmri_lm_number(q, "q", lower = 0, integer = TRUE)
  if (is.null(p)) {
    p <- switch(struct, iid = NULL, ar1 = 1L, ar2 = 2L, arp = NULL)
  } else {
    p <- .fmri_lm_number(p, "p", lower = 1, integer = TRUE)
  }
  if (struct == "arp" && is.null(p)) {
    stop("`p` must be supplied when `struct = 'arp'`.", call. = FALSE)
  }
  iter_gls <- .fmri_lm_number(iter_gls, "iter_gls", lower = 0, integer = TRUE)
  voxelwise <- .fmri_lm_flag(voxelwise, "voxelwise")
  exact_first <- .fmri_lm_flag(exact_first, "exact_first")
  shrink_c0 <- .fmri_lm_number(shrink_c0, "shrink_c0", lower = 0,
                               strictly = TRUE, integer = TRUE)
  if (!is.null(censor)) {
    valid <- (is.character(censor) && length(censor) == 1L && identical(censor, "auto")) ||
      is.numeric(censor) || is.logical(censor)
    if (!valid || anyNA(censor)) {
      stop("`censor` must be NULL, 'auto', integer indices, or a logical mask.",
           call. = FALSE)
    }
  }
  if (q > 0L) {
    if (!identical(pooling, "run")) {
      stop("MA terms currently require `pooling = 'run'`.", call. = FALSE)
    }
    if (voxelwise) {
      stop("MA terms do not yet support voxelwise temporal covariance.",
           call. = FALSE)
    }
    if (!is.null(censor)) {
      stop("MA terms do not yet support censoring because estimation across gaps is biased.",
           call. = FALSE)
    }
    if (iter_gls < 1L) {
      stop("MA terms require `iter_gls >= 1`.", call. = FALSE)
    }
  }
  .new_fmri_lm_spec(list(
    struct = struct, p = p, q = q, iter_gls = iter_gls,
    pooling = pooling, shared_estimator = shared_estimator,
    parcels = parcels, voxelwise = voxelwise,
    exact_first = exact_first, censor = censor, shrink_c0 = shrink_c0,
    global = identical(pooling, "global"),
    by_cluster = identical(pooling, "parcel")
  ), "fmri_lm_noise_spec")
}

#' Robust-estimation specification
#'
#' @param type Robust loss, or `"none"` for ordinary least squares.
#' @param k_huber Tuning constant for the Huber loss.
#' @param c_tukey Tuning constant for Tukey's bisquare loss.
#' @param max_iter Maximum number of robust reweighting iterations.
#' @param scale_scope Scope at which robust scale is estimated.
#' @param reestimate_phi Re-estimate temporal-noise parameters after robust
#'   reweighting.
#' @export
robust_spec <- function(type = c("none", "huber", "bisquare"),
                        k_huber = 1.345, c_tukey = 4.685,
                        max_iter = 2L,
                        scale_scope = c("run", "global", "voxel"),
                        reestimate_phi = FALSE) {
  type <- match.arg(type)
  scale_scope <- match.arg(scale_scope)
  k_huber <- .fmri_lm_number(k_huber, "k_huber", lower = 0, strictly = TRUE)
  c_tukey <- .fmri_lm_number(c_tukey, "c_tukey", lower = 0, strictly = TRUE)
  max_iter <- .fmri_lm_number(max_iter, "max_iter", lower = 1, integer = TRUE)
  reestimate_phi <- .fmri_lm_flag(reestimate_phi, "reestimate_phi")
  .new_fmri_lm_spec(list(
    type = type, k_huber = k_huber, c_tukey = c_tukey,
    max_iter = max_iter, scale_scope = scale_scope,
    reestimate_phi = reestimate_phi
  ), "fmri_lm_robust_spec")
}

#' Variance and reference-distribution specification
#'
#' @param method Covariance estimator: model-based, HAC, or sandwich.
#' @param max_lag Maximum HAC lag, or `"auto"` for a data-dependent lag.
#' @param taper HAC lag-window taper.
#' @param debias Apply the estimator's finite-sample correction.
#' @param df Reference degrees-of-freedom method.
#' @export
variance_spec <- function(method = c("model", "hac", "sandwich"),
                          max_lag = "auto",
                          taper = c("none", "tukey", "parzen"),
                          debias = TRUE,
                          df = c("residual", "satterthwaite")) {
  method <- match.arg(method)
  taper <- match.arg(taper)
  df <- match.arg(df)
  debias <- .fmri_lm_flag(debias, "debias")
  if (is.character(max_lag)) {
    max_lag <- .fmri_lm_choice(max_lag, "auto", "max_lag")
  } else {
    max_lag <- .fmri_lm_number(max_lag, "max_lag", lower = 0, integer = TRUE)
  }
  .new_fmri_lm_spec(list(method = method, max_lag = max_lag, taper = taper,
                         debias = debias, df = df),
                    "fmri_lm_variance_spec")
}

#' Volume-weight specification
#'
#' @param method Volume-weight construction rule.
#' @param threshold Positive threshold or tuning constant for the selected rule.
#' @param values Optional explicit finite volume weights.
#' @export
weights_spec <- function(method = c("none", "inverse_squared", "soft_threshold", "tukey"),
                         threshold = 1.5, values = NULL) {
  method <- match.arg(method)
  threshold <- .fmri_lm_number(threshold, "threshold", lower = 0, strictly = TRUE)
  if (!is.null(values) && (!is.numeric(values) || any(!is.finite(values)))) {
    stop("`values` must be NULL or a finite numeric vector.", call. = FALSE)
  }
  .new_fmri_lm_spec(list(method = method, threshold = threshold, values = values,
                         enabled = method != "none", weights = values),
                    "fmri_lm_weights_spec")
}

#' Nuisance-projection specification
#'
#' @param method Projection rule.
#' @param nuisance_mask Optional logical or numeric selector for nuisance
#'   features.
#' @param nuisance_matrix Optional explicit nuisance-feature matrix.
#' @param lambda Nonnegative regularization strength, `"auto"`, or `"gcv"`.
#' @param warn_redundant Warn when projection overlaps modeled nuisance terms.
#' @export
projection_spec <- function(method = c("none", "soft_subspace"),
                            nuisance_mask = NULL, nuisance_matrix = NULL,
                            lambda = "auto", warn_redundant = TRUE) {
  method <- match.arg(method)
  warn_redundant <- .fmri_lm_flag(warn_redundant, "warn_redundant")
  if (is.character(lambda)) {
    lambda <- .fmri_lm_choice(lambda, c("auto", "gcv"), "lambda")
  } else {
    lambda <- .fmri_lm_number(lambda, "lambda", lower = 0)
  }
  if (method == "soft_subspace" && is.null(nuisance_mask) && is.null(nuisance_matrix)) {
    stop("Soft-subspace projection requires `nuisance_mask` or `nuisance_matrix`.",
         call. = FALSE)
  }
  .new_fmri_lm_spec(list(
    method = method, nuisance_mask = nuisance_mask,
    nuisance_matrix = nuisance_matrix, lambda = lambda,
    warn_redundant = warn_redundant, enabled = method != "none"
  ), "fmri_lm_projection_spec")
}

#' Mechanical execution specification
#'
#' @param voxel_chunks Number of voxel partitions for joint fitting.
#' @param backend Matrix or formula/reference implementation.
#' @param parallel Parallelization dimension.
#' @param progress Display progress.
#' @export
compute_spec <- function(voxel_chunks = 10L,
                         backend = c("matrix", "reference"),
                         parallel = c("none", "voxels", "chunks"),
                         progress = FALSE) {
  backend <- match.arg(backend)
  parallel <- match.arg(parallel)
  voxel_chunks <- .fmri_lm_number(voxel_chunks, "voxel_chunks", lower = 1,
                                  integer = TRUE)
  progress <- .fmri_lm_flag(progress, "progress")
  .new_fmri_lm_spec(list(
    voxel_chunks = voxel_chunks, backend = backend,
    parallel = parallel, progress = progress
  ), "fmri_lm_compute_spec")
}

#' @export
print.fmri_lm_spec <- function(x, ...) {
  cat(sprintf("<%s>\n", class(x)[1L]))
  str(unclass(x), give.attr = FALSE)
  invisible(x)
}

.fmri_lm_robust_enabled <- function(x) {
  type <- if (is.list(x)) x$type else x
  !is.null(type) && !identical(type, FALSE) && !identical(type, "none") &&
    !identical(type, "FALSE")
}
