.fmri_lm_noise_from_legacy <- function(ar_options) {
  if (inherits(ar_options, "fmri_lm_noise_spec")) return(ar_options)
  ar_options <- ar_options %||% list()
  .fmri_lm_check_keys(
    ar_options,
    c("struct", "p", "q", "iter_gls", "global", "voxelwise", "exact_first",
      "shared_estimator",
      "censor", "by_cluster", "order", "shrink_c0", "cor_struct", "iter"),
    "ar_options"
  )
  struct <- ar_options$struct %||% ar_options$cor_struct
  order <- ar_options$order
  if (is.null(struct) && !is.null(order)) {
    order <- .fmri_lm_number(order, "ar_options$order", lower = 0, integer = TRUE)
    struct <- if (order == 0L) "iid" else if (order == 1L) "ar1" else if (order == 2L) "ar2" else "arp"
  }
  struct <- struct %||% "iid"
  if (identical(struct, "none")) struct <- "iid"
  p <- ar_options$p
  if (is.null(p) && identical(struct, "arp")) p <- order
  global <- ar_options$global %||% FALSE
  by_cluster <- ar_options$by_cluster %||% FALSE
  global <- .fmri_lm_flag(global, "ar_options$global")
  by_cluster <- .fmri_lm_flag(by_cluster, "ar_options$by_cluster")
  if (global && by_cluster) {
    stop("`ar_options$global` and `ar_options$by_cluster` cannot both be TRUE.",
         call. = FALSE)
  }
  noise_spec(
    struct = struct,
    p = p,
    q = ar_options$q %||% 0L,
    iter_gls = ar_options$iter_gls %||% ar_options$iter %||% 1L,
    pooling = if (by_cluster) "parcel" else if (global) "global" else "run",
    shared_estimator = ar_options$shared_estimator %||% "pooled_acvf",
    voxelwise = ar_options$voxelwise %||% FALSE,
    exact_first = ar_options$exact_first %||% FALSE,
    censor = ar_options$censor,
    shrink_c0 = ar_options$shrink_c0 %||% 100L
  )
}

.fmri_lm_robust_from_legacy <- function(robust_options) {
  if (inherits(robust_options, "fmri_lm_robust_spec")) return(robust_options)
  robust_options <- robust_options %||% list()
  .fmri_lm_check_keys(
    robust_options,
    c("type", "k_huber", "c_tukey", "max_iter", "scale_scope", "reestimate_phi"),
    "robust_options"
  )
  type <- robust_options$type %||% "none"
  if (identical(type, FALSE) || identical(type, "FALSE")) type <- "none"
  if (identical(type, TRUE)) type <- "huber"
  scale_scope <- robust_options$scale_scope %||% "run"
  if (identical(scale_scope, "local")) scale_scope <- "voxel"
  robust_spec(
    type = type,
    k_huber = robust_options$k_huber %||% 1.345,
    c_tukey = robust_options$c_tukey %||% 4.685,
    max_iter = robust_options$max_iter %||% 2L,
    scale_scope = scale_scope,
    reestimate_phi = robust_options$reestimate_phi %||% FALSE
  )
}

.fmri_lm_weights_from_legacy <- function(options) {
  if (inherits(options, "fmri_lm_weights_spec")) return(options)
  options <- options %||% list()
  .fmri_lm_check_keys(options, c("enabled", "method", "threshold", "weights"),
                      "volume_weights_options")
  enabled <- options$enabled %||% FALSE
  enabled <- .fmri_lm_flag(enabled, "volume_weights_options$enabled")
  method <- if (enabled) options$method %||% "inverse_squared" else "none"
  weights_spec(method = method, threshold = options$threshold %||% 1.5,
               values = options$weights)
}

.fmri_lm_projection_from_legacy <- function(options) {
  if (inherits(options, "fmri_lm_projection_spec")) return(options)
  options <- options %||% list()
  .fmri_lm_check_keys(
    options,
    c("enabled", "nuisance_mask", "nuisance_matrix", "lambda", "warn_redundant"),
    "soft_subspace_options"
  )
  enabled <- options$enabled %||% FALSE
  enabled <- .fmri_lm_flag(enabled, "soft_subspace_options$enabled")
  projection_spec(
    method = if (enabled) "soft_subspace" else "none",
    nuisance_mask = options$nuisance_mask,
    nuisance_matrix = options$nuisance_matrix,
    lambda = options$lambda %||% "auto",
    warn_redundant = options$warn_redundant %||% TRUE
  )
}

#' Configuration for fmri_lm fitting
#'
#' `fmri_lm_control()` creates an `fmri_lm_control` object collecting all
#' options for robust and autoregressive modelling. It validates inputs and
#' applies defaults so downstream functions receive a single structured list.
#'
#' For common use cases, `fmri_lm()` provides convenience parameters that
#' are easier to use than these detailed option lists:
#' \itemize{
#'   \item \code{volume_weights = TRUE} enables volume weighting with defaults
#'   \item \code{volume_weights = "tukey"} enables with Tukey method
#'   \item \code{nuisance_projection = N} enables soft projection with matrix N
#'   \item \code{nuisance_projection = "mask.nii"} enables with mask file
#' }
#'
#' Use the `*_options` lists below only when you need fine-grained control.
#'
#' @param robust_options list of robust fitting options. See Details.
#' @param ar_options list of autoregressive modelling options. See Details.
#' @param volume_weights_options list of volume weighting options. See Details.
#'   For simple cases, use the \code{volume_weights} parameter in \code{fmri_lm()} instead.
#' @param soft_subspace_options list of soft subspace projection options. See Details.
#'   For simple cases, use the \code{nuisance_projection} parameter in \code{fmri_lm()} instead.
#' @return An object of class `fmri_lm_control`.
#' @details
#' `robust_options` may contain:
#'   * `type` (`FALSE`, "huber", "bisquare")
#'   * `k_huber`
#'   * `c_tukey`
#'   * `max_iter`
#'   * `scale_scope` ("run", "global", "voxel")
#'   * `reestimate_phi` (logical)
#'
#' `ar_options` may contain:
#'   * `struct` ("iid", "ar1", "ar2", "arp")
#'   * `p` (order for "arp")
#'   * `q` (moving-average order; positive values request ARMA inference)
#'   * `iter_gls` (maximum GLS iterations; pure design-corrected AR uses one)
#'   * `global` (logical, use global phi)
#'   * `shared_estimator` ("pooled_acvf" or experimental "mean_series")
#'   * `voxelwise` (logical)
#'   * `exact_first` (logical)
#'   * `censor` (integer vector of timepoints to exclude from AR estimation,
#'     logical vector where TRUE = censored, or "auto" to extract from dataset)
#'
#' `volume_weights_options` may contain:
#'   * `enabled` (logical, whether to compute and apply volume weights)
#'   * `method` ("inverse_squared", "soft_threshold", "tukey")
#'   * `threshold` (numeric, DVARS threshold for weighting)
#'   * `weights` (optional pre-computed weight vector)
#'
#' `soft_subspace_options` may contain:
#'   * `enabled` (logical, whether to apply soft subspace projection)
#'   * `nuisance_mask` (path to NIfTI mask or logical vector)
#'   * `nuisance_matrix` (pre-computed nuisance timeseries matrix)
#'   * `lambda` (numeric, "auto", or "gcv")
#'   * `warn_redundant` (logical, warn if baseline has nuisance terms)
#'
#' This list controls soft subspace preprocessing. It is separate from the
#' built-in fast `fmri_lm()` engines. To use reduced-rank-regression GLS with
#' conditional or block-bootstrap standard errors, call `fmri_lm()` with
#' `engine = "rrr_gls"` and pass reduced-rank options in `engine_args`. To use
#' the sketched GLM path, call `fmri_lm()` with `engine = "latent_sketch"` and
#' pass a `lowrank_control()` object via `lowrank`.
#'
#' When `fmri_lm()` is called with the convenience argument
#' `nuisance_projection`, `enabled` is set automatically. When constructing a
#' `soft_subspace_options` list directly, set `enabled = TRUE` yourself.
#'
#' @noRd
.fmri_lm_control_legacy <- function(estimation = NULL,
                            noise = NULL,
                            robust = NULL,
                            variance = NULL,
                            weights = NULL,
                            projection = NULL,
                            robust_options = NULL,
                            ar_options = NULL,
                            volume_weights_options = NULL,
                            soft_subspace_options = NULL,
                            na_action = c("error", "propagate")) {
  na_action <- match.arg(na_action)

  conflicts <- c(
    noise = !is.null(noise) && !is.null(ar_options),
    robust = !is.null(robust) && !is.null(robust_options),
    weights = !is.null(weights) && !is.null(volume_weights_options),
    projection = !is.null(projection) && !is.null(soft_subspace_options)
  )
  if (any(conflicts)) {
    stop("Canonical specs cannot be combined with their legacy option lists: ",
         paste(names(conflicts)[conflicts], collapse = ", "), ".", call. = FALSE)
  }

  estimation <- estimation %||% estimation_spec()
  noise <- noise %||% .fmri_lm_noise_from_legacy(ar_options)
  robust <- robust %||% .fmri_lm_robust_from_legacy(robust_options)
  variance <- variance %||% variance_spec()
  weights <- weights %||% .fmri_lm_weights_from_legacy(volume_weights_options)
  projection <- projection %||% .fmri_lm_projection_from_legacy(soft_subspace_options)

  expected <- list(
    estimation = "fmri_lm_estimation_spec", noise = "fmri_lm_noise_spec",
    robust = "fmri_lm_robust_spec", variance = "fmri_lm_variance_spec",
    weights = "fmri_lm_weights_spec", projection = "fmri_lm_projection_spec"
  )
  values <- list(estimation = estimation, noise = noise, robust = robust,
                 variance = variance, weights = weights, projection = projection)
  for (nm in names(expected)) {
    if (!inherits(values[[nm]], expected[[nm]])) {
      stop(sprintf("`%s` must be created by `%s_spec()`.", nm,
                   sub("fmri_lm_", "", sub("_spec$", "", expected[[nm]]))),
           call. = FALSE)
    }
  }

  if (noise$q > 0L && .fmri_lm_robust_enabled(robust)) {
    stop("MA terms cannot yet be combined with robust fitting.", call. = FALSE)
  }
  if (noise$q > 0L && !identical(estimation$scope, "runwise_meta")) {
    stop("MA terms currently require `estimation_spec('runwise_meta')`.",
         call. = FALSE)
  }

  structure(
    c(values, list(
      # Temporary compatibility aliases used by the existing fitting kernels.
      ar = noise,
      volume_weights = weights,
      soft_subspace = projection,
      na_action = na_action,
      schema_version = 2L
    )),
    class = "fmri_lm_control"
  )
}

#' Configure an fMRI linear model
#'
#' `fmri_lm_control()` is the single statistical configuration boundary for
#' [fmri_lm()]. Each section is a validated typed specification; mechanical
#' choices such as chunking and parallelism belong in [compute_spec()].
#'
#' @param estimation An [estimation_spec()].
#' @param noise A [noise_spec()].
#' @param robust A [robust_spec()].
#' @param variance A [variance_spec()].
#' @param weights A [weights_spec()].
#' @param projection A [projection_spec()].
#' @param na_action Either `"error"` or `"propagate"`.
#' @param ... Transitional legacy option lists retained for one compatibility
#'   window. New code should use the typed sections above.
#' @return An object of class `fmri_lm_control`.
#' @export
fmri_lm_control <- function(estimation = NULL,
                            noise = NULL,
                            robust = NULL,
                            variance = NULL,
                            weights = NULL,
                            projection = NULL,
                            na_action = c("error", "propagate"), ...) {
  legacy_dots <- list(...)
  if (length(legacy_dots) &&
      (is.null(names(legacy_dots)) || any(!nzchar(names(legacy_dots))))) {
    stop("Every legacy control argument must be named.", call. = FALSE)
  }
  do.call(
    .fmri_lm_control_legacy,
    c(list(estimation = estimation, noise = noise, robust = robust,
           variance = variance, weights = weights, projection = projection,
           na_action = na_action), legacy_dots)
  )
}

# Validate combinations whose support depends on the built-in fitting path.
# Keep these checks out of fmri_lm_control(): registered engines receive the
# same lossless typed specification and may support a broader combination.
.validate_builtin_voxelwise_noise <- function(cfg) {
  noise <- cfg$noise %||% cfg$ar
  if (!isTRUE(noise$voxelwise)) {
    return(invisible(cfg))
  }

  if (identical(noise$struct, "iid")) {
    stop("Built-in voxelwise temporal covariance requires an AR structure.",
         call. = FALSE)
  }
  if (!identical(cfg$estimation$scope, "runwise_meta")) {
    stop(
      "Built-in voxelwise temporal covariance currently requires ",
      "`estimation_spec('runwise_meta')`; joint fitting estimates one shared ",
      "AR model.",
      call. = FALSE
    )
  }
  if (!identical(noise$pooling, "run")) {
    stop("Built-in voxelwise temporal covariance currently requires `pooling = 'run'`.",
         call. = FALSE)
  }
  if (!is.null(noise$censor)) {
    stop("Built-in voxelwise temporal covariance does not yet support censoring.",
         call. = FALSE)
  }
  if (!identical(noise$iter_gls, 1L)) {
    stop("Built-in voxelwise temporal covariance currently requires `iter_gls = 1`.",
         call. = FALSE)
  }
  if (isTRUE(cfg$robust$reestimate_phi)) {
    stop(
      "Built-in voxelwise temporal covariance does not yet support robust ",
      "AR re-estimation (`reestimate_phi = TRUE`).",
      call. = FALSE
    )
  }
  if (isTRUE(cfg$weights$enabled) || isTRUE(cfg$projection$enabled)) {
    stop(
      "Built-in voxelwise temporal covariance does not yet support volume ",
      "weighting or soft-subspace projection.",
      call. = FALSE
    )
  }

  invisible(cfg)
}

#' @export
print.fmri_lm_control <- function(x, ...) {
  cat("<fmri_lm_control>\n")
  str(list(estimation = x$estimation, noise = x$noise,
           robust = x$robust, variance = x$variance,
           weights = x$weights, projection = x$projection,
           na_action = x$na_action),
      give.attr = FALSE)
  invisible(x)
}
