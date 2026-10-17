#' Fit a Linear Regression Model for fMRI Data Analysis
#'
#' `fmri_lm` is a generic for fitting fMRI regression models. The
#' default interface accepts a model formula and dataset. An
#' alternative method can be used with a preconstructed
#' \code{fmri_model} object that already contains the design and data.
#'
#' @param formula A model formula describing the event structure or an
#'   \code{fmri_model} object.
#' @param ... Additional method arguments. Recognized engine arguments include
#'   `engine`, `engine_args`, and `lowrank`; see Details.
#' @return An object of class \code{fmri_lm}.
#' @export
fmri_lm <- function(formula, ...) {
  UseMethod("fmri_lm")
}

#' @keywords internal
#' @noRd
.fmri_lm_extract_engine_context <- function(dots) {
  engine <- dots$engine
  lowrank <- dots$lowrank
  engine_cfg <- dots$cfg
  engine_args <- dots$engine_args

  dots$engine <- NULL
  dots$lowrank <- NULL
  dots$cfg <- NULL
  dots$engine_args <- NULL

  if (!is.null(engine_args) && !is.list(engine_args)) {
    engine_args <- list(engine_args)
  }

  if (!is.null(engine) && !is.null(dots[[engine]])) {
    add_args <- dots[[engine]]
    if (!is.null(add_args) && !is.list(add_args)) {
      add_args <- list(add_args)
    }
    engine_args <- if (is.null(engine_args)) add_args else utils::modifyList(engine_args, add_args)
    dots[[engine]] <- NULL
  }

  if (is.null(engine_args)) {
    engine_args <- list()
  }

  if (length(dots) > 0L) {
    dot_names <- names(dots)
    if (is.null(dot_names)) {
      dot_names <- rep("<unnamed>", length(dots))
    }
    dot_names[dot_names == ""] <- "<unnamed>"
    stop("Unexpected arguments: ", paste(dot_names, collapse = ", "), call. = FALSE)
  }

  list(
    engine = engine,
    lowrank = lowrank,
    engine_cfg = engine_cfg,
    engine_args = engine_args
  )
}

#' @keywords internal
#' @noRd
.fmri_lm_normalize_robust_options <- function(robust,
                                              robust_options = NULL,
                                              robust_psi = NULL,
                                              robust_max_iter = NULL,
                                              robust_scale_scope = NULL) {
  robust_options <- robust_options %||% list()

  # `type` can be requested three ways. Previously the first non-NULL of
  # (robust, robust_psi, robust_options$type) won by position, which made
  # `robust_psi` unreachable (robust defaults to NULL only since this fix; it
  # used to default to FALSE, so its branch always fired first) and let an
  # explicit `robust = FALSE` be overridden by `robust_options$type`. Requests
  # are now collected and a genuine disagreement is an error.
  canon <- function(x, name) {
    if (length(x) != 1L || is.na(x)) {
      stop(sprintf("`%s` must be a single non-missing robust setting.", name),
           call. = FALSE)
    }
    if (is.logical(x)) {
      return(if (isTRUE(x)) "huber" else "none")
    }
    if (!is.character(x)) {
      stop(sprintf("`%s` must be TRUE, FALSE, 'huber', or 'bisquare'.", name),
           call. = FALSE)
    }
    value <- tolower(x)
    if (!value %in% c("none", "false", "huber", "bisquare")) {
      stop(sprintf("`%s` must be TRUE, FALSE, 'huber', or 'bisquare'.", name),
           call. = FALSE)
    }
    if (value == "false") "none" else value
  }
  requested <- c(
    if (!is.null(robust))                  stats::setNames(canon(robust, "robust"), "robust"),
    if (!is.null(robust_psi))              stats::setNames(canon(robust_psi, "robust_psi"), "robust_psi"),
    if ("type" %in% names(robust_options)) stats::setNames(canon(robust_options$type, "robust_options$type"), "robust_options$type")
  )
  if (length(unique(requested)) > 1L) {
    stop("Conflicting robust settings: ",
         paste(sprintf("`%s` = %s", names(requested), sQuote(requested)), collapse = ", "),
         ". Supply only one.", call. = FALSE)
  }
  if (length(requested)) {
    robust_options$type <- if (requested[[1L]] == "none") FALSE else requested[[1L]]
  } else {
    robust_options$type <- FALSE
  }

  if (!is.null(robust_max_iter) && !("max_iter" %in% names(robust_options))) {
    robust_options$max_iter <- robust_max_iter
  }
  if (!is.null(robust_scale_scope) && !("scale_scope" %in% names(robust_options))) {
    robust_options$scale_scope <- robust_scale_scope
  }

  robust_options
}

#' @keywords internal
#' @noRd
.fmri_lm_normalize_ar_options <- function(ar_options = NULL,
                                          ar_voxelwise = NULL,
                                          cor_struct = NULL,
                                          cor_iter = NULL,
                                          cor_global = NULL,
                                          ar1_exact_first = NULL,
                                          ar_p = NULL) {
  ar_options <- ar_options %||% list()

  # Each shorthand is an alias for one `ar_options` key. Previously the list won
  # by position and the shorthand was dropped without a word -- which also made
  # the documented precedence for `ar_voxelwise` ("overrides ar_options$voxelwise")
  # exactly backwards. A shorthand that agrees with the list is accepted; one that
  # disagrees is an error rather than a silent choice.
  aliases <- list(
    cor_struct      = list(key = "struct",      value = cor_struct),
    cor_iter        = list(key = "iter_gls",    value = cor_iter),
    cor_global      = list(key = "global",      value = cor_global),
    ar1_exact_first = list(key = "exact_first", value = ar1_exact_first),
    ar_p            = list(key = "p",           value = ar_p),
    ar_voxelwise    = list(key = "voxelwise",   value = ar_voxelwise)
  )
  for (nm in names(aliases)) {
    a <- aliases[[nm]]
    if (is.null(a$value)) next
    if (a$key %in% names(ar_options)) {
      if (!isTRUE(all.equal(ar_options[[a$key]], a$value))) {
        stop(sprintf(
          "Conflicting AR settings: `%s` = %s and `ar_options$%s` = %s. Supply only one.",
          nm, sQuote(format(a$value)), a$key, sQuote(format(ar_options[[a$key]]))),
          call. = FALSE)
      }
    } else {
      ar_options[[a$key]] <- a$value
    }
  }
  if (!("voxelwise" %in% names(ar_options))) {
    ar_options$voxelwise <- FALSE
  }

  if (!is.null(ar_options$order) && is.null(ar_options$struct)) {
    ar_order_tmp <- as.integer(ar_options$order[1])
    if (!is.finite(ar_order_tmp) || ar_order_tmp <= 0L) {
      ar_options$struct <- "iid"
    } else if (ar_order_tmp == 1L) {
      ar_options$struct <- "ar1"
    } else if (ar_order_tmp == 2L) {
      ar_options$struct <- "ar2"
    } else {
      ar_options$struct <- "arp"
      ar_options$p <- ar_options$p %||% ar_order_tmp
    }
  }
  ar_options$order <- NULL

  ar_options
}

#' @keywords internal
#' @noRd
.fmri_lm_normalize_preprocessing_options <- function(volume_weights_options = NULL,
                                                     soft_subspace_options = NULL,
                                                     volume_weights = NULL,
                                                     nuisance_projection = NULL) {
  if (!is.null(volume_weights) && !identical(volume_weights, FALSE)) {
    if (is.null(volume_weights_options)) {
      volume_weights_options <- list()
    }
    volume_weights_options$enabled <- TRUE
    if (is.character(volume_weights)) {
      volume_weights_options$method <- volume_weights
    }
  }

  if (!is.null(nuisance_projection)) {
    if (is.null(soft_subspace_options)) {
      soft_subspace_options <- list()
    }
    soft_subspace_options$enabled <- TRUE
    if (is.matrix(nuisance_projection)) {
      soft_subspace_options$nuisance_matrix <- nuisance_projection
    } else if (is.character(nuisance_projection)) {
      soft_subspace_options$nuisance_mask <- nuisance_projection
    }
  }

  list(
    volume_weights_options = volume_weights_options,
    soft_subspace_options = soft_subspace_options
  )
}

#' @keywords internal
#' @noRd
.fmri_lm_build_config <- function(robust,
                                  control = NULL,
                                  robust_options = NULL,
                                  robust_psi = NULL,
                                  robust_max_iter = NULL,
                                  robust_scale_scope = NULL,
                                  ar_options = NULL,
                                  ar_voxelwise = NULL,
                                  cor_struct = NULL,
                                  cor_iter = NULL,
                                  cor_global = NULL,
                                  ar1_exact_first = NULL,
                                  ar_p = NULL,
                                  volume_weights_options = NULL,
                                  soft_subspace_options = NULL,
                                  volume_weights = NULL,
                                  nuisance_projection = NULL,
                                  engine_cfg = NULL) {
  supplied <- c(
    robust = !is.null(robust), robust_options = !is.null(robust_options),
    robust_psi = !is.null(robust_psi), robust_max_iter = !is.null(robust_max_iter),
    robust_scale_scope = !is.null(robust_scale_scope),
    ar_options = !is.null(ar_options), ar_voxelwise = !is.null(ar_voxelwise),
    cor_struct = !is.null(cor_struct), cor_iter = !is.null(cor_iter),
    cor_global = !is.null(cor_global), ar1_exact_first = !is.null(ar1_exact_first),
    ar_p = !is.null(ar_p),
    volume_weights_options = !is.null(volume_weights_options),
    soft_subspace_options = !is.null(soft_subspace_options),
    volume_weights = !is.null(volume_weights),
    nuisance_projection = !is.null(nuisance_projection)
  )
  also_supplied <- paste0("`", names(supplied)[supplied], "`", collapse = ", ")
  if (!is.null(control)) {
    if (!inherits(control, "fmri_lm_control")) {
      stop("`control` must be an object created by `fmri_lm_control()`.", call. = FALSE)
    }
    if (!is.null(engine_cfg)) {
      stop("`control` and deprecated `cfg` cannot both be supplied.", call. = FALSE)
    }
    if (any(supplied)) {
      stop("`control` cannot be combined with legacy configuration arguments. Also supplied: ",
           also_supplied, ".", call. = FALSE)
    }
    return(control)
  }
  # A supplied `cfg` used to replace the config wholesale, silently discarding
  # every other configuration argument -- so `cor_struct = "ar2"` alongside a cfg
  # naming "iid" ran iid without a word. Refuse the ambiguity instead.
  if (!is.null(engine_cfg)) {
    if (!inherits(engine_cfg, "fmri_lm_control")) {
      stop("`cfg` must be an object created by `fmri_lm_control()`.", call. = FALSE)
    }
    if (any(supplied)) {
      stop("`cfg` cannot be combined with other configuration arguments. ",
           "Also supplied: ", also_supplied,
           ". Put every setting in `cfg`, or drop `cfg` and use the individual arguments.",
           call. = FALSE)
    }
    return(engine_cfg)
  }

  robust_options <- .fmri_lm_normalize_robust_options(
    robust = robust,
    robust_options = robust_options,
    robust_psi = robust_psi,
    robust_max_iter = robust_max_iter,
    robust_scale_scope = robust_scale_scope
  )

  ar_options <- .fmri_lm_normalize_ar_options(
    ar_options = ar_options,
    ar_voxelwise = ar_voxelwise,
    cor_struct = cor_struct,
    cor_iter = cor_iter,
    cor_global = cor_global,
    ar1_exact_first = ar1_exact_first,
    ar_p = ar_p
  )

  preprocessing <- .fmri_lm_normalize_preprocessing_options(
    volume_weights_options = volume_weights_options,
    soft_subspace_options = soft_subspace_options,
    volume_weights = volume_weights,
    nuisance_projection = nuisance_projection
  )

  .fmri_lm_control_legacy(
    robust_options = robust_options,
    ar_options = ar_options,
    volume_weights_options = preprocessing$volume_weights_options,
    soft_subspace_options = preprocessing$soft_subspace_options
  )
}

#' Resolve the statistical-scope/execution boundary
#'
#' With a canonical `compute`, the estimation scope comes from `cfg`. The
#' legacy execution arguments instead set the scope (`strategy`) and are
#' translated to a `compute_spec()`.
#' @return List with the `cfg` and `compute` the fit runs with.
#' @keywords internal
#' @noRd
.fmri_lm_resolve_execution <- function(cfg, compute = NULL,
                                       strategy, nchunks, use_fast_path,
                                       progress, parallel_voxels,
                                       parallel_chunks,
                                       legacy_supplied = logical()) {
  if (!is.null(compute)) {
    if (!inherits(compute, "fmri_lm_compute_spec")) {
      stop("`compute` must be an object created by `compute_spec()`.", call. = FALSE)
    }
    if (any(legacy_supplied)) {
      stop("`compute` cannot be combined with legacy execution arguments: ",
           paste0("`", names(legacy_supplied)[legacy_supplied], "`", collapse = ", "),
           ".", call. = FALSE)
    }
    return(list(cfg = cfg, compute = compute))
  }

  strategy <- match.arg(strategy, c("runwise", "chunkwise"))
  cfg$estimation <- estimation_spec(
    if (identical(strategy, "runwise")) "runwise_meta" else "joint"
  )
  compute <- compute_spec(
    voxel_chunks = as.integer(nchunks),
    backend = if (isTRUE(use_fast_path)) "matrix" else "reference",
    parallel = if (isTRUE(parallel_voxels)) "voxels" else if (isTRUE(parallel_chunks)) "chunks" else "none",
    progress = progress
  )
  list(cfg = cfg, compute = compute)
}

#' @keywords internal
#' @noRd
.fmri_lm_attach_config_metadata <- function(fit, requested_cfg, executed_cfg = requested_cfg) {
  stopifnot(inherits(requested_cfg, "fmri_lm_control"))
  stopifnot(inherits(executed_cfg, "fmri_lm_control"))

  attr(fit, "requested_config") <- requested_cfg
  attr(fit, "executed_config") <- executed_cfg
  attr(fit, "config") <- executed_cfg
  attr(fit, "requested_control") <- requested_cfg
  attr(fit, "executed_control") <- executed_cfg
  fit
}

#' @keywords internal
#' @noRd
.fmri_lm_normalize_engine_name <- function(engine) {
  if (identical(engine, "sketch")) {
    "latent_sketch"
  } else {
    engine
  }
}

#' @keywords internal
#' @noRd
.fmri_lm_resolve_engine_spec <- function(engine) {
  engine_name <- .fmri_lm_normalize_engine_name(engine)
  spec <- get_engine_spec(engine_name, include_functions = TRUE)
  if (is.null(spec)) {
    stop(sprintf("Unknown engine '%s'.", engine), call. = FALSE)
  }
  spec
}

#' Warn about execution arguments an engine silently ignores
#'
#' `.fmri_lm_dispatch_engine()` is called without `strategy`, `nchunks`,
#' `use_fast_path`, `progress`, or the `parallel_*` flags, so supplying them
#' alongside `engine` had no effect and said nothing. Detection is via
#' `match.call()` rather than `missing()` because these arguments have non-NULL
#' defaults and some are reassigned by `match.arg()` before dispatch.
#'
#' @keywords internal
#' @noRd
.fmri_lm_warn_engine_ignores <- function(call, engine) {
  ignored <- c("strategy", "nchunks", "use_fast_path", "progress",
               "parallel_voxels", "parallel_chunks")
  supplied <- intersect(ignored, names(call)[-1L])
  if (length(supplied)) {
    warning(sprintf(
      "engine = \"%s\" ignores %s; %s control the built-in fitter only.",
      engine, paste0("`", supplied, "`", collapse = ", "),
      if (length(supplied) == 1L) "it controls" else "they control"),
      call. = FALSE)
  }
  invisible(supplied)
}

#' @keywords internal
#' @noRd
.fmri_lm_dispatch_engine <- function(model, dataset, engine, lowrank = NULL, cfg, engine_args = list()) {
  if (!is.null(lowrank) && is.null(engine_args$lowrank)) {
    engine_args$lowrank <- lowrank
  }

  spec <- .fmri_lm_resolve_engine_spec(engine = engine)
  # Engines default to one joint fit. A control object built without an
  # explicit `estimation` carries fmri_lm()'s built-in runwise default, which
  # is not a request made of the engine (#247).
  if (!.fmri_lm_estimation_supplied(cfg)) {
    cfg$estimation <- estimation_spec("joint")
  }
  .validate_engine_capabilities(spec$name, cfg, spec$capabilities)
  .validate_engine_context(spec$name, model, dataset, engine_args, cfg, spec$capabilities)
  executed_cfg <- .derive_engine_execution_config(cfg, spec$capabilities)

  if (is.function(spec$preflight)) {
    spec$preflight(model, dataset, engine_args, executed_cfg)
  }

  fit <- spec$fit(model, dataset, engine_args, executed_cfg)
  if (!inherits(fit, "fmri_lm")) {
    stop(sprintf("Engine '%s' must return an object of class 'fmri_lm'", spec$name), call. = FALSE)
  }

  fit$result <- .fmri_lm_finalize_result(
    result = fit$result, cfg = executed_cfg,
    compute = attr(fit, "compute"), model = model, dataset = dataset,
    strategy = spec$strategy, engine = spec$name
  )

  fit <- .fmri_lm_attach_config_metadata(fit, requested_cfg = cfg, executed_cfg = executed_cfg)
  if (is.null(attr(fit, "strategy"))) {
    attr(fit, "strategy") <- spec$strategy
  }
  if (is.null(attr(fit, "engine"))) {
    attr(fit, "engine") <- spec$name
  }

  fit
}

#' @noRd
#' @param block The model formula for block structure.
#' @param baseline_model (Optional) A \code{baseline_model} object. Default is \code{NULL}.
#' @param dataset An \code{fmri_frame} object containing the time-series data.
#' @param durations A vector of event durations. Default is \code{0}.
#' @param drop_empty Logical. Whether to remove factor levels with zero size. Default is \code{TRUE}.
#' @param robust Logical or character. Either \code{FALSE} (no robust fitting), 
#'   \code{TRUE} (use Huber), or one of \code{"huber"} or \code{"bisquare"}. Default is \code{FALSE}.
#' @param robust_options List of robust fitting options. See Details.
#' @param ar_options List of autoregressive modeling options. See Details.
#' @param strategy The data splitting strategy, either \code{"runwise"} or \code{"chunkwise"}. Default is \code{"runwise"}.
#' @param nchunks Number of data chunks when strategy is \code{"chunkwise"}.
#'   This controls memory partitioning; chunks are processed sequentially unless
#'   \code{parallel_chunks = TRUE}. Default is \code{10}.
#' @param use_fast_path Logical. If \code{TRUE} (the default), use the fast
#'   matrix engine, which supports OLS, AR whitening, robust, and preprocessing.
#'   \code{FALSE} selects the formula/lm reference engine (no robust or full
#'   preprocessing support); it is retained mainly as a parity oracle.
#' @param progress Logical. Whether to display a progress bar during model fitting. Default is \code{FALSE}.
#' @param ar_voxelwise Logical. Estimate AR parameters voxel-wise. Shorthand for
#'   \code{ar_options$voxelwise}; supplying both with different values is an error.
#' @param parallel_voxels Logical. Parallelize across voxels where supported;
#'   this does not control chunkwise execution.
#' @param parallel_chunks Logical. For \code{strategy = "chunkwise"}, process
#'   chunks with \code{future.apply::future_lapply()} using the active
#'   \code{future} plan. Default is \code{FALSE}.
#' @param cor_struct Character. Shorthand for \code{ar_options$struct} (e.g., "ar1", "ar2", "arp").
#' @param cor_iter Integer. Shorthand for \code{ar_options$iter_gls}.
#' @param cor_global Logical. Shorthand for \code{ar_options$global}.
#' @param ar1_exact_first Logical. Shorthand for \code{ar_options$exact_first}.
#' @param ar_p Integer. Shorthand for \code{ar_options$p}.
#' @param robust_psi Character. Shorthand for \code{robust_options$type} (e.g., "huber", "bisquare").
#' @param robust_max_iter Integer. Shorthand for \code{robust_options$max_iter}.
#' @param robust_scale_scope Character. Shorthand for \code{robust_options$scale_scope} ("run", "global", or "voxel").
#' @param volume_weights_options List of volume weighting options. See \code{\link{fmri_lm_control}}.
#' @param soft_subspace_options List of soft subspace projection options. See \code{\link{fmri_lm_control}}.
#' @param volume_weights Logical or character. Simple interface for volume weighting:
#'   \itemize{
#'     \item \code{TRUE}: Enable with default method ("inverse_squared")
#'     \item \code{"inverse_squared"}, \code{"soft_threshold"}, \code{"tukey"}: Enable with specific method
#'     \item \code{FALSE} or \code{NULL}: Disable (default)
#'   }
#'   For fine-grained control, use \code{volume_weights_options} instead.
#' @param nuisance_projection Matrix, character path, or NULL. Simple interface for
#'   soft subspace projection:
#'   \itemize{
#'     \item Matrix: Use as nuisance timeseries directly
#'     \item Character: Path to NIfTI mask for WM/CSF voxels
#'     \item \code{NULL}: Disable (default)
#'   }
#'   For fine-grained control (lambda selection, warnings), use \code{soft_subspace_options}.
#' @return A fitted linear regression model for fMRI data analysis.
#' 
#' @details
#' \code{robust_options} may contain:
#' \itemize{
#'   \item \code{type}: Character or logical. Type of robust fitting (\code{FALSE}, \code{"huber"}, \code{"bisquare"})
#'   \item \code{k_huber}: Numeric. Tuning constant for Huber's psi (default: 1.345)
#'   \item \code{c_tukey}: Numeric. Tuning constant for Tukey's bisquare psi (default: 4.685)
#'   \item \code{max_iter}: Integer. Maximum IRLS iterations (default: 2)
#'   \item \code{scale_scope}: Character. Scope for scale estimation (\code{"run"} or \code{"global"})
#'   \item \code{reestimate_phi}: Logical. Whether to re-estimate AR parameters after robust fitting
#' }
#' 
#' \code{ar_options} may contain:
#' \itemize{
#'   \item \code{struct}: Character. Correlation structure (\code{"iid"}, \code{"ar1"}, \code{"ar2"}, \code{"arp"})
#'   \item \code{p}: Integer. AR order when \code{struct = "arp"}
#'   \item \code{q}: Nonnegative MA order. Positive values fit ARMA(p, q)
#'   through the built-in matrix backend and currently require runwise
#'   meta-estimation, run pooling, model-based variance, no censoring, and no
#'   robust or voxelwise fitting.
#'   \item \code{iter_gls}: Integer. Maximum GLS iterations (default: 1).
#'   Pure AR uses one initial OLS estimate; iterative refinement applies to ARMA.
#'   \item \code{global}: Logical. Use global AR coefficients (default: FALSE)
#'   \item \code{shared_estimator}: Spatial estimator for shared AR:
#'   \code{"pooled_acvf"} (default) or experimental \code{"mean_series"}
#'   \item \code{voxelwise}: Logical. Estimate AR parameters voxel-wise
#'   (default: FALSE). The built-in path currently requires runwise
#'   meta-estimation, run pooling, one GLS iteration, no censoring, and no
#'   volume weighting, soft-subspace projection, or robust AR re-estimation.
#'   \item \code{exact_first}: Logical. Apply exact AR(1) scaling to first sample (default: FALSE)
#' }
#'
#' With shared built-in AR estimation, the default `shared_estimator =
#' "pooled_acvf"` pools voxel residual autocovariances and targets a typical
#' voxel covariance. The experimental `shared_estimator = "mean_series"`
#' instead estimates from the cross-voxel mean OLS residual series and targets
#' the coherent spatial component, which can be substantially more
#' autocorrelated than an individual voxel. OLS projection commutes with this
#' averaging, so both shared estimators retain the matching design correction.
#' Design-corrected AR estimates are formed once from initial OLS residuals and
#' held fixed for the GLS solve; later GLS residuals do not share the OLS
#' residual operator.
#'
#' Built-in fast engines are selected through `...`:
#' \itemize{
#'   \item \code{engine = "latent_sketch"} uses the sketched GLM path. Configure
#'   it with \code{lowrank = lowrank_control(...)}. The alias
#'   \code{engine = "sketch"} is normalized to \code{"latent_sketch"}.
#'   \item \code{engine = "rrr_gls"} uses the reduced-rank-regression GLS path.
#'   Configure it with \code{engine_args = list(...)}. This engine supports
#'   shared AR whitening from \code{ar_options} and inference for event/task
#'   parameters only; nuisance or baseline terms may be estimated, but
#'   standard-error and contrast inference is restricted to event/task
#'   coefficients.
#' }
#'
#' For \code{engine = "rrr_gls"}, \code{engine_args} may contain:
#' \itemize{
#'   \item \code{rank}: Positive integer rank. If \code{NULL} and
#'   \code{rank_mode = "fixed"}, the engine uses the full available task rank.
#'   \item \code{rank_mode}: One of \code{"fixed"}, \code{"energy"}, or
#'   \code{"rss_budget"} (default: \code{"fixed"}).
#'   \item \code{energy_keep}: Fraction of task singular-value energy to retain
#'   when \code{rank_mode = "energy"} (default: \code{0.99}).
#'   \item \code{rss_budget}: Non-negative residual-sum-of-squares budget used
#'   when \code{rank_mode = "rss_budget"}.
#'   \item \code{se_mode}: Standard-error mode, either \code{"conditional"} or
#'   \code{"bootstrap"} (default: \code{"conditional"}). Bootstrap standard
#'   errors resample task-subspace residuals in blocks.
#'   \item \code{bootstrap_n}: Number of bootstrap resamples when
#'   \code{se_mode = "bootstrap"} (default: \code{200}).
#'   \item \code{bootstrap_block_size}: Positive integer block size for
#'   bootstrap resampling. If \code{NULL}, a block size of 1 is used.
#'   \item \code{bootstrap_seed}: Optional integer seed for reproducible
#'   bootstrap resampling.
#'   \item \code{contrast_policy}: How to handle post-hoc contrasts that include
#'   non-event coefficients: \code{"warn_drop"}, \code{"drop"}, or
#'   \code{"error"} (default: \code{"warn_drop"}).
#' }
#'
#' The aliases \code{energy} and \code{nboot} are accepted as legacy shorthand
#' for \code{energy_keep} and \code{bootstrap_n}, respectively. Prefer the
#' canonical names in new code.
#' 
#' @export
#' @seealso \code{\link{matrix_frame}}, \code{\link{fmri_lm_fit}}, \code{\link{fmri_lm_control}}
#' @examples
#' 
#' facedes <- subset(read.table(system.file("extdata", "face_design.txt", package = "fmrireg"), 
#' header=TRUE), face_gen != "n/a")
#' facedes$face_gen <- droplevels(factor(facedes$face_gen))
#' sframe <- sampling_frame(rep(430/2,6), TR=2)
#' ev <- event_model(onset ~ hrf(face_gen, basis="gaussian"), data=facedes, 
#' block= ~ run, sampling_frame=sframe)
#' globonsets <- fmrihrf::global_onsets(sframe, facedes$onset, facedes$run)
#' reg1_signal <- regressor(globonsets[facedes$face_gen == "male"], hrf=fmrihrf::HRF_GAUSSIAN)
#' reg2_signal <- regressor(globonsets[facedes$face_gen == "female"], hrf=fmrihrf::HRF_GAUSSIAN)
#' time <- samples(sframe, global=TRUE)
#' y1 <- fmrihrf::evaluate(reg1_signal, time)*1.5
#' y2 <- fmrihrf::evaluate(reg2_signal, time)*3.0
#' y <- y1+y2
#' ys1 <- y + rnorm(length(y), sd=.02)
#' ys2 <- y + rnorm(length(y), sd=.02)
#' 
#' h <<- gen_hrf(fmrihrf::hrf_bspline, N=7, span=25)
#' dset <- matrix_frame(cbind(ys1,ys2), TR=2, 
#'                        run_length=fmrihrf::blocklens(sframe), 
#'                        event_table=facedes)
#' flm <- fmri_lm(onset ~ hrf(face_gen, 
#'                            basis=gen_hrf(fmrihrf::hrf_bspline, N=7, span=25)), 
#'                block = ~ run, 
#'                strategy="chunkwise", nchunks=1, dataset=dset)
#'
#' \dontrun{
#' fit_rrr <- fmri_lm(
#'   onset ~ hrf(condition), block = ~ run, dataset = dset,
#'   ar_options = list(struct = "ar1"),
#'   engine = "rrr_gls",
#'   engine_args = list(
#'     rank_mode = "energy",
#'     energy_keep = 0.99,
#'     se_mode = "bootstrap",
#'     bootstrap_n = 200L,
#'     bootstrap_block_size = 24L,
#'     bootstrap_seed = 1L
#'   )
#' )
#' standard_error(fit_rrr, type = "estimates")
#' }
#' 
.fmri_lm_formula_legacy <- function(formula, block, baseline_model = NULL, dataset,
                                    durations = 0, drop_empty = TRUE, ...) {
  formula <- .fmrireg_inject_registered_bases(formula)
  assert_that(is.formula(formula), msg = "'formula' must be a formula")
  assert_that(is.formula(block), msg = "'block' must be a formula")
  assert_that(inherits(dataset, "fmri_frame"), msg = "'dataset' must be an 'fmri_frame'")
  assert_that(is.numeric(durations), msg = "'durations' must be numeric")
  assert_that(is.logical(drop_empty), msg = "'drop_empty' must be logical")

  model <- create_fmri_model(formula, block, baseline_model, dataset,
                             durations = durations, drop_empty = drop_empty)
  .fmri_lm_model_legacy(model, dataset = dataset, ...)
}

#' @keywords internal
#' @noRd
.fmri_lm_model_legacy <- function(formula, dataset = NULL,
                               control = NULL, compute = NULL,
                               robust = NULL, robust_options = NULL,
                               ar_options = NULL,
                               volume_weights_options = NULL, soft_subspace_options = NULL,
                               strategy = c("runwise", "chunkwise"), nchunks = 10,
                               use_fast_path = TRUE, progress = FALSE,
                               ar_voxelwise = NULL, parallel_voxels = FALSE,
                               cor_struct = NULL, cor_iter = NULL,
                               cor_global = NULL, ar1_exact_first = NULL,
                               ar_p = NULL,
                               robust_psi = NULL, robust_max_iter = NULL,
                               robust_scale_scope = NULL,
                               volume_weights = NULL, nuisance_projection = NULL,
                               parallel_chunks = FALSE,
                               ...) {
  legacy_execution_supplied <- c(
    strategy = !missing(strategy), nchunks = !missing(nchunks),
    use_fast_path = !missing(use_fast_path), progress = !missing(progress),
    parallel_voxels = !missing(parallel_voxels),
    parallel_chunks = !missing(parallel_chunks)
  )
  engine_ctx <- .fmri_lm_extract_engine_context(list(...))
  assert_that(inherits(formula, "fmri_model"))
  if (!is.null(robust)) {
    .fmri_lm_normalize_robust_options(robust = robust)
  }
  assert_that(
    is.logical(parallel_chunks) && length(parallel_chunks) == 1L && !is.na(parallel_chunks),
    msg = "'parallel_chunks' must be TRUE or FALSE"
  )

  dataset <- dataset %||% formula$dataset %||% attr(formula, "dataset")
  if (is.null(dataset)) {
    stop("No dataset found in 'formula' and none supplied.")
  }
  assert_that(inherits(dataset, "fmri_frame"))

  cfg <- .fmri_lm_build_config(
    robust = robust,
    control = control,
    robust_options = robust_options,
    robust_psi = robust_psi,
    robust_max_iter = robust_max_iter,
    robust_scale_scope = robust_scale_scope,
    ar_options = ar_options,
    ar_voxelwise = ar_voxelwise,
    cor_struct = cor_struct,
    cor_iter = cor_iter,
    cor_global = cor_global,
    ar1_exact_first = ar1_exact_first,
    ar_p = ar_p,
    volume_weights_options = volume_weights_options,
    soft_subspace_options = soft_subspace_options,
    volume_weights = volume_weights,
    nuisance_projection = nuisance_projection,
    engine_cfg = engine_ctx$engine_cfg
  )

  if (!is.null(control) && is.null(compute)) compute <- compute_spec()
  execution <- .fmri_lm_resolve_execution(
    cfg = cfg, compute = compute, strategy = strategy, nchunks = nchunks,
    use_fast_path = use_fast_path, progress = progress,
    parallel_voxels = parallel_voxels, parallel_chunks = parallel_chunks,
    legacy_supplied = legacy_execution_supplied
  )
  cfg <- execution$cfg
  compute <- execution$compute

  if (!is.null(engine_ctx$engine)) {
    # Without a control object the scope came from the legacy `strategy`
    # default, which is not a request made of the engine.
    if (is.null(control)) cfg$estimation <- estimation_spec("joint")
    .fmri_lm_warn_engine_ignores(match.call(), engine_ctx$engine)
    fit <- .fmri_lm_dispatch_engine(
      model = formula,
      dataset = dataset,
      engine = engine_ctx$engine,
      lowrank = engine_ctx$lowrank,
      cfg = cfg,
      engine_args = engine_ctx$engine_args
    )
    attr(fit, "compute") <- compute
    return(fit)
  }

  fmri_lm_fit(formula, dataset, cfg, compute = compute)
}

#' Fit an fMRI linear model from a formula
#'
#' The public fitting boundary separates statistical choices in [fmri_lm_control()]
#' from mechanical execution choices in [compute_spec()]. Optional engines use
#' the same control object and receive their private arguments through
#' `engine_args`.
#'
#' @param control A validated [fmri_lm_control()]. The default fits runs
#'   separately and pools them by fixed-effect meta-analysis
#'   (`estimation_spec("runwise_meta")`). A control object built without an
#'   `estimation` argument keeps this runwise default for the built-in fitter;
#'   pass `fmri_lm_control(estimation = estimation_spec("joint"))` to fit one
#'   joint model across runs. HAC, sandwich, and Satterthwaite inference
#'   (see [variance_spec()]) are implemented only for a joint fit; when such a
#'   `variance` is requested without an `estimation`, [fmri_lm_control()]
#'   selects the joint scope and says so in a message. When `engine` is
#'   supplied and `control` does not specify `estimation` (or is omitted), the
#'   engine receives `estimation_spec("joint")`.
#' @param compute A validated [compute_spec()].
#' @param block Formula describing run/block structure.
#' @param baseline_model Optional baseline/nuisance model.
#' @param dataset An `fmri_frame`. For an `fmri_model` method this may be
#'   omitted when the model already owns its dataset.
#' @param durations Event durations passed to model construction.
#' @param drop_empty Remove empty factor levels during model construction.
#' @param engine Optional registered engine name.
#' @param engine_args Named list passed only to `engine`.
#' @param lowrank Optional [lowrank_control()] for the latent-sketch engine.
#' @param ... Transitional legacy arguments retained for one compatibility
#'   window. New code must express statistical choices in `control` and
#'   execution choices in `compute`; unknown arguments are rejected.
#' @rdname fmri_lm
#' @export
fmri_lm.formula <- function(formula, block, baseline_model = NULL, dataset,
                            durations = 0, drop_empty = TRUE,
                            control = fmri_lm_control(
                              estimation = estimation_spec("runwise_meta")
                            ),
                            compute = compute_spec(),
                            engine = NULL, engine_args = list(),
                            lowrank = NULL, ...) {
  control_missing <- missing(control)
  compute_missing <- missing(compute)
  .fmri_lm_call_front(
    .fmri_lm_formula_legacy,
    list(formula = formula, block = block, baseline_model = baseline_model,
         dataset = dataset, durations = durations, drop_empty = drop_empty),
    control = control, compute = compute,
    control_missing = control_missing, compute_missing = compute_missing,
    engine = engine, engine_args = engine_args, lowrank = lowrank,
    legacy_dots = list(...)
  )
}

#' @rdname fmri_lm
#' @export
fmri_lm.fmri_model <- function(formula, dataset = NULL,
                               control = fmri_lm_control(
                                 estimation = estimation_spec("runwise_meta")
                               ),
                               compute = compute_spec(),
                               engine = NULL, engine_args = list(),
                               lowrank = NULL, ...) {
  control_missing <- missing(control)
  compute_missing <- missing(compute)
  .fmri_lm_call_front(
    .fmri_lm_model_legacy,
    list(formula = formula, dataset = dataset),
    control = control, compute = compute,
    control_missing = control_missing, compute_missing = compute_missing,
    engine = engine, engine_args = engine_args, lowrank = lowrank,
    legacy_dots = list(...)
  )
}

#' Route an fmri_lm() method call through the legacy-argument adapter
#'
#' Keeps the migration adapter isolated behind the canonical public formals:
#' it can be deleted as one unit after the deprecation window. Defaults for
#' `control`/`compute` are dropped when legacy arguments are supplied so the
#' two vocabularies do not collide.
#' @keywords internal
#' @noRd
.fmri_lm_call_front <- function(front, args, control, compute,
                                control_missing, compute_missing,
                                engine, engine_args, lowrank, legacy_dots) {
  if (!is.list(engine_args)) {
    stop("`engine_args` must be a named list.", call. = FALSE)
  }
  if (length(engine_args) && (is.null(names(engine_args)) || any(!nzchar(names(engine_args))))) {
    stop("Every element of `engine_args` must be named.", call. = FALSE)
  }
  if (control_missing && !is.null(engine) && !length(legacy_dots)) {
    control <- fmri_lm_control(estimation = estimation_spec("joint"))
  }
  call_args <- c(
    args,
    list(control = if (control_missing && length(legacy_dots)) NULL else control,
         compute = if (compute_missing && length(legacy_dots)) NULL else compute),
    legacy_dots,
    list(engine = engine, engine_args = engine_args, lowrank = lowrank)
  )
  do.call(front, call_args)
}


#' Warn when a dataset's censor column cannot influence the fit
#'
#' Censoring only ever feeds AR estimation and whitening; it never drops
#' volumes from the regression itself. Under an iid noise model it therefore
#' has no effect at all, which is silent and easy to mistake for censoring
#' having been applied.
#'
#' @keywords internal
#' @noRd
.warn_if_censor_inert <- function(dataset, cfg) {
  spec <- cfg$ar$censor
  if (is.character(spec) && identical(spec, "none")) {
    return(invisible(NULL))
  }
  has_censor <- !is.null(spec) && !identical(spec, "auto")
  if (!has_censor) {
    has_censor <- !is.null(tryCatch(.dset_censor(dataset), error = function(e) NULL))
  }
  if (!has_censor) {
    return(invisible(NULL))
  }
  if (identical(cfg$ar$struct %||% "iid", "iid")) {
    warning(
      "Censoring was supplied but has no effect under an iid noise model: ",
      "censored volumes feed AR estimation only, they are not dropped from ",
      "the regression. Set `ar_options = list(struct = \"ar1\")` (or another ",
      "AR structure) to use it, or `censor = \"none\"` to silence this.",
      call. = FALSE
    )
  }
  invisible(NULL)
}

#' Fit an fMRI Linear Regression Model with a Specified Fitting Strategy
#'
#' This function fits an fMRI linear regression model using the specified
#' \code{fmri_model} object and dataset. The estimation scope in \code{cfg}
#' selects the strategy (\code{"runwise_meta"} fits runs separately and pools
#' them; \code{"joint"} fits one chunkwise model) and \code{compute} the
#' mechanics. It is primarily an internal function used by \code{fmri_lm}.
#'
#' @param fmrimod An \code{fmri_model} object.
#' @param dataset An \code{fmri_frame} object containing the time-series data.
#' @param cfg An \code{fmri_lm_control} object containing all fitting options. See \code{\link{fmri_lm_control}}.
#' @param compute A \code{\link{compute_spec}}: voxel chunks, matrix or
#'   reference backend, parallelism and progress.
#' @return A fitted fMRI linear regression model with the specified fitting strategy.
#' @keywords internal
#' @seealso \code{\link{fmri_lm}}, \code{\link{fmri_model}}, \code{\link{matrix_frame}}
fmri_lm_fit <- function(fmrimod, dataset, cfg, compute = compute_spec()) {
  # Validate config before inspecting path-dependent combinations.
  assert_that(inherits(cfg, "fmri_lm_control"), msg = "'cfg' must be an 'fmri_lm_control' object")
  assert_that(inherits(compute, "fmri_lm_compute_spec"),
              msg = "'compute' must be an object created by `compute_spec()`")
  .validate_builtin_voxelwise_noise(cfg)

  if (!identical(cfg$estimation$scope, "joint") &&
      (!identical(cfg$variance$method, "model") ||
       identical(cfg$variance$df, "satterthwaite"))) {
    stop("HAC, sandwich, and Satterthwaite inference currently require ",
         "`estimation_spec(scope = 'joint')`; runwise meta-estimation must ",
         "combine run-level covariance explicitly.", call. = FALSE)
  }
  
  # Error checking
  assert_that(inherits(fmrimod, "fmri_model"), msg = "'fmrimod' must be an 'fmri_model' object")
  assert_that(inherits(dataset, "fmri_frame"), msg = "'dataset' must be an 'fmri_frame' object")
  .warn_if_censor_inert(dataset, cfg)

  strategy <- if (identical(cfg$estimation$scope, "runwise_meta")) "runwise" else "chunkwise"
  nchunks <- compute$voxel_chunks
  use_fast_path <- identical(compute$backend, "matrix")
  progress <- compute$progress
  parallel_voxels <- identical(compute$parallel, "voxels")
  parallel_chunks <- identical(compute$parallel, "chunks")
  if (strategy == "runwise" && parallel_chunks) {
    warning("'parallel_chunks' only applies when strategy = 'chunkwise'.", call. = FALSE)
  }

  # Check for redundant nuisance regressors when soft subspace projection is enabled
  if (cfg$soft_subspace$enabled && cfg$soft_subspace$warn_redundant) {
    .check_redundant_nuisance(fmrimod$baseline_model, cfg$soft_subspace)
  }

  contrast_prep <- prepare_fmri_lm_contrasts(fmrimod)
  processed_conlist <- contrast_prep$processed
  standard_path_conlist <- contrast_prep$standard

  # Global AR and global robust scale both start from per-run OLS residuals;
  # one pass over the runs serves either or both.
  phi_global <- NULL
  sigma_global <- NULL
  need_global_ar <- cfg$ar$global && cfg$ar$struct != "iid"
  need_global_scale <- .fmri_lm_robust_enabled(cfg$robust) &&
    cfg$robust$scale_scope == "global"
  if (need_global_ar || need_global_scale) {
    run_chunks <- .dset_run_chunks(dataset)
    resid_chunks <- vector("list", length(run_chunks))
    design_chunks <- vector("list", length(run_chunks))
    row_med_chunks <- vector("list", length(run_chunks))
    for (ri in seq_along(run_chunks)) {
      rch <- run_chunks[[ri]]
      X_run <- .run_design(fmrimod, rch$chunk_num)
      Y_run <- as.matrix(.chunk_data(rch))
      ols <- .fast_lm_matrix(X_run, Y_run, .fast_preproject(X_run), return_fitted = TRUE)
      resid_run <- Y_run - ols$fitted
      if (need_global_ar) {
        resid_chunks[[ri]] <- resid_run
        design_chunks[[ri]] <- X_run
      }
      if (need_global_scale) {
        row_med_chunks[[ri]] <- matrixStats::rowMedians(abs(resid_run))
      }
    }

    if (need_global_ar) {
      resid_mat <- do.call(rbind, resid_chunks)
      # Each run was fitted with its own coefficient vector. The exact combined
      # residual operator is therefore block diagonal; row-binding the designs
      # would instead describe one coefficient vector shared across runs.
      design_mat <- .block_diagonal_design(design_chunks)
      run_lengths <- vapply(resid_chunks, nrow, integer(1))
      run_indices <- split(
        seq_len(sum(run_lengths)),
        rep(seq_along(run_lengths), run_lengths)
      )
      phi_global <- .estimate_shared_ar_parameters(
        resid_mat, .target_ar_order(cfg$ar) %||% 0L, cfg$ar,
        run_indices = run_indices,
        design = design_mat
      )
      cfg$ar$iter_gls <- 1L
    }

    if (need_global_scale) {
      row_med_all <- unlist(row_med_chunks, use.names = FALSE)
      sigma_global <- 1.4826 * median(row_med_all)
      if (sigma_global <= .Machine$double.eps) sigma_global <- .Machine$double.eps
    }
  }

  result <- switch(strategy,
                   "runwise" = runwise_lm(dataset, fmrimod, standard_path_conlist, # Pass full objects
                                         cfg = cfg, use_fast_path = use_fast_path,
                                         progress = progress,
                                         phi_fixed = phi_global,
                                         sigma_fixed = sigma_global,
                                         parallel_voxels = parallel_voxels
                                         ),
                   "chunkwise" = {
                    if (.dset_is_latent(dataset)) {
                      if (isTRUE(parallel_chunks)) {
                        warning(
                          "'parallel_chunks' is not currently supported for latent frames; ",
                          "using sequential chunkwise fitting.",
                          call. = FALSE
                        )
                      }
                      # latent frames have no fast matrix engine (they route
                      # through multiresponse_arma); force the slow path
                      # regardless of the global default. Intentional.
                      chunkwise_lm_latent(dataset, fmrimod, standard_path_conlist, # Pass full objects
                                   nchunks, cfg, verbose = FALSE, use_fast_path = FALSE,
                                   progress = progress,
                                   phi_fixed = phi_global,
                                   sigma_fixed = sigma_global
                                   ) # Do not pass use_fast_path
                    } else {
                      chunkwise_lm(dataset, fmrimod, standard_path_conlist, # Pass full objects
                                   nchunks, cfg = cfg, use_fast_path = use_fast_path,
                                   progress = progress,
                                   parallel_chunks = parallel_chunks,
                                   phi_fixed = phi_global,
                                   sigma_fixed = sigma_global
                                   )
                    }
                  })

  result <- .fmri_lm_finalize_result(
    result = result, cfg = cfg, compute = compute,
    model = fmrimod, dataset = dataset, strategy = strategy, engine = "builtin"
  )
  
  ret <- list(
    result = result,
    model = fmrimod,
    strategy = strategy,
    bcons = processed_conlist,
    dataset = dataset,
    ar_coef = result$ar_coef,
    ma_coef = result$ma_coef
  )
  
  class(ret) <- "fmri_lm"
  
  ret <- .fmri_lm_attach_config_metadata(ret, requested_cfg = cfg, executed_cfg = cfg)
  attr(ret, "strategy") <- strategy
  attr(ret, "compute") <- compute
  
  return(ret)
}


#' Compute Fitted Hemodynamic Response Functions for an fmri_lm Object
#'
#' This method computes the fitted hemodynamic response functions (HRFs) for an \code{fmri_lm} object.
#'
#' @param x An \code{fmri_lm} object for which the fitted HRFs should be computed.
#' @param sample_at A numeric vector of time points at which the HRFs should be sampled. Default is \code{seq(0, 24, by = 1)}.
#' @param ... Additional arguments (currently unused).
#' @return A list where each element corresponds to an event term in the \code{fmri_lm} object. Each element contains:
#' \describe{
#'   \item{\code{pred}}{A matrix of predicted HRF values.}
#'   \item{\code{design}}{A tibble containing the design matrix for the HRFs.}
#' }
#' @export
fitted_hrf.fmri_lm <- function(x, sample_at = seq(0, 24, by = 1), ...) {
  # Error checking
  assert_that(inherits(x, "fmri_lm"), msg = "'x' must be an 'fmri_lm' object")
  assert_that(is.numeric(sample_at), msg = "'sample_at' must be numeric")
  assert_that(all(is.finite(sample_at)), msg = "'sample_at' must contain only finite values")
  
  eterms <- terms(x$model$event_model)
  
  # If no event terms, return empty list
  if (is.null(eterms) || length(eterms) == 0) {
    return(list())
  }

  # Feature terms are sampled series, not trial HRFs. Reconstructing a
  # canonical HRF shape from their coefficients is the wrong object.
  keep <- !vapply(eterms, function(t) inherits(t, "feature_term"), logical(1))
  eterms <- eterms[keep]
  if (length(eterms) == 0L) {
    return(list())
  }

  # Pull coefficients directly from model results (voxels x coefficients)
  betas <- x$result$betas$data[[1]]$estimate[[1]]
  if (!is.matrix(betas)) {
    betas <- as.matrix(betas)
  }

  # Column mapping for each event term
  col_indices <- attr(x$model$event_model$design_matrix, "col_indices")
  if (is.null(col_indices) || length(col_indices) == 0L) {
    stop("Event model design matrix is missing 'col_indices' metadata.", call. = FALSE)
  }

  term_names <- names(eterms)
  sample_n <- length(sample_at)
  nvox <- nrow(betas)
  
  # For each event term, compute its indices based on the term's structure
  pred <- lapply(seq_along(eterms), function(i) {
    eterm <- eterms[[i]]
    term_name <- term_names[[i]]
    if (is.null(term_name) || !nzchar(term_name)) {
      term_name <- paste0("term_", i)
    }

    ind <- col_indices[[term_name]]
    if (is.null(ind) || length(ind) == 0L) {
      stop(sprintf("No coefficient indices found for event term '%s'.", term_name), call. = FALSE)
    }
    if (any(!is.finite(ind))) {
      stop(sprintf("Coefficient indices for term '%s' contain non-finite values.", term_name), call. = FALSE)
    }
    ind <- as.integer(ind)
    if (any(ind < 1L) || any(ind > ncol(betas))) {
      stop(
        sprintf(
          "Coefficient indices for term '%s' are out of bounds (max index %d, available %d).",
          term_name, max(ind), ncol(betas)
        ),
        call. = FALSE
      )
    }
    
    # Get the HRF specification (stored as an attribute in fmridesign)
    hrf_spec <- attr(eterm, "hrfspec")
    if (is.null(hrf_spec) && !is.null(eterm$hrfspec)) {
      # Backward-compatibility: some versions may store as a list element
      hrf_spec <- eterm$hrfspec
    }
    if (is.null(hrf_spec) && !is.null(eterm$hrf) && is.function(eterm$hrf)) {
      hrf_spec <- list(hrf = eterm$hrf)
    }
    
    # Fallback HRF if spec is missing
    hrf <- if (!is.null(hrf_spec) && !is.null(hrf_spec$hrf)) hrf_spec$hrf else fmrihrf::HRF_SPMG1
    if (!is.function(hrf)) {
      stop(sprintf("HRF specification for term '%s' is not callable.", term_name), call. = FALSE)
    }
    
    # Get the conditions (cells) for this term
    excond <- cells(eterm, exclude_basis = TRUE)
    ncond <- nrow(excond)

    excond_expanded <- excond[rep(seq_len(ncond), each = sample_n), , drop = FALSE]
    design <- tibble::add_column(
      tibble::as_tibble(excond_expanded),
      time = rep(sample_at, ncond),
      .before = 1L
    )

    if (sample_n == 0L || ncond == 0L) {
      return(list(pred = matrix(numeric(0), nrow = 0L, ncol = nvox), design = design))
    }

    # Create the HRF basis matrix at sample points
    G <- as.matrix(hrf(sample_at))
    if (nrow(G) != sample_n) {
      stop(
        sprintf(
          "HRF basis for term '%s' returned %d rows for %d sample points.",
          term_name, nrow(G), sample_n
        ),
        call. = FALSE
      )
    }

    nbasis_term <- ncol(G)
    expected_cols <- ncond * nbasis_term
    if (length(ind) != expected_cols) {
      stop(
        sprintf(
          "Term '%s' expects %d coefficient columns (%d conditions x %d basis), but found %d.",
          term_name, expected_cols, ncond, nbasis_term, length(ind)
        ),
        call. = FALSE
      )
    }

    # Select term coefficients and apply basis per condition (avoids large block-diagonal allocation).
    B <- t(betas[, ind, drop = FALSE])
    yh_blocks <- vector("list", ncond)
    for (j in seq_len(ncond)) {
      idx <- ((j - 1L) * nbasis_term + 1L):(j * nbasis_term)
      yh_blocks[[j]] <- G %*% B[idx, , drop = FALSE]
    }
    yh <- do.call(rbind, yh_blocks)

    list(pred = as.matrix(yh), design = design)
  })
  
  # Set names from event terms
  names(pred) <- term_names
  
  return(pred)
}




pull_stat <- function(x, type, element) {
  if (type == "betas") {
    ret <- x$result$betas$data[[1]][[element]][[1]]
    
    # Check bounds and filter valid indices
    max_col <- ncol(ret)
    valid_event_indices <- x$result$event_indices[x$result$event_indices <= max_col]
    
    if (length(valid_event_indices) == 0) {
      warning("No valid event indices found in pull_stat. Using all available columns.")
      valid_event_indices <- 1:max_col
    }
    
    ret <- ret[, valid_event_indices, drop = FALSE]
    
    # Use the actual column names from the design matrix instead of conditions()
    # This avoids duplicate names when multiple terms have the same variables
    dm <- design_matrix(x$model)
    if (!is.null(dm) && ncol(dm) >= max(valid_event_indices)) {
      actual_colnames <- colnames(dm)[valid_event_indices]
      colnames(ret) <- actual_colnames
    } else {
      # Fallback: use conditions but make them unique
      condition_names <- conditions(x$model$event_model)[1:length(valid_event_indices)]
      colnames(ret) <- make.names(condition_names, unique = TRUE)
    }
    
    # Ensure tibble output for consistency with original behavior
    res <- suppressMessages(tibble::as_tibble(ret, .name_repair = "check_unique"))
  } else if (type %in% c("contrasts", "F")) {
    con_type <- if (type == "F") "Fcontrast" else "contrast"
    ret <- x$result$contrasts %>% dplyr::filter(type == con_type)
    if (nrow(ret) == 0) {
      stop(if (type == "F") "No F contrasts for this model." else "No simple contrasts for this model.")
    }
    cnames <- ret$name
    # Name inputs before binding to avoid temporary name-repair messages.
    out <- lapply(ret$data, function(x) x[[element]]) %>% stats::setNames(cnames) %>% dplyr::bind_cols()
    names(out) <- cnames
    out
  } else {
    stop("Invalid type specified. Must be 'betas', 'contrasts', or 'F'.")
  }
}




#' Fit Linear Model Contrasts
#'
#' This function computes contrasts and beta statistics for a fitted linear model.
#'
#' @param fit A fitted linear model object.
#' @param conlist A list of contrast matrices.
#' @param fcon A list of F-contrasts.
#' @param vnames Variable names corresponding to the model coefficients.
#' @param se Logical. Whether to compute standard errors. Default is \code{TRUE}.
#' @return A list containing contrasts, beta statistics, and the fitted model.
#' @keywords internal
fit_lm_contrasts <- function(fit, conlist, fcon, vnames, se = TRUE) {
  # One (X'X)^-1 for the model, shared by every contrast and the beta table.
  cov.unscaled <- .lm_cov_unscaled(fit, varnames = vnames)
  conres <- if (!is.null(conlist)) {
    ret <- lapply(conlist, function(con) {
      # Extract colind from the contrast object's attributes
      colind <- attr(con, "colind")
      if (is.null(colind)) {
        stop(sprintf(
          "Missing colind attribute for requested contrast '%s'.",
          con$name %||% "unnamed"
        ), call. = FALSE)
      }
      estimate_contrast(con, fit, colind, cov.unscaled = cov.unscaled)
    })
    names(ret) <- sapply(conlist, function(x) x$name %||% "unnamed")
    ret
  } else {
    list()
  }
  
  bstats <- beta_stats(fit, vnames, se = se, cov.unscaled = cov.unscaled)
  list(contrasts = conres, bstats = bstats, fit = fit)
}




#' Fit Multiresponse Linear Model
#'
#' This function fits a linear model to multiple responses in an fMRI dataset.
#'
#' @param form The formula used to define the linear model.
#' @param data_env The environment containing the data to be used in the linear model.
#' @param conlist The list of contrasts used in the analysis.
#' @param vnames The names of the variables used in the linear model.
#' @param fcon The F-contrasts used in the analysis.
#' @param modmat The model matrix (default is \code{NULL}, which will calculate the model matrix using the formula).
#' @return A list containing the results from the multiresponse linear model analysis.
#' @keywords internal
multiresponse_lm <- function(form, data_env, conlist, vnames, fcon, modmat = NULL) {
  lm_fit <- if (is.null(modmat)) {
    lm(as.formula(form), data = data_env)
  } else {
    lm.fit(modmat, data_env$.y)
  }
  
  # Use the actual column names from the model matrix instead of vnames
  # This ensures the dimensions match correctly
  actual_vnames <- if (is.null(modmat)) {
    names(coef(lm_fit))
  } else {
    colnames(modmat)
  }
  
  fit_lm_contrasts(lm_fit, conlist, fcon, actual_vnames)
}





unpack_chunkwise <- function(cres, event_indices, baseline_indices) {
  # --- Beta Processing (Seems OK) ---
  cbetas <- lapply(cres, function(x) x$bstats) %>% dplyr::bind_rows()
  dat_beta <- cbetas$data %>% dplyr::bind_rows()
  
  # Check validity (assuming estimate column now correctly holds matrices)
  valid_estimates_idx <- sapply(dat_beta$estimate, function(x) !is.null(x) && is.matrix(x) && nrow(x) > 0)
  if (!any(valid_estimates_idx)) {
      stop("No valid beta estimates found across chunks in unpack_chunkwise.")
  }
  dat_beta_valid <- dat_beta[valid_estimates_idx, , drop = FALSE]

  # Concatenate beta results across chunks
  estimate_beta <- do.call(rbind, dat_beta_valid$estimate)
  se_beta <- do.call(rbind, dat_beta_valid$se)
  stat_beta <- do.call(rbind, dat_beta_valid$stat)
  prob_beta <- do.call(rbind, dat_beta_valid$prob)
  sigma_beta <- do.call(c, dat_beta_valid$sigma)
  
  # Re-package combined beta results
  cbetas_out <- dplyr::tibble(
    type = cbetas$type[1],
    stat_type = cbetas$stat_type[1],
    df.residual = cbetas$df.residual[1],
    conmat = list(NULL),
    colind = list(NULL),
    data = list(
      dplyr::tibble(
        estimate = list(estimate_beta),   
        se = list(se_beta),               
        stat = list(stat_beta),            
        prob = list(prob_beta),            
        sigma = list(sigma_beta)          
      )
    )
  )

  # --- Contrast Processing --- 
  ncon <- if (length(cres) > 0 && !is.null(cres[[1]]$contrasts) && length(cres[[1]]$contrasts) > 0) {
      length(cres[[1]]$contrasts)
  } else { 0 }

  if (ncon > 0) {
    contab <- lapply(cres, function(x) { 
        # Ensure contrasts is a list of tibbles, even if only one contrast
        cons <- x$contrasts 
        if (!is.list(cons)) cons <- list(cons) # Handle single contrast case if needed
        if (length(cons) > 0 && !is.null(names(cons))) { # Ensure names exist
             dplyr::bind_rows(cons, .id = "contrast_internal_name") # Requires names
        } else if (length(cons) > 0) {
             # Fallback if names are missing, might need adjustment based on actual structure
             warning("Contrast list per chunk lacks names, attempting bind_rows without .id")
             dplyr::bind_rows(cons)
        } else {
             dplyr::tibble() # Return empty tibble for chunks with no contrasts
        }
    }) %>% dplyr::bind_rows() # Bind results from all chunks

    # Check if contab is empty after binding
    if (nrow(contab) == 0) {
        con <- empty_contrast_table()
    } else {
        # Group by original contrast name and type
        # Use 'name' column if it exists, otherwise fallback might be needed
        grouping_vars <- intersect(c("name", "type"), names(contab))
        if (length(grouping_vars) == 0) stop("Cannot group contrasts: 'name' or 'type' column missing.")
        
        gsplit <- contab %>% dplyr::group_by(dplyr::across(dplyr::all_of(grouping_vars))) %>% dplyr::group_split()

        # Process each contrast group (combine results across chunks)
        con <- lapply(gsplit, function(g) {
            dat <- g$data %>% dplyr::bind_rows() 
            
            # Both paths now produce the same structure, so no conditional logic needed.
            # Simply assign the vectors directly.
            estimate_full <- dat$estimate
            se_full <- dat$se
            stat_full <- dat$stat
            prob_full <- dat$prob
            sigma_full <- if ("sigma" %in% names(dat)) dat$sigma else NULL
            
            # Re-package combined data using the canonical contrast schema:
            # one row per voxel and one atomic column per statistic.
            combined_data_tibble <- dplyr::tibble(
                estimate = estimate_full,
                se = se_full,
                stat = stat_full,
                prob = prob_full
            )
            if (!is.null(sigma_full)) {
                combined_data_tibble$sigma = sigma_full
            }

            # Take metadata from the first chunk's entry for this contrast
            g %>% dplyr::select(-data) %>% dplyr::slice_head() %>% 
                dplyr::mutate(data = list(combined_data_tibble))
                
        }) %>% dplyr::bind_rows() 
    }
  } else {
    con <- empty_contrast_table()
  }

  # --- DEBUG FINAL CONTRAST TIBBLE ---
  # message("Structure of final 'con' tibble before returning from unpack_chunkwise:")
  # print(str(con))
  # --- END DEBUG ---

  list(
    betas = cbetas_out,
    contrasts = con,
    event_indices = event_indices,
    baseline_indices = baseline_indices
  )
}




#' Print an fmri_lm_result object
#'
#' Provides a colorful and informative printout.
#'
#' @param x An fmri_lm_result object.
#' @param ... Additional arguments (unused).
#' @export
#' @rdname print
print.fmri_lm <- function(x, ...) {
  # optional: check if crayon is installed
  if (!requireNamespace("crayon", quietly = TRUE)) {
    # fallback to standard cat if crayon is missing
    cat("fmri_lm_result object (install 'crayon' for color)\n\n")
    
    cat("Model formula:\n",
        as.character(x$model$event_model$model_spec$formula), "\n")
    cat("Strategy: ", x$strategy, "\n")
    cat("Baseline parameters: ",
        ncol(design_matrix(x$model$baseline_model)), "\n")
    cat("Design parameters: ",
        ncol(design_matrix(x$model$event_model)), "\n")
    cat("Contrasts: ",
        paste(names(x$bcons), collapse = ", "), "\n\n")
    return(invisible(x))
  }
  
  # If we do have crayon, let's color it up:
  cat(crayon::blue$bold("\n==================================\n"))
  cat(crayon::blue$bold("        fmri_lm_result          \n"))
  cat(crayon::blue$bold("==================================\n\n"))
  
  # Print the model formula
  cat(crayon::green("Model formula:\n  "))
  cat(crayon::silver(as.character(x$model$event_model$model_spec$formula)), "\n\n")
  
  # Print strategy
  cat(crayon::green("Fitting strategy:  "))
  cat(crayon::silver(x$strategy), "\n\n")
  
  # Some stats about baseline, design, and contrasts
  bdim <- ncol(design_matrix(x$model$baseline_model))
  ddim <- ncol(design_matrix(x$model$event_model))
  
  cat(crayon::green("Baseline parameters: "), crayon::silver(bdim), "\n")
  cat(crayon::green("Design parameters:   "), crayon::silver(ddim), "\n")
  
  # If you have some # of simple contrasts
  c_names <- names(x$bcons)
  if (length(c_names) > 0) {
    cat(crayon::green("Contrasts:          "), crayon::silver(paste(c_names, collapse = ", ")), "\n\n")
  } else {
    cat(crayon::green("Contrasts:          "), crayon::silver("None\n\n"))
  }
  
  cat(crayon::yellow("Use coef(...), stats(...), etc. to extract results.\n\n"))
  
  invisible(x)
}
  
    
    
