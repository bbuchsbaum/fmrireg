#' @method get_formula fmri_model
#' @rdname get_formula
#' @export
get_formula.fmri_model <- function(x,...) {
  assert_that(inherits(x, "fmri_model"), msg = "'x' must be an 'fmri_model' object")
  term_names <- names(terms(x))
  form <- paste(".y ~", paste(term_names, collapse = " + "), "-1")
  return(as.formula(form))
}
#' Estimate fixed-order AR coefficients via fmriAR adapter with legacy fallback
#' @keywords internal
#' @noRd
.estimate_ar_parameters_routed <- function(residuals_vec, ar_order, run_indices = NULL, censor = NULL) {
  if (is.null(ar_order) || ar_order <= 0L) {
    return(numeric(0))
  }

  resid_mat <- as.matrix(residuals_vec)
  if (ncol(resid_mat) == 0L || nrow(resid_mat) <= ar_order) {
    return(rep(0, ar_order))
  }

  ar_cfg <- switch(as.integer(ar_order),
    `1` = list(struct = "ar1"),
    `2` = list(struct = "ar2"),
    list(struct = "arp", p = as.integer(ar_order))
  )

  plan <- tryCatch(
    .estimate_ar_via_fmriAR(
      residuals = resid_mat,
      cfg = ar_cfg,
      run_indices = run_indices,
      censor = censor
    ),
    error = function(e) NULL
  )

  phi <- NULL
  if (!is.null(plan)) {
    if (!is.null(plan$phi) && length(plan$phi) > 0L) {
      phi <- as.numeric(plan$phi[[1]])
    } else if (!is.null(plan$phi_by_parcel) && length(plan$phi_by_parcel) > 0L) {
      phi <- as.numeric(plan$phi_by_parcel[[1]])
    }
  }

  if (is.null(phi) || !length(phi)) {
    if (!isTRUE(getOption("fmrireg.ar.nonvoxel_legacy_fallback", FALSE))) {
      phi <- rep(0, ar_order)
    } else {
      phi <- estimate_ar_parameters(drop(resid_mat[, 1]), ar_order, censor = censor)
    }
  }

  phi <- as.numeric(phi)
  if (length(phi) < ar_order) {
    phi <- c(phi, rep(0, ar_order - length(phi)))
  } else if (length(phi) > ar_order) {
    phi <- phi[seq_len(ar_order)]
  }
  phi[!is.finite(phi)] <- 0
  phi
}

#' Fast row-wise robust regression for a single run
#'
#' Wrapper around robust_iterative_fitter for backward compatibility.
#' This function implements an IRLS algorithm using Huber or Tukey bisquare 
#' weights on the time-point residuals.
#'
#' @param X Design matrix (time points \eqn{\times} predictors)
#' @param Y Data matrix (time points \eqn{\times} voxels)
#' @param proj Preprojection list from \code{.fast_preproject(X)}
#' @param psi Psi function for weighting, either \code{"huber"} or
#'   \code{"bisquare"}
#' @param k_huber Tuning constant for Huber weights
#' @param c_tukey Tuning constant for Tukey bisquare weights
#' @param max_it Maximum number of IRLS iterations
#' @param sigma_fixed Optional fixed sigma value (for global scale estimation)
#' @keywords internal
#' @noRd
fast_rlm_run <- function(X, Y, proj,
                         psi = c("huber", "bisquare"),
                         k_huber = 1.345,
                         c_tukey = 4.685,
                         max_it = 2L,
                         sigma_fixed = NULL) {

  psi <- match.arg(psi)

  # Validate inputs: robust path should not accept NA/Inf
  if (anyNA(X) || any(!is.finite(X)) || anyNA(Y) || any(!is.finite(Y))) {
    stop("fast_rlm_run: X/Y contain NA/Inf values", call. = FALSE)
  }


  
  # Create initial GLM context
  if (missing(proj) || is.null(proj)) {
    proj <- .fast_preproject(X)
  }
  
  initial_ctx <- glm_context(X = X, Y = Y, proj = proj)
  
  # Create robust options
  cfg_robust_options <- list(
    type = psi,
    k_huber = k_huber,
    c_tukey = c_tukey,
    max_iter = max_it,
    scale_scope = if (is.null(sigma_fixed)) "run" else "global"
  )
  
  # Call robust_iterative_fitter
  result <- robust_iterative_fitter(
    initial_glm_ctx = initial_ctx,
    cfg_robust_options = cfg_robust_options,
    X_orig_for_resid = X,
    sigma_fixed = sigma_fixed
  )
  
  # Calculate standard errors
  se_beta <- sqrt(diag(result$XtWXi_final)) * result$sigma_robust_scale_final
  
  # Return in the expected format
  list(
    betas = result$betas_robust,
    se = se_beta,
    sigma = result$sigma_robust_scale_final,
    dfres = result$dfres,
    XtXinv = result$XtWXi_final,
    weights = result$robust_weights_final
  )
}

#' Create an fMRI Model
#'
#' This function creates an \code{fmri_model} by combining an event model and a baseline model.
#' If a baseline model is not provided, a default one is created based on the dataset.
#'
#' @param formula The model formula for experimental events.
#' @param block The model formula for block structure.
#' @param baseline_model (Optional) A \code{baseline_model} object. If \code{NULL}, a default baseline model is created.
#' @param dataset An \code{fmri_dataset} containing the event table and sampling frame.
#' @param drop_empty Logical. Whether to remove factor levels with zero size. Default is \code{TRUE}.
#' @param durations A vector of event durations. Default is \code{0}.
#' @return An \code{fmri_model} object.
#' @keywords internal
create_fmri_model <- function(formula, block, baseline_model = NULL, dataset, drop_empty = TRUE, durations = 0) {
  assert_that(is.formula(formula), msg = "'formula' must be a formula")
  formula <- .fmrireg_inject_registered_bases(formula)
  assert_that(is.formula(block), msg = "'block' must be a formula")
  assert_that(inherits(dataset, "fmri_dataset"), msg = "'dataset' must be an 'fmri_dataset'")
  assert_that(is.numeric(durations), msg = "'durations' must be numeric")
  
  if (is.null(baseline_model)) {
    baseline_model <- baseline_model(
      basis = "bs",
      degree = max(ceiling(median(fmrihrf::blocklens(dataset$sampling_frame)) / 100), 3),
      sframe = dataset$sampling_frame
    )
  } else {
    assert_that(inherits(baseline_model, "baseline_model"),
                msg = "'baseline_model' must have class 'baseline_model'")
  }

  ev_model <- event_model(
    formula_or_list = formula,
    block = block,
    data = dataset$event_table,
    sampling_frame = dataset$sampling_frame,
    drop_empty = drop_empty,
    durations = durations
  )
  
  fmri_model(ev_model, baseline_model, dataset)
}



#' Fit a Linear Regression Model for fMRI Data Analysis
#'
#' `fmri_lm` is a generic for fitting fMRI regression models. The
#' default interface accepts a model formula and dataset. An
#' alternative method can be used with a preconstructed
#' \code{fmri_model} object that already contains the design and data.
#'
#' @param formula A model formula describing the event structure or an
#'   \code{fmri_model} object.
#' @param ... Additional arguments passed to the chosen method.
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
  engine_ar_options <- dots$ar_options
  engine_robust_options <- dots$robust_options
  engine_cfg <- dots$cfg
  engine_args <- dots$engine_args

  dots$engine <- NULL
  dots$lowrank <- NULL
  dots$ar_options <- NULL
  dots$robust_options <- NULL
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
    engine_ar_options = engine_ar_options,
    engine_robust_options = engine_robust_options,
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
                                              robust_scale_scope = NULL,
                                              engine_robust_options = NULL) {
  robust_type <- if (is.logical(robust)) {
    if (robust) "huber" else FALSE
  } else {
    robust
  }

  robust_options <- robust_options %||% list()
  if (!is.null(engine_robust_options)) {
    robust_options <- utils::modifyList(robust_options, engine_robust_options)
  }
  if (!is.null(robust_type) && !("type" %in% names(robust_options))) {
    robust_options$type <- robust_type
  }
  if (!is.null(robust_psi) && !("type" %in% names(robust_options))) {
    robust_options$type <- robust_psi
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
                                          ar_voxelwise = FALSE,
                                          cor_struct = NULL,
                                          cor_iter = NULL,
                                          cor_global = NULL,
                                          ar1_exact_first = NULL,
                                          ar_p = NULL,
                                          engine_ar_options = NULL) {
  ar_options <- ar_options %||% list()
  if (!is.null(engine_ar_options)) {
    ar_options <- utils::modifyList(ar_options, engine_ar_options)
  }

  if (!is.null(cor_struct) && !("struct" %in% names(ar_options))) {
    ar_options$struct <- cor_struct
  }
  if (!is.null(cor_iter) && !("iter_gls" %in% names(ar_options))) {
    ar_options$iter_gls <- cor_iter
  }
  if (!is.null(cor_global) && !("global" %in% names(ar_options))) {
    ar_options$global <- cor_global
  }
  if (!is.null(ar1_exact_first) && !("exact_first" %in% names(ar_options))) {
    ar_options$exact_first <- ar1_exact_first
  }
  if (!is.null(ar_p) && !("p" %in% names(ar_options))) {
    ar_options$p <- ar_p
  }
  if (!("voxelwise" %in% names(ar_options))) {
    ar_options$voxelwise <- ar_voxelwise
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
                                  robust_options = NULL,
                                  robust_psi = NULL,
                                  robust_max_iter = NULL,
                                  robust_scale_scope = NULL,
                                  ar_options = NULL,
                                  ar_voxelwise = FALSE,
                                  cor_struct = NULL,
                                  cor_iter = NULL,
                                  cor_global = NULL,
                                  ar1_exact_first = NULL,
                                  ar_p = NULL,
                                  volume_weights_options = NULL,
                                  soft_subspace_options = NULL,
                                  volume_weights = NULL,
                                  nuisance_projection = NULL,
                                  engine_ar_options = NULL,
                                  engine_robust_options = NULL,
                                  engine_cfg = NULL) {
  robust_options <- .fmri_lm_normalize_robust_options(
    robust = robust,
    robust_options = robust_options,
    robust_psi = robust_psi,
    robust_max_iter = robust_max_iter,
    robust_scale_scope = robust_scale_scope,
    engine_robust_options = engine_robust_options
  )

  ar_options <- .fmri_lm_normalize_ar_options(
    ar_options = ar_options,
    ar_voxelwise = ar_voxelwise,
    cor_struct = cor_struct,
    cor_iter = cor_iter,
    cor_global = cor_global,
    ar1_exact_first = ar1_exact_first,
    ar_p = ar_p,
    engine_ar_options = engine_ar_options
  )

  preprocessing <- .fmri_lm_normalize_preprocessing_options(
    volume_weights_options = volume_weights_options,
    soft_subspace_options = soft_subspace_options,
    volume_weights = volume_weights,
    nuisance_projection = nuisance_projection
  )

  cfg <- if (!is.null(engine_cfg) && inherits(engine_cfg, "fmri_lm_config")) {
    engine_cfg
  } else {
    fmri_lm_control(
      robust_options = robust_options,
      ar_options = ar_options,
      volume_weights_options = preprocessing$volume_weights_options,
      soft_subspace_options = preprocessing$soft_subspace_options
    )
  }

  if (!is.null(engine_cfg) && inherits(engine_cfg, "fmri_lm_config")) {
    cfg$robust <- engine_cfg$robust
    cfg$ar <- utils::modifyList(cfg$ar, engine_cfg$ar)
  }

  cfg
}

#' @keywords internal
#' @noRd
.fmri_lm_attach_config_metadata <- function(fit, requested_cfg, executed_cfg = requested_cfg) {
  stopifnot(inherits(requested_cfg, "fmri_lm_config"))
  stopifnot(inherits(executed_cfg, "fmri_lm_config"))

  attr(fit, "requested_config") <- requested_cfg
  attr(fit, "executed_config") <- executed_cfg
  attr(fit, "config") <- executed_cfg
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

#' @keywords internal
#' @noRd
.fmri_lm_dispatch_engine <- function(model, dataset, engine, lowrank = NULL, cfg, engine_args = list()) {
  if (!is.null(lowrank) && is.null(engine_args$lowrank)) {
    engine_args$lowrank <- lowrank
  }

  spec <- .fmri_lm_resolve_engine_spec(engine = engine)
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

  fit <- .fmri_lm_attach_config_metadata(fit, requested_cfg = cfg, executed_cfg = executed_cfg)
  if (is.null(attr(fit, "strategy"))) {
    attr(fit, "strategy") <- spec$strategy
  }
  if (is.null(attr(fit, "engine"))) {
    attr(fit, "engine") <- spec$name
  }

  fit
}

#' @rdname fmri_lm
#' @export
#' @param block The model formula for block structure.
#' @param baseline_model (Optional) A \code{baseline_model} object. Default is \code{NULL}.
#' @param dataset An \code{fmri_dataset} object containing the time-series data.
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
#' @param use_fast_path Logical. If \code{TRUE}, use matrix-based computation for speed. Default is \code{FALSE}.
#' @param progress Logical. Whether to display a progress bar during model fitting. Default is \code{FALSE}.
#' @param ar_voxelwise Logical. Estimate AR parameters voxel-wise (overrides \code{ar_options$voxelwise}).
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
#'   \item \code{iter_gls}: Integer. Number of GLS iterations (default: 1)
#'   \item \code{global}: Logical. Use global AR coefficients (default: FALSE)
#'   \item \code{voxelwise}: Logical. Estimate AR parameters voxel-wise (default: FALSE)
#'   \item \code{exact_first}: Logical. Apply exact AR(1) scaling to first sample (default: FALSE)
#' }
#' 
#' @export
#' @seealso \code{\link{fmri_dataset}}, \code{\link{fmri_lm_fit}}, \code{\link{fmri_lm_control}}
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
#' dset <- matrix_dataset(cbind(ys1,ys2), TR=2, 
#'                        run_length=fmrihrf::blocklens(sframe), 
#'                        event_table=facedes)
#' flm <- fmri_lm(onset ~ hrf(face_gen, 
#'                            basis=gen_hrf(fmrihrf::hrf_bspline, N=7, span=25)), 
#'                block = ~ run, 
#'                strategy="chunkwise", nchunks=1, dataset=dset)
#' 
fmri_lm.formula <- function(formula, block, baseline_model = NULL, dataset, durations = 0, drop_empty = TRUE,
                         robust = FALSE, robust_options = NULL, ar_options = NULL,
                         volume_weights_options = NULL, soft_subspace_options = NULL,
                         strategy = c("runwise", "chunkwise"), nchunks = 10, use_fast_path = FALSE, progress = FALSE,
                         ar_voxelwise = FALSE,
                         parallel_voxels = FALSE,
                    # Individual AR parameters for backward compatibility
                    cor_struct = NULL, cor_iter = NULL, cor_global = NULL,
                    ar1_exact_first = NULL, ar_p = NULL,
                    # Individual robust parameters for backward compatibility
                    robust_psi = NULL, robust_max_iter = NULL, robust_scale_scope = NULL,
                    # Convenience parameters for preprocessing features
                    volume_weights = NULL, nuisance_projection = NULL,
                    parallel_chunks = FALSE,
                    ...) {
  engine_ctx <- .fmri_lm_extract_engine_context(list(...))
  engine <- engine_ctx$engine
  lowrank <- engine_ctx$lowrank
  engine_ar_options <- engine_ctx$engine_ar_options
  engine_robust_options <- engine_ctx$engine_robust_options
  engine_cfg <- engine_ctx$engine_cfg
  engine_args <- engine_ctx$engine_args

  strategy <- match.arg(strategy)
  formula <- .fmrireg_inject_registered_bases(formula)
  
  # Error checking
  assert_that(is.formula(formula), msg = "'formula' must be a formula")
  assert_that(is.formula(block), msg = "'block' must be a formula")
  assert_that(inherits(dataset, "fmri_dataset"), msg = "'dataset' must be an 'fmri_dataset'")
  assert_that(is.numeric(durations), msg = "'durations' must be numeric")
  assert_that(is.logical(drop_empty), msg = "'drop_empty' must be logical")
  assert_that(is.logical(robust) || robust %in% c("huber", "bisquare"), 
              msg = "'robust' must be logical or one of 'huber', 'bisquare'")
  assert_that(is.logical(use_fast_path), msg = "'use_fast_path' must be logical")
  assert_that(is.logical(ar_voxelwise), msg = "'ar_voxelwise' must be logical")
  assert_that(
    is.logical(parallel_chunks) && length(parallel_chunks) == 1L && !is.na(parallel_chunks),
    msg = "'parallel_chunks' must be TRUE or FALSE"
  )
  if (strategy == "chunkwise") {
    assert_that(is.numeric(nchunks) && nchunks > 0, msg = "'nchunks' must be a positive number")
  }

  # Create configuration object
  cfg <- .fmri_lm_build_config(
    robust = robust,
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
    engine_ar_options = engine_ar_options,
    engine_robust_options = engine_robust_options,
    engine_cfg = engine_cfg
  )
  
  model <- create_fmri_model(formula, block, baseline_model, dataset, durations = durations, drop_empty = drop_empty)

  if (!is.null(engine)) {
    return(.fmri_lm_dispatch_engine(
      model = model,
      dataset = dataset,
      engine = engine,
      lowrank = lowrank,
      cfg = cfg,
      engine_args = engine_args
    ))
  }
  
  # Pass configuration object down
  # Note: We don't pass ... here because all parameters have been processed
  # and included in the cfg object
  ret <- fmri_lm_fit(model, dataset, strategy, cfg, nchunks,
                     use_fast_path = use_fast_path, progress = progress,
                     parallel_voxels = parallel_voxels,
                     parallel_chunks = parallel_chunks)
  return(ret)
}

#' @rdname fmri_lm
#' @export
fmri_lm.fmri_model <- function(formula, dataset = NULL,
                               robust = FALSE, robust_options = NULL,
                               ar_options = NULL,
                               volume_weights_options = NULL, soft_subspace_options = NULL,
                               strategy = c("runwise", "chunkwise"), nchunks = 10,
                               use_fast_path = FALSE, progress = FALSE,
                               ar_voxelwise = FALSE, parallel_voxels = FALSE,
                               cor_struct = NULL, cor_iter = NULL,
                               cor_global = NULL, ar1_exact_first = NULL,
                               ar_p = NULL,
                               robust_psi = NULL, robust_max_iter = NULL,
                               robust_scale_scope = NULL,
                               volume_weights = NULL, nuisance_projection = NULL,
                               parallel_chunks = FALSE,
                               ...) {
  engine_ctx <- .fmri_lm_extract_engine_context(list(...))
  engine <- engine_ctx$engine
  lowrank <- engine_ctx$lowrank
  engine_ar_options <- engine_ctx$engine_ar_options
  engine_robust_options <- engine_ctx$engine_robust_options
  engine_cfg <- engine_ctx$engine_cfg
  engine_args <- engine_ctx$engine_args
  strategy <- match.arg(strategy)
  assert_that(inherits(formula, "fmri_model"))
  assert_that(
    is.logical(parallel_chunks) && length(parallel_chunks) == 1L && !is.na(parallel_chunks),
    msg = "'parallel_chunks' must be TRUE or FALSE"
  )

  dataset <- dataset %||% formula$dataset %||% attr(formula, "dataset")
  if (is.null(dataset)) {
    stop("No dataset found in 'formula' and none supplied.")
  }
  assert_that(inherits(dataset, "fmri_dataset"))

  cfg <- .fmri_lm_build_config(
    robust = robust,
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
    engine_ar_options = engine_ar_options,
    engine_robust_options = engine_robust_options,
    engine_cfg = engine_cfg
  )

  if (!is.null(engine)) {
    return(.fmri_lm_dispatch_engine(
      model = formula,
      dataset = dataset,
      engine = engine,
      lowrank = lowrank,
      cfg = cfg,
      engine_args = engine_args
    ))
  }

  ret <- fmri_lm_fit(formula, dataset, strategy, cfg, nchunks,
                     use_fast_path = use_fast_path, progress = progress,
                     parallel_voxels = parallel_voxels,
                     parallel_chunks = parallel_chunks)
  return(ret)
}


#' Fit an fMRI Linear Regression Model with a Specified Fitting Strategy
#'
#' This function fits an fMRI linear regression model using the specified \code{fmri_model} object, dataset,
#' and data splitting strategy (either \code{"runwise"} or \code{"chunkwise"}). It is primarily an internal function
#' used by the \code{fmri_lm} function.
#'
#' @param fmrimod An \code{fmri_model} object.
#' @param dataset An \code{fmri_dataset} object containing the time-series data.
#' @param strategy The data splitting strategy, either \code{"runwise"} or \code{"chunkwise"}. Default is \code{"runwise"}.
#' @param cfg An \code{fmri_lm_config} object containing all fitting options. See \code{\link{fmri_lm_control}}.
#' @param nchunks Number of data chunks when strategy is \code{"chunkwise"}.
#'   This controls memory partitioning; chunks are processed sequentially unless
#'   \code{parallel_chunks = TRUE}. Default is \code{10}.
#' @param use_fast_path Logical. If \code{TRUE}, use matrix-based computation for speed. Default is \code{FALSE}.
#' @param progress Logical. Whether to display a progress bar during model fitting. Default is \code{FALSE}.
#' @param parallel_voxels Logical. If TRUE, voxelwise AR processing within runs
#'   is parallelised using `future.apply`. Default is \code{FALSE}.
#' @param parallel_chunks Logical. If TRUE and \code{strategy = "chunkwise"},
#'   process chunks with \code{future.apply::future_lapply()} using the active
#'   \code{future} plan. Default is \code{FALSE}.
#' @param ... Additional arguments.
#' @return A fitted fMRI linear regression model with the specified fitting strategy.
#' @keywords internal
#' @seealso \code{\link{fmri_lm}}, \code{\link{fmri_model}}, \code{\link{fmri_dataset}}
fmri_lm_fit <- function(fmrimod, dataset, strategy = c("runwise", "chunkwise"),
                        cfg, nchunks = 10, use_fast_path = FALSE, progress = FALSE,
                        parallel_voxels = FALSE, parallel_chunks = FALSE, ...) {
  strategy <- match.arg(strategy)
  
  # Validate config
  assert_that(inherits(cfg, "fmri_lm_config"), msg = "'cfg' must be an 'fmri_lm_config' object")
  
  # Error checking
  assert_that(inherits(fmrimod, "fmri_model"), msg = "'fmrimod' must be an 'fmri_model' object")
  assert_that(inherits(dataset, "fmri_dataset"), msg = "'dataset' must be an 'fmri_dataset' object")
  assert_that(is.logical(use_fast_path), msg = "'use_fast_path' must be logical")
  assert_that(
    is.logical(parallel_chunks) && length(parallel_chunks) == 1L && !is.na(parallel_chunks),
    msg = "'parallel_chunks' must be TRUE or FALSE"
  )
  if (strategy == "chunkwise") {
    assert_that(is.numeric(nchunks) && nchunks > 0, msg = "'nchunks' must be a positive number")
  } else if (isTRUE(parallel_chunks)) {
    warning("'parallel_chunks' only applies when strategy = 'chunkwise'.", call. = FALSE)
  }

  # Check for redundant nuisance regressors when soft subspace projection is enabled
  if (cfg$soft_subspace$enabled && cfg$soft_subspace$warn_redundant) {
    .check_redundant_nuisance(fmrimod$baseline_model, cfg$soft_subspace)
  }

  contrast_prep <- prepare_fmri_lm_contrasts(fmrimod)
  processed_conlist <- contrast_prep$processed
  standard_path_conlist <- contrast_prep$standard
  contrast_objects <- processed_conlist
  
  # Pass the full processed contrast objects list down.
  # The fitting function (chunkwise/runwise) will decide whether to use the objects (slow path)
  # or extract weights (fast path).
  phi_global <- NULL
  sigma_global <- NULL
  if (cfg$ar$global && cfg$ar$struct != "iid") {
    ar_order <- switch(cfg$ar$struct,
                       ar1 = 1L,
                       ar2 = 2L,
                       arp = cfg$ar$p)

    chunk_iter <- exec_strategy("runwise")(dataset)
    run_chunks <- collect_chunks(chunk_iter)
    
    form <- get_formula(fmrimod)
    resid_chunks <- vector("list", length(run_chunks))
    for (ri in seq_along(run_chunks)) {
      rch <- run_chunks[[ri]]
      tmats_run <- term_matrices(fmrimod, rch$chunk_num)
      data_env_run <- list2env(tmats_run)
      n_time_run <- nrow(tmats_run[[1]])
      data_env_run[[".y"]] <- rep(0, n_time_run)
      X_run <- model.matrix(form, data_env_run)
      proj_run <- .fast_preproject(X_run)
      Y_run <- as.matrix(rch$data)
      ols <- .fast_lm_matrix(X_run, Y_run, proj_run, return_fitted = TRUE)
      resid_chunks[[ri]] <- rowMeans(Y_run - ols$fitted)
    }
    resid_vec <- unlist(resid_chunks, use.names = FALSE)
    phi_global <- .estimate_ar_parameters_routed(resid_vec, ar_order)
    cfg$ar$iter_gls <- 1L
  }

  if (cfg$robust$type != FALSE && cfg$robust$scale_scope == "global") {
    chunk_iter <- exec_strategy("runwise")(dataset)
    run_chunks <- collect_chunks(chunk_iter)
    form <- get_formula(fmrimod)
    row_med_chunks <- vector("list", length(run_chunks))
    for (ri in seq_along(run_chunks)) {
      rch <- run_chunks[[ri]]
      tmats_run <- term_matrices(fmrimod, rch$chunk_num)
      data_env_run <- list2env(tmats_run)
      n_time_run <- nrow(tmats_run[[1]])
      data_env_run[[".y"]] <- rep(0, n_time_run)
      X_run <- model.matrix(form, data_env_run)
      proj_run <- .fast_preproject(X_run)
      Y_run <- as.matrix(rch$data)
      ols <- .fast_lm_matrix(X_run, Y_run, proj_run, return_fitted = TRUE)
      row_med_chunks[[ri]] <- matrixStats::rowMedians(abs(Y_run - ols$fitted))
    }
    row_med_all <- unlist(row_med_chunks, use.names = FALSE)
    sigma_global <- 1.4826 * median(row_med_all)
    if (sigma_global <= .Machine$double.eps) sigma_global <- .Machine$double.eps
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
                    if (inherits(dataset, "latent_dataset")) {
                      if (isTRUE(parallel_chunks)) {
                        warning(
                          "'parallel_chunks' is not currently supported for latent_dataset; ",
                          "using sequential chunkwise fitting.",
                          call. = FALSE
                        )
                      }
                      chunkwise_lm(dataset, fmrimod, standard_path_conlist, # Pass full objects
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
  
  ret <- list(
    result = result,
    model = fmrimod,
    strategy = strategy,
    bcons = processed_conlist,
    dataset = dataset,
    ar_coef = result$ar_coef
  )
  
  class(ret) <- "fmri_lm"
  
  ret <- .fmri_lm_attach_config_metadata(ret, requested_cfg = cfg, executed_cfg = cfg)
  attr(ret, "strategy") <- strategy
  
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



#' Reshape Coefficient Data
#'
#' This function reshapes coefficient data from wide to long format and merges it with design information.
#'
#' @param df A data frame containing coefficient estimates.
#' @param des A data frame containing design information.
#' @param measure The name of the value column in the reshaped data. Default is \code{"value"}.
#' @return A data frame in long format with merged design information.
#' @keywords internal
#' @autoglobal
reshape_coef <- function(df, des, measure = "value") {
  # Create a unique identifier for each row
  df <- df %>% dplyr::mutate(row_id = dplyr::row_number())
  
  des <- des %>% dplyr::mutate(key = do.call(paste, c(.[, colnames(des)], sep = ":")))
  
  colnames(df)[-ncol(df)] <- des$key  # assign new column names excluding the last column (row_id)
  
  df_long <- df %>%
    tidyr::pivot_longer(-row_id, names_to = "col_name", values_to = measure)
  
  # Match the long dataframe with the design dataframe
  df_long <- dplyr::left_join(df_long, des, by = c("col_name" = "key"))
  
  return(df_long)
}


#' Extract Statistical Measures from an fmri_lm Object
#'
#' This function extracts statistical measures (e.g., estimates, standard errors) from an \code{fmri_lm} object.
#'
#' @param x An \code{fmri_lm} object.
#' @param type The type of statistic to extract: \code{"betas"}, \code{"contrasts"}, or \code{"F"}.
#' @param element The specific element to extract, such as \code{"estimate"}, \code{"se"}, \code{"stat"}, or \code{"prob"}.
#' @return A tibble containing the requested statistical measures.
#' @keywords internal
pull_stat_revised <- function(x, type, element) {
  if (type == "betas") {
    # Ensure we access the matrix correctly from the list structure
    beta_matrix <- x$result$betas$data[[1]]$estimate[[1]]
    ret <- beta_matrix[, x$result$event_indices, drop = FALSE]
    colnames(ret) <- conditions(x$model$event_model)
    suppressMessages(tibble::as_tibble(ret, .name_repair = "check_unique"))
  } else if (type == "contrasts") {
    ret <- x$result$contrasts %>% dplyr::filter(type == "contrast")
    if (nrow(ret) == 0) {
      stop("No simple contrasts for this model.")
    }
    cnames <- ret$name
    # Extract the specific element (e.g., estimate), which is a list(vector)
    # Then extract the vector itself (element [[1]]) before binding
    out <- lapply(ret$data, function(inner_tibble) inner_tibble[[element]][[1]]) %>% 
             dplyr::bind_cols()
    names(out) <- cnames
    out
  } else if (type == "F") {
    ret <- x$result$contrasts %>% dplyr::filter(type == "Fcontrast")
    if (nrow(ret) == 0) {
      stop("No F contrasts for this model.")
    }
    cnames <- ret$name
    # Extract the specific element (e.g., estimate), which is list(vector)
    # Then extract the vector itself (element [[1]]) before binding
    out <- lapply(ret$data, function(inner_tibble) inner_tibble[[element]][[1]]) %>% 
             dplyr::bind_cols()
    names(out) <- cnames
    out
  } else {
    stop("Invalid type specified. Must be 'betas', 'contrasts', or 'F'.")
  }
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
  } else if (type == "contrasts") {
    ret <- x$result$contrasts %>% dplyr::filter(type == "contrast")
    if (nrow(ret) == 0) {
      stop("No simple contrasts for this model.")
    }
    cnames <- ret$name
    out <- lapply(ret$data, function(x) x[[element]]) %>% dplyr::bind_cols()
    names(out) <- cnames
    out
  } else if (type == "F") {
    ret <- x$result$contrasts %>% dplyr::filter(type == "Fcontrast")
    if (nrow(ret) == 0) {
      stop("No F contrasts for this model.")
    }
    cnames <- ret$name
    out <- lapply(ret$data, function(x) x[[element]]) %>% dplyr::bind_cols()
    names(out) <- cnames
    out
  } else {
    stop("Invalid type specified. Must be 'betas', 'contrasts', or 'F'.")
  }
}

#' @method coef fmri_lm
#' @export
coef.fmri_lm <- function(object, type = c("betas", "contrasts"), include_baseline = FALSE, recon = FALSE, ...) {
  type <- match.arg(type)
  
  if (type == "contrasts") {
    # Contrast handling remains the same
    res <- pull_stat(object, "contrasts", "estimate")
  } else if (type == "betas") {
    # Get all beta estimates first
    all_betas <- object$result$betas$data[[1]]$estimate[[1]]
    
    if (include_baseline) {
      # Return all betas, ensure correct names from the full design matrix
      res <- all_betas
      dm_colnames <- colnames(design_matrix(object$model))
      if (!is.null(dm_colnames)) {
        if (length(dm_colnames) == ncol(res)) {
          colnames(res) <- dm_colnames
        } else {
          beta_colind <- tryCatch(object$result$betas$colind[[1]], error = function(e) NULL)
          if (!is.null(beta_colind) &&
              length(beta_colind) == ncol(res) &&
              max(beta_colind) <= length(dm_colnames)) {
            colnames(res) <- dm_colnames[beta_colind]
          } else {
            colnames(res) <- make.names(paste0("beta_", seq_len(ncol(res))), unique = TRUE)
          }
        }
      }
      # Convert back to tibble for consistency if needed, though matrix might be better here
      # res <- as_tibble(res)
    } else {
      # Default: return only event betas
      # Check bounds and filter valid indices
      max_col <- ncol(all_betas)
      valid_event_indices <- object$result$event_indices[object$result$event_indices <= max_col]
      
      if (length(valid_event_indices) == 0) {
        warning("No valid event indices found in coef.fmri_lm. Using all available columns.")
        valid_event_indices <- 1:max_col
      }
      
      res <- all_betas[, valid_event_indices, drop = FALSE]
      
      # Use the actual column names from the design matrix instead of conditions()
      # This avoids duplicate names when multiple terms have the same variables
      dm <- design_matrix(object$model)
      if (!is.null(dm) && ncol(dm) >= max(valid_event_indices)) {
        actual_colnames <- colnames(dm)[valid_event_indices]
        colnames(res) <- actual_colnames
      } else {
        # Fallback: use conditions but make them unique
        condition_names <- conditions(object$model$event_model)[1:length(valid_event_indices)]
        colnames(res) <- make.names(condition_names, unique = TRUE)
      }
      
      # Return as matrix - transpose to get conditions x voxels
      res <- t(res)
    }
  } else {
    # Should not happen due to match.arg, but defensive coding
    stop("Invalid type specified.")
  }
  
  # Reconstruction functionality can be added here if necessary (applies to the 'res' matrix/tibble)
  # if (recon && inherits(object$dataset, "fmri_dataset")) { ... }
  
  return(res)
}

#' @method stats fmri_lm
#' @rdname stats
#' @export
stats.fmri_lm <- function(x, type = c("estimates", "contrasts", "F"), ...) {
  type <- match.arg(type)
  if (type == "estimates") {
    pull_stat(x, "betas", "stat")
  } else if (type == "contrasts") {
    pull_stat(x, "contrasts", "stat")
  } else if (type == "F") {
    pull_stat(x, "F", "stat")
  }
}

#' @method standard_error fmri_lm
#' @rdname standard_error
#' @export
standard_error.fmri_lm <- function(x, type = c("estimates", "contrasts"),...) {
  type <- match.arg(type)
  if (type == "estimates") {
    pull_stat(x, "betas", "se")
  } else if (type == "contrasts") {
    pull_stat(x, "contrasts", "se")
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
  conres <- if (!is.null(conlist)) {
    ret <- lapply(conlist, function(con) {
      # Extract colind from the contrast object's attributes
      colind <- attr(con, "colind")
      if (is.null(colind)) {
        warning(paste("Missing colind attribute for contrast:", con$name %||% "unnamed"))
        return(NULL) # Skip this contrast
      }
      estimate_contrast(con, fit, colind)
    })
    # Filter out NULL results
    ret <- ret[!sapply(ret, is.null)]
    names(ret) <- sapply(conlist[!sapply(ret, is.null)], function(x) x$name %||% "unnamed")
    ret
  } else {
    list()
  }
  
  bstats <- beta_stats(fit, vnames, se = se)
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
        con <- dplyr::tibble()
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
            
            # Re-package combined data for this contrast
            combined_data_tibble <- dplyr::tibble(
                estimate = list(estimate_full), 
                se = list(se_full),             
                stat = list(stat_full),          
                prob = list(prob_full)           
            )
            if (!is.null(sigma_full)) {
                combined_data_tibble$sigma = list(sigma_full)
            }

            # Take metadata from the first chunk's entry for this contrast
            g %>% dplyr::select(-data) %>% dplyr::slice_head() %>% 
                dplyr::mutate(data = list(combined_data_tibble))
                
        }) %>% dplyr::bind_rows() 
    }
  } else {
    con <- dplyr::tibble() # Return empty tibble if no contrasts
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




#' Perform Chunkwise Linear Modeling on fMRI Dataset
#'
#' This function performs a chunkwise linear model analysis on an fMRI dataset,
#' splitting the dataset into chunks and running the linear model on each chunk.
#'
#' @param x An \code{fmri_dataset} object.
#' @param model The \code{fmri_model} used for the analysis.
#' @param contrast_objects The list of full contrast objects.
#' @param nchunks The number of chunks to divide the dataset into.
#' @param cfg An \code{fmri_lm_config} object containing all fitting options.
#' @param verbose Logical. Whether to display progress messages (default is \code{FALSE}).
#' @param use_fast_path Logical. If \code{TRUE}, use matrix-based computation for speed. Default is \code{FALSE}.
#' @param progress Logical. Display a progress bar for chunk processing. Default is \code{FALSE}.
#' @param parallel_chunks Logical. If \code{TRUE}, process chunks with
#'   \code{future.apply::future_lapply()} using the active future plan.
#' @param phi_fixed Optional fixed AR parameters.
#' @param sigma_fixed Optional fixed robust scale estimate.
#' @return A list containing the unpacked chunkwise results.
#' @keywords internal
chunkwise_lm.fmri_dataset_old <- function(x, model, contrast_objects, nchunks, cfg,
                                      verbose = FALSE, use_fast_path = FALSE, progress = FALSE,
                                      parallel_chunks = FALSE,
                                      phi_fixed = NULL,
                                      sigma_fixed = NULL, ...) {
  # Legacy shim retained for callers that still reference the historical helper.
  chunkwise_lm.fmri_dataset(
    x = x,
    model = model,
    contrast_objects = contrast_objects,
    nchunks = nchunks,
    cfg = cfg,
    verbose = verbose,
    use_fast_path = use_fast_path,
    progress = progress,
    parallel_chunks = parallel_chunks,
    phi_fixed = phi_fixed,
    sigma_fixed = sigma_fixed,
    ...
  )
}



#' Perform Runwise Linear Modeling on fMRI Dataset
#'
#' This function performs a runwise linear model analysis on an fMRI dataset by
#' running the linear model for each data run and combining the results.
#'
#' @param dset An \code{fmri_dataset} object.
#' @param model The \code{fmri_model} used for the analysis.
#' @param contrast_objects The list of full contrast objects.
#' @param cfg An \code{fmri_lm_config} object containing all fitting options.
#' @param verbose Logical. Whether to display progress messages (default is \code{FALSE}).
#' @param use_fast_path Logical. Whether to use fast path computation (default is \code{FALSE}).
#' @param progress Logical. Display a progress bar for run processing. Default is \code{FALSE}.
#' @param phi_fixed Optional fixed AR parameters.
#' @param sigma_fixed Optional fixed robust scale estimate.
#' @param parallel_voxels Logical. If TRUE, process voxels in parallel using
#'   `future.apply`. Default is \code{FALSE}.
#' @return A list containing the combined results from runwise linear model analysis.
#' @keywords internal
#' @autoglobal
runwise_lm <- function(dset, model, contrast_objects, cfg, verbose = FALSE,
                       use_fast_path = FALSE, progress = FALSE,
                       phi_fixed = NULL,
                       sigma_fixed = NULL,
                       parallel_voxels = FALSE) {
  runwise_lm_impl(
    dset = dset,
    model = model,
    contrast_objects = contrast_objects,
    cfg = cfg,
    verbose = verbose,
    use_fast_path = use_fast_path,
    progress = progress,
    phi_fixed = phi_fixed,
    sigma_fixed = sigma_fixed,
    parallel_voxels = parallel_voxels
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
  
    
    
