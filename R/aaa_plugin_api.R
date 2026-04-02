#' Plugin registration and helpers for fmrireg
#'
.fmrireg_engine_registry <- new.env(parent = emptyenv())
.fmrireg_basis_registry <- new.env(parent = emptyenv())

#' Register a plugin engine for `fmri_lm`
#'
#' @param name Character scalar identifier advertised to users (e.g. "friman").
#' @param fit Function invoked as `fit(model, dataset, args, cfg)` and expected
#'   to return an `fmri_lm` object.
#' @param preflight Optional function invoked before fitting; receives the same
#'   arguments as `fit` and can signal errors early.
#' @param capabilities Optional named list describing engine support for global
#'   `fmri_lm()` options. Recognized fields currently include `robust`,
#'   `preprocessing`, `ar_voxelwise`, `ar_by_cluster`, plus contextual rules
#'   such as `requires_event_regressors`, `requires_parcels_for_by_cluster`,
#'   and `forbid_by_cluster_dataset_classes`.
#' @return Invisibly, `TRUE`.
#' @export
register_engine <- function(name, fit, preflight = NULL, capabilities = list()) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  if (!is.function(fit)) {
    stop("`fit` must be a function", call. = FALSE)
  }
  if (!is.null(preflight) && !is.function(preflight)) {
    stop("`preflight` must be NULL or a function", call. = FALSE)
  }
  if (is.null(capabilities)) {
    capabilities <- list()
  }
  .fmrireg_engine_registry[[name]] <- list(
    fit = fit,
    preflight = preflight,
    capabilities = capabilities
  )
  invisible(TRUE)
}

#' @keywords internal
get_engine <- function(name) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    return(NULL)
  }
  .fmrireg_engine_registry[[name]]
}

#' @keywords internal
.builtin_engine_aliases <- function(name) {
  if (identical(name, "latent_sketch")) {
    "sketch"
  } else {
    character()
  }
}

#' @keywords internal
.builtin_engine_source <- function(name) {
  if (name %in% c("latent_sketch", "rrr_gls")) {
    "builtin"
  } else {
    "plugin"
  }
}

#' @keywords internal
.new_engine_spec <- function(name, registration, include_functions = FALSE) {
  spec <- list(
    name = name,
    source = .builtin_engine_source(name),
    aliases = .builtin_engine_aliases(name),
    strategy = if (identical(name, "latent_sketch")) "sketch" else "engine",
    capabilities = .normalize_engine_capabilities(registration$capabilities),
    has_preflight = is.function(registration$preflight)
  )

  if (isTRUE(include_functions)) {
    spec$fit <- registration$fit
    spec$preflight <- registration$preflight
  }

  class(spec) <- c("fmrireg_engine_spec", "list")
  spec
}

#' @keywords internal
get_engine_spec <- function(name, include_functions = FALSE) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    return(NULL)
  }

  registration <- get_engine(name)
  if (is.null(registration)) {
    return(NULL)
  }

  .new_engine_spec(name, registration, include_functions = include_functions)
}

#' @keywords internal
list_engine_specs <- function(include_functions = FALSE) {
  engine_names <- sort(ls(.fmrireg_engine_registry, all.names = TRUE))
  specs <- lapply(engine_names, get_engine_spec, include_functions = include_functions)
  names(specs) <- engine_names
  specs
}

#' Inspect a Registered Engine Specification
#'
#' Returns a read-only description of a registered engine, including its
#' normalized capabilities, source (`"builtin"` vs `"plugin"`), aliases, and
#' dispatch strategy. This is intended for extension authors and diagnostic
#' tooling; it does not expose the underlying fit/preflight functions.
#'
#' @param name Engine name, such as `"rrr_gls"` or `"latent_sketch"`.
#' @return A list of class `fmrireg_engine_spec`, or `NULL` if the engine is not
#'   registered.
#' @export
engine_spec <- function(name) {
  get_engine_spec(name, include_functions = FALSE)
}

#' List Registered Engine Specifications
#'
#' Returns the read-only specifications for all currently registered engines.
#'
#' @return A named list of `fmrireg_engine_spec` objects.
#' @export
engine_specs <- function() {
  list_engine_specs(include_functions = FALSE)
}

#' @export
print.fmrireg_engine_spec <- function(x, ...) {
  stopifnot(inherits(x, "fmrireg_engine_spec"))

  aliases <- if (length(x$aliases) > 0L) {
    paste(x$aliases, collapse = ", ")
  } else {
    "<none>"
  }

  caps <- x$capabilities %||% list()
  capability_summary <- c(
    paste0("robust=", caps$robust),
    paste0("preprocessing=", caps$preprocessing),
    paste0("ar_voxelwise=", caps$ar_voxelwise),
    paste0("ar_by_cluster=", caps$ar_by_cluster)
  )

  cat("<fmrireg_engine_spec>\n", sep = "")
  cat("name: ", x$name, "\n", sep = "")
  cat("source: ", x$source, " | strategy: ", x$strategy, "\n", sep = "")
  cat("aliases: ", aliases, "\n", sep = "")
  cat("capabilities: ", paste(capability_summary, collapse = ", "), "\n", sep = "")

  if (isTRUE(caps$requires_event_regressors)) {
    cat("requires: event regressors\n", sep = "")
  }
  if (isTRUE(caps$requires_parcels_for_by_cluster)) {
    cat("requires: parcels when by_cluster = TRUE\n", sep = "")
  }
  if (length(caps$forbid_by_cluster_dataset_classes) > 0L) {
    cat(
      "forbids by_cluster for: ",
      paste(caps$forbid_by_cluster_dataset_classes, collapse = ", "),
      "\n",
      sep = ""
    )
  }

  invisible(x)
}

#' @keywords internal
.normalize_engine_capabilities <- function(capabilities = NULL) {
  defaults <- list(
    robust = TRUE,
    preprocessing = TRUE,
    ar_voxelwise = TRUE,
    ar_by_cluster = TRUE,
    requires_event_regressors = FALSE,
    requires_parcels_for_by_cluster = FALSE,
    forbid_by_cluster_dataset_classes = character()
  )

  capabilities <- capabilities %||% list()
  utils::modifyList(defaults, capabilities)
}

#' @keywords internal
.validate_engine_capabilities <- function(engine, cfg, capabilities = NULL) {
  stopifnot(is.character(engine), length(engine) == 1L, nzchar(engine))
  stopifnot(inherits(cfg, "fmri_lm_config"))

  caps <- .normalize_engine_capabilities(capabilities)

  if (identical(caps$robust, FALSE) && !identical(cfg$robust$type %||% FALSE, FALSE)) {
    stop(sprintf("%s does not support robust fitting; set robust = FALSE", engine), call. = FALSE)
  }

  preprocessing_enabled <- isTRUE(cfg$volume_weights$enabled) || isTRUE(cfg$soft_subspace$enabled)
  if (identical(caps$preprocessing, FALSE) && preprocessing_enabled) {
    stop(sprintf("%s does not support volume_weights or soft_subspace preprocessing", engine), call. = FALSE)
  }

  if (identical(caps$ar_voxelwise, FALSE) && isTRUE(cfg$ar$voxelwise)) {
    stop(sprintf("%s supports only shared (non-voxelwise) temporal covariance", engine), call. = FALSE)
  }

  if (identical(caps$ar_by_cluster, FALSE) && isTRUE(cfg$ar$by_cluster)) {
    stop(sprintf("%s does not support parcel-specific AR whitening", engine), call. = FALSE)
  }

  invisible(caps)
}

#' @keywords internal
.validate_engine_context <- function(engine, model, dataset, args, cfg, capabilities = NULL) {
  stopifnot(is.character(engine), length(engine) == 1L, nzchar(engine))
  stopifnot(inherits(model, "fmri_model"))
  stopifnot(inherits(cfg, "fmri_lm_config"))

  caps <- .normalize_engine_capabilities(capabilities)
  args <- args %||% list()

  if (isTRUE(caps$requires_event_regressors)) {
    tmats <- term_matrices(model)
    event_indices <- attr(tmats, "event_term_indices")
    if (is.null(event_indices) || length(event_indices) == 0L) {
      stop(sprintf("%s engine requires at least one event/task regressor", engine), call. = FALSE)
    }
  }

  if (isTRUE(cfg$ar$by_cluster) && isTRUE(caps$requires_parcels_for_by_cluster)) {
    lowrank <- args$lowrank %||% list()
    if (is.null(lowrank$parcels)) {
      stop(sprintf("%s requires `lowrank$parcels` when ar_options$by_cluster = TRUE", engine), call. = FALSE)
    }
  }

  if (isTRUE(cfg$ar$by_cluster) && length(caps$forbid_by_cluster_dataset_classes) > 0L) {
    matches <- caps$forbid_by_cluster_dataset_classes[
      vapply(caps$forbid_by_cluster_dataset_classes, function(cls) inherits(dataset, cls), logical(1))
    ]
    if (length(matches) > 0L) {
      stop(
        sprintf("%s does not support by_cluster AR whitening for %s inputs", engine, matches[[1L]]),
        call. = FALSE
      )
    }
  }

  invisible(caps)
}

#' @keywords internal
.derive_engine_execution_config <- function(cfg, capabilities = NULL) {
  stopifnot(inherits(cfg, "fmri_lm_config"))
  caps <- .normalize_engine_capabilities(capabilities)

  executed <- unserialize(serialize(cfg, NULL))

  if (identical(caps$robust, FALSE)) {
    executed$robust <- fmri_lm_control()$robust
  }

  if (identical(caps$preprocessing, FALSE)) {
    defaults <- fmri_lm_control()
    executed$volume_weights <- defaults$volume_weights
    executed$soft_subspace <- defaults$soft_subspace
  }

  if (identical(caps$ar_voxelwise, FALSE)) {
    executed$ar$voxelwise <- FALSE
  }

  if (identical(caps$ar_by_cluster, FALSE) && "by_cluster" %in% names(executed$ar)) {
    executed$ar$by_cluster <- FALSE
  }

  class(executed) <- class(cfg)
  executed
}

#' Register an HRF/basis constructor
#'
#' @param name Character scalar used in formulas (e.g. `basis = "friman2"`).
#' @param constructor Function returning an object understood by `fmrihrf`
#'   (typically an `HRF`). Additional arguments from the original `hrf()` call
#'   are forwarded when the basis is constructed.
#' @return Invisibly, `TRUE`.
#' @export
register_basis <- function(name, constructor) {
  stopifnot(is.character(name), length(name) == 1L, nzchar(name))
  if (!is.function(constructor)) {
    stop("`constructor` must be a function", call. = FALSE)
  }
  .fmrireg_basis_registry[[name]] <- constructor
  invisible(TRUE)
}

#' @keywords internal
get_basis <- function(name) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    return(NULL)
  }
  .fmrireg_basis_registry[[name]]
}

#' Resolve a registered basis function by name
#'
#' This function is used internally by the formula processing system to lazily
#' resolve basis names (e.g., "cca2", "cca3", "spmg1") into actual basis objects.
#' It's part of the plugin API that allows extension packages to register custom
#' basis functions.
#'
#' @param name Character string naming a registered basis function
#' @param ... Additional arguments passed to the basis constructor
#' @return An HRF basis object
#' @export
#' @examples
#' \dontrun{
#' # Resolve the cca3 basis with specific parameters
#' basis <- resolve_basis("cca3", span = 30, TR = 2)
#' }
resolve_basis <- function(name, ...) {
  ctor <- get_basis(name)
  if (is.null(ctor)) {
    stop(sprintf("Unknown basis '%s'. Did you forget to register it?", name), call. = FALSE)
  }
  ctor(...)
}

#' @keywords internal
.fmrireg_resolve_basis <- function(name, ...) {
  # Deprecated: kept for backward compatibility
  resolve_basis(name, ...)
}

.fmrireg_replace_basis_expr <- function(expr, env) {
  if (!is.call(expr)) {
    return(expr)
  }
  parts <- as.list(expr)
  if (length(parts) == 0L) {
    return(expr)
  }
  if (length(parts) > 1L) {
    for (i in seq(2L, length(parts))) {
      parts[[i]] <- .fmrireg_replace_basis_expr(parts[[i]], env)
    }
  }
  call_name <- as.character(parts[[1L]])[1L]
  if (!nzchar(call_name) || !(call_name %in% c("hrf", "trialwise"))) {
    return(as.call(parts))
  }
  arg_names <- names(parts)
  if (is.null(arg_names)) {
    arg_names <- rep("", length(parts))
  }
  basis_idx <- which(arg_names == "basis")
  if (length(basis_idx) != 1L) {
    return(as.call(parts))
  }
  basis_expr <- parts[[basis_idx]]
  basis_name <- NULL
  if (is.character(basis_expr) && length(basis_expr) == 1L) {
    basis_name <- basis_expr
  } else {
    basis_val <- try(eval(basis_expr, envir = env), silent = TRUE)
    if (!inherits(basis_val, "try-error") && is.character(basis_val) && length(basis_val) == 1L) {
      basis_name <- basis_val
    } else {
      return(as.call(parts))
    }
  }
  if (is.null(get_basis(basis_name))) {
    return(as.call(parts))
  }
  extra_idx <- which(seq_along(parts) != 1L & seq_along(parts) != basis_idx & nzchar(arg_names))
  extra_args <- parts[extra_idx]
  extra_names <- arg_names[extra_idx]
  if (length(extra_args) > 0L) {
    names(extra_args) <- extra_names
  }
  # Create fmrireg::resolve_basis call for proper namespace access
  resolve_call <- call("::", quote(fmrireg), quote(resolve_basis))
  basis_call_parts <- c(list(resolve_call), list(basis_name), extra_args)
  names(basis_call_parts) <- c("", "", extra_names)
  parts[[basis_idx]] <- as.call(basis_call_parts)
  if (length(extra_idx) > 0L) {
    for (idx in sort(extra_idx, decreasing = TRUE)) {
      parts[[idx]] <- NULL
    }
  }
  as.call(parts)
}

#' @keywords internal
.fmrireg_inject_registered_bases <- function(formula) {
  if (!inherits(formula, "formula")) {
    return(formula)
  }
  env <- environment(formula)
  out <- formula
  rhs_idx <- length(out)
  env_eval <- env
  if (is.null(env_eval)) {
    env_eval <- globalenv()
  }
  out[[rhs_idx]] <- .fmrireg_replace_basis_expr(out[[rhs_idx]], env_eval)
  if (!is.null(env)) {
    environment(out) <- env
  }
  out
}

#' @keywords internal
prepare_fmri_lm_contrasts <- function(fmrimod) {
  processed_conlist <- list()
  contrast_info <- contrast_weights(fmrimod$event_model)
  col_indices <- attr(fmrimod$event_model$design_matrix, "col_indices")

  if (length(contrast_info) > 0 && !is.null(col_indices)) {
    for (contrast_name in names(contrast_info)) {
      con_spec <- contrast_info[[contrast_name]]

      term_name <- trimws(contrast_name)
      if (grepl("#", term_name, fixed = TRUE)) {
        term_name <- sub("#.*$", "", term_name)
      } else {
        term_name <- sub("\\..*$", "", term_name)
      }

      colind <- col_indices[[term_name]]
      if (is.null(colind)) {
        warning(sprintf("Contrast '%s' refers to term '%s' but col_indices are missing.", contrast_name, term_name))
        next
      }
      if (length(colind) == 0L) {
        warning(sprintf("No column indices found for term '%s'.", term_name))
        next
      }

      if (inherits(con_spec, "contrast") || inherits(con_spec, "Fcontrast")) {
        attr(con_spec$weights, "colind") <- colind
        attr(con_spec, "colind") <- colind
        processed_conlist[[contrast_name]] <- con_spec
      } else {
        warning(sprintf("Item '%s' is not a contrast or Fcontrast object.", contrast_name))
      }
    }
  }

  simple_conlist <- Filter(function(x) inherits(x, "contrast"), processed_conlist)
  fconlist <- Filter(function(x) inherits(x, "Fcontrast"), processed_conlist)

  list(
    processed = processed_conlist,
    simple = simple_conlist,
    f = fconlist,
    standard = c(simple_conlist, fconlist)
  )
}

#' Fit the core GLM on a transformed time-series matrix
#'
#' @param model An `fmri_model` describing the design.
#' @param Y Numeric matrix with `nrow(Y)` time points and columns matching voxels.
#' @param cfg Optional `fmri_lm_config`; defaults to `fmri_lm_control()`.
#' @param dataset Optional dataset backing the model. Defaults to `model$dataset` when available.
#' @param strategy Character label recorded on the returned object. Defaults to "external".
#' @param engine Character label indicating the engine that produced the fit. Defaults to "external".
#' @return An object of class `fmri_lm`.
#' @export
fit_glm_on_transformed_series <- function(model, Y, cfg = NULL, dataset = NULL,
                                          strategy = "external", engine = "external") {
  if (!inherits(model, "fmri_model")) {
    stop("`model` must be an 'fmri_model'", call. = FALSE)
  }
  if (is.null(cfg)) {
    cfg <- fmri_lm_control()
  } else if (!inherits(cfg, "fmri_lm_config")) {
    stop("`cfg` must inherit from 'fmri_lm_config'", call. = FALSE)
  }
  if (!is.matrix(Y)) {
    Y <- as.matrix(Y)
  }
  mode(Y) <- "numeric"

  X <- as.matrix(design_matrix(model))
  if (nrow(X) != nrow(Y)) {
    stop("Row mismatch between design matrix and response matrix", call. = FALSE)
  }

  dataset <- dataset %||% model$dataset

  contrast_prep <- prepare_fmri_lm_contrasts(model)
  processed_conlist <- contrast_prep$processed

  tmats <- term_matrices(model)
  event_indices <- attr(tmats, "event_term_indices")
  baseline_indices <- attr(tmats, "baseline_term_indices")
  vnames <- attr(tmats, "varnames")

  proj <- .fast_preproject(X)
  ctx <- glm_context(X = X, Y = Y, proj = proj)
  fit <- solve_glm_core(ctx)

  sigma_vec <- sqrt(fit$sigma2)
  dfres <- proj$dfres

  bstats <- beta_stats_matrix(fit$betas, proj$XtXinv, sigma_vec, dfres, vnames)

  simple_conlist <- contrast_prep$simple
  fconlist <- contrast_prep$f
  simple_weights <- lapply(simple_conlist, `[[`, "weights")
  if (length(simple_weights)) {
    names(simple_weights) <- names(simple_conlist)
  }
  fcon_weights <- lapply(fconlist, `[[`, "weights")
  if (length(fcon_weights)) {
    names(fcon_weights) <- names(fconlist)
  }

  contrast_results <- fit_lm_contrasts_fast(
    fit$betas,
    fit$sigma2,
    proj$XtXinv,
    simple_weights,
    fcon_weights,
    dfres
  )

  combined_contrasts <- if (length(contrast_results) > 0L) {
    dplyr::bind_rows(contrast_results)
  } else {
    tibble::tibble()
  }

  result <- list(
    contrasts = combined_contrasts,
    betas = bstats,
    event_indices = event_indices,
    baseline_indices = baseline_indices,
    cov.unscaled = proj$XtXinv,
    sigma = sigma_vec,
    rss = fit$rss,
    rdf = dfres,
    resvar = fit$rss / dfres,
    ar_coef = NULL
  )

  ret <- list(
    result = result,
    model = model,
    strategy = strategy,
    bcons = processed_conlist,
    dataset = dataset,
    ar_coef = result$ar_coef
  )
  class(ret) <- "fmri_lm"
  ret <- .fmri_lm_attach_config_metadata(ret, requested_cfg = cfg, executed_cfg = cfg)
  attr(ret, "strategy") <- strategy
  attr(ret, "engine") <- engine
  ret
}

#' Fit GLM with full config (AR/robust) on a transformed series
#'
#' Runs the same integrated solver used by `fmri_lm`, honoring AR/robust options
#' from `cfg`, but on an externally provided response matrix `Y` (T×V). This is
#' intended for engines that transform the time-series before inference.
#'
#' @param model An `fmri_model` describing the design.
#' @param Y Numeric matrix with `nrow(Y)` time points and columns matching voxels.
#' @param cfg Optional `fmri_lm_config`; defaults to `fmri_lm_control()`.
#' @param dataset Optional dataset backing the model. Defaults to `model$dataset` when available.
#' @param strategy Character label recorded on the returned object. Defaults to "external".
#' @param engine Character label indicating the engine that produced the fit. Defaults to "external".
#' @return An object of class `fmri_lm`.
#' @export
fit_glm_with_config <- function(model, Y, cfg = NULL, dataset = NULL,
                                strategy = "external", engine = "external") {
  if (!inherits(model, "fmri_model")) stop("`model` must be an 'fmri_model'", call. = FALSE)
  if (is.null(cfg)) cfg <- fmri_lm_control()
  if (!inherits(cfg, "fmri_lm_config")) stop("`cfg` must inherit from 'fmri_lm_config'", call. = FALSE)
  if (!is.matrix(Y)) Y <- as.matrix(Y)
  mode(Y) <- "numeric"

  X <- as.matrix(design_matrix(model))
  if (nrow(X) != nrow(Y)) stop("Row mismatch between design matrix and response matrix", call. = FALSE)
  dataset <- dataset %||% model$dataset

  contrast_prep <- prepare_fmri_lm_contrasts(model)
  processed_conlist <- contrast_prep$processed

  tmats <- term_matrices(model)
  event_indices <- attr(tmats, "event_term_indices")
  baseline_indices <- attr(tmats, "baseline_term_indices")
  vnames <- attr(tmats, "varnames")

  # Derive run indices from sampling frame when available
  run_indices <- NULL
  sframe <- tryCatch(model$event_model$sampling_frame, error = function(e) NULL)
  if (!is.null(sframe)) {
    runs <- tryCatch(fmrihrf::blockids(sframe), error = function(e) NULL)
    if (is.numeric(runs) && length(runs) == nrow(X)) {
      split_idx <- split(seq_len(nrow(X)), runs)
      run_indices <- lapply(split_idx, as.integer)
    }
  }

  fit <- solve_integrated_glm(X, Y, cfg, run_indices = run_indices)

  sigma_vec <- sqrt(fit$sigma2)
  dfres <- fit$dfres %||% (nrow(X) - ncol(X))
  XtXinv <- fit$XtXinv %||% .fast_preproject(X)$XtXinv

  bstats <- beta_stats_matrix(fit$betas, XtXinv, sigma_vec, dfres, vnames)

  simple_conlist <- contrast_prep$simple
  fconlist <- contrast_prep$f
  simple_weights <- lapply(simple_conlist, `[[`, "weights")
  if (length(simple_weights)) names(simple_weights) <- names(simple_conlist)
  fcon_weights <- lapply(fconlist, `[[`, "weights")
  if (length(fcon_weights)) names(fcon_weights) <- names(fconlist)

  contrast_results <- fit_lm_contrasts_fast(
    fit$betas,
    fit$sigma2,
    XtXinv,
    simple_weights,
    fcon_weights,
    dfres
  )

  combined_contrasts <- if (length(contrast_results) > 0L) dplyr::bind_rows(contrast_results) else tibble::tibble()

  result <- list(
    contrasts = combined_contrasts,
    betas = bstats,
    event_indices = event_indices,
    baseline_indices = baseline_indices,
    cov.unscaled = XtXinv,
    sigma = sigma_vec,
    rss = fit$rss %||% (fit$sigma2 * dfres),
    rdf = dfres,
    resvar = (fit$rss %||% (fit$sigma2 * dfres)) / dfres,
    ar_coef = fit$ar_coef %||% NULL
  )

  ret <- list(
    result = result,
    model = model,
    strategy = strategy,
    bcons = processed_conlist,
    dataset = dataset,
    ar_coef = result$ar_coef
  )
  class(ret) <- "fmri_lm"
  ret <- .fmri_lm_attach_config_metadata(ret, requested_cfg = cfg, executed_cfg = cfg)
  attr(ret, "strategy") <- strategy
  attr(ret, "engine") <- engine
  ret
}

#' Fit GLM from sufficient statistics
#'
#' Computes OLS/GLS-equivalent estimates from cross-products without materializing
#' the transformed series. Engines can stream `XtX`, `XtS`, and `StS` and call this
#' to obtain a standard `fmri_lm` object. AR/robust are not estimated here; this is
#' a low-level helper for OLS-equivalent inference from suffstats.
#'
#' @param model An `fmri_model` describing the design.
#' @param XtX p×p cross-product of the design matrix.
#' @param XtS p×V cross-product of design with data.
#' @param StS length-V vector of sum of squares per voxel.
#' @param df Residual degrees of freedom.
#' @param cfg Optional `fmri_lm_config`; used for metadata only.
#' @param dataset Optional dataset backing the model.
#' @param strategy Character label for the returned object.
#' @param engine Character label for the returned object.
#' @return An object of class `fmri_lm`.
#' @export
fit_glm_from_suffstats <- function(model, XtX, XtS, StS, df,
                                   cfg = NULL, dataset = NULL,
                                   strategy = "external", engine = "external") {
  if (!inherits(model, "fmri_model")) stop("`model` must be an 'fmri_model'", call. = FALSE)
  if (is.null(cfg)) cfg <- fmri_lm_control()
  if (!inherits(cfg, "fmri_lm_config")) stop("`cfg` must inherit from 'fmri_lm_config'", call. = FALSE)

  XtX <- as.matrix(XtX); XtS <- as.matrix(XtS); StS <- as.numeric(StS)
  p <- nrow(XtX); V <- ncol(XtS)
  stopifnot(ncol(XtX) == p, nrow(XtS) == p, length(StS) == V)

  # Invert XtX robustly
  XtXinv <- tryCatch(chol2inv(chol(XtX)), error = function(e) solve(XtX))
  betas <- XtXinv %*% XtS
  SSE <- StS - colSums(betas * XtS)
  sigma2 <- pmax(SSE / df, .Machine$double.eps)
  sigma_vec <- sqrt(sigma2)

  contrast_prep <- prepare_fmri_lm_contrasts(model)
  processed_conlist <- contrast_prep$processed
  tmats <- term_matrices(model)
  event_indices <- attr(tmats, "event_term_indices")
  baseline_indices <- attr(tmats, "baseline_term_indices")
  vnames <- attr(tmats, "varnames")

  bstats <- beta_stats_matrix(betas, XtXinv, sigma_vec, df, vnames)
  simple_conlist <- contrast_prep$simple
  fconlist <- contrast_prep$f
  simple_weights <- lapply(simple_conlist, `[[`, "weights"); if (length(simple_weights)) names(simple_weights) <- names(simple_conlist)
  fcon_weights <- lapply(fconlist, `[[`, "weights"); if (length(fcon_weights)) names(fcon_weights) <- names(fconlist)
  contrast_results <- fit_lm_contrasts_fast(betas, sigma2, XtXinv, simple_weights, fcon_weights, df)
  combined_contrasts <- if (length(contrast_results) > 0L) dplyr::bind_rows(contrast_results) else tibble::tibble()

  result <- list(
    contrasts = combined_contrasts,
    betas = bstats,
    event_indices = event_indices,
    baseline_indices = baseline_indices,
    cov.unscaled = XtXinv,
    sigma = sigma_vec,
    rss = SSE,
    rdf = df,
    resvar = SSE / df,
    ar_coef = NULL
  )
  dataset <- dataset %||% model$dataset
  ret <- list(result = result, model = model, strategy = strategy, bcons = processed_conlist, dataset = dataset, ar_coef = result$ar_coef)
  class(ret) <- "fmri_lm"
  ret <- .fmri_lm_attach_config_metadata(ret, requested_cfg = cfg, executed_cfg = cfg)
  attr(ret, "strategy") <- strategy
  attr(ret, "engine") <- engine
  ret
}
