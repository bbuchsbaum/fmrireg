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
#' @param capabilities Optional named list describing the engine (for future use).
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
  attr(ret, "config") <- cfg
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
  attr(ret, "config") <- cfg
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
  attr(ret, "config") <- cfg
  attr(ret, "strategy") <- strategy
  attr(ret, "engine") <- engine
  ret
}
