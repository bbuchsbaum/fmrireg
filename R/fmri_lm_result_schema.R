# Versioned result and inference schema for fmri_lm().

#' The canonical empty contrast table
#'
#' Single source of truth for "this fit has no contrasts". Every estimation
#' path must use it so that zero requested contrasts stays zero rows: a
#' synthesised placeholder row is indistinguishable from a real contrast to
#' accessors ([coef()], [tidy()], `coef_names()`) and to [write_results()],
#' which already short-circuits correctly on an empty table.
#'
#' The columns are the ones every contrast table carries, so `filter(ct, type
#' == "contrast")`, `ct$name` and `bind_rows()` work on an empty table -- a
#' zero-column `tibble()` fails those with "object 'type' not found". It is
#' deliberately not a full schema: `package_tcontrast_result()` adds
#' `df.residual` and `unpack_chunkwise()` adds `contrast_internal_name`, so a
#' populated table from those paths is wider. Consumers must therefore still
#' branch on `nrow()`, never on `names()`.
#'
#' @keywords internal
#' @noRd
empty_contrast_table <- function() {
  dplyr::tibble(
    type = character(0),
    name = character(0),
    stat_type = character(0),
    conmat = list(),
    colind = list(),
    data = list()
  )
}

.fmri_lm_beta_payload <- function(result) {
  if (is.null(result$betas) || !nrow(result$betas)) return(NULL)
  data <- result$betas$data[[1L]]
  if (is.null(data) || !nrow(data)) return(NULL)
  data
}

.fmri_lm_result_nvox <- function(result) {
  payload <- .fmri_lm_beta_payload(result)
  if (!is.null(payload) && length(payload$estimate)) {
    estimate <- payload$estimate[[1L]]
    return(if (is.matrix(estimate)) nrow(estimate) else length(estimate))
  }
  length(result$sigma %||% result$resvar %||% result$rss)
}

.fmri_lm_nominal_df <- function(result) {
  df <- result$rdf
  if (is.null(df) && !is.null(result$betas) && nrow(result$betas)) {
    df <- result$betas$df.residual[[1L]]
  }
  as.numeric(df %||% NA_real_)
}

.fmri_lm_df_from_weights <- function(weights, nominal, nvox) {
  if (is.null(weights) || !is.finite(nominal)) return(rep(nominal, nvox))
  adjust <- function(w) {
    if (is.null(w)) return(nominal)
    w <- as.numeric(w)
    w <- w[is.finite(w)]
    if (!length(w)) return(1)
    max(nominal * sum(w) / length(w), 1)
  }
  if (is.numeric(weights)) return(rep(adjust(weights), nvox))
  if (is.list(weights) && length(weights) == nvox &&
      !any(vapply(weights, is.list, logical(1)))) {
    return(vapply(weights, adjust, numeric(1)))
  }
  # Runwise robust weights are nested by run and do not admit a single
  # observation-level weight vector after fixed-effect pooling. The individual
  # run statistics have already used their retained weights; report nominal
  # pooled df here rather than inventing a second adjustment.
  rep(nominal, nvox)
}

.new_fmri_variance_model <- function(method, covariance, covariance_scope,
                                     scale, df_nominal, df_inference,
                                     standard_errors, metadata = list()) {
  structure(list(
    method = method,
    covariance = covariance,
    covariance_scope = covariance_scope,
    scale = scale,
    df_nominal = df_nominal,
    df_inference = df_inference,
    standard_errors = standard_errors,
    metadata = metadata
  ), class = "fmri_lm_variance_model")
}

.fmri_lm_finalize_result <- function(result, cfg, compute = NULL,
                                     model = NULL, dataset = NULL, strategy = NULL,
                                     engine = "builtin") {
  if (!is.list(result)) stop("An fmri_lm result payload must be a list.", call. = FALSE)
  if (!inherits(cfg, "fmri_lm_control")) {
    stop("Result finalization requires an `fmri_lm_control`.", call. = FALSE)
  }

  nvox <- .fmri_lm_result_nvox(result)
  nominal <- .fmri_lm_nominal_df(result)
  df_inference <- .fmri_lm_df_from_weights(
    result$robust_weights, nominal = nominal, nvox = nvox
  )
  if (inherits(result$variance_model, "fmri_lm_variance_model") &&
      length(result$variance_model$df_inference)) {
    df_inference <- as.numeric(result$variance_model$df_inference)
  }
  status <- result$voxel_status
  if (is.null(status) || length(status) != nvox) status <- rep("ok", nvox)

  covariance <- result$covariance_model_basis %||%
    result$covariance_by_voxel %||%
    result$covariance_by_run %||%
    result$cov.unscaled
  covariance_scope <- if (!is.null(result$covariance_by_voxel)) {
    "voxel"
  } else if (!is.null(result$covariance_by_run)) {
    "run"
  } else if (!is.null(covariance)) {
    "shared"
  } else {
    "summary"
  }
  payload <- .fmri_lm_beta_payload(result)
  standard_errors <- if (!is.null(payload) && length(payload$se)) payload$se[[1L]] else NULL
  scale <- result$resvar %||% if (!is.null(result$sigma)) result$sigma^2 else NULL
  if (is.null(scale) && !is.null(payload) && length(payload$sigma)) {
    scale <- as.numeric(payload$sigma[[1L]])^2
  }
  if (is.null(result$rdf)) result$rdf <- nominal
  if (is.null(result$resvar) && !is.null(scale)) result$resvar <- scale
  if (is.null(result$sigma) && !is.null(scale)) result$sigma <- sqrt(pmax(scale, 0))

  needs_context <- !identical(cfg$variance$method, "model") ||
    identical(cfg$variance$df, "satterthwaite")
  variance_fit <- NULL
  if (needs_context) {
    if (!identical(cfg$estimation$scope, "joint")) {
      stop("HAC, sandwich, and Satterthwaite inference currently require ",
           "`estimation_spec(scope = 'joint')`; runwise meta-estimation must ",
           "combine run-level covariance explicitly.", call. = FALSE)
    }
    context <- result$inference_context
    if (is.null(context)) {
      can_reconstruct <- !is.null(model) && !is.null(dataset) &&
        identical(cfg$noise$struct, "iid") &&
        !.fmri_lm_robust_enabled(cfg$robust) &&
        !isTRUE(cfg$weights$enabled) && !isTRUE(cfg$projection$enabled)
      if (!can_reconstruct) {
        stop("The selected variance method requires the fitted-domain design ",
             "and residuals, but this engine did not return an inference context.",
             call. = FALSE)
      }
      X <- design_matrix(model)
      Y <- as.matrix(fmridataset::get_data_matrix(dataset))
      estimate <- as.matrix(payload$estimate[[1L]])
      run_chunks <- collect_chunks(exec_strategy("runwise")(dataset))
      context <- list(
        X = X,
        residuals = Y - X %*% t(estimate),
        runs = lapply(run_chunks, `[[`, "row_ind"),
        censor = NULL
      )
    }
    variance_fit <- .fmri_lm_variance_from_context(context, cfg$variance)
    df_inference <- variance_fit$df

    full_covariance <- NULL
    if (!identical(cfg$variance$method, "model")) {
      full_covariance <- variance_fit$covariance
      result <- .fmri_lm_update_inference(result, full_covariance, df_inference)
      covariance <- full_covariance
      covariance_scope <- "voxel"
      scale <- NULL
    } else if (identical(covariance_scope, "shared") && !is.null(scale)) {
      if (length(scale) == 1L) scale <- rep(scale, nvox)
      full_covariance <- lapply(seq_len(nvox), function(v) covariance * scale[v])
      result <- .fmri_lm_update_inference(result, full_covariance, df_inference)
    }
    result$inference_context <- NULL
  }

  # Inference may have been regenerated above; capture the authoritative
  # standard errors after that update rather than retaining the model-path
  # values that entered the finalizer.
  payload <- .fmri_lm_beta_payload(result)
  standard_errors <- if (!is.null(payload) && length(payload$se)) payload$se[[1L]] else NULL

  result$schema_version <- 2L
  result$voxel_status <- status
  result$df <- list(nominal = nominal, inference = df_inference,
                    method = cfg$variance$df)
  result$fit_state <- list(
    robust_weights = result$robust_weights,
    ar_parameters = result$ar_coef,
    ma_parameters = result$ma_coef,
    temporal_diagnostics = result$temporal_diagnostics,
    voxel_status = status
  )
  if (is.null(result$variance_model)) {
    result$variance_model <- .new_fmri_variance_model(
      method = cfg$variance$method,
      covariance = covariance,
      covariance_scope = covariance_scope,
      scale = scale,
      df_nominal = nominal,
      df_inference = df_inference,
      standard_errors = standard_errors,
      metadata = list(
        estimator = engine,
        estimation_scope = cfg$estimation$scope,
        noise = cfg$noise,
        robust = cfg$robust,
        df_method = cfg$variance$df,
        taper = cfg$variance$taper,
        compute = compute,
        strategy = strategy,
        model_based_factorization = identical(cfg$variance$method, "model") &&
          covariance_scope %in% c("shared", "voxel"),
        selected_max_lag = variance_fit$max_lag %||% NULL
      )
    )
  }
  result
}

#' Extract the structured variance model
#'
#' @param x A fitted `fmri_lm` object.
#' @return An `fmri_lm_variance_model` object.
#' @export
variance_model <- function(x) {
  if (!inherits(x, "fmri_lm")) stop("`x` must be an fmri_lm fit.", call. = FALSE)
  x$result$variance_model
}

#' @export
print.fmri_lm_variance_model <- function(x, ...) {
  cat("<fmri_lm_variance_model>\n")
  cat("  method:", x$method, "\n")
  cat("  covariance scope:", x$covariance_scope, "\n")
  cat("  nominal df:", paste(unique(x$df_nominal), collapse = ", "), "\n")
  cat("  inference df:", paste(utils::head(unique(x$df_inference), 5L), collapse = ", "), "\n")
  invisible(x)
}
