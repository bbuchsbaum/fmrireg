# Smooth condition-level HRF estimation

.validate_estimate_hrf_formula <- function(form) {
  if (!inherits(form, "formula") || length(form) != 3L) {
    stop("form must be a two-sided event-model formula", call. = FALSE)
  }
  rhs_names <- all.names(form[[3L]], functions = TRUE, unique = TRUE)
  if ("trialwise" %in% rhs_names) {
    stop(
      "estimate_hrf() estimates condition-level hrf() curves; trialwise() terms are not supported.",
      call. = FALSE
    )
  }
  if (!("hrf" %in% rhs_names)) {
    stop("form must contain at least one condition-level hrf() term", call. = FALSE)
  }
  invisible(TRUE)
}

.order_hrf_event_data <- function(data, form, block) {
  onset_variables <- all.vars(form[[2L]])
  if (length(onset_variables) != 1L || !(onset_variables %in% names(data))) {
    stop("The left-hand side of form must name one onset column in the event table", call. = FALSE)
  }
  if (!inherits(block, "formula") || length(block) != 2L) {
    stop("block must be a one-sided formula", call. = FALSE)
  }
  block_frame <- stats::model.frame(
    block,
    data = data,
    na.action = stats::na.pass,
    drop.unused.levels = FALSE
  )
  ordering <- do.call(
    order,
    c(as.list(block_frame), list(data[[onset_variables]], na.last = TRUE))
  )
  data[ordering, , drop = FALSE]
}

.new_hrf_basis_spec <- function(type, k, span) {
  degree <- switch(type, bspline = 3L, tent = 1L)
  full_df <- k + 2L
  n_internal <- full_df - degree - 1L
  knots <- if (n_internal > 0L) {
    seq(0, span, length.out = n_internal + 2L)[-c(1L, n_internal + 2L)]
  } else {
    numeric()
  }
  list(type = type, k = k, span = span, degree = degree, knots = knots)
}

.evaluate_hrf_basis <- function(time, spec) {
  time <- as.numeric(time)
  out <- matrix(0, nrow = length(time), ncol = spec$k)
  inside <- is.finite(time) & time >= 0 & time <= spec$span
  if (any(inside)) {
    full <- splines::bs(
      time[inside],
      knots = spec$knots,
      degree = spec$degree,
      intercept = TRUE,
      Boundary.knots = c(0, spec$span)
    )
    if (ncol(full) != spec$k + 2L) {
      stop("Internal HRF basis dimension mismatch", call. = FALSE)
    }
    out[inside, ] <- full[, -c(1L, ncol(full)), drop = FALSE]
  }
  colnames(out) <- sprintf("basis_%02d", seq_len(spec$k))
  out
}

.as_estimation_hrf_basis <- function(spec) {
  basis_spec <- spec
  basis_fun <- function(t) .evaluate_hrf_basis(t, basis_spec)
  fmrihrf::as_hrf(
    basis_fun,
    name = paste0("estimate_hrf_", spec$type),
    nbasis = spec$k,
    span = spec$span,
    params = list(.estimator = "estimate_hrf")
  )
}

.hrf_curve_map <- function(event_model, design, k) {
  design_names <- colnames(design)
  if (is.null(design_names) || any(!nzchar(design_names))) {
    stop("The event design must have non-empty column names", call. = FALSE)
  }

  colmap <- tryCatch(
    as.data.frame(fmridesign::design_colmap(event_model)),
    error = function(e) NULL
  )
  has_semantic_map <- !is.null(colmap) &&
    all(c("col", "term_tag", "condition", "basis_ix") %in% names(colmap)) &&
    nrow(colmap) == ncol(design) &&
    setequal(as.integer(colmap$col), seq_len(ncol(design)))

  if (has_semantic_map) {
    colmap <- colmap[order(as.integer(colmap$col)), , drop = FALSE]
    term_by_column <- as.character(colmap$term_tag)
    condition_by_column <- as.character(colmap$condition)
    curve_id <- paste(term_by_column, condition_by_column, sep = "_")
    basis_index <- as.integer(colmap$basis_ix)
  } else {
    suffix <- "_b([0-9]+)$"
    matched <- grepl(suffix, design_names)
    if (!all(matched)) {
      stop(
        "Could not map all event-design columns to HRF basis components",
        call. = FALSE
      )
    }
    curve_id <- sub(suffix, "", design_names)
    basis_index <- as.integer(sub("^.*_b([0-9]+)$", "\\1", design_names))
    term_by_column <- sub("_.*$", "", curve_id)
    condition_by_column <- ifelse(
      startsWith(curve_id, paste0(term_by_column, "_")),
      substring(curve_id, nchar(term_by_column) + 2L),
      curve_id
    )
  }
  curves <- unique(curve_id)
  raw_indices <- lapply(curves, function(curve) which(curve_id == curve))

  for (i in seq_along(raw_indices)) {
    observed <- sort(basis_index[raw_indices[[i]]])
    if (!identical(observed, seq_len(k))) {
      stop(
        sprintf(
          "Curve '%s' does not contain exactly basis components 1 through %d",
          curves[[i]], k
        ),
        call. = FALSE
      )
    }
    raw_indices[[i]] <- raw_indices[[i]][order(basis_index[raw_indices[[i]]])]
  }

  order_index <- unlist(raw_indices, use.names = FALSE)
  reordered <- design[, order_index, drop = FALSE]
  colnames(reordered) <- unlist(lapply(seq_along(curves), function(i) {
    sprintf("%s_b%02d", curves[[i]], seq_len(k))
  }), use.names = FALSE)
  ordered_indices <- split(
    seq_len(ncol(reordered)),
    rep(seq_along(curves), each = k)
  )

  term <- character(length(curves))
  condition <- character(length(curves))
  for (i in seq_along(curves)) {
    first_column <- raw_indices[[i]][1L]
    term[[i]] <- term_by_column[[first_column]]
    condition[[i]] <- condition_by_column[[first_column]]
  }
  display_curve <- condition
  duplicated_condition <- duplicated(display_curve) |
    duplicated(display_curve, fromLast = TRUE)
  display_curve[duplicated_condition] <- paste(
    term[duplicated_condition],
    condition[duplicated_condition],
    sep = ":"
  )
  display_curve <- make.unique(display_curve, sep = "#")

  list(
    design = reordered,
    indices = unname(ordered_indices),
    info = tibble::tibble(curve = display_curve, term = term, condition = condition)
  )
}

.prepare_hrf_estimation <- function(form, fixed, block, dataset, basemod, basis_spec) {
  hrf_basis <- .as_estimation_hrf_basis(basis_spec)
  estimation_formula <- inject_basis(form, hrf_basis, fun_names = "hrf")
  event_data <- .order_hrf_event_data(dataset$event_table, form, block)
  event_mod <- event_model(
    estimation_formula,
    data = event_data,
    block = block,
    sampling_frame = dataset$sampling_frame
  )
  event_design_raw <- as.matrix(design_matrix(event_mod))
  curve_map <- .hrf_curve_map(event_mod, event_design_raw, basis_spec$k)

  baseline_mod <- if (is.null(basemod)) {
    baseline_model("constant", sframe = dataset$sampling_frame)
  } else {
    basemod
  }
  baseline_design <- as.matrix(design_matrix(baseline_mod))

  fixed_mod <- NULL
  fixed_design <- NULL
  if (!is.null(fixed)) {
    if (!inherits(fixed, "formula")) {
      stop("fixed must be NULL or an event-model formula", call. = FALSE)
    }
    fixed_data <- .order_hrf_event_data(dataset$event_table, fixed, block)
    fixed_mod <- event_model(
      fixed,
      data = fixed_data,
      block = block,
      sampling_frame = dataset$sampling_frame
    )
    fixed_design <- as.matrix(design_matrix(fixed_mod))
  }

  nuisance_design <- if (is.null(fixed_design)) {
    baseline_design
  } else {
    cbind(baseline_design, fixed_design)
  }

  list(
    event_design = curve_map$design,
    nuisance_design = nuisance_design,
    curve_indices = curve_map$indices,
    curve_info = curve_map$info,
    event_model = event_mod,
    fixed_model = fixed_mod,
    baseline_model = baseline_mod,
    formula = estimation_formula
  )
}

.partial_out_hrf_nuisance <- function(event_design, response, nuisance_design) {
  if (is.null(nuisance_design) || ncol(nuisance_design) == 0L) {
    return(list(
      event = event_design,
      response = response,
      nuisance_rank = 0L,
      nuisance_columns = 0L
    ))
  }

  nuisance_qr <- qr(nuisance_design, LAPACK = FALSE)
  nuisance_rank <- nuisance_qr$rank
  if (nuisance_rank < ncol(nuisance_design)) {
    warning(
      sprintf(
        "Nuisance design is rank deficient: retained %d of %d independent directions.",
        nuisance_rank, ncol(nuisance_design)
      ),
      call. = FALSE
    )
  }

  Q <- if (nuisance_rank > 0L) {
    qr.Q(nuisance_qr, complete = FALSE)[, seq_len(nuisance_rank), drop = FALSE]
  } else {
    matrix(0, nrow = nrow(event_design), ncol = 0L)
  }
  residualize <- function(x) {
    if (ncol(Q) == 0L) x else x - Q %*% crossprod(Q, x)
  }

  list(
    event = residualize(event_design),
    response = residualize(response),
    nuisance_rank = nuisance_rank,
    nuisance_columns = ncol(nuisance_design)
  )
}

.hrf_smoothness_penalty <- function(k, n_curves) {
  difference_order <- if (k >= 3L) 2L else 1L
  D <- diff(diag(k), differences = difference_order)
  block <- crossprod(D)
  kronecker(diag(n_curves), block)
}

.solve_hrf_system <- function(XtX, XtY, penalty, lambda) {
  A <- XtX + lambda * penalty
  factor <- tryCatch(chol(A), error = function(e) NULL)
  if (is.null(factor)) {
    return(NULL)
  }
  inverse <- chol2inv(factor)
  coefficients <- inverse %*% XtY
  edf <- sum(diag(inverse %*% XtX))
  list(coefficients = coefficients, inverse = inverse, edf = edf)
}

.select_hrf_lambda <- function(
    XtX, XtY, response_ss, penalty, n_effective,
    lambda, lambda_grid, progress) {
  if (is.character(lambda)) {
    if (length(lambda) != 1L || !identical(tolower(lambda), "gcv")) {
      stop("lambda must be a non-negative number or 'gcv'", call. = FALSE)
    }
    candidates <- sort(unique(as.numeric(lambda_grid)))
  } else {
    if (!is.numeric(lambda) || length(lambda) != 1L ||
        !is.finite(lambda) || lambda < 0) {
      stop("lambda must be a non-negative number or 'gcv'", call. = FALSE)
    }
    candidates <- as.numeric(lambda)
  }
  if (!length(candidates) || any(!is.finite(candidates)) || any(candidates < 0)) {
    stop("lambda_grid must contain finite non-negative values", call. = FALSE)
  }

  scores <- rep(Inf, length(candidates))
  edf <- rep(NA_real_, length(candidates))
  progress_id <- NULL
  if (isTRUE(progress) && length(candidates) > 1L) {
    progress_id <- cli::cli_progress_bar(
      "Selecting HRF smoothness",
      total = length(candidates)
    )
    on.exit(cli::cli_progress_done(id = progress_id), add = TRUE)
  }

  for (i in seq_along(candidates)) {
    solved <- .solve_hrf_system(XtX, XtY, penalty, candidates[[i]])
    if (!is.null(solved)) {
      beta <- solved$coefficients
      fitted_ss <- colSums(beta * (XtX %*% beta))
      cross_term <- colSums(beta * XtY)
      rss <- pmax(response_ss - 2 * cross_term + fitted_ss, 0)
      denominator <- 1 - solved$edf / n_effective
      if (is.finite(denominator) && denominator > 0) {
        relative_rss <- rss / pmax(response_ss, .Machine$double.eps)
        scores[[i]] <- mean(relative_rss) / denominator^2
        edf[[i]] <- solved$edf
      }
    }
    if (!is.null(progress_id)) cli::cli_progress_update(id = progress_id)
  }

  if (!any(is.finite(scores))) {
    stop("No smoothing candidate produced an identifiable HRF fit", call. = FALSE)
  }
  best_score <- min(scores, na.rm = TRUE)
  tied <- which(is.finite(scores) & scores <= best_score * (1 + 1e-10))
  best <- tied[[length(tied)]]

  list(
    lambda = candidates[[best]],
    table = tibble::tibble(lambda = candidates, score = scores, edf = edf)
  )
}

#' Estimate smooth condition-level hemodynamic response functions
#'
#' `estimate_hrf()` estimates condition-level HRF curves directly from an fMRI
#' dataset. It builds an event-aligned spline basis, removes baseline and fixed
#' nuisance effects once, and solves every voxel together with one penalized
#' multiresponse linear system. This replaces the former voxel-by-voxel GAM
#' implementation, whose predictor did not represent post-stimulus time.
#'
#' A single smoothing parameter is shared across voxels. With `lambda = "gcv"`,
#' it is selected by scale-normalized generalized cross-validation so high-
#' variance voxels do not dominate the choice. Curves are constrained to zero at
#' the beginning and end of `rsam`. Standard errors and confidence intervals use
#' the fitted penalized-linear-model covariance under independent, homoscedastic
#' time-point errors; they are not autocorrelation-robust.
#'
#' @param form A two-sided event-model formula containing one or more `hrf()`
#'   terms. `trialwise()` terms are rejected because this function estimates
#'   condition-level curves.
#' @param fixed Optional event-model formula whose design is treated as nuisance.
#' @param block Formula identifying acquisition runs or blocks.
#' @param dataset An `fmri_dataset`.
#' @param bs Deprecated legacy basis selector. Values from the former GAM API are
#'   mapped to `basis = "bspline"`.
#' @param rsam Strictly increasing, finite post-stimulus times beginning at zero.
#' @param basemod Optional baseline model. The default is a constant baseline.
#' @param k Number of free HRF basis functions per curve. At least four for the
#'   default cubic B-spline basis.
#' @param fx Deprecated legacy selector. `fx = TRUE` maps to `lambda = 0` when
#'   `lambda` is omitted.
#' @param progress Show progress while scanning the GCV grid.
#' @param basis Either `"bspline"` (cubic) or `"tent"` (piecewise linear).
#' @param lambda Non-negative smoothing strength or `"gcv"` for shared automatic
#'   selection.
#' @param lambda_grid Candidate smoothing strengths used for GCV.
#' @param ci_level Confidence level in `(0, 1)`, or `NULL` to omit intervals.
#'
#' @return An object of class `fmri_hrf_estimate`. Its `estimate` and
#'   `std.error` arrays have dimensions time by curve by voxel. The object also
#'   contains labeled curve metadata, basis coefficients, smoothing diagnostics,
#'   and the designs used for estimation. Use `tidy()`, `predict()`, `coef()`, or
#'   `as.matrix()` for common downstream representations.
#'
#' @seealso [fitted_hrf()], [event_model()], [baseline_model()]
#' @examples
#' set.seed(18)
#' n <- 80L
#' events <- data.frame(
#'   onset = seq(6, 62, by = 8),
#'   condition = factor(rep(c("A", "B"), 4)),
#'   run = 1L
#' )
#' dataset <- fmridataset::matrix_dataset(
#'   matrix(rnorm(n * 2), nrow = n),
#'   TR = 1,
#'   run_length = n,
#'   event_table = events
#' )
#' fit <- estimate_hrf(
#'   onset ~ hrf(condition),
#'   block = ~run,
#'   dataset = dataset,
#'   rsam = 0:12,
#'   k = 6,
#'   lambda = 1
#' )
#' tidy(fit, voxel = 1)
#' @export
estimate_hrf <- function(
    form,
    fixed = NULL,
    block,
    dataset,
    bs = NULL,
    rsam = seq(0, 20, by = 1),
    basemod = NULL,
    k = 8L,
    fx = NULL,
    progress = FALSE,
    basis = c("bspline", "tent"),
    lambda = "gcv",
    lambda_grid = c(0, 10^seq(-4, 4, length.out = 25L)),
    ci_level = 0.95) {
  lambda_missing <- missing(lambda)
  .validate_estimate_hrf_formula(form)
  if (!inherits(dataset, "fmri_dataset")) {
    stop("dataset must inherit from 'fmri_dataset'", call. = FALSE)
  }

  if (!is.null(bs)) {
    bs_value <- match.arg(as.character(bs)[[1L]], c("tp", "ts", "cr", "ps", "bspline", "tent"))
    warning("bs is deprecated; use basis instead", call. = FALSE)
    basis <- if (identical(bs_value, "tent")) "tent" else "bspline"
  } else {
    basis <- match.arg(basis)
  }
  if (!is.null(fx)) {
    if (!is.logical(fx) || length(fx) != 1L || is.na(fx)) {
      stop("fx must be NULL, TRUE, or FALSE", call. = FALSE)
    }
    warning("fx is deprecated; use lambda = 0 for an unpenalized fit", call. = FALSE)
    if (isTRUE(fx) && lambda_missing) lambda <- 0
  }

  rsam <- as.numeric(rsam)
  if (length(rsam) < 2L || any(!is.finite(rsam))) {
    stop("rsam must contain at least two finite values", call. = FALSE)
  }
  if (abs(rsam[[1L]]) > sqrt(.Machine$double.eps)) {
    stop("rsam must start at zero", call. = FALSE)
  }
  if (any(diff(rsam) <= 0)) {
    stop("rsam must be strictly increasing", call. = FALSE)
  }

  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k != as.integer(k)) {
    stop("k must be a single integer", call. = FALSE)
  }
  k <- as.integer(k)
  minimum_k <- if (identical(basis, "bspline")) 4L else 2L
  if (k < minimum_k) {
    stop(sprintf("k must be at least %d for basis = '%s'", minimum_k, basis), call. = FALSE)
  }
  if (!is.null(ci_level) &&
      (!is.numeric(ci_level) || length(ci_level) != 1L ||
       !is.finite(ci_level) || ci_level <= 0 || ci_level >= 1)) {
    stop("ci_level must be NULL or a number strictly between zero and one", call. = FALSE)
  }

  span <- max(rsam)
  basis_spec <- .new_hrf_basis_spec(basis, k, span)
  prepared <- .prepare_hrf_estimation(form, fixed, block, dataset, basemod, basis_spec)

  response <- as.matrix(fmridataset::get_data_matrix(dataset))
  if (nrow(response) != nrow(prepared$event_design)) {
    stop("Dataset and HRF event design have different numbers of time points", call. = FALSE)
  }
  if (any(!is.finite(response))) {
    stop("Dataset contains non-finite response values", call. = FALSE)
  }
  if (ncol(response) < 1L) {
    stop("Dataset contains no response columns", call. = FALSE)
  }

  partial <- .partial_out_hrf_nuisance(
    prepared$event_design,
    response,
    prepared$nuisance_design
  )
  X <- partial$event
  Y <- partial$response
  event_rank <- qr(X, LAPACK = FALSE)$rank
  if (event_rank < ncol(X)) {
    stop(
      sprintf(
        "Residualized HRF design is rank deficient: rank %d for %d coefficients; add events, reduce k, or simplify form.",
        event_rank, ncol(X)
      ),
      call. = FALSE
    )
  }

  n_effective <- nrow(X) - partial$nuisance_rank
  if (n_effective <= ncol(X)) {
    stop("Too few residual degrees of freedom for HRF estimation", call. = FALSE)
  }
  XtX <- crossprod(X)
  XtY <- crossprod(X, Y)
  response_ss <- colSums(Y^2)
  penalty <- .hrf_smoothness_penalty(k, nrow(prepared$curve_info))

  selected <- .select_hrf_lambda(
    XtX = XtX,
    XtY = XtY,
    response_ss = response_ss,
    penalty = penalty,
    n_effective = n_effective,
    lambda = lambda,
    lambda_grid = lambda_grid,
    progress = progress
  )
  solved <- .solve_hrf_system(XtX, XtY, penalty, selected$lambda)
  if (is.null(solved)) {
    stop("Selected smoothing strength did not produce an identifiable HRF fit", call. = FALSE)
  }

  coefficients <- solved$coefficients
  rownames(coefficients) <- colnames(prepared$event_design)
  voxel_names <- colnames(response)
  if (is.null(voxel_names) || any(!nzchar(voxel_names))) {
    voxel_names <- paste0("voxel_", seq_len(ncol(response)))
  }
  colnames(coefficients) <- voxel_names

  fitted_ss <- colSums(coefficients * (XtX %*% coefficients))
  cross_term <- colSums(coefficients * XtY)
  rss <- pmax(response_ss - 2 * cross_term + fitted_ss, 0)
  df_residual <- n_effective - solved$edf
  if (!is.finite(df_residual) || df_residual <= 0) {
    stop("The fitted HRF model has no residual degrees of freedom", call. = FALSE)
  }
  sigma2 <- rss / df_residual
  covariance_unscaled <- solved$inverse %*% XtX %*% solved$inverse

  basis_at_time <- .evaluate_hrf_basis(rsam, basis_spec)
  curve_names <- prepared$curve_info$curve
  estimate <- array(
    NA_real_,
    dim = c(length(rsam), length(curve_names), ncol(response)),
    dimnames = list(time = as.character(rsam), curve = curve_names, voxel = voxel_names)
  )
  std_error <- estimate
  for (curve in seq_along(prepared$curve_indices)) {
    index <- prepared$curve_indices[[curve]]
    estimate[, curve, ] <- basis_at_time %*% coefficients[index, , drop = FALSE]
    curve_covariance <- covariance_unscaled[index, index, drop = FALSE]
    variance_scale <- rowSums((basis_at_time %*% curve_covariance) * basis_at_time)
    std_error[, curve, ] <- sqrt(pmax(outer(variance_scale, sigma2), 0))
  }

  lower <- upper <- NULL
  if (!is.null(ci_level)) {
    critical <- stats::qt(1 - (1 - ci_level) / 2, df = df_residual)
    lower <- estimate - critical * std_error
    upper <- estimate + critical * std_error
  }

  structure(
    list(
      estimate = estimate,
      std.error = std_error,
      lower = lower,
      upper = upper,
      time = rsam,
      curves = curve_names,
      voxels = voxel_names,
      curve_info = prepared$curve_info,
      coefficients = coefficients,
      sigma = sqrt(sigma2),
      df.residual = df_residual,
      edf = solved$edf,
      lambda = selected$lambda,
      gcv = selected$table,
      basis = basis,
      basis_spec = basis_spec,
      basis_at_time = basis_at_time,
      span = span,
      ci_level = ci_level,
      covariance_unscaled = covariance_unscaled,
      event_model = prepared$event_model,
      fixed_model = prepared$fixed_model,
      baseline_model = prepared$baseline_model,
      event_design = prepared$event_design,
      nuisance_design = prepared$nuisance_design,
      formula = form,
      fixed = fixed,
      block = block,
      diagnostics = list(
        algorithm = "shared penalized multiresponse solve",
        n_time = nrow(response),
        n_voxels = ncol(response),
        n_curves = length(curve_names),
        event_rank = event_rank,
        event_columns = ncol(X),
        nuisance_rank = partial$nuisance_rank,
        nuisance_columns = partial$nuisance_columns,
        condition_number = kappa(XtX)
      ),
      call = match.call()
    ),
    class = "fmri_hrf_estimate"
  )
}

#' Methods for smooth HRF estimates
#'
#' @param x,object An `fmri_hrf_estimate` object.
#' @param curve A curve name or one-based curve index. Required by `as.matrix()`
#'   when the object contains multiple curves.
#' @param what Return estimated curves or their standard errors.
#' @param newdata Numeric post-stimulus times within the fitted span.
#' @param se.fit Return prediction standard errors with fitted curves.
#' @param ... Unused.
#' @name fmri_hrf_estimate-methods
NULL

#' @rdname fmri_hrf_estimate-methods
#' @export
print.fmri_hrf_estimate <- function(x, ...) {
  cat("<fmri_hrf_estimate>\n")
  cat(sprintf("  %d curve%s x %d voxel%s at %d time points\n",
              length(x$curves), if (length(x$curves) == 1L) "" else "s",
              length(x$voxels), if (length(x$voxels) == 1L) "" else "s",
              length(x$time)))
  cat(sprintf("  basis: %s (%d functions per curve)\n", x$basis, x$basis_spec$k))
  cat(sprintf("  lambda: %.6g; effective df: %.3f; residual df: %.3f\n",
              x$lambda, x$edf, x$df.residual))
  invisible(x)
}

#' @rdname fmri_hrf_estimate-methods
#' @export
coef.fmri_hrf_estimate <- function(object, ...) {
  object$coefficients
}

#' @rdname fmri_hrf_estimate-methods
#' @export
as.matrix.fmri_hrf_estimate <- function(x, curve = NULL,
                                        what = c("estimate", "std.error"), ...) {
  what <- match.arg(what)
  if (is.null(curve)) {
    if (length(x$curves) != 1L) {
      stop("curve must be supplied when more than one HRF curve was estimated", call. = FALSE)
    }
    curve <- x$curves[[1L]]
  }
  if (is.numeric(curve)) {
    if (length(curve) != 1L || is.na(curve) || curve < 1L || curve > length(x$curves)) {
      stop("curve index is out of bounds", call. = FALSE)
    }
    curve <- x$curves[[as.integer(curve)]]
  }
  if (!is.character(curve) || length(curve) != 1L || !(curve %in% x$curves)) {
    stop("curve must identify one fitted HRF curve", call. = FALSE)
  }
  values <- if (identical(what, "estimate")) x$estimate else x$std.error
  out <- values[, curve, , drop = FALSE]
  dim(out) <- c(length(x$time), length(x$voxels))
  dimnames(out) <- list(time = as.character(x$time), voxel = x$voxels)
  out
}

#' @rdname fmri_hrf_estimate-methods
#' @export
predict.fmri_hrf_estimate <- function(object, newdata = object$time,
                                      se.fit = FALSE, ...) {
  newdata <- as.numeric(newdata)
  if (!length(newdata) || any(!is.finite(newdata)) ||
      any(newdata < 0) || any(newdata > object$span)) {
    stop(sprintf("newdata must contain finite times within [0, %g]", object$span), call. = FALSE)
  }
  basis <- .evaluate_hrf_basis(newdata, object$basis_spec)
  fit <- array(
    NA_real_,
    dim = c(length(newdata), length(object$curves), length(object$voxels)),
    dimnames = list(time = as.character(newdata), curve = object$curves, voxel = object$voxels)
  )
  se <- fit
  sigma2 <- object$sigma^2
  for (curve in seq_along(object$curves)) {
    index <- seq.int((curve - 1L) * object$basis_spec$k + 1L,
                     curve * object$basis_spec$k)
    fit[, curve, ] <- basis %*% object$coefficients[index, , drop = FALSE]
    curve_covariance <- object$covariance_unscaled[index, index, drop = FALSE]
    variance_scale <- rowSums((basis %*% curve_covariance) * basis)
    se[, curve, ] <- sqrt(pmax(outer(variance_scale, sigma2), 0))
  }
  if (isTRUE(se.fit)) list(fit = fit, se.fit = se, df = object$df.residual) else fit
}

#' Tidy a smooth HRF estimate
#'
#' @param x An `fmri_hrf_estimate` object.
#' @param curve Optional curve names to retain.
#' @param voxel Optional voxel names or indices to retain.
#' @param ... Unused.
#' @return A tibble with one row per time, curve, and voxel.
#' @export
tidy.fmri_hrf_estimate <- function(x, curve = NULL, voxel = NULL, ...) {
  curve_index <- if (is.null(curve)) {
    seq_along(x$curves)
  } else if (is.numeric(curve)) {
    as.integer(curve)
  } else {
    match(curve, x$curves)
  }
  voxel_index <- if (is.null(voxel)) {
    seq_along(x$voxels)
  } else if (is.numeric(voxel)) {
    as.integer(voxel)
  } else {
    match(voxel, x$voxels)
  }
  if (anyNA(curve_index) || any(curve_index < 1L) || any(curve_index > length(x$curves))) {
    stop("curve contains unknown or out-of-range values", call. = FALSE)
  }
  if (anyNA(voxel_index) || any(voxel_index < 1L) || any(voxel_index > length(x$voxels))) {
    stop("voxel contains unknown or out-of-range values", call. = FALSE)
  }

  estimate <- x$estimate[, curve_index, voxel_index, drop = FALSE]
  std_error <- x$std.error[, curve_index, voxel_index, drop = FALSE]
  lower <- if (is.null(x$lower)) array(NA_real_, dim(estimate)) else x$lower[, curve_index, voxel_index, drop = FALSE]
  upper <- if (is.null(x$upper)) array(NA_real_, dim(estimate)) else x$upper[, curve_index, voxel_index, drop = FALSE]
  grid <- expand.grid(
    time = x$time,
    curve = x$curves[curve_index],
    voxel = x$voxels[voxel_index],
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  curve_rows <- match(grid$curve, x$curve_info$curve)

  tibble::tibble(
    time = as.numeric(grid$time),
    curve = as.character(grid$curve),
    term = x$curve_info$term[curve_rows],
    condition = x$curve_info$condition[curve_rows],
    voxel = as.character(grid$voxel),
    estimate = as.numeric(estimate),
    std.error = as.numeric(std_error),
    lower = as.numeric(lower),
    upper = as.numeric(upper)
  )
}
