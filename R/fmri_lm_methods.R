# S3 Methods for fmri_lm Objects
# Methods for extracting results and information from fitted fmri_lm objects

#' Reshape Coefficients
#' 
#' @keywords internal
#' @noRd
reshape_coef <- function(df, des, measure = "value") {
  nrun <- length(levels(des$blockids))
  run_order <- order(des$blockids)
  ret <- matrix(t(as.matrix(df[, run_order, drop = FALSE])),
                nrow = nrun, byrow = TRUE)
  as.data.frame(ret)
}

## pull_stat_revised/pull_stat are defined in R/fmrilm.R
## Remove duplicate definitions here to ensure a single source of truth.

#' Extract Coefficients from an fmri_lm Fit
#'
#' Returns estimated coefficients from a fitted \code{fmri_lm} model. Every
#' form of the result has one orientation: \strong{voxels x terms}, with one
#' row per voxel (or latent component) and one column per coefficient or
#' contrast, and the term names on the column margin. This matches
#' \code{\link{stats}()}, \code{\link{standard_error}()} and
#' \code{\link{p_values}()}, so results from the accessor family line up
#' element-wise (for example \code{coef(fit) / standard_error(fit)}).
#'
#' @param object A fitted \code{fmri_lm} object.
#' @param type \code{"betas"} (default) for regression coefficients, or
#'   \code{"contrasts"} for the estimates of the simple contrasts defined in
#'   the model.
#' @param include_baseline Logical. For \code{type = "betas"}, also return the
#'   baseline (drift, block intercept and nuisance) coefficients, in design
#'   matrix column order. The default, \code{FALSE}, returns event
#'   coefficients only.
#' @param recon Logical; used by latent-space fits to reconstruct voxel-space
#'   coefficients. Ignored for ordinary fits.
#' @param ... Unused.
#' @return For \code{type = "betas"}, a numeric matrix with one row per voxel
#'   and one column per coefficient, columns named after the design matrix
#'   (see \code{\link{coef_names}()}). For \code{type = "contrasts"}, a tibble
#'   with one row per voxel and one column per contrast.
#' @section Orientation change in fmrireg 0.2.0:
#'   Before 0.2.0 the default call (\code{type = "betas"},
#'   \code{include_baseline = FALSE}) returned the transpose, terms x voxels,
#'   while the other forms were voxels x terms. Code written for the old
#'   default that indexed \code{coef(fit)[term, ]} must now use
#'   \code{coef(fit)[, term]}.
#' @seealso \code{\link{coef_names}()}, \code{\link{stats}()} (t-statistics),
#'   \code{\link{standard_error}()}, \code{\link{p_values}()},
#'   \code{\link{coef_image}()}.
#' @examples
#' X <- matrix(rnorm(50 * 3), 50, 3)
#' edata <- data.frame(
#'   condition = factor(c("A", "B", "A", "B")),
#'   onsets = c(1, 12, 25, 38),
#'   run = 1
#' )
#' dset <- matrix_frame(X, TR = 2, run_length = 50, event_table = edata)
#' fit <- fmri_lm(onsets ~ hrf(condition), block = ~run, dataset = dset)
#' b <- coef(fit) # 3 voxels x 2 event coefficients
#' dim(b)
#' colnames(b)
#' dim(coef(fit, include_baseline = TRUE))
#' @method coef fmri_lm
#' @export
coef.fmri_lm <- function(object, type = c("betas", "contrasts"), include_baseline = FALSE, recon = FALSE, ...) {
  type <- match.arg(type)

  if (type == "contrasts") {
    return(pull_stat(object, "contrasts", "estimate"))
  }

  all_betas <- as.matrix(object$result$betas$data[[1]]$estimate[[1]])
  dm_colnames <- tryCatch(colnames(design_matrix(object$model)), error = function(e) NULL)

  if (include_baseline) {
    res <- all_betas
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
  } else {
    max_col <- ncol(all_betas)
    valid_event_indices <- object$result$event_indices[object$result$event_indices <= max_col]

    if (length(valid_event_indices) == 0) {
      warning("No valid event indices found in coef.fmri_lm. Using all available columns.")
      valid_event_indices <- seq_len(max_col)
    }

    res <- all_betas[, valid_event_indices, drop = FALSE]

    # Use the design-matrix column names rather than conditions(): they are
    # unique even when several terms share variables.
    if (!is.null(dm_colnames) && length(dm_colnames) >= max(valid_event_indices)) {
      colnames(res) <- dm_colnames[valid_event_indices]
    } else {
      condition_names <- conditions(object$model$event_model)[seq_along(valid_event_indices)]
      colnames(res) <- make.names(condition_names, unique = TRUE)
    }
  }

  # Voxels x terms (#227); rows are unnamed voxels, columns are terms.
  rownames(res) <- NULL
  res
}

#' Resolve the parameter-family selector of the fmri_lm accessors
#'
#' \code{"betas"} is the canonical name of the event-coefficient family; it
#' matches \code{coef(type = "betas")}. \code{"estimates"} is accepted as a
#' legacy synonym. In \code{stats()} it is deprecated, because
#' \code{stats(fit, "estimates")} reads as "the estimates" but returns
#' t-statistics (bbuchsbaum/fmrireg#217).
#' @keywords internal
#' @noRd
.accessor_family <- function(type, allowed, caller, deprecate_estimates = FALSE) {
  type <- match.arg(type, c(allowed, "estimates"))
  if (identical(type, "estimates")) {
    if (deprecate_estimates &&
        !isTRUE(getOption("fmrireg.suppress_deprecation", FALSE))) {
      rlang::warn(
        c(
          sprintf("`%s(type = \"estimates\")` is deprecated; use `type = \"betas\"`.", caller),
          i = sprintf(
            "`%s()` returns t-statistics for the betas, not the estimates. Use `coef()` for the estimates.",
            caller
          )
        ),
        class = c("fmrireg_deprecated_estimates_type", "deprecatedWarning")
      )
    }
    type <- "betas"
  }
  type
}

#' @method stats fmri_lm
#' @rdname stats
#' @export
stats.fmri_lm <- function(x, type = c("betas", "contrasts", "F"), ...) {
  type <- .accessor_family(type[1], c("betas", "contrasts", "F"), "stats",
                           deprecate_estimates = TRUE)
  pull_stat(x, type, "stat")
}

#' @method p_values fmri_lm
#' @rdname p_values
#' @export
p_values.fmri_lm <- function(x, type = c("betas", "contrasts"), ...) {
  type <- .accessor_family(type[1], c("betas", "contrasts"), "p_values")
  pull_stat(x, type, "prob")
}

#' @method standard_error fmri_lm
#' @rdname standard_error
#' @export
standard_error.fmri_lm <- function(x, type = c("betas", "contrasts"), ...) {
  type <- .accessor_family(type[1], c("betas", "contrasts"), "standard_error")
  pull_stat(x, type, "se")
}


#' @noRd
#' @keywords internal
.tidy_stat_block <- function(tbl, block = c("betas", "contrasts")) {
  block <- match.arg(block)
  if (!inherits(tbl, "tbl_df") || !"data" %in% names(tbl)) {
    stop("Expected tibble with nested data column", call. = FALSE)
  }

  nested <- tbl$data[[1]]
  if (nrow(nested) == 0) {
    return(tibble::tibble())
  }

  # sigma is solver-specific metadata, not a statistic needed by tidy().
  # In particular, runwise meta-estimation returns the four inferential
  # matrices below without a sigma column.
  key_cols <- c("estimate", "se", "stat", "prob")

  if (!all(key_cols %in% names(nested))) {
    stop("Unexpected structure in nested statistics table", call. = FALSE)
  }

  cols <- lapply(key_cols, function(nm) {
    values <- nested[[nm]]
    if (length(values) == 0L) {
      return(NULL)
    }

    # Beta statistics are normally stored as a one-element list-column of
    # matrices. Contrast statistics are normally atomic vectors, but some
    # backends use the same list-column representation. Normalize both forms
    # without discarding matrix dimensions.
    if (is.list(values)) {
      if (length(values) != 1L || is.null(values[[1L]])) {
        return(NULL)
      }
      values[[1L]]
    } else {
      values
    }
  })
  names(cols) <- key_cols
  cols
}

#' @noRd
#' @keywords internal
.tidy_fmri_stats <- function(model, type = c("estimates", "contrasts")) {
  type <- match.arg(type)
  block <- if (type == "estimates") model$result$betas else model$result$contrasts

  if (is.null(block) || nrow(block) == 0) {
    return(tibble::tibble())
  }

  mats <- .tidy_stat_block(block, block = if (type == "estimates") "betas" else "contrasts")
  if (length(mats) == 0 || any(vapply(mats, is.null, logical(1)))) {
    return(tibble::tibble())
  }
  estimate_mat <- as.matrix(mats$estimate)
  se_mat <- as.matrix(mats$se)
  stat_mat <- as.matrix(mats$stat)
  prob_mat <- as.matrix(mats$prob)

  stat_dims <- lapply(
    list(estimate = estimate_mat, se = se_mat, stat = stat_mat, prob = prob_mat),
    dim
  )
  if (!all(vapply(stat_dims[-1L], identical, logical(1), stat_dims[[1L]]))) {
    stop("Estimate, SE, statistic, and p-value layouts do not agree", call. = FALSE)
  }

  add_names <- function(target, reference) {
    if (is.null(colnames(target)) && !is.null(colnames(reference))) {
      colnames(target) <- colnames(reference)
    }
    target
  }

  stat_mat <- add_names(stat_mat, estimate_mat)
  se_mat <- add_names(se_mat, estimate_mat)
  prob_mat <- add_names(prob_mat, estimate_mat)
  estimate_mat <- add_names(estimate_mat, stat_mat)

  if (is.null(colnames(estimate_mat))) {
    design_names <- tryCatch(
      colnames(design_matrix(model$model)),
      error = function(e) NULL
    )
    colind <- block$colind[[1L]] %||% NULL
    if (!is.null(colind) && length(colind) == ncol(estimate_mat) &&
        length(design_names) >= max(colind)) {
      design_names <- design_names[colind]
    }
    if (length(design_names) == ncol(estimate_mat)) {
      colnames(estimate_mat) <- design_names
    } else {
      colnames(estimate_mat) <- paste0("term", seq_len(ncol(estimate_mat)))
    }
  }
  se_mat <- add_names(se_mat, estimate_mat)
  stat_mat <- add_names(stat_mat, estimate_mat)
  prob_mat <- add_names(prob_mat, estimate_mat)

  term_names <- colnames(estimate_mat)
  if (type == "contrasts" && !is.null(block$conmat) && length(block$conmat) > 0) {
    cm <- block$conmat[[1]]
    if (!is.null(colnames(cm))) {
      term_names <- colnames(cm)
    }
  }

  estimate_df <- as.data.frame(estimate_mat)
  se_df <- as.data.frame(se_mat)
  stat_df <- as.data.frame(stat_mat)
  prob_df <- as.data.frame(prob_mat)
  colnames(estimate_df) <- term_names
  colnames(se_df) <- term_names
  colnames(stat_df) <- term_names
  colnames(prob_df) <- term_names

  n_vox <- nrow(estimate_df)
  term_names <- colnames(estimate_df)

  result <- tibble::tibble(
    voxel = rep(seq_len(n_vox), each = length(term_names)),
    term = rep(term_names, times = n_vox),
    # Labels enumerate all terms within a voxel. R matrices flatten by column,
    # so transpose first to preserve that row-major semantic order.
    estimate = as.vector(t(as.matrix(estimate_df))),
    std_error = as.vector(t(as.matrix(se_df))),
    statistic = as.vector(t(as.matrix(stat_df))),
    p_value = as.vector(t(as.matrix(prob_df)))
  )

  inference_df <- model$result$df$inference %||% NULL
  if (!is.null(inference_df)) {
    if (length(inference_df) == 1L) inference_df <- rep(inference_df, n_vox)
    if (length(inference_df) == n_vox) {
      result$df_inference <- rep(as.numeric(inference_df), each = length(term_names))
    }
  }
  if (!is.null(model$result$df$nominal)) {
    result$df_residual <- rep(as.numeric(model$result$df$nominal)[1L], nrow(result))
  } else if (!is.null(block$df.residual)) {
    result$df_residual <- rep(block$df.residual[1], nrow(result))
  }

  result
}

#' @export
tidy.fmri_lm <- function(x, type = c("estimates", "contrasts"), ...) {
  stats_tbl <- .tidy_fmri_stats(x, match.arg(type))
  clean_condition <- function(label) {
    label <- gsub("^conditioncondition_condition\\.", "", label)
    label <- gsub("^condition_condition\\.", "", label)
    gsub("\\.", " ", label)
  }
  if (nrow(stats_tbl)) {
    stats_tbl$term <- clean_condition(stats_tbl$term)
  }
  stats_tbl
}

#' @method print fmri_lm
#' @export
print.fmri_lm <- function(x, ...) {
  cli::cli_h1("fMRI Linear Model Results")
  
  # Model info
  cli::cli_h2("Model Information")
  cli::cli_ul()
  cli::cli_li("Dataset: {.field {class(x$dataset)[1]}}")
  scope <- attr(x, "requested_control")$estimation$scope %||%
    attr(x, "config")$estimation$scope %||% attr(x, "strategy")
  cli::cli_li("Estimation scope: {.field {scope}}")
  
  # Design info
  n_events <- length(x$result$event_indices)
  n_baseline <- length(x$result$baseline_indices)
  cli::cli_li("Parameters: {.val {n_events}} event + {.val {n_baseline}} baseline")
  
  # Data dimensions
  beta_mat <- x$result$betas$data[[1]]$estimate[[1]]
  n_voxels <- nrow(beta_mat)
  cli::cli_li("Voxels analyzed: {.val {n_voxels}}")
  
  # Degrees of freedom
  df_values <- as.numeric(x$result$df$inference %||% x$result$betas$df.residual[1])
  df_values <- df_values[is.finite(df_values)]
  if (length(df_values)) {
    df_label <- if (length(unique(df_values)) == 1L) {
      format(df_values[[1L]], digits = 5L)
    } else {
      sprintf("%s to %s", format(min(df_values), digits = 5L),
              format(max(df_values), digits = 5L))
    }
    df_method <- x$result$df$method %||% "residual"
    cli::cli_li("Inference df ({.field {df_method}}): {.val {df_label}}")
  }
  cli::cli_end()
  
  # Contrasts info if available
  if (!is.null(x$result$contrasts) && nrow(x$result$contrasts) > 0) {
    cli::cli_h2("Contrasts")
    cli::cli_ul()
    
    # Simple contrasts
    simple_cons <- x$result$contrasts %>% dplyr::filter(type == "contrast")
    if (nrow(simple_cons) > 0) {
      cli::cli_li("Simple contrasts: {.val {nrow(simple_cons)}}")
      con_names <- paste(simple_cons$name, collapse = ", ")
      cli::cli_text("  {.emph {con_names}}")
    }
    
    # F contrasts  
    f_cons <- x$result$contrasts %>% dplyr::filter(type == "Fcontrast")
    if (nrow(f_cons) > 0) {
      cli::cli_li("F-contrasts: {.val {nrow(f_cons)}}")
      fcon_names <- paste(f_cons$name, collapse = ", ")
      cli::cli_text("  {.emph {fcon_names}}")
    }
    cli::cli_end()
  }
  
  # Config info if available
  if (!is.null(attr(x, "config"))) {
    cfg <- attr(x, "config")
    cli::cli_h2("Model Configuration")
    cli::cli_ul()
    
    # AR info
    noise <- cfg$noise %||% cfg$ar
    if (noise$struct != "iid") {
      cli::cli_li("AR structure: {.field {noise$struct}}")
      cli::cli_li("Noise pooling: {.field {noise$pooling %||% 'run'}}")
      if (noise$voxelwise) {
        cli::cli_li("AR estimation: {.emph voxelwise}")
      } else {
        cli::cli_li("Shared AR estimator: {.field {noise$shared_estimator %||% 'pooled_acvf'}}")
      }
    }
    
    # Robust info
    if (.fmri_lm_robust_enabled(cfg$robust)) {
      cli::cli_li("Robust method: {.field {cfg$robust$type}}")
      cli::cli_li("Robust tuning: {.val {cfg$robust$c_tukey}}")
    }
    cli::cli_li("Variance: {.field {cfg$variance$method}}")
    cli::cli_li("Reference df: {.field {cfg$variance$df}}")
    cli::cli_end()
  }
  
  cli::cli_rule()
  cli::cli_text("{.emph Use coef(), stats(), or standard_error() to extract results}")

  invisible(x)
}

#' @rdname coef_names
#' @param type Which set of names to return: \code{"estimates"} (default) for
#'   event regressor names, \code{"contrasts"} for simple contrast names,
#'   \code{"F"} for F-contrast names, or \code{"all"} for a named list of all
#'   three.
#' @examples
#' # Create a small example
#' X <- matrix(rnorm(50 * 4), 50, 4)
#' edata <- data.frame(
#'   condition = factor(c("A", "B", "A", "B")),
#'   onsets = c(1, 12, 25, 38),
#'   run = c(1, 1, 1, 1)
#' )
#' dset <- matrix_frame(X, TR = 2, run_length = 50,
#'                                     event_table = edata)
#' fit <- fmri_lm(onsets ~ hrf(condition), block = ~run, dataset = dset)
#' coef_names(fit)
#' @method coef_names fmri_lm
#' @export
coef_names.fmri_lm <- function(x, type = c("estimates", "contrasts", "F", "all"), ...) {
  type <- match.arg(type)

  get_estimate_names <- function() {
    dm <- design_matrix(x$model)
    ei <- x$result$event_indices
    max_col <- if (!is.null(dm)) ncol(dm) else max(ei)
    valid_idx <- ei[ei <= max_col]
    if (!is.null(dm) && length(valid_idx) > 0) {
      colnames(dm)[valid_idx]
    } else {
      conds <- conditions(x$model$event_model)
      make.names(conds[seq_along(valid_idx)], unique = TRUE)
    }
  }

  get_contrast_names <- function() {
    ct <- x$result$contrasts
    if (is.null(ct) || nrow(ct) == 0) return(character(0))
    simple <- ct[ct$type == "contrast", , drop = FALSE]
    if (nrow(simple) == 0) return(character(0))
    simple$name
  }

  get_f_names <- function() {
    ct <- x$result$contrasts
    if (is.null(ct) || nrow(ct) == 0) return(character(0))
    fcons <- ct[ct$type == "Fcontrast", , drop = FALSE]
    if (nrow(fcons) == 0) return(character(0))
    fcons$name
  }

  switch(type,
    estimates = get_estimate_names(),
    contrasts = get_contrast_names(),
    "F" = get_f_names(),
    all = list(
      estimates = get_estimate_names(),
      contrasts = get_contrast_names(),
      "F" = get_f_names()
    )
  )
}

#' @rdname coef_image
#' @param statistic For \code{fmri_lm} objects: one of \code{"estimate"},
#'   \code{"se"}, \code{"tstat"}, or \code{"prob"}.
#' @param type For \code{fmri_lm} objects: which coefficient set to index into:
#'   \code{"estimates"} (default), \code{"contrasts"}, or \code{"F"}.
#' @param ... Additional arguments (currently unused).
#' @examples
#' # Create a small example
#' X <- matrix(rnorm(50 * 4), 50, 4)
#' edata <- data.frame(
#'   condition = factor(c("A", "B", "A", "B")),
#'   onsets = c(1, 12, 25, 38),
#'   run = c(1, 1, 1, 1)
#' )
#' dset <- matrix_frame(X, TR = 2, run_length = 50,
#'                                     event_table = edata)
#' fit <- fmri_lm(onsets ~ hrf(condition), block = ~run, dataset = dset)
#' # Get coefficient estimates as a numeric vector
#' coef_image(fit, coef = 1)
#' @method coef_image fmri_lm
#' @export
coef_image.fmri_lm <- function(object, coef = 1,
                                statistic = c("estimate", "se", "tstat", "prob"),
                                type = c("estimates", "contrasts", "F"),
                                ...) {
  statistic <- match.arg(statistic)
  type <- match.arg(type)

  # Map statistic to the internal element names used by pull_stat
  element <- switch(statistic,
    estimate = "estimate",
    se       = "se",
    tstat    = "stat",
    prob     = "prob"
  )

  # ---- resolve coef to a column index ----
  available <- coef_names(object, type = type)
  if (length(available) == 0) {
    stop("No coefficients of type '", type, "' available in this model.",
         call. = FALSE)
  }

  if (is.character(coef)) {
    idx <- match(coef, available)
    if (is.na(idx)) {
      stop("Coefficient '", coef, "' not found. Available names: ",
           paste(available, collapse = ", "), call. = FALSE)
    }
  } else {
    idx <- as.integer(coef)
    if (idx < 1 || idx > length(available)) {
      stop("Coefficient index ", idx, " out of range [1, ",
           length(available), "].", call. = FALSE)
    }
  }

  # ---- extract the values vector for the requested coefficient ----
  if (type == "estimates") {
    mat <- object$result$betas$data[[1]][[element]][[1]]
    ei <- object$result$event_indices
    valid_idx <- ei[ei <= ncol(mat)]
    values <- mat[, valid_idx[idx]]
  } else if (type == "contrasts") {
    ct <- object$result$contrasts
    simple <- ct[ct$type == "contrast", , drop = FALSE]
    values <- as.vector(simple$data[[idx]][[element]])
  } else {
    # F-contrasts
    ct <- object$result$contrasts
    fcons <- ct[ct$type == "Fcontrast", , drop = FALSE]
    values <- as.vector(fcons$data[[idx]][[element]])
  }

  # ---- reconstruct spatial image if possible ----
  tryCatch({
    spatial <- .fmri_dataset_mask_space(object$dataset, "coefficient image reconstruction")
    sp <- spatial$space
    mask_idx <- which(spatial$mask_array)
    vol_array <- array(NA_real_, dim(sp))
    vol_array[mask_idx] <- values
    neuroim2::NeuroVol(vol_array, sp)
  }, error = function(e) {
    # Non-spatial dataset: return raw vector with informative name
    names(values) <- NULL
    values
  })
}

#' @rdname coef_images
#' @param statistic For \code{fmri_lm} objects: one of \code{"estimate"},
#'   \code{"se"}, \code{"tstat"}, or \code{"prob"}.
#' @param type For \code{fmri_lm} objects: which coefficient set to extract:
#'   \code{"estimates"} (default), \code{"contrasts"}, or \code{"F"}.
#' @param coefs Optional character vector of coefficient names (a subset of
#'   \code{coef_names(object, type = type)}) to extract. Defaults to all
#'   coefficients of the requested \code{type}.
#' @examples
#' # Create a small example
#' X <- matrix(rnorm(50 * 4), 50, 4)
#' edata <- data.frame(
#'   condition = factor(c("A", "B", "A", "B")),
#'   onsets = c(1, 12, 25, 38),
#'   run = c(1, 1, 1, 1)
#' )
#' dset <- matrix_frame(X, TR = 2, run_length = 50,
#'                                     event_table = edata)
#' fit <- fmri_lm(onsets ~ hrf(condition), block = ~run, dataset = dset)
#' # Named list of one volume per event regressor
#' imgs <- coef_images(fit, statistic = "estimate", type = "estimates")
#' names(imgs)
#' @method coef_images fmri_lm
#' @export
coef_images.fmri_lm <- function(object,
                                statistic = c("estimate", "se", "tstat", "prob"),
                                type = c("estimates", "contrasts", "F"),
                                coefs = NULL,
                                ...) {
  statistic <- match.arg(statistic)
  type <- match.arg(type)

  available <- coef_names(object, type = type)
  if (length(available) == 0) {
    return(stats::setNames(list(), character(0)))
  }

  if (is.null(coefs)) {
    coefs <- available
  } else {
    missing <- setdiff(coefs, available)
    if (length(missing) > 0) {
      stop("Coefficient(s) not found for type '", type, "': ",
           paste(missing, collapse = ", "),
           ". Available names: ", paste(available, collapse = ", "),
           call. = FALSE)
    }
  }

  images <- lapply(coefs, function(nm) {
    coef_image(object, coef = nm, statistic = statistic, type = type)
  })
  names(images) <- coefs
  images
}
