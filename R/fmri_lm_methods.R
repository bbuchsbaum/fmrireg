# S3 Methods for fmri_lm Objects
# Methods for extracting results and information from fitted fmri_lm objects


## Deprecated stub: fitted_hrf.fmri_lm
## The real implementation lives in R/fmrilm.R. The placeholder
## implementation previously here has been removed to avoid clobbering
## the robust method and to prevent duplicate S3 registrations.
## If needed for reference, see git history prior to this change.

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
      # Only assign column names up to the number of columns in the beta matrix
      dm_colnames <- colnames(design_matrix(object$model))
      if (ncol(res) <= length(dm_colnames)) {
        colnames(res) <- dm_colnames[1:ncol(res)]
      } else {
        # This shouldn't happen, but be defensive
        colnames(res) <- paste0("beta_", 1:ncol(res))
      }
      # Return as matrix - voxels x predictors
      # res <- t(res)  # Removed transpose to maintain standard orientation
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
      
      # Return as matrix - voxels x predictors
      # res <- t(res)  # Removed transpose to maintain standard orientation
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
  
  element <- "stat"
  
  if (type == "estimates") {
    pull_stat(x, "betas", element)
  } else {
    pull_stat(x, type, element)
  }
}

#' @method p_values fmri_lm
#' @rdname p_values
#' @export
p_values.fmri_lm <- function(x, type = c("estimates", "contrasts"), ...) {
  type <- match.arg(type)
  
  element <- "prob"
  
  if (type == "estimates") {
    pull_stat(x, "betas", element)
  } else {
    pull_stat(x, type, element)
  }
}

#' @method standard_error fmri_lm
#' @rdname standard_error
#' @export
standard_error.fmri_lm <- function(x, type = c("estimates", "contrasts"),...) {
  type <- match.arg(type)
  
  element <- "se"
  
  if (type == "estimates") {
    pull_stat(x, "betas", element)
  } else {
    pull_stat(x, type, element)
  }
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

  key_cols <- switch(block,
    betas = c("estimate", "se", "stat", "prob", "sigma"),
    contrasts = c("estimate", "se", "stat", "prob", "sigma")
  )

  if (!all(key_cols %in% names(nested))) {
    stop("Unexpected structure in nested statistics table", call. = FALSE)
  }

  cols <- lapply(key_cols, function(nm) {
    mats <- nested[[nm]]
    if (length(mats) == 0 || is.null(mats[[1]])) {
      NULL
    } else {
      mats[[1]]
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
    colnames(estimate_mat) <- paste0("term", seq_len(ncol(estimate_mat)))
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
    estimate = as.vector(as.matrix(estimate_df)),
    std_error = as.vector(as.matrix(se_df)),
    statistic = as.vector(as.matrix(stat_df)),
    p_value = as.vector(as.matrix(prob_df))
  )

  if (!is.null(block$df.residual)) {
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
  cli::cli_li("Strategy: {.field {attr(x, 'strategy')}}")
  
  # Design info
  n_events <- length(x$result$event_indices)
  n_baseline <- length(x$result$baseline_indices)
  cli::cli_li("Parameters: {.val {n_events}} event + {.val {n_baseline}} baseline")
  
  # Data dimensions
  beta_mat <- x$result$betas$data[[1]]$estimate[[1]]
  n_voxels <- nrow(beta_mat)
  cli::cli_li("Voxels analyzed: {.val {n_voxels}}")
  
  # Degrees of freedom
  df_resid <- x$result$betas$df.residual[1]
  cli::cli_li("Residual df: {.val {df_resid}}")
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
    if (cfg$ar$struct != "iid") {
      cli::cli_li("AR structure: {.field {cfg$ar$struct}}")
      if (cfg$ar$global) cli::cli_li("AR scope: {.emph global}")
      if (cfg$ar$voxelwise) cli::cli_li("AR estimation: {.emph voxelwise}")
    }
    
    # Robust info
    if (cfg$robust$type != FALSE) {
      cli::cli_li("Robust method: {.field {cfg$robust$type}}")
      cli::cli_li("Robust tuning: {.val {cfg$robust$c_tukey}}")
    }
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
#' dset <- fmridataset::matrix_dataset(X, TR = 2, run_length = 50,
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
#' dset <- fmridataset::matrix_dataset(X, TR = 2, run_length = 50,
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
    values <- simple$data[[idx]][[element]][[1]]
  } else {
    # F-contrasts
    ct <- object$result$contrasts
    fcons <- ct[ct$type == "Fcontrast", , drop = FALSE]
    values <- fcons$data[[idx]][[element]][[1]]
  }

  # ---- reconstruct spatial image if possible ----
  tryCatch({
    mask <- fmridataset::get_mask(object$dataset)
    sp <- neuroim2::space(mask)
    mask_idx <- which(as.logical(mask))
    vol_array <- array(NA_real_, dim(sp))
    vol_array[mask_idx] <- values
    neuroim2::NeuroVol(vol_array, sp)
  }, error = function(e) {
    # Non-spatial dataset: return raw vector with informative name
    names(values) <- NULL
    values
  })
}
