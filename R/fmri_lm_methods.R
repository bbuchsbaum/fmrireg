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
