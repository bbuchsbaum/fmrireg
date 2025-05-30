# S3 Methods for fmri_lm Objects
# Methods for extracting results and information from fitted fmri_lm objects


#' Extract HRF from Fitted Model
#'
#' @param x An fmri_lm object
#' @param sample_at Time points at which to sample the HRF
#' @param ... Additional arguments
#' @return A data frame with the fitted HRF
#' @export
fitted_hrf.fmri_lm <- function(x, sample_at = seq(0, 24, by = 1), ...) {
  # Create time vector
  time <- sample_at
  
  # Get the HRF from the model
  hrf_obj <- getHRF(x$model$event_model)
  
  # Check if it's a basis HRF
  if (inherits(hrf_obj, "hrfbasis")) {
    # Get basis functions
    basis_mat <- predict(hrf_obj, time)
    
    # Get coefficients for the first event term
    # This is simplified - in practice might want to specify which term
    event_indices <- x$result$event_indices
    coefs <- coef(x, type = "betas", include_baseline = FALSE)
    
    # Take first few coefficients corresponding to basis
    nbasis <- ncol(basis_mat)
    if (ncol(coefs) >= nbasis) {
      basis_coefs <- as.matrix(coefs[1, 1:nbasis])
      fitted_values <- basis_mat %*% t(basis_coefs)
    } else {
      warning("Not enough coefficients for basis functions")
      fitted_values <- rep(0, length(time))
    }
  } else {
    # Standard HRF - just evaluate it
    fitted_values <- evaluate(hrf_obj, time)
  }
  
  data.frame(
    time = time,
    value = as.vector(fitted_values)
  )
}

#' Reshape Coefficients
#' 
#' @keywords internal
#' @noRd
reshape_coef <- function(df, des, measure = "value") {
  nvox <- nrow(df)
  ncon <- ncol(df)
  nrun <- length(levels(des$blockids))
  ret <- matrix(0, nvox * nrun, ncon / nrun)
  
  runids <- sort(unique(des$blockids))
  
  currow <- 1
  for (i in 1:length(runids)) {
    runs <- which(des$blockids == runids[i])
    ind <- seq(currow, currow + nvox - 1)
    currow <- currow + nvox
    for (j in 1:length(runs)) {
      ret[ind, j] <- df[, runs[j]]
    }
  }
  
  as.data.frame(ret)
}

#' Pull Statistics from fmri_lm Results
#'
#' @keywords internal
#' @noRd
pull_stat_revised <- function(x, type, element) {
  if (type == "betas") {
    # Extract the beta data from the tibble structure
    beta_data <- x$result$betas$data[[1]]
    ret <- beta_data[[element]][[1]]
    
    # Get valid event indices
    max_col <- ncol(ret)
    valid_event_indices <- x$result$event_indices[x$result$event_indices <= max_col]
    
    if (length(valid_event_indices) == 0) {
      warning("No valid event indices found. Using all available columns.")
      valid_event_indices <- 1:max_col
    }
    
    # Extract only event-related betas
    ret <- ret[, valid_event_indices, drop = FALSE]
    
    # Set column names
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

#' Pull Statistics (Legacy Version)
#'
#' @keywords internal
#' @noRd
pull_stat <- function(x, type, element) {
  # Redirect to revised version
  pull_stat_revised(x, type, element)
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
      colnames(res) <- colnames(design_matrix(object$model))
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
      
      # Ensure tibble output for consistency with original behavior
      res <- suppressMessages(tibble::as_tibble(res, .name_repair = "check_unique"))
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
