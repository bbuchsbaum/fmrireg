###############################################################################
## fmrimodel.R
##
## This file creates the overall fMRI model by combining an event model 
## (describing experimental events) and a baseline model (modeling drift, 
## nuisance, and block effects).
##
## The file provides two main functions:
##   - create_fmri_model(): Creates an fMRI model from a formula, block formula,
##     (optionally) a baseline_model, and an fmri_dataset.
##   - fmri_model(): Combines an event_model and a baseline_model into an fmri_model.
##
## Additional functions build design matrices, compute contrasts, print, and plot
## the overall fMRI model.
###############################################################################

## ============================================================================
## Section 1: fMRI Model Construction Functions
## ============================================================================

#' Create an fMRI Model
#'
#' This function creates an fMRI model consisting of an event model and a baseline model.
#'
#' @param formula The model formula for experimental events.
#' @param block The model formula for block structure.
#' @param baseline_model (Optional) A \code{baseline_model} object. Default is \code{NULL}.
#' @param dataset An \code{fmri_dataset} object containing the time-series data.
#' @param drop_empty Logical. Whether to remove factor levels with zero size. Default is \code{TRUE}.
#' @param durations A vector of event durations. Default is \code{0}.
#' @return An \code{fmri_model} object.
#' @export
#' @examples
#' \dontrun{
#' # Assuming you have an fmri_dataset object named ds and a formula for events:
#' fmri_mod <- create_fmri_model(formula = onset ~ hrf(x) + hrf(y),
#'                               block = ~ run,
#'                               dataset = ds,
#'                               drop_empty = TRUE,
#'                               durations = rep(0, nrow(ds$event_table)))
#' }
create_fmri_model <- function(formula, block, baseline_model = NULL, dataset, drop_empty = TRUE, durations = 0) {
  assert_that(is.formula(formula), msg = "'formula' must be a formula")
  assert_that(is.formula(block), msg = "'block' must be a formula")
  assert_that(inherits(dataset, "fmri_dataset"), msg = "'dataset' must be an 'fmri_dataset'")
  assert_that(is.numeric(durations), msg = "'durations' must be numeric")
  
  # Replicate durations if a single value is provided.
  if (length(durations) == 1) {
    durations <- rep(durations, nrow(dataset$event_table))
  }
  
  # Resolve conflict: use a temporary variable to hold the baseline model.
  if (is.null(baseline_model)) {
    base_model_obj <- baseline_model(
      basis = "bs",
      degree = max(ceiling(median(fmrihrf::blocklens(dataset$sampling_frame)) / 100), 3),
      sframe = dataset$sampling_frame
    )
  } else {
    assert_that(inherits(baseline_model, "baseline_model"),
                msg = "'baseline_model' must have class 'baseline_model'")
    base_model_obj <- baseline_model
  }
  
  ev_model <- event_model(
    x = formula,
    block = block,
    data = dataset$event_table,
    sampling_frame = dataset$sampling_frame,
    drop_empty = drop_empty,
    durations = durations
  )
  
  fmri_model(ev_model, base_model_obj, dataset)
}


#' Construct an fMRI Regression Model
#'
#' This function constructs an fMRI regression model consisting of an event model
#' and a baseline model. The resulting model can be used for the analysis of fMRI data.
#'
#' @param event_model An object of class "event_model" representing the event-related part of the fMRI regression model.
#' @param baseline_model An object of class "baseline_model" representing the baseline-related part of the fMRI regression model.
#' @param dataset An \code{fmri_dataset} used to build the model.
#' @return An object of class \code{fmri_model} containing the event and baseline models along with the dataset.
#' @export
#' @seealso event_model, baseline_model
fmri_model <- function(event_model, baseline_model, dataset) {
  assert_that(inherits(event_model, "event_model"))
  assert_that(inherits(baseline_model, "baseline_model"))
  assert_that(inherits(dataset, "fmri_dataset"))

  fmodel <- list(event_model = event_model,
                 baseline_model = baseline_model,
                 dataset = dataset)
  class(fmodel) <- "fmri_model"
  fmodel
}


#' (Internal) Prediction Matrix
#'
#' This function is intended to compute a prediction matrix for the model.
#' (Currently a stub.)
#'
#' @param x An fmri_model object.
#' @return (Not implemented)
#' @keywords internal
#' @noRd
prediction_matrix <- function(x) {
  stop("not implemented")
}


## ============================================================================
## Section 2: Design Matrix and Environment for fMRI Models
## ============================================================================

#' Design Matrix for fMRI Models
#' 
#' Extract the combined design matrix from an fMRI model containing both event and baseline terms.
#' 
#' @param x An fmri_model object
#' @param blockid Optional numeric vector specifying which blocks/runs to include
#' @param ... Additional arguments (not used)
#' @return A tibble containing the combined design matrix with event and baseline terms
#' @method design_matrix fmri_model
#' @export
#' @importFrom tibble as_tibble
design_matrix.fmri_model <- function(x, blockid = NULL, ...) {
  suppressMessages(
    tibble::as_tibble(
      cbind(
        design_matrix(x$event_model, blockid),
        design_matrix(x$baseline_model, blockid)
      ),
      .name_repair = "check_unique"
    )
  )
}


#' @importFrom tibble as_tibble
#' @keywords internal
#' @noRd
design_env.fmri_model <- function(x, blockid = NULL) {
  stop("Not implemented")
}


## ============================================================================
## Section 3: Accessor Functions for fMRI Models
## ============================================================================

#' @export
terms.fmri_model <- function(x, ...) {
  c(terms(x$event_model), terms(x$baseline_model))
}

#' @export
#' @autoglobal
cells.fmri_model <- function(x, ...) {
  c1 <- cells(x$event_model) %>% dplyr::mutate(type = "event")
  c2 <- cells(x$baseline_model) %>% dplyr::mutate(type = "baseline")
  rbind(c1, c2) %>% dplyr::relocate(index, type)
}

#' @export
blocklens.fmri_model <- function(x, ...) {
  fmrihrf::blocklens(x$event_model)
}

#' @export
event_terms.fmri_model <- function(x, ...) {
  terms(x$event_model)
}

#' @export
baseline_terms.fmri_model <- function(x, ...) {
  terms(x$baseline_model)
}

#' @export
contrast_weights.fmri_model <- function(x, ...) {
  contrast_weights(x$event_model, ...)
}

#' @export
conditions.fmri_model <- function(x, ...) {
  unlist(lapply(terms(x), function(t) conditions(t)), use.names = FALSE)
}

#' @export
conditions.baseline_model <- function(x, ...) {
  unlist(lapply(terms(x), function(t) conditions(t)), use.names = FALSE)
}


## ============================================================================
## Section 4: Plot and Print Methods for fMRI Models
## ============================================================================

#' @export
plot.fmri_model <- function(x, ...) {
  with_package("cowplot")
  p1 <- plot(x$event_model) + ggplot2::ggtitle("Event Model")
  p2 <- plot(x$baseline_model) + ggplot2::ggtitle("Baseline Model")
  cowplot::plot_grid(p1, p2, nrow = 2, align = "h")
}

#' @export
#' @rdname print
print.fmri_model <- function(x, ...) {
  # Header with fancy border
  cat("\n=============================================")
  cat("\n             fMRI Model                     ")
  cat("\n=============================================")
  
  # Event Model Section
  cat("\n Event Model                                ")
  cat("\n---------------------------------------------")
  cat("\n║ Formula:", crayon::cyan(Reduce(paste, deparse(x$event_model$model_spec$formula))))
  
  # Event Model Summary
  cat("\n║ Summary:")
  cat("\n║   • Terms:", crayon::yellow(length(terms(x$event_model))))
  cat("\n║   • Events:", crayon::yellow(nrow(x$event_model$model_spec$event_table)))
  cat("\n║   • Design Columns:", crayon::yellow(length(conditions(x$event_model))))
  cat("\n║   • Blocks:", crayon::yellow(length(unique(x$event_model$blockids))))
  
  # Baseline Model Section (if present)
  if (!is.null(x$baseline_model)) {
    cat("\n---------------------------------------------")
    cat("\n║ Baseline Model                           ║")
    cat("\n║ Components:")
    
    # Drift term info
    if (!is.null(x$baseline_model$drift_term)) {
      drift_name <- x$baseline_model$drift_term$varname
      basis_type <- x$baseline_model$drift_spec$basis
      degree <- x$baseline_model$drift_spec$degree
      drift_cols <- ncol(design_matrix(x$baseline_model$drift_term))
      cat("\n║   • Drift:", crayon::magenta(drift_name))
      cat("\n║     - Type:", crayon::blue(basis_type))
      cat("\n║     - Degree:", crayon::blue(degree))
      cat("\n║     - Columns:", crayon::yellow(drift_cols))
    }
    
    # Block term info
    if (!is.null(x$baseline_model$block_term)) {
      const_cols <- ncol(design_matrix(x$baseline_model$block_term))
      cat("\n║   • Block Terms:", crayon::yellow(const_cols), "columns")
    }
    
    # Nuisance term info
    if (!is.null(x$baseline_model$nuisance_term)) {
      nuis_cols <- ncol(design_matrix(x$baseline_model$nuisance_term))
      cat("\n║   • Nuisance Terms:", crayon::yellow(nuis_cols), "columns")
    }
  }
  
  # Total Model Summary
  cat("\n---------------------------------------------")
  cat("\n║ Total Model                              ║")
  total_cols <- ncol(design_matrix(x))
  cat("\n║   • Total Design Columns:", crayon::yellow(total_cols))
  
  # Footer
  cat("\n╚══════════════════════════════════════════╝\n")
}

#' correlation_map.fmri_model
#'
#' @description
#' Generates a correlation heatmap of the columns in an \code{fmri_model}'s combined
#' event+baseline design matrix.
#'
#' @param x An \code{fmri_model}.
#' @param method Correlation method (e.g., "pearson", "spearman").
#' @param half_matrix Logical; if TRUE, display only the lower triangle of the matrix.
#' @param absolute_limits Logical; if TRUE, set color scale limits from -1 to 1.
#' @param ... Additional arguments passed to internal plotting functions.
#' @export
correlation_map.fmri_model <- function(x,
                                       method          = c("pearson", "spearman"),
                                       half_matrix     = FALSE,
                                       absolute_limits = TRUE,
                                       ...) {
  DM <- as.matrix(design_matrix(x))
  .correlation_map_common(DM, method=method, half_matrix=half_matrix,
                          absolute_limits=absolute_limits, ...)
}


#' Heatmap visualization of the combined fmri_model design matrix
#'
#' @description
#' Produces a single heatmap of *all* columns in the design matrix from an
#' \code{fmri_model} object, which merges both the event_model and baseline_model
#' regressors. Rows are scans; columns are regressors.  
#' Optionally draws horizontal lines between blocks (runs), and rotates x‐axis
#' labels diagonally for readability.
#'
#' @param x An \code{fmri_model} object.
#' @param block_separators Logical; if \code{TRUE}, draw white horizontal lines between blocks.
#' @param rotate_x_text Logical; if \code{TRUE}, rotate x-axis labels by 45 degrees.
#' @param fill_midpoint Numeric or \code{NULL}; if not \code{NULL}, passed to
#'   \code{\link[ggplot2]{scale_fill_gradient2}} to center the color scale (e.g. \code{fill_midpoint=0}).
#' @param fill_limits Numeric vector of length 2 or \code{NULL}; passed to the fill scale
#'   \code{limits=} argument. This can clip or expand the color range.
#' @param ... Additional arguments passed to \code{\link[ggplot2]{geom_tile}}.
#'
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#' @return A ggplot2 plot object.
#' @export
design_map.fmri_model <- function(x,
                                  block_separators = TRUE,
                                  rotate_x_text    = TRUE,
                                  fill_midpoint    = NULL,
                                  fill_limits      = NULL,
                                  ...) {
  # 1) Extract full design matrix (event + baseline)
  DM <- design_matrix(x)
  n_scans <- nrow(DM)
  
  # 2) Convert to a long data frame for ggplot
  df_long <- tibble::as_tibble(DM, .name_repair = "unique")
  df_long$scan_number <- seq_len(n_scans)
  df_long <- tidyr::pivot_longer(
    df_long,
    cols      = -scan_number,
    names_to  = "Regressor",
    values_to = "Value"
  )
  
  # 3) Construct the ggplot with tile geometry
  plt <- ggplot(df_long, aes(x = Regressor, y = scan_number, fill = Value)) +
    geom_tile(...)
  
  # 4) Reverse the y-axis so scan_number=1 is at the top
  plt <- plt + scale_y_reverse()
  
  # 5) Select a color scale
  #    - If fill_midpoint != NULL => scale_fill_gradient2 with that midpoint
  #    - Otherwise => a 3-color gradient
  if (is.null(fill_midpoint)) {
    plt <- plt + scale_fill_gradientn(
      colours = c("navy", "white", "firebrick"),
      limits  = fill_limits
    )
  } else {
    plt <- plt + scale_fill_gradient2(
      midpoint = fill_midpoint,
      low      = "navy",
      mid      = "white",
      high     = "firebrick",
      limits   = fill_limits
    )
  }
  
  # 6) Optionally draw block-separators if we have block information
  #    The fmri_model inherits block info from x$event_model (and baseline).
  #    We'll just rely on x$event_model$blockids for run boundaries
  if (!is.null(x$event_model$blockids) && block_separators) {
    block_ids  <- x$event_model$blockids
    run_info   <- rle(block_ids)
    row_breaks <- cumsum(run_info$lengths)
    num_cols   <- ncol(DM)
    
    # Draw white lines at each boundary
    for (rb in row_breaks[-length(row_breaks)]) {
      plt <- plt + 
        annotate("segment",
                 x    = 0.5, 
                 xend = num_cols + 0.5,
                 y    = rb + 0.5,
                 yend = rb + 0.5,
                 color = "white", size = 1)
    }
  }
  
  # 7) Apply some theming
  plt <- plt + 
    theme_minimal(base_size = 14) +
    labs(x = "Regressors", y = "Scan Number", fill = "Value") +
    theme(
      panel.grid  = element_blank(),
      axis.text.x = if (rotate_x_text) element_text(angle = 45, hjust = 1) else element_text()
    )
  
  plt
}
