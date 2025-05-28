#' @title Model Utilities for fmri_model Objects
#' @description Functions for creating and manipulating fmri_model objects
#' @keywords internal
#' @importFrom assertthat assert_that
#' @importFrom stats as.formula

#' Get the formula representation of an fMRI model
#'
#' This function extracts the formula from an \code{fmri_model} object.
#'
#' @return A formula representing the model.
#' @param x An \code{fmri_model} object.
#' @param ... Additional arguments.
#' @rdname get_formula
#' @export
get_formula.fmri_model <- function(x,...) {
  assert_that(inherits(x, "fmri_model"), msg = "'x' must be an 'fmri_model' object")
  term_names <- names(terms(x))
  form <- paste(".y ~", paste(term_names, collapse = " + "), "-1")
  return(as.formula(form))
}

#' Extract Term Matrices from an fMRI Model
#'
#' This function extracts the term matrices from an \code{fmri_model}, which consists of event-related terms
#' and baseline-related terms. These matrices are used to build the design matrix in fMRI data analysis.
#'
#' @param x An \code{fmri_model} object containing the event and baseline models.
#' @param blocknum (Optional) A numeric vector specifying the block numbers to include. Defaults to all blocks.
#' @param ... Additional arguments
#' @return A named list of term matrices, with event terms followed by baseline terms.
#'         Attributes \code{"event_term_indices"} and \code{"baseline_term_indices"} store the indices of event and baseline terms,
#'         \code{"blocknum"} stores the block numbers, and \code{"varnames"} stores the variable names.
#' @export
#' @rdname term_matrices
term_matrices.fmri_model <- function(x, blocknum = NULL,...) {
  assert_that(inherits(x, "fmri_model"), msg = "'x' must be an 'fmri_model' object")
  
  if (is.null(blocknum)) {
    blocknum <- sort(unique(x$event_model$blockids))
  }
  
  # Get the full convolved design matrix from the event model
  event_dm <- design_matrix(x$event_model, blockid = blocknum)
  
  # Get the baseline design matrix
  baseline_dm <- design_matrix(x$baseline_model, blockid = blocknum)
  
  # Extract individual term matrices using the col_indices attribute
  col_indices <- attr(x$event_model$design_matrix, "col_indices")
  if (is.null(col_indices)) {
    stop("Event model design matrix missing 'col_indices' attribute needed to extract individual term matrices.")
  }
  
  # Extract event term matrices from the full convolved design matrix
  eterms <- lapply(names(col_indices), function(term_name) {
    indices <- col_indices[[term_name]]
    as.matrix(event_dm[, indices, drop = FALSE])
  })
  names(eterms) <- names(col_indices)
  
  # Extract baseline term matrices (baseline terms are simpler, one per term)
  bterms <- lapply(baseline_terms(x), function(term) as.matrix(design_matrix(term, blockid = blocknum)))
  
  # Compute indices for event and baseline terms
  num_event_cols <- ncol(event_dm)
  num_baseline_cols <- ncol(baseline_dm)
  
  eterm_indices <- 1:num_event_cols
  bterm_indices <- (num_event_cols + 1):(num_event_cols + num_baseline_cols)
  
  # Combine term matrices
  term_matrices <- c(eterms, bterms)
  names(term_matrices) <- names(terms(x))
  
  # Collect variable names
  vnames <- c(colnames(event_dm), colnames(baseline_dm))
  
  # Set attributes
  attr(term_matrices, "event_term_indices") <- eterm_indices
  attr(term_matrices, "baseline_term_indices") <- bterm_indices
  attr(term_matrices, "blocknum") <- blocknum
  attr(term_matrices, "varnames") <- vnames
  
  return(term_matrices)
}

#' Create an fmri_model from a Formula
#'
#' @description
#' This function creates an \code{fmri_model} object from a formula, block specification,
#' and dataset. It's a convenience function that combines event and baseline models.
#'
#' @param formula A formula specifying the event model
#' @param block A factor vector indicating the block structure
#' @param baseline_model An optional baseline model. If NULL, a default polynomial baseline is created
#' @param dataset The fmri_dataset object
#' @param drop_empty Logical, whether to drop empty blocks
#' @param durations Duration of events (default 0 for instantaneous events)
#' @return An fmri_model object
#' @export
create_fmri_model <- function(formula, block, baseline_model = NULL, dataset, drop_empty = TRUE, durations = 0) {
  # Handle block variable
  if (is.character(block)) {
    block <- factor(block)
  }
  
  if (!is.factor(block)) {
    block <- factor(block)
  }
  
  # Drop empty levels if requested
  if (drop_empty) {
    block <- droplevels(block)
  }
  
  # Create sampling frame
  sframe <- sampling_frame(block, dataset$TR)
  
  # Create baseline model if not provided
  if (is.null(baseline_model)) {
    baseline_model <- baseline_model(sframe, degree = 3, intercept = TRUE)
  }
  
  # Create event model
  event_model <- event_model(formula, block = sframe, durations = durations)
  
  # Combine into fmri_model
  fmri_model(event_model, baseline_model)
}