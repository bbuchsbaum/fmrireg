#' @title Model Utilities for fmri_model Objects
#' @description Functions for creating and manipulating fmri_model objects
#' @keywords internal
#' @importFrom assertthat assert_that
#' @importFrom stats as.formula

#' @method get_formula fmri_model
#' @rdname get_formula
#' @examples
#' fm <- fmrireg:::.demo_fmri_model()
#' get_formula(fm)
#' @export
get_formula.fmri_model <- function(x,...) {
  assert_that(inherits(x, "fmri_model"), msg = "'x' must be an 'fmri_model' object")
  term_names <- names(terms(x))
  form <- paste(".y ~", paste(term_names, collapse = " + "), "-1")
  return(as.formula(form))
}

#' Extract term matrices from fmri_model
#'
#' Extract design matrices for individual terms from an fmri_model object.
#'
#' @param x An fmri_model object
#' @param blocknum Optional vector of block numbers to extract matrices for
#' @param ... Additional arguments (currently unused)
#' @return A list of matrices, one for each term in the model
#' @examples
#' fm <- fmrireg:::.demo_fmri_model()
#' term_matrices(fm)
#' @method term_matrices fmri_model
#' @export
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
#' @return An fmri_model object that also stores the \code{dataset}
#' @export
create_fmri_model <- function(formula, block, baseline_model = NULL, dataset, drop_empty = TRUE, durations = 0) {
  # Handle block variable
  formula <- .fmrireg_inject_registered_bases(formula)
  if (inherits(block, "formula")) {
    # Evaluate formula against event_table
    block_var <- all.vars(block)
    if (length(block_var) != 1) {
      stop("Block formula must specify exactly one variable")
    }
    if (!block_var %in% names(dataset$event_table)) {
      stop(sprintf("Block variable '%s' not found in event_table", block_var))
    }
    block <- dataset$event_table[[block_var]]
  }
  
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
  
  # Create sampling frame using run_length from dataset
  if (!is.null(dataset$sampling_frame)) {
    sframe <- dataset$sampling_frame
  } else {
    # Fallback to using run_length if available
    if (!is.null(dataset$run_length)) {
      sframe <- sampling_frame(dataset$run_length, dataset$TR)
    } else {
      stop("Dataset must have either sampling_frame or run_length")
    }
  }
  
  # Create baseline model if not provided
  if (is.null(baseline_model)) {
    baseline_model <- baseline_model(basis = "poly", degree = 3, sframe = sframe, intercept = "global")
  }
  
  # Create event model
  # For event_model, we need to pass the event table data along with block assignments
  event_data <- dataset$event_table
  event_model <- event_model(formula, data = event_data, block = block, 
                           sampling_frame = sframe, durations = durations)
  
  # Combine into fmri_model
  fmri_model(event_model, baseline_model, dataset)
}
