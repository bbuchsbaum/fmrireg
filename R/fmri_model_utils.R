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
  
  all_runs <- .fmri_model_run_ids(x)
  if (is.null(blocknum)) {
    blocknum <- all_runs
  }
  # Request the whole design without `blockid` when every run is selected:
  # subsetting by all runs is not equivalent to the full design for terms
  # that span runs (e.g. a global intercept), which some fmridesign versions
  # replicate once per requested run.
  dm_blockid <- .full_design_blockid(blocknum, all_runs)
  
  # Get the full convolved design matrix from the event model
  event_dm <- design_matrix(x$event_model, blockid = dm_blockid)
  
  # Get the baseline design matrix
  baseline_dm <- design_matrix(x$baseline_model, blockid = dm_blockid)
  
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
  bterms <- lapply(baseline_terms(x), function(term) as.matrix(design_matrix(term, blockid = dm_blockid)))
  
  # Compute indices for event and baseline terms
  num_event_cols <- ncol(event_dm)
  num_baseline_cols <- ncol(baseline_dm)
  
  # The term matrices are what the solvers stack, so they must account for
  # exactly the columns of the design they index.
  bterm_width <- sum(vapply(bterms, ncol, integer(1)))
  if (bterm_width != num_baseline_cols) {
    stop(sprintf(
      "Baseline term matrices have %d columns but the baseline design has %d; cannot index the design consistently.",
      bterm_width, num_baseline_cols
    ), call. = FALSE)
  }

  eterm_indices <- seq_len(num_event_cols)
  bterm_indices <- num_event_cols + seq_len(num_baseline_cols)
  
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

#' Sorted unique run identifiers of an fmri_model
#' @keywords internal
#' @noRd
.fmri_model_run_ids <- function(x) {
  # Runs come from the sampling frame (one id per scan). The event model's
  # blockids are per event and omit runs without events, so "every run"
  # judged from them would miss those runs and re-request the design by id.
  bids <- fmrihrf::blockids(x$event_model$sampling_frame)
  if (length(bids) == 0L) {
    bids <- x$event_model$blockids
  }
  sort(unique(bids))
}

#' Normalise a run selection for fmridesign's design_matrix()
#'
#' Returns `NULL` (meaning "the whole design") when `blockid` selects every
#' run, and `blockid` unchanged otherwise.
#' @keywords internal
#' @noRd
.full_design_blockid <- function(blockid, all_runs) {
  if (is.null(blockid)) {
    return(NULL)
  }
  if (setequal(blockid, all_runs)) NULL else blockid
}

