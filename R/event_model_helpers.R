# Internal Helper Functions for Event Model Construction

#' @importFrom tibble tibble
#' @importFrom rlang is_formula f_rhs f_lhs f_env eval_tidy syms sym is_list is_call call_name call_args %||% empty_env %||% is_expression as_label expr quo_get_expr env_bury
#' @import assertthat
#' @importFrom purrr map map_dfr
#' @importFrom cli cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom stats setNames approxfun quantile dnorm dgamma poly model.matrix terms reformulate
#' @importFrom utils head combn modifyList

# Helper to recursively find and evaluate hrf/trialwise calls in formula RHS
#' @keywords internal
#' @noRd
find_and_eval_hrf_calls <- function(expr, data, f_env) {
  
  if (rlang::is_call(expr)) {
    fun_name <- rlang::call_name(expr)
    
    # Base case: Found a target function call
    if (fun_name %in% c("hrf", "trialwise", "afni_hrf", "afni_trialwise", "covariate")) {
      # Evaluate the call in the formula env, using data as a mask
      eval_env <- rlang::env_bury(f_env %||% rlang::empty_env(), !!!data)
      evaluated_spec <- try(rlang::eval_tidy(expr, env = eval_env), silent = TRUE)
      
      if (inherits(evaluated_spec, "try-error")) {
          stop(sprintf("Failed to evaluate term '%s': %s", 
                       rlang::as_label(expr), attr(evaluated_spec, "condition")$message), call. = FALSE)
      }
      if (!inherits(evaluated_spec, "hrfspec")) {
          stop(sprintf("Term '%s' did not evaluate to an 'hrfspec' object (result class: %s)",
                       rlang::as_label(expr), class(evaluated_spec)[1]), call. = FALSE)
      }
      return(list(evaluated_spec))
    }
    
    # Recursive case: Traverse arguments of other calls (like + or *)
    # Collect results from each argument and concatenate while preserving
    # list structure (so each hrfspec remains intact)
    results_list <- lapply(rlang::call_args(expr), find_and_eval_hrf_calls, data, f_env)
    results <- do.call(c, results_list)
    return(results)
    
  } else {
    # Ignore symbols, constants, other non-calls that aren't target functions
    return(list())
  }
}


#' Parse fmri model formula (Internal Helper for EM-2)
#'
#' Takes a model formula (e.g., `onset ~ hrf(cond) + hrf(mod, basis="spmg3")`)
#' and extracts the onset variable, evaluates the `hrf()` calls on the RHS
#' to create `hrfspec` objects, and returns these components along with the
#' formula's environment.
#'
#' @param formula The model formula.
#' @param data The data frame providing context for evaluation.
#' @return A list containing: 
#'  `spec_tbl` (tibble with term_id, spec), 
#'  `onsets` (numeric vector), 
#'  `formula_env` (the formula environment).
#' @keywords internal
#' @noRd
parse_event_formula <- function(formula, data) {
  stopifnot(rlang::is_formula(formula))
  formula_env <- rlang::f_env(formula)
  
  # --- Extract and Evaluate LHS (Onsets) --- 
  lhs_expr <- rlang::f_lhs(formula)
  if (is.null(lhs_expr)) stop("Model formula must have a left-hand side (LHS) specifying the onset variable.", call.=FALSE)
  
  # Evaluate LHS in the data environment (masked)
  onsets <- try(rlang::eval_tidy(lhs_expr, data = data, env = formula_env %||% rlang::empty_env()), silent = TRUE)
  if (inherits(onsets, "try-error")) {
       stop(sprintf("Failed to evaluate onset variable '%s' from formula LHS: %s", 
                    rlang::as_label(lhs_expr), attr(onsets, "condition")$message), call. = FALSE)
  }
  if (!is.numeric(onsets)) {
       stop(sprintf("Onset variable '%s' extracted from formula LHS is not numeric (class: %s).", 
                    rlang::as_label(lhs_expr), class(onsets)[1]), call. = FALSE)
  }
  
  # --- Parse and Evaluate RHS (hrfspec terms) --- 
  rhs_expr <- rlang::f_rhs(formula)
  
  # Recursively find and evaluate hrfspec calls
  hrfspec_list <- find_and_eval_hrf_calls(rhs_expr, data, formula_env)
  
  if (length(hrfspec_list) == 0) {
      stop("No valid hrf(...) or similar terms found on the right-hand side of the formula.", call.=FALSE)
  }
  
  # Create the spec table
  spec_tbl <- tibble::tibble(
    term_id = seq_along(hrfspec_list),
    spec = hrfspec_list
  )
  
  # Return components
  list(
    spec_tbl = spec_tbl,
    onsets = onsets,
    formula_env = formula_env
  )
}

#' Construct an event_term object from an hrfspec (Internal Helper for EM-4)
#'
#' Takes a single hrf specification and the parsed model spec containing shared
#' timing info and evaluation context. Evaluates the variables/expressions 
#' in the hrfspec and calls the public `event_term` constructor.
#'
#' @param hrfspec An `hrfspec` object.
#' @param model_spec The list returned by `parse_event_model`.
#' @return An `event_term` object with the `hrfspec` attached as an attribute.
#' @keywords internal
#' @noRd
construct_event_term <- function(hrfspec, model_spec) {
  stopifnot(inherits(hrfspec, "hrfspec"))
  stopifnot(is.list(model_spec))
  
  vars_quos <- hrfspec$vars # List of quosures
  data_env <- model_spec$data
  f_env <- model_spec$formula_env %||% rlang::empty_env()
  
  # Create safe evaluation environment with onsets available for trialwise()
  eval_env <- rlang::env_bury(f_env, !!!data_env, onsets = model_spec$onsets)
  
  evaluated_vars <- list()
  var_names <- character(length(vars_quos))
  
  for (i in seq_along(vars_quos)) {
    quo <- vars_quos[[i]]
    expr <- rlang::quo_get_expr(quo)
    var_label <- rlang::as_label(expr)
    
    val <- try(rlang::eval_tidy(quo, data = data_env, env = eval_env), silent = TRUE)
    
    # Special handling for trialwise expressions - use eval() instead of eval_tidy()
    # This is needed because eval_tidy() doesn't properly handle the onsets variable
    # in the evaluation environment for expressions like .trial_factor(length(onsets))
    if (!inherits(val, 'try-error') && grepl("trial_factor", var_label)) {
      val <- try(eval(rlang::quo_get_expr(quo), envir = eval_env), silent = TRUE)
    }
    
    if (inherits(val, 'try-error')) {
        stop(sprintf("Failed to evaluate variable/expression '%s' for term '%s': %s", 
                     var_label, hrfspec$name %||% "UnnamedTerm", 
                     attr(val, "condition")$message), call. = FALSE)
    }
    
    # Store the evaluated value directly
    evaluated_vars[[i]] <- val
    
    # Determine the name to use for this variable within the event_term list
    # Use the basis object name if available, otherwise the label
    var_names[i] <- if (inherits(val, 'ParametricBasis') && !is.null(val$name)) {
      # Use the pre-defined name from the basis object (e.g., "z_RT", "poly_RT", "RT1_RT2")
      # This name might not be globally unique yet, but identifies the component
      val$name 
    } else {
      # Use the sanitized label for simple variables or basis without names
      make.names(var_label) 
    }
  }
  
  # Ensure unique names *within this specific term's evaluated variable list*
  # These names are internal to the event_term object and NOT the final design matrix colnames.
  names(evaluated_vars) <- make.unique(var_names)
  
  # Evaluate subset expression if present
  subset_result <- NULL
  if (!is.null(hrfspec$subset)) {
    subset_result <- try(rlang::eval_tidy(hrfspec$subset, data = data_env, env = eval_env), silent = TRUE)
    if (inherits(subset_result, 'try-error')) {
      stop(sprintf("Failed to evaluate subset expression for term '%s': %s", 
                   hrfspec$name %||% "UnnamedTerm", 
                   attr(subset_result, "condition")$message), call. = FALSE)
    }
    if (!is.logical(subset_result)) {
      stop(sprintf("Subset expression for term '%s' must evaluate to a logical vector, got: %s", 
                   hrfspec$name %||% "UnnamedTerm", class(subset_result)[1]), call. = FALSE)
    }
  }
  
  # Call the public event_term constructor
  # Use hrfspec durations if specified, otherwise fall back to model_spec durations
  durations_to_use <- if (!is.null(hrfspec$durations)) hrfspec$durations else model_spec$durations
  
  et <- event_term(evlist = evaluated_vars, 
                   onsets = model_spec$onsets, 
                   blockids = model_spec$blockids, 
                   durations = durations_to_use,
                   subset = subset_result)
                   
  attr(et, "hrfspec") <- hrfspec
  
  # TODO (EM-19): Add check here
  
  et
}

#' Build event model design matrix (Stage 3 of Pipeline)
#'
#' Takes the list of realised event terms, convolves each one, 
#' cbinds the results, and attaches term span/index attributes.
#'
#' @param terms A named list of `event_term` objects.
#' @param sampling_frame The fmri sampling_frame object.
#' @param precision Numeric precision for convolution.
#' @param parallel Logical indicating whether to use parallel processing (requires `future.apply`).
#' @return A tibble representing the full event model design matrix, with 
#'         `term_spans` and `col_indices` attributes.
#' @keywords internal
#' @noRd
build_event_model_design_matrix <- function(terms, sampling_frame, precision, parallel = FALSE) {

  if (length(terms) == 0) {
    warning("No terms provided to build design matrix.", call. = FALSE)
    return(tibble::tibble(.rows = sum(fmrihrf::blocklens(sampling_frame))))
  }

  # Convolution function for a single term
  convolve_one_term <- function(term) {
      # Check if this is a covariate term (already "convolved" - no HRF needed)
      if (inherits(term, "covariate_convolved_term")) {
          # Covariate terms don't need convolution, just return their design matrix
          dm <- as.matrix(term$design_matrix)
          # Ensure proper column names
          if (is.null(colnames(dm))) {
              colnames(dm) <- paste0("cov_", seq_len(ncol(dm)))
          }
          return(dm)
      }
      
      # Check if this is already a convolved term (AFNI only now)
      if (inherits(term, "afni_hrf_convolved_term") || inherits(term, "afni_trialwise_convolved_term")) {
          # AFNI terms are already "convolved" and should return NULL
          # since they are meant to be handled by AFNI's 3dDeconvolve, not R's convolution
          warning("AFNI terms are not fully supported in the current event_model pipeline. ", 
                  "They are intended for use with AFNI's 3dDeconvolve.", call. = FALSE)
          # Return NULL to indicate this term should not contribute columns to the design matrix
          return(NULL)
      }
      
      # For regular event_term objects, check for hrfspec and convolve
      hrfspec <- attr(term, "hrfspec")
      if (is.null(hrfspec) || is.null(hrfspec$hrf)) {
          stop(sprintf("Cannot convolve term '%s': missing attached hrfspec or HRF object within spec.", 
                       attr(term, "term_tag") %||% term$varname %||% "Unnamed"), call. = FALSE)
      }
      
      # Regular convolution for event_term objects
      convolve.event_term(term, hrf = hrfspec$hrf, sampling_frame = sampling_frame, 
                          precision = precision, drop.empty = TRUE)
  }

  # Apply convolution (potentially in parallel)
  term_matrices <- if (parallel) {
      warning("Parallel processing not fully implemented yet, using sequential lapply.", call. = FALSE)
      lapply(terms, convolve_one_term)
  } else {
      lapply(terms, convolve_one_term)
  }

  # Filter out NULL matrices (from AFNI terms that don't contribute to design matrix)
  non_null_indices <- which(!sapply(term_matrices, is.null))
  term_matrices_filtered <- term_matrices[non_null_indices]
  terms_filtered <- terms[non_null_indices]
  
  # Handle case where all terms are NULL (all AFNI terms)
  if (length(term_matrices_filtered) == 0) {
    warning("All terms are AFNI terms that don't contribute to the design matrix. Returning empty design matrix.", call. = FALSE)
    empty_dm <- tibble::tibble(.rows = sum(fmrihrf::blocklens(sampling_frame)))
    attr(empty_dm, "term_spans") <- integer(0)
    attr(empty_dm, "col_indices") <- list()
    attr(empty_dm, "colnames_final") <- TRUE
    return(empty_dm)
  }

  # Calculate col_indices based on term_ncols before combining matrices
  term_ncols <- vapply(term_matrices_filtered, ncol, integer(1))
  term_spans <- cumsum(term_ncols)
  col_indices <- vector("list", length(terms))
  start_idx <- 1
  for (i in seq_along(terms)) {
      if (i %in% non_null_indices) {
        # Find position in filtered list
        filtered_pos <- which(non_null_indices == i)
        end_idx <- if (filtered_pos == 1) term_ncols[1] else term_spans[filtered_pos]
        col_indices[[i]] <- seq.int(start_idx, end_idx)
        start_idx <- end_idx + 1
      } else {
        # AFNI term - no columns
        col_indices[[i]] <- integer(0)
      }
  }
  names(col_indices) <- names(terms) # Names are now the term_tags

  # --- Combine convolved term matrices --- 
  # Store original names before potential modification by cbind/tibble
  original_colnames_list <- lapply(term_matrices_filtered, colnames)
  # Perform cbind
  dm_mat <- do.call(cbind, term_matrices_filtered)
  
  # Assign combined list of original names initially proposed by convolve.event_term
  current_colnames <- unlist(original_colnames_list)
  
  # Ensure global uniqueness using make.names before tibble conversion
  # This respects the naming_refactoring.md guarantee for unique final names.
  unique_colnames <- make.names(current_colnames, unique = TRUE)
  
  # Warn if names were actually changed by make.names (indicates initial clash)
  if (!identical(current_colnames, unique_colnames)) {
      changed_indices <- which(current_colnames != unique_colnames)
      
      # Filter out changes that are just empty strings being converted to valid names
      # These are not real duplicates, just invalid names being fixed
      real_duplicate_indices <- changed_indices[current_colnames[changed_indices] != ""]
      
      # Only warn if there are actual meaningful duplicates (not just empty string fixes)
      if (length(real_duplicate_indices) > 0) {
          # Provide a more informative warning about which names were auto-suffixed
          warning_message <- paste0(
              "Duplicate column names detected and automatically resolved by appending suffixes (e.g., '.1', '.2'). ",
              "This typically occurs when multiple terms generate identical column name proposals (e.g., from multiple hrf(Ident(X,Y)) calls where 'X' is repeated, or if term_tags clash in unforeseen ways).",
              " Review model specification if this is unexpected. Changed names include:\n",
              paste0("  Original: '", current_colnames[real_duplicate_indices], "' -> Became: '", unique_colnames[real_duplicate_indices], "'", collapse = "\n")
          )
          warning(warning_message, call. = FALSE)
      }
  }
  colnames(dm_mat) <- unique_colnames
  
  # Convert to tibble - .name_repair = "minimal" is now safe as names are unique and valid R names.
  dm <- tibble::as_tibble(dm_mat, .name_repair = "minimal") 

  # Check column count consistency
  if (ncol(dm) != term_spans[length(term_spans)]) {
      warning("Column mismatch after combining term matrices. Check names/cbind.", call. = FALSE)
  }
  
  # --- Attach Attributes --- 
  # Attach term_spans: vector for formula interface, matrix for list interface
  interface_type <- attr(terms, "interface") %||% "formula"
  if (interface_type == "list") {
    starts <- c(1, head(term_spans, -1) + 1)
    spans_mat <- cbind(starts, term_spans)
    attr(dm, "term_spans") <- spans_mat
  } else {
    attr(dm, "term_spans") <- term_spans
  }
  attr(dm, "col_indices") <- col_indices
  # Add flag indicating names are final according to new grammar
  attr(dm, "colnames_final") <- TRUE
  
  dm
}


#' Realise event terms from parsed specification (Stage 2 of Pipeline)
#'
#' Takes the parsed model specification and iterates through each term spec,
#' calling `construct_event_term` to create the actual `event_term` objects.
#' Handles optional progress display.
#'
#' @param parsed_spec Output list from `parse_event_model`.
#' @param sampling_frame The fmri sampling_frame object.
#' @param drop_empty Logical indicating whether to drop empty events (passed to term construction).
#' @param progress Logical indicating whether to show a progress bar.
#' @return A named list of `event_term` objects.
#' @keywords internal
#' @noRd
realise_event_terms <- function(parsed_spec, sampling_frame, drop_empty = TRUE, progress = FALSE) {

  n_terms <- nrow(parsed_spec$spec_tbl)
  term_specs <- parsed_spec$spec_tbl$spec

  if (progress) {
    pb <- cli::cli_progress_bar("Realising event terms", total = n_terms, clear = FALSE)
    on.exit(cli::cli_progress_done(id = pb), add = TRUE)
  }

  # Create a complete model_spec that includes sampling_frame for construct methods
  model_spec <- parsed_spec
  model_spec$sampling_frame <- sampling_frame

  term_list <- vector("list", n_terms)
  term_tags <- character(n_terms)
  term_uids <- character(n_terms)

  for (i in 1:n_terms) {
    spec <- term_specs[[i]]

    # 1. Generate UID (t01, t02, ...)
    term_uid <- sprintf("t%02d", i)
    spec$uid <- term_uid # Store on spec
    term_uids[i] <- term_uid # Collect UIDs

    # 2. Generate Term Tag using helper
    # Pass *previously generated* tags for uniqueness check (use actual tags collected so far)
    # Note: existing_tags should contain the actual tags, not placeholders like NA
    term_tag_value <- make_term_tag(spec, existing_tags = stats::na.omit(term_tags[seq_len(i - 1)]))
    # ^^^ Get the value (could be NULL or a character tag)

    spec$term_tag <- term_tag_value # Store actual value (NULL or tag) on spec

    # Store a placeholder (NA) in the character vector if tag is NULL,
    # otherwise store the actual tag.
    term_tags[i] <- if (is.null(term_tag_value)) NA_character_ else term_tag_value
    
    # 3. Construct the event_term object
    # Use the generic construct method which will dispatch to the appropriate method
    # (construct.afni_hrfspec, construct.hrfspec, etc.) based on the spec type
    # Pass the complete model_spec that includes sampling_frame
    term_list[[i]] <- construct(spec, model_spec)

    # 4. Store UID and *actual* Term Tag (NULL or character) attribute on the event_term object itself
    if (!is.null(term_list[[i]])) {
        attr(term_list[[i]], "uid") <- term_uid
        attr(term_list[[i]], "term_tag") <- term_tag_value # Assign NULL or tag here
    } else {
        # This might happen if construct_event_term fails or returns NULL
        warning(paste("Term", i, "construction resulted in NULL. Cannot set attributes."))
    }

    # Update progress
    if (progress) {
      cli::cli_progress_update(id = pb)
    }
  }

  # --- Post-processing --- 

  # Assign names to the list, handling NAs from NULL term tags.
  final_list_names <- term_tags
  na_indices <- which(is.na(final_list_names))
  if (length(na_indices) > 0) {
      # Generate default names for NA elements, e.g., based on index or variable names
      for(idx in na_indices) {
         # Try to get a fallback name from the term object's varname attribute if available
         term_obj <- term_list[[idx]]
         fallback_name <- if (!is.null(term_obj) && !is.null(term_obj$varname)) {
                              term_obj$varname # Use varname (e.g., "Ident.RT1..RT2.") as list name
                          } else {
                              paste0("term_", idx) # Generic fallback
                          } 
         final_list_names[idx] <- fallback_name
      }
      message(sprintf("Generated default list names for %d terms with NULL term_tags (likely Ident).", length(na_indices)))
  }
  
  # Ensure final list names are unique
  names(term_list) <- make.unique(final_list_names) 
  
  # Tag input interface type for later builder
  attr(term_list, "interface") <- parsed_spec$interface

  term_list
}

#' Parse event model inputs (Stage 1 of Pipeline)
#'
#' Internal helper to normalize formula/list inputs, validate/extract shared
#' onset, block, duration info, and create a standardized spec table.
#'
#' @param formula_or_list Formula or list of hrfspec objects.
#' @param data The data frame containing event variables.
#' @param block Formula or vector specifying blocks/runs.
#' @param durations Numeric vector or scalar for event durations.
#' @return A list containing: spec_tbl, onsets, durations, blockids, data, formula_env.
#' @keywords internal
#' @noRd
parse_event_model <- function(formula_or_list, data, block, durations = 0) {
  
  stopifnot(inherits(data, "data.frame"))
  
  spec_tbl <- NULL
  onsets <- NULL
  formula_env <- NULL
  
  if (rlang::is_formula(formula_or_list)) {
    parsed_form <- parse_event_formula(formula_or_list, data)
    spec_tbl <- parsed_form$spec_tbl
    onsets <- parsed_form$onsets
    formula_env <- parsed_form$formula_env
    stopifnot(is.numeric(onsets))
    if(length(onsets) != nrow(data)) {
         stop(sprintf("Length of extracted onset variable (%d) != nrow(data) (%d)",
                      length(onsets), nrow(data)), call.=FALSE)
    }
    
  } else if (is.list(formula_or_list)) {
    stopifnot(all(vapply(formula_or_list, inherits, TRUE, "hrfspec")), 
              length(formula_or_list) > 0)
              
    spec_tbl <- tibble::tibble(
        term_id = seq_along(formula_or_list),
        spec = formula_or_list
    )
    
    if (!"onset" %in% names(data)) {
        stop("When providing a list of hrfspecs, `data` must contain an 'onset' column. If named differently, please use the formula interface.", call.=FALSE)
    }
    onsets <- data$onset
    stopifnot(is.numeric(onsets))
    if(length(onsets) != nrow(data)) {
         stop(sprintf("Length of 'onset' column (%d) != nrow(data) (%d)",
                      length(onsets), nrow(data)), call.=FALSE)
    }
    formula_env <- rlang::empty_env() 
    
  } else {
    stop("`formula_or_list` must be a formula or a list of hrfspec objects.", call.=FALSE)
  }
  
  stopifnot(all(!is.na(onsets)))

  blockids_raw <- NULL
  if (rlang::is_formula(block)) {
    block_rhs <- rlang::f_rhs(block)
    if (identical(block_rhs, 1L) || identical(block_rhs, 1)) {
      blockids_raw <- rep(1, nrow(data))
    } else {
      blockids_raw <- rlang::eval_tidy(block_rhs, data)
    }
  } else {
    blockids_raw <- block
  }
  
  if (is.factor(blockids_raw)) {
    blockids_raw <- as.integer(as.character(blockids_raw))
  }
  
  stopifnot(length(blockids_raw) == length(onsets),
            all(!is.na(blockids_raw))) 
            
  # Use stopifnot for non-decreasing check for better trace
  if (is.unsorted(blockids_raw, strictly = FALSE)) {
       stop("'blockids' must be non-decreasing.", call. = FALSE)
  }
  
  # Canonicalize blockids
  blockids <- as.integer(match(blockids_raw, unique(blockids_raw)))
  
  if (is.null(durations)) durations <- 0 
  durations_proc <- recycle_or_error(durations, length(onsets), "durations")
  stopifnot(all(!is.na(durations_proc)))

  list(
    spec_tbl    = spec_tbl,
    onsets      = onsets,
    durations   = durations_proc,
    blockids    = blockids, # Use canonicalized IDs
    data        = data, 
    formula_env = formula_env,
    interface   = if (rlang::is_formula(formula_or_list)) "formula" else "list"
  )
} 