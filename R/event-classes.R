#' @include utils-internal.R
# Internal unified constructor for event sequences
# 
# Creates an 'event-sequence' object (many rows, one variable).
# This is the core internal function that standardizes various input types 
# into a common representation: a list with class c("event", "event_seq").
# It stores the event payload in `$value` (always a numeric matrix) and 
# additional metadata (like factor levels or basis objects) in `$meta`.
# Public-facing wrappers like `event_factor`, `event_variable` call this function.
#
# @param value Event values (factor, character, numeric vector, matrix, ParametricBasis).
# @param name Name of the event variable (will be sanitized).
# @param onsets Numeric vector of event onsets (seconds).
# @param blockids Numeric vector of block IDs (non-decreasing integers).
# @param durations Numeric vector of event durations (seconds), or a scalar.
# @param subset Optional logical vector indicating which events to keep. If
#   provided, the vector must be the same length as `onsets` and contain no `NA`
#   values.
# 
# @return An S3 object of class `event` and `event_seq`, containing:
#   \item{varname}{Sanitized variable name.} 
#   \item{onsets}{Numeric vector of onsets (subsetted).} 
#   \item{durations}{Numeric vector of durations (subsetted and recycled).}
#   \item{blockids}{Numeric vector of block IDs (subsetted).}
#   \item{value}{Numeric matrix storing the event payload. 
#                 For factors, this contains integer codes (1-based).
#                 For basis objects, this is the basis matrix (`basis$y`).
#                 For numeric vectors/matrices, it's the numeric data.} 
#   \item{continuous}{Logical flag. TRUE for numeric, matrix, basis; FALSE for factor/character.} 
#   \item{meta}{A list containing metadata, or NULL. 
#               For factors: `meta = list(levels = c("level1", "level2", ...))`. 
#               For basis objects: `meta = list(basis = <ParametricBasis object>)`.} 
# 
# @name event
# @keywords internal
# @noRd
event <- function(value, name, onsets, blockids, durations = 0, subset = NULL) {
  
  # --- Add NA check for input value --- 
  if (!inherits(value, "ParametricBasis") && anyNA(value)) { # Don't check inside basis objects
      warning(sprintf("NA values detected in input `value` for event variable '%s'. Coercion or filtering might occur.", name), call. = FALSE)
  }
  # -------------------------------------
  
  n_initial <- length(onsets)
  if (!is.null(subset)) {
      assertthat::assert_that(length(subset) == n_initial,
        msg = sprintf("subset length (%d) must match onsets length (%d)",
                      length(subset), n_initial))
      assertthat::assert_that(!anyNA(subset),
        msg = "subset cannot contain NA values")
  } else {
      subset <- rep(TRUE, n_initial)
  }
  
  # --------------------------------------------
  # 1. Validate / recycle / pre-process with the existing helper
  # TODO: Review/refactor/integrate .checkEVArgs itself later (Ticket EV-4)
  meta_checked <- .checkEVArgs(name, value, onsets, blockids, durations)
  
  # Overwrite with validated/processed objects from the helper
  name      <- meta_checked$varname # Use sanitized name
  value     <- meta_checked$value
  onsets    <- meta_checked$onsets
  blockids  <- meta_checked$blockids
  durations <- meta_checked$durations # Recycled by helper
  # --------------------------------------------
  
  # Apply subsetting safely AFTER validation and recycling
  keep_indices <- which(subset)
  onsets_sub    <- onsets[keep_indices]
  blockids_sub  <- blockids[keep_indices]
  durations_sub <- durations[keep_indices]
  n_final <- length(onsets_sub)
  
  # Determine continuous flag based on *original* input type
  # Basis, numeric, matrix are continuous; factor, character are not.
  is_continuous <- !(is.factor(value) || is.character(value))
  
  # Process value based on type, applying subsetting.
  # Use if/else if/else for conditional logic
  
  if (inherits(value, "ParametricBasis")) {
      basis_subset <- sub_basis(value, subset)
      switch_result <- list(value = basis_subset$y, meta = list(basis = basis_subset))
  } else if (is.factor(value) || is.character(value)) {
      # Store original levels - handle both factors and character vectors
      original_levels <- if (is.factor(value)) {
          levels(value)
      } else {
          # For character vectors, get unique values as levels
          unique(as.character(value))
      }
      # Create factor 'f' using subsetted values BUT preserving original levels
      f <- factor(as.character(value)[keep_indices], levels = original_levels)
      if (length(f) == 0) {
          val_payload <- matrix(numeric(0), nrow = 0, ncol = 1L)
      } else {
          val_payload <- matrix(as.integer(f), ncol = 1L)
      }
      # Store the ORIGINAL levels in meta$levels
      switch_result <- list(value = val_payload, meta = list(levels = original_levels))
  } else if (is.matrix(value)) {
      switch_result <- list(value = value[keep_indices, , drop = FALSE], meta = NULL)
  } else if (is.numeric(value) && (is.vector(value) || length(dim(value)) == 1)) {
      switch_result <- list(value = matrix(value[keep_indices], ncol = 1L), meta = NULL)
  } else {
      # Default case (should be caught by .checkEVArgs, but defensive)
      stop(paste("Unsupported `value` type:", class(value)[1], "for variable:", name))
  }
  
  # Assign results from the conditional logic
  val_mat   <- switch_result$value
  meta_list <- switch_result$meta

  # Ensure payload is a matrix
  if (!is.matrix(val_mat)) {
      stop("Internal error: Processed event value is not a matrix.") 
  }
  
  # Ensure colnames exist
  if (is.null(colnames(val_mat))) {
      if (ncol(val_mat) == 1 && !is.null(name)) {
          colnames(val_mat) <- name # Use varname for single column
      } else {
           # Use simple sequence for multi-column fallback
          colnames(val_mat) <- paste0("V", seq_len(ncol(val_mat)))
      } 
  }
  
  # Final length check (simplified as .checkEVArgs ensures alignment before subset)
  stopifnot(nrow(val_mat) == n_final)
  
  # Construct the final list structure
  structure(list(
      varname    = name, 
      onsets     = onsets_sub, 
      durations  = durations_sub, 
      blockids   = blockids_sub,
      value      = val_mat,        # Standardized matrix payload
      continuous = is_continuous,  # Robust flag based on input type
      meta       = meta_list       # List containing levels or basis object, or NULL
      ),
    class = c("event", "event_seq")
  )
} 


## ============================================================================
## Section: Public Event Constructor Wrappers
## ============================================================================

#' Create a categorical event sequence from a factor.
#' 
#' This is a user-facing wrapper around the internal `event()` constructor,
#' specifically for creating categorical event sequences from factors or characters.
#'
#' @param fac A factor or something coercible to a factor.
#' @param name Name of the event variable.
#' @param onsets Numeric vector of event onsets (seconds).
#' @param blockids Numeric vector of block IDs.
#' @param durations Numeric vector of event durations (seconds), or a scalar.
#' @param subset Optional logical vector indicating which events to keep. If
#'   provided, the vector must match `onsets` in length and contain no `NA`
#'   values.
#'
#' Column names are sanitized using `.sanitizeName()` if provided. If
#' column names are missing or not unique, deterministic feature suffixes
#' (`f01`, `f02`, ...) are generated instead. The resulting names are
#' returned by `levels()` for the event object.
#'
#' @return An S3 object of class `event` and `event_seq`.
#'
#' @examples
#' efac <- event_factor(factor(c("a", "b", "c", "a", "b", "c")), "abc", 
#'         onsets = seq(1, 100, length.out = 6))
#' print(efac)
#' levels(efac)
#'
#' @seealso \code{\link{event_model}}, \code{\link{event_variable}}, \code{\link{event_matrix}}, \code{\link{event_basis}}
#' @export 
event_factor <- function(fac, name, onsets, blockids = 1, durations = 0, subset = NULL) {
  # Convert to factor early if needed, but event() handles characters too.
  if (!is.factor(fac) && !is.character(fac)) {
    warning("Input `fac` is not a factor or character, attempting to convert.")
    fac <- factor(as.character(fac))
  }
  
  # Call the unified internal constructor
  event(value = fac, 
        name = name, 
        onsets = onsets, 
        blockids = blockids, 
        durations = durations, 
        subset = subset)
}        

#' Create a continuous event sequence from a numeric vector.
#'
#' This is a user-facing wrapper around the internal `event()` constructor,
#' specifically for creating continuous event sequences from numeric vectors.
#' 
#' @param vec Numeric vector representing continuous event values.
#' @param name Name of the event variable.
#' @param onsets Numeric vector of event onsets (seconds).
#' @param blockids Numeric vector of block IDs.
#' @param durations Numeric vector of event durations (seconds), or a scalar.
#' @param subset Optional logical vector indicating which events to keep. If
#'   provided, the vector must match `onsets` in length and contain no `NA`
#'   values.
#'
#' @return An S3 object of class `event` and `event_seq`.
#'
#' @examples
#' evar <- event_variable(c(1, 2, 3, 4, 5, 6), "example_var", onsets = seq(1, 100, length.out = 6))
#' print(evar)
#' is_continuous(evar)
#'
#' @seealso \code{\link{event_factor}}
#' @export
event_variable <- function(vec, name, onsets, blockids = 1, durations = 0, subset = NULL) {
  if (!is.numeric(vec) || !(is.vector(vec) || length(dim(vec)) <= 1)) {
      stop("`vec` must be a numeric vector.")
  }
  if (is.factor(vec)) {
    stop("Cannot create an event_variable from a factor, use 'event_factor'.")
  }

  # Call the unified internal constructor
  event(value = vec, 
        name = name, 
        onsets = onsets, 
        blockids = blockids, 
        durations = durations, 
        subset = subset)
}       

#' Create a continuous event set from a matrix.
#'
#' This is a user-facing wrapper around the internal `event()` constructor,
#' specifically for creating continuous event sequences from numeric matrices.
#' 
#' @param mat A numeric matrix of continuous event values (one row per event).
#' @param name Name of the event variable.
#' @param onsets Numeric vector of event onsets (seconds).
#' @param blockids Numeric vector of block IDs.
#' @param durations Numeric vector of event durations (seconds), or a scalar.
#' @param subset Optional logical vector indicating which events to keep. If
#'   provided, the vector must match `onsets` in length and contain no `NA`
#'   values.
#'
#' If `mat` has column names and more than one column, those names are
#' sanitized using `.sanitizeName()` before being stored. The sanitized
#' column names are returned by `levels()` for the resulting event object.
#'
#' @return An S3 object of class `event` and `event_seq`.
#'
#' @examples
#' mat <- matrix(rnorm(20), 10, 2, dimnames=list(NULL, c("Val1", "Val2")))
#' onsets <- seq(1, 100, length.out = 10)
#' durations <- rep(1, 10)
#' blockids <- rep(1, 10)
#' eset <- event_matrix(mat, "eset", onsets, blockids, durations)
#' print(eset)
#' columns(eset) # Alias for levels
#'
#' @export
event_matrix <- function(mat, name, onsets, blockids = 1, durations = 0, subset = NULL) {
  if (!is.matrix(mat) || !is.numeric(mat)) {
      stop("`mat` must be a numeric matrix.")
  }
  assert_that(nrow(mat) == length(onsets),
              msg = sprintf("Length mismatch for '%s': nrow(mat)=%d, length(onsets)=%d",
                          name, nrow(mat), length(onsets)))
  # Ensure deterministic column names
  cn <- colnames(mat)
  if (is.null(cn) || anyDuplicated(cn)) {
      colnames(mat) <- feature_suffix(seq_len(ncol(mat)), ncol(mat))
  } else {
      colnames(mat) <- .sanitizeName(cn)
  }
  
  # Call the unified internal constructor
  event(value = mat, 
        name = name, 
        onsets = onsets, 
        blockids = blockids, 
        durations = durations, 
        subset = subset)
}

#' Create an event set from a ParametricBasis object.
#'
#' This is a user-facing wrapper around the internal `event()` constructor,
#' specifically for creating event sequences modulated by a basis set.
#' 
#' @param basis A `ParametricBasis` object (e.g., from `BSpline`, `PolynomialBasis`).
#' @param name Optional name for the event variable. If NULL, uses `basis$name`.
#' @param onsets Numeric vector of event onsets (seconds).
#' @param blockids Numeric vector of block IDs.
#' @param durations Numeric vector of event durations (seconds), or a scalar.
#' @param subset Optional logical vector indicating which events to keep. If
#'   provided, the vector must match `onsets` in length and contain no `NA`
#'   values.
#'
#' @return An S3 object of class `event` and `event_seq`.
#'
#' @import assertthat
#' @examples
#' basis <- BSpline(1:21, 3)
#' onsets <- seq(0, 20, length.out = 21)
#' blockids <- rep(1, length(onsets))
#' ebasis <- event_basis(basis, onsets=onsets, blockids=blockids)
#' print(ebasis)
#' levels(ebasis)
#'
#' @export
event_basis <- function(basis, name = NULL, onsets, blockids = 1, durations = 0, subset = NULL) {
  assertthat::assert_that(inherits(basis, "ParametricBasis"))
  
  # Use basis$name if name is not provided
  if (is.null(name)) {
      name <- basis$name
      if (is.null(name)) { # Fallback if basis name is also NULL
          warning("No name provided and basis$name is NULL, using default name 'basis'")
          name <- "basis"
      }
  }
  
  # Call the unified internal constructor
  event(value = basis, 
        name = name, 
        onsets = onsets, 
        blockids = blockids, 
        durations = durations, 
        subset = subset)
} 

## ============================================================================
## Section: Unified S3 Methods for 'event' class
## ============================================================================

#' Extract Levels from fmrireg Objects
#'
#' @description
#' Extract levels from various fmrireg objects. These methods extend the base R 
#' \code{\link[base]{levels}} generic to work with fmrireg-specific classes.
#'
#' @param x An object from which to extract levels. Can be:
#'   \itemize{
#'     \item An \code{event} object - returns factor levels or column names
#'     \item A \code{Scale} object - returns the variable name
#'     \item A \code{ScaleWithin} object - returns the variable name
#'     \item A \code{RobustScale} object - returns the variable name
#'   }
#' @return A character vector of levels or names, depending on the object type:
#'   \itemize{
#'     \item For categorical events: the factor levels
#'     \item For continuous events: the column names (matrices) or variable name (vectors)
#'     \item For scale objects: the variable name being scaled
#'   }
#' @examples
#' # Create a categorical event
#' fac_event <- event_factor(
#'   factor(c("A", "B", "A", "B")),
#'   name = "condition",
#'   onsets = c(1, 10, 20, 30),
#'   blockids = rep(1, 4)
#' )
#' levels(fac_event)  # Returns: c("A", "B")
#' 
#' # Create a continuous event
#' cont_event <- event_variable(
#'   c(1.2, 0.8, 1.5, 0.9),
#'   name = "reaction_time",
#'   onsets = c(1, 10, 20, 30),
#'   blockids = rep(1, 4)
#' )
#' levels(cont_event)  # Returns: "reaction_time"
#' @name levels.event
#' @aliases levels.Scale levels.ScaleWithin levels.RobustScale
#' @method levels event
#' @export
levels.event <- function(x) {
  if (!is.null(x$meta$basis)) {
    # Use levels method for the basis object itself
    levels(x$meta$basis) 
  } else if (!is.null(x$meta$levels)) {
    # Use stored levels for factors
    x$meta$levels
  } else {
    # Default to column names of the value matrix
    colnames(x$value)
  }
}

#' @describeIn levels.event Alias for levels.event
#' @export
columns.event <- levels.event

#' @method is_continuous event
#' @rdname is_continuous
#' @export
is_continuous.event <- function(x) {
  x$continuous
}

#' @method is_categorical event
#' @rdname is_categorical
#' @export
is_categorical.event <- function(x) {
  !is_continuous(x)
}

#' @method cells event
#' @rdname cells
#' @export
#' @importFrom tibble tibble
cells.event <- function(x, drop.empty = TRUE, ...) {
  var_name <- x$varname # Get the variable name
  
  if (is_categorical(x)) {
    lvl <- levels(x)
    if (length(lvl) == 0) { # Handle factor with no levels
        # Use dynamic name for the empty tibble column
        out <- tibble::tibble(!!var_name := factor(character(0)))
        attr(out, "count") <- integer(0)
        return(out)
    }
    counts <- tabulate(match(x$value[, 1L], seq_along(lvl)), nbins = length(lvl))
    names(counts) <- lvl
    
    if (drop.empty) {
        keep <- counts > 0
        lvl_out <- lvl[keep]
        counts_out <- counts[keep]
    } else {
        lvl_out <- lvl
        counts_out <- counts
    }
    # Use dynamic name var_name for the tibble column
    out <- tibble::tibble(!!var_name := factor(lvl_out, levels = lvl)) 
    attr(out, "count") <- counts_out
    
  } else {
    # Continuous -> single pseudo-cell representing the variable
    # Use dynamic name var_name for the tibble column
    out <- tibble::tibble(!!var_name := var_name)
    count_val <- length(x$onsets)
    attr(out, "count") <- count_val
    names(attr(out, "count")) <- var_name # Name the count
     # Handle drop.empty for zero-event continuous case
    if (drop.empty && count_val == 0) {
          out <- out[0, , drop=FALSE]
          attr(out, "count") <- integer(0)
          names(attr(out, "count")) <- character(0)
      }
  }
  out
}

#' @method elements event
#' @rdname elements
#' @export
elements.event <- function(x, what = c("values", "labels"), transformed = TRUE, ...) {
  
  what <- match.arg(what)
  var_name_sanitized <- .sanitizeName(x$varname)
  
  element_data <- NULL
  
  if (what == "values") {
      # --- Handle VALUES --- 
      
      # Check cache first
      cache_val_attr <- paste0("_elements_cache_", what) # Unique cache attribute name
      cached_data <- attr(x, cache_val_attr)
      if (!is.null(cached_data)) {
          # Return cached data wrapped in a named list
          return(stats::setNames(list(cached_data), var_name_sanitized))
      }
      
      # Handle basis 'transformed' case (although !transformed is tricky)
      if (!is.null(x$meta$basis) && !transformed) {
          warning("'transformed = FALSE' for basis elements is not reliably supported, returning the transformed basis matrix.")
          element_data <- x$value
      } else {
           element_data <- x$value # This is always an N x K matrix (or 0xK)
      }
      
      # Store in cache before returning
      attr(x, cache_val_attr) <- element_data
      
  } else {
      # --- Handle LABELS --- 
      # Return descriptive names/levels for EACH event instance
      n_events <- length(x$onsets)
      lvls <- levels.event(x) # Get the names/levels (length K)
      
      if (n_events == 0) {
          # Handle empty event sequence consistently
          if (!x$continuous && !is.null(x$meta$levels)) {
              element_data <- factor(character(0), levels = lvls)
          } else if (length(lvls) <= 1) { # Single name/level
              element_data <- character(0)
          } else { # Multi-column name/level
              element_data <- matrix(character(0), nrow=0, ncol=length(lvls), dimnames=list(NULL, lvls))
          }
      } else if (!x$continuous && !is.null(x$meta$levels)) {
          # Categorical: Reconstruct factor values as a factor vector (length = n_events)
          element_data <- factor(x$value[, 1L], levels = seq_along(lvls), labels = lvls)
      } else if (length(lvls) == 1) {
          # Single continuous variable or basis: repeat name n_events times (length = n_events)
          element_data <- rep(lvls, n_events)
      } else {
          # Multi-column continuous (matrix/basis): create matrix repeating colnames (N x K)
          element_data <- matrix(rep(lvls, each = n_events), 
                                 nrow = n_events, 
                                 ncol = length(lvls), 
                                 dimnames = list(NULL, lvls))
      } 
  }
  
  # Return the data directly, not wrapped in a list
  element_data
}

#' Print event objects
#'
#' Provides a concise summary of an event object using cli.
#'
#' @param x An event object.
#' @param ... Additional arguments (unused).
#' @import cli
#' @export
#' @rdname print
print.event <- function(x, ...) {
  nevents <- length(x$onsets)
  type <- if (is_continuous(x)) "Continuous" else "Categorical"
  
  cli::cli_h1("Event Sequence: {.field {x$varname}}")
  
  cli::cli_div(theme = list(span.info = list(color = "blue")))
  cli::cli_text("{.info * Type:} {type}")
  
  # Display Levels (categorical) or Columns (continuous)
  lvls <- levels(x) # Uses levels.event S3 method
  if (!is_continuous(x)) {
    cli::cli_text("{.info * Levels:} {paste(lvls, collapse = ', ')}")
  } else {
    cli::cli_text("{.info * Columns:} {paste(lvls, collapse = ', ')}")
  }
  
  cli::cli_text("{.info * Events:} {nevents}")
  
  if (nevents > 0) {
    cli::cli_h2("Timing")
    onset_range <- range(x$onsets, na.rm = TRUE)
    dur_range <- range(x$durations, na.rm = TRUE)
    onset_range_str <- sprintf("%.2f - %.2f sec", onset_range[1], onset_range[2])
    dur_range_str <- sprintf("%.2f - %.2f sec", dur_range[1], dur_range[2])
    cli::cli_text("{.info * Onset Range:} {onset_range_str}")
    cli::cli_text("{.info • Duration Range:} {dur_range_str}")
    
    cli::cli_h2("Blocks")
    blocks_table <- table(x$blockids)
    nblocks <- length(blocks_table)
    cli::cli_text("{.info • Number of Blocks:} {nblocks}")
    max_show_blocks <- 10
    blocks_display <- if(nblocks > max_show_blocks) {
                          paste(c(names(blocks_table)[1:max_show_blocks], "..."), collapse = ", ")
                      } else {
                          paste(names(blocks_table), collapse = ", ")
                      }
    cli::cli_text("{.info • Block IDs:} {blocks_display}")
  } else {
      cli::cli_alert_info("Event sequence is empty.")
  }
  
  # Display Value Range for continuous (non-basis) variables
  if (is_continuous(x) && is.null(x$meta$basis)) {
      cli::cli_h2("Values")
      value_range <- range(x$value, na.rm = TRUE)
      value_range_str <- sprintf("%.2f - %.2f", value_range[1], value_range[2])
       cli::cli_text("{.info • Value Range:} {value_range_str}")
  }
  
  cli::cli_end()
  
  invisible(x)
}

## ============================================================================
## Section: Get Formatted Labels for Single Event
## ============================================================================

#' Get Formatted Labels for a Single Event
#'
#' Returns a character vector of formatted labels for an event object,
#' using the `Variable[Level]` style for categorical events, 
#' `Variable[Index]` for multi-column continuous events, or just
#' `Variable` for single continuous events.
#' Useful for getting consistent labels for individual event components.
#' This is distinct from `levels()` which returns the raw level names or column names.
#' Relies on the internal `.level_vector` helper function.
#'
#' @param x An object of class `event`.
#' @param ... Additional arguments (unused).
#' @return A character vector of formatted labels, or `character(0)` if not applicable.
#' @export
#' @examples
#' fac <- factor(rep(c("A", "B"), 3))
#' onsets <- 1:6
#' ev_fac <- event_factor(fac, "Condition", onsets)
#' labels(ev_fac) # Should return c("Condition[A]", "Condition[B]")
#' 
#' vals <- 1:6
#' ev_num <- event_variable(vals, "Modulator", onsets)
#' labels(ev_num) # Should return "Modulator"
#' 
#' mat <- matrix(1:12, 6, 2)
#' colnames(mat) <- c("C1", "C2")
#' ev_mat <- event_matrix(mat, "MatrixVar", onsets)
#' labels(ev_mat) # Should return c("MatrixVar[1]", "MatrixVar[2]") 
#'
labels.event <- function(x, ...) {
  # Directly call the internal helper (now in utils-internal.R)
  # Assumes utils-internal.R functions are available in package namespace
  .level_vector(x)
}

## ============================================================================
## Section: End of File
