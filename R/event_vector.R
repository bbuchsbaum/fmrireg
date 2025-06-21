###############################################################################
# EVENT_VECTOR.R
#
# This file contains helper routines for cleaning names, checking ordering,
# constructing event model terms (and various event types), extracting design
# matrices, computing contrasts, and printing event-related objects for fMRI.
#
###############################################################################

## ============================================================================
## Section 1: Helper Functions (Moved to R/utils-internal.R)
## ============================================================================

# Removed .sanitizeName
# Removed is.increasing
# Removed is.strictly.increasing
# Removed .checkEVArgs

## ============================================================================
## Section 2: Event Term Construction
## ============================================================================

#' Create an event model term from a named list of variables.
#'
#' Generates an `event_term` object which represents the combination of one 
#' or more event sequences (e.g., a factor crossed with a numeric modulator).
#' It takes a list of variables (factors, numeric vectors, matrices, basis objects)
#' along with shared onsets, block IDs, and durations.
#' It uses the `EV` factory internally to create standardized `event` objects for each variable.
#'
#' @param evlist A named list of variables (factors, numeric, matrices, ParametricBasis objects).
#'        The names are used as variable identifiers within the term.
#' @param onsets Numeric vector of onset times (in seconds).
#' @param blockids Numeric vector of block numbers (non-decreasing integers).
#' @param durations Numeric vector of event durations (seconds, default is 0). 
#'        Can be scalar (recycled) or vector matching length of `onsets`.
#' @param subset Optional logical vector indicating which events to retain (applied before processing).
#'
#' @return A list object with class `c("event_term", "event_seq")`. Contains:
#'   \item{varname}{Concatenated variable names from `evlist`.} 
#'   \item{events}{A named list of the processed `event` objects.} 
#'   \item{subset}{The `subset` vector used.} 
#'   \item{event_table}{A tibble representing the combinations of descriptive levels/names 
#'                    for each event in the term, constructed using `elements(..., values=FALSE)`.} 
#'   \item{onsets}{Numeric vector of onsets (after processing/subsetting).} 
#'   \item{blockids}{Numeric vector of block IDs (after processing/subsetting).} 
#'   \item{durations}{Numeric vector of durations (after processing/subsetting).} 
#'
#' @examples 
#' x1 <- factor(rep(letters[1:3], 10))
#' x2 <- factor(rep(1:3, each = 10))
#' onsets <- seq(1, 100, length.out = 30)
#' blockids <- rep(1:3, each = 10)
#' 
#' eterm <- event_term(list(Condition = x1, Group = x2),
#'                     onsets = onsets,
#'                     blockids = blockids)
#' print(eterm)
#' head(event_table(eterm))
#' levels(eterm)
#' head(design_matrix(eterm))
#'
#' @export
#' @import assertthat
#' @importFrom tibble as_tibble tibble
event_term <- function(evlist, onsets, blockids, durations = 0, subset = NULL) {
  
  # Convert blockids to numeric if they are factors.
  if (is.factor(blockids)) {
    blockids <- as.numeric(as.character(blockids))
  }
  
  # Ensure blockids are non-decreasing.
  # Use base R equivalent of is.increasing
  assert_that(!is.unsorted(blockids), msg = "'blockids' must consist of non-decreasing integers")
  
  # NOTE: subset handling is now primarily inside event(), 
  # but event_term might still need its own subset logic if it uses it before calling EV?
  # Currently, it passes subset down to EV.
  if (is.null(subset)) { 
    subset <- rep(TRUE, length(onsets)) 
  }
  
  # Duration recycling is handled inside event() via .checkEVArgs
  # if (length(durations) == 1) {
  #   durations <- rep(durations, length(onsets))
  # }
  
  vnames <- names(evlist)
  onlen <- length(onsets)
  
  # Basic check on input lengths before calling EV factory
  getlen <- function(v) {
    if (inherits(v, "event")) length(v$onsets)
    else if (is.matrix(v)) nrow(v)
    else if (inherits(v, "ParametricBasis")) nrow(v$y) # Check basis matrix dim
    else length(v)
  }
  for(i in seq_along(evlist)) {
      assert_that(getlen(evlist[[i]]) == onlen, 
                  msg=sprintf("Length mismatch between onsets (%d) and variable '%s' (%d)", 
                              onlen, vnames[i], getlen(evlist[[i]])))
  }
  
  # Create event objects by dispatching to the appropriate public wrapper
  evs <- lapply(seq_along(evlist), function(i) {
    vals <- evlist[[i]]
    vname_i <- vnames[i]
    
    # Type checking and dispatching (replaces EV factory logic)
    if (inherits(vals, "event")) {
        # If it's already an event object, just use it as-is
        # TODO: Consider whether to re-apply subset/durations/etc from current call
        vals
    } else if (inherits(vals, "ParametricBasis")) {
        event_basis(basis = vals, name = vname_i, onsets = onsets, blockids = blockids, durations = durations, subset = subset)
    } else if (is.factor(vals) || is.character(vals)) {
        event_factor(fac = vals, name = vname_i, onsets = onsets, blockids = blockids, durations = durations, subset = subset)
    } else if (is.matrix(vals)) {
        # Handle single-column matrices as vectors for event_variable?
        # No, event_matrix handles matrices directly now.
        event_matrix(mat = vals, name = vname_i, onsets = onsets, blockids = blockids, durations = durations, subset = subset)
    } else if (is.numeric(vals) && (is.vector(vals) || length(dim(vals)) <= 1)) {
         # Check specifically for numeric vectors (or 1D arrays)
        event_variable(vec = vals, name = vname_i, onsets = onsets, blockids = blockids, durations = durations, subset = subset)
    } else {
        stop(sprintf("Unsupported value type '%s' for variable '%s' in event_term", class(vals)[1], vname_i))
    }
  })
  
  names(evs) <- sapply(evs, function(ev) ev$varname)
  pterms <- unlist(lapply(evs, function(ev) ev$varname))
  
  # Check if any event object creation resulted in zero events after subsetting
  if (length(evs) > 0 && length(evs[[1]]$onsets) == 0) {
      warning(sprintf("Event term '%s' resulted in zero events after subsetting/processing.", 
                      paste(vnames, collapse=":")))
      # Return an empty structure? Or let downstream fail? Let downstream handle for now.
      # It should have empty onsets, durations, blockids, value matrix. 
  }
  
  # Rebuild event_table based on the *actual* content of the event objects
  # event() now handles the internal structure, so we use elements()
  # This relies on elements.event(..., values=FALSE) returning the appropriate factor levels/names
  
  # Get descriptive elements (levels/names) for each event
  descriptive_elements <- elements(evs[[1]], values = FALSE) # Need a way to call elements on list 'evs'? 
  # No, elements.event_term works on the term object after it's built.
  # Let's reconstruct the table *after* building the initial list. 

  varname <- paste(sapply(evs, function(x) x$varname), collapse = ":")
  
  # Create the list structure first
  ret <- list(varname = varname, 
              events = evs, 
              subset = subset, # Keep original subset for reference?
              # Placeholder for event_table, rebuild below
              event_table = NULL, 
              # Use onsets/blockids/durations from the *first* processed event object
              # Assumes they are consistent across all events after processing (checked by event())
              onsets = if(length(evs) > 0) evs[[1]]$onsets else numeric(0), 
              blockids = if(length(evs) > 0) evs[[1]]$blockids else numeric(0), 
              durations = if(length(evs) > 0) evs[[1]]$durations else numeric(0))
  
  class(ret) <- c("event_term", "event_seq")
  
  # Now build the event_table using elements()
  # Get descriptive elements (levels/names) for the term
  descriptive_elements_list <- elements(ret, what = "labels") # Explicitly request labels
  # Combine into a tibble directly from the list
  etab <- try(tibble::as_tibble(descriptive_elements_list), silent = TRUE)
  if (inherits(etab, "try-error")) {
      warning("Failed to create event_table for term: ", varname)
      etab <- tibble::tibble()
  }
  ret$event_table <- etab
  
  ret
}

#' @export
event_table.event_term <- function(x) x$event_table

## ============================================================================
## Section 3: EV Factory and Event Constructors (REMOVE EV Factory)
## ============================================================================

# Removed EV factory function. Logic is now inlined in event_term.

# Removed event_factor, event_variable, event_matrix, event_basis wrappers.
# They are now located in R/event-classes.R

## ============================================================================
## Section 4: Levels and Formula Methods
## ============================================================================

#' @method formula event_term
#' @export
formula.event_term <- function(x, ...) {
  # NOTE: This uses parent_terms.event_term which still exists below
  #       It might need adjustment if parent_terms logic changes.
  as.formula(paste("~ (", paste(parent_terms(x), collapse = ":"), "-1)"))
}

#' @noRd
.vector_of_labels <- function(ev) {
  # Helper to generate the vector of condition labels for a single event.
  # Uses levels.event() for the actual levels/column names.
  lvls <- levels(ev) # Get levels/colnames from levels.event
  
  if (is_continuous(ev) && length(lvls) > 1) {
    # Continuous multi-column (matrix/basis): Use index 1:ncol
    vapply(seq_along(lvls), 
           \(k) .label_component(ev, k), 
           character(1))
  } else if (is_categorical(ev)) {
    # Categorical: Use actual levels
    vapply(lvls, \(lvl) .label_component(ev, lvl), character(1))
  } else {
    # Single continuous variable (vector): Just the variable name
    .label_component(ev)
  }
}

## ============================================================================
## Section 5: Cells Extraction
## ============================================================================

#' @method cells event_term
#' @rdname cells
#' @export
cells.event_term <- function(x, drop.empty = TRUE, ...) {
  ## ----------------------------------------------------------------
  ## 0. fast cache ---------------------------------------------------
  # Use fixed name as levels rarely change post-construction
  cache_attr_name <- "..cells" 
  if (!is.null(cached <- attr(x, cache_attr_name))) {
    cnt <- attr(cached, "count")
    # Need to handle potential NULL count if cache is invalid
    if(is.null(cnt)) {
        warning("Invalid cache detected for cells.event_term, recomputing.")
    } else {
        return(if (drop.empty) cached[cnt > 0, , drop = FALSE] else cached)
    }
  }

  ## ----------------------------------------------------------------
  ## 1. categorical events only -------------------------------------
  # Use Filter and Negate for conciseness
  cats <- Filter(Negate(is_continuous), x$events)

  if (length(cats) == 0) {                  # no factors ⇒ one big cell
    # Use a more descriptive name if needed, maybe based on varname
    # Consistent with cells.event: use the (first/only) varname
    var_name_cont <- if (length(x$events) > 0) x$events[[1]]$varname else "all_events"
    out <- tibble::tibble(!!var_name_cont := var_name_cont) # Use varname for column
    count_val <- length(x$onsets)
    attr(out, "count") <- count_val
    # Assign name to the count attribute
    names(attr(out, "count")) <- var_name_cont 
    
    attr(x, cache_attr_name) <- out # Cache the result
    return(out)
  }

  ## ----------------------------------------------------------------
  ## 2. observed combination counts ---------------------------------
  # Reconstruct factors from internal representation
  obs_list <- lapply(cats, \(ev) {
      # Add checks for valid meta$levels and value structure
      if (is.null(ev$meta$levels) || !is.matrix(ev$value) || ncol(ev$value) != 1) {
           stop(paste("Invalid internal structure for categorical event:", ev$varname), call.=FALSE)
      }
      factor(ev$value[, 1L], levels = seq_along(ev$meta$levels), labels = ev$meta$levels)
  })
  # Create data frame of observed factor combinations
  obs <- do.call(data.frame, c(obs_list, stringsAsFactors = FALSE))
  
  # Use table() for efficient contingency counting (if few factors)
  # For many factors, alternative might be needed 
  # tbl <- as.data.frame.matrix(table(obs)) # This creates wide format, not needed directly
  
  # Generate the grid of all possible level combinations
  levels_list <- lapply(cats, levels)
  grid <- expand.grid(levels_list, stringsAsFactors = FALSE) # Use FALSE then convert relevant columns
  # Ensure column names match original variable names
  colnames(grid) <- names(cats)
  # Convert grid columns to factors matching the levels in obs_list
  for(i in seq_along(grid)){
      grid[[i]] <- factor(grid[[i]], levels=levels(obs_list[[i]]))
  }

  ## match() much faster than join for simple cases ------------------
  # Create unique string keys for observed and grid rows
  # Use sep that's unlikely to appear in levels
  key_obs  <- do.call(paste, c(obs, sep = "\001")) 
  key_grid <- do.call(paste, c(grid, sep = "\001"))
  # Count occurrences by matching observed keys to grid keys
  count    <- tabulate(match(key_obs, key_grid), nbins = nrow(grid))

  ## ----------------------------------------------------------------
  ## 3. assemble result ---------------------------------------------
  out <- tibble::as_tibble(grid) # Convert final grid to tibble

  # --- Add names to the count vector ---
  # Generate names for the counts based on the grid rows (factor level combinations)
  if (nrow(grid) > 0) {
      # Use a separator consistent with how interactions might be named elsewhere
      cell_names <- apply(grid, 1, paste, collapse = ":")
      # Ensure counts has names before attaching
      if (length(count) == length(cell_names)) {
         names(count) <- cell_names
      } else {
         # This shouldn't happen if logic above is correct, but add a warning
         warning("Mismatch between number of cells and counts calculated in cells.event_term")
      }
  } # Else count is likely integer(0) and doesn't need names
  # --- End adding names ---

  attr(out, "count") <- count # Attach the now named count vector

  attr(x, cache_attr_name) <- out # Cache the result (with named counts)

  # Filter based on drop.empty
  if (drop.empty) {
      keep_idx <- count > 0
      filtered_out <- out[keep_idx, , drop = FALSE]
      # Ensure the count attribute on the filtered output is also correct
      attr(filtered_out, "count") <- count[keep_idx]
      filtered_out
  } else {
       out
  }
}

#' @noRd
.event_set <- function(x, exclude_basis = FALSE) {
  evtab <- event_table(x)
  
  evset <- if (fmrihrf::nbasis(x) > 1 & !exclude_basis) {
    ncond <- fmrihrf::nbasis(x)
    # Construct a zero-padded string for basis labels.
    zstr <- paste0(rep("0", ceiling(log10(ncond + 1e-6))), collapse = "")
    
          evlist <- c(list(factor(paste("basis", zstr, 1:fmrihrf::nbasis(x), sep = ""))), cells(x$evterm))
    names(evlist) <- c("basis", parent_terms(x$evterm))
    evlist <- lapply(evlist, levels)
    ret <- expand.grid(evlist, stringsAsFactors = TRUE)
    ret[c(2:length(ret), 1)]
  } else {
    cells(x$evterm)
  }
}

#' @export
#' @rdname cells
cells.covariate_convolved_term <- function(x, ...) {
  unique(event_table(x))
}

#' @export
#' @importFrom stringr str_trim
cells.convolved_term <- function(x, exclude_basis = FALSE, ...) {
  evtab <- event_table(x)
  evset <- .event_set(x, exclude_basis = exclude_basis)
  
  strels <- apply(apply(evtab, 2, stringr::str_trim), 1, paste, collapse = ":")
  
  strlevs <- if (nrow(evset) > 1) {
    apply(apply(evset, 2, stringr::str_trim), 1, paste, collapse = ":")
  } else {
    as.character(evset[1, 1])
  }
  
  attr(evset, "rownames") <- strlevs
  
  counts <- if (exclude_basis) {
    rep(attr(cells(x$evterm), "count"), each = 1)
  } else {
          rep(attr(cells(x$evterm), "count"), each = fmrihrf::nbasis(x))
  }
  
  ret <- evset[counts > 0, , drop = FALSE]
  attr(ret, "count") <- counts[counts > 0]
  ret
}

## ============================================================================
## Section 6: Conditions and Parent Terms
## ============================================================================

#' @method conditions event_term
#' @rdname conditions
#' @export
#' @importFrom tibble tibble
conditions.event_term <- function(x, drop.empty = TRUE, expand_basis = FALSE, ...) {
  
  # --- Caching --- 
  opts_key <- paste(drop.empty, expand_basis, sep="|") # Keep drop.empty in key for now
  cached_val <- attr(x, "..conds")
  cached_opts <- attr(x, "..conds_opts")
  
  # --- RE-ENABLE CACHE --- 
  if (!is.null(cached_val) && !is.null(cached_opts) && identical(cached_opts, opts_key)) {
    return(cached_val)
  }
  # message("--- conditions.event_term: Cache bypassed/missed, recalculating ---") # Keep commented out
  # --- END RE-ENABLE ---
  
  # --- Shortcut for single continuous event with one column --- 
  if (length(x$events) == 1 && is_continuous(x$events[[1]])) {
    cols <- try(columns(x$events[[1]]), silent=TRUE)
    if (!inherits(cols, "try-error") && length(cols) == 1) {
        base_cond_tags <- cols 
        if (expand_basis) {
            hrfspec <- attr(x, "hrfspec")
            nb <- if (!is.null(hrfspec) && !is.null(hrfspec$hrf)) fmrihrf::nbasis(hrfspec$hrf) else 1L
            final_cond_tags <- add_basis(base_cond_tags, nb)
        } else {
            final_cond_tags <- base_cond_tags
        }
        attr(x, "..conds") <- final_cond_tags
        attr(x, "..conds_opts") <- opts_key
        return(final_cond_tags)
    }
  }
  
  # --- Generate Tokens for Each Component --- 
  comp_tokens_list <- lapply(x$events, function(ev) {
      if (is_categorical(ev)) {
          levs <- levels(ev) 
          # Handle case where factor might have no levels after subsetting? levels() should return character(0)
          if (length(levs) == 0) return(character(0))
          level_token(ev$varname, levs)
      } else {
          columns(ev) 
      }
  })
  
  # Filter out components that returned empty tokens (e.g., factors with no levels)
  comp_tokens_list <- Filter(function(tk) length(tk) > 0, comp_tokens_list)
  
  if (length(comp_tokens_list) == 0) { # If ALL components became empty
       final_out <- character(0)
       attr(x, "..conds") <- final_out
       attr(x, "..conds_opts") <- opts_key
       return(final_out)
  }
  
  # --- Combine Tokens using expand.grid and make_cond_tag --- 
  names(comp_tokens_list) <- names(Filter(function(tk) length(tk) > 0, x$events)) # Match names to filtered tokens
  full_grid <- expand.grid(comp_tokens_list, stringsAsFactors = FALSE)
  base_cond_tags_all <- apply(full_grid, 1, make_cond_tag)
  
  # --- REMOVED drop.empty LOGIC BLOCK --- 
  # The logic relying on cells() was flawed for mixed continuous/categorical terms.
  # Dropping based on actual matrix rank deficiency is handled by design_matrix() / model.matrix().
  base_cond_tags_final <- base_cond_tags_all
  
  # --- Handle expand_basis --- 
  if (expand_basis) {
      hrfspec <- attr(x, "hrfspec")
      nb <- if (!is.null(hrfspec) && !is.null(hrfspec$hrf)) fmrihrf::nbasis(hrfspec$hrf) else 1L
      final_cond_tags <- add_basis(base_cond_tags_final, nb)
  } else {
      final_cond_tags <- base_cond_tags_final
  }
  
  # --- Cache and Return --- 
  final_out <- as.vector(final_cond_tags)
  attr(x, "..conds") <- final_out
  attr(x, "..conds_opts") <- opts_key
  
  return(final_out)
}

#' @noRd
parent_terms.event_term <- function(x) unlist(lapply(x$events, function(ev) ev$varname))

## ============================================================================
## Section 7: Continuous/Categorical Checks and Elements Extraction
## ============================================================================

#' @export
is_continuous.event_term <- function(x) all(sapply(x$events, function(x) is_continuous(x)))

#' @export
is_categorical.event_term <- function(x) !is_continuous(x)

#' @method elements event_term
#' @rdname elements
#' @export
elements.event_term <- function(x, what = c("values", "labels"), ...) {
  # Ensure 'what' is determined correctly, respecting 'values' if passed via ...
  dots <- list(...)
  if (!missing(what)) {
      what <- match.arg(what)
  } else if (!is.null(dots$values)) {
      what <- if (dots$values) "values" else "labels"
  } else {
      what <- "values" # Default if neither 'what' nor 'values' is specified
  }
  
  # Pass the determined 'what' argument and any other arguments down
  # lapply will now create a named list directly as elements.event returns vectors/matrices
  els <- lapply(x$events, elements, what = what, ...)
  
  # Sanitize names of the resulting list (which should already be named by lapply)
  # If lapply didn't preserve names from x$events, this is needed.
  names(els) <- vapply(names(x$events), .sanitizeName, character(1))
  els
}

## ============================================================================
## Section 8: Onsets and Block IDs
## ============================================================================

#' @export
onsets.convolved_term <- function(x) {
  onsets(x$evterm)
}

#' @export
onsets.event_term <- function(x) {
  x$onsets
}

#' @export
blockids.event_term <- function(x) {
  x$blockids
}

#' @export
blockids.convolved_term <- function(x) {
  fmrihrf::blockids(x$evterm)
}

## ============================================================================
## Section 9: Splitting Onsets
## ============================================================================

#' @method split_onsets event_term
#' @rdname split_onsets
#' @export
split_onsets.event_term <- function(x, sframe, global = FALSE, blocksplit = FALSE, ...) {
  # Get categorical events.
  facs <- x$events[!sapply(x$events, is_continuous)]
  
  if (length(facs) == 0) {
    ons <- if (global) {
      fmrihrf::global_onsets(sframe, onsets(x), fmrihrf::blockids(x))
    } else {
      onsets(x)
    }
    return(list(split(ons, fmrihrf::blockids(x))))
  }
  
  # For categorical events, construct a crossed factor.
  facs <- lapply(facs, function(fac) unlist(elements(fac)))
  
  f <- function(...) {
    interaction(..., drop = TRUE, sep = ":")
  }
  
  cfac <- try(do.call(f, facs))
  # If error, a more informative error message might be warranted.
  
  ret <- if (global) {
    split(fmrihrf::global_onsets(sframe, onsets(x), fmrihrf::blockids(x)), cfac)
  } else {
    split(onsets(x), cfac)
  }
  
  if (blocksplit) {
    bsplit <- split(fmrihrf::blockids(x), cfac)
    ret <- lapply(seq_along(ret), function(i) {
      split(ret[[i]], bsplit[[i]])
    })
  }
  
  names(ret) <- longnames(x)
  ret
}

## ============================================================================
## Section 10: Convolution and Regressor Generation
## ============================================================================

#' Convolve HRF with Design Matrix.
#'
#' Convolves a HRF with a design matrix (one column per condition) to produce a
#' list of regressors.
#'
#' @param hrf A function representing the HRF.
#' @param dmat Design matrix (with named columns).
#' @param globons Numeric vector of global onsets.
#' @param durations Numeric vector of event durations.
#' @param summate Logical; if TRUE, summate the convolved HRF (default: TRUE).
#'
#' @return A list of regressors (one for each column).
#' @export
convolve_design <- function(hrf, dmat, globons, durations, summate = TRUE) {
  cond.names <- names(dmat)
  
  # Remove rows with NA values.
  if (any(is.na(dmat)) || any(is.na(globons))) {
    keep <- apply(dmat, 1, function(vals) all(!is.na(vals)))
    keep[is.na(globons)] <- FALSE
    dmat <- dmat[keep, ]
    durations <- durations[keep]
    globons <- globons[keep]
  }
  
  reglist <- purrr::map(1:ncol(dmat), function(i) {
    amp <- dmat[, i][[1]]
    nonzero <- which(amp != 0)
    if (length(nonzero) == 0) {
      # Call Reg directly to create an empty regressor object
      fmrihrf::regressor(onsets = numeric(0), hrf = hrf, amplitude = 0)
    } else {
      # Call Reg directly here as well
      fmrihrf::regressor(onsets = globons[nonzero], hrf = hrf, amplitude = amp[nonzero], duration = durations[nonzero], summate = summate)
    }
  })
  
  reglist
}

#' Extract regressors for an event term
#'
#' Convolve the event-term design matrix with an HRF and return the
#' resulting regressors.
#'
#' @rdname regressors
#' @param hrf HRF function
#' @param sampling_frame sampling_frame object
#' @param summate Logical; sum HRF responses
#' @param drop.empty Logical; drop empty conditions
#' @export
regressors.event_term <- function(x, hrf, sampling_frame, summate = FALSE, drop.empty = TRUE, ...) {
  globons <- fmrihrf::global_onsets(sampling_frame, x$onsets, x$blockids)
  durations <- x$durations
  blockids <- x$blockids
  nimages <- sum(fmrihrf::blocklens(sampling_frame))
  cnames <- conditions(x)
  dmat <- design_matrix(x, drop.empty)
  ncond <- ncol(dmat)
  
  reg <- convolve_design(hrf, dmat, globons, durations, summate = summate)
  names(reg) <- colnames(dmat)
  reg
}

#' @method convolve event_term
#' @rdname convolve
#' @export
#' @importFrom tibble as_tibble
#' @importFrom dplyr %>% group_by select do ungroup
convolve.event_term <- function(x, hrf, sampling_frame, drop.empty = TRUE, 
                                summate = TRUE, precision = 0.3, ...) {
  # Check for term_tag attribute (should have been added in realise_event_terms)
  term_tag <- attr(x, "term_tag")
  # --- REMOVED FALLBACK LOGIC FOR term_tag ---
  # If term_tag is NULL (e.g., for Ident()-only terms), make_column_names will handle it 
  # by not prepending a term_tag, resulting in direct variable names.
  # if (is.null(term_tag)) {
  #     warning(sprintf("Missing 'term_tag' attribute on event_term '%s'. Using $varname as fallback.", x$varname %||% "Unnamed"), call.=FALSE)
  #     term_tag <- x$varname %||% "UnnamedTerm"
  # }
  
  # --- Basic Setup --- 
  globons <- fmrihrf::global_onsets(sampling_frame, x$onsets, x$blockids)
  durations <- x$durations
  blockids <- x$blockids
  nimages <- sum(fmrihrf::blocklens(sampling_frame))
  
  # --- Get Unconvolved Design Matrix (dmat) --- 
  # design_matrix handles dropping empty/rank-deficient columns based on drop.empty
  dmat <- design_matrix(x, drop.empty = drop.empty)
  
  # --- Get Base Condition Names from the *actual* matrix columns --- 
  # This ensures names match the columns being convolved
  base_cnames <- colnames(dmat)
  
  # Check if dmat became empty after dropping
  if (ncol(dmat) == 0 || nrow(dmat) == 0) {
      warning(sprintf("Design matrix for term '%s' became empty after dropping. Convolution will result in an empty matrix.", term_tag), call.=FALSE)
      # Proceed to generate names for an empty matrix
      base_cnames <- character(0) # Use empty names
  }
  
  # --- Convolution per Block --- 
  block_ids <- unique(blockids)
  # Precompute global sample times once to avoid buggy blockids handling
  global_samples <- fmrihrf::samples(sampling_frame, global = TRUE)
  sample_blockids <- fmrihrf::blockids(sampling_frame)
  cmat_list <- lapply(block_ids, function(bid) {
    idx <- which(blockids == bid)
    # Ensure we subset the correct dmat based on drop.empty consistency
    dblock <- dmat[idx, , drop = FALSE] 
    globons_block <- globons[idx]
    durations_block <- durations[idx]
    
    # Skip block if dblock is empty (e.g., no events for this term in this block)
    if(nrow(dblock) == 0 || ncol(dblock) == 0) return(NULL) 
    
    reg <- convolve_design(hrf, dblock, globons_block, durations_block, summate = summate)
    
    # FIXED METHOD: Evaluate against all global samples but extract only block-specific rows
    # This preserves temporal accuracy while maintaining correct dimensions
    sam_all <- fmrihrf::samples(sampling_frame, global = TRUE)
    
    # Evaluate regressors against ALL samples (preserves temporal accuracy)
    full_block_mat <- do.call(cbind, lapply(reg, function(r) fmrihrf::evaluate(r, sam_all, precision = precision)))
    
    # Extract only the rows for this specific block (maintains correct dimensions)
    block_lengths <- fmrihrf::blocklens(sampling_frame)
    if (as.integer(bid) == 1) {
      block_start <- 1
      block_end <- block_lengths[1]
    } else {
      block_start <- sum(block_lengths[1:(as.integer(bid)-1)]) + 1
      block_end <- sum(block_lengths[1:as.integer(bid)])
    }
    block_mat <- full_block_mat[block_start:block_end, , drop = FALSE]
    block_mat
  })
  
  # Remove NULLs (from blocks with no events/cols) and rbind
  cmat_list <- Filter(Negate(is.null), cmat_list)
  
  # --- Generate Final Column Names --- 
  nb <- fmrihrf::nbasis(hrf)
  # Use the base_cnames derived directly from the dmat that was convolved
  cn <- make_column_names(term_tag, base_cnames, nb)
  
  # Handle case where convolution results in an empty matrix
  if (length(cmat_list) == 0) {
      warning(sprintf("Convolution resulted in an empty matrix for term '%s\'.\n  Returning tibble with correct names but 0 rows.", term_tag), call.=FALSE)
      # Return empty tibble with correct names and 0 rows
      return(tibble::as_tibble(matrix(numeric(0), nrow=0, ncol=length(cn)), 
                               .name_repair="minimal", .names_minimal = cn))
  }
  cmat <- do.call(rbind, cmat_list)
  
  # Handle add_sum flag if present (set by trialwise)
  if (isTRUE(attr(x, "add_sum"))) {
    if (ncol(cmat) > 0) { # Ensure there are columns to average
      mean_col <- matrix(rowMeans(cmat, na.rm = TRUE), ncol = 1)
      mean_col_name <- make.names(paste0(attr(x, "add_sum_label") %||% term_tag, "_mean"))
      colnames(mean_col) <- mean_col_name
      cmat <- cbind(cmat, mean_col)
      # Update column names to include the new mean column
      cn <- c(cn, mean_col_name)
    } else {
      warning(sprintf("Cannot add sum column for term '%s': no base columns generated.", term_tag))
    }
  }
  
  # Assign names, checking for length consistency
  if (length(cn) == ncol(cmat)) {
    colnames(cmat) <- cn
  } else {
      warning(sprintf("Final column name count (%d) mismatch with convolved matrix columns (%d) for term '%s'. Using generic names.", 
                      length(cn), ncol(cmat), term_tag), call. = FALSE)
      colnames(cmat) <- make.names(paste0("col_", seq_len(ncol(cmat))), unique=TRUE)
  }
  
  # --- Optional Debug Validation --- 
  if (getOption("fmrireg.debug", FALSE)) {
     if (exists("is_valid_heading", mode="function")){
        stopifnot(all(is_valid_heading(colnames(cmat))))
     } else {
        warning("fmrireg.debug=TRUE: is_valid_heading helper not found for validation.")
     }
  }
  
  # --- Return Result --- 
  suppressMessages(tibble::as_tibble(cmat, .name_repair = "minimal"))
}

## ============================================================================
## Section 11: F-Contrast Computation
## ============================================================================

#' @export
Fcontrasts.event_term <- function(x, max_inter = 4L, ...) {

  ## --- helpers -------------------------------------------------------------
  .is_cat <- function(ev) !is_continuous(ev)
  .Dmat   <- function(n) {
      if (n < 2) stop("Need at least 2 levels for contrasts.")
      con <- contr.sum(n)
      colnames(con) <- paste0("c", 1:(n - 1))
      con
  }
  .Cvec   <- function(n) matrix(1, nrow = n, ncol = 1)

  ## --- preparation ---------------------------------------------------------
  evs_cat <- Filter(.is_cat, x$events)
  if (!length(evs_cat)) stop("No categorical variables found in term '", x$varname, "' for Fcontrasts.", call.=FALSE)

  C <- lapply(evs_cat, function(ev) .Cvec(length(levels(ev))))
  D <- lapply(evs_cat, function(ev) .Dmat(length(levels(ev))))
  names(C) <- names(D) <- names(evs_cat)

  # --- Get expected row names in Kronecker order ---------------------------
  cat_levels_list <- lapply(evs_cat, levels)
  # --- build Cartesian product of levels (Kronecker order) ------------------ 
  lvl_grid        <- do.call(expand.grid, cat_levels_list)
  cat_cond_names  <- apply(lvl_grid, 1, paste, collapse = ":")
  expected_rows   <- nrow(lvl_grid)
  # ---------------------------------------------------------------------- 

  ## --- Compute main effects matrices (without rownames yet) ---------------
  main <- Map(function(i) {
      mat_list <- C
      mat_list[[i]] <- D[[i]] 
      Reduce(kronecker, mat_list)
  }, seq_along(D)) |> 
    stats::setNames(names(evs_cat))

  ## --- Compute interaction effects matrices (without rownames yet) -------
  final_contrasts_list <- if (length(D) > 1 && length(D) <= max_inter) {
      inter <- unlist(lapply(2:length(D), function(k) {
          combn(length(D), k, simplify = FALSE, FUN = function(ix) {
              mat_list <- C
              mat_list[ix] <- D[ix] 
              M <- Reduce(kronecker, mat_list)
              attr(M, "name") <- paste(names(evs_cat)[ix], collapse=":")
              M
          })
      }), recursive = FALSE)
      names(inter) <- vapply(inter, attr, "", "name")
      c(main, inter)
  } else {
      main
  }
  
  ## --- Assign correct categorical rownames to all matrices -------
  final_contrasts_named <- lapply(seq_along(final_contrasts_list), function(i) {
       M <- final_contrasts_list[[i]]
       mat_name <- names(final_contrasts_list)[i]
       if (!is.matrix(M)) { 
            warning(paste("Skipping rownames for invalid matrix in Fcontrasts list element:", mat_name))
            return(M)
       }
       # Compare nrow(M) to expected rows from CATEGORICAL grid/interaction
       if (nrow(M) == expected_rows) {
            # Use dimnames[[1]] <- assignment (should be correct now)
            dimnames(M)[[1]] <- cat_cond_names
       } else {
           # This warning should be less likely now, but keep for safety
           warning(paste("Dimension mismatch for contrast '", mat_name, "': expected ", 
                         expected_rows, " rows (from categorical interaction), but matrix has ", nrow(M), ". Rownames not assigned."))
       }
       M # Return matrix (modified in place)
  })
  
  names(final_contrasts_named) <- names(final_contrasts_list)
  
  final_contrasts_named
}

#' Retrieve contrast definitions for an event term
#'
#' This accessor returns the list of contrast specifications attached to the
#' term's originating hrfspec, if any.
#'
#' @param x An `event_term` object.
#' @param ... Unused.
#' @return A list of contrast specifications or `NULL` when none are defined.
#' @export
contrasts.event_term <- function(x, ...) {
  hrfspec <- attr(x, "hrfspec")
  if (is.null(hrfspec)) return(NULL)
  hrfspec$contrasts
}

## ============================================================================
## Section 12: Design Matrix Construction
## ============================================================================

#' @method design_matrix event_term
#' @rdname design_matrix
#' @export
#' @importFrom tibble as_tibble
design_matrix.event_term <- function(x, drop.empty = TRUE, ...) {

  # --- Special case: term contains only one continuous event --- 
  # This includes single numeric variables and multi-column basis functions.
  # Bypass model.matrix and return the value matrix directly.
  if (is_continuous(x) && length(x$events) == 1) {
    ev      <- x$events[[1]]
    # Directly use the value matrix (N x K)
    out_mat <- ev$value 
    # Ensure it's a data frame before naming
    out_df <- as.data.frame(out_mat) # Convert matrix to data frame
    
    # Use columns() to get the base condition tags for naming intermediate matrix
    # These tags represent the *parts* of the final name (e.g., "01", "02" for Poly)
    # Do NOT sanitize with make.names here, as these are intermediate component names.
    intermediate_cond_tags <- try(columns(ev), silent = TRUE)

    if (inherits(intermediate_cond_tags, "try-error") || length(intermediate_cond_tags) != ncol(out_df)) {
      warning(sprintf("Failed to get valid condition tags via columns() or column count mismatch for event term '%s' (varname: '%s'). Using generic V# names for intermediate matrix.", 
                      x$varname %||% "UnnamedTerm", ev$varname %||% "UnnamedEventVar"), call. = FALSE)
      # Fallback to generic, sanitized names if columns() fails or gives wrong number
      cnames <- make.names(paste0("V", seq_len(ncol(out_df))), unique = TRUE) 
    } else {
      # Use the raw condition tags (e.g., "01", "02") directly as intermediate names.
      # These are not necessarily valid full R names yet but are components.
      cnames <- intermediate_cond_tags 
    }
    
    names(out_df) <- cnames
    # Return as a tibble
    return(tibble::as_tibble(out_df, .name_repair = "minimal"))
  }
  
  ## ----------------------------------------------------------------
  ## 1. Build the "data" data-frame that model.matrix() needs
  ## ----------------------------------------------------------------
  #   * For categorical events → factor column (1 per event)
  #   * For continuous events  → numeric column(s) (one per matrix column)
  
  # Helper to extract appropriate column(s) from an event object
  build_cols <- function(ev) {
    if (is_categorical(ev)) {
      # Categorical: return data frame with factor column named ev$varname
      fac <- factor(ev$value[, 1L],
                    levels = seq_along(ev$meta$levels),
                    labels = ev$meta$levels)
      df_out <- data.frame(fac, check.names=FALSE)
      colnames(df_out) <- ev$varname 
      df_out
    } else {
      # Continuous/Basis: return data frame with a single matrix column named ev$varname
      mat_col <- ev$value # This is the N x K matrix
      df_out <- data.frame(I(mat_col)) # Use I() to store matrix in one column
      colnames(df_out) <- ev$varname # Name the column containing the matrix
      df_out
    }
  }
  
  # Combine columns from all events into a single data frame
  # cbind should now handle the mix of factor columns and matrix columns
  df_list <- lapply(x$events, build_cols)
  df <- try(do.call(cbind, df_list), silent=TRUE)
  if (inherits(df, "try-error")) {
      stop("Failed to construct data frame for model.matrix from event term: ", x$varname, 
           "\n  Original error: ", attr(df, "condition")$message)
  }
  
  # === DEBUG PRINT ===
  #message(sprintf("Term: %s, nrow(df): %d, length(x$onsets): %d", x$varname, nrow(df), length(x$onsets)))
  # === END DEBUG ===
  
  if (nrow(df) != length(x$onsets)) {
      stop("Internal error: Row mismatch when building data frame for event term: ", x$varname)
  }

  ## ----------------------------------------------------------------
  ## 2. model.matrix() with the term's formula
  ## ----------------------------------------------------------------
  # Get the formula (e.g., ~ Condition:Modulator - 1)
  form <- formula(x)
  
  # Special case: Check for single-level factors which cause model.matrix to fail
  # If any factor column has only one level, handle it specially
  has_single_level_factor <- FALSE
  for (col_name in colnames(df)) {
    col_data <- df[[col_name]]
    if (is.factor(col_data) && nlevels(col_data) <= 1) {
      has_single_level_factor <- TRUE
      break
    }
  }
  
  if (has_single_level_factor) {
    # For single-level factors, create a simple matrix of ones
    # This represents the constant effect of that single level
    n_rows <- nrow(df)
    mm <- matrix(1, nrow = n_rows, ncol = 1)
    # Use the condition name from conditions() for naming
    expected_colnames_raw <- conditions(x, drop.empty = FALSE)
    if (length(expected_colnames_raw) > 0) {
      colnames(mm) <- make.names(expected_colnames_raw[1], unique = TRUE)
    } else {
      colnames(mm) <- make.names(x$varname, unique = TRUE)
    }
  } else {
    # Normal case: use model.matrix
    # Ensure the data frame column names match what the formula expects
    # (formula uses original names, df uses sanitized/numbered names)
    # model.matrix should handle this via the data=df argument.
    
    mm <- try(model.matrix(form, data = df), silent=TRUE)
    if (inherits(mm, "try-error")) {
         stop("Failed to create model matrix for event term: ", x$varname, 
              "\n  Formula was: ", deparse(form),
              "\n  Data frame columns: ", paste(colnames(df), collapse=", "),
              "\n  Original error: ", attr(mm, "condition")$message)
    }
    
    # --- SET COLNAMES using conditions() as the single source of truth --- 
    # Get potentially *all* condition names first (before drop.empty)
    expected_colnames_raw <- conditions(x, drop.empty = FALSE)
    
    # Sanitize the raw names using make.names
    expected_colnames <- make.names(expected_colnames_raw, unique = TRUE)
    
    # Basic check: Does the number of columns match?
    # model.matrix might produce fewer columns if rank-deficient, 
    # but conditions() should produce the full set based on expand.grid.
    # This mismatch needs careful handling. For now, assume model.matrix is correct
    # in terms of *which* columns are estimable, and use conditions() to name them.
    # If ncol(mm) < length(expected_colnames), it implies model.matrix dropped some.
    
    # TODO: How to robustly map expected_colnames to the columns present in mm?
    # This is tricky. model.matrix column names (e.g., ConditionB:Modulator1)
    # don't directly map to conditions() output (e.g., Condition[B]:Modulator[1]).
    # For now, we ASSUME the order is the same and the number of columns matches 
    # *if the design is full rank*. If not, this naming will be wrong.
    # A more robust solution might involve parsing model.matrix colnames or attributes.
    # Let's proceed with the direct assignment as per the reviewer's suggestion, 
    # but acknowledge this fragility.
    if (ncol(mm) != length(expected_colnames)) {
        warning(sprintf("Column count mismatch for '%s': model.matrix (%d) vs conditions (%d). Naming may be incorrect due to rank deficiency.",
                        x$varname, ncol(mm), length(expected_colnames_raw)), call. = FALSE)
        # Attempt to name the existing columns anyway, hoping the order matches
        # This might fail if length(expected_colnames) is shorter, though unlikely.
        colnames(mm) <- expected_colnames[1:ncol(mm)] 
    } else {
        colnames(mm) <- expected_colnames
    }
  }
  
  ## ----------------------------------------------------------------
  ## 3. Drop empty columns if requested (optional)
  ## ----------------------------------------------------------------
  # model.matrix might return fewer columns than expected if interactions 
  # lead to rank deficiency. drop.empty applies to the *output* matrix.
  if (isTRUE(drop.empty)) {
      # Check for intercept columns and constant columns
      # Intercept columns are named "(Intercept)" 
      # Constant columns have zero variance but non-zero values (e.g., all ones)
      is_intercept <- (colnames(mm) == "(Intercept)")
      
      # Calculate variance and check for all-zero columns
      col_vars <- apply(mm, 2, var, na.rm = TRUE)
      col_all_zero <- colSums(abs(mm), na.rm = TRUE) == 0
      
      # A column should be kept if:
      # 1. It's an intercept column, OR
      # 2. It has non-zero variance (not constant), OR  
      # 3. It's a constant non-zero column (zero variance but not all zeros)
      is_constant_nonzero <- (col_vars < 1e-8 | is.na(col_vars)) & !col_all_zero
      keep_cols <- which(is_intercept | col_vars > 1e-8 | is.na(col_vars) | is_constant_nonzero)
      
      if (length(keep_cols) < ncol(mm)){
          # message("Dropping empty columns: ", paste(colnames(mm)[! (1:ncol(mm)) %in% keep_cols], collapse=", "))
          mm <- mm[, sort(keep_cols), drop = FALSE]
      }
  }
  
  # Optional: Handle NAs by zero-filling (historical behavior)
  # mm[is.na(mm)] <- 0 

  tibble::as_tibble(mm, .name_repair = "check_unique")
}

## ============================================================================
## Section 13: Print Methods
## ============================================================================

#' Print fmri_term objects.
#'
#' @param x An fmri_term object.
#' @param ... Additional arguments.
#' @export
#' @rdname print
print.fmri_term <- function(x, ...) {
  cat("fmri_term: ", class(x)[[1]], "\n")
  cat("  Term Name: ", x$varname, "\n")
  cat("  Num Rows: ", nrow(design_matrix(x)), "\n")
  cat("  Num Columns: ", ncol(design_matrix(x)), "\n")
}

#' Print convolved_term objects.
#'
#' @param x A convolved_term object.
#' @param ... Additional arguments.
#' @export
#' @rdname print
print.convolved_term <- function(x, ...) {
  cat("fmri_term: ", class(x)[[1]], "\n")
  cat("  Term Name: ", x$varname, "\n")
  cat("  Formula:  ", as.character(formula(x$evterm)), "\n")
  cat("  Num Events: ", nrow(x$evterm$event_table), "\n")
  cat("  Num Rows: ", nrow(design_matrix(x)), "\n")
  cat("  Num Columns: ", ncol(design_matrix(x)), "\n")
  cat("  Conditions: ", conditions(x), "\n")
  cat("  Term Types: ", paste(purrr::map_chr(x$evterm$events, ~ class(.)[[1]])), "\n")
}

#' Print afni_hrf_convolved_term objects.
#'
#' @param x An afni_hrf_convolved_term object.
#' @param ... Additional arguments.
#' @export
#' @rdname print
print.afni_hrf_convolved_term <- function(x, ...) {
  cat("fmri_term: ", class(x)[[1]], "\n")
  cat("  Term Name: ", x$varname, "\n")
  cat("  Formula:  ", as.character(formula(x$evterm)), "\n")
  cat("  Num Events: ", nrow(x$evterm$event_table), "\n")
  cat("  Conditions: ", conditions(x), "\n")
  cat("  Term Types: ", paste(purrr::map_chr(x$evterm$events, ~ class(.)[[1]])), "\n")
}

#' Print event_term objects
#'
#' Provides a concise summary of an event_term object using cli.
#'
#' @param x An event_term object.
#' @param ... Additional arguments (unused).
#' @import cli
#' @export
#' @rdname print
print.event_term <- function(x, ...) {
  nevents <- length(x$onsets)
  nvars <- length(x$events)
  
  cli::cli_h1("Event Term: {.field {x$varname}}")
  
  cli::cli_div(theme = list(span.info = list(color = "blue")))
  cli::cli_text("{.info • Number of Events:} {nevents}")
  cli::cli_text("{.info • Variables:} {paste(names(x$events), collapse = ", ")}")

  cli::cli_h2("Variable Types")
  if (nvars > 0) {
    for (name in names(x$events)) {
      # Use is_continuous generic method
      type <- if (is_continuous(x$events[[name]])) "Continuous" else "Categorical"
      # Use {.field {name}} for safer interpolation of the variable name
      cli::cli_text("{.info  • {.field {name}}:} {type}")
    }
  } else {
     cli::cli_text(" (No variables in term)")
  }
  
  if (nevents > 0) {
    cli::cli_h2("Timing")
    onset_range <- range(x$onsets, na.rm = TRUE)
    dur_range <- range(x$durations, na.rm = TRUE)
    # Evaluate sprintf outside cli::cli_text to avoid interpolation issues
    onset_range_str <- sprintf("%.2f - %.2f sec", onset_range[1], onset_range[2])
    dur_range_str <- sprintf("%.2f - %.2f sec", dur_range[1], dur_range[2])
    cli::cli_text("{.info • Onset Range:} {onset_range_str}")
    cli::cli_text("{.info • Duration Range:} {dur_range_str}")
    
    cli::cli_h2("Blocks")
    blocks_table <- table(x$blockids)
    nblocks <- length(blocks_table)
    cli::cli_text("{.info • Number of Blocks:} {nblocks}")
    # Truncate long block lists
    max_show_blocks <- 10
    blocks_display <- if(nblocks > max_show_blocks) {
                          paste(c(names(blocks_table)[1:max_show_blocks], "..."), collapse = ", ")
                      } else {
                          paste(names(blocks_table), collapse = ", ")
                      }
    cli::cli_text("{.info • Block IDs:} {blocks_display}")
    events_per_block_display <- if(nblocks > max_show_blocks) {
                                     paste(c(blocks_table[1:max_show_blocks], "..."), collapse = ", ")
                                 } else {
                                     paste(blocks_table, collapse = ", ")
                                 }
    cli::cli_text("{.info • Events per Block:} {events_per_block_display}")
  } else {
      cli::cli_alert_info("Event term is empty.")
  }
  cli::cli_end()
  
  invisible(x)
}

## ============================================================================
## Section 14: End of File
###############################################################################
# End of event_vector.R
