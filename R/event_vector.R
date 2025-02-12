###############################################################################
# EVENT_VECTOR.R
#
# This file contains helper routines for cleaning names, checking ordering,
# constructing event model terms (and various event types), extracting design
# matrices, computing contrasts, and printing event‐related objects for fMRI.
#
###############################################################################

## ============================================================================
## Section 1: Helper Functions
## ============================================================================

#' Sanitize a variable name.
#'
#' This function replaces problematic characters (such as colons, spaces,
#' parentheses, and commas) to produce a valid name. For potential performance
#' improvements, caching of sanitized names could be implemented.
#'
#' @param name A character string.
#' @return A sanitized version of the input name.
#' @keywords internal
.sanitizeName <- function(name) {
  name <- gsub(":", ".", name)
  name <- gsub(" ", "", name)
  name <- gsub("\\)$", "", name)
  name <- gsub("[\\(\\)]", "_", name, perl = TRUE)
  name <- gsub(",", "_", name)
  name <- gsub("\\.$", "", name)
  name
}

#' Check if a vector is non-decreasing.
#'
#' @param vec A numeric vector.
#' @return TRUE if the differences are all non-negative.
#' @keywords internal
is.increasing <- function(vec) {
  all(diff(vec) >= 0)
}

#' Check if a vector is strictly increasing.
#'
#' @param vec A numeric vector.
#' @return TRUE if all differences are strictly positive.
#' @keywords internal
is.strictly.increasing <- function(vec) {
  all(diff(vec) > 0)
}

#' Validate event arguments.
#'
#' Ensures that the event onsets, values, block IDs, and durations are
#' consistent. In particular, it asserts that onsets have no NA values and
#' are strictly increasing within each block.
#'
#' @param name The name of the event.
#' @param vals The event values.
#' @param onsets Numeric vector of event onsets.
#' @param blockids Numeric vector of block IDs.
#' @param durations Numeric vector of durations (or a scalar).
#' @return A list of validated event parameters.
#' @keywords internal
.checkEVArgs <- function(name, vals, onsets, blockids, durations = NULL) {
  assert_that(length(onsets) == length(vals))
  # Ensure onsets are not NA.
  assert_that(all(!is.na(onsets)))
  
  # Split onsets by block and verify that they are strictly increasing.
  sons <- split(onsets, blockids)
  for (ons in sons) {
    assert_that(is.strictly.increasing(ons))
  }
  
  assert_that(is.increasing(blockids))
  
  
  # Replicate durations if necessary.
  if (is.null(durations) || length(durations) == 1) {
    durations <- rep(durations, length(onsets))
  }
  
  assert_that(length(durations) == length(vals))
  assert_that(length(blockids) == length(vals))
  
  list(varname = name, value = vals, onsets = onsets, durations = durations, blockids = blockids)
}

## ============================================================================
## Section 2: Event Term Construction
## ============================================================================

#' Create an event model term from a named list of variables.
#'
#' Generates an event model term from a list of named variables along with their
#' onsets, block IDs, and durations. Optionally, a subset of onsets can be retained.
#'
#' @param evlist A list of named variables.
#' @param onsets Numeric vector of onset times (in seconds).
#' @param blockids Numeric vector of block numbers.
#' @param durations Numeric vector of event durations (default is 1).
#' @param subset Logical vector indicating which onsets to retain (default is NULL).
#'
#' @return A list with components:
#'   - varname: concatenated variable names,
#'   - events: a list of event objects,
#'   - subset: the retained onsets,
#'   - event_table: a tibble with event information,
#'   - onsets, blockids, durations.
#'
#' @examples 
#' x1 <- factor(rep(letters[1:3], 10))
#' x2 <- factor(rep(1:3, each = 10))
#' eterm <- event_term(list(x1 = x1, x2 = x2),
#'                     onsets = seq(1, 100, length.out = 30),
#'                     blockids = rep(1, 30))
#'
#' @export
event_term <- function(evlist, onsets, blockids, durations = 1, subset = NULL) {
  
  # Convert blockids to numeric if they are factors.
  if (is.factor(blockids)) {
    blockids <- as.numeric(as.character(blockids))
  }
  
  # Ensure blockids are non-decreasing.
  assert_that(is.increasing(blockids), msg = "'blockids' must consist of strictly increasing integers")
  
  if (is.null(subset)) { 
    subset <- rep(TRUE, length(onsets)) 
  }
  
  if (length(durations) == 1) {
    durations <- rep(durations, length(onsets))
  }
  
  vnames <- names(evlist)
  onlen <- length(onsets)
  # Create event objects using the EV factory function.
  
  getlen <- function(v) {
    if (is.matrix(v)) {
      nrow(v)
    } else {
      length(v)
    }
  }
  
  
  
  evs <- lapply(seq_along(evlist), function(i) {
    #assertthat::assert_that(getlen(evlist[[i]]) == onlen, msg="all event variables must have the same length")
    EV(evlist[[i]], vnames[i], onsets = onsets, blockids = blockids, durations = durations, subset = subset)
  })
  
  names(evs) <- sapply(evs, function(ev) ev$varname)
  pterms <- unlist(lapply(evs, function(ev) ev$varname))
  
  len <- sum(subset)
  
  # Build event_table: For continuous events, repeat the sanitized name; otherwise, use the event values.
  etab <- suppressMessages(tibble::as_tibble(lapply(pterms, function(termname) {
    if (is_continuous(evs[[termname]])) {
      rep(.sanitizeName(termname), len)
    } else {
      evs[[termname]]$value
    }			
  }), .name_repair = "check_unique"))
  
  names(etab) <- sapply(pterms, .sanitizeName)
  varname <- paste(sapply(evs, function(x) x$varname), collapse = ":")
  
  ret <- list(varname = varname, 
              events = evs, 
              subset = subset, 
              event_table = etab, 
              onsets = evs[[1]]$onsets, 
              blockids = evs[[1]]$blockids, 
              durations = evs[[1]]$durations)
  class(ret) <- c("event_term", "event_seq")
  ret
}

#' @export
event_table.event_term <- function(x) x$event_table

## ============================================================================
## Section 3: EV Factory and Event Constructors
## ============================================================================

#' EV
#'
#' Factory function for creating event objects (e.g. event_factor, event_variable,
#' event_basis, event_matrix) based on the type of input values.
#'
#' @param vals Event values.
#' @param name Name of the event variable.
#' @param onsets Numeric vector of event onsets.
#' @param blockids Numeric vector of block IDs.
#' @param durations Numeric vector of event durations (or scalar).
#' @param subset Logical vector indicating which events to keep.
#'
#' @return An event object.
#'
#' @examples
#' ev_fac <- EV(factor(c("A", "B", "C")), "fac", onsets = c(1, 10, 20), blockids = rep(1, 3))
#'
#' @keywords internal
#' @export
EV <- function(vals, name, onsets, blockids, durations = 1, subset = rep(TRUE, length(onsets))) {
  if (length(durations) == 1) {
    durations <- rep(durations, length(onsets))
  }
  
  if (is.matrix(vals) && NCOL(vals) == 1) {
    vals <- vals[, 1, drop = TRUE]
  }
  
  ## Eagerly apply subset.
  if (inherits(vals, "ParametricBasis")) {
    event_basis(vals, onsets, blockids, durations, subset)	
  } else if (is.factor(vals) || is.character(vals)) {
    vals <- factor(as.character(vals)[subset])
    event_factor(vals, name, onsets[subset], blockids[subset], durations[subset])
  } else if (is.matrix(vals)) {
    event_matrix(vals[subset, ], name, onsets[subset], blockids[subset], durations[subset])
  } else if (is.numeric(vals) & is.vector(vals)) {
    event_variable(vals[subset], name, onsets[subset], blockids[subset], durations[subset])
  } else {
    stop(paste("cannot create event_seq from type:", typeof(vals)))
  }
}

#' Create a categorical event sequence from a factor.
#'
#' @param fac A factor representing the categorical event sequence.
#' @param name Name of the event sequence.
#' @param onsets Numeric vector of onsets.
#' @param blockids Numeric vector of block IDs (default: rep(1, length(fac))).
#' @param durations Numeric vector of durations (default: rep(0, length(fac))).
#'
#' @return An object of class "event_factor" and "event_seq".
#'
#' @examples
#' efac <- event_factor(factor(c("a", "b", "c", "a", "b", "c")), "abc", onsets = seq(1, 100, length.out = 6))
#'
#' @seealso \code{\link{event_model}}
#' @export 
event_factor <- function(fac, name, onsets, blockids = rep(1, length(fac)), durations = rep(0, length(fac))) {
  if (!is.factor(fac)) {
    warning("argument 'fac' is not a factor, converting to factor")
    fac <- factor(as.character(fac))
  }
  
  ret <- .checkEVArgs(name, fac, onsets, blockids, durations)
  ret$continuous <- FALSE
  class(ret) <- c("event_factor", "event_seq")
  ret
}        

#' Create a continuous event sequence from a numeric vector.
#'
#' @param vec Numeric vector representing continuous event values.
#' @param name Name of the event sequence.
#' @param onsets Numeric vector of onsets.
#' @param blockids Numeric vector of block IDs (default: 1).
#' @param durations Numeric vector of durations (default: 0).
#'
#' @return An object of class "event_variable" and "event_seq".
#'
#' @examples
#' evar <- event_variable(c(1, 2, 3, 4, 5, 6), "example_var", onsets = seq(1, 100, length.out = 6))
#'
#' @seealso \code{\link{event_factor}}
#' @export
event_variable <- function(vec, name, onsets, blockids = 1, durations = 0) {
  stopifnot(is.vector(vec))
  
  if (is.factor(vec)) {
    stop("cannot create an event_variable from a factor, use 'event_factor'.")
  }
  
  if (length(durations) == 1) {
    durations <- rep(durations, length(onsets))
  }
  if (length(blockids) == 1) {
    blockids <- rep(blockids, length(onsets))
  }
  
  ret <- .checkEVArgs(name, vec, onsets, blockids, durations)
  ret$continuous <- TRUE
  class(ret) <- c("event_variable", "event_seq")
  ret
}       

#' Create a continuous event set from a matrix.
#'
#' @param mat A matrix of continuous event values (one row per event).
#' @param name Name for the event set.
#' @param onsets Numeric vector of onsets.
#' @param blockids Numeric vector of block IDs (default: rep(1, ncol(mat))).
#' @param durations Numeric vector of durations (default: NULL).
#'
#' @return An object of class "event_matrix" and "event_seq".
#'
#' @examples
#' mat <- matrix(rnorm(200), 100, 2)
#' onsets <- seq(1, 1000, length.out = 100)
#' durations <- rep(1, 100)
#' blockids <- rep(1, 100)
#' eset <- event_matrix(mat, "eset", onsets, durations, blockids)
#'
#' @export
event_matrix <- function(mat, name, onsets, blockids = rep(1, ncol(mat)), durations = NULL) {
  stopifnot(is.matrix(mat))
  if (is.null(durations)) {
    durations <- rep(0, nrow(mat))
  }
  ret <- .checkEVArgs(name, as.vector(mat[, 1]), onsets, blockids, durations)
  ret$continuous <- TRUE
  
  if (is.null(colnames(mat))) {
    colnames(mat) <- 1:NCOL(mat)
  }
  
  ret$value <- mat
  class(ret) <- c("event_matrix", "event_seq")
  ret
}

#' Create an event set from a ParametricBasis object.
#'
#' @param basis A ParametricBasis object.
#' @param onsets Numeric vector of onsets.
#' @param blockids Numeric vector of block IDs (default: 1).
#' @param durations Numeric vector of durations (default: 0).
#' @param subset Logical vector for subsetting (default: rep(TRUE, length(onsets))).
#'
#' @return An object of class "event_basis" and "event_seq".
#'
#' @import assertthat
#' @examples
#' basis <- BSpline(1:21, 3)
#' onsets <- seq(0, 20, length.out = 21)
#' blockids <- rep(1, length(onsets))
#' ebasis <- event_basis(basis, onsets, blockids)
#'
#' @export
event_basis <- function(basis, onsets, blockids = 1, durations = 0, subset = rep(TRUE, length(onsets))) {
  assertthat::assert_that(inherits(basis, "ParametricBasis"))
  
  if (any(!subset)) {
    basis <- sub_basis(basis, subset)
  }
  
  if (length(durations) == 1) {
    durations <- rep(durations, length(onsets))
  }
  
  ret <- .checkEVArgs(basis$name, basis$y[, 1], onsets[subset], blockids[subset], durations[subset])
  ret$value <- basis$y
  ret$continuous <- TRUE
  ret$basis <- basis
  class(ret) <- c("event_basis", "event_seq")
  ret
}

## ============================================================================
## Section 4: Levels and Formula Methods
## ============================================================================

#' @export
levels.event_factor <- function(x) levels(x$value)

#' @export
levels.event_variable <- function(x) x$varname

#' @export
levels.event_matrix <- function(x) colnames(x$value)

#' @export
levels.event_set <- function(x) colnames(x$value)

#' @export
levels.event_basis <- function(x) levels(x$basis)

#' Build a formula for an event term.
#'
#' @export
formula.event_term <- function(x, ...) {
  as.formula(paste("~ (", paste(parent_terms(x), collapse = ":"), "-1)"))
}

#' @export
levels.event_term <- function(x) {
  facs <- x$events[!sapply(x$events, is_continuous)]
  if (length(facs) == 1) {
    levels(facs[[1]])
  } else {
    facs <- lapply(facs, function(x) x$value)
    f <- function(...) interaction(..., drop = TRUE, sep = ":")
    levels(do.call(f, facs))
  }
}

## ============================================================================
## Section 5: Cells Extraction
## ============================================================================

#' Retrieve cells of an event_factor object.
#'
#' @param x An event_factor object.
#' @param drop.empty Logical; if TRUE, empty cells are removed (not implemented).
#' @param ... Additional arguments.
#'
#' @return A list of data frames containing the cells.
#'
#' @export
#' @rdname cells
cells.event_factor <- function(x, drop.empty = TRUE, ...) {
  etab <- data.frame(onsets = x$onsets, durations = x$durations, blockids = x$blockids)
  split(etab, x$value)
}

#' Retrieve cells of an event_term object.
#'
#' @param x An event_term object.
#' @param drop.empty Logical; if TRUE, remove cells with no events.
#' @param ... Additional arguments.
#'
#' @return A tibble containing the cells with a "count" attribute.
#'
#' @export
#' @rdname cells
#' @examples
#' evlist <- list(fac1 = factor(c("A", "B", "A", "B")), 
#'                fac2 = factor(c("1", "1", "2", "2")))
#' eterm <- event_term(evlist, onsets = 1:4, blockids = rep(1, 4))
#' cells(eterm)
cells.event_term <- function(x, drop.empty = TRUE, ...) {
  evtab <- x$event_table
  evset <- suppressMessages(tibble::as_tibble(expand.grid(lapply(x$events, levels)), .name_repair = "check_unique"))
  
  which.cat <- which(!sapply(x$events, is_continuous))
  
  if (length(which.cat) > 0) {
    evs <- suppressMessages(tibble::as_tibble(lapply(evset[, which.cat], as.character), .name_repair = "check_unique"))
    evt <- suppressMessages(tibble::as_tibble(lapply(evtab[, which.cat], as.character), .name_repair = "check_unique"))
    
    counts <- apply(evs, 1, function(row1) {
      sum(apply(evt, 1, function(row2) { all(as.character(row1) == as.character(row2)) }))
    })
    
    if (drop.empty) {
      evset <- evset[counts > 0, , drop = FALSE]
      attr(evset, "count") <- counts[counts > 0]
    } else {
      attr(evset, "count") <- counts
    }
  } else {
    attr(evset, "count") <- nrow(evtab)
  }
  
  evset
}

#' @noRd
.event_set <- function(x, exclude_basis = FALSE) {
  evtab <- event_table(x)
  
  evset <- if (nbasis(x) > 1 & !exclude_basis) {
    ncond <- nbasis(x)
    # Construct a zero-padded string for basis labels.
    zstr <- paste0(rep("0", ceiling(log10(ncond + 1e-6))), collapse = "")
    
    evlist <- c(list(factor(paste("basis", zstr, 1:nbasis(x), sep = ""))), cells(x$evterm))
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
    rep(attr(cells(x$evterm), "count"), each = nbasis(x))
  }
  
  ret <- evset[counts > 0, , drop = FALSE]
  attr(ret, "count") <- counts[counts > 0]
  ret
}

## ============================================================================
## Section 6: Conditions and Parent Terms
## ============================================================================

#' @export
conditions.fmri_term <- function(x, ...) {
  colnames(design_matrix(x))
}

#' @export
#' @family conditions
conditions.convolved_term <- function(x, ...) {
  colnames(design_matrix(x))
}

#' @export
#' @family conditions
conditions.afni_hrf_convolved_term <- function(x, ...) {
  conditions(x$evterm)
}

#' @export
#' @family conditions
conditions.afni_trialwise_convolved_term <- function(x, ...) {
  conditions(x$evterm)
}

#' Extract conditions from an event_term.
#'
#' Constructs a vector of condition labels by combining the levels of
#' the categorical events.
#'
#' @param x An event_term object.
#' @param drop.empty Logical; if TRUE, conditions with no cells are dropped (default is TRUE).
#' @param ... Additional arguments.
#' @return A character vector of condition labels.
#' @export
conditions.event_term <- function(x, drop.empty = TRUE, ...) {
  .cells <- cells(x, drop.empty = drop.empty)
  pterms <- parent_terms(x)
  levs <- apply(.cells, 1, paste, collapse = ":")
  
  splitlevs <- strsplit(levs, ":")
  ret <- lapply(seq_along(pterms), function(i) {
    lev <- sapply(splitlevs, "[[", i)
    term <- pterms[[i]]
    if (length(levels(x$events[[i]])) > 1) {
      paste(.sanitizeName(pterms[i]), "[", lev, "]", sep = "")
    } else {
      .sanitizeName(pterms[i])
    }
  })
  
  do.call(function(...) paste(..., sep = ":"), ret)
}

#' @noRd
columns.event_term <- function(x) as.vector(unlist(lapply(x$events, columns)))

#' @noRd
columns.event_seq <- function(x) x$varname

#' @noRd
columns.event_matrix <- function(x) paste0(.sanitizeName(x$varname), ".", levels(x))

#' @noRd
columns.event_set <- function(x) paste0(.sanitizeName(x$varname), ".", levels(x))

#' @noRd
columns.event_basis <- function(x) columns(x$basis)

#' @noRd
parent_terms.event_term <- function(x) unlist(lapply(x$events, function(ev) ev$varname))

## ============================================================================
## Section 7: Continuous/Categorical Checks and Elements Extraction
## ============================================================================

#' @export
is_continuous.event_seq <- function(x) x$continuous

#' @export
is_continuous.event_factor <- function(x) FALSE

#' @export
is_categorical.event_seq <- function(x) TRUE

#' @export
is_continuous.event_term <- function(x) all(sapply(x$events, function(x) is_continuous(x)))

#' @export
is_categorical.event_term <- function(x) !is_continuous(x)

#' @export
is_categorical.event_seq <- function(x) !x$continuous

#' @export
elements.event_matrix <- function(x, values = TRUE, ...) {
  if (values) {
    ret <- x$value
    colnames(ret) <- colnames(x)
    ret <- list(ret)
    names(ret) <- .sanitizeName(x$varname)
    ret
  } else {
    N <- length(x$onsets)
    vnames <- colnames(x)
    res <- lapply(vnames, function(el) rep(el, N))
    mat <- do.call(cbind, res)
    colnames(mat) <- vnames			
    ret <- list(mat)
    names(ret) <- .sanitizeName(x$varname)
    ret			
  }
}

#' @export
elements.event_seq <- function(x, values = TRUE, ...) {
  if (values) {
    ret <- list(x$value)
    names(ret) <- x$varname
    ret
  } else {
    ret <- list(rep(x$varname, length(x)))
    names(ret) <- x$varname
    ret
  }
}

#' @export
elements.event_basis <- function(x, values = TRUE, transformed = TRUE, ...) {
  if (values && !transformed) {
    x$value$x				
  } else if (values) {
    ret <- x$basis$y
    colnames(ret) <- columns(x)
    n <- .sanitizeName(x$varname)
    ret <- list(ret)
    names(ret) <- n
    ret
  } else {
    N <- length(x)
    vnames <- columns(x)
    res <- lapply(vnames, function(el) rep(el, N))
    mat <- do.call(cbind, res)
    colnames(mat) <- vnames			
    ret <- list(mat)
    names(ret) <- .sanitizeName(x$varname)
    ret		
  }
}

#' @export
elements.event_term <- function(x, values = TRUE, ...) {
  els <- lapply(x$events, elements, values = values)
  n <- sapply(names(els), function(nam) .sanitizeName(nam))
  names(els) <- as.vector(n)
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
  blockids(x$evterm)
}

## ============================================================================
## Section 9: Splitting Onsets
## ============================================================================

#' Split onsets of an event_term object.
#'
#' Splits the onsets based on factor levels or block IDs.
#'
#' @param x An event_term object.
#' @param sframe A data frame representing the sampling frame.
#' @param global Logical; if TRUE, use global onsets (default: FALSE).
#' @param blocksplit Logical; if TRUE, split onsets by block IDs (default: FALSE).
#' @param ... Additional arguments.
#'
#' @return A list of numeric vectors for each factor level or block.
#'
#' @export
split_onsets.event_term <- function(x, sframe, global = FALSE, blocksplit = FALSE, ...) {
  # Get categorical events.
  facs <- x$events[!sapply(x$events, is_continuous)]
  
  if (length(facs) == 0) {
    ons <- if (global) {
      global_onsets(sframe, onsets(x), blockids(x))
    } else {
      onsets(x)
    }
    return(list(split(ons, blockids(x))))
  }
  
  # For categorical events, construct a crossed factor.
  facs <- lapply(facs, function(fac) unlist(elements(fac)))
  
  f <- function(...) {
    interaction(..., drop = TRUE, sep = ":")
  }
  
  cfac <- try(do.call(f, facs))
  # If error, a more informative error message might be warranted.
  
  ret <- if (global) {
    split(global_onsets(sframe, onsets(x), blockids(x)), cfac)
  } else {
    split(onsets(x), cfac)
  }
  
  if (blocksplit) {
    bsplit <- split(blockids(x), cfac)
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
      null_regressor(hrf)
    } else {
      regressor(globons[nonzero], hrf, amplitude = amp[nonzero], duration = durations[nonzero], summate = summate)
    }
  })
  
  reglist
}

#' @export
regressors.event_term <- function(x, hrf, sampling_frame, summate = FALSE, drop.empty = TRUE) {
  globons <- global_onsets(sampling_frame, x$onsets, x$blockids)
  durations <- x$durations
  blockids <- x$blockids
  nimages <- sum(sampling_frame$blocklens)
  cnames <- conditions(x)
  dmat <- design_matrix(x, drop.empty)
  ncond <- ncol(dmat)
  
  reg <- convolve_design(hrf, dmat, globons, durations, summate = summate)
  names(reg) <- colnames(dmat)
  reg
}

#' Convolve an event-related design matrix with an HRF.
#'
#' This function takes an event-related design matrix and convolves it with a given
#' HRF to produce a new design matrix suitable for fMRI analysis.
#'
#' @importFrom tibble as_tibble
#' @importFrom dplyr %>% group_by select do ungroup
#' @autoglobal
#' @export
#' @param x A data frame with the design matrix.
#' @param hrf The HRF function.
#' @param sampling_frame Data frame specifying the sampling frame.
#' @param drop.empty Logical; if TRUE, remove empty rows.
#' @param summate Logical; if TRUE, sum the convolved HRF.
#' @param precision Numeric; the convolution precision (default: 0.3).
#' @param ... Additional arguments.
#'
#' @return A tibble of the convolved design matrix.
convolve.event_term <- function(x, hrf, sampling_frame, drop.empty = TRUE, 
                                summate = TRUE, precision = 0.3, ...) {
  globons <- global_onsets(sampling_frame, x$onsets, x$blockids)
  durations <- x$durations
  blockids <- x$blockids
  nimages <- sum(sampling_frame$blocklens)
  cnames <- conditions(x)
  dmat <- design_matrix(x, drop.empty)
  ncond <- ncol(dmat)
  
  # Group by block IDs and process each group separately.
  cmat <- dmat |>
    dplyr::mutate(.blockids = blockids, .globons = globons, .durations = durations) |>
    dplyr::group_by(.blockids) |>
    dplyr::do({
      d <- dplyr::select(., 1:ncond)
      reg <- convolve_design(hrf, d, .$.globons, .$.durations, summate = summate)
      sam <- samples(sampling_frame, blockids = as.integer(as.character(.$.blockids[1])), global = TRUE)
      ## NOTE: This loop could be parallelized if necessary.
      ret <- do.call(cbind, lapply(seq_along(reg), function(ri) {
        vname <- paste0("v", ri)
        evaluate(reg[[ri]], sam, precision = precision)
      }))
      ret <- suppressMessages(tibble::as_tibble(ret, .name_repair = "minimal"))
      names(ret) <- paste0("v", 1:length(reg))
      ret
    })
  
  cmat <- cmat %>% dplyr::ungroup() %>% dplyr::select(-.blockids)
  
  if (nbasis(hrf) > 1) {
    blevs <- paste("[", 1:nbasis(hrf), "]", sep = "")
    cnames <- unlist(lapply(cnames, function(prefix) paste(prefix, ":basis", blevs, sep = "")))
  }
  
  colnames(cmat) <- cnames
  suppressMessages(tibble::as_tibble(cmat, .name_repair = "check_unique"))
}

## ============================================================================
## Section 11: F-Contrast Computation
## ============================================================================

#' Compute F-contrasts for an event_term.
#'
#' Computes F-contrast matrices for main effects and interactions based on the
#' categorical event factors. The function uses kronecker products to combine
#' difference matrices (Dlist) and constant vectors (Clist).
#'
#' @param x An event_term object.
#' @param ... Additional arguments.
#'
#' @return A list of contrast matrices.
#' @export
Fcontrasts.event_term <- function(x, ...) {
  cellcount <- attr(cells(x, drop.empty = FALSE), "count")
  if (any(cellcount) == 0) {
    stop("Currently cannot compute Fcontrasts for non-orthogonal design.")
  }
  
  which_cat <- which(sapply(x$events, function(obj) is_categorical(obj)))
  assert_that(length(which_cat) > 0, msg = "Fcontrasts cannot be computed for terms with no categorical variables")
  
  # Clist: For each categorical event, a vector of ones (constant effect).
  Clist <- lapply(x$events[which_cat], function(ev) rep(1, length(levels(ev))))
  # Dlist: For each categorical event, a difference matrix that captures contrasts.
  Dlist <- lapply(x$events[which_cat], function(ev) t(-diff(diag(length(levels(ev))))))
  
  nfac <- length(Clist)
  valid_cells <- cellcount > 0
  
  # Compute main effects contrasts by combining (via kronecker product)
  # the difference matrix for one factor and constant vectors for others.
  main_effects <- lapply(seq(from = length(Clist), to = 1), function(i) {
    Dcon <- Dlist[[i]]
    Cs <- Clist[-i]
    mats <- vector(nfac, mode = "list")
    mats[[i]] <- Dcon
    mats[setdiff(seq_len(nfac), i)] <- Cs
    ret <- Reduce(kronecker, rev(mats))
    # If some cells are invalid, adjust via SVD or scaling.
    if (!all(valid_cells)) {
      ret <- ret[valid_cells, , drop = FALSE]
      if (ncol(ret) > 1) {
        ret <- svd(ret)$u
      } else {
        ret <- scale(ret, center = TRUE, scale = FALSE)
      }
    }
    row.names(ret) <- conditions(x)
    ret
  })
  
  names(main_effects) <- rev(names(x$events)[which_cat])
  
  # If more than one categorical factor, compute interactions.
  if (length(which_cat) > 1 && all(valid_cells)) {
    interactions <- vector(length(Clist) - 1, mode = "list")
    for (i in seq(from = length(Clist), to = 2)) {
      icomb <- combn(nfac, i)
      ret <- lapply(1:ncol(icomb), function(j) {
        ind <- icomb[, j]
        mats <- vector(nfac, mode = "list")
        mats[ind] <- Dlist[ind]
        if (length(ind) < nfac) {
          mats[-ind] <- Clist[-ind]
        }
        cmat <- Reduce(kronecker, mats)
        row.names(cmat) <- conditions(x)
        cmat
      })
      cnames <- apply(icomb, 2, function(i) paste0(names(x$events)[which_cat][i], collapse = ":"))
      names(ret) <- cnames
      interactions[[i - 1]] <- ret
    }
    
    return(c(main_effects, unlist(interactions, recursive = FALSE)))
  } else {
    main_effects
  }
}

## ============================================================================
## Section 12: Design Matrix Construction
## ============================================================================

#' Construct a design matrix for an event_term.
#'
#' This function creates a design matrix from an event_term object by first
#' populating a local environment with event values and then applying model.matrix
#' to the formula derived from the event_term.
#'
#' @param x An event_term object.
#' @param drop.empty Logical; if TRUE, columns with no events are removed.
#' @param ... Additional arguments.
#'
#' @return A tibble representing the design matrix.
#' @export
design_matrix.event_term <- function(x, drop.empty = TRUE, ...) {
  # Create a new environment to store event elements.
  locenv <- new.env()
  pterms <- purrr::map_chr(parent_terms(x), .sanitizeName)
  
  # Populate the environment with each event's elements.
  for (ev in x$events) {
    vname <- .sanitizeName(ev$varname)
    els <- elements(ev, values = TRUE)
    lapply(names(els), function(n) assign(n, els[[n]], envir = locenv))
  }
  
  # --- MINIMAL PATCH: Make all columns from x$event_table visible in locenv ---
  # This ensures that if the formula references something like `modulator` 
  # (a column of event_table), it's available for model.matrix().
  #browser()
  for (nm in names(x$event_table)) {
    if (!exists(nm, envir = locenv)) {
      assign(nm, x$event_table[[nm]], envir = locenv)
    }
  }
  #browser()
  
  # Create a data frame of event elements.
  els <- as.data.frame(elements(x))
  #print(els)
  nas <- try(apply(els, 1, function(vals) any(is.na(vals))))
  counts <- attr(cells(x, drop.empty = FALSE), "count")
  
  # If the event term consists of a single constant factor, return a constant column.
  mat <- if (ncol(els) == 1 && is.factor(els[, 1]) && length(levels(els[, 1])) == 1) {
    cbind(rep(1, NROW(els))) 
  } else { 		
    out <- try(model.matrix(formula(x), data = locenv))
    if (inherits(out, "try-error")) {
      browser()
    }
    out
  }
  
  rmat <- mat  # Optionally, one might multiply by a subset indicator.
  
  # Replace rows with NA if any are detected.
  if (any(nas)) {
    rmat <- matrix(0, nrow(x$event_table), length(conditions(x, drop.empty = FALSE)))
    rmat[!nas, ] <- mat
    rmat[nas, ] <- NA				
  }
  
  # Optionally remove columns with zero event counts.
  if (any(counts == 0) && (length(conditions(x, drop = FALSE)) == length(counts)) && drop.empty) {
    rmat <- rmat[, !(counts == 0), drop = FALSE]
    colnames(rmat) <- conditions(x, drop = TRUE)
  } else {
    colnames(rmat) <- conditions(x, drop = FALSE)			
  }
  
  suppressMessages(tibble::as_tibble(rmat, .name_repair = "check_unique"))
}

## ============================================================================
## Section 13: Print Methods
## ============================================================================

#' Print fmri_term objects.
#'
#' @param x An fmri_term object.
#' @param ... Additional arguments.
#' @export
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
print.afni_hrf_convolved_term <- function(x, ...) {
  cat("fmri_term: ", class(x)[[1]], "\n")
  cat("  Term Name: ", x$varname, "\n")
  cat("  Formula:  ", as.character(formula(x$evterm)), "\n")
  cat("  Num Events: ", nrow(x$evterm$event_table), "\n")
  cat("  Conditions: ", conditions(x), "\n")
  cat("  Term Types: ", paste(purrr::map_chr(x$evterm$events, ~ class(.)[[1]])), "\n")
}

#' Print event_factor objects.
#'
#' @param x An event_factor object.
#' @param ... Additional arguments.
#' @export
print.event_factor <- function(x, ...) {
  cat("\n═══ Event Factor ═══\n")
  cat("\n📋 Basic Info:\n")
  cat(crayon::blue("  • Name:"), x$varname, "\n")
  cat(crayon::blue("  • Events:"), length(x$value), "\n")
  cat(crayon::blue("  • Levels:"), paste(levels(x$value), collapse = ", "), "\n")
  cat("\n⏱️  Timing:\n")
  cat(crayon::blue("  • Duration range:"), sprintf("%.2f - %.2f seconds", min(x$durations), max(x$durations)), "\n")
  cat(crayon::blue("  • Onset range:"), sprintf("%.2f - %.2f seconds", min(x$onsets), max(x$onsets)), "\n")
  cat("\n🔳 Blocks:\n")
  blocks <- table(x$blockids)
  cat(crayon::blue("  • Number of blocks:"), length(blocks), "\n")
  cat(crayon::blue("  • Events per block:"), paste(blocks, collapse = ", "), "\n\n")
}

#' Print event_variable objects.
#'
#' @param x An event_variable object.
#' @param ... Additional arguments.
#' @export
print.event_variable <- function(x, ...) {
  cat("\n═══ Event Variable ═══\n")
  cat("\n📊 Variable Info:\n")
  cat(crayon::blue("  • Name:"), x$varname, "\n")
  cat(crayon::blue("  • Events:"), length(x$value), "\n")
  cat(crayon::blue("  • Range:"), sprintf("%.2f - %.2f", min(x$value), max(x$value)), "\n")
  cat("\n⏱️  Timing:\n")
  cat(crayon::blue("  • Duration range:"), sprintf("%.2f - %.2f seconds", min(x$durations), max(x$durations)), "\n")
  cat(crayon::blue("  • Onset range:"), sprintf("%.2f - %.2f seconds", min(x$onsets), max(x$onsets)), "\n")
  cat("\n🔳 Blocks:\n")
  blocks <- table(x$blockids)
  cat(crayon::blue("  • Number of blocks:"), length(blocks), "\n")
  cat(crayon::blue("  • Events per block:"), paste(blocks, collapse = ", "), "\n\n")
}

#' Print event_matrix objects.
#'
#' @param x An event_matrix object.
#' @param ... Additional arguments.
#' @export
print.event_matrix <- function(x, ...) {
  cat("\n═══ Event Matrix ═══\n")
  cat("\n📊 Matrix Info:\n")
  cat(crayon::blue("  • Name:"), x$varname, "\n")
  cat(crayon::blue("  • Dimensions:"), paste(dim(x$value), collapse = " × "), "\n")
  cat(crayon::blue("  • Column names:"), paste(colnames(x$value), collapse = ", "), "\n")
  cat("\n⏱️  Timing:\n")
  cat(crayon::blue("  • Duration range:"), sprintf("%.2f - %.2f seconds", min(x$durations), max(x$durations)), "\n")
  cat(crayon::blue("  • Onset range:"), sprintf("%.2f - %.2f seconds", min(x$onsets), max(x$onsets)), "\n")
  cat("\n🔳 Blocks:\n")
  blocks <- table(x$blockids)
  cat(crayon::blue("  • Number of blocks:"), length(blocks), "\n")
  cat(crayon::blue("  • Events per block:"), paste(blocks, collapse = ", "), "\n")
  cat("\n📈 Values:\n")
  ranges <- apply(x$value, 2, range)
  cat(crayon::blue("  • Ranges per column:\n"))
  for (i in 1:ncol(ranges)) {
    cat(sprintf("    %s: %.2f - %.2f\n", colnames(x$value)[i], ranges[1, i], ranges[2, i]))
  }
  cat("\n")
}

#' Print event_term objects.
#'
#' @param x An event_term object.
#' @param ... Additional arguments.
#' @export
print.event_term <- function(x, ...) {
  cat("\n═══ Event Term ═══\n")
  cat("\n📋 Term Info:\n")
  cat(crayon::blue("  • Name:"), x$varname, "\n")
  cat(crayon::blue("  • Number of events:"), nrow(x$event_table), "\n")
  cat(crayon::blue("  • Variables:"), paste(names(x$events), collapse = ", "), "\n")
  cat("\n📊 Variable Types:\n")
  for (name in names(x$events)) {
    type <- class(x$events[[name]])[1]
    cat(sprintf("  • %s: %s\n", name, type))
  }
  cat("\n⏱️  Timing:\n")
  cat(crayon::blue("  • Duration range:"), sprintf("%.2f - %.2f seconds", min(x$durations), max(x$durations)), "\n")
  cat(crayon::blue("  • Onset range:"), sprintf("%.2f - %.2f seconds", min(x$onsets), max(x$onsets)), "\n")
  cat("\n🔳 Blocks:\n")
  blocks <- table(x$blockids)
  cat(crayon::blue("  • Number of blocks:"), length(blocks), "\n")
  cat(crayon::blue("  • Events per block:"), paste(blocks, collapse = ", "), "\n\n")
}

## ============================================================================
## Section 14: End of File
###############################################################################
# End of event_vector.R