
#' @keywords internal
#' @noRd
.sanitizeName <- function(name) {
  name <- gsub(":", ".", name)
  name <- gsub(" ", "", name)
  name <- gsub("\\)$", "", name)
  name <- gsub("[\\(\\)]", "_", name, perl=TRUE)
  name <- gsub(",", "_", name)
  name <- gsub("\\.$", "", name)
  name
}

#' @keywords internal
#' @noRd
is.increasing <- function(vec) {
  all(diff(vec) >= 0)
}

#' @keywords internal
#' @noRd
is.strictly.increasing <- function(vec) {
  all(diff(vec) > 0)
}

#' @keywords internal
#' @import assertthat
#' @noRd
.checkEVArgs <- function(name, vals, onsets, blockids, durations=NULL) {

  assert_that(length(onsets) == length(vals))
  
  ## no NA onsets allowed
  assert_that(all(!is.na(onsets)))
  
  sons <- split(onsets, blockids)

  for (ons in sons) {
    assert_that(is.strictly.increasing(ons))
  }
  
  if (is.null(durations) || length(durations) == 1) {
    durations <- rep(durations, length(onsets))
  }
  
  assert_that(length(durations) == length(vals))
  assert_that(length(blockids) == length(vals))
  
  list(varname=name, value=vals, onsets=onsets, durations=durations, blockids=blockids)
}



#' Create an event model term from a named list of variables.
#'
#' This function generates an event model term from a list of named variables,
#' along with their onsets, block IDs, and durations. Optionally, a subset of
#' onsets can be retained.
#'
#' @param evlist A list of named variables.
#' @param onsets A vector of onset times for the experimental events in seconds.
#' @param blockids A vector of block numbers associated with each onset.
#' @param durations A vector of event durations (default is 1).
#' @param subset A logical vector indicating the subset of onsets to retain (default is NULL).
#' 
#' @return A list containing the following components:
#'   - varname: A character string representing the variable names, concatenated with colons.
#'   - events: A list of event variables.
#'   - subset: A logical vector indicating the retained onsets.
#'   - event_table: A tibble containing event information.
#'   - onsets: A vector of onset times.
#'   - blockids: A vector of block numbers.
#'   - durations: A vector of event durations.
#' 
#' @examples 
#' x1 <- factor(rep(letters[1:3], 10))
#' x2 <- factor(rep(1:3, each=10))
#' eterm <- event_term(list(x1=x1,x2=x2), onsets=seq(1,100,length.out=30), 
#'                     blockids=rep(1,30))
#' 
#' x1 <- rnorm(30)
#' x2 <- factor(rep(1:3, each=10))
#' eterm <- event_term(list(x1=x1,x2=x2), onsets=seq(1,100,length.out=30), 
#'                     blockids=rep(1,30), subset=x1>0)
#'
#' @export
event_term <- function(evlist, onsets, blockids, durations = 1, subset=NULL) {
  
  assert_that(is.increasing(blockids), msg="'blockids' must consist of strictly increasing integers")
            
  if (is.null(subset)) { subset=rep(TRUE, length(onsets)) }
  
  if (length(durations) == 1) {
    durations <- rep(durations, length(onsets))
  }
  
  vnames <- names(evlist)
  evs <- lapply(1:length(evlist), function(i) EV(evlist[[i]], vnames[i], 
                                                 onsets=onsets, 
                                                 blockids=blockids, 
                                                 durations=durations,
                                                 subset=subset))
  
  names(evs) <- sapply(evs, function(ev) ev$varname)

  pterms <- unlist(lapply(evs, function(ev) ev$varname))
  
  len <- sum(subset)
  

  etab <- suppressMessages(tibble::as_tibble(lapply(pterms, function(termname) {
    if (is_continuous(evs[[termname]])) {
      rep(.sanitizeName(termname), len)
    } else {
      evs[[termname]]$value
    }			
  }), .name_repair="check_unique"))
  
  names(etab) <- sapply(pterms, .sanitizeName)
  varname <- paste(sapply(evs, function(x) x$varname), collapse=":")
  
  ret <- list(varname=varname, 
              events=evs, 
              subset=subset, 
              event_table=etab, 
              onsets=evs[[1]]$onsets, 
              blockids=evs[[1]]$blockids, 
              durations=evs[[1]]$durations)
  class(ret) <- c("event_term", "event_seq")
  ret
}


#' @export
event_table.event_term <- function(x) x$event_table


#' EV
#' 
#' factory function for creating 'event' types: event_factor, event_variable, event_basis, event_matrix.
#' 
#' @param vals the event values
#' @param name the name of the event variable
#' @param onsets the event onsets.
#' @param blockids the block ids associated with each event (must be non-decreasing)
#' @param durations the duration of each event.
#' @param subset a \code{logical} vector indicating the subset of events to keep
#' 
#' @examples 
#' 
#' ev_fac <- EV(factor(c("A", "B", "C")), "fac", onsets=c(1,10,20), 
#' blockids=rep(1,3))
#' 
#' ev_fac2 <- EV(factor(c("A", "B", "C")), "fac", onsets=c(1,10,20), 
#' blockids=rep(1,3), subset=c(TRUE, TRUE, FALSE))
#' 
#' ev_numeric <- EV(c(1,2,3), "fac", onsets=c(1,10,20), 
#' blockids=rep(1,3))
#' @keywords internal
#' @noRd
EV <- function(vals, name, onsets, blockids, durations = 1, subset=rep(TRUE,length(onsets))) {
  
  if (length(durations) == 1) {
    durations <- rep(durations, length(onsets))
  }
  
  if (is.matrix(vals) && NCOL(vals) == 1) {
    vals <- vals[, 1, drop=TRUE]
  }
  
  ## subset is eagerly applied here.
  
  if (inherits(vals, "ParametricBasis")) {
    event_basis(vals, onsets, blockids, durations,subset)	
  } else if (is.factor(vals) || is.character(vals)) {
    vals <- factor(as.character(vals)[subset])
    event_factor(vals, name, onsets[subset], blockids[subset], durations[subset])
  }else if (is.numeric(vals)) {
    event_variable(vals[subset], name, onsets[subset], blockids[subset], durations[subset])
  } else if (is.matrix(vals)) {
    event_matrix(vals[subset,], name, onsets[subset], blockids[subset], durations[subset])
  } else {
    stop(paste("cannot create event_seq from type: ", typeof(vals)))
  }
  
}

#' Create a categorical event sequence from a factor
#'
#' This function generates a categorical event sequence object from a given factor.
#' It can be used to create event sequences for categorical predictors in fMRI data analyses.
#'
#' @param fac A factor representing the categorical event sequence
#' @param name A character string providing a name for the event sequence
#' @param onsets A numeric vector of onsets for each event in the sequence
#' @param blockids A numeric vector of block identifiers for each event in the sequence (default: rep(1, length(fac)))
#' @param durations A numeric vector of durations for each event in the sequence (default: rep(0, length(fac)))
#'
#' @return An object representing the categorical event sequence, with class "event_factor" and "event_seq"
#'
#' @examples 
#' efac <- event_factor(factor(c("a", "b", "c", "a", "b", "c")), "abc", onsets=seq(1, 100, length.out=6))
#'
#' @seealso \code{\link{EV}}, \code{\link{event_model}}
#' @export 
event_factor <- function(fac, name, onsets, blockids=rep(1,length(fac)), durations=rep(0, length(fac))) {
  if (!is.factor(fac)) {
    warning("argument 'fac' is not a factor, converting to factor")
    fac <- factor(as.character(fac))
  }
  
  ret <- .checkEVArgs(name, fac, onsets, blockids, durations)
  ret$continuous = FALSE
  #ret$split_onsets <- split(efac$onsets, efac$value)
  class(ret) <- c("event_factor", "event_seq")
  ret
}        

#' Create a continuous valued event sequence from a numeric vector
#'
#' This function generates a continuous valued event sequence object from a given numeric vector.
#' It can be used to create event sequences for continuous predictors in fMRI data analyses.
#'
#' @param vec A numeric vector representing the continuous event sequence values
#' @param name A character string providing a name for the event sequence
#' @param onsets A numeric vector of onsets for each event in the sequence
#' @param blockids A numeric vector of block identifiers for each event in the sequence (default: 1)
#' @param durations A numeric vector of durations for each event in the sequence (default: NULL)
#'
#' @return An object representing the continuous valued event sequence, with class "event_variable" and "event_seq"
#'
#' @examples 
#' evar <- event_variable(c(1, 2, 3, 4, 5, 6), "example_var", onsets=seq(1, 100, length.out=6))
#'
#' @seealso \code{\link{EV}}, \code{\link{event_model}}
#' @export
event_variable <- function(vec, name, onsets, blockids=1, durations=NULL) {
  stopifnot(is.vector(vec))
  
  if (is.factor(vec)) {
    stop("cannot create an event_variable from a factor, use 'event_factor'.")
  }
  
  ret <- .checkEVArgs(name, vec, onsets, blockids, durations)
  ret$continuous <- TRUE
  class(ret) <- c("event_variable", "event_seq")
  ret
  
}       

#' Create a continuous valued event set from a matrix
#'
#' This function generates a continuous valued event set object from a given matrix.
#' It is useful for creating event sequences with multiple continuous predictors for fMRI data analyses.
#'
#' @param mat A matrix of continuous event sequence values, with one row per event
#' @param name A character string providing a name for the event set
#' @param onsets A numeric vector of onsets for each event in the sequence
#' @param blockids A numeric vector of block identifiers for each event in the sequence (default: 1)
#' @param durations A numeric vector of durations for each event in the sequence (default: NULL)
#'
#' @return An object representing the continuous valued event set, with class "event_matrix" and "event_seq"
#'
#' @examples 
#' mat <- matrix(rnorm(200), 100, 2)
#' onsets <- seq(1, 1000, length.out=100)
#' durations <- rep(1, 100)
#' blockids <- rep(1, 100)
#'
#' eset <- event_matrix(mat, "eset", onsets, durations, blockids)
#'
#' @export
event_matrix <- function(mat, name, onsets, blockids=rep(1, ncol(mat)), durations=NULL) {
#  browser()
  stopifnot(is.matrix(mat))
  
  ret <- .checkEVArgs(name, as.vector(mat[,1]), onsets, blockids, durations)
  ret$continuous <- TRUE
  
  if (is.null(colnames(mat))) {
    colnames(mat) <- 1:NCOL(mat)
  }
  
  ret$value <- mat
  
  class(ret) <- c("event_matrix", "event_seq")
  ret
}


#' Create an event set from a ParametricBasis object
#'
#' This function generates an event set object from a given basis object of type `ParametricBasis`.
#' It is useful for creating event sequences based on basis functions for fMRI data analyses.
#'
#' @param basis A ParametricBasis object containing the basis functions
#' @param onsets A numeric vector of onsets for each event in the sequence
#' @param blockids A numeric vector of block identifiers for each event in the sequence (default: 1)
#' @param durations A numeric vector of durations for each event in the sequence (default: NULL)
#' @param subset A logical vector indicating a subset of the basis object to use (default: TRUE for all elements)
#'
#' @return An object representing the event set based on the basis functions, with class "event_basis" and "event_seq"
#'
#' @import assertthat
#' @examples 
#' # Create a ParametricBasis object
#' basis <- ParametricBasis("Gamma", shape = 6, rate = 0.9)
#' onsets <- seq(0, 20, length.out = 21)
#' blockids <- rep(1, length(onsets))
#'
#' # Generate an event_basis object
#' ebasis <- event_basis(basis, onsets, blockids)
#'
#' @export
event_basis <- function(basis, onsets, blockids=1, durations=NULL, subset=rep(TRUE, length(onsets))) {
  assertthat::assert_that(inherits(basis, "ParametricBasis"))
  
  if (any(!subset)) {
    basis <- sub_basis(basis, subset)
  }

  ret <- .checkEVArgs(basis$name, basis$y[,1], onsets[subset], blockids[subset], durations[subset])
  #basis$y <- basis$y[subset,]
  ret$value <- basis$y
  ret$continuous <- TRUE
  ret$basis <- basis
  class(ret) <- c("event_basis", "event_seq")
  ret
}


#' @export
#' @rdname levels
levels.event_factor <- function(x) levels(x$value) 

#' @export
#' @rdname levels
levels.event_variable <- function(x) x$varname 

#' @export
#' @rdname levels
levels.event_matrix <- function(x) colnames(x$value) 

#' @export
#' @rdname levels
levels.event_set <- function(x) colnames(x$value) 

#' @export
levels.event_basis <- function(x) levels(x$basis)

#' @export
formula.event_term <- function(x, ...) as.formula(paste("~ ", "(", paste(parent_terms(x), collapse=":"), 
                                                   "-1", ")"))

#' @export
#' @rdname levels
levels.event_term <- function(x) {
  facs <- x$events[!sapply(x$events, is_continuous)]
  if (length(facs) == 1) {
    levels(facs[[1]])
  } else {
    facs <- lapply(facs, function(x) x$value)
    f <- function(...) interaction(..., drop=TRUE, sep=":")
    levels(do.call(f, facs))
  }
}


#' Retrieve cells of an event_factor object.
#'
#' This function extracts cells from an event_factor object, optionally removing
#' empty cells. Note that the removal of empty cells is not implemented.
#'
#' @param x An event_factor object.
#' @param drop.empty Logical. If TRUE, empty cells will be removed (not implemented).
#' @param ... Additional arguments to be passed to the function.
#'
#' @return A list of data frames containing the cells of the event_factor object.
#'
#' @export
#' @rdname cells
cells.event_factor <- function(x, drop.empty=TRUE,...) {
  etab <- data.frame(onsets=x$onsets, durations=x$durations, blockids=x$blockids)
  split(etab, x$value)
}


#' Retrieve cells of an event_term object.
#'
#' This function extracts cells from an event_term object, optionally removing
#' empty cells.
#'
#' @param x An event_term object.
#' @param drop.empty Logical. If TRUE, empty cells will be removed.
#' @param ... Additional arguments to be passed to the function.
#'
#' @return A tibble containing the cells of the event_term object, with an
#'   additional "count" attribute.
#'
#' @export
#' @rdname cells
#' @examples
#' 
#' evlist <- list(fac1=factor(c("A", "B", "A", "B")), 
#'                fac2=factor(c("1", "1", "2", "2")))
#' eterm <- event_term(evlist,onsets=1:4, blockids=rep(1,4))
#' cells(eterm)
cells.event_term <- function(x, drop.empty=TRUE,...) {
  evtab <- x$event_table
  evset <- suppressMessages(tibble::as_tibble(expand.grid(lapply(x$events, levels)), .name_repair="check_unique"))
  
  which.cat <- which(!sapply(x$events, is_continuous))
  
  if (length(which.cat) > 0) {
    evs <- suppressMessages(tibble::as_tibble(lapply(evset[,which.cat], as.character), .name_repair="check_unique"))
    evt <- suppressMessages(tibble::as_tibble(lapply(evtab[,which.cat], as.character), .name_repair="check_unique"))
  
    counts <- apply(evs, 1, function(row1) {
      sum(apply(evt, 1, function(row2) {										
        all(as.character(row1) == as.character(row2))
      }))
    })
  
    if (drop.empty) {
      evset <- evset[counts > 0,,drop=F]
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
#' @keywords internal
.event_set <- function(x, exclude_basis=FALSE) {
  evtab <- event_table(x)
  
  evset <- if (nbasis(x) > 1 & !exclude_basis) {
    evlist <- c(list(factor(paste("basis", 1:nbasis(x), sep = ""))), cells(x$evterm))
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
cells.covariate_convolved_term <- function(x,...) {
  unique(event_table(x))
}

#' @export
#' @importFrom stringr str_trim
cells.convolved_term <- function(x, exclude_basis=FALSE,...) {
  evtab <- event_table(x)
  evset <- .event_set(x, exclude_basis=exclude_basis)
  
  strels <- apply(apply(evtab, 2, stringr::str_trim), 1, paste, collapse = ":")
  
  strlevs <- if (nrow(evset) > 1) {
    apply(apply(evset, 2, stringr::str_trim), 1, paste, collapse = ":")
  } else {
    as.character(evset[1,1])
  }
  
  attr(evset, "rownames") <- strlevs
  
  counts <- if (exclude_basis) {
    rep(attr(cells(x$evterm), "count"), each = 1)
  } else {
    rep(attr(cells(x$evterm), "count"), each = nbasis(x))
  }
  
  ret <- evset[counts > 0, , drop = F]
  attr(ret, "count") <- counts[counts > 0]
  ret
  
}

#' @export
#' @rdname conditions
conditions.fmri_term <- function(x, ...) {
  colnames(design_matrix(x))
}

#' @export
#' @rdname conditions
conditions.convolved_term <- function(x,...) {
  colnames(design_matrix(x))
}

#' @export
#' @rdname conditions
conditions.afni_hrf_convolved_term <- function(x,...) {
  conditions(x$evterm)
}

#' @export
#' @rdname conditions
conditions.afni_trialwise_convolved_term <- function(x,...) {
  conditions(x$evterm)
}


#' @export
#' @rdname conditions
conditions.event_term <- function(x, drop.empty=TRUE, ...) {
  
  .cells <- cells(x, drop.empty=drop.empty)
  pterms <- parent_terms(x)
  levs <- apply(.cells, 1, paste, collapse=":")
  
  splitlevs <- strsplit(levs, ":")
  ret <- lapply(1:length(pterms), function(i) {
    lev <- sapply(splitlevs, "[[", i)
    term <- pterms[[i]]
    
    if (length(levels(x$events[[i]])) > 1) {
      paste(.sanitizeName(pterms[i]), "[", lev, "]", sep="")
    } else {
      .sanitizeName(pterms[i])
    }
  })
  
  do.call(function(...) paste(..., sep=":"), ret)
}


## TODO is columns ever used? do we need this function

#' @export
#' @rdname columns
columns.event_term <- function(x) as.vector(unlist(lapply(x$events, columns)))

#' @export
#' @rdname columns
columns.event_seq <- function(x) x$varname

#' @export
#' @rdname columns
columns.event_matrix <- function(x) paste0(.sanitizeName(x$varname), ".", levels(x))

#' @export
#' @rdname columns
columns.event_set <- function(x) paste0(.sanitizeName(x$varname), ".", levels(x))

#' @export
#' @rdname columns
columns.event_basis <- function(x) columns(x$basis)


#' @export
#' @rdname parent_terms
parent_terms.event_term <- function(x) unlist(lapply(x$events, function(ev) ev$varname))


#' @export
is_continuous.event_seq <- function(x) x$continuous

#' @export
is_continuous.event_factor <- function(x) x$continuous

#' @export
is_categorical.event_seq <- function(x) !x$continuous

#' @export
is_continuous.event_term <- function(x) all(sapply(x$events, function(x) is_continuous(x)))

#' @export
is_categorical.event_term <- function(x) !is_continuous(x)


#' @export
is_categorical.event_seq <- function(x) !x$continuous


#' @export
elements.event_matrix <- function(x, values=TRUE, ...) {
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
    #ret <- list(rep(varname(x), length(x)))
    ret <- list(rep(x$varname, length(x)))
    names(ret) <- x$varname
    ret
  }
}

#' @export
elements.event_basis <- function(x, values=TRUE, transformed=TRUE, ...) {
  if (values && !transformed) {
    x$value$x				
  } else if (values) {
    #browser()
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
elements.event_term <- function(x, values=TRUE, ...) {
  els <- lapply(x$events, elements, values=values)
  n <- sapply(names(els), function(nam) .sanitizeName(nam))
  names(els) <- as.vector(n)
  els
}

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

#' Split onsets of an event_term object
#'
#' This function splits the onsets of an event_term object based on its factor levels or block identifiers.
#' It is useful for processing fMRI data when analyzing event-related designs with multiple conditions or blocks.
#'
#' @param x An event_term object
#' @param sframe A data frame representing the sampling frame
#' @param global A logical value indicating whether to use global onsets (default: FALSE)
#' @param blocksplit A logical value indicating whether to split onsets by block identifiers (default: FALSE)
#' @param ... Additional arguments passed to other methods
#'
#' @return A list of numeric vectors representing the split onsets for each factor level or block identifier
#'
#' @export
split_onsets.event_term <- function(x, sframe, global=FALSE,blocksplit=FALSE, ...) {
  ### need to check for 0 factors
  facs <- x$events[!sapply(x$events, is_continuous)]
  
  if (length(facs) == 0) {
    ons <- if (global) {
      global_onsets(sframe, onsets(x), blockids(x))
    } else {
      onsets(x)
    }
    return(list(split(ons, blockids(x))))
    
  }
  
  facs <- lapply(facs, function(fac) unlist(elements(fac)))
  
  f <- function(...) {
    interaction(..., drop=TRUE, sep=":")
  }
            
  cfac <- try(do.call(f, facs))
  #if (inherits(cfac, "try-error")) {
  #  stop("could not construct crossed factor for event_term `x`")
  #}
            
  ret <- if (global) {
    split(global_onsets(sframe, onsets(x), blockids(x)), cfac)
  } else {
    split(onsets(x), cfac)
  }
  
  if (blocksplit) {
    bsplit <- split(blockids(x), cfac)
    ret <- lapply(1:length(ret), function(i) {
      split(ret[[i]], bsplit[[i]])
    })
  }
  
  names(ret) <- longnames(x)
  ret
}




#' Convolve HRF with Design Matrix
#'
#' This function convolves a given hemodynamic response function (HRF) with the design matrix of an fMRI study.
#' It is useful for modeling the expected BOLD signal in response to experimental conditions.
#'
#' @param hrf A numeric vector representing the hemodynamic response function
#' @param dmat A design matrix with columns representing different experimental conditions
#' @param globons A numeric vector of global onsets for each event
#' @param durations A numeric vector of event durations
#' @param summate A logical value indicating whether to summate the convolved HRF (default: TRUE)
#'
#' @return A list of regressors, one for each column in the design matrix.
#' @export
convolve_design <- function(hrf, dmat, globons, durations, summate=TRUE) {
  cond.names <- names(dmat)
  
  if (any(is.na(dmat)) || any(is.na(globons))) {
    keep <- apply(dmat, 1, function(vals) all(!is.na(vals)))
    keep[is.na(globons)] <- FALSE
    dmat <- dmat[keep,]
    durations <- durations[keep]
    globons <- globons[keep]
  } 
  
  reglist <- purrr::map(1:ncol(dmat), function(i) {
    amp <- dmat[,i][[1]]
    nonzero <- which(amp != 0)
    ## issue with single_trial_regressor
    if (length(nonzero) == 0) {
      null_regressor(hrf)
      ## scaling issue with single_trial_regressor
      #} #else if (length(nonzero) == 1) {
      #single_trial_regressor(globons[nonzero], hrf, amplitude=amp[nonzero], duration=durations[nonzero])
    } else {
      regressor(globons[nonzero], hrf, amplitude=amp[nonzero], duration=durations[nonzero], summate=summate)
    }
  })
  
  reglist
  
}

#' Convolve an event-related design matrix with an HRF.
#'
#' This function takes an event-related design matrix and convolves it with
#' a specified Hemodynamic Response Function (HRF) to create a new design matrix
#' suitable for fMRI analysis. It also supports additional arguments for
#' flexibility and customization.
#'
#' @importFrom tibble as_tibble
#' @importFrom magrittr %>%
#' @importFrom dplyr group_by select do ungroup
#' @autoglobal
#' @export
#'
#' @param x A data frame containing the input design matrix.
#' @param hrf A Hemodynamic Response Function to convolve the design matrix with.
#' @param sampling_frame A data frame specifying the sampling frame for the analysis.
#' @param drop.empty Logical. If TRUE, empty rows in the design matrix will be removed.
#' @param summate Logical. If TRUE, the convolved design matrix will be summed.
#' @param precision Numeric. The desired precision for the calculations.
#' @param ... Additional arguments to be passed to the function.
#'
#' @return A convolved design matrix, in tibble format.
convolve.event_term <- function(x, hrf, sampling_frame, drop.empty=TRUE, 
                                summate=TRUE, precision=.3,...) {
  
  
  globons <- global_onsets(sampling_frame, x$onsets, x$blockids)
  
  durations <- x$durations
  blockids <- x$blockids
  
  nimages <- sum(sampling_frame$blocklens)
  
  cnames <- conditions(x)
  dmat <- design_matrix(x, drop.empty)
  ncond <- ncol(dmat)

  #print("convolving")
 
  cmat <- dmat |> dplyr::mutate(.blockids=blockids, .globons=globons, .durations=durations) |> 
    dplyr::group_by(.blockids) %>%
    dplyr::do({
      d <- dplyr::select(., 1:ncond)
      reg <- convolve_design(hrf, d, .$.globons, .$.durations, summate=summate)
      sam <- samples(sampling_frame, blockids=as.integer(as.character(.$.blockids[1])), global=TRUE)
      
      ## TODO could bee parallelized
      ## ret <- do.call(rbind, furrr::future_map(reg, function(r) evaluate(r, sam))) 
      ret <- do.call(cbind, lapply(seq_along(reg), function(ri) {
        vname <- paste0("v", ri)
        evaluate(reg[[ri]], sam, precision=precision) 
        
        #%>% select({{vname}} := value)
        #names(tmp) <- paste0("v", ri)
        #suppressMessages(as_tibble(tmp))
      })) 
      
      ret <- suppressMessages(tibble::as_tibble(ret, .name_repair="minimal"))
      names(ret) <- paste0("v", 1:length(reg))
      ret
  }) 
  
  cmat <- cmat %>% dplyr::ungroup() %>% dplyr::select(-.blockids)
  
  if (nbasis(hrf) > 1) {
    blevs <- paste("[", 1:nbasis(hrf), "]", sep="")
    cnames <- unlist(lapply(cnames, function(prefix) paste(prefix, ":basis", blevs, sep="")))
  } 
            
  
  colnames(cmat) <- cnames
  suppressMessages(tibble::as_tibble(cmat, .name_repair="check_unique"))
  
  
}

#' Compute F-contrasts for Event Term
#'
#' This function computes F-contrasts for an event term object, considering main effects and interactions.
#'
#' @param x An event term object
#' @param ... Additional arguments passed to the function
#'
#' @return A list of contrast matrices for main effects and interactions
#' @export
Fcontrasts.event_term <- function(x,...) {
  cellcount <- attr(cells(x, drop.empty=FALSE), "count")
  if (any(cellcount) == 0) {
    stop("currently cannot compute Fcontrasts for non-orthogonal design.")
  }
  ##browser()
  ## TODO check for no empty cells, otherwise everything fails
  which_cat <- which(sapply(x$events, function(obj) is_categorical(obj)))
  assert_that(length(which_cat) > 0, msg="Fcontrasts cannot be computed for terms with no categorical variables")
  ## factors comprising this term
  pterms <- parent_terms(x)[which_cat]
  evs <- x$events[which_cat]
  cond <- conditions(x)
  facnames <- names(evs)
  
  Clist <- lapply(evs, function(ev) rep(1, length(levels(ev))))
  Dlist <- lapply(evs, function(ev) t(-diff(diag(length(levels(ev))))))
  
  nfac <- length(Clist)
  
  valid_cells <- cellcount > 0

  main_effects <- lapply(length(Clist):1, function(i) {
    #print(i)
    Dcon <- Dlist[[i]]
    Cs <- Clist[-i]
    mats <- vector(nfac, mode="list")
    mats[[i]] <- Dcon
    mats[seq(1, nfac)[-i]] <- Cs
    ret <- Reduce(kronecker, rev(mats))
    if (!(all(valid_cells))) {
      ret <- ret[valid_cells,,drop=FALSE]
      if (ncol(ret) > 1) {
        ret <- svd(ret)$u
      } else {
        ret <- scale(ret, center=TRUE, scale=FALSE)
      }
    }
    row.names(ret) <- cond
    ret
  })
  
  names(main_effects) <- rev(facnames)
  
  if (length(facnames) > 1 && all(valid_cells)) {
    interactions <- vector(length(Clist)-1, mode="list")
    for (i in length(Clist):2) {
      icomb <- combn(nfac, i)
      ret <- lapply(1:ncol(icomb), function(j) {
        ind <- icomb[,j]
        mats <- vector(nfac, mode="list")
        mats[ind] <- Dlist[ind]
        if (length(ind) < nfac) {
          mats[-ind] <- Clist[-ind]
        }
      
        cmat <- Reduce(kronecker, mats)
        row.names(cmat) <- cond
        cmat
      })
      cnames <- apply(icomb, 2, function(i) paste0(facnames[i], collapse=":"))
      names(ret) <- cnames
      interactions[[i-1]] <- ret
    }
    
    return(c(main_effects, unlist(interactions, recursive=FALSE)))
  } else {
    main_effects
  }
  
}
                 

#' @importFrom tibble as_tibble
#' @importFrom purrr map_chr
#' @export
design_matrix.event_term <- function(x, drop.empty=TRUE,...) {

  locenv <- new.env()
  pterms <- map_chr(parent_terms(x), .sanitizeName)	
  
  for (ev in x$events) {
    vname <- .sanitizeName(ev$varname)
    els <- elements(ev, values=TRUE)
    #if (length(els) == 1 && is.matrix(els[[1]])) {
    ## TODO remove special casing for matrix elements 
    #  mat <- els[[1]]
    #  for (i in 1:ncol(mat)) {
    #    assign(columns(ev)[i], mat[,i],envir=locenv)
    #  }
    #} else {
      lapply(names(els), function(n) assign(n, els[[n]],envir=locenv))
    #}
  }
  
  
  els <- as.data.frame(elements(x))
  nas <- try(apply(els,1, function(vals) any(is.na(vals))))
  counts <- attr(cells(x, drop=FALSE), "count")
  
  #print(ncol(els))
  #browser()
  mat <- if (ncol(els) == 1 && is.factor(els[,1]) && length(levels(els[,1])) == 1) {
    ## a 1 level term
    cbind(rep(1, NROW(els))) 
  } else { 		
    model.matrix(formula(x), data=locenv) 
  }
  

  ### multiply design matrix by subset
  rmat <- mat #* x$subset
  
  #remove rows with NAS
  if (any(nas)) {
    rmat <- matrix(0, nrow(x$event_table), length(conditions(x, drop.empty=FALSE)))
    rmat[!nas,] <- mat
    rmat[nas,] <- NA				
  } 
  

  # remove columns with no events (postpone this to later stage?) 
  if (any(counts == 0) && (length(conditions(x, drop=F)) == length(counts)) && drop.empty) {
    rmat <- rmat[, !(counts==0), drop=FALSE]
    colnames(rmat) <- conditions(x, drop=T)
  } else {
    colnames(rmat) <- conditions(x, drop=F)			
  }
  
  suppressMessages(tibble::as_tibble(rmat, .name_repair="check_unique"))
}



#' @export
print.event_term <- function(x, ...) {
  cat("event_term", "\n")
  cat("  ", "Term Name: ", x$varname, "\n")
  cat("  ", "Formula:  ", as.character(formula(x)), "\n")
  cat("  ", "Num Events: ", nrow(x$event_table), "\n")
  cat("  ", "Term Types: ", paste(map_chr(x$events, ~ class(.)[[1]])))
  cat("\n")
}

#' @export
print.fmri_term <- function(x,...) {
  cat("fmri_term: ", class(x)[[1]], "\n")
  cat("  ", "Term Name: ", x$varname, "\n")
  cat("  ", "Num Rows: ", nrow(design_matrix(x)), "\n")
  cat("  ", "Num Columns: ", ncol(design_matrix(x)), "\n")
}

#' @export
print.convolved_term <- function(x,...) {
  cat("fmri_term: ", class(x)[[1]], "\n")
  cat("  ", "Term Name: ", x$varname, "\n")
  cat("  ", "Formula:  ", as.character(formula(x$evterm)), "\n")
  cat("  ", "Num Events: ", nrow(x$evterm$event_table), "\n")
  cat("  ", "Num Rows: ", nrow(design_matrix(x)), "\n")
  cat("  ", "Num Columns: ", ncol(design_matrix(x)), "\n")
  cat("  ", "Conditions: ", conditions(x), "\n")
  cat("  ", "Term Types: ", paste(map_chr(x$evterm$events, ~ class(.)[[1]])))
  cat("\n")
}

#' @export
print.afni_hrf_convolved_term <- function(x,...) {
  cat("fmri_term: ", class(x)[[1]], "\n")
  cat("  ", "Term Name: ", x$varname, "\n")
  cat("  ", "Formula:  ", as.character(formula(x$evterm)), "\n")
  cat("  ", "Num Events: ", nrow(x$evterm$event_table), "\n")
  cat("  ", "Conditions: ", conditions(x), "\n")
  cat("  ", "Term Types: ", paste(map_chr(x$evterm$events, ~ class(.)[[1]])))
  cat("\n")
}




