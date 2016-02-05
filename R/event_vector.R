#' @import checkmate
NULL

.sanitizeName <- function(name) {
  name <- gsub(":", ".", name)
  name <- gsub(" ", "", name)
  name <- gsub("[\\(\\)]", ".", name, perl=TRUE)
  name <- gsub(",", ".", name)
  name <- gsub("\\.$", "", name)
  name
}

is.increasing <- function(vec) {
  all(diff(vec) >= 0)
}

is.strictly.increasing <- function(vec) {
  all(diff(vec) > 0)
}


#' @import assertthat
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



#' event_term
#' 
#' Create a event model term from a named list of variables.
#' @param evlist
#' @param onsets the onset times from the experimental events in seconds
#' @param blockids the block number associated with each onset
#' @param durations
#' @param subset
#' @export
#' @rdname event_term-class
event_term <- function(evlist, onsets, blockids, durations = 1, subset=NULL) {
  assert_that(is.increasing(blockids))
  
  vnames <- names(evlist)
  evs <- lapply(1:length(evlist), function(i) EV(evlist[[i]], vnames[i], onsets=onsets, blockids=blockids, durations=durations))
  names(evs) <- sapply(evs, function(ev) ev$varname)
  
  if (is.null(subset)) { subset=rep(TRUE, length(evs[[1]]$onsets)) }
  
  pterms <- unlist(lapply(evs, function(ev) ev$varname))
  
  etab <- as.data.frame(lapply(pterms, function(termname) {
    if (isContinuous(evs[[termname]])) {
      rep(.sanitizeName(termname), length(onsets))
    } else {
      evs[[termname]]$value
    }			
  }))
  
  names(etab) <- sapply(pterms, .sanitizeName)
  varname <- paste(sapply(evs, function(x) x$varname), collapse=":")
  
  ret <- list(varname=varname, events=evs, subset=subset, eventTable=etab, onsets=evs[[1]]$onsets, blockids=evs[[1]]$blockids, durations=evs[[1]]$durations)
  class(ret) <- c("event_term", "event_seq")
  ret
}

event_table.event_term <- function(x) x$eventTable


#' EV
#' 
#' factory function for creating 'event' types: event_factor, event_variable, event_basis, event_matrix.
#' 
#' @param vals the event values
#' @param name the name of the event variable
#' @param onsets the event onsets.
#' @param blockids the block ids associated with each event (must be non-decreasing)
#' @export
EV <- function(vals, name, onsets, blockids, durations = 1) {
  
  if (length(durations) == 1) {
    durations <- rep(durations, length(onsets))
  }
  
  if (is.matrix(vals) && NCOL(vals) == 1) {
    vals <- vals[, 1, drop=TRUE]
  }
  
  if (inherits(vals, "ParametricBasis")) {
    event_basis(vals, onsets, blockids, durations)	
  } else if (is.vector(vals)) {
    event_variable(vals, name, onsets, blockids, durations)
  } else if (is.matrix(vals)) {
    event_matrix(vals, name, onsets, blockids, durations)
  } else if (is.factor(vals)) {
    event_factor(vals, name, onsets, blockids, durations)
  } else {
    stop(paste("cannot create event_seq from type: ", typeof(vals)))
  }
  
}

#' event_factor
#' 
#' Create an categorical event sequence from a \code{factor} 
#' 
#' @param fac a factor
#' @param name the name for the factor
#' @param onsets vector of event onsets in seconds
#' @param blockids block index variable
#' @param durations the durations in seconds of the onsets
#' @export
event_factor <- function(fac, name, onsets, blockids=1, durations=NULL) {
  if (!is.factor(fac)) {
    warning("argument 'fac' is not a factor, converting to factor")
    fac <- as.factor(factor())
  }
  
  ret <- .checkEVArgs(name, fac, onsets, blockids, durations)
  ret$continuous = FALSE
  class(ret) <- c("event_factor", "event_seq")
  ret
}        

#' event_variable
#' 
#' Create a continuous valued event sequence from a \code{numeric} vector.
#' @param name the name of the variable
#' @param onsets the event onsets in seconds
#' @param blockids the index of the block/scan in which the event occurs
#' @param durations the durations of each event in seconds
#' @export
event_variable <- function(vec, name, onsets, blockids=1, durations=NULL) {
  stopifnot(is.vector(vec))
  
  if (is.factor(vec)) {
    stop("cannot create an event_variable from a factor, use EventFactor.")
  }
  
  ret <- .checkEVArgs(name, vec, onsets, blockids, durations)
  ret$continuous <- TRUE
  class(ret) <- c("event_variable", "event_seq")
  ret
  
}       

#' event_matrix
#' 
#' Create a continuous valued event set from a \code{matrix}
#' 
#' @param mat a matrix of values, one row per event, indicating the amplitude/intensity of each event.
#' @param name the name of the variable
#' @param onsets the event onsets in seconds
#' @param durations the durations of each event in seconds
#' @param blockids the index of the block/scan in which the event occurs
#' @examples 
#' 
#' mat <- matrix(rnorm(200), 100, 2)
#' onsets <- seq(1, 1000, length.out=100)
#' durations <- rep(1, 100)
#' blockids <- rep(1, 100)
#' 
#' eset <- event_matrix(mat, "eset", onsets,durations,blockids)
#' 
#' @export
event_matrix <- function(mat, name, onsets, durations=NULL, blockids=1 ) {
  stopifnot(is.matrix(mat))
  
  ret <- .checkEVArgs(name, as.vector(mat[,1]), onsets, blockids, durations)
  ret$continuous <- TRUE
  
  if (is.null(colnames(mat))) {
    colnames(mat) <- 1:NCOL(mat)
  }
  
  class(ret) <- c("event_matrix", "event_seq")
  ret
}


#' event_basis
#' 
#' Create a event set from a basis object of type \code{\linkS4class{ParametricBasis}}. 
#' @param basis
#' @param onsets
#' @param blockids
#' @param durations
#' @import assertthat
#' @export
event_basis <- function(basis, onsets, blockids=1, durations=NULL) {
  assertthat::assert_that(inherits(basis, "ParametricBasis"))
  ret <- .checkEVArgs(basis$name, basis$x, onsets, blockids, durations)
  ret$continuous <- TRUE
  ret$basis <- basis
  class(ret) <- c("event_basis", "event_seq")
  ret
}


#' @export
levels.event_factor <- function(x) levels(x$value) 

#' @export
levels.event_variable <- function(x) x$varname 


levels.event_matrix <- function(x) colnames(x$value) 


#' @export
levels.event_set <- function(x) colnames(x$value) 

#' @export
levels.event_basis <- function(x) seq(1, ncol(x$basis$y))

#' @export
formula.event_term <- function(x) as.formula(paste("~ ", "(", paste(parentTerms(x), collapse=":"), "-1", ")"))

levels.event_term <- function(x) {
  facs <- x$events[!sapply(x$events, isContinuous)]
  if (length(facs) == 1) {
    levels(facs[[1]])
  } else {
    facs <- lapply(facs, function(x) x$value)
    f <- function(...) interaction(..., drop=TRUE, sep=":")
    levels(do.call(f, facs))
  }
}

#' @export
cells.event_factor <- function(x, drop.empty=TRUE) {
  etab <- data.frame(onsets=x$onsets, durations=x$durations, blockids=x$blockids)
  split(etab, x$value)
}

#' @export
cells.event_term <- function(x, drop.empty=TRUE) {
  evtab <- x$eventTable
  evset <- expand.grid(lapply(x$events, levels))
  which.cat <- which(!sapply(x$events, isContinuous))
  
  counts <- apply(evset[,which.cat,drop=F], 1, function(row1) {
    sum(apply(evtab[x$subset,which.cat,drop=F], 1, function(row2) {										
      all(row1 == row2)
    }))
  })
  
  if (drop.empty) {
    evset <- evset[counts > 0,,drop=F]
    attr(evset, "count") <- counts[counts > 0]
  } else {
    attr(evset, "count") <- counts
  }
  evset
}

.event_set <- function(x) {
  evtab <- event_table(x)
  
  evset <- if (nbasis(x) > 1) {
    evlist <- c(list(factor(paste("basis", 1:nbasis(x), sep = ""))), cells(x@eventTerm))
    names(evlist) <- c("basis", parentTerms(x@eventTerm))
    evlist <- lapply(evlist, levels)
    ret <- expand.grid(evlist, stringsAsFactors = TRUE)
    ret[c(2:length(ret), 1)]
  } else {
    cells(x$evterm)
  }
  
}

cells.convolved_term <- function(x) {
  evtab <- event_table(x)
  evset <- .event_set(x)
  
  strels <- apply(apply(evtab, 2, str_trim), 1, paste, collapse = ":")
  strlevs <- apply(apply(evset, 2, str_trim), 1, paste, collapse = ":")
  row.names(evset) <- strlevs
  counts <- rep(attr(cells(x$evterm), "count"), each = nbasis(x))
  
  ret <- evset[counts > 0, , drop = F]
  attr(ret, "count") <- counts[counts > 0]
  ret
  
}

conditions.fmri_term <- function(x) {
  colnames(design_matrix(x))
}


#' @export
conditions.event_term <- function(x, drop.empty=TRUE) {
  
  .cells <- cells(x, drop.empty=drop.empty)
  pterms <- parentTerms(x)
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

#' @export
columns.event_term <- function(x) as.vector(unlist(lapply(x$events, columns)))

#' @export
columns.event_seq <- function(x) x$varname

columns.event_matrix <- function(x) paste0(.sanitizeName(x$varname), ".", levels(x))


#' @export
columns.event_set <- function(x) paste0(.sanitizeName(x$varname), ".", levels(x))

#' @export
columns.event_basis <- function(x) paste0(.sanitizeName(x$varname), ".", levels(x))


#' @export
parentTerms.event_term <- function(x) unlist(lapply(x$events, function(ev) ev$varname))

#' @export
isContinuous.event_seq <- function(x) x$continuous



#' @export
elements.event_matrix <- function(x, values=TRUE) {
  if (values) {
    ret <- x$value
    colnames(ret) <- colnames(x)
    ret <- list(ret)
    names(ret) <- .sanitizeName(varname(x))
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
elements.event_seq <- function(x, values = TRUE) {
  if (values) {
    ret <- list(x$value)
    names(ret) <- x$varname
    ret
  } else {
    ret <- list(rep(varname(x), length(x)))
    names(ret) <- varname(x)
    ret
  }
}

#' @export
elements.event_basis <- function(x, values=TRUE, transformed=TRUE) {
  if (values && !transformed) {
    x$value$x				
  } else if (values) {
    ret <- x$basis$y
    colnames(ret) <- colnames(x)
    n <- .sanitizeName(x$varname)
    ret <- list(ret)
    names(ret) <- n
    ret
  } else {
    N <- length(x)
    vnames <- colnames(x)
    res <- lapply(vnames, function(el) rep(el, N))
    mat <- do.call(cbind, res)
    colnames(mat) <- vnames			
    ret <- list(mat)
    names(ret) <- .sanitizeName(varname(x))
    ret		
  }
}


#' @export
elements.event_term <- function(x, values=TRUE) {
  els <- lapply(x$events, elements, values=values)
  n <- sapply(names(els), function(nam) .sanitizeName(nam))
  names(els) <- as.vector(n)
  els
}


#' @export
convolve_design <- function(hrf, dmat, globons, durations) {
  cond.names <- names(dmat)
  keep <- if (any(is.na(dmat)) || any(is.na(globons))) {
    ret <- apply(dmat, 1, function(vals) all(!is.na(vals)))
    ret[is.na(globons)] <- FALSE
  } else {
    TRUE
  }
  
  
  lapply(1:NCOL(dmat), function(i) {
    regressor(globons[keep], hrf, amplitude=dmat[keep,i], duration=durations[keep])
  })
  
}

#' @export
convolve.event_term <- function(x, hrf, samplingFrame, drop.empty=TRUE) {

  globons <- global_onsets(samplingFrame, x$onsets, x$blockids)
  
  nimages <- sum(samplingFrame$blocklens)
  samples <- seq(samplingFrame$startTime, length.out=nimages, by=samplingFrame$TR)
  
  dmat <- as.data.frame(design_matrix(x, drop.empty))
  
  blockids <- factor(x$blockids)
  split.dmat <- split(dmat, blockids)
  split.ons <- split(globons, blockids)
  split.durations <- split(x$durations, blockids)
  split.samples <- split(samples, rep(1:length(samplingFrame$blocklens), samplingFrame$blocklens))
  
  reglist <- lapply(1:length(split.dmat), function(i) {
    reg <- convolve_design(hrf, split.dmat[[i]], split.ons[[i]], split.durations[[i]])
    do.call(cbind, lapply(reg, function(r) evaluate(r, split.samples[[i]])))
  })
  
  ret <- do.call(rbind, reglist)
  
  cnames <- conditions(x)
 
  if (nbasis(hrf) > 1) {
    blevs <- paste("[", 1:nbasis(hrf), "]", sep="")
    cnames <- unlist(lapply(cnames, function(prefix) paste(prefix, ":basis", blevs, sep="")))
  } 
            
  colnames(ret) <- cnames
  as.data.frame(ret)
  
  #lapply(reglist, function(reg) evaluate(reg, )
  
}

#' @export
design_matrix.event_term <- function(x, drop.empty=TRUE) {
  locenv <- new.env()
  pterms <- sapply(parentTerms(x), .sanitizeName)	
  
  for (ev in x$events) {
    vname <- .sanitizeName(ev$varname)
    els <- elements(ev, values=TRUE)
    lapply(names(els), function(n) assign(n, els[[n]],envir=locenv))			
  }
  
  els <- as.data.frame(elements(x))
  nas <- try(apply(els,1, function(vals) any(is.na(vals))))
  counts <- attr(cells(x, drop=FALSE), "count")
  
  mat <- if (ncol(els) == 1 && is.factor(els[,1]) && length(levels(els[,1])) == 1) {
    ## a 1 level term
    cbind(rep(1, NROW(els))) 
  } else { 		
    model.matrix(formula(x), data=locenv) 
  }
  
  ### multiply design matrix by subset
  rmat <- mat * x$subset
  
  #remove rows with NAS
  if (any(nas)) {
    rmat <- matrix(0, nrow(x$eventTable), length(conditions(x, drop.empty)))
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
  
  rmat
}



#' @export
print.event_term <- function(object) {
  cat("event_term", "\n")
  cat("  ", "Term Name: ", object$varname, "\n")
  cat("  ", "Formula:  ", as.character(formula(object)), "\n")
  cat("  ", "Num Events: ", nrow(object$eventTable), "\n")
  cat("  ", "Term Types: ", paste(sapply(object$events, function(ev) class(ev)[[1]])))
}

#' @export
print.fmri_term <- function(object) {
  cat("fmri_term", "\n")
  cat("  ", "Term Name: ", object$varname, "\n")
  cat("  ", "Num Events: ", nrow(design_matrix(object)), "\n")
  cat("  ", "Num Columns: ", ncol(design_matrix(object)), "\n")
}

#' @export
print.convolved_term <- function(object) {
  cat("fmri_term", "\n")
  cat("  ", "Term Name: ", object$varname, "\n")
  cat("  ", "Formula:  ", as.character(formula(object$evterm)), "\n")
  cat("  ", "Num Events: ", nrow(design_matrix(object)), "\n")
  cat("  ", "Num Columns: ", ncol(design_matrix(object)), "\n\n")
  cat("  ", "Conditions: ", conditions(object), "\n\n")
  cat("  ", "Term Types: ", paste(sapply(object$evterm$events, function(ev) class(ev)[[1]])))
}




