.sanitizeName <- function(name) {
  name <- gsub(":", ".", name)
  name <- gsub(" ", "", name)
  name <- gsub("[\\(\\)]", ".", name, perl=TRUE)
  name <- gsub(",", ".", name)
  name <- gsub("\\.$", "", name)
  name
}

.checkEVArgs <- function(name, vals, onsets, blockids, durations=NULL) {
  stopifnot(length(onsets) == length(vals))
  
  if (is.null(durations) || length(durations) == 1) {
    durations <- rep(durations, length(onsets))
  }
  
  stopifnot(length(durations) == length(vals))
  stopifnot(length(blockids) == length(vals))
  
  list(varname=name, value=vals, onsets=onsets, durations=durations, blockids=blockids)
}

#' @export
#' @rdname event_term-class
event_term <- function(evlist, onsets, blockids, durations = 1, subset=NULL) {
  vnames <- names(evlist)
  evs <- lapply(1:length(evlist), function(i) EV(evlist[[i]], vnames[i], onsets=onsets, blockids=blockids, durations=durations))
  names(evs) <- sapply(evs, function(ev) ev$varname)
  
  if (is.null(subset)) {
    subset=rep(TRUE, length(evs[[1]]$onsets))
  }
  
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
    EventVariable(vals, name, onsets, blockids, durations)
  } else if (is.matrix(vals)) {
    event_set(vals, name, onsets, blockids, durations)
  } else if (is.factor(vals)) {
    EventFactor(vals, name, onsets, blockids, durations)
  } else {
    stop(paste("cannot create event_seq from type: ", typeof(vals)))
  }
  
}

#' @description Create an categorical event sequence from a \code{factor} 
#' @export
EventFactor <- function(fac, name, onsets, blockids=1, durations=NULL) {
  if (!is.factor(fac)) {
    warning("argument 'fac' is not a factor, converting to factor")
    fac <- as.factor(factor())
  }
  
  ret <- .checkEVArgs(name, fac, onsets, blockids, durations)
  ret$continuous = FALSE
  class(ret) <- c("EventFactor", "event_seq")
  ret
}        

#' @description Create a continuous valued event sequence from a \code{numeric} vector.
#' @export
EventVariable <- function(vec, name, onsets, blockids=1, durations=NULL) {
  stopifnot(is.vector(vec))
  
  if (is.factor(vec)) {
    stop("cannot create an EventVariable from a factor, use EventFactor.")
  }
  
  ret <- .checkEVArgs(name, vec, onsets, blockids, durations)
  ret$continuous <- TRUE
  class(ret) <- c("EventVariable", "event_seq")
  ret
  
}       

#' @description Create a continuous valued event set from a \code{matrix}
#' @export
event_set <- function(mat, name, onsets, durations=NULL, blockids=1 ) {
  stopifnot(is.matrix(mat))
  
  ret <- .checkEVArgs(name, mat[,1], onsets, blockids, durations)
  ret$continuous <- TRUE
  
  if (is.null(colnames(mat))) {
    colnames(mat) <- 1:NCOL(mat)
  }
  
  class(ret) <- c("event_set", "event_seq")
  ret
}

#' @description Create a event set from a basis object of type \code{\linkS4Class{ParametricBasis}}. 
#' @export
event_basis <- function(basis, onsets, blockids=1, durations=NULL) {
  stopifnot(inherits(basis, "ParametricBasis"))
  ret <- .checkEVArgs(basis$name, basis$x, onsets, blockids, durations)
  ret$continuous <- TRUE
  ret$basis <- basis
  class(ret) <- c("event_basis", "event_seq")
  ret
}



levels.EventFactor <- function(x) levels(x$value) 

levels.EventVariable <- function(x) x$varname 

levels.event_set <- function(x) colnames(x$value) 

levels.event_basis <- function(x) seq(1, ncol(x$basis$y))

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
cells.EventFactor <- function(x, drop.empty=TRUE) {
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
    attr(evset, "count") <- counts[counts > 0]
  }
  evset
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

colnames.event_term <- function(x) as.vector(unlist(lapply(x$events, colnames)))
colnames.event_seq <- function(x) x$varname
colnames.event_set <- function(x) paste0(.sanitizeName(x$varname), ".", levels(x))
colnames.event_basis <- function(x) paste0(.sanitizeName(x$varname), ".", levels(x))

parentTerms.event_term <- function(x) unlist(lapply(x$events, function(ev) ev$varname))

isContinuous.event_seq <- function(x) x$continuous



elements.event_set <- function(x, values=TRUE) {
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

elements.event_term <- function(x, values=TRUE) {
  els <- lapply(x$events, elements, values=values)
  n <- sapply(names(els), function(nam) .sanitizeName(nam))
  names(els) <- as.vector(n)
  els
}

convolveDesign <- function(hrf, dmat, globons, durations) {
  cond.names <- names(dmat)
  keep <- if (any(is.na(dmat)) || any(is.na(globons))) {
    ret <- apply(dmat, 1, function(vals) all(!is.na(vals)))
    ret[is.na(globons)] <- FALSE
  } else {
    TRUE
  }
  
  
  lapply(1:NCOL(dmat), function(i) {
    Regressor(globons[keep], hrf, amplitude=dmat[keep,i], duration=durations[keep])
  })
  
}

convolve.event_term <- function(x, hrf, samplingFrame, drop.empty=TRUE) {
  globons <- globalOnsets(samplingFrame, x$onsets)
  
  nimages <- sum(samplingFrame$blocklens)
  samples <- seq(samplingFrame$startTime, length.out=nimages, by=samplingFrame$TR)
  
  dmat <- as.data.frame(designMatrix(x, drop.empty))
  
  blockids <- factor(x$blockids)
  split.dmat <- split(dmat, blockids)
  split.ons <- split(globons, blockids)
  split.durations <- split(x$durations, blockids)
  split.samples <- split(samples, rep(1:length(samplingFrame$blocklens), samplingFrame$blocklens))
  
  reglist <- lapply(1:length(split.dmat), function(i) {
    convolveDesign(hrf, split.dmat[[i]], split.ons[[i]], split.durations[[i]])
  })
  
}

designMatrix.event_term <- function(x, drop.empty=TRUE) {
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
  #if (any(counts == 0) && (length(conditions(x, drop=F)) == length(counts)) && drop.unused.levels) {
  #  rmat <- rmat[, !(counts==0), drop=FALSE]
  #  colnames(rmat) <- conditions(x, drop=T)
  #} else {
  #  colnames(rmat) <- conditions(x, drop=F)			
  #}
  
  rmat
}


print.event_term <- function(object) {
  cat("Event Term", "\n")
  cat(" ", "Term Name: ", object$varname, "\n")
  cat(" ", "Formula:  ", as.character(formula(object)), "\n")
  cat(" ", "Num Events: ", nrow(object$eventTable), "\n")
  cat(" ", "Term Types: ", paste(sapply(object$events, function(ev) class(ev)[[1]])))
}
