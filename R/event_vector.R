#' @import checkmate
NULL

.sanitizeName <- function(name) {
  name <- gsub(":", ".", name)
  name <- gsub(" ", "", name)
  name <- gsub("[\\(\\)]", ".", name, perl=TRUE)
  name <- gsub(",", "_", name)
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
#' 
#' @param evlist a list of named variables
#' @param onsets the onset times from the experimental events in seconds
#' @param blockids the block number associated with each onset
#' @param durations
#' @param subset
#' @export
#' @rdname event_term-class
event_term <- function(evlist, onsets, blockids, durations = 1, subset=NULL) {
  assert_that(is.increasing(blockids))
            
  if (is.null(subset)) { subset=rep(TRUE, length(onsets)) }
  
  if (length(durations) == 1) {
    durations <- rep(durations, length(onsets))
  }
  
  vnames <- names(evlist)
  evs <- lapply(1:length(evlist), function(i) EV(evlist[[i]][subset], vnames[i], 
                                                 onsets=onsets[subset], 
                                                 blockids=blockids[subset], 
                                                 durations=durations[subset]))
  
  names(evs) <- sapply(evs, function(ev) ev$varname)
  
  pterms <- unlist(lapply(evs, function(ev) ev$varname))
  
  len <- sum(subset)
  etab <- tibble::as_data_frame(lapply(pterms, function(termname) {
    if (is_continuous(evs[[termname]])) {
      rep(.sanitizeName(termname), len)
    } else {
      evs[[termname]]$value
    }			
  }))
  
  names(etab) <- sapply(pterms, .sanitizeName)
  varname <- paste(sapply(evs, function(x) x$varname), collapse=":")
  ret <- list(varname=varname, events=evs, subset=subset, event_table=etab, 
              onsets=evs[[1]]$onsets, 
              blockids=evs[[1]]$blockids, 
              durations=evs[[1]]$durations)
  class(ret) <- c("event_term", "event_seq")
  ret
}

event_table.event_term <- function(x) x$event_table


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
#' 
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

#' @export
levels.event_matrix <- function(x) colnames(x$value) 

#' @export
levels.event_set <- function(x) colnames(x$value) 

#' @export
levels.event_basis <- function(x) seq(1, ncol(x$basis$y))

#' @export
formula.event_term <- function(x) as.formula(paste("~ ", "(", paste(parent_terms(x), collapse=":"), "-1", ")"))

#' @export
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

#' @export
cells.event_factor <- function(x, drop.empty=TRUE) {
  etab <- data.frame(onsets=x$onsets, durations=x$durations, blockids=x$blockids)
  split(etab, x$value)
}

#' @export
cells.event_term <- function(x, drop.empty=TRUE) {
  evtab <- x$event_table
  evset <- tibble::as_tibble(expand.grid(lapply(x$events, levels)))
  which.cat <- which(!sapply(x$events, is_continuous))
  
  evs <- tibble::as_tibble(lapply(evset[,which.cat], as.character))
  evt <- tibble::as_tibble(lapply(evtab[,which.cat], as.character))
  
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
  evset
}

.event_set <- function(x) {
  evtab <- event_table(x)
  
  evset <- if (nbasis(x) > 1) {
    evlist <- c(list(factor(paste("basis", 1:nbasis(x), sep = ""))), cells(x@eventTerm))
    names(evlist) <- c("basis", parent_terms(x@eventTerm))
    evlist <- lapply(evlist, levels)
    ret <- expand.grid(evlist, stringsAsFactors = TRUE)
    ret[c(2:length(ret), 1)]
  } else {
    cells(x$evterm)
  }
  
}

#' @export
cells.convolved_term <- function(x) {
  evtab <- event_table(x)
  evset <- .event_set(x)
  
  strels <- apply(apply(evtab, 2, str_trim), 1, paste, collapse = ":")
  strlevs <- if (nrow(evset) > 1) {
    apply(apply(evset, 2, str_trim), 1, paste, collapse = ":")
  } else {
    as.character(evset[1,1])
  }
  
  attr(evset, "rownames") <- strlevs
  counts <- rep(attr(cells(x$evterm), "count"), each = nbasis(x))
  
  ret <- evset[counts > 0, , drop = F]
  attr(ret, "count") <- counts[counts > 0]
  ret
  
}

#' @export
conditions.fmri_term <- function(x) {
  colnames(design_matrix(x))
}

#' @export
conditions.convolved_term <- function(x) {
  colnames(x$design_matrix)
}


#' @export
conditions.event_term <- function(x, drop.empty=TRUE) {
  
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
parent_terms.event_term <- function(x) unlist(lapply(x$events, function(ev) ev$varname))

#' @export
is_continuous.event_seq <- function(x) x$continuous

#' @export
is_categorical.event_seq <- function(x) !x$continuous

#' @export
is_continuous.event_term <- function(x) all(sapply(x$events, function(x) is_continuous(x)))

#' @export
is_categorical.event_term <- function(x) !is_continuous(x)

#' @export
is_categorical.event_seq <- function(x) !x$continuous


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
    ret
  } else {
    rep(TRUE, nrow(dmat))
  }
  

  lapply(1:NCOL(dmat), function(i) {
    regressor(globons[keep], hrf, amplitude=unlist(dmat[keep,i]), duration=durations[keep])
  })
  
}

#' @importFrom tibble as_tibble
#' @importFrom dplyr group_by select do ungroup
#' @export
convolve.event_term <- function(x, hrf, sframe, drop.empty=TRUE) {
  globons <- global_onsets(sframe, x$onsets, x$blockids)
  durations <- x$durations
  blockids <- x$blockids
  
  nimages <- sum(sframe$blocklens)
  
  cnames <- conditions(x)
  
  dmat <- design_matrix(x, drop.empty)
  ncond <- ncol(dmat)

  cmat <- dmat %>% dplyr::mutate(.blockids=blockids, .globons=globons, .durations=durations) %>% 
    dplyr::group_by(.blockids) %>%
    dplyr::do({
      d <- dplyr::select(., 1:ncond)
      reg <- convolve_design(hrf, d, .$.globons, .$.durations)
      sam <- samples(sframe, blockids=as.integer(as.character(.$.blockids[1])), global=TRUE)
      ret <- do.call(cbind, lapply(reg, function(r) evaluate(r, sam)))
      tibble::as_tibble(ret)
  }) %>% dplyr::ungroup() %>% dplyr::select(-.blockids)
  
 
  if (nbasis(hrf) > 1) {
    blevs <- paste("[", 1:nbasis(hrf), "]", sep="")
    cnames <- unlist(lapply(cnames, function(prefix) paste(prefix, ":basis", blevs, sep="")))
  } 
            
  colnames(cmat) <- cnames
  tibble::as_tibble(cmat)
  
  #lapply(reglist, function(reg) evaluate(reg, )
  
}


Fcontrasts.event_term <- function(x) {
  which_cat <- which(sapply(x$events, function(obj) is_categorical(obj)))
  assert_that(length(which_cat) > 0)
  pterms <- parent_terms(x)[which_cat]
  evs <- x$events[which_cat]
  cond <- conditions(x)
  facnames <- names(evs)
  
  Clist <- lapply(evs, function(ev) rep(1, length(levels(ev))))
  Dlist <- lapply(evs, function(ev) t(-diff(diag(length(levels(ev))))))
  
  nfac <- length(Clist)

  main_effects <- lapply(length(Clist):1, function(i) {
    Dcon <- Dlist[[i]]
    Cs <- Clist[-i]
    mats <- vector(nfac, mode="list")
    mats[[i]] <- Dcon
    mats[seq(1, nfac)[-i]] <- Cs
    ret <- Reduce(kronecker, mats)
    row.names(ret) <- cond
    ret
  })
  
  names(main_effects) <- facnames
  
  if (length(facnames) > 1) {
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
#' @export
design_matrix.event_term <- function(x, drop.empty=TRUE) {
  locenv <- new.env()
  pterms <- sapply(parent_terms(x), .sanitizeName)	
  
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
  rmat <- mat #* x$subset
  
  #remove rows with NAS
  if (any(nas)) {
    rmat <- matrix(0, nrow(x$event_table), length(conditions(x, drop.empty)))
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
  
  tibble::as_tibble(rmat)
}



#' @export
print.event_term <- function(object) {
  cat("event_term", "\n")
  cat("  ", "Term Name: ", object$varname, "\n")
  cat("  ", "Formula:  ", as.character(formula(object)), "\n")
  cat("  ", "Num Events: ", nrow(object$event_table), "\n")
  cat("  ", "Term Types: ", paste(sapply(object$events, function(ev) class(ev)[[1]])))
}

#' @export
print.fmri_term <- function(object) {
  cat("fmri_term: ", class(object)[[1]], "\n")
  cat("  ", "Term Name: ", object$varname, "\n")
  cat("  ", "Num Rows: ", nrow(design_matrix(object)), "\n")
  cat("  ", "Num Columns: ", ncol(design_matrix(object)), "\n")
}

#' @export
print.convolved_term <- function(object) {
  cat("fmri_term: ", class(object)[[1]], "\n")
  cat("  ", "Term Name: ", object$varname, "\n")
  cat("  ", "Formula:  ", as.character(formula(object$evterm)), "\n")
  cat("  ", "Num Events: ", nrow(object$evterm$event_table), "\n")
  cat("  ", "Num Rows: ", nrow(design_matrix(object)), "\n")
  cat("  ", "Num Columns: ", ncol(design_matrix(object)), "\n")
  cat("  ", "Conditions: ", conditions(object), "\n")
  cat("  ", "Term Types: ", paste(sapply(object$evterm$events, function(ev) class(ev)[[1]])))
}




