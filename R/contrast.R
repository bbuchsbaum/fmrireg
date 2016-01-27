
#' @export
contrast <- function(A, B=NULL, where=TRUE, split_by=NULL) {
  ret <- list(A=substitute(A),
              B=substitute(B),
              where=substitute(where),
              split_by=substitute(split_by)
  )
  
  class(ret) <- c("contrast_spec", "list")
  ret
}

poly_contrast <- function(A, where=TRUE, degree=1, value_map=NULL) {
  ret <- list(
    A=substitute(A),
    where=substitute(where),
    degree=degree,
    value_map=value_map)
  
  class(ret) <- c("poly_contrast_spec", "list")
  ret
}


contrast_weights.poly_contrast_spec <- function(x, term) {
 
  term.cells <- cells(term)
  row.names(term.cells) <- longnames(term)
 
  keep <- eval(x$where, envir=term.cells, enclos=parent.frame())	
  reduced.term.cells <- subset(term.cells, keep)
  
  keepA <- eval(x$A, envir=term.cells, enclos=parent.frame())
  
  vals <- if (is.null(x$value_map)) {
    as.numeric(as.character(reduced.term.cells[keepA,]))
  } else {
    unlist(x$value_map[as.character(reduced.term.cells[keepA,])])
  }
  

  
  weights <- matrix(0, NROW(term.cells), x$degree)	
  pvals <- stats::poly(vals, degree=x$degree)
  row.names(weights) <- row.names(term.cells)
  colnames(weights) <- paste("poly", 1:x$degree, sep="")
  
  weights[keep, ] <- pvals
  weights
  
}



# 
# .computeWeights <- function(cmod, pat) {
#   patmat <- patternMatch(cmod, pat)
#   if (sum(patmat) <= 0) {
#     stop(paste("failed to match any conditions with pattern: ", pat))
#   }
#   
#   cvec <- ifelse(patmat, 1, 0)
#   cvec/sum(cvec)	
# }


makeWeights <- function(keepA, keepB=NULL) {
  weights <- matrix(0, length(keepA), 1)
  numA <- sum(keepA)
  
  weights[keepA,1] <- rep(1/numA, numA)
  if (!is.null(keepB)) {
    numB <- sum(keepB)
    weights[keepB,1] <- -rep(1/numB, numB)
  }
  
  weights
}

contrast_weights.contrast_spec <- function(x, term) {
  term.cells <- cells(term)
  row.names(term.cells) <- longnames(term)
  
  count <- attr(term.cells, "count")		
  term.cells <- subset(term.cells, count > 0)
  
  keep <- eval(x$where, envir=term.cells, enclos=parent.frame())	
  keepA <- eval(x$A, envir=term.cells, enclos=parent.frame())
  keepB <- if (is.null(x$B)) NULL else eval(x$B, envir=term.cells, enclos=parent.frame())
  split_fac <- eval(x$split_by, envir=term.cells, enclos=parent.frame())
  
  weights <- if (!is.null(split.fac)) {
    split.levs <- levels(x$split.fac)
    do.call(cbind,lapply(split.levs, function(lev) {					
      keepC <- split.fac == lev 
      if (is.null(keepB)) {
        makeWeights(keep & keepA & keepC)
      } else {
        makeWeights(keep & keepA & keepC, keep & keepB & keepC)
      }
    }))
  } else { 
    if (is.null(keepB)) {
      makeWeights(keep & keepA)
    } else {
      makeWeights(keep & keepA, keep & keepB) 
    }
  }
    
  row.names(weights) <- row.names(term.cells)
  weights
}



