
#' contrast_set
#' 
#' @export
#' @import assertthat
contrast_set <- function(...) {
  ret <- list(...)
  assertthat::assert_that(all(sapply(ret, inherits, "contrast_spec")))
  class(ret) <- c("contrast_set", "list")
  ret
}


#' pairwise_contrasts
#' @param levels
#' @param where the subset over which the contrast is computed
#' @export
pairwise_contrasts <- function(levels, where=TRUE) {
  cbns <- combn(length(levels),2)
  ret <- lapply(1:ncol(cbns), function(i) {
    lev1 <- levels[cbns[1,i]]
    lev2 <- levels[cbns[2,i]]
    pair_contrast(as.formula(paste("~", lev1)), as.formula(paste("~", lev2)), where=where, name=paste0("con_", lev1, "_", lev2))
  })
  
  do.call(contrast_set, ret)
}
  


#' pair_contrast
#' 
#' @param A the first \code{formula} expression in the contrast
#' @param B the second \code{formula} expression in the contrast
#' @param name the name of the contrast
#' @param where the subset over which the contrast is computed
#' @param split_by split the contrast by a second factor variable
#' @export
pair_contrast <- function(A, B, name, where=TRUE, split_by=NULL) {
  assert_that(lazyeval::is_formula(A))
  assert_that(lazyeval::is_formula(B))
  ret <- list(A=A,
              B=B,
              where=where,
              split_by=split_by,
              name=name)
    
  class(ret) <- c("contrast_spec", "list")
  ret
}


#' unit_contrast
#' 
#' @param A the contrast expression as a formula
#' @param name the name of the contrast
#' @param where the subset of conditions to apply contrast to
#' @export
unit_contrast <- function(A, name=NULL, where=TRUE, split_by=NULL) {
  assert_that(lazyeval::is_formula(A)) 
  
  if (is.null(name)) {
    name <- as.character(lazyeval::f_rhs(A))
  }
    
  structure(
          list(A=A,
              B=NULL,
              where=lazyeval::expr_find(where),
              split_by=substitute(split_by),
              name=name),
          class=c("unit_contrast_spec", "contrast_spec", "list")
  )
  
}

#' @export
contrast_weights.unit_contrast_spec <- function(x, term) {
  term.cells <- cells(term)
  #row.names(term.cells) <- longnames(term)
  
  keep <- eval(x$where, envir=term.cells, enclos=parent.frame())	
  reduced.term.cells <- subset(term.cells, keep)
  
  mm <- model.matrix(~ reduced.term.cells[[1]] -1)
  colnames(mm) <- levels(reduced.term.cells[[1]])
  
  cvec <- rep(1, nrow(reduced.term.cells))/nrow(reduced.term.cells)
  #cvec <- lazyeval::f_eval_rhs(x$A, data=as.data.frame(mm))
   
  weights <- matrix(0, NROW(term.cells), 1)	
  row.names(weights) <- longnames(term)
  
  weights[keep, ] <- cvec
  
  ret <- list(
    name=x$name,
    weights=weights,
    contrast_spec=x)
  
  class(ret) <- c("unit_contrast", "contrast", "list")
  ret
}


#' @export
poly_contrast <- function(A, name, where=TRUE, degree=1, value_map=NULL) {
  assert_that(lazyeval::is_formula(A))
  ret <- list(
    A=A,
    B=NULL,
    where=substitute(where),
    degree=degree,
    value_map=value_map,
    name=name)
  
  class(ret) <- c("poly_contrast_spec", "contrast_spec", "list")
  ret
}



#' @export
contrast_weights.poly_contrast_spec <- function(x, term) {
 
  term.cells <- cells(term)
  row.names(term.cells) <- longnames(term)
 
  keep <- eval(x$where, envir=term.cells, enclos=parent.frame())	
  reduced.term.cells <- subset(term.cells, keep)
  
 
  vals <- lazyeval::f_eval_rhs(x$A, data=reduced.term.cells)
  
  vals <- if (is.null(x$value_map)) {
    as.numeric(as.character(vals))
  } else {
    unlist(x$value_map[as.character(vals)])
  }
 
  weights <- matrix(0, NROW(term.cells), x$degree)	
  pvals <- stats::poly(vals, degree=x$degree)
  row.names(weights) <- row.names(term.cells)
  colnames(weights) <- paste("poly", 1:x$degree, sep="")
  
  weights[keep, ] <- pvals
  
  ret <- list(
    name=x$name,
    weights=weights,
    condnames=longnames(term),
    contrast_spec=x)
  
  class(ret) <- c("poly_contrast", "contrast", "list")
  ret
  
}


#' @keywords internal
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

#' @export
contrast_weights.contrast_spec <- function(x, term) {
  
  term.cells <- cells(term)
 
  count <- attr(term.cells, "count")		
  term.cells <- subset(term.cells, count > 0)
  keep <- eval(x$where, envir=term.cells, enclos=parent.frame())	
  keepA <- lazyeval::f_eval_rhs(x$A, data=term.cells)
  keepB <- if (is.null(x$B)) NULL else lazyeval::f_eval_rhs(x$B, data=term.cells)
  
  split_fac <- if (!is.null(x$split_fac)) {
    eval(x$split_by, envir=term.cells, enclos=parent.frame())
  }
  
  weights <- if (!is.null(split_fac)) {
    split.levs <- levels(x$split_fac)
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
    
  row.names(weights) <- longnames(term)
  
  ret <- list(
    name=x$name,
    weights=weights,
    condnames=longnames(term),
    contrast_spec=x)
  
  class(ret) <- c("contrast", "list")
  ret  
}

#' @importFrom gmodels estimable
#' @export
estcon.contrast <- function(x, fit, indices) {
  wts <- numeric(length(fit$assign))
  wts[indices] <- x$weights
  
  gmodels::estimable(fit, wts)
}

#' @export
print.contrast_set <- function(x) {
  for (con in x) {
    print(con)
    cat("\n")
  }
}

#' @export
print.contrast_spec <- function(x) {
  cat("contrast:", x$name, "\n")
  cat(" A: ", Reduce(paste, deparse(x$A)), "\n")
  if (!is.null(x$B))
    cat(" B: ", Reduce(paste, deparse(x$B)), "\n")
  if (x$where[1] != TRUE && length(x$where) > 1)
    cat(" where: ", x$where, "\n")
  

}

#' @export
print.contrast <- function(x) {
  print(x$contrast_spec)
  cat(" weights: ", x$weights, "\n")
  cat(" conditions: ", row.names(x$weights))
}

#' @export
print.poly_contrast_spec <- function(x) {
  cat("poly contrast", "\n")
  cat(" A: ", Reduce(paste, deparse(x$A)), "\n")
  cat(" degree: ", x$degree, "\n")
  cat(" values: ", unlist(x$value_map), "\n")
}


