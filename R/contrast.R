
#' contrast
#' 
#' @param form
#' @param name
#' @param where
contrast <- function(form, name, where=NULL) {
  assert_that(lazyeval::is_formula(form))
  if (!is.null(where)) {
    assert_that(lazyeval::is_formula(where))
  }
  ret <- list(A=form,
              B=NULL,
              where=where,
              name=name)
  
  class(ret) <- c("contrast_formula_spec", "contrast_spec", "list")
  ret
  
}


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
#' @export
pair_contrast <- function(A, B, name, where=NULL) {
  assert_that(lazyeval::is_formula(A))
  assert_that(lazyeval::is_formula(B))
  
  if (!is.null(where)) {
    assert_that(lazyeval::is_formula(where))
  }
  
  ret <- list(A=A,
              B=B,
              where=where,
              name=name)
    
  class(ret) <- c("pair_contrast_spec", "contrast_spec", "list")
  ret
}


#' unit_contrast
#' 
#' @param A the contrast expression as a formula
#' @param name the name of the contrast
#' @param where the subset of conditions to apply contrast to
#' @export
unit_contrast <- function(A, name, where=NULL) {
  assert_that(lazyeval::is_formula(A)) 
  
  if (!is.null(where)) {
    assert_that(lazyeval::is_formula(where))
  }
    
  structure(
          list(A=A,
               B=NULL,
               where=where,
               name=name),
          class=c("unit_contrast_spec", "contrast_spec", "list")
  )
  
}

#' @export
contrast_weights.unit_contrast_spec <- function(x, term) {
  term.cells <- cells(term)
  
  if (!is.null(x$where)) {
    keep <- lazyeval::f_eval_rhs(x$where, data=term.cells)	
    reduced.term.cells <- subset(term.cells, keep)
  } else {
    keep <- rep(TRUE, nrow(term.cells))
    reduced.term.cells <- term.cells
  }
  
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

`-.contrast_spec` <- function(e1, e2, ...){
  assert_that(inherits(e2, "contrast_spec"))
  structure(list(
    name=paste0(e1$name, ":", e2$name),
    con1=e1,
    con2=e2),
    class="contrast_diff_spec"
  )
}

contrast_weights.contrast_diff_spec <- function(x, term) {
  wts1 <- contrast_weights(x$con1, term)
  wts2 <- contrast_weights(x$con2, term)

  ret <- structure(
    list(
      name=x$name,
      weights=wts1$weights - wts2$weights,
      condnames=longnames(term),
      contrast_spec=x),
    class="contrast")
  
  class(ret) <- c("contrast", "list")
  ret  
}

#' @export
contrast_weights.pair_contrast_spec <- function(x, term) {
  term.cells <- cells(term)
  
  count <- attr(term.cells, "count")		
  term.cells <- subset(term.cells, count > 0)
  keep <- if (!is.null(x$where)) {
    eval(x$where, envir=term.cells, enclos=parent.frame())	
  } else {
    rep(TRUE, nrow(term.cells))
  }
  
  keepA <- lazyeval::f_eval_rhs(x$A, data=term.cells)
  keepB <- if (is.null(x$B)) NULL else lazyeval::f_eval_rhs(x$B, data=term.cells)

  weights <- if (is.null(keepB)) {
    makeWeights(keep & keepA)
  } else {
    makeWeights(keep & keepA, keep & keepB) 
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

#' @export
contrast_weights.contrast_formula_spec <- function(x, term) {
 
  term.cells <- cells(term)
  cform <- as.formula(paste("~", paste0(names(term.cells), collapse=":"), "-1"))
  condnames <- shortnames(term)
  count <- attr(term.cells, "count")		
  term.cells <- subset(term.cells, count > 0)
  
  if (!is.null(x$where)) {
    keep <- lazyeval::f_eval_rhs(x$where, data=term.cells)
    assert_that(sum(keep) > 0)
    term.cells <- term.cells[keep,]
  } else {
    keep <- rep(TRUE, nrow(term.cells))
  }

   
  A <- as.formula(paste("~", gsub(":", ".", deparse(lazyeval::f_rhs(x$A)))))
  modmat <- tibble::as_tibble(model.matrix(cform,data=term.cells))
  names(modmat) <- gsub(":", ".", condnames)
  
  weights <- matrix(0, NROW(term.cells), 1)
  weights[keep,1] <-  as.vector(lazyeval::f_eval_rhs(A, data=modmat))
  row.names(weights) <- row.names(term.cells)
  
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
  if (!is.null(x$where))
    cat(" where: ", Reduce(paste, deparse(x$where)), "\n")
  

}

#' @export
print.contrast <- function(x) {
  print(x$contrast_spec)
  cat(" weights: ", x$weights, "\n")
  cat(" conditions: ", x$condnames)
}

#' @export
print.poly_contrast_spec <- function(x) {
  cat("poly contrast", "\n")
  cat(" A: ", Reduce(paste, deparse(x$A)), "\n")
  cat(" degree: ", x$degree, "\n")
  cat(" values: ", unlist(x$value_map), "\n")
}


