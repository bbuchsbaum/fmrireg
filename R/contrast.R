## TODO

## when contrasts are misspecified or have spelling errors, errors need to be infromative.


#' Contrast Specification
#'
#' @description
#' Define a linear contrast using a formula expression.
#'
#' @param form A formula describing the contrast.
#' @param name A character label for the contrast.
#' @param where An expression defining the subset over which the contrast is applied (default: NULL).
#'
#' @return A list containing the contrast specification.
#'
#' @examples
#' # A minus B contrast
#' contrast(~ A - B, name="A_B")
#'
#' @export
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

#' Create a Set of Contrasts
#'
#' @description
#' Construct a list of contrast_spec objects.
#'
#' @param ... A variable-length list of contrast_spec objects.
#'
#' @return A list of contrast_spec objects with class "contrast_set".
#'
#' @examples
#' c1 <- contrast(~ A - B, name="A_B")
#' c2 <- contrast(~ B - C, name="B_C")
#' contrast_set(c1,c2)
#'
#' @export
#' @import assertthat
#' @importFrom purrr map_lgl
contrast_set <- function(...) {
  ret <- list(...)
  assertthat::assert_that(all(map_lgl(ret, inherits, "contrast_spec")))
  class(ret) <- c("contrast_set", "list")
  ret
}

#' Pairwise Contrasts
#'
#' @description
#' Construct pairwise contrasts for all combinations of factor levels.
#'
#' @param levels A vector of factor levels to be compared.
#' @param where An optional formula specifying the subset over which the contrast is computed.
#'
#' @return A contrast_set object containing pairwise contrasts for all combinations of factor levels.
#'
#' @examples
#' pairwise_contrasts(c("A", "B", "C"))
#'
#' @export
#' @importFrom utils combn
pairwise_contrasts <- function(levels, where=NULL) {
  if (!is.null(where)) {
    assert_that(lazyeval::is_formula(where))
  }
  
  cbns <- combn(length(levels),2)
  ret <- lapply(1:ncol(cbns), function(i) {
    lev1 <- levels[cbns[1,i]]
    lev2 <- levels[cbns[2,i]]
    pair_contrast(as.formula(paste("~", lev1)), as.formula(paste("~", lev2)), where=where, name=paste0("con_", lev1, "_", lev2))
  })
  
  do.call(contrast_set, ret)
}
  
#' One Against All Contrast
#'
#' @description
#' Construct contrasts comparing each factor level against the average of the other levels.
#'
#' @param levels A vector of factor levels to be compared.
#' @param facname A character string specifying the name of the factor containing the supplied levels.
#' @param where An optional formula specifying the subset over which the contrast is computed.
#'
#' @return A contrast_set object containing contrasts comparing each factor level against the average of the other levels.
#'
#' @examples
#' fac <- factor(rep(c("A", "B", "C"), 2))
#' con <- one_against_all_contrast(levels(fac), "fac")
#'
#' @export
one_against_all_contrast <- function(levels, facname, where=NULL) {
  if (!is.null(where)) {
    assert_that(lazyeval::is_formula(where))
  }
  
  ret <- lapply(1:length(levels), function(i) {
    lev1 <- levels[i]
    levother <- levels[-i]
    pair_contrast(as.formula(paste("~", facname, " == ", paste0('"', lev1, '"'))), 
                             as.formula(paste0("~", facname, "!= ", paste0('"', lev1, '"'))), 
                  where=where, name=paste0("con_", lev1, "_vs_", "other"))
  })
  
  do.call(contrast_set, ret)
  
}

#' Pair Contrast
#'
#' @description
#' Construct a sum-to-zero contrast between two logical expressions.
#'
#' @param A A formula representing the first logical expression in the contrast.
#' @param B A formula representing the second logical expression in the contrast.
#' @param name A character string specifying the name of the contrast (mandatory).
#' @param where An optional formula specifying the subset over which the contrast is computed.
#'
#' @return A pair_contrast_spec object containing the sum-to-zero contrast between the two logical expressions.
#'
#' @examples
#' # A hypothetical experiment with a factor 'category' that takes on values of 'face' and 'scene'
#' pair_contrast(~ category == "face", ~ category == "scene", name="face_scene")
#'
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


#' Unit Contrast
#'
#' @description
#' Construct a contrast that sums to 1 and is used to define contrasts against the baseline.
#'
#' @param A A formula representing the contrast expression.
#' @param name A character string specifying the name of the contrast.
#' @param where An optional formula specifying the subset of conditions to apply the contrast to.
#'
#' @return A unit_contrast_spec object containing the contrast that sums to 1.
#'
#' @examples
#' con <- unit_contrast(~ Face, name="Main_face")
#'
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

#' Unit Contrast Weights
#'
#' @description
#' Compute the contrast weights for a unit_contrast_spec object.
#'
#' @param x A unit_contrast_spec object.
#' @param term A term object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing the term, name, weights, condition names, and contrast specification.
#'
#' @export
contrast_weights.unit_contrast_spec <- function(x, term,...) {
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
    term=term,
    name=x$name,
    weights=weights,
    condnames=longnames(term),
    contrast_spec=x)
  
  class(ret) <- c("unit_contrast", "contrast", "list")
  ret
}

#' Polynomial Contrast
#'
#' @description
#' Create a polynomial contrast evaluated over a set of (orderable) factor levels.
#'
#' @param A A formula describing the contrast.
#' @param name A character string representing the name of the contrast.
#' @param where An optional formula defining the subset over which the contrast is computed.
#' @param degree An integer representing the degree of the polynomial (default: 1).
#' @param value_map An optional list that maps between levels of a factor and a numeric value.
#'
#' @return A list containing the polynomial contrast specification.
#'
#' @examples
#' pcon <- poly_contrast(~ time, name="poly_time_3", degree=3)
#'
#' @export
poly_contrast <- function(A, name, where=NULL, degree=1, value_map=NULL) {
  assert_that(lazyeval::is_formula(A))
  
  if (!is.null(where)) {
    assert_that(lazyeval::is_formula(where))
  }

  ret <- list(
    A=A,
    B=NULL,
    where=where,
    degree=degree,
    value_map=value_map,
    name=name)
  
  class(ret) <- c("poly_contrast_spec", "contrast_spec", "list")
  ret
}



#' Polynomial Contrast Weights
#'
#' @description
#' Compute the contrast weights for a polynomial contrast specification.
#'
#' @param x A poly_contrast_spec object.
#' @param term A term object.
#' @param ... Additional arguments (currently not used).
#'
#' @return A list containing the polynomial contrast weights and related information.
#'
#' @export
contrast_weights.poly_contrast_spec <- function(x, term,...) {
  term.cells <- cells(term)
  term.cells <- term.cells %>% dplyr::mutate(rowname=longnames(term))
  #row.names(term.cells) <- longnames(term)
 
  #keep <- eval(x$where, envir=term.cells, enclos=parent.frame())	
  if (!is.null(x$where)) {
    keep <- lazyeval::f_eval_rhs(x$where, data=term.cells)
    assert_that(sum(keep) > 0)
    #term.cells <- term.cells[keep,]
  } else {
    keep <- rep(TRUE, nrow(term.cells))
  }
  
  
  reduced.term.cells <- subset(term.cells, keep)
  
  vals <- lazyeval::f_eval_rhs(x$A, data=reduced.term.cells)
  
  vals <- if (is.null(x$value_map)) {
    as.numeric(as.character(vals))
  } else {
    unlist(x$value_map[as.character(vals)])
  }
 
  weights <- matrix(0, NROW(term.cells), x$degree)	
  pvals <- stats::poly(vals, degree=x$degree)
  #row.names(weights) <- row.names(term.cells)
  row.names(weights) <- term.cells$rowname
  colnames(weights) <- paste("poly", 1:x$degree, sep="")
  
  weights[keep, ] <- pvals
  
  ret <- list(
    term=term,
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
`-.contrast_spec` <- function(e1, e2, ...){
  assert_that(inherits(e2, "contrast_spec"))
  structure(list(
    name=paste0(e1$name, ":", e2$name),
    con1=e1,
    con2=e2),
    class=c("contrast_diff_spec", "contrast_spec", "list")
  )
}

#' @export
contrast_weights.contrast_diff_spec <- function(x, term,...) {
  wts1 <- contrast_weights(x$con1, term)
  wts2 <- contrast_weights(x$con2, term)

  ret <- structure(
    list(
      term=term,
      name=x$name,
      weights=wts1$weights - wts2$weights,
      condnames=longnames(term),
      contrast_spec=x),
    class=c("contrast_diff", "contrast")
  )

  ret  
}

#' @export
contrast_weights.pair_contrast_spec <- function(x, term,...) {
  term.cells <- cells(term)
  
  count <- attr(term.cells, "count")		
  term.cells <- subset(term.cells, count > 0)
  
  keep <- if (!is.null(x$where)) {
    lazyeval::f_eval_rhs(x$where, data=term.cells)
    #eval(x$where, envir=term.cells, enclos=parent.frame())	
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
    term=term,
    name=x$name,
    weights=weights,
    condnames=longnames(term),
    contrast_spec=x)
  
  class(ret) <- c("contrast", "list")
  ret  
}

#' @export
contrast_weights.contrast_formula_spec <- function(x, term,...) {

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

  frm <- paste0(gsub(":", ".", deparse(lazyeval::f_rhs(x$A))), collapse="")
  A <- as.formula(paste("~", frm))
  
  if (is_continuous.event_term(term)) {
    ## hack to handle contrasts with continuous terms
    facs <- !sapply(term$evterm$events, is_continuous)
    term.cells[,!facs] <- 1
  } 
  
  modmat <- suppressMessages(tibble::as_tibble(model.matrix(cform,data=term.cells), 
                                               .name_repair="check_unique"))
  names(modmat) <- gsub(":", ".", condnames)
  weights <- matrix(0, NROW(term.cells), 1)
  weights[keep,1] <-  as.vector(lazyeval::f_eval_rhs(A, data=modmat))
  row.names(weights) <- row.names(term.cells)
  
  ret <- list(
    term=term,
    name=x$name,
    weights=weights,
    condnames=longnames(term),
    contrast_spec=x)
  
  
  class(ret) <- c("contrast", "list")
  ret  
}

#' convert a contrast to an AFNI 'GLT' 
#' 
#' @param x the contrast to convert
#' @param ... extra args
#' @keywords internal
to_glt <- function(x, ...) UseMethod("to_glt")


#' @export
to_glt.contrast <- function(x,...) {
  if (is.matrix(x$weights) && ncol(x$weights) > 1) {
    glts <- lapply(1:ncol(x$weights), function(i) {
      paste0(signif(x$weights[,i],4), "*", x$condnames, collapse=" ")
    })
    
    ret <- list(glt_str=glts,
                name=paste0("GLT_", x$name, "_", 1:ncol(x$weights)),
                con=x)
    
    class(ret) <- "glt_contrast_list"
    ret
  } else {
    glt <- paste0(signif(x$weights,4), "*", x$condnames, collapse=" ")
    ret <- list(glt_str=glt,
       name=paste0("GLT_", x$name),
       con=x)
  
    class(ret) <- "glt_contrast"
    ret
  }
}


#' write a GLT file to disk
#' 
#' @param x the object
#' @param fname the file name to write to
#' @keywords internal
write_glt <- function(x, fname) UseMethod("write_glt")

#' @export
write_glt.glt_contrast <- function(x, fname=NULL) {
  con <- if (is.null(fname)) {
    file(fname, "w")
  } else {
    file(paste0(x$name, ".txt"), "w")
  }
  
  writeLines(x$glt_str, con=con)
  close(con)
}

#' @importFrom gmodels estimable
#' @export
estcon.contrast <- function(x, fit, indices, ...) {
  wts <- numeric(length(fit$assign))
  wts[indices] <- x$weights
  
  gmodels::estimable(fit, wts)
}

#' @export
print.contrast_set <- function(x,...) {
  for (con in x) {
    print(con)
    cat("\n")
  }
}

#' @export
print.contrast_spec <- function(x,...) {
  cat("contrast:", x$name, "\n")
  cat(" A: ", Reduce(paste, deparse(x$A)), "\n")
  if (!is.null(x$B))
    cat(" B: ", Reduce(paste, deparse(x$B)), "\n")
  if (!is.null(x$where))
    cat(" where: ", Reduce(paste, deparse(x$where)), "\n")
  

}

#' @export
print.contrast <- function(x,...) {
  print(x$contrast_spec)
  cat(" term: ", x$term$varname, "\n")
  cat(" weights: ", "\n")
  print(x$weights)
  cat(" conditions: ", x$condnames)
  
}

#' @export
print.poly_contrast_spec <- function(x,...) {
  cat("poly contrast:", "\n")
  cat(" A: ", Reduce(paste, deparse(x$A)), "\n")
  cat(" degree: ", x$degree, "\n")
  if (!is.null(x$where)) {
    cat(" where: ", deparse(x$where), "\n")
  }
  
  if (!is.null(x$value_map)) {
    cat(" values: ", unlist(x$value_map), "\n")
  }
}

#' @export
print.contrast_diff_spec <- function(x,...) {
  cat("contrast difference:", "\n")
  cat("  ", x$con1$name, "-", x$con2$name, "\n")
}



