#' @keywords internal
#' @noRd
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
#' @param facname The name of the factor variable (column name in the design) these levels belong to.
#' @param where An optional formula specifying the subset over which the contrast is computed.
#' @param name_prefix A character string to prefix the generated contrast names (default: "con").
#'
#' @return A contrast_set object containing pairwise contrasts for all combinations of factor levels.
#'
#' @examples
#' # Assuming 'my_factor' is a column name
#' pairwise_contrasts(c("A", "B", "C"), facname = "my_factor")
#' pairwise_contrasts(c("A", "B", "C"), facname = "my_factor", name_prefix = "pair")
#'
#' @export
#' @importFrom utils combn
pairwise_contrasts <- function(levels, facname, where=NULL, name_prefix = "con") {
  assert_that(is.character(facname), length(facname) == 1, msg = "'facname' must be a single string.")
  if (!is.null(where)) {
    assert_that(lazyeval::is_formula(where))
  }
  
  if (length(levels) < 2) {
    stop("pairwise_contrasts requires at least two levels.")
  }
  
  cbns <- combn(length(levels), 2)
  ret <- lapply(1:ncol(cbns), function(i) {
    lev1 <- levels[cbns[1,i]]
    lev2 <- levels[cbns[2,i]]
    # Construct formulas using the factor name
    formula_A <- as.formula(paste("~", facname, "==", paste0('"', lev1, '"')))
    formula_B <- as.formula(paste("~", facname, "==", paste0('"', lev2, '"')))
    pair_contrast(formula_A, formula_B, 
                  where=where, name=paste0(name_prefix, "_", lev1, "_", lev2))
  })
  
  do.call(contrast_set, ret)
}




#' Pair Contrast
#'
#' @description
#' Construct a sum-to-zero contrast between two logical expressions. This function is
#' particularly useful for comparing specific conditions or combinations of conditions.
#'
#' @param A A formula representing the first logical expression in the contrast.
#' @param B A formula representing the second logical expression in the contrast.
#' @param name A character string specifying the name of the contrast (mandatory).
#' @param where An optional formula specifying the subset over which the contrast is computed.
#'
#' @return A pair_contrast_spec object containing:
#'   \item{A}{First logical expression}
#'   \item{B}{Second logical expression}
#'   \item{where}{Subsetting formula (if provided)}
#'   \item{name}{Contrast name}
#'
#' @details
#' The contrast is constructed as (A - B), where A and B are logical expressions that
#' evaluate to TRUE/FALSE for each observation. The resulting contrast weights sum to zero.
#'
#' @examples
#' # Compare faces vs scenes
#' pair_contrast(~ category == "face", ~ category == "scene", name = "face_vs_scene")
#'
#' # Compare with subsetting
#' pair_contrast(~ category == "face", ~ category == "scene",
#'              name = "face_vs_scene_block1",
#'              where = ~ block == 1)
#'
#' @seealso
#' \code{\link{pairwise_contrasts}} for all pairwise comparisons,
#' \code{\link{contrast_set}} for creating sets of contrasts
#'
#' @export
pair_contrast <- function(A, B, name, where = NULL) {
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


#' One-way Contrast
#'
#' @description
#' Create a one-way contrast specification
#'
#' @param A A formula specifying the contrast
#' @param name The name of the contrast
#' @param where Optional environment for evaluating the formula
#' @return A oneway_contrast_spec object that can be used to generate contrast weights
#' @examples
#' # Create a one-way contrast for a factor 'basis'
#' con <- oneway_contrast(~ basis, name = "Main_basis")
#'
#' # Create a one-way contrast with a specific environment
#' con <- oneway_contrast(~ basis, name = "Main_basis",
#'                       where = new.env())
#'
#' @seealso \code{\link{interaction_contrast}} for testing interactions,
#'          \code{\link{pair_contrast}} for pairwise comparisons
#' @export
oneway_contrast <- function(A, name, where = NULL) {
  assert_that(lazyeval::is_formula(A)) 
  
  if (!is.null(where)) {
    assert_that(lazyeval::is_formula(where))
  }
  
  structure(
    list(A=A,
         B=NULL,
         where=where,
         name=name),
    class=c("oneway_contrast_spec", "contrast_spec", "list")
  )
}

#' Interaction Contrast
#'
#' @description
#' Create an interaction contrast specification
#'
#' @param A A formula specifying the interaction contrast
#' @param name The name of the contrast
#' @param where Optional environment for evaluating the formula
#' @return An interaction_contrast_spec object containing the specification for
#'         generating interaction contrast weights
#' @examples
#' # Create an interaction contrast for factors A and B
#' con <- interaction_contrast(~ A * B, name = "A_by_B")
#'
#' # Create an interaction contrast with a specific environment
#' con <- interaction_contrast(~ A * B, name = "A_by_B",
#'                           where = new.env())
#'
#' @seealso \code{\link{oneway_contrast}} for main effects,
#'          \code{\link{pair_contrast}} for pairwise comparisons
#' @export
interaction_contrast <- function(A, name, where = NULL) {
  assert_that(lazyeval::is_formula(A)) 
  
  if (!is.null(where)) {
    assert_that(lazyeval::is_formula(where))
  }
  
  
  
  structure(
    list(A=A,
         B=NULL,
         where=where,
         name=name),
    class=c("interaction_contrast_spec", "contrast_spec", "list")
  )
}


#' Polynomial Contrast
#'
#' @description
#' Create polynomial contrasts for testing trends across ordered factor levels. This is
#' particularly useful for analyzing factors with a natural ordering (e.g., time, dose).
#'
#' @param A A formula specifying the ordered factor.
#' @param name A character string identifying the contrast.
#' @param where An optional formula for subsetting the data.
#' @param degree An integer specifying the degree of the polynomial (default: 1).
#' @param value_map An optional list mapping factor levels to numeric values.
#'
#' @return A poly_contrast_spec object containing the specification for generating
#'   polynomial contrast weights.
#'
#' @details
#' The function creates orthogonal polynomial contrasts up to the specified degree.
#' These contrasts can test for linear, quadratic, cubic, and higher-order trends
#' in the data. The value_map parameter allows for non-uniform spacing between levels.
#'
#' @examples
#' # Linear trend across time points
#' pcon <- poly_contrast(~ time, name = "linear_time", degree = 1)
#'
#' # Cubic trend with custom spacing
#' pcon <- poly_contrast(~ dose, name = "dose_cubic",
#'                      degree = 3,
#'                      value_map = list("low" = 0, "med" = 2, "high" = 5))
#'
#' @seealso
#' \code{\link{oneway_contrast}} for categorical contrasts,
#' \code{\link{interaction_contrast}} for interaction effects
#'
#' @export
poly_contrast <- function(A, name, where = NULL, degree = 1, value_map = NULL) {
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


#' One-way Contrast Weights
#'
#' @description
#' Compute the contrast weights for an oneway_contrast_spec object.
#'
#' @param x An oneway_contrast_spec object.
#' @param term A term object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing the term, name, weights, condition names, and contrast specification.
#'
#' @export
contrast_weights.oneway_contrast_spec <- function(x, term,...) {
  term.cells <- cells(term)
  count <- attr(cells(term, drop.empty=FALSE), "count")
  
  if (!is.null(x$where)) {
    keep <- lazyeval::f_eval_rhs(x$where, data=term.cells)	
    reduced.term.cells <- subset(term.cells, keep)
  } else {
    keep <- rep(TRUE, nrow(term.cells))
    reduced.term.cells <- term.cells
  }
  

  
  fac <- all.vars(rlang::f_rhs(x$A))
  if (length(fac) > 1) {
    warning("one-way contrast has more than one factor specified, taking first one only.")
  }
  
  cmat <- generate_main_effect_contrast(reduced.term.cells, fac[1])
  #fac <- as.factor(lazyeval::f_eval_rhs(x$A, reduced.term.cells))
  #contrasts(fac) <- contr.sum(levels(fac))
  #excon <- model.matrix(~ fac - 1)
  #excon <- contrasts(fac)[as.numeric(fac), ]
  
  weights <- matrix(0, NROW(term.cells), ncol(cmat))	
  row.names(weights) <- longnames(term)
  
  weights[keep, ] <- cmat
  
  ret <- list(
    term=term,
    name=x$name,
    weights=weights,
    condnames=longnames(term),
    contrast_spec=x)
  
  class(ret) <- c("oneway_contrast", "Fcontrast", "contrast", "list")
  ret
}

#' Interaction Contrast Weights
#'
#' @description
#' Compute the contrast weights for an interaction_contrast_spec object.
#'
#' @param x An interaction_contrast_spec object.
#' @param term A term object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing the term, name, weights, condition names, and contrast specification.
#'
#' @export
contrast_weights.interaction_contrast_spec <- function(x, term,...) {
  term.cells <- cells(term)
  if (!is.null(x$where)) {
    keep <- lazyeval::f_eval_rhs(x$where, data=term.cells)	
    reduced.term.cells <- subset(term.cells, keep)
  } else {
    keep <- rep(TRUE, nrow(term.cells))
    reduced.term.cells <- term.cells
  }
  
  factors <- all.vars(rlang::f_rhs(x$A))
  cmat <- generate_interaction_contrast(reduced.term.cells, factors)
  
  weights <- matrix(0, NROW(term.cells), ncol(cmat))	
  row.names(weights) <- longnames(term)
  
  weights[keep, ] <- cmat
  
  ret <- list(
    term=term,
    name=x$name,
    weights=weights,
    condnames=longnames(term),
    contrast_spec=x)
  
  class(ret) <- c("interaction_contrast", "Fcontrast",  "contrast", "list")
  ret
  
  
}

#' Polynomial Contrast Weights
#'
#' @description
#' Compute the contrast weights for a poly_contrast_spec object.
#'
#' @param x A poly_contrast_spec object.
#' @param term A term object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing the term, name, weights, condition names, and contrast specification.
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

#' Pair Contrast Weights
#'
#' @description
#' Compute the contrast weights for a pair_contrast_spec object.
#'
#' @param x A pair_contrast_spec object.
#' @param term A term object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing the term, name, weights, condition names, and contrast specification.
#'
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

#' Contrast Formula Weights
#'
#' @description
#' Compute the contrast weights for a contrast_formula_spec object.
#'
#' @param x A contrast_formula_spec object.
#' @param term A term object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing the term, name, weights, condition names, and contrast specification.
#'
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
  
  frm <- paste0(deparse(lazyeval::f_rhs(x$A)), collapse="")
  frm <- paste0(gsub(":", ".", deparse(lazyeval::f_rhs(x$A))), collapse="")
  A <- as.formula(paste("~", frm))
  
  if (is_continuous(term$evterm)) {
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

#' Contrast Difference Weights
#'
#' @description
#' Compute the contrast weights for a contrast_diff_spec object.
#'
#' @param x A contrast_diff_spec object.
#' @param term A term object.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing the term, name, weights, condition names, and contrast specification.
#'
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

#' Convert Contrast to GLT
#'
#' @description
#' Convert a contrast to an AFNI 'GLT' format.
#'
#' @param x The contrast to convert.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list containing the GLT string, name, and contrast specification.
#'
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

#' Write GLT to File
#'
#' @description
#' Write a generic function for writing GLT contrast to a file.
#'
#' @param x The GLT contrast object
#' @param ... Additional arguments passed to methods
#' @export
write_glt <- function(x, ...) {
  UseMethod("write_glt")
}

#' Write GLT to File
#'
#' @description
#' Write a GLT contrast to a file.
#'
#' @param x The GLT contrast to write.
#' @param fname The file name to write to.
#'
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

#' Estimated Contrast
#'
#' @description
#' Compute the estimated contrast for a given fit and indices.
#'
#' @param x The contrast to estimate.
#' @param fit The fit object.
#' @param indices The indices to use.
#' @param ... Additional arguments (currently unused).
#'
#' @return The estimated contrast.
#'
#' @noRd
#' @keywords internal
estcon.contrast <- function(x, fit, indices, ...) {
  wts <- numeric(length(fit$assign))
  wts[indices] <- x$weights
  
  gmodels::estimable(fit, wts)
}

#' Print Contrast Set
#'
#' @description
#' Print a contrast set.
#'
#' @param x The contrast set to print.
#' @param ... Additional arguments (currently unused).
#'
#' @export
print.contrast_set <- function(x, ...) {
  n_contrasts <- length(x)
  
  # Header
  cat("\n═══ Contrast Set ═══\n")
  
  # Summary
  cat("\n Overview:\n")
  cat("  • Number of contrasts:", n_contrasts, "\n")
  
  # Group contrasts by type
  types <- sapply(x, function(con) class(con)[1])
  type_table <- table(types)
  if (length(type_table) > 0) {
    cat("  • Types of contrasts:\n")
    for (type in names(type_table)) {
      cat("    -", type, ":", type_table[type], "\n")
    }
  }
  
  # List all contrasts
  cat("\n🔍 Individual Contrasts:\n")
  for (i in seq_along(x)) {
    cat("\n[", i, "] ", x[[i]]$name, " (", class(x[[i]])[1], ")\n", sep="")
    cat("    Formula: ")
    if (!is.null(x[[i]]$A)) cat(deparse(x[[i]]$A))
    if (!is.null(x[[i]]$B)) cat(" vs ", deparse(x[[i]]$B))
    cat("\n")
    if (!is.null(x[[i]]$where)) {
      cat("    Subset: ", deparse(x[[i]]$where), "\n")
    }
  }
  
  cat("\n")
  invisible(x)
}

#' Print Contrast Specification
#'
#' @description
#' Print a contrast specification.
#'
#' @param x The contrast specification to print.
#' @param ... Additional arguments (currently unused).
#'
#' @export
print.contrast_spec <- function(x,...) {
  cat("contrast:", x$name, "\n")
  cat(" A: ", Reduce(paste, deparse(x$A)), "\n")
  if (!is.null(x$B))
    cat(" B: ", Reduce(paste, deparse(x$B)), "\n")
  if (!is.null(x$where))
    cat(" where: ", Reduce(paste, deparse(x$where)), "\n")
  

}

#' Print Contrast
#'
#' @description
#' Print a contrast.
#'
#' @param x The contrast to print.
#' @param ... Additional arguments (currently unused).
#'
#' @export
print.contrast <- function(x,...) {
  print(x$contrast_spec)
  cat(" term: ", x$term$varname, "\n")
  cat(" weights: ", "\n")
  print(x$weights)
  cat(" conditions: ", x$condnames)
  
}

#' Print Polynomial Contrast Specification
#'
#' @description
#' Print a polynomial contrast specification.
#'
#' @param x The polynomial contrast specification to print.
#' @param ... Additional arguments (currently unused).
#'
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

#' Print Contrast Difference Specification
#'
#' @description
#' Print a contrast difference specification.
#'
#' @param x The contrast difference specification to print.
#' @param ... Additional arguments (currently unused).
#'
#' @export
print.contrast_diff_spec <- function(x,...) {
  cat("contrast difference:", "\n")
  cat("  ", x$con1$name, "-", x$con2$name, "\n")
}

#' plot_contrasts
#'
#' @description
#' Generic function for plotting contrasts.
#'
#' @param x Object containing contrast information
#' @param ... Additional arguments passed to methods
#' @export
plot_contrasts <- function(x, ...) {
  UseMethod("plot_contrasts")
}

#' plot_contrasts.event_model
#'
#' @description
#' Produces a heatmap of all contrasts defined for an \code{event_model}.
#' Rows = each contrast (or column of an F-contrast), columns = each regressor in
#' the full design matrix, and the fill color = the contrast weight.
#'
#' @param x An \code{event_model} with (lazily) defined contrasts.
#' @param absolute_limits Logical; if \code{TRUE}, the color scale is fixed at (-1,1).
#'   If \code{FALSE}, the range is set to (min, max) of the weights.
#' @param rotate_x_text Logical; if \code{TRUE}, rotate x-axis labels for readability.
#' @param ... Further arguments passed to \code{geom_tile}, e.g. \code{color="grey80"}.
#'
#' @return A \code{ggplot2} object (a heatmap).
#' @import ggplot2
#' @export
plot_contrasts.event_model <- function(
    x,
    absolute_limits = FALSE,
    rotate_x_text   = TRUE,
    scale_mode      = c("auto", "diverging", "one_sided"),
    coord_fixed     = TRUE,
    ...
) {
  # 1) Extract all the design-matrix column names
  dm <- design_matrix(x)
  regressor_names <- colnames(dm)
  
  # 2) Gather contrast weights (the nested list by term, then by contrast)
  cws <- contrast_weights(x)
  
  # Flatten everything into one big matrix of contrast weights
  big_mat  <- NULL
  rownames <- character(0)
  
  add_contrast_row <- function(vec, row_name) {
    if (is.null(big_mat)) {
      big_mat <<- matrix(vec, nrow = 1)
      colnames(big_mat) <<- regressor_names
      rownames <<- row_name
    } else {
      big_mat <<- rbind(big_mat, vec)
      rownames <<- c(rownames, row_name)
    }
  }
  
  for (term_nm in names(cws)) {
    term_level <- cws[[term_nm]]
    if (is.null(term_level)) next
    
    for (contrast_nm in names(term_level)) {
      cw_obj <- term_level[[contrast_nm]]
      # By default, we store offset_weights in cw_obj$offset_weights
      W <- cw_obj$offset_weights
      if (is.null(W)) next  # skip if no offset_weights
      
      # (#designCols x #contrastCols)
      ncols <- ncol(W)
      for (k in seq_len(ncols)) {
        this_col <- W[, k]
        if (ncols > 1) {
          row_label <- paste0(contrast_nm, "_component", k)
        } else {
          row_label <- contrast_nm
        }
        add_contrast_row(this_col, row_label)
      }
    }
  }
  
  if (is.null(big_mat)) {
    stop("No contrasts found in this event_model.")
  }
  
  rownames(big_mat) <- rownames
  
  # 3) Convert big_mat to a long data frame
  df_long <- as.data.frame(big_mat, check.names = FALSE)
  df_long$ContrastName <- rownames(big_mat)
  df_long <- reshape2::melt(
    df_long,
    id.vars      = "ContrastName",
    variable.name = "Regressor",
    value.name    = "Weight"
  )
  
  # 4) Build the ggplot
  plt <- ggplot2::ggplot(
    df_long,
    ggplot2::aes(
      x    = ReorderFactor(Regressor),
      y    = ReorderFactor(ContrastName, reverse = TRUE),
      fill = Weight
    )
  ) +
    ggplot2::geom_tile(...)
  
  # Decide on color scale
  scale_mode <- match.arg(scale_mode)
  wmin <- min(df_long$Weight, na.rm = TRUE)
  wmax <- max(df_long$Weight, na.rm = TRUE)
  
  # If user set absolute_limits=TRUE, we might forcibly use [-1,1] or [0,1]
  # but let's allow scale_mode to override as well.
  if (scale_mode == "diverging") {
    # Diverging scale, centered on 0
    lim_low  <- if (absolute_limits) -1 else wmin
    lim_high <- if (absolute_limits)  1 else wmax
    midpt <- 0
    
    plt <- plt + ggplot2::scale_fill_gradient2(
      limits   = c(lim_low, lim_high),
      midpoint = midpt,
      low      = "blue",
      mid      = "white",
      high     = "red"
    )
    
  } else if (scale_mode == "one_sided") {
    # One-sided scale for 0..1 or 0..something
    lim_low  <- if (absolute_limits) 0 else wmin
    lim_high <- if (absolute_limits) 1 else wmax
    
    plt <- plt + ggplot2::scale_fill_gradient(
      limits = c(lim_low, lim_high),
      low    = "white",
      high   = "red"
    )
    
  } else {
    # scale_mode == "auto"
    # If we detect any negative weight, do diverging; otherwise do one-sided
    if (wmin < 0) {
      # diverging
      lim_low  <- if (absolute_limits) -1 else wmin
      lim_high <- if (absolute_limits)  1 else wmax
      plt <- plt + ggplot2::scale_fill_gradient2(
        limits   = c(lim_low, lim_high),
        midpoint = 0,
        low      = "blue",
        mid      = "white",
        high     = "red"
      )
    } else {
      # one-sided
      lim_low  <- if (absolute_limits) 0 else wmin
      lim_high <- if (absolute_limits) 1 else wmax
      plt <- plt + ggplot2::scale_fill_gradient(
        limits = c(lim_low, lim_high),
        low    = "white",
        high   = "red"
      )
    }
  }
  
  # 5) Theming
  plt <- plt +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::labs(
      x    = "Regressor",
      y    = "Contrast",
      fill = "Weight"
    ) +
    ggplot2::theme(
      panel.grid  = ggplot2::element_blank(),
      axis.ticks  = ggplot2::element_blank()
    )
  
  if (rotate_x_text) {
    plt <- plt + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  }
  
  # 6) Optionally fix the coordinate ratio to keep tiles square
  if (coord_fixed) {
    plt <- plt +
      ggplot2::scale_x_discrete(expand = c(0, 0)) +
      ggplot2::scale_y_discrete(expand = c(0, 0)) +
      ggplot2::coord_fixed()
  }
  
  plt
}

#' A small utility to preserve factor order in ggplot
#' 
#' Makes a factor from a character vector but preserves the order of appearance.
#' If `reverse=TRUE`, it reverses that order.
#' @keywords internal
#' @noRd
ReorderFactor <- function(x, reverse=FALSE) {
  levs <- unique(as.character(x))
  if (reverse) levs <- rev(levs)
  factor(x, levels=levs)
}

#' Contrast Weights for a Contrast Set
#'
#' @description
#' Compute the contrast weights for each contrast specification within a contrast_set object.
#'
#' @param x A contrast_set object (a list of contrast_spec objects).
#' @param term A term object against which weights should be computed.
#' @param ... Additional arguments passed to individual contrast_weights methods.
#'
#' @return A named list where each element is the result of calling contrast_weights 
#'         on the corresponding contrast_spec in the set. The list names are the 
#'         names of the individual contrasts.
#'
#' @export
#' @importFrom purrr map set_names
contrast_weights.contrast_set <- function(x, term, ...) {
  # Ensure x is a list (contrast_set inherits from list)
  if (!is.list(x)) {
    stop("Input 'x' must be a contrast_set (list).")
  }
  
  # Iterate through each contrast spec in the set
  results_list <- purrr::map(x, function(contrast_spec) {
    # Check if the element is actually a contrast_spec
    if (!inherits(contrast_spec, "contrast_spec")) {
      warning(paste("Element", contrast_spec$name, "is not a contrast_spec object, skipping."))
      return(NULL)
    }
    # Compute weights for the individual contrast spec
    contrast_weights(contrast_spec, term, ...)
  })
  
  # Filter out any NULL results (if any elements weren't contrast_spec)
  results_list <- results_list[!sapply(results_list, is.null)]
  
  # Set the names of the results list based on the names of the contrasts
  contrast_names <- purrr::map_chr(results_list, "name")
  results_list <- purrr::set_names(results_list, contrast_names)
  
  return(results_list)
}



