

extract_terms <- function(formula, data) {
  if (!inherits(formula, "terms")) {
    terms(formula, data = data)
  } else {
    formula
  }	
}

extract_covariates <- function(.terms, variables, resp, etab) {
  vars <- attr(.terms, "variables") 
  varnames <- sapply(vars, deparse, width.cutoff = 500)[-1]
  ind.vars <- varnames[-resp] 
  orig.covar.names <- ind.vars[which(sapply(variables[-resp], function(obj) is.numeric(obj) || .is.parametric.basis(obj)))]
  new.covar.names <- names(events(etab))[sapply(events(etab), isContinuous)]
  covar.names <- as.list(orig.covar.names)
  names(covar.names) <- new.covar.names
  covar.names
}

is.parametric.basis <- function(obj) { inherits(obj, "ParametricBasis") }

extract_variables <- function(.terms, data) {
  env <- environment(.terms)
  varnames <- sapply(attr(.terms, "variables") , deparse, width.cutoff = 500)[-1]
  variables <- eval(attr(.terms, "variables"), data, env)  
  names(variables) <- varnames
  variables
  
}



fmri_model <- function(formula, event_table, basis=HRF_SPMG1, durations=0, blockids, blocklens, TR, aux_data=data.frame(), drop_empty=TRUE) {
  stopifnot(inherits(formula, "formula"))
  
  vterms <- extract_terms(formula, event_table)
  resp <- attr(vterms, "response")
  assert_that(resp > 0)
  
  assert_that(all(blocklens>0))
  assert_that(length(TR) == 1)
  assert_that(TR > 0)
  
  assert_that(is.data.frame(aux_data))
  
  if (is.null(resp)) {
    stop("need to provide onset vector on left side of formula, e.g. Onsets ~  a + b")
  }
  
  if (missing(durations)) {
    durations <- rep(0, nrow(event_table))
  }
  
  variables <- extract_variables(vterms, event_table)
  lhs <- variables[[resp]]
  
  #assert_that(lhs %in% names(event_table))
  
  rhs <- variables[(resp+1):length(variables)]
  vclass <- sapply(rhs, class)
  
  model_spec <- list(formula=formula, event_table=event_table, onsets=lhs, varspec=rhs, varclass=vclass,
                     durations=durations, blocklens=blocklens, blockids=blockids, TR=TR, drop_empty=drop_empty,
                     aux_data=aux_data)
  
  class(model_spec) <- c("model_spec", "list")
  
  fmodel <- construct(model_spec)
  fmodel
}



parse_term <- function(vars, ttype) {
  dim <- length(vars) # number of variables
  term <- deparse(vars[[1]],backtick=TRUE) # first covariate
  if (dim>1) # then deal with further covariates
    for (i in 2:dim) term[i]<-deparse(vars[[i]],backtick=TRUE)
  for (i in 1:dim) term[i] <- attr(terms(reformulate(term[i])),"term.labels")
  
  full.call<-paste(ttype, "(",term[1],sep="")
  if (dim > 1) for (i in 2:dim) full.call<-paste(full.call,",",term[i],sep="")
  label <- paste(full.call,")",sep="")   # label for parameters of this term  
  
  list(term=term, label=label)
}

design_matrix.model_spec <- function(x) {
  termlist <- lapply(x$varspec, function(m) construct(m,x))
  ret <- lapply(termlist, design_matrix)
  vnames <- unlist(lapply(ret, names))
  dmat <- as.data.frame(do.call(cbind, ret))
  names(dmat) <- vnames
  dmat
}

design_matrix.fmri_model <- function(x) {
  ret <- lapply(x$terms, design_matrix)
  vnames <- unlist(lapply(ret, names))
  dmat <- as.data.frame(do.call(cbind, ret))
  names(dmat) <- vnames
  dmat
}

construct.model_spec <- function(x) {
  terms <- lapply(x$varspec, function(m) construct(m,x))
  
  ret <- list(
    terms=terms,
    model_spec=x
  )
  
  class(ret) <- c("fmri_model", "list")
  ret
}

terms.fmri_model <- function(x) {
  x$terms
}


#' @export
print.fmri_model <- function(object) {
  cat("fmri_model", "\n")
  cat(" ", "Formula:  ", as.character(object$model_spec$forumula), "\n")
  cat(" ", "Num Events: ", nrow(object$model_spec$event_table), "\n")
  #cat(" ", "Term Types: ", paste(sapply(object$events, function(ev) class(ev)[[1]])))
}


  

#' nuisance
#' 
#' a 'nuisance' term that consists of an arbitrary numeric matrix with the same number of rows as image time points.
#' 
#' @export
#' @param x a \code{matrix} 
#' @return a class of type \code{nuisancespec}
nuisance <- function(x) {
  varname <- substitute(x)
  
  ret <- list(
    name=varname
  )
  
  class(ret) <- "nuisancespec"
  ret
}






#' baseline
#' 
#' A matrix of polynomial regressors for modeling low-frequency drift in fmri time series.
#' @importFrom splines bs ns
#' @param number of polynomial terms for each image block
#' @param basis the type of polynomial basis.
#' @param name the name of the term
#' @export
baseline <- function(degree=5, basis=c("bs", "poly", "ns")[1], name=paste0("Baseline_", basis, "_", degree)) {
  bfun <- switch(basis,
                 bs=bs,
                 ns=ns,
                 poly=poly)
  
  ret <- list(
    degree=degree,
    fun=bfun,
    name=name
  )
  
  class(ret) <- c("baselinespec", "nuisancespec")
  ret
}


#' @export
block <- function(x) {
  varname <- substitute(x)
  ret <- list(
    name=varname
  )
  
  class(ret) <- "blockspec"
  ret
  
}


.hrf_parse <- function(..., prefix=NULL, basis=HRF_SPMG1, nbasis=1) {
  vars <- as.list(substitute(list(...)))[-1] 
  parsed <- parse_term(vars, "hrf")
  term <- parsed$term
  label <- parsed$label
  
  basis <- if (is.character(basis)) {
    getHRF(basis, nbasis=nbasis)
  } else if (is.function(basis)) {
    test <- basis(1:10)
    HRF(basis, name="custom_hrf", nbasis=ncol(test), ...)
  } else if (inherits(basis, "HRF")) {
    basis
  } else {
    stop("invalid basis function: must be 1) character string indicating hrf type, e.g. 'gamma' 2) a function or 3) an object of class 'HRF': ", basis)
  }
  
  
  varnames <- if (!is.null(prefix)) {
    paste0(prefix, "_", term)
  } else {
    term
  }
  
  termname <- paste0(varnames, collapse="::")
  
  list(vars=vars, parsed=parsed, term=term, label=label, basis=basis, varnames=varnames, termname=termname)
}
  
  
#' @export
trialwise <- function(..., basis=HRF_SPMG1, onsets=NULL, durations=NULL, prefix=NULL, subset=NULL, precision=.2, nbasis=1,contrasts=list()) {
 
  parsed <- .hrf_parse(..., prefix=prefix, basis=basis, nbasis=nbasis)
  
  ret <- list(
    name=parsed$termname,
    varnames=parsed$varnames,
    vars=parsed$term,
    label=parsed$label,
    hrf=parsed$basis,
    onsets=onsets,
    durations=durations,
    prefix=prefix,
    subset=substitute(subset),
    precision=precision,
    contrasts=contrasts)
  
  class(ret) <- c("trialwisespec", "hrfspec", "list")
  ret
}


#' hrf
#' 
#' hemodynamic regressor specification function
#' 
#' This function is to be used in formulas for fitting fucntions, e.g. onsets ~ hrf(fac1,fac2) ...
#' 
#' 
#' @param ... the variable names
#' @param basis the impulse response function.
#' @param onsets optional onsets override. If missing, onsets will be taken from global model specification duration evaluation.
#' @param durations optional durations override. If missing, onsets will be taken from global model specification during evlauation.
#' @param prefix a character string that is prepended to the variables names and used to identify the term.
#' @param subset
#' @param precision 
#' @param nbasis number of basis functions -- only used for hemodynamic response functions (e.g. bspline) that take a variable number of bases.
#' @param contrasts one or more \code{contrastspec} objects created with the \code{contrast} function. 
#' If multiple contrasts are required, then these should be wrapped in a \code{list}.
#' @export
hrf <- function(..., basis="spmg1", onsets=NULL, durations=NULL, prefix=NULL, subset=NULL, precision=.2, nbasis=1, contrasts=NULL) {
  vars <- as.list(substitute(list(...)))[-1] 
  parsed <- parse_term(vars, "hrf")
  term <- parsed$term
  label <- parsed$label
  
  basis <- if (is.character(basis)) {
    getHRF(basis, nbasis=nbasis)
  } else if (is.function(basis)) {
    test <- basis(1:10)
    HRF(basis, name="custom_hrf", nbasis=ncol(test), ...)
  } else if (inherits(basis, "HRF")) {
    basis
  } else {
    stop("invalid basis function: must be 1) character string indicating hrf type, e.g. 'gamma' 2) a function or 3) an object of class 'HRF': ", basis)
  }
  
  varnames <- if (!is.null(prefix)) {
    paste0(prefix, "_", term)
  } else {
    term
  }
  
  termname <- paste0(varnames, collapse="::")
  
  
  ret <- list(
    name=termname,
    varnames=varnames,
    vars=term,
    label=label,
    hrf=basis,
    onsets=onsets,
    durations=durations,
    prefix=prefix,
    subset=substitute(subset),
    precision=precision,
    contrasts=contrasts)
  
  class(ret) <- c("hrfspec", "list")
  ret
}

#' @export
construct.blockspec <- function(x, model_spec) {
  blockids <- base::eval(parse(text=x$name), envir=model_spec$aux_data, enclos=parent.frame())
  blockids <- as.factor(blockids)
  mat <- model.matrix(~ blockids - 1)
  colnames(mat) <- paste0(x$name, "_", levels(blockids))
  matrix_term(x$name, mat)
}

#' @export
construct.nuisancespec <- function(x, model_spec) {
  mat <- base::eval(parse(text=x$varname), envir=model_spec$aux_data, enclos=parent.frame())
  matrix_term(x$varname, mat)
}  

#' @export  
construct.baselinespec <- function(x, model_spec) {
  ret <- lapply(model_spec$blocklens, function(bl) cbind(rep(1, bl), x$fun(seq(1, bl), x$degree)))
  mat <- matrix(0, sum(model_spec$blocklens), (x$degree+1)*length(model_spec$blocklens))
  
  for (i in seq_along(ret)) {
    rowstart <- sum(model_spec$blocklens[1:i]) - model_spec$blocklens[1] + 1
    rowend <- rowstart + model_spec$blocklens[i] -1
    colstart <- (x$degree+1)*(i-1) + 1
    colend <- colstart + x$degree
    mat[rowstart:rowend, colstart:colend] <- ret[[i]]
  }
  

  cnames <- apply(expand.grid(paste0("Baseline#",1:length(model_spec$blocklens)), paste0("Run#", 1:(x$degree+1))), 1, paste, collapse="_")
  colnames(mat) <- cnames
  matrix_term(x$name, mat)	
}

#' @export
construct.trialwisespec <- function(x, model_spec) {
  ## compied almost verbatim from construct.hrfspec
  onsets <- if (!is.null(x$onsets)) x$onsets else model_spec$onsets
  durations <- if (!is.null(x$durations)) x$durations else model_spec$durations
  varlist <- lapply(seq_along(x$vars), function(i) {
    base::eval(parse(text=x$vars[[i]]), envir=model_spec$event_table, enclos=parent.frame())
  })
  
  
  ## syntheticlly adds '+trial_index+' variable
  trial_index <- factor(seq(1, length(onsets)))
  varlist <- c(varlist, list(trial_index))

  names(varlist) <- c(x$varnames, "trial_index")
  
  subs <- if (!is.null(x$subset)) base::eval(x$subset, envir=model_spec$event_table, enclos=parent.frame()) else rep(TRUE, length(onsets))
  
  et <- event_term(varlist, onsets, model_spec$blockids, durations, subs)
  sframe <- sampling_frame(model_spec$blocklens, model_spec$TR, model_spec$TR/2, x$precision)
  cterm <- convolve(et, x$hrf, sframe)
  ret <- list(
    evterm=et,
    design_matrix=cterm,
    sampling_frame=sframe,
    contrasts=x$contrasts,
    hrfspec=x
  )
  
  class(ret) <- c("trialwise_convolved_term", "convolved_term", "fmri_term", "list") 
  ret
}

#' @export
construct.hrfspec <- function(x, model_spec) {
  onsets <- if (!is.null(x$onsets)) x$onsets else model_spec$onsets
  durations <- if (!is.null(x$durations)) x$durations else model_spec$durations
 
  varlist <- lapply(seq_along(x$vars), function(i) {
    base::eval(parse(text=x$vars[[i]]), envir=model_spec$event_table, enclos=parent.frame())
  })
  
  names(varlist) <- x$varnames
  subs <- if (!is.null(x$subset)) base::eval(x$subset, envir=model_spec$event_table, enclos=parent.frame()) else rep(TRUE, length(onsets))
  
  et <- event_term(varlist, onsets, model_spec$blockids, durations, subs)
  sframe <- sampling_frame(model_spec$blocklens, model_spec$TR, model_spec$TR/2, x$precision)
  cterm <- convolve(et, x$hrf, sframe)
   
  ret <- list(
    evterm=et,
    design_matrix=as.data.frame(cterm),
    sampling_frame=sframe,
    hrfspec=x,
    contrasts=x$contrasts
  )
  
  class(ret) <- c("convolved_term", "fmri_term", "list") 
  ret
}

#' @export
design_matrix.convolved_term <- function(x) {
  x$design_matrix
}


#' @export
matrix_term <- function(varname, mat) {
  stopifnot(is.matrix(mat))
  ret <- list(varname=varname, design_matrix=mat)
  class(ret) <- c("matrix_term", "fmri_term", "list")
  ret
}

#' @export
design_matrix.matrix_term <- function(x,...) {
  if (is.null(colnames(x$design_matrix))) {
    cnames <- paste0(x$varname, "_", 1:ncol(x$design_matrix))
    colnames(x$design_matrix) <- cnames
  }
  
  dmat <- as.data.frame(x$design_matrix)
  dmat
}

#' @export
event_table.convolved_term <- function(x) event_table(x$evterm)

#' @export
nbasis.convolved_term <- function(x) nbasis(x$hrf)

#' @export
longnames.convolved_term <- function(x) {
  # ignores exclude.basis
  term.cells <- cells(x)
  # ignores exclude.basis
  apply(sapply(1:ncol(term.cells), 
    function(i) {
      paste(names(term.cells)[i], "#", term.cells[,i], sep="")
  }), 1, paste, collapse=":")

}
  

  
  
