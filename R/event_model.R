

extractTerms <- function(formula, data) {
  if (!inherits(formula, "terms")) {
    terms(formula, data = data)
  } else {
    formula
  }	
}



extractCovariates <- function(.terms, variables, resp, etab) {
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

extractVariables <- function(.terms, data) {
  
  env <- environment(.terms)
  varnames <- sapply(attr(.terms, "variables") , deparse, width.cutoff = 500)[-1]
  variables <- eval(attr(.terms, "variables"), data, env)  
  names(variables) <- varnames
  variables
  
}

createEventTerms <- function(.terms, variables, resp, etab, facnames, expmat) {
  covar.names <- extractCovariates(.terms, variables, resp, etab)
  facterm <- do.call(EventTerm, lapply(facnames, function(fac) etab[[fac]]))
  var.terms <- if(length(covar.names) > 0) .extractVarTerms(covar.names, facnames, expmat, etab) else NULL
  c(facterm,var.terms)
}



fmrilm <- function(formula, event_table, voxel_data, basis=HRF.SPMG1, drop.empty=TRUE) {
  stopifnot(inherits(formula, "formula"))
  
  vterms <- extractTerms(formula, event_table)
  resp <- attr(vterms, "response")
  stopifnot(resp > 0)
  
  if (is.null(resp)) {
    stop("need to provide onset vector on left side of formula, e.g. Onsets ~  a + b")
  }
  
  variables <- extractVariables(vterms, event_table)
  lhs <- variables[[resp]]
  
  rhs <- variables[(resp+1):length(variables)]
  vclass <- sapply(rhs, class)
  
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
    varname=varname
  )
  
  class(ret) <- "nuisancespec"
  ret
}

<<<<<<< HEAD

#' @importFrom splines bs ns
=======
#' baseline
#' 
#' A matrix of polynomial regressors for modleing low-frequency drift in fmri time series.
#' @importFrom splines bs, ns
#' @param number of polynomial terms for each image block
#' @param basis the type of polynomial basis.
#' @param name the name of the term
>>>>>>> 852d8b5091ae7f8c6c09adf3bdf3731924c18e72
#' @export
baseline <- function(degree=5, basis=c("bs", "poly", "ns")[1], name=paste0("Baseline_", basis, "_", degree)) {
  bfun <- switch(basis,
                 bs=bs,
                 ns=ns,
                 poly=poly)
  
  ret <- list(
    degree=degree,
    fun=bfun,
    varname=name
  )
  
  class(ret) <- c("baselinespec", "nuisancespec")
  ret
}

<<<<<<< HEAD

#' @export
=======
#' hrf
#' 
#' hemodynamic response specification
#' @param ... the variable names
#' @param basis the impulse response function.
#' @param onsets optional onsets override. If missing, onsets will be taken from overall model specification.
#' @param durations optional durations override. If missing, onsets will be taken from overall model specification.
#' @param prefix
#' @param subset
#' @param precision
>>>>>>> 852d8b5091ae7f8c6c09adf3bdf3731924c18e72
hrf <- function(..., basis=HRF.GAMMA, onsets=NULL, durations=NULL, prefix=NULL, subset=NULL, precision=.2) {
  vars <- as.list(substitute(list(...)))[-1] 
  parsed <- parse_term(vars, "hrf")
  term <- parsed$term
  label <- parsed$label
  
  basis <- if (is.character(basis)) {
    getHRF(basis)
  } else if (is.function(basis)) {
    test <- basis(1:10)
    HRF(basis, name="custom_hrf", nbasis=ncol(test))
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
  
  ret <- list(
    varnames=varnames,
    vars=term,
    label=label,
    hrf=basis,
    onsets=onsets,
    durations=durations,
    prefix=prefix,
    subset=substitute(subset))
  
  class(ret) <- "hrfspec"
  ret
}


#' @export
nuisancespec.construct <- function(x, model_spec) {
  mat <- base::eval(parse(text=x$varname), envir=model_spec$raw_table, enclos=parent.frame())
  matrix_term(x$varname, mat)
}  

#' @export  
baselinespec.construct <- function(x, model_spec) {
  ret <- lapply(model_spec$blocklens, function(bl) cbind(rep(1, bl), bfun(seq(1, bl), x$degree)))
  mat <- matrix(0, sum(model_spec$blocklens), (x$degree+1)*length(model_spec$blocklens))
  for (i in seq_along(ret)) {
    rowstart <- sum(model_spec$blocklens[1:i]) - model_spec$blocklens[1] + 1
    rowend <- rowstart + model_spec$blocklens[i] -1
    colstart <- (x$degree+1)*(i-1) + 1
    colend <- colstart + x$degree
    mat[rowstart:rowend, colstart:colend] <- ret[[i]]
  }
  
  cnames <- apply(expand.grid(paste("Baseline",1:length(blocklens),sep=""), paste("Run", 1:(x$degree+1), sep="")), 1, paste, collapse="_")
  colnames(mat) <- cnames
  matrix_term(x$varname, mat)	
}

#' @export
hrfspec.construct <- function(x, model_spec) {
  onsets <- if (!is.null(x$onsets)) {
    x$onsets
  }	else {
    model_spec$onsets
  }	
  
  durations <- if (!is.null(x$durations)) {
    x$durations
  } else {
    model_spec$durations
  }
 
 
  varlist <- lapply(seq_along(x$vars), function(i) {
    base::eval(parse(text=x$vars[[i]]), envir=model_spec$event_table, enclos=parent.frame())
  })
  
  names(varlist) <- x$varnames
  
  subs <- if (is.null(x$subset)) {
    NULL
  } else {
    base::eval(x$subset, envir=model_spec$event_table, enclos=parent.frame())
  }
  
  et <- event_term(varlist, onsets, model_spec$blockids, durations, subs)
  #sframe <- sampling_frame(blocklens, TR, startTime, x$precision)
  #cterm <- convolve(et, x$hrf, sframe)
  
  ret <- list(
    evterm=et,
    hrfspec=x
  )
  
  class(ret) <- c("convolvable_term", "fmri_term", "list") 
  ret
}
  

  
  
