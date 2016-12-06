
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

#' a block variable, which is constant over the span of a scanning run
#' @value an instance of a class \code{blockspec}
#' @export
block <- function(x) {
  varname <- substitute(x)
  pterm <- parse_term(as.list(substitute(x)), "block")
  
  ret <- list(
    name=varname,
    label=pterm$label
  )
  
  class(ret) <- "blockspec"
  ret
  
}

#' @export
trialwise <- function(..., basis=HRF_SPMG1, onsets=NULL, durations=NULL, 
                      prefix=NULL, subset=NULL, precision=.2, nbasis=1,contrasts=list()) {
  
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

