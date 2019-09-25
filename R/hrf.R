#' @import splines
NULL

#' Construct an \code{HRF} instance 
#' 
#' \code{HRF} takes a faw function f(t) and returns an \code{HRF} instance
#' 
#' @param hrf a function mapping from time --> signal
#' @param lag optional lag in seconds
#' @param width optional block width in seconds
#' @param precision sampling precision in seconds
#' @param ... extra parameters for the \code{hrf} function
#' 
#' @return an instance of type \code{HRF} inheriting from \code{function}
#' 
#' @examples 
#' 
#' ## generate and hrf using spmg canonical hrf, a lag of 3, and a width of 2.
#' grf <- gen_hrf(hrf_spmg1, lag=3, width=2)
#' grf(0:20)
#' 
#' hg <- purrr::partial(hrf_gaussian, mean=4, sd=1)
#' grf <- gen_hrf(hg, lag=1, width=2)
#' 
#' vals <- grf(0:20)
#' @export
gen_hrf <- function(hrf, lag=0, width=NULL, precision=.1, half_life=Inf, summate=TRUE, normalize=FALSE, name="gen_hrf", ...) {
  .orig <- list(...)
  
  if (!is.null(width)) {
    assert_that(width > 0)
    hrf <- gen_hrf_blocked(hrf, width=width, precision=precision, summate=summate, normalize=normalize, summate=summate)
  }
  
  if (lag > 0) {
    #force(hrf)
    hrf <- gen_hrf_lagged(hrf, lag=lag)
  }
  
  f <- if (length(.orig) > 0) {
    ret <- function(t) {
      do.call(hrf, c(list(t), .orig))
    }
    attr(ret, "params") <- .orig
    ret
  } else {
    hrf
  }
  
  vals <- f(0:2)
  
  nb <- if (is.vector(vals)) {
    1
  } else if (is.matrix(vals)) {
    ncol(vals)
  } else {
    stop("gen_hrf: constructed hrf is invalid")
  }
  HRF(f, name=name, nbasis=nb)
}


#' Generate an HRF basis set
#' 
#' \code{gen_hrf_set} construct an hrf basis set from a one or more component functions.
#' This function is used to create arbitrary sets of basis functions for fMRI modeling.
#' 
#' @param ... a list of functions f(t) mapping from time to amplitude 
#' @examples 
#' 
#' hrf1 <- hrf_spmg1 %>% gen_hrf(lag=0)
#' hrf2 <- hrf_spmg1 %>% gen_hrf(lag=3)
#' hrf3 <- hrf_spmg1 %>% gen_hrf(lag=6)
#' 
#' hset <- gen_hrf_set(hrf1,hrf2,hrf3)
#' @export
gen_hrf_set <- function(..., name="hrf_set") {
  hrflist <- list(...)
  assertthat::assert_that(all(sapply(hrflist, is.function)))
  f <- function(t) {
    do.call("cbind", lapply(hrflist, function(fun) fun(t)))
  }
  
  HRF(f, name=name, nbasis=length(hrflist))
}


#' HRF constructor function
#' 
#' a class used to represent a hemodynamic response function.
#' 
#' @param fun hemodynamic response function mapping from time --> BOLD response
#' @param name the name of the function
#' @param nbasis the number of basis functions, e.g. the columnar dimension of the hrf.
#' @param param_names the names of the parameters
#' @examples 
#' 
#' hrf <- HRF(hrf_gamma, "gamma", nbasis=1, param_names=c("shape", "rate"))
#' resp <- evaluate(hrf, seq(0,24,by=1))
#' 
#' @export
#' @rdname HRF
HRF <- function(fun, name, nbasis=1, param_names=NULL) {
  vals <- fun(0:30)
  
  if (nbasis == 1) {
    peak <- max(vals)
  } else {
    peak <- max(apply(vals, 2, max))
  }
  
  scale_factor <- 1/peak
  
  
  structure(fun, name=name, 
            nbasis=as.integer(nbasis), 
            param_names=param_names, 
            scale_factor=scale_factor, 
            class=c("HRF", "function"))
  
}

#' @inheritParams HRF
#' @export
AFNI_HRF <- function(name, nbasis, params) {
  structure(name, 
            nbasis=as.integer(nbasis), 
            params=params, 
            class=c("AFNI_HRF"))
  
}


#' @export
as.character.AFNI_HRF <- function(x) {
  paste(x, "\\(", paste(attr(x, "params"), collapse=","), "\\)", sep="")
}



#' @importFrom numDeriv grad
makeDeriv <- function(HRF, n=1) {
  if (n == 1) {
    function(t) numDeriv::grad(HRF, t)
  } else {
    Recall(function(t) numDeriv::grad(HRF,t), n-1)
  }
}


#' gen_hrf_lagged
#' 
#' @param hrf the underlying hrf function to shift
#' @param lag the lag/delay in seconds
#' @param ... extra args supplied to \code{hrf} function
#' @examples 
#' hrf_lag5 <- gen_hrf_lagged(HRF_SPMG1, lag=5)
#' hrf_lag5(0:20)
#' @export
gen_hrf_lagged <- function(hrf, lag=2,...) {
  force(hrf)
  
  if (length(lag)>1) {
   
    do.call(gen_hrf_set, lapply(lag, function(l) gen_hrf_lagged(hrf, l)))
  } else {
    function(t) {
      hrf(t-lag,...)
    }
  }
}

#' @export
#' @inheritParams gen_hrf_lagged
hrf_lagged <- gen_hrf_lagged


#' gen_hrf_blocked
#' 
#' @param hrf the hemodynmaic response function
#' @param width the width of the block
#' @param precision the sampling resolution
#' @param half_life the half_life of the exponential decay function (used to model response attenuation)
#' @param summate whether to allow each impulse response function to "add" up.
#' @param normalize rescale so that the peak of the output is 1.
#' @importFrom purrr partial
#' @export
gen_hrf_blocked <- function(hrf=hrf_gaussian, width=5, precision=.1, half_life=Inf, summate=TRUE, normalize=FALSE, ...) {
  force(hrf)
  purrr::partial(convolve_block, hrf=hrf, width=width, precision=precision, half_life=half_life, summate=summate, normalize=normalize, ...)
}

#' @export
#' @inheritParams gen_hrf_blocked
hrf_blocked <- gen_hrf_blocked


#' a convolve hemodynamic response with a block duration
#' 
#' apply a hemodynamic response function with times \code{t} and duration \code{width}
#' 
#' @param t time in seconds
#' @param hrf the underlying hemodynamic response function
#' @param width the fixed width of the response in seconds.
#' @param precision the sampling precision in seconds
#' @param half_life the half_life of the exponential decay function (used to model attenutation)
#' @param summate whether to allow each impulse response function to "add" up.
#' @param normalize rescale so that the peak of the output is 1.
#' @export
convolve_block <- function(t, hrf=hrf_gaussian, width=5, precision=.1, half_life=Inf, summate=TRUE, normalize=FALSE, ...) {
 
  hmat <- sapply(seq(0, width, by=precision), function(offset) {
    hrf(t-offset,...) * exp(-offset/half_life)
  })
  
  ret <- if (summate) {
    rowSums(hmat)
  } else {
    #r <- range(hmat[,1])
    #apply(hmat,1,function(vals) vals[which.max(abs(vals))])
    apply(hmat,1,function(vals) vals[which.max(vals)])
  }
  
  if (normalize) {
    ret <- ret/max(abs(ret))
  } 
  
  ret
}
  
  
#' hrf_time
#' 
#' @param t time in seconds
#' @param maxt the maximum time point in domain
#' @export
hrf_time <- function(t, maxt) {
  ifelse(t > 0 & t < maxt, t, 0)
}

#' hrf_bspline
#' 
#' @param t a vector of times
#' @param span the temporal window over which the basis sets spans
#' @param N the number of basis functions
#' @param degree the degree of the spline
#' @examples 
#' 
#' hrfb <- hrf_bspline(seq(0,20,by=.5), N=4, degree=2)
#' @export
#' @importFrom splines bs
hrf_bspline <- function(t, span=20, N=5, degree=3) {
	
	ord <- 1 + degree
	nIknots <- N - ord + 1
	if (nIknots < 0) {
		nIknots <- 0
		#warning("'df' was too small; have used  ", ord - (1 - intercept))
	}
	
	knots <- if (nIknots > 0) {
				knots <- seq.int(from = 0, to = 1, length.out = nIknots + 2)[-c(1, nIknots + 2)]
				stats::quantile(0:span, knots)
			} else {
				0
			}
	
	if (any(t < 0)) {
		t[t < 0] <- 0
	}
	
	if(any(t > span)) {
		t[t > span] <- 0
	}
	
	bs(t, df=N, knots=knots, degree=degree, Boundary.knots=c(0,span))
}


#' hrf_gamma
#' 
#' A hemodynamic response function using the gamma density function
#' 
#' @param t time
#' @param shape the shape parameter for gamma pdf
#' @param rate the rate parameter for gamma pdf
#' @export
hrf_gamma <- function(t, shape=6, rate=1) {
  dgamma(t, shape=shape, rate=rate)
}

#' hrf_gaussian
#' 
#' A hemodynamic response function using the gamma density function
#' 
#' @param t time
#' @param mean the mean of Gaussian pdf
#' @param sd the standard deviation of Gaussian pdf
#' @export
hrf_gaussian <- function(t, mean=6, sd=2) {
	dnorm(t, mean=mean, sd=sd)
}

#' hrf_spmg1
#' 
#' A hemodynamic response function based on the SPM canonical double gamma parameterzation.
#' 
#' @param t time
#' @export
hrf_spmg1 <- function(t, P1=5, P2=15) {
  A1=.00833
  A2=1.274527e-13

	ifelse(t < 0, 0, exp(-t)*(A1*t^P1 - A2*t^P2))
	
}

#' @export
#' @rdname HRF
HRF_GAMMA <- HRF(hrf_gamma, "gamma", param_names=c("shape", "rate"))

#' @export
#' @rdname HRF
HRF_GAUSSIAN <- HRF(hrf_gaussian, "gaussian", param_names=c("mean", "sd"))

#' @export
#' @rdname HRF
HRF_BSPLINE <- HRF(gen_hrf(hrf_bspline), "bspline", nbasis=5)

#' @export
#' @rdname HRF
HRF_SPMG1 <- HRF(hrf_spmg1, 
                 "SPMG1", param_names=c("A1", "A2"))

#' @export
#' @rdname HRF
HRF_SPMG2 <- HRF(gen_hrf_set(hrf_spmg1, makeDeriv(hrf_spmg1)), 
                 "SPMG2", nbasis=2, param_names=c("A1", "A2"))

#' @export
#' @rdname HRF
HRF_SPMG3 <- HRF(gen_hrf_set(hrf_spmg1, makeDeriv(hrf_spmg1), makeDeriv(makeDeriv(hrf_spmg1))), 
                 "SPMG3", nbasis=3, param_names=c("A1", "A2"))


#' evaluate
#' 
#' @param amplitude the scaling value for the event
#' @param duration the duration of the event
#' @param precision the temporal resolution used for computing summed responses when duration > 0 
#' @param summate whether the HRF response increases its amplitude as a function of stimulus duration.
#' @param normalize scale output so that peak is 1.
#' @rdname evaluate
#' @export
#' @examples 
#' 
#' hrf1 <- evaluate(HRF_SPMG1, grid=seq(0,20,by=1.5), duration=2, precision=1)
#' 
#' # the same, except now turn off temporal summation.
#' hrf2 <- evaluate(HRF_SPMG1, grid=seq(0,20,by=1.5), duration=2, precision=1,summate=FALSE)
#' 
evaluate.HRF <- function(x, grid, amplitude=1, duration=0, precision=.2, summate=TRUE, normalize=FALSE) {
  if (duration < precision) {
    x(grid)*amplitude*attr(x, "scale_factor")       
  } else if (nbasis(x) == 1) {
    samples <- seq(0, duration, by=precision)
    sfac <- attr(x, "scale_factor") 
    
    hmat <- sapply(samples, function(offset) {
      x(grid-offset)*amplitude*sfac
    })
    
    ret <- if (summate) {
      rowSums(hmat)
    } else {
      #rowMeans(hmat)
      apply(hmat,1,function(vals) vals[which.max(vals)])
    }
    
    if (normalize) {
      ret <- ret/max(abs(ret))
    }
    
    ret
  } else {
    sfac <- attr(x, "scale_factor")
    Reduce("+", lapply(seq(0, duration, by=precision), function(offset) {
      x(grid-offset)*amplitude*sfac   
    }))
  }
}


#' @export
evaluate.hrfspec <- function(x, grid, amplitude=1, duration=0, precision=.1) {
  evaluate(x$hrf, grid,amplitude, duration, precision)
}


#' @export
nbasis.HRF <- function(x) attr(x, "nbasis")

#' getHRF
#' 
#' @param name the name of the hrf function
#' @param nbasis the numbe rof basis functions (if relevant)
#' @export
getHRF <- function(name=c("gamma", "spmg1", "spmg2", "spmg3", "bspline", "gaussian", "tent", "bs"), nbasis=5, lag=0, ...) {
  name <- match.arg(name)
	
	hrf <- switch(name,
			gamma=HRF(gen_hrf_lagged(hrf_gamma,lag=lag),name="gamma"),
			gaussian=HRF(gen_hrf_lagged(HRF_GAUSSIAN,lag=lag), name="gaussian"),
			spmg1=HRF(gen_hrf_lagged(hrf_spmg1,lag=lag), name="spmg1", nbasis=1),
			spmg2=HRF(gen_hrf_lagged(HRF_SPMG2,lag=lag), name="spmg2", nbasis=2),
			spmg3=HRF(gen_hrf_lagged(HRF_SPMG3,lag=lag), name="spmg3", nbasis=3),
			tent=HRF(gen_hrf_lagged(hrf_bspline, lag=lag, N=nbasis,degree=1,...), name="bspline", nbasis=nbasis),
			bs=HRF(gen_hrf_lagged(hrf_bspline, lag=lag, N=nbasis,...), name="bspline", nbasis=nbasis),
			bspline=HRF(gen_hrf_lagged(hrf_bspline, lag=lag, N=nbasis,...), name="bspline", nbasis=nbasis))
	
	if (is.null(hrf)) {
		stop("could not find create hrf named: ", name)
	}
	
	hrf
}






#' hrf
#' 
#' hemodynamic regressor specification function for model formulas.
#' 
#' This function is to be used in formulas for fitting fucntions, e.g. onsets ~ hrf(fac1,fac2) ...
#' 
#' 
#' @param ... the variable names, all of which must be present in the enclosing environment (e.g. an \code{event_model} object)
#' @param basis the impulse response function or the name of a pre-supplied function, one of: "gamma", "spmg1", "spmg2", "spmg3", "bspline", "gaussian".
#' @param onsets optional onsets override. If missing, onsets will be taken from the \code{event_model}
#' @param durations optional durations override. If missing, onsets will be taken from the \code{event_model}
#' @param prefix a character string that is prepended to the variable names and used to identify the term. 
#'               Can be used to disambiguate two \code{hrf} terms with the same variable(s) but different onsets or basis functions.
#' @param subset an expression indicating the subset of 'onsets' to keep
#' @param precision sampling precision in seconds
#' @param nbasis number of basis functions -- only used for hemodynamic response functions (e.g. bspline) that take a variable number of bases.
#' @param contrasts one or more \code{contrast_spec} objects created with the \code{contrast} function. 
#' If multiple contrasts are required, then these should be wrapped in a \code{list} or \code{contrast_set}.
#' @param id a  unique \code{character} identifier used to refer to term, otherwise will be determined from variable names.
#' @param lag a temporal offset in seconds which is added to onset before convolution
#' 
#' @examples 
#' 
#' ## 'hrf' is typically used in the context of \code{formula}s.
#' hspec <- hrf(x)
#' hspec2 <- hrf(x, basis="gamma")
#' hspec3 <- hrf(x, basis="bs", nbasis=4)
#' 
#' form <- onsets ~ hrf(x) + hrf(y) + hrf(x,y)
#' 
#' @export
hrf <- function(..., basis="spmg1", onsets=NULL, durations=NULL, prefix=NULL, subset=NULL, precision=.2, 
                nbasis=1, contrasts=NULL, id=NULL, lag=0, summate=TRUE) {
  
  vars <- as.list(substitute(list(...)))[-1] 
  parsed <- parse_term(vars, "hrf")
  term <- parsed$term
  label <- parsed$label
  
  if (!is.numeric(lag) || length(lag) > 1) {
    stop("hrf: 'lag' must be a numeric scalar")
  }
  
  basis <- if (is.character(basis)) {
    getHRF(basis, nbasis=nbasis, lag=lag)
  } else if (inherits(basis, "HRF")) {
    if (lag > 0) {
      HRF(gen_hrf_lagged(basis, lag=lag), name=basis$name, nbasis=basis$nbasis)
    } else {
      basis
    }
    
  } else if (is.function(basis)) {
    test <- basis(1:10)
    HRF(gen_hrf_lagged(basis,lag=lag), name="custom_hrf", nbasis=ncol(test))
  } else {
    stop("invalid basis function: must be 1) character string indicating hrf type, e.g. 'gamma' 2) a function or 3) an object of class 'HRF': ", basis)
  }
  
  
  varnames <- if (!is.null(prefix)) {
    paste0(prefix, "_", term)
  } else {
    term
  }
  
  termname <- paste0(varnames, collapse="::")
  
  if (is.null(id)) {
    id <- termname
  }  
  
  cset <- if (inherits(contrasts, "contrast_spec")) {
    contrast_set(con1=contrasts)
  } else if (inherits(contrasts, "contrast_set")) {
    contrasts
  } 
  
  ret <- list(
    name=termname, ## hrf(x,y), where termname = "x::y"
    id=id, ## hrf(x), id by default is "x::y"
    varnames=varnames, ## list of all variables (e.g. list(x,y))
    vars=term, ## list of unparsed vars
    label=label, ## "hrf(x)" the full expression
    hrf=basis,
    onsets=onsets,
    durations=durations,
    prefix=prefix,
    subset=substitute(subset),
    precision=precision,
    lag=lag,
    contrasts=cset,
    summate=summate)
  
  class(ret) <- c("hrfspec", "list")
  ret
}

#' @keywords internal
construct_event_term <- function(x, model_spec) {
  ## TODO what if we are missing a block id?
  onsets <- if (!is.null(x$onsets)) x$onsets else model_spec$onsets
  durations <- if (!is.null(x$durations)) x$durations else model_spec$durations
  
  varlist <- lapply(seq_along(x$vars), function(i) {
    base::eval(parse(text=x$vars[[i]]), envir=model_spec$event_table, enclos=parent.frame())
  })
  
  names(varlist) <- x$varnames
  
  subs <- if (!is.null(x$subset)) {
    base::eval(x$subset, envir=model_spec$event_table, enclos=parent.frame()) 
  } else {
    rep(TRUE, length(onsets))
  }
  
  et <- event_term(varlist, onsets, model_spec$blockids, durations, subs)
}

#' @export
construct.hrfspec <- function(x, model_spec) {
  
  et <- construct_event_term(x,model_spec)
  
  cterm <- convolve(et, x$hrf, model_spec$sampling_frame, summate=x$summate)
  
  ret <- list(
    varname=et$varname,
    evterm=et,
    design_matrix=as.data.frame(cterm),
    sampling_frame=model_spec$sampling_frame,
    hrfspec=x,
    contrasts=x$contrasts,
    id=if(!is.null(x$id)) x$id else et$varname
  )
  
  class(ret) <- c("convolved_term", "fmri_term", "list") 
  ret
}


.hrf_parse <- function(..., prefix=NULL, basis=HRF_SPMG1, nbasis=1, lag=0, termsep="::") {
  vars <- as.list(substitute(list(...)))[-1] 
  browser()
  if (length(vars) > 0) {
    parsed <- parse_term(vars, "hrf")
    term <- parsed$term
    label <- parsed$label
  } else {
    stop("hrf: must have at least one variable")
  }
  
  basis <- if (is.character(basis)) {
    getHRF(basis, nbasis=nbasis, lag=lag)
  } else if (inherits(basis, "HRF")) {
    if (lag > 0) {
      HRF(gen_hrf_lagged(basis, lag=lag), name=basis$name)
    } else {
      basis
    }
    
  } else if (is.function(basis)) {
    test <- basis(1:10)
    HRF(gen_hrf_lagged(basis,lag=lag,...), name="custom_hrf", nbasis=ncol(test))
  } else {
    stop("invalid basis function: must be 1) character string indicating hrf type, e.g. 'gamma' 2) a function or 3) an object of class 'HRF': ", basis)
  }
  
  
  varnames <- if (!is.null(prefix)) {
    paste0(prefix, "_", parsed$term)
  } else {
    parsed$term
    
  }
  
  termname <- paste0(varnames, collapse=termsep)
  
  list(vars=vars, parsed=parsed, term=term, label=label, basis=basis, varnames=varnames, termname=termname)
}


# construct_additive_event_term <- function(x, model_spec) {
#   ## TODO what if we are missing a block id?
#   onsets <- if (!is.null(x$onsets)) x$onsets else model_spec$onsets
#   durations <- if (!is.null(x$durations)) x$durations else model_spec$durations
#   
#   varlist <- lapply(seq_along(x$vars), function(i) {
#     base::eval(parse(text=x$vars[[i]]), envir=model_spec$event_table, enclos=parent.frame())
#   })
#   
#   names(varlist) <- x$varnames
#   
#   subs <- if (!is.null(x$subset)) {
#     base::eval(x$subset, envir=model_spec$event_table, enclos=parent.frame()) 
#   } else {
#     rep(TRUE, length(onsets))
#   }
#   
#   mat <- do.call(cbind, varlist)
#   vlist <- list(mat)
#   names(vlist) <- x$name
#   
#   et <- event_term(vlist, onsets, model_spec$blockids, durations, subs)
# }


# hrf_add <- function(..., basis=HRF_SPMG1, onsets=NULL, durations=NULL,
#                     prefix=NULL, subset=NULL, precision=.2, nbasis=1,contrasts=list(), id=NULL) {
#   parsed <- .hrf_parse(..., prefix=prefix, basis=basis, nbasis=nbasis, termsep="+")
# 
#   if (is.null(id)) {
#     id <- parsed$termname
#   }
# 
# 
# 
#   ret <- list(
#     name=parsed$termname,
#     varnames=parsed$varnames,
#     vars=parsed$term,
#     label=parsed$label,
#     hrf=parsed$basis,
#     onsets=onsets,
#     durations=durations,
#     prefix=prefix,
#     subset=substitute(subset),
#     precision=precision,
#     contrasts=contrasts)
# 
#   class(ret) <- c("hrf_add_spec", "hrfspec", "list")
#   ret
# }
# 
# # @export
# construct.hrf_add_spec <- function(x, model_spec) {
#   et <- construct_additive_event_term(x, model_spec)
#   cterm <- convolve(et, x$hrf, model_spec$sampling_frame, summate=x$summate)
# 
#   ret <- list(
#     varname=et$varname,
#     evterm=et,
#     design_matrix=cterm,
#     sampling_frame=model_spec$sampling_frame,
#     contrasts=x$contrasts,
#     hrfspec=x,
#     id=x$id
#   )
# 
#   class(ret) <- c("convolved_term", "fmri_term", "list")
#   ret
# }


#' trialwise
#' 
#' This function is to be used in formulas for fitting functions, e.g. onsets ~ trialwise() ...
#' 
#' @inheritParams hrf
#' @examples 
#' 
#' 
#' 
#' @export
trialwise <- function(basis=HRF_SPMG1, onsets=NULL, durations=NULL, 
                      prefix=NULL, subset=NULL, precision=.2,
                      contrasts=list(), id=NULL, add_sum=FALSE) {
  
  termname = "trialwise"
  
  if (is.null(id)) {
    id <- termname
  }  
  
  if (is.function(basis) && !inherits(basis, "HRF")) {
    basis <- gen_hrf(basis)
  } else if (is.character(basis)) {
    basis <- getHRF(basis)
  } 
  
  assert_that(inherits(basis, "HRF"))
  
  ret <- list(
    name=termname,
    varnames=list("trialwise"),
    label="trialwise",
    hrf=basis,
    onsets=onsets,
    durations=durations,
    prefix=prefix,
    subset=substitute(subset),
    precision=precision,
    contrasts=contrasts,
    add_sum=add_sum)
  
  class(ret) <- c("trialwisespec", "hrfspec", "list")
  ret
}

#' @export
construct.trialwisespec <- function(x, model_spec) {
  
  ## compied almost verbatim from construct.hrfspec
  onsets <- if (!is.null(x$onsets)) x$onsets else model_spec$onsets
  durations <- if (!is.null(x$durations)) x$durations else model_spec$durations
  
  trial_index <- factor(seq(1, length(onsets)))
  
  varlist <- list(trial_index)
  names(varlist) <- x$varname
  
  subs <- if (!is.null(x$subset)) base::eval(x$subset, envir=model_spec$event_table, enclos=parent.frame()) else rep(TRUE, length(onsets))
  
  et <- event_term(varlist, onsets, model_spec$blockids, durations, subs)
  cterm <- convolve(et, x$hrf, model_spec$sampling_frame)
  
  ret <- list(
    varname=et$varname,
    evterm=et,
    design_matrix=cterm,
    sampling_frame=model_spec$sampling_frame,
    contrasts=x$contrasts,
    hrfspec=x,
    id=x$id
  )
  
  class(ret) <- c("trialwise_convolved_term", "convolved_term", "fmri_term", "list") 
  ret
}


#' construct an native AFNI hrf specification for '3dDeconvolve' with the 'stim_times' argument.
#' 
#' @inheritParams hrf
#' @examples 
#' 
#' @export
afni_hrf <- function(..., basis=c("spmg1", "block", "dmblock",           
                                  "tent",   "csplin", "poly",  "sin",        
                                  "gam", "spmg2", "spmg3", "wav"), 
                                  onsets=NULL, durations=NULL, prefix=NULL, subset=NULL, 
                                  nbasis=1, contrasts=NULL, id=NULL, start=NULL, stop=NULL) {
  
  ## TODO cryptic error message when argument is mispelled and is then added to ...
  basis <- match.arg(basis)
  
  vars <- as.list(substitute(list(...)))[-1] 
  parsed <- parse_term(vars, "afni_hrf")
  term <- parsed$term
  label <- parsed$label
  
  hrf <- if (!is.null(durations)) {
    assert_that(length(durations) == 1, msg="afni_hrf does not currently accept variable durations")
    get_AFNI_HRF(basis, nbasis=nbasis, duration=durations[1], b=start, c=stop)
  } else {
    get_AFNI_HRF(basis, nbasis=nbasis, b=start, c=stop)
  }
  
  
  varnames <- if (!is.null(prefix)) {
    paste0(prefix, "_", term)
  } else {
    term
  }
  
  termname <- paste0(varnames, collapse="::")
  
  if (is.null(id)) {
    id <- termname
  }  
  
  cset <- if (inherits(contrasts, "contrast_spec")) {
    contrast_set(con1=contrasts)
  } else if (inherits(contrasts, "contrast_set")) {
    contrasts
  } 
  
  ret <- list(
    name=termname,
    id=id,
    varnames=varnames,
    vars=term,
    label=label,
    hrf=hrf,
    onsets=onsets,
    durations=durations,
    prefix=prefix,
    subset=substitute(subset),
    lag=lag,
    contrasts=cset)
  
  class(ret) <- c("afni_hrfspec", "hrfspec", "list")
  ret
  
}

#' construct an native AFNI hrf specification for '3dDeconvolve' and indifdually modulated events using the 'stim_times_IM' argument.
#' 
#' @inheritParams hrf
#' @examples 
#' 
#' @export
afni_trialwise <- function(label, basis=c("spmg1", "block", "dmblock", "gam", "wav"), 
                     onsets=NULL, durations=NULL, subset=NULL, 
                      id=NULL, start=NULL, stop=NULL) {
  
  ## TODO cryptic error message when argument is mispelled and is then added to ...
  basis <- match.arg(basis)
  
  hrf <- if (!is.null(durations)) {
    assert_that(length(durations) == 1, msg="afni_trialwise does not currently accept variable durations")
    get_AFNI_HRF(basis, nbasis=1, duration=durations[1], b=start, c=stop)
  } else {
    get_AFNI_HRF(basis, nbasis=1, b=start, c=stop)
  }
  
  
  if (is.null(id)) {
    id <- label
  }  
  
  ret <- list(
    name=label,
    id=id,
    hrf=hrf,
    onsets=onsets,
    durations=durations,
    subset=substitute(subset))
  
  class(ret) <- c("afni_trialwise_hrfspec", "hrfspec", "list")
  ret
  
}

#' @export
construct.afni_hrfspec <- function(x, model_spec) {
  
  et <- construct_event_term(x, model_spec)
  
  ## do not convolve an afni term
  ##cterm <- convolve(et, x$hrf, model_spec$sampling_frame, summate=x$summate)
  
  ret <- list(
    varname=et$varname,
    evterm=et,
    sampling_frame=model_spec$sampling_frame,
    hrfspec=x,
    contrasts=x$contrasts,
    id=if(!is.null(x$id)) x$id else et$varname
  )
  
  class(ret) <- c("afni_hrf_convolved_term", "convolved_term", "fmri_term", "list") 
  ret
}



AFNI_SPMG1 <- function(d=1) AFNI_HRF(name="SPMG1", nbasis=as.integer(1), params=list(d=d)) 
AFNI_SPMG2 <- function(d=1) AFNI_HRF(name="SPMG2", nbasis=as.integer(2), params=list(d=d))
AFNI_SPMG3 <- function(d=1) AFNI_HRF(name="SPMG3", nbasis=as.integer(3), params=list(d=d))
AFNI_BLOCK <- function(d=1,p=1) AFNI_HRF(name="BLOCK", nbasis=as.integer(1), params=list(d=d,p=p))
AFNI_dmBLOCK <- function(d=1,p=1) AFNI_HRF(name="dmBLOCK", nbasis=as.integer(1), params=list(d=d,p=p))

AFNI_TENT <- function(b=0,c=18, n=10) AFNI_HRF(name="TENT", nbasis=as.integer(n), params=list(b=b,c=c,n=n))
AFNI_CSPLIN <- function(b=0,c=18, n=6) AFNI_HRF(name="CSPLIN", nbasis=as.integer(n), params=list(b=b,c=c,n=n))
AFNI_POLY <- function(b=0,c=18, n=10) AFNI_HRF(name="POLY", nbasis=as.integer(n), params=list(b=b,c=c,n=n))
AFNI_SIN <- function(b=0,c=18, n=10) AFNI_HRF(name="SIN", nbasis=as.integer(n), params=list(b=b,c=c,n=n))
AFNI_GAM <- function(p=8.6,q=.547) AFNI_HRF(name="GAM", nbasis=as.integer(1), params=list(p=p,q=q))
AFNI_WAV <- function(d=1) AFNI_HRF(name="WAV", nbasis=as.integer(1), params=list(d=1))


get_AFNI_HRF <- function(name, nbasis=1, duration=1, b=0, c=18) {
  hrf <- switch(name,
                gamma=AFNI_GAM(),
                spmg1=AFNI_SPMG1(d=duration),
                spmg2=AFNI_SPMG2(d=duration),
                spmg3=AFNI_SPMG3(d=duration),
                csplin=AFNI_CSPLIN(b=b,c=c, n=nbasis),
                poly=AFNI_POLY(b=b,c=c, n=nbasis),
                sine=AFNI_SIN(b=b,c=c,n=nbasis),
                wav=AFNI_WAV(),
                block=AFNI_BLOCK(d=duration),
                dmblock=AFNI_dmBLOCK())
  
  if (is.null(hrf)) {
    stop("could not find afni hrf named: ", name)
  }
  
  hrf
  
}


# hrf.logit <- function(t, a1=1, T1, T2, T3, D1, D2, D3) {
#    a2 <- a1 * (((inv.logit(-T3)/D3) - (inv.logit(-T1)/D1)) / ((inv.logit(-T3)/D3) + (inv.logit(-T2/D2))))
#    
#    a3 <- abs(a2) - abs(a1)
#    
#    a1*inv.logit((t-T1)/D1) + a2*inv.logit((t-T2)/D2) + a3*inv.logit((t-T3)/D3)
# }
# 
# hrf.logit2 <- function(t, a1, a2, a3, T1, T2, T3, D1, D2, D3) {
#   a1*inv.logit((t-T1)/D1) + a2*inv.logit((t-T2)/D2) + a3*inv.logit((t-T3)/D3)
# }
# 
# 
# getfun <- function(time) {
#   TIME <- time
#   ret <- function(par) {
#     hrf.logit2(TIME, par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9])
#   }
#   return(ret)
# }
# 
# minimize.fun <- function(yvals, fun, par) {
#   ypred <- fun(par)
#   return(sum((yvals-ypred)^2))
# }
# 
# get.minimizer <- function(yfun, yvals) {
#   YFUN <- yfun
#   YVALS <- yvals
#   ret <- function(par) {
#     ypred <- YFUN(par)
#     return(sum((YVALS-ypred)^2))
#   }
#   
#   return(ret)
# }
# 
# 
# shift.HRF <- function(HRF, shift) {
#   localShift <- shift
#   function(t) {
#     HRF(t+localShift)
#   }
# }
# 
# 
# 
# makeBlock <- function(HRF, duration) {
#   d1 <- duration
#   if (duration < 2) {
#     stop("duration must be greater than 1")
#   }
#   
#   funlist <- c(HRF, lapply(seq(-1, -(duration-1)), function(i) shift.HRF(HRF, i)))
#   function(t) {
#     ret <- numeric(length(t))
#     for (i in 1:length(t)) {
#       ret[i] <- sum(unlist(lapply(funlist, function(fun) fun(t[i]))))
#     }
#     ret
#     
#   }
# }
# 
# 
# makeBlockHRF <- function(eventOnset, duration, HRF) {
#   
#   localHRF <- HRF
#   localOnset <- eventOnset
#   onsets <- seq(eventOnset, eventOnset+duration, 1)
#   funlist <- lapply(onsets, function(onset) makeEventHRF(onset, localHRF))
#   
#   function(t) {
#     ret <- lapply(funlist, function(fun) fun(t)) 
#     Reduce("+", ret)
#   }
# }    
# 
# 
# makeEventHRF <- function(eventOnset, HRF, amp=1) {
#   
#   localHRF <- HRF
#   localOnset <- eventOnset
#   localAmp <- amp
#   function(t) {
#     
#     localHRF(t-localOnset)*amp
#     #for (i in 1:length(t)) {
#     #  ret[i] <- localHRF(t[i]-localOnset)*amp
#     #}
#     
#   }
# }
# 
# .makeEventHRF2 <- function(eventOnset, HRF, amp=1) {
#   
#   localHRF <- HRF
#   localOnset <- eventOnset
#   localAmp <- amp
#   function(t) {
#     ret <- numeric(length(t))
#     for (i in 1:length(t)) {
#       if (t[i] < localOnset) {
#         ret[i] <- 0
#       } else {
#         ret[i] <- localHRF(t[i]-localOnset)*amp
#       }
#     }
#     ret
#   }
# }
# 
# .makeEventHRF3 <- function(eventOnset, HRF, amp=1) {
#   
#   localHRF <- HRF
#   localOnset <- eventOnset
#   localAmp <- amp
#   function(t) {
#     if (t < localOnset) {
#       0
#     } else {
#       localHRF(t-localOnset)*amp
#     }
#     
#   }
# }


    
      
      
