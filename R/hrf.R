#' @import splines
NULL

#' gen_hrf
#' 
#' @param hrf a function mapping from time --> signal
#' @param lag optional lag in seconds
#' @param width optional block width in seconds
#' @param precision sampling precision in seconds
#' @param ... extra parameters for the \code{hrf} function
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
#' @export
gen_hrf <- function(hrf, lag=0, width=NULL, precision=.1, ...) {
  .orig <- list(...)
  
  if (!is.null(width)) {
    assert_that(width > 0)
    hrf <- gen_hrf_blocked(hrf, width=width, precision=precision)
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
  
  f
}


#' gen_hrf_set
#' 
#' @param ... a list of hrf functions
#' @export
gen_hrf_set <- function(...) {
  hrflist <- list(...)
  function(t) {
    do.call("cbind", lapply(hrflist, function(fun) fun(t)))
  }
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
  function(t) {
    hrf(t-lag,...)
  }
}


#' gen_hrf_blocked
#' 
#' @param hrf the hemodynmaic response function
#' @param width the width of the block
#' @param precision the sampling resolution
#' @importFrom purrr partial
#' @export
gen_hrf_blocked <- function(hrf=hrf_gaussian, width=5, precision=.1,...) {
  force(hrf)
  purrr::partial(hrf_block, hrf=hrf, width=width, precision=precision,...)
}

#' hrf_block
#' 
#' @param t time in seconds
#' @param hrf the underlying hemodynamic response function
#' @param width the fixed width of the response in seconds.
#' @param precision the sampling precision in seconds
#' @export
hrf_block <- function(t, hrf=hrf_gaussian, width=5, precision=.1,...) {
  Reduce("+", lapply(seq(0, width, by=precision), function(offset) {
    hrf(t-offset,...)
  }))
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
#' @rdname evaluate
#' @export
evaluate.HRF <- function(x, grid, amplitude=1, duration=0, precision=.1) {
  if (duration < precision) {
    x(grid)*amplitude*attr(x, "scale_factor")       
  } else if (nbasis(x) == 1) {
    rowSums(sapply(seq(0, duration, by=precision), function(offset) {
      x(grid-offset)*amplitude*attr(x, "scale_factor") 
    }))
    
    # TODO FIX ME
    #rowSums(map_dbl(seq(0, duration, by=precision), function(offset) {
    #            x(grid-offset)*amplitude*attr(x, "scale_factor") 
    #}))
  } else {
    Reduce("+", lapply(seq(0, duration, by=precision), function(offset) {
      x(grid-offset)*amplitude*attr(x, "scale_factor")   
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
                nbasis=1, contrasts=NULL, id=NULL, lag=0) {
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
    HRF(gen_hrf_lagged(basis,lag=lag,...), name="custom_hrf", nbasis=ncol(test))
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
    name=termname,
    id=id,
    varnames=varnames,
    vars=term,
    label=label,
    hrf=basis,
    onsets=onsets,
    durations=durations,
    prefix=prefix,
    subset=substitute(subset),
    precision=precision,
    lag=lag,
    contrasts=cset)
  
  class(ret) <- c("hrfspec", "list")
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
  
  subs <- if (!is.null(x$subset)) {
    base::eval(x$subset, envir=model_spec$event_table, enclos=parent.frame()) 
  } else {
    rep(TRUE, length(onsets))
  }
  
  et <- event_term(varlist, onsets, model_spec$blockids, durations, subs)
  
  cterm <- convolve(et, x$hrf, model_spec$sampling_frame)
  
  ret <- list(
    varname=et$varname,
    evterm=et,
    design_matrix=as.data.frame(cterm),
    sampling_frame=model_spec$sampling_frame,
    hrfspec=x,
    contrasts=x$contrasts,
    id=x$id
  )
  
  class(ret) <- c("convolved_term", "fmri_term", "list") 
  ret
}


.hrf_parse <- function(..., prefix=NULL, basis=HRF_SPMG1, nbasis=1, lag=0) {
  vars <- as.list(substitute(list(...)))[-1] 
  parsed <- parse_term(vars, "hrf")
  term <- parsed$term
  label <- parsed$label
  
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
    paste0(prefix, "_", term)
  } else {
    term
  }
  
  termname <- paste0(varnames, collapse="::")
  
  list(vars=vars, parsed=parsed, term=term, label=label, basis=basis, varnames=varnames, termname=termname)
}



#' trialwise
#' 
#' This function is to be used in formulas for fitting functions, e.g. onsets ~ trialwise(fac1) ...
#' 
#' @inheritParams hrf
#' @examples 
#' 
#' 
#' ## trialwise can be used with a factor with a single level because it splits each element in to separate regressor
#' fac <- factor(rep(1, 10))
#' twise <- trialwise(fac)
#' @export
trialwise <- function(..., basis=HRF_SPMG1, onsets=NULL, durations=NULL, 
                      prefix=NULL, subset=NULL, precision=.2, nbasis=1,contrasts=list(), id=NULL) {
  
  parsed <- .hrf_parse(..., prefix=prefix, basis=basis, nbasis=nbasis)
  
  if (is.null(id)) {
    id <- parsed$termname
  }  
  
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
  #sframe <- sampling_frame(model_spec$blocklens, model_spec$TR, model_spec$TR/2, x$precision)
  
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






# hrf.logit <- function(t, a1=1, T1, T2, T3, D1, D2, D3) {
#   a2 <- a1 * (((inv.logit(-T3)/D3) - (inv.logit(-T1)/D1)) / ((inv.logit(-T3)/D3) + (inv.logit(-T2/D2))))
#   
#   a3 <- abs(a2) - abs(a1)
#   
#   a1*inv.logit((t-T1)/D1) + a2*inv.logit((t-T2)/D2) + a3*inv.logit((t-T3)/D3)
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


    
      
      
