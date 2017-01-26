#' @import splines
NULL

#' helper function
#' @export
createHRF <- function(HRF, ...) {
  .orig <- list(...)
  
  f <- if (length(.orig) > 0) {
    ret <- function(t) {
      do.call(HRF, c(list(t), .orig))
    }
    attr(ret, "params") <- .orig
    ret
  } else {
    HRF
  }
  
  f
}

#' HRF constructor function
#' @rdname HRF-class
#' @param fun hemodynamic response function mapping from time --> BOLD response
#' @param name the name of the function
#' @param nbasis the number of basis, e.g. the columnar dimension of the response.
#' @export
HRF <- function(fun, name, nbasis=1, param_names=NULL) {
  ret <- list(hrf=fun, name=name, nbasis=as.integer(nbasis), param_names=param_names)
  
  class(ret) <- "HRF"
  ret
}


createHRFSet <- function(...) {
	hrflist <- list(...)
  function(t) {
    do.call("cbind", lapply(hrflist, function(fun) fun(t)))
  }
}

#' @importFrom numDeriv grad
makeDeriv <- function(HRF, n=1) {
  if (n == 1) {
    function(t) numDeriv::grad(HRF, t)
  } else {
    Recall(function(t) numDeriv::grad(HRF,t), n-1)
  }
}

#' @export
hrf_time <- function(t, maxt) {
  ifelse(t > 0 & t < maxt, t, 0)
}

#' @export
#' @importFrom splines bs
hrf_bspline <- function(t, width=20, N=5, degree=3) {
	
	ord <- 1 + degree
	nIknots <- N - ord + 1
	if (nIknots < 0) {
		nIknots <- 0
		#warning("'df' was too small; have used  ", ord - (1 - intercept))
	}
	
	knots <- if (nIknots > 0) {
				knots <- seq.int(from = 0, to = 1, length.out = nIknots + 2)[-c(1, nIknots + 2)]
				stats::quantile(0:width, knots)
			} else {
				0
			}
	
	if (any(t < 0)) {
		t[t < 0] <- 0
	}
	
	if(any(t > width)) {
		t[t > width] <- 0
	}
	
	bs(t, df=N, knots=knots, degree=degree, Boundary.knots=c(0,width))
}

#' @export
hrf_gamma <- function(t, shape=6, rate=1) {
  dgamma(t, shape=shape, rate=rate)
}

#' @export
hrf_gaussian <- function(t, mean=6, sd=2) {
	dnorm(t, mean=mean, sd=sd)
}

#' @export
hrf_spmg1 <- function(t, A1=.00833, A2=1.274527e-13, P1=5, P2=15) {
	ifelse(t < 0, 0, exp(-t)*(A1*t^P1 - A2*t^P2))
	
}

#' @export
HRF_GAMMA <- HRF(hrf_gamma, "gamma", param_names=c("shape", "rate"))

#' @export
HRF_GAUSSIAN <- HRF(hrf_gaussian, "gaussian", param_names=c("mean", "sd"))

#' @export
HRF_BSPLINE <- HRF(createHRF(hrf_bspline), "bspline", 5)

#' @export
HRF_SPMG1 <- HRF(hrf_spmg1, 
                 "SPMG1", param_names=c("A1", "A2"))

#' @export
HRF_SPMG2 <- HRF(createHRFSet(hrf_spmg1, makeDeriv(hrf_spmg1)), 
                 "SPMG2", nbasis=2, param_names=c("A1", "A2"))

#' @export
HRF_SPMG3 <- HRF(createHRFSet(hrf_spmg1, makeDeriv(hrf_spmg1), makeDeriv(makeDeriv(hrf_spmg1))), 
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
    x$hrf(grid)*amplitude       
  } else if (nbasis(x) == 1) {
    rowSums(sapply(seq(0, duration, by=precision), function(offset) {
                x$hrf(grid-offset)*amplitude
    }))
  } else {
    Reduce("+", lapply(seq(0, duration, by=precision), function(offset) {
      x$hrf(grid-offset)*amplitude
    }))
  }
}

evaluate.hrfspec <- function(x, grid, amplitude=1, duration=0, precision=.1) {
  evaluate(x$hrf, grid,amplitude, duration, precision)
}

nbasis.HRF <- function(x) x$nbasis

#' @export
getHRF <- function(name=c("gamma", "spmg1", "spmg2", "spmg3", "bspline"), nbasis=5) {
	
	hrf <- switch(name,
			gamma=HRF_GAMMA,
			gaussian=HRF_GAUSSIAN,
			spmg1=HRF_SPMG1,
			spmg2=HRF_SPMG2,
			spmg3=HRF_SPMG3,
			bspline=HRF(createHRF(hrf_bspline, N=nbasis), "bspline", nbasis))
	
	if (is.null(hrf)) {
		stop("could not find hrf named: ", name)
	}
	
	hrf
}

hrf.logit <- function(t, a1=1, T1, T2, T3, D1, D2, D3) {
  a2 <- a1 * (((inv.logit(-T3)/D3) - (inv.logit(-T1)/D1)) / ((inv.logit(-T3)/D3) + (inv.logit(-T2/D2))))
 
  a3 <- abs(a2) - abs(a1)
  
  a1*inv.logit((t-T1)/D1) + a2*inv.logit((t-T2)/D2) + a3*inv.logit((t-T3)/D3)
}

hrf.logit2 <- function(t, a1, a2, a3, T1, T2, T3, D1, D2, D3) {
  a1*inv.logit((t-T1)/D1) + a2*inv.logit((t-T2)/D2) + a3*inv.logit((t-T3)/D3)
}


getfun <- function(time) {
  TIME <- time
  ret <- function(par) {
    hrf.logit2(TIME, par[1], par[2], par[3], par[4], par[5], par[6], par[7], par[8], par[9])
  }
  return(ret)
}

minimize.fun <- function(yvals, fun, par) {
  ypred <- fun(par)
  return(sum((yvals-ypred)^2))
}

get.minimizer <- function(yfun, yvals) {
  YFUN <- yfun
  YVALS <- yvals
  ret <- function(par) {
    ypred <- YFUN(par)
    return(sum((YVALS-ypred)^2))
  }

  return(ret)
}


shift.HRF <- function(HRF, shift) {
  localShift <- shift
  function(t) {
    HRF(t+localShift)
  }
}



makeBlock <- function(HRF, duration) {
  d1 <- duration
  if (duration < 2) {
    stop("duration must be greater than 1")
  }

  funlist <- c(HRF, lapply(seq(-1, -(duration-1)), function(i) shift.HRF(HRF, i)))
  function(t) {
    ret <- numeric(length(t))
    for (i in 1:length(t)) {
      ret[i] <- sum(unlist(lapply(funlist, function(fun) fun(t[i]))))
    }
    ret
    
  }
}


makeBlockHRF <- function(eventOnset, duration, HRF) {
  
  localHRF <- HRF
  localOnset <- eventOnset
  onsets <- seq(eventOnset, eventOnset+duration, 1)
  funlist <- lapply(onsets, function(onset) makeEventHRF(onset, localHRF))
  
  function(t) {
    ret <- lapply(funlist, function(fun) fun(t)) 
    Reduce("+", ret)
  }
}    
  
  
makeEventHRF <- function(eventOnset, HRF, amp=1) {

  localHRF <- HRF
  localOnset <- eventOnset
  localAmp <- amp
  function(t) {
 
    localHRF(t-localOnset)*amp
    #for (i in 1:length(t)) {
    #  ret[i] <- localHRF(t[i]-localOnset)*amp
    #}
 
  }
}

.makeEventHRF2 <- function(eventOnset, HRF, amp=1) {

  localHRF <- HRF
  localOnset <- eventOnset
  localAmp <- amp
  function(t) {
    ret <- numeric(length(t))
    for (i in 1:length(t)) {
      if (t[i] < localOnset) {
        ret[i] <- 0
      } else {
        ret[i] <- localHRF(t[i]-localOnset)*amp
      }
    }
    ret
  }
}

.makeEventHRF3 <- function(eventOnset, HRF, amp=1) {

  localHRF <- HRF
  localOnset <- eventOnset
  localAmp <- amp
  function(t) {
    if (t < localOnset) {
      0
    } else {
      localHRF(t-localOnset)*amp
    }

  }
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
#' @param durations optional durations override. If missing, onsets will be taken from global model specification during evaluation.
#' @param prefix a character string that is prepended to the variables names and used to identify the term.
#' @param subset
#' @param precision 
#' @param nbasis number of basis functions -- only used for hemodynamic response functions (e.g. bspline) that take a variable number of bases.
#' @param contrasts one or more \code{contrastspec} objects created with the \code{contrast} function. 
#' If multiple contrasts are required, then these should be wrapped in a \code{list}.
#' @param id a  unique \code{character} identifier used to refer to term.
#' @export
hrf <- function(..., basis="spmg1", onsets=NULL, durations=NULL, prefix=NULL, subset=NULL, precision=.2, 
                nbasis=1, contrasts=NULL, id=NULL) {
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
  
  if (is.null(id)) {
    id <- termname
  }  
  
  cset <- if (inherits(contrasts, "contrast_spec")) {
    #vname <- deparse(substitute(contrasts))
    #eval(parse(text=paste0("contrast_set(", vname, "=contrasts)")))
    contrast_set(con1=contrasts)
  } else if (inherits(contrasts, "contrast_set")) {
    contrasts
  } #else if (!is.null(contrasts)) {
  ## try creating a contrast
  #vname <- deparse(substitute(contrasts))
  #eval(parse(text=paste0("contrast_set(", vname, "=contrasts)")))
  #contrast_set(con1)
  #}
  
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
  subs <- if (!is.null(x$subset)) base::eval(x$subset, envir=model_spec$event_table, enclos=parent.frame()) else rep(TRUE, length(onsets))
  
  et <- event_term(varlist, onsets, model_spec$blockids, durations, subs)
  #sframe <- sampling_frame(model_spec$sampling_frame$blocklens, model_spec$TR, model_spec$sampling_frame$TR/2, x$precision)
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







    
      
      
