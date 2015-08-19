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
HRF <- function(fun, name, nbasis=1) {
  ret <- list(hrf=fun, name=name, nbasis=as.integer(nbasis))
  
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
hrf.time <- function(t, maxt) {
  ifelse(t > 0 & t < maxt, t, 0)
}

#' @export
#' @importFrom splines bs
hrf.bspline <- function(t, width=20, N=5, degree=3) {
	
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
hrf.gamma <- function(t, shape=6, rate=1) {
  dgamma(t, shape=shape, rate=rate)
}

#' @export
hrf.gaussian <- function(t, mean=6, sd=2) {
	dnorm(t, mean=mean, sd=sd)
}

#' @export
hrf.spmg1 <- function(t, A1=.00833, A2=1.274527e-13, P1=5, P2=15) {
	ifelse(t < 0, 0, exp(-t)*(A1*t^P1 - A2*t^P2))
	
}

#' @export
HRF.GAMMA <- HRF(hrf.gamma, "gamma")

#' @export
HRF.GAUSSIAN <- HRF(hrf.gaussian, "gaussian")

#' @export
HRF.BSPLINE <- HRF(createHRF(hrf.bspline), "bspline", 5)

#' @export
HRF.SPMG1 <- HRF(hrf.spmg1, "SPMG1")

#' @export
HRF.SPMG2 <- HRF(createHRFSet(hrf.spmg1, makeDeriv(hrf.spmg1)), "SPMG2")

#' @export
HRF.SPMG3 <- HRF(createHRFSet(hrf.spmg1, makeDeriv(hrf.spmg1), makeDeriv(makeDeriv(hrf.spmg1))), "SPMG3")

#' evaluate an HRF for a single event along a set of points
#' @param grid the sampling grid, e.g. the points in time at which the HRF will be evaluated.
#' @param amplitude the scaling value for the event
#' @param duration the duration of the event
#' @param precision the temporal resolution used for computing summed responses when duration > 0 
#' @export
evaluate.HRF <- function(x, grid, amplitude=1, duration=0, precision=.1) {
  if (duration < precision) {
    x$hrf(grid)*amplitude       
  } else {
    rowSums(sapply(seq(0, duration, by=precision), function(offset) {
                x$hrf(grid-offset)*amplitude
    }))
  }
}

nbasis.HRF <- function(x) x$nbasis

#' @export
getHRF <- function(name=c("gamma", "spmg1", "spmg2", "spmg3", "bspline"), ...) {
	
	hrf <- switch(name,
			gamma=create.HRF(hrf.gamma, ...),
			gaussian=create.HRF(hrf.gaussian, ...),
			spmg1=HRF.SPMG1,
			spmg2=HRF.SPMG2,
			spmg3=HRF.SPMG3,
			bspline=create.HRF(hrf.bspline, ...))		
	
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
    #browser()
    
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


    
      



    






    
      
      
