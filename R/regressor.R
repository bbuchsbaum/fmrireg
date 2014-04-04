#' @importFrom RANN nn2
#' @include AllClass.R
NULL

#' Regressor constructor function
#' @param onset the event onsets in seconds
#' @param hrf a hemodynamic response function
#' @param duration duration of events (defualt is 0)
#' @param amplitude scaling vector (default is 1)
#' @param span the tmporal window of the response function (default is 20)
Regressor <- function(onsets, hrf, duration=0, amplitude=1, span=20) {
  if (length(duration) == 1) {
    duration = ConstantVector(as.vector(duration), length(onsets))
  }
  if (length(amplitude) == 1) {
    amplitude = ConstantVector(as.vector(amplitude), length(onsets))
  }
  
  new("Regressor", onsets=onsets,hrf=hrf, duration=duration,amplitude=amplitude,span=span)  
}

#' extract terms from formula
extractTerms <- function(formula, data) {
  if (!inherits(formula, "terms")) {
    terms(formula, data = data)
  } else {
    formula
  }	
}

dots <- function(...) {
  eval(substitute(alist(...)))
}

RegressorTerm <- function(..., data, duration=0, span=0) {
  varlist <- dots(...)
  
  vframe <- as.data.frame(lapply(varlist, function(v) eval(v, data, enclos=parent.frame())))
  vnames <- sapply(varlist, deparse)
  names(vframe) <- vnames
  vframe
  #varlist <- substitute(...)
  #varlist
  #varnames <- parse(text=match.call())
  #varnames <- as.list(varnames)[2:length(varnames)]
  #print(varnames)
  #print(varlist)
  
  ## make it all factors and then require a "by" argument
  ## by could by by=Poly(x,5))
  
  
}



#' @param repTime amount of time between each sample on the grid
#' @export
setMethod(f="evaluate", signature=signature(x = "Regressor", grid="numeric"),
          function (x, grid, repTime) {
         
            dspan <- x@span/repTime          
            outmat <- matrix(0, length(grid), length(x@onsets))
            nidx <- apply(RANN::nn2(matrix(grid), matrix(x@onsets), k=2)$nn.idx, 1, min)
            
            valid <- sapply(x@onsets, function(o) o > grid[1] && o < grid[length(grid)])
            valid.ons <- x@onsets[valid]
            nidx <- nidx[valid]
            end <- length(grid)
            
            for (i in seq_along(valid.ons)) { 
              grid.idx <- seq(nidx[i], min(nidx[i] + dspan, end))             
              relOns <- grid[grid.idx] - valid.ons[i]        
              outmat[grid.idx,i] <- x@hrf(relOns)                        
            }
            rowSums(outmat)
          })