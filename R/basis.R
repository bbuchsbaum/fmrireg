

dctbasis <- function(n, p=n, const=FALSE) {
  m <- 1:n
  ret <- do.call(cbind, lapply(2:p, function(k) {
    (2/n)^(1/2) * cos(((2*m-1)*k*pi)/(2*n))
  }))
  
  if (const) {
    ret <- cbind((1/n)^.5, ret)
  }
  
  ret
}

#osplinebasis <- function() {
#  
#}


#' sub_basis
#' 
#' subset a parametric basis regressor
#' 
#' 
#' @param x the object
#' @param subset the subset
sub_basis <-  function(x, subset) UseMethod("sub_basis")


#' Ident
#' 
#' A basis that applies identity transform to a set of raw variables.
#' 
#' @param ... a list of variable names
#' @return an instance of class \code{Ident} extending \code{ParametricBasis}
#' @export 
Ident <- function(...) {
  mc <- match.call()
  vlist <- list(...)

  varnames <- as.list(substitute(list(...)))[-1]
  varnames <- unlist(lapply(varnames, as.character))
  names(vlist) <- varnames
  y <- do.call(cbind, vlist)
  
  ret <- list(x=vlist, y=y, varnames=varnames, name=paste(varnames, collapse="_"))
  class(ret) <- c("Ident", "ParametricBasis")
  ret
}


#' Polynomila basis
#' 
#' Orthogonal polynomial expansion of a linear term based on \code{\link[stats]{poly}}
#' 
#' @param x a numeric vector at which to evaluate the polynomial. Missing values are not allowed in x.
#' @param degree the degree of the polynomial. Must be less than the number of unique points.
#' @return an instance of class \code{Poly} extending \code{ParametricBasis}
#' 
#' @seealso \link[stats]{poly}
#' @export 
Poly <- function(x, degree) {
  mc <- match.call()
  pres <- poly(x,degree)
  n <- paste0("poly", "_", as.character(mc[["x"]]), "_", degree)
  ret <- list(x=x,y=pres,fun="poly",argname=as.character(mc[["x"]]), name=n, degree=degree)
  class(ret) <- c("Poly", "ParametricBasis")
  ret
}



#' B-spline basis
#' 
#' Generate the B-spline basis matrix for a polynomial spline.
#' 
#' @param x a numeric vector at which to evaluate the spline. Missing values are not allowed in x
#' @param degree the degree of the piecewise polynomial
#' @importFrom splines bs
#' @seealso \link[splines]{bs}
#' @export
BSpline <- function(x, degree) {
  mc <- match.call()
  
  pres <- bs(x,degree)
  n <- paste0("bs", "_", as.character(mc[["x"]]), "_", degree)[1]
  ret <- list(x=x,y=pres,fun="bs",argname=as.character(mc[["x"]])[1], name=n, degree=degree)
  class(ret) <- c("BSpline", "ParametricBasis")
  
  ret
}


#' @export
predict.Poly <- function(object,newdata,...) {
  predict(object$y, newdata)
}

sub_basis.Poly <- function(x, subset) {
  ret <- list(x=x$x[subset],
       y=x$y[subset,],
       fun=x$fun,
       argname=x$argname,
       name=x$name,
       degree=x$degree)
  class(ret) <- c("Poly", "ParametricBasis")
  ret
}

sub_basis.BSpline <- function(x, subset) {
  ret <- list(x=x$x[subset],
              y=x$y[subset,],
              fun=x$fun,
              argname=x$argname,
              name=x$name,
              degree=x$degree)
  class(ret) <- c("Poly", "ParametricBasis")
  ret
}

sub_basis.Ident <- function(x, subset) {
  vlist <- lapply(x$vlist, function(v) v[subset])
  ret <- list(vlist=vlist, varnames=x$varnames, name="Ident")
  class(ret) <- c("Ident", "ParametricBasis")
  ret
}

#' @export
predict.BSpline <- function(object,newdata,...) {
  predict(object$y, newdata)
}

#' @export
predict.Ident <- function(object,newdata,...) {
  ret <- as.data.frame(do.call(rbind, lapply(object$varnames, function(v) base::eval(v, newdata))))
  names(ret) <- object$varnames
  ret
}

#' @export
levels.Ident <- function(x) {
  x$varnames
}

#' @export
levels.BSpline <- function(x) {
  seq(1, x$degree)
}

#' @export
levels.Poly <- function(x) {
  seq(1,x$degree)
}

#' @export
columns.Poly <- function(x) {
  paste0(x$name, ".", seq(1, x$degree))
}

#' @export
columns.BSpline <- function(x) {
  paste0(x$name, ".", seq(1, x$degree))
}

#' @export
columns.Ident <- function(x) {
  x$varnames
}


