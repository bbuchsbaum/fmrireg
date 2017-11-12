


#' Poly
#' 
#' @param x a numeric vector at which to evaluate the polynomial. Missing values are not allowed in x
#' @param degree the degree of the polynomial. Must be less than the number of unique points.
#' @return an instance of class \code{Poly} extending \code{ParametricBasis}
#' @export 
Poly <- function(x, degree) {
  mc <- match.call()
  pres <- poly(x,degree)
  n <- paste0("poly", "_", as.character(mc[["x"]]), "_", degree)
  ret <- list(x=x,y=pres,fun="poly",argname=as.character(mc[["x"]]), name=n, degree=degree)
  class(ret) <- c("Poly", "ParametricBasis")
  ret
}



#' BSpline
#' 
#' @param x a numeric vector at which to evaluate the spline. Missing values are not allowed in x
#' @param degree the degree of the piecewise polynomial
#' @importFrom splines bs
#' @export
BSpline <- function(x, degree) {
  mc <- match.call()
  
  pres <- poly(x,degree)
  n <- paste0("bs", "_", as.character(mc[["x"]]), "_", degree)[1]
  ret <- list(x=x,y=pres,fun="bs",argname=as.character(mc[["x"]])[1], name=n, degree=degree)
  class(ret) <- c("BSpline", "ParametricBasis")
  
  ret
}

#' @export
predict.Poly <- function(object,newdata) {
  predict(object$y, newdata)
}

#' @export
predict.BSpline <- function(object,newdata) {
  predict(object$y, newdata)
}
