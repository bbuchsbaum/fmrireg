
<<<<<<< HEAD

#' @export
=======
#' Poly
#' 
#' @param x a numeric vector at which to evaluate the polynomial. Missing values are not allowed in x
#' @param degree the degree of the polynomial. Must be less than the number of unique points.
#' @export polynomial basis wrapping \code{poly}
>>>>>>> 852d8b5091ae7f8c6c09adf3bdf3731924c18e72
Poly <- function(x, degree) {
  mc <- match.call()
  pres <- poly(x,degree)
  n <- paste0("poly", "_", as.character(mc[["x"]]), "_", degree)
  ret <- list(x=x,y=pres,fun="poly",argname=as.character(mc[["x"]]), name=n, degree=degree)
  class(ret) <- c("Poly", "ParametricBasis")
  ret
}

<<<<<<< HEAD
#' @export
#' @importFrom splines bs
=======
#' BSpline
#' 
#' @param x a numeric vector at which to evaluate the spline. Missing values are not allowed in x
#' @param degree the degree of the piecewise polynomial
>>>>>>> 852d8b5091ae7f8c6c09adf3bdf3731924c18e72
BSpline <- function(x, degree) {
  mc <- match.call()
  
  pres <- poly(x,degree)
  n <- paste0("bs", "_", as.character(mc[["x"]]), "_", degree)[1]
  ret <- list(x=x,y=pres,fun="bs",argname=as.character(mc[["x"]])[1], name=n, degree=degree)
  class(ret) <- c("BSpline", "ParametricBasis")
  
  #print(paste("mc:", as.character(mc["x"])))
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
