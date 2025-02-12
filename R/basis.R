#' @noRd 
#' @keywords internal
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
#' @export
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


#' Polynomial basis
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

# Poly <- function(x, degree) {
#   mc <- match.call()
#   
#   # If x is already numeric, do the polynomial expansion now (old behavior)
#   if (is.numeric(x)) {
#     pres <- stats::poly(x, degree)
#   } else {
#     # If x is not numeric (likely a symbol like `modulator`),
#     # defer actual polynomial expansion until sub_basis.Poly() or similar.
#     pres <- NULL
#   }
#   
#   n <- paste0("poly_", as.character(mc[["x"]]), "_", degree)
#   
#   ret <- list(
#     x       = x,       # might be a numeric vector or a symbol
#     y       = pres,    # either the poly(...) matrix OR NULL
#     fun     = "poly",
#     argname = as.character(mc[["x"]]),
#     name    = n,
#     degree  = degree,
#     call    = mc       # store original call for lazy eval
#   )
#   
#   class(ret) <- c("Poly", "ParametricBasis")
#   ret
# }

#' Standardized basis
#' 
#' Standardize a numeric vector by centering and scaling, handling NAs appropriately
#' 
#' @param x a numeric vector to standardize. Missing values are allowed and will be replaced with 0 after standardization.
#' @return an instance of class \code{Standardized} extending \code{ParametricBasis}
#' @export 
Standardized <- function(x) {
  mc <- match.call()
  
  # Remove NAs for computing mean and sd
  x_clean <- x[!is.na(x)]
  x_mean <- mean(x_clean)
  x_sd <- sd(x_clean)
  
  # Standardize, replacing NAs with 0
  standardized <- (x - x_mean) / x_sd
  standardized[is.na(standardized)] <- 0
  
  n <- paste0("std", "_", as.character(mc[["x"]]))
  ret <- list(x=x, y=matrix(standardized, ncol=1), 
              mean=x_mean, sd=x_sd,
              fun="standardized", 
              argname=as.character(mc[["x"]]), 
              name=n)
  class(ret) <- c("Standardized", "ParametricBasis")
  ret
}

#' @export
predict.Standardized <- function(object, newdata, ...) {
  # Standardize new data using stored mean and sd
  standardized <- (newdata - object$mean) / object$sd
  # Replace NAs with 0 to match training behavior
  standardized[is.na(standardized)] <- 0
  matrix(standardized, ncol=1)
}

#' @export
sub_basis.Standardized <- function(x, subset) {
  ret <- list(x=x$x[subset],
              y=x$y[subset,],
              fun=x$fun,
              argname=x$argname,
              name=x$name,
              mean=x$mean,
              sd=x$sd)
  class(ret) <- c("Standardized", "ParametricBasis")
  ret
}

#' @export
levels.Standardized <- function(x) {
  x$name
}

#' @export
columns.Standardized <- function(x) {
  paste0(x$name)
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
#' @return an \code{BSpline} list instance
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

#' @export
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

#' @export
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

#' @export
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

#' @noRd
columns.Poly <- function(x) {
  paste0(x$name, ".", seq(1, x$degree))
}

#' @noRd
columns.BSpline <- function(x) {
  paste0(x$name, ".", seq(1, x$degree))
}

#' @noRd
columns.Ident <- function(x) {
  x$varnames
}


