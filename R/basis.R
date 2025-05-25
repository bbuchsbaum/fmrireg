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
  vlist <- list(...)

  varnames <- as.list(substitute(list(...)))[-1]
  varnames <- unlist(lapply(varnames, as.character))
  names(vlist) <- varnames
  y <- do.call(cbind, vlist)
  
  ret <- list(vlist=vlist, y=y, varnames=varnames, name=paste(varnames, collapse="_"))
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
  argname <- as.character(mc[["x"]])
  pres <- poly(x,degree)
  n <- paste0("poly_", argname)
  ret <- list(x=x,y=pres,fun="poly",argname=argname, 
              name=n,
              degree=degree)
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
#' Standardize a numeric vector by centering and scaling, handling NAs appropriately.
#' If the computed standard deviation is \code{NA} or zero, a small constant
#' (\code{1e-6}) is used instead to avoid division by zero.
#' The returned basis matrix has one column with this standardized name.
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
  if (is.na(x_sd) || x_sd == 0) x_sd <- 1e-6
  
  # Standardize, replacing NAs with 0
  standardized <- (x - x_mean) / x_sd
  standardized[is.na(standardized)] <- 0
  
  n <- paste0("std", "_", as.character(mc[["x"]]))
  ret <- list(x=x, y=matrix(standardized, ncol=1, dimnames=list(NULL, n)),
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
  sd_val <- object$sd
  if (is.na(sd_val) || sd_val == 0) sd_val <- 1e-6
  standardized <- (newdata - object$mean) / sd_val
  # Replace NAs with 0 to match training behavior
  standardized[is.na(standardized)] <- 0
  matrix(standardized, ncol=1, dimnames=list(NULL, object$name))
}

#' @export
sub_basis.Standardized <- function(x, subset) {
  ret <- list(x=x$x[subset],
              y=x$y[subset, , drop = FALSE],
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
  x$argname 
}

#' @export
columns.Standardized <- function(x) {
  continuous_token(x$name)
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
  argname <- as.character(mc[["x"]])[1]
  pres <- bs(x, degree = degree)
  n <- paste0("bs_", argname)
  ret <- list(x=x, y=pres, fun="bs", argname=argname, 
              name=n,
              degree=degree)
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
       degree=x$degree)
  class(ret) <- c("Poly", "ParametricBasis")
  ret
}

#' @export
sub_basis.BSpline <- function(x, subset) {
  ret <- list(x=x$x[subset],
              y=x$y[subset, , drop=FALSE],
              fun=x$fun,
              argname=x$argname,
              degree=x$degree)
  class(ret) <- c("BSpline", "ParametricBasis")
  ret
}

#' @export
sub_basis.Ident <- function(x, subset) {
  ret <- list(y = x$y[subset, , drop = FALSE], 
              varnames = x$varnames, 
              fun = "ident")
  class(ret) <- c("Ident", "ParametricBasis")
  ret
}

#' @export
predict.BSpline <- function(object, newdata, ...) {
  # Rebuild using bs() and stored attributes
  splines::bs(newdata, degree = object$degree,
     knots = attr(object$y, "knots"),
     Boundary.knots = attr(object$y, "Boundary.knots"))
}

#' @export
predict.Ident <- function(object,newdata,...) {
  if (!is.matrix(newdata) && !is.data.frame(newdata)) {
      stop("newdata for predict.Ident should be matrix or data.frame containing necessary columns")
  }
  return(as.matrix(newdata[, object$varnames, drop=FALSE]))
}

#' @export
levels.Ident <- function(x) {
  x$varnames
}

#' @export
levels.BSpline <- function(x) {
  ncols <- ncol(x$y)
  if (is.null(ncols) || ncols == 0) {
    return(character(0))
  }
  indices <- seq_len(ncols)
  zeropad(indices, ncols) # Pad based on the number of columns
}

#' @export
levels.Poly <- function(x) {
  ncols <- ncol(x$y) 
  if (is.null(ncols) || ncols == 0) {
    return(character(0))
  }
  indices <- seq_len(ncols) # Or seq_len(x$degree)
  zeropad(indices, x$degree) 
}

#' @export
columns.Poly <- function(x) {
  # Use degree for padding width
  indices <- seq_len(x$degree)
  padded_indices <- zeropad(indices, x$degree) # Enforces min width 
  # Return only the padded index as the condition tag component
  continuous_token(padded_indices)
}

#' @export
columns.BSpline <- function(x) {
  # Use degree for padding width
  indices <- seq_len(x$degree)
  padded_indices <- zeropad(indices, x$degree) # Enforces min width 
  # Return only the padded index as the condition tag component
  continuous_token(padded_indices)
}

#' @export
columns.Ident <- function(x) {
  continuous_token(x$varnames)
}


#' @export
nbasis.BSpline <- function(x) {
  x$degree
}

#' @export
nbasis.Poly <- function(x) {
  x$degree
}

#' Z-score (global) basis
#'
#' @param x numeric vector (NAs allowed)
#' @return object of class c("Scale","ParametricBasis")
#' @export
Scale <- function(x) {
  mc <- match.call()
  varname <- as.character(mc[["x"]])
  mu <- mean(x, na.rm = TRUE)
  sd_ <- stats::sd(x, na.rm = TRUE)
  if (is.na(sd_) || sd_ == 0) sd_ <- 1e-6 # Add guard for sd_ == 0 or NA
  z  <- (x - mu)/sd_
  z[is.na(z)] <- 0
  
  final_name <- paste0("z_", varname)
  ret <- list(x = x,
              y = matrix(z, ncol = 1, dimnames=list(NULL, final_name)),
              mean = mu,
              sd   = sd_,
              fun  = "scale",
              argname = varname,
              name = final_name)
  class(ret) <- c("Scale","ParametricBasis")
  ret
}

#' @export
#' @name predict.ParametricBasis
#' @title Predict Method for ParametricBasis Objects
#' @rdname predict.ParametricBasis
predict.Scale <- function(object, newdata, ...) {
  z <- (newdata - object$mean)/object$sd
  z[is.na(z)] <- 0
  matrix(z, ncol = 1, dimnames=list(NULL, object$name))
}

#' @export
#' @rdname sub_basis
sub_basis.Scale <- function(x, subset) {
  ret <- x
  ret$x <- x$x[subset]
  # Perform subsetting
  subsetted_y <- x$y[subset, , drop = FALSE]
  # Ensure the result is a matrix, specifically N x 1 if it was originally
  if (!is.matrix(subsetted_y)) {
    # This case indicates an unexpected dimension drop.
    # Reform it into an N x 1 matrix.
    # The number of rows should be length(ret$x) if subsetting was consistent.
    # warning("Dimension dropped unexpectedly during Standardized/Scale basis subsetting. Re-casting to matrix.", call.=FALSE)
    ret$y <- matrix(subsetted_y, ncol = 1, dimnames = list(NULL, colnames(x$y)))
  } else {
    ret$y <- subsetted_y
  }
  class(ret) <- class(x)
  ret
}

#' @export
#' @name levels.ParametricBasis
#' @title Get Levels/Names for ParametricBasis Objects
#' @rdname levels.ParametricBasis
levels.Scale  <- function(x) x$argname

#' @export
#' @name columns.ParametricBasis
#' @title Get Column Name Information for ParametricBasis Objects
#' @rdname columns.ParametricBasis
columns.Scale <- function(x) {
  continuous_token(x$name)
}


#' Z-score within groups
#' 
#' @param x numeric vector
#' @param g grouping factor / character / integer of same length as x
#' @export
ScaleWithin <- function(x, g) {
  mc <- match.call()
  varname <- as.character(mc[["x"]])
  grpname <- as.character(mc[["g"]])
  stopifnot(length(x) == length(g))
  g <- as.factor(g)
  
  # pre-compute means/sds per group (ignoring NAs)
  stats <- tapply(seq_along(x), g, function(idx) {
    mu <- mean(x[idx], na.rm = TRUE)
    sd_ <- stats::sd(x[idx], na.rm = TRUE)
    # Handle groups with sd=0 or only one non-NA value
    if (is.na(sd_) || sd_ == 0) sd_ <- 1e-6 # Use small value instead of 0 or NA
    c(mean = mu, sd = sd_)
  })
  mus <- sapply(stats, `[[`, "mean")
  sds <- sapply(stats, `[[`, "sd")
  
  z <- mapply(function(val, grp) {
              grp_name <- as.character(grp) # Ensure factor level name is used
              mu <- mus[grp_name]; sd_ <- sds[grp_name]
              if (is.na(val) || is.na(mu) || is.na(sd_)) 0
              else (val - mu) / sd_
            },
            x, g, SIMPLIFY = TRUE)
  
  final_name <- paste0("z_", varname, "_by_", grpname)
  ret <- list(x = x,
              group = g,
              y = matrix(z, ncol = 1, dimnames=list(NULL, final_name)),
              means = mus,
              sds   = sds,
              fun   = "scale_within",
              argname = varname,
              grpname = grpname,
              name = final_name)
  class(ret) <- c("ScaleWithin","ParametricBasis")
  ret
}

#' @export
#' @rdname predict.ParametricBasis
predict.ScaleWithin <- function(object, newdata, newgroup, ...) {
  stopifnot(length(newdata) == length(newgroup))
  newgroup <- as.factor(newgroup)
  
  z <- mapply(function(val, grp) {
      grp_name <- as.character(grp) # Ensure factor level name used for indexing
      mu <- object$means[grp_name]; sd_ <- object$sds[grp_name]
      if (is.na(val) || is.na(mu) || is.na(sd_)) 0 # Handle missing group stats in newdata
      else (val - mu)/sd_
    },
    newdata, newgroup, SIMPLIFY = TRUE)
  matrix(z, ncol = 1, dimnames=list(NULL, object$name))
}

#' @export
#' @rdname sub_basis
sub_basis.ScaleWithin <- function(x, subset) {
  ret <- x
  ret$x     <- x$x[subset]
  ret$group <- x$group[subset]
  ret$y     <- x$y[subset,,drop = FALSE]
  class(ret) <- class(x)
  ret
}

#' @export
#' @rdname levels.ParametricBasis
levels.ScaleWithin  <- function(x) x$argname

#' @export
#' @rdname columns.ParametricBasis
columns.ScaleWithin <- function(x) {
  continuous_token(x$name)
}


#' Robust Scaling (Median/MAD)
#'
#' @param x numeric vector (NAs allowed)
#' @return object of class c("RobustScale","ParametricBasis")
#' @export
RobustScale <- function(x) {
  mc <- match.call()
  varname <- as.character(mc[["x"]])
  med <- stats::median(x, na.rm = TRUE)
  mad_ <- stats::mad(x, na.rm = TRUE)  # Default constant=1.4826
  if (is.na(mad_) || mad_ == 0) mad_ <- 1e-6 # Handle zero MAD
  z <- (x - med)/mad_
  z[is.na(z)] <- 0
  
  final_name <- paste0("robz_", varname)
  ret <- list(x = x,
              y = matrix(z,  ncol = 1, dimnames=list(NULL, final_name)),
              median = med,
              mad    = mad_,
              fun = "robust_scale",
              argname = varname,
              name = final_name)
  class(ret) <- c("RobustScale","ParametricBasis")
  ret
}

#' @export
#' @rdname predict.ParametricBasis
predict.RobustScale <- function(object, newdata, ...) {
  z <- (newdata - object$median)/object$mad
  z[is.na(z)] <- 0
  matrix(z, ncol = 1, dimnames=list(NULL, object$name))
}

#' @export
#' @rdname sub_basis
# Use the same logic as Scale
sub_basis.RobustScale <- sub_basis.Scale

#' @export
#' @rdname levels.ParametricBasis
levels.RobustScale  <- function(x) x$argname

#' @export
#' @rdname columns.ParametricBasis
columns.RobustScale <- columns.Scale

# Add nbasis methods -----

#' @export
nbasis.Scale        <- function(x) 1L

#' @export
nbasis.ScaleWithin  <- function(x) 1L

#' @export
nbasis.RobustScale  <- function(x) 1L

#' @export
nbasis.Standardized <- function(x) 1L

#' @export
nbasis.Ident        <- function(x) ncol(x$y)

