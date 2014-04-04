#' constructor function for \code{ConstantVector}
#' @param constant the constant value
#' @param length the length of the vector
ConstantVector <- function(value, length) {
  if (length <= 0) {
    stop("length must be greater than 0")
  }
  new("ConstantVector", .Data=value, length=as.integer(length))
}

#' @export
setMethod(f="[", signature=signature(x = "ConstantVector", i = "numeric", j = "missing", drop = "ANY"),
          function (x, i, j, drop) {
            if (any(i > x@length)) {
              stop(paste("index range", max(i), "is outside vector bounds", x@length))
            } else {
              rep(x@.Data[1], length(i))     
            }
          })

#' @export
setMethod(f="[[", signature=signature(x = "ConstantVector", i = "numeric", j = "missing"),
          function (x, i, j) {
            if (any(i > x@length)) {
              stop(paste("index range", max(i), "is outside than vector bounds", x@length))
            } else {
              rep(x@.Data[1], length(i))        
            }
          })




              