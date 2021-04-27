#' @keywords internal
#' @importFrom forecast auto.arima
#' @importFrom forecast Arima
auto_whiten <- function(Y, modelmat, method=c("auto", "ar1", "ar2", "arma")) {
  method <- match.arg(method)
  out <- foreach( i = seq(1,ncol(Y)), .combine=cbind) %dopar% {
    y0 <- Y[,i]
    lm.1 = lm(y0 ~ modelmat)
    yresid <- resid(lm.1)
    
    afit <- if (method == "auto") {
      auto.arima(yresid, max.d=0, max.D=0, seasonal=FALSE)
    } else if (method == "ar1") {
      Arima(yresid, order=c(1,0,0), seasonal=FALSE)
    } else if (method == "ar2") {
      Arima(yresid, order=c(2,0,0), seasonal=FALSE)
    } else if (method== "arma") {
      Arima(yresid, order=c(1,0,1), seasonal=FALSE)
    }
    
    fitted(lm.1) + resid(afit)
    #list(phi=afit$model$phi, theta=afit$model$theta)
  }
  
  out
}