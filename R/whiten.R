#' @keywords internal
auto_whiten <- function(Y, modelmat) {
  
  out <- foreach( i = seq(1,ncol(Y)), .combine=cbind) %dopar% {
    y0 <- Y[,i]
    lm.1 = lm(y0 ~ modelmat)
    yresid <- resid(lm.1)
    afit <- auto.arima(yresid, seasonal=FALSE)
    y1 <- fitted(lm.1) + resid(afit)
    #list(phi=afit$model$phi, theta=afit$model$theta)
  }
  
  out
}