#' @keywords internal
#' @importFrom forecast auto.arima
#' @importFrom forecast Arima
auto_whiten <- function(Y, modelmat, method=c("auto", "ar1", "ar2", "arma")) {
  method <- match.arg(method)
  out <- foreach( i = seq(1,ncol(Y)), .combine=cbind) %dopar% {
    y0 <- Y[,i]
    lm.1 = lm(y0 ~ modelmat-1)
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
    
    
    
    #browser()
    
    fitted(lm.1) + resid(afit)
    #list(phi=afit$model$phi, theta=afit$model$theta)
  }
  
  out
}

sq_inv_ar1 <- function(phi,n) {
  phi_hat_vect <- rep( - phi, (n - 1))
  A <- bandSparse(n, k = c(1, 0),
             diagonals = list(phi_hat_vect,
                              c(sqrt(1 - phi^2),
                                rep(1, (n - 1)))
             ))
  
}



sq_inv_arma <- function(phi, theta, n) {
  #acf_theo_hat <- zapsmall(ARMAacf(ar = phi, ma = theta,lag.max = (n - 1)))
  acf_theo_hat <- ARMAacf(ar = phi, ma = theta,lag.max = (n - 1))
  #ind <- which(acf_theo_hat != 0)
  #acf_theo_hat <- sparseVector(acf_theo_hat[ind], i=ind, length=length(acf_theo_hat))
  
  psi_hat <- ARMAtoMA(ar = phi,ma = theta, 100)
  variance_hat <- 1 + sum(psi_hat ^ 2)
  Sigma_hat <- Matrix::toeplitz(acf_theo_hat) * variance_hat
  S <- Matrix(Sigma_hat, sparse=TRUE)
  A <- Matrix(round(solve(chol(Sigma_hat)),digits = 4))
  A
}


