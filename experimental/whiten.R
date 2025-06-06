#' @keywords internal
#' @autoglobal
#' @noRd
auto_whiten <- function(Y, modelmat, method=c("auto", "ar1", "ar2", "arma")) {
  with_package("forecast")
  method <- match.arg(method)
  out <- foreach( i = seq(1,ncol(Y)), .combine=cbind) %dopar% {
    y0 <- Y[,i]
    lm.1 = lm(y0 ~ modelmat-1)
    yresid <- resid(lm.1)
    
    afit <- if (method == "auto") {
      forecast::auto.arima(yresid, max.d=0, max.D=0, seasonal=FALSE)
    } else if (method == "ar1") {
      forecast::Arima(yresid, order=c(1,0,0), seasonal=FALSE)
    } else if (method == "ar2") {
      forecast::Arima(yresid, order=c(2,0,0), seasonal=FALSE)
    } else if (method== "arma") {
      forecast::Arima(yresid, order=c(1,0,1), seasonal=FALSE)
    }
    
    
    
    fitted(lm.1) + resid(afit)
    #list(phi=afit$model$phi, theta=afit$model$theta)
  }
  
  out
}

#' @keywords internal
#' @noRd
sq_inv_ar1_by_run <- function(phi, n) {
  fun <- memoise::memoise(sq_inv_ar1)
  out <- lapply(n, function(len) fun(phi, len))
  Matrix::bdiag(out)
}


#' @keywords internal
#' @noRd
sq_inv_ar1 <- function(phi, n) {
  phi_hat_vect <- rep(-phi, (n - 1))
  A <- Matrix::bandSparse(n,
                  k = c(1, 0),
                  diagonals = list(phi_hat_vect,
                                   c(sqrt(1 - phi ^ 2),
                                     rep(1, (
                                       n - 1
                                     )))))
  
  
}

#' @keywords internal
#' @noRd
sq_inv_arma_by_run <- function(phi, theta, n) {
  fun <- memoise::memoise(sq_inv_arma)
  out <- lapply(n, function(len) fun(phi, theta, len))
  Matrix::bdiag(out)
}

## not so slow if n is smaller
#' @keywords internal
#' @noRd
sq_inv_arma <- function(phi, theta, n) {
  #acf_theo_hat <- zapsmall(ARMAacf(ar = phi, ma = theta,lag.max = (n - 1)), digits=8)
  acf_theo_hat <- ARMAacf(ar = phi, ma = theta,lag.max = (n - 1))
  #ind <- which(acf_theo_hat != 0)
  #acf_theo_hat <- sparseVector(acf_theo_hat[ind], i=ind, length=length(acf_theo_hat))
  
  psi_hat <- ARMAtoMA(ar = phi,ma = theta, 10)
  variance_hat <- 1 + sum(psi_hat ^ 2)
  Sigma_hat <- Matrix::toeplitz(acf_theo_hat) * variance_hat
  S <- Matrix::Matrix(Sigma_hat, sparse=TRUE)
  A <- Matrix::Matrix(round(solve(chol(Sigma_hat)),digits = 4))
  A
}

# gen_signal <- function(n) {
#   y <- rnorm(n)
#   fitted(locfit(y ~ lp(seq(1,length(y)), nn=.03)))
#  
# }
# 
# n=800
# y=gen_signal(n)
# afit = auto.arima(y, seasonal=FALSE, max.q=1)
# phi <- afit$coef[1:afit$arma[1]]
# 
# nar <- afit$arma[1]
# nam <- afit$arma[2]
# theta <- afit$coef[(nar+1):(nar+1):nam]

