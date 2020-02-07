library(rrpack)
library(pls)
library(RGCCA)
library(microbenchmark)

Y <- matrix(rnorm(1000*100), 1000, 100)
X <- matrix(rnorm(1000*100), 1000, 100)

dat <- list(Y=Y, X=X)

ncomp=3
pls1 <- microbenchmark(plsr(Y ~ X, ncomp=ncomp, data=dat))
cca1 <- microbenchmark(rgcca(A= list(Y, X),
      C = matrix(c(0, 1, 1, 0), 2, 2),
      tau = c(1, 1), ncomp=c(ncomp,ncomp), scale=FALSE, verbose=FALSE))



Q <- cca1$a[[1]]
Qt <- corpcor::pseudoinverse(Q)
lambda <- matrix(0, ncomp,ncomp)
diag(lambda) <- sqrt(cca1$AVE$AVE_inner_model)
#lambda = matrix(cca1$AVE$AVE_X[[2]], 1,1)
Bpls <- coef(pls1)[,,1]
Bcca <- cca1$a[[2]] %*% lambda %*% Qt

Ypred <- scale(X,center=TRUE) %*% Bcca
Ycen <- scale(Y,center=TRUE, scale=FALSE)
sum((Ycen - Ypred)^2)
sum(Ycen^2)

Ypred2 <- scale(X,scale=FALSE) %*% Bpls
sum((Ycen - Ypred2)^2)
