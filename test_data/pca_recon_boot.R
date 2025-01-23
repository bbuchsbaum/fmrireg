library(mvtnorm)
library(multivarious)
library(boot)

X <- matrix(rnorm(80 * 10), 80, 10)
X <- scale(X, center=TRUE, scale=FALSE)
pcres <- pca(X, ncomp = 10, preproc=pass())
ndvar=5

y0 <- X[,1]

## design matrix
Xd <- rmvnorm(nrow(X), sigma = diag(seq(ndvar, 1)))
Xd <- scale(Xd, center=TRUE)

set.seed(4237)
lm.boot <- boot(y0, function(dat, ind) {
  coef(lm(dat[ind] ~ Xd[ind,] - 1))[1]
}, R=50)



se.boot.orig <- sd(lm.boot$t)

lm.orig <- lm(y0 ~ Xd-1)


### lm in pca space
fits <- lapply(1:length(sdev(pcres)), function(i) {
  y <- scores(pcres)[,i]
  lm(y ~ Xd -1)
})

bvals <- unlist(lapply(fits, function(x) coef(x)[1]))

## reconstructed betas
b_lds <- coef(pcres)
brec <- (t(as.matrix(bvals)) %*% t(b_lds))[1]
abs(brec[1] - coef(lm.orig)[1]) < .1

gen_seboot <- function(R) {
  print(R)
  
  fits <- lapply(1:length(sdev(pcres)), function(i) {
    y <- scores(pcres)[,i]
    lm(y ~ Xd - 1)
  })
  
  boot_fits <- lapply(1:500, function(i) {
    ind <- sample(1:nrow(Xd), replace=TRUE)
    tmp <- scores(pcres)[,1]
    f0 <- lm(tmp ~ Xd[ind,] - 1)
    Qr <- qr.lm(f0)
    covu <- chol2inv(Qr$qr)
    ret <- lapply(1:length(sdev(pcres)), function(j) {
      y <- scores(pcres)[,j]
      f <- lm(y[ind] ~ Xd[ind,] - 1)
      
      list(beta=coef(f)[1], res=resid(f), covu=covu)
    })
    
    betas <- sapply(ret, "[[", "beta")
    resid <- do.call(cbind, lapply(ret, "[[", "res"))
    
    list(betas=betas, R=cov(resid), covu=covu)
  })
  # boot_fits <- lapply(1:length(sdev(pcres)), function(i) {
  #   y <- scores(pcres)[,i]
  #   lm.boot <- boot(y, function(dat, ind) {
  #     coef(lm(dat[ind] ~ Xd[ind,] - 1))[1]
  #   }, R=R)
  #   
  # })
  betas <- do.call(rbind, lapply(boot_fits, "[[", "betas"))
  seboot <- apply(betas,2,function(vals) sd(vals))
  varboot <- apply(betas,2,function(vals) var(vals))
  
  pooledR <- Reduce("+", lapply(boot_fits, "[[", "R"))/(length(boot_fits))
  pooledRCorrected <- pooledR * (nrow(X)-1)/lm.orig$df.residual
  Qr <- qr.lm(lm.orig)
  covu <- chol2inv(Qr$qr)
  
  pooledCov <- Reduce("+", lapply(boot_fits, "[[", "covu"))/(length(boot_fits))
  bootrec1 <- sqrt(b_lds[1,] %*% (pooledR * covu[1,1]) %*% b_lds[1,])
  bootrec2 <- sqrt(b_lds[1,] %*% (pooledRCorrected * covu[1,1]) %*% b_lds[1,])
  bootrec3 <- sqrt(b_lds[1,] %*% (pooledRCorrected * pooledCov[1,1]) %*% b_lds[1,])
  bootrec4 <- sqrt(b_lds[1,] %*% (pooledR * pooledCov[1,1]) %*% b_lds[1,])
  bootrec5 <- (b_lds[1,] * b_lds[1,]) %*% seboot
  bootrec6 <- sqrt((b_lds[1,] * b_lds[1,]) %*% (varboot))
  sebootrec <- sd(betas %*% b_lds[1,])
  ### bingo
  sqrt(b_lds[1,] %*% cov(betas) %*% b_lds[1,])
  sqrt(rowSums((b_lds %*% (cov(betas)) * b_lds)))
  ### bingo
  Qr <- qr.lm(fits[[1]])
  covu <- chol2inv(Qr$qr)
  
  seboot <- sapply(boot_fits, function(x) sd(x$t))
  
  R <- do.call(cbind, lapply(fits, resid))
  covR <- cov(R) * (nrow(R) - 1)/lm.orig$df.residual
  #covR <- cov(R) 
  
  bootrec1 <- sqrt(b_lds[1,] %*% (covR * covu[1,1]) %*% b_lds[1,])
  
  bootrec2 <- (b_lds[1,] * b_lds[1,]) %*% seboot
  

  bdist <- do.call(cbind, lapply(boot_fits, function(x) x$t))
  sebootrec <- sd(betas %*% b_lds[1,])
  
}



