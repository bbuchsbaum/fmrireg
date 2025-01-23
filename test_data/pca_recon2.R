library(mvtnorm)
library(multivarious)

X <- matrix(rnorm(100 * 10), 100, 10)
X <- scale(X, center=TRUE, scale=FALSE)
pcres <- pca(X, ncomp = 10, preproc=pass())


y0 <- X[,1]

## design matrix
Xd <- rmvnorm(nrow(X), sigma = diag(seq(10, 1)))
Xd <- scale(Xd, center=TRUE)

lm.orig <- lm(y0 ~ Xd-1)


### lm in pca space
fits <- lapply(1:ndvar, function(i) {
  y <- scores(pcres)[,i]
  lm(y ~ Xd -1)
})

## reconstructed betas
b_lds <- coef(pcres)
brec <- (t(bvals) %*% t(b_lds))[1]
abs(brec[1] - coef(lm.orig)[1]) < .1

### sparse PCA
run <- function() {
  
  ## 50 voxels, 30 bases
  basis = matrix(rnorm(50 * 30), 50, 30)
  
  ## 50 voxels, 15 time points
  X <- matrix(rnorm(50 * 15), 50, 15)
  
  basis <- scale(basis, scale = FALSE)
  X <- scale(X, scale = FALSE)
  
  reg <- multivarious::regress(basis, X, method="enet", alpha=0, lambda=0)
  Xrecon <- reconstruct(reg)
  
  Xrt <- t(Xrecon)
  Xrt2 = reg$coefficients %*%  t(reg$s)
  
  Xd <- rmvnorm(ncol(X), sigma = diag(seq(5, 1)))
  
  ## time series of 1st voxel
  y0 <- Xrt[,1]
  y0 <- X[1,]
  lm.recon <- lm(y0 ~ Xd - 1)
  
  ### lm in pca space
  fits <- lapply(1:ncol(basis), function(i) {
    y <- reg$coefficients[, i]
    lm(y ~ Xd - 1)
  })
  
  #b_lds <- inverse_projection(reg)
  #b_lds <- t(coef(reg))
  bvals <- sapply(fits, function(fit)
    coef(fit)[1])
  
  brec <- (t(bvals) %*% t(scores(reg)))
  c(brec[1,1], coef(lm.recon)[1])
}




