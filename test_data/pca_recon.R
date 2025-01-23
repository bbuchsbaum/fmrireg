


library(multivarious)
library(mvtnorm)


runit <- function(nvar=20, ndvar=5) {
  X <- matrix(rnorm(nvar * 1000), 1000, nvar)
  pcres <- pca(X, ncomp = nvar)
  recon <- reconstruct(pcres, comp = 1:nvar)
  recon <- scale(recon, center = TRUE, scale = FALSE)
  
  y0 <- recon[, 1]
  
  Xd <- rmvnorm(nrow(X), sigma = diag(seq(ndvar, 1)))
  Xd <- scale(Xd, center=TRUE)
  #Xd <- matrix(rnorm(nrow(X) * ndvar), nrow(X), ndvar)

  lm.orig <- lm(y0 ~ Xd-1)
  
  fits <- lapply(1:nvar, function(i) {
    y <- multivarious::scores(pcres)[,i]
    lm(y ~ Xd -1)
  })
  
  R <- do.call(cbind, lapply(fits, fitted))
  
  se_lm <- function(fit) sqrt(diag(vcov(fit)))
  
  se_lm2 <- function(fit, dofpen=0) {
    Qr <- qr.lm(fit)
    cov.unscaled <- chol2inv(Qr$qr)
    rss <- colSums(as.matrix(fit$residuals^2))
    rdf <- fit$df.residual + dofpen
    
    resvar <- rss/rdf
    resvar <- rss
    sigma <- sqrt(resvar)
    
    vc <- sapply(length(fit$coefficients), function(i) {
      vcv <- cov.unscaled * sigma[i]^2
      #vcv <- cov.unscaled * sigma^2
      sqrt(diag(vcv))
    })
    list(se=vc, sigma=sigma, covu=cov.unscaled)
  }
  
  b_lds <- coef(pcres)[, 1:nvar,drop=FALSE]
  bvals <- sapply(fits, function(fit) coef(fit)[1])
 
  brec <- (t(as.matrix(bvals)) %*% t(b_lds))[1]
  
  n = nrow(X)
  #serrs <- sapply(fits, function(fit) c(se_lm2(fit,length(fits) + 1/length(fits))[1]))
  serrs <- sapply(1:length(fits), function(i) c(se_lm(fits[[i]])))
  scovs <- lapply(fits, function(fit) vcov(fit))
  
  #dr <- lm.orig$df.residual
  #fac=dr/(dr - 1/2)
  #fac=1
  #serrs2 = sqrt((serrs * n/(n-(2)))^2)
  #serrs2 = sqrt((serrs ^ 2))
  
  #brec <- t(bvals) %*% t(b_lds)
  
  #fac=sqrt(9 + .013)/(sqrt(9))
  
  
  whitened_se <- function(fits, whiten=TRUE) {
    #Qr <- stats:::qr.lm(lm.orig)
    Qr <- qr.lm(fits[[1]])
    cov.unscaled <- chol2inv(Qr$qr)
    
    tmp = do.call(cbind, lapply(fits, resid))
    if (whiten) {
      w=whiteningMatrix(cor(tmp), method="Cholesky")
      tmp <- tmp %*% w
    }
    
    rss <- colSums(as.matrix(tmp^2))
    rdf <- lm.orig$df.residual 
    resvar <- rss/rdf
    sigma <- sqrt(resvar)
    
    vc <- sapply(1:ncol(tmp), function(i) {
      vcv <- cov.unscaled * sigma[i]^2
      #vcv <- cov.unscaled * sigma^2
      sqrt(diag(vcv))
    })
    ##
    ##list(se=se, cov=cov(tmp) * (nrow(tmp) - 1)/lm.orig$df.residual, covu=cov.unscaled)
    list(se=1, cov=cov(tmp) * (nrow(tmp) - 1)/lm.orig$df.residual, covu=cov.unscaled)
  }
  
  se_no  <- whitened_se(fits, FALSE)
  ##serec_w <- sqrt(b_lds[1, ] %*% (diag(se_w^2)) %*% b_lds[1,])
  ##serec <- sqrt(b_lds[1, ] %*% (diag(serrs^2)) %*% b_lds[1,])[1]
  
  ###### CORRECT RECONSTRUCTION of SE
  for (i in 1:5) {
    b_lds[1,] %*% (se_no[[2]] * se_no[[3]][i,i]) %*% b_lds[1,]
  }
  serec2 <- sqrt(b_lds[1, ] %*% (se_no[[2]] * se_no[[3]][2,2]) %*% b_lds[1,])
  
  ###################################
  
  serr <- se_lm(lm.orig)[1]
  c(serr, serec, serec2)
}





