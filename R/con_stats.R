fit_Ftests <- function(object) {
  w <- object$weights
  ssr <- if (is.null(w)) {
    apply(object$residuals, 2, function(vals) sum(vals^2))
  } else {
    apply(object$residuals, 2, function(vals) sum((vals^2 *w)))
  }
  
  mss <- if (is.null(w)) {
    apply(object$fitted.values, 2, function(vals) sum(vals^2))
  } else {
    apply(object$fitted.values, 2, function(vals) sum(vals^2 * w))
  }
  
  #if (ssr < 1e-10 * mss) 
  #  warning("ANOVA F-tests on an essentially perfect fit are unreliable")
  
  dfr <- df.residual(object)
  p <- object$rank
  
  p1 <- 1L:p
  #comp <- object$effects[p1]
  asgn <- object$assign[stats:::qr.lm(object)$pivot][p1]
  
  
  #nmeffects <- c("(Intercept)", attr(object$terms, "term.labels"))
  #tlabels <- nmeffects[1 + unique(asgn)]
  
  df <- c(lengths(split(asgn, asgn)), dfr)
  
  I.p <- diag(nrow(coefficients(object)))
  nterms <- length(unique(asgn))
  hmat <- lapply(1:nterms, function(i) {
    subs <- which(asgn == i)
    hyp.matrix <- I.p[subs, , drop=FALSE]
    hyp.matrix[!apply(hyp.matrix, 1, function(x) all(x == 0)), , drop = FALSE]
  })
  
  ret <- lapply(seq_along(ssr), function(i) {
    comp <- object$effects[,i]
    ss <- c(unlist(lapply(split(comp^2, asgn), sum)), ssr[i])
    ms <- ss/df
    f <- ms/(ssr[i]/dfr)
    
    P <- pf(f, df, dfr, lower.tail = FALSE)
    list(F=f, P=p)
  })
  
  FMat <- do.call(rbind, lapply(ret, "[[", "F"))
  PMat <- do.call(rbind, lapply(ret, "[[", "P"))
  
  list(F=FMat, P=PMat)
  
  
}

beta_stats <- function(lmfit) {
  Qr <- stats:::qr.lm(lmfit)
  cov.unscaled <- chol2inv(Qr$qr)
  betamat <- lmfit$coefficients
  
  
  rss <- colSums(lmfit$residuals^2)
  rdf <- lmfit$df.residual
  resvar <- rss/rdf
  sigma <- sqrt(resvar)
  
  vc <- sapply(1:ncol(betamat), function(i) {
    vcv <- cov.unscaled * sigma[i]^2
    sqrt(diag(vcv))
  })
  
  prob <- 2 * (1 - pt(abs(betamat/vc), lmfit$df.residual))
  tstat <- betamat/vc
  return(
    list(
      estimate=function() betamat,
      stat=function() betamat/vc,
      se=function() vc,
      prob=function() 2 * (1 - pt(abs(betamat/vc), lmfit$df.residual)),
      stat_type="tstat")
  )
  
}


fit_Fcontrasts <- function(lmfit, conmat, colind) {
  Qr <- stats:::qr.lm(lmfit)
  cov.unscaled <- chol2inv(Qr$qr[colind,colind,drop=FALSE])
  betamat <- lmfit$coefficients[colind,]
  
  df <- lmfit$df.residual
  
  rss <- colSums(lmfit$residuals^2)
  rdf <- lmfit$df.residual
  resvar <- rss/rdf
  sigma <- sqrt(resvar)
  sigma2 <- sigma^2
  
  #msr <- summary.lm(reg)$sigma  # == SSE / (n-p)
  r <- nrow(conmat)
  
  
  #                        -1
  #     (Cb - d)' ( C (X'X)   C' ) (Cb - d) / r
  # F = ---------------------------------------
  #                 SSE / (n-p)
  #
  
  cm <- solve((conmat %*% cov.unscaled %*% t(conmat)))
  
  Fstat <- sapply(1:ncol(betamat), function(i) {
    b <- betamat[,i]
    cb <- conmat %*% b
    Fstat <- t(cb) %*% cm %*% (cb) / r / sigma2[i]
  })
  
  
  return(
    list(
      estimate=function() betamat,
      se=function() sigma2,
      stat=function() Fstat,
      prob=function() 1-pf(Fstat,r,rdf),
      stat_type="Fstat")
  )
  
}

#' fit_contrasts
#' 
#' @param lmfit the \code{lm} object
#' @param conmat the contrast \code{matrix} or contrast \code{vector}
#' @param colind the subset column indices in the design associated with the contrast. 
#' @export
fit_contrasts <- function(lmfit, conmat, colind) {
  if (!is.matrix(conmat) && is.numeric(conmat)) {
    conmat <- matrix(conmat, 1, length(conmat))
  }
  
  Qr <- stats:::qr.lm(lmfit)
  
  p1 <- 1:lmfit$rank
  cov.unscaled <- chol2inv(Qr$qr[colind,colind,drop=FALSE])
  betamat <- lmfit$coefficients[colind,]
  ct <- conmat %*% betamat
  
  rss <- colSums(lmfit$residuals^2)
  rdf <- lmfit$df.residual
  resvar <- rss/rdf
  sigma <- sqrt(resvar)
  
  vc <- sapply(1:ncol(betamat), function(i) {
    vcv <- cov.unscaled * sigma[i]^2
    sqrt(diag(conmat %*% vcv %*% t(conmat)))
  })
  
  
  return(
    list(
      estimate=function() ct,
      se=function() vc,
      stat=function() ct/vc,
      prob=function() 2 * (1 - pt(abs(ct/vc), lmfit$df.residual)),
      stat_type="tstat")
  )
  
}