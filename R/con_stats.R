#' @keywords internal
qr.lm <- getFromNamespace("qr.lm", "stats")

#' @keywords internal
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
  asgn <- object$assign[qr.lm(object)$pivot][p1]
  
  
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


old_beta_stats <- function(lmfit, varnames, se=TRUE) {
  cfs <- coef(lmfit)
  betamat <- if (is.vector(cfs)) {
    as.matrix(coef(lmfit))
  } else {
    betamat <- cfs
  }

  if (se) {
    Qr <- qr.lm(lmfit)
    cov.unscaled <- chol2inv(Qr$qr)
    rss <- colSums(as.matrix(lmfit$residuals^2))
  
    rdf <- lmfit$df.residual
    resvar <- rss/rdf
    sigma <- sqrt(resvar)

    vc <- sapply(1:ncol(betamat), function(i) {
      vcv <- cov.unscaled * sigma[i]^2
      sqrt(diag(vcv))
    })
  
  
    prob <- 2 * (1 - pt(abs(betamat/vc), lmfit$df.residual))
    tstat <- betamat/vc
  
    vc <- t(vc)
    colnames(vc) <- varnames
    
    betamat <- t(betamat)
    vc <- t(vc)
  
  
    list(
        estimate=function() betamat,
        stat=function() betamat/vc,
        se=function() vc,
        prob=function() 2 * (1 - pt(abs(betamat/vc), lmfit$df.residual)),
        stat_type="tstat"
    )
  } else {
    betamat <- t(betamat)
    colnames(betamat) <- varnames
    list(
      estimate=function() betamat,
      stat_type="effects"
    )
  }
}

beta_stats <- function(lmfit, varnames, se=TRUE) {
  cfs <- coef(lmfit)
  
  betamat <- if (is.vector(cfs)) {
    as.matrix(coef(lmfit))
  } else {
    betamat <- cfs
  }
  
  if (se) {
    Qr <- qr.lm(lmfit)
    cov.unscaled <- chol2inv(Qr$qr)
    rss <- colSums(as.matrix(lmfit$residuals^2))
    
    rdf <- lmfit$df.residual
    resvar <- rss/rdf
    sigma <- sqrt(resvar)
    
    vc <- sapply(1:ncol(betamat), function(i) {
      vcv <- cov.unscaled * sigma[i]^2
      sqrt(diag(vcv))
    })
    
    
    #prob <- 2 * (1 - pt(abs(betamat/vc), lmfit$df.residual))
    #tstat <- betamat/vc
    
    vc <- t(vc)
    colnames(vc) <- varnames
    
    betamat <- t(betamat)
    colnames(betamat) <- varnames
    #vc <- t(vc)
  
    list(
      estimate=betamat,
      stat=betamat/vc,
      se=vc,
      prob=2 * (1 - pt(abs(betamat/vc), lmfit$df.residual)),
      stat_type="tstat"
    )
  } else {
    betamat <- t(betamat)
    colnames(betamat) <- varnames
    list(
      estimate=betamat,
      stat_type="effects"
    )
  }
}


fit_Fcontrasts <- function(lmfit, conmat, colind) {
  #browser()
  Qr <- qr.lm(lmfit)

  cov.unscaled <- try(chol2inv(Qr$qr))
  
  cmat <- matrix(0, nrow(conmat), ncol(cov.unscaled))
  cmat[,colind] <- conmat
  
  cfs <- coef(lmfit)
  
  betamat <- if (is.vector(cfs)) {
    as.matrix(coef(lmfit))
  } else {
    cfs
  }
  
  df <- lmfit$df.residual
  
  rss <- colSums(as.matrix(lmfit$residuals^2))
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
  
  cm <- solve((cmat %*% cov.unscaled %*% t(cmat)))
  
  #Fstat <- map_dbl(1:ncol(betamat), function(i) {
  #  b <- betamat[,i]
  #  cb <- cmat %*% b
  #  Fstat <- t(cb) %*% cm %*% (cb) / r / sigma2[i]
  #})
  
  estimate <- map_dbl(1:ncol(betamat), function(i) {
    b <- betamat[,i]
    cb <- cmat %*% b
    t(cb) %*% cm %*% (cb) / r 
  })
  
  se <- sigma2
  Fstat <- estimate/se
  
  return(
    list(
      estimate=estimate,
      se=sigma2,
      stat=Fstat,
      prob=1-pf(Fstat,r,rdf),
      stat_type="Fstat")
  )
  
}


old_fit_Fcontrasts <- function(lmfit, conmat, colind) {
  #browser()
  Qr <- qr.lm(lmfit)
  
  cov.unscaled <- try(chol2inv(Qr$qr))
  
  cmat <- matrix(0, nrow(conmat), ncol(cov.unscaled))
  cmat[,colind] <- conmat
  
  cfs <- coef(lmfit)
  
  betamat <- if (is.vector(cfs)) {
    as.matrix(coef(lmfit))
  } else {
    cfs
  }
  
  df <- lmfit$df.residual
  
  rss <- colSums(as.matrix(lmfit$residuals^2))
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
  
  cm <- solve((cmat %*% cov.unscaled %*% t(cmat)))
  
  Fstat <- map_dbl(1:ncol(betamat), function(i) {
    b <- betamat[,i]
    cb <- cmat %*% b
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
#' @param se whetheer to compute standard errors, t-statistics, and p-values
#' @importFrom purrr map_dbl
fit_contrasts <- function(lmfit, conmat, colind, se=TRUE) {
  if (!is.matrix(conmat) && is.numeric(conmat)) {
    conmat <- t(matrix(conmat, 1, length(conmat)))
  }
  
  cfs <- coef(lmfit)
  
  if (is.vector(cfs)) {
    betamat <- as.matrix(coef(lmfit))
  } else {
    betamat <- cfs
  }
  
  ncoef <- nrow(betamat)
  cmat <- matrix(0, ncoef, 1)
  cmat[colind,] <- conmat
  
  ct <- as.vector(t(cmat) %*% betamat)
  rss <- colSums(as.matrix(lmfit$residuals^2))
  rdf <- lmfit$df.residual
  resvar <- rss/rdf
  sigma <- sqrt(resvar)
  
  p1 <- 1:lmfit$rank
  
  if (se) {
    Qr <- qr.lm(lmfit)
    cov.unscaled <- try(chol2inv(Qr$qr))
  
    if (inherits(cov.unscaled, "try-error")) {
      stop("fit_contrasts: error computing contrast covariance")
    }
  
    vc <- map_dbl(1:ncol(betamat), function(i) {
      vcv <- cov.unscaled * sigma[i]^2
      sqrt(diag(t(cmat) %*% vcv %*% cmat))
    })
  
  
    return(
      list(
        conmat=cmat,
        sigma=sigma,
        df.residual=lmfit$df.residual,
        estimate=ct,
        se=vc,
        stat=ct/vc,
        prob=2 * (1 - pt(abs(ct/vc), lmfit$df.residual)),
        stat_type="tstat")
    )
  } else {
    return(
      list(
        conmat=cmat,
        sigma=sigma,
        df.residual=lmfit$df.residual,
        estimate=ct,
        stat_type="effects")
    )
  }
}


old_fit_contrasts <- function(lmfit, conmat, colind, se=TRUE) {
  if (!is.matrix(conmat) && is.numeric(conmat)) {
    conmat <- t(matrix(conmat, 1, length(conmat)))
  }
  
  ncoef <- nrow(coef(lmfit))
  cmat <- matrix(0, ncoef, 1)
  cmat[colind,] <- conmat
  
  cfs <- coef(lmfit)
  
  if (is.vector(cfs)) {
    betamat <- as.matrix(coef(lmfit))
  } else {
    betamat <- cfs
  }
  
  ct <- as.vector(t(cmat) %*% betamat)
  rss <- colSums(as.matrix(lmfit$residuals^2))
  rdf <- lmfit$df.residual
  resvar <- rss/rdf
  sigma <- sqrt(resvar)
  
  p1 <- 1:lmfit$rank
  
  if (se) {
    Qr <- qr.lm(lmfit)
    cov.unscaled <- try(chol2inv(Qr$qr))
    
    if (inherits(cov.unscaled, "try-error")) {
      stop("fit_contrasts: error computing contrast covariance")
    }
    
    vc <- map_dbl(1:ncol(betamat), function(i) {
      vcv <- cov.unscaled * sigma[i]^2
      sqrt(diag(t(cmat) %*% vcv %*% cmat))
    })
    
    
    return(
      list(
        conmat=cmat,
        sigma=sigma,
        df.residual=lmfit$df.residual,
        estimate=function() ct,
        se=function() vc,
        stat=function() ct/vc,
        prob=function() 2 * (1 - pt(abs(ct/vc), lmfit$df.residual)),
        stat_type="tstat")
    )
  } else {
    return(
      list(
        conmat=cmat,
        sigma=sigma,
        df.residual=lmfit$df.residual,
        estimate=function() ct,
        stat_type="effects")
    )
  }
}