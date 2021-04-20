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
# 
# ### from vcov package -- suppodely faster
# Vcov.lm = function(object, ...) {
#   if (p <- object$rank) {
#     p1 = seq_len(p)
#     rss = if (is.null(w <- object$weights)) { 
#       sum(object$residuals^2)
#     } else {
#       sum(w * object$residuals^2)
#     }
#     covmat = rss * chol2inv(object$qr$qr[p1, p1, drop = FALSE])/
#       object$df.residual
#     nm = names(object$coefficients)
#     dimnames(covmat) = list(nm, nm)
#     return(covmat)
#   } else return(numeric(0))
# }
# 
# Vcov.glm = function(object, dispersion = NULL, ...) {
#   if (p <- object$rank) {
#     if (is.null(dispersion)) {
#       dispersion = if (object$family$family %in% c('poisson', 'binomial')) {
#         1
#       } else {
#         df_r = object$df.residual
#         if (df_r) {
#           if (any(!object$weights)) 
#             warning('observations with zero weight not',
#                     'used for calculating dispersion')
#           w = object$weights
#           idx = w > 0
#           sum(w[idx] * object$residuals[idx]^2)/df_r
#         } else NaN
#       }
#     }
#     p1 = seq_len(p)
#     nm <- names(object$coefficients[object$qr$pivot[p1]])
#     covmat = dispersion * chol2inv(object$qr$qr[p1, p1, drop = FALSE])
#     dimnames(covmat) = list(nm, nm)
#     return(covmat)
#   } else return(numeric(0))
# }

beta_stats <- function(lmfit, varnames) {
  Qr <- stats:::qr.lm(lmfit)
  cov.unscaled <- chol2inv(Qr$qr)
  
  cfs <- coef(lmfit)
  
  betamat <- if (is.vector(cfs)) {
    as.matrix(coef(lmfit))
  } else {
    betamat <- cfs
  }
 
  
  rss <- colSums(as.matrix(lmfit$residuals^2))
  
  rdf <- lmfit$df.residual
 
  resvar <- rss/rdf
  sigma <- sqrt(resvar)

  vc <- sapply(1:ncol(betamat), function(i) {
    vcv <- cov.unscaled * sigma[i]^2
    sqrt(diag(vcv))
  })
  
  #sqrt(diag(vcov(lmfit))
  
  prob <- 2 * (1 - pt(abs(betamat/vc), lmfit$df.residual))
  tstat <- betamat/vc
  
  betamat <- t(betamat)
  vc <- t(vc)
  
  colnames(betamat) <- varnames
  colnames(vc) <- varnames
  
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
  #browser()
  Qr <- stats:::qr.lm(lmfit)
  
  ## subset covariance matrix by column indices of the contrast matrix associated with term in questions
  #cov.unscaled <- try(chol2inv(Qr$qr[colind,colind,drop=FALSE]))
  cov.unscaled <- try(chol2inv(Qr$qr))
  
  cmat <- matrix(0, nrow(conmat), ncol(cov.unscaled))
  cmat[,colind] <- conmat
  
  cfs <- coef(lmfit)
  
  betamat <- if (is.vector(cfs)) {
    as.matrix(coef(lmfit))
  } else {
    cfs
  }
  
  #betamat <- lmfit$coefficients[colind,]
  
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
  
  ## conmat is same dimensions as the number of levels contained in the term
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
#' @importFrom purrr map_dbl
fit_contrasts <- function(lmfit, conmat, colind) {
  
  if (!is.matrix(conmat) && is.numeric(conmat)) {
    conmat <- t(matrix(conmat, 1, length(conmat)))
  }
 
  
  Qr <- stats:::qr.lm(lmfit)
  
  p1 <- 1:lmfit$rank
  #cov.unscaled <- try(chol2inv(Qr$qr[colind,colind,drop=FALSE]))
  cov.unscaled <- try(chol2inv(Qr$qr))
  cmat <- matrix(0, nrow(cov.unscaled), 1)
  cmat[colind,] <- conmat
  
  
  if (inherits(cov.unscaled, "try-error")) {
    stop("fit_contrasts: error computing contrast covariance")
    #browser()
  }
 
  #print(length(colind))
  
  cfs <- coef(lmfit)
  if (is.vector(cfs)) {
    #as.matrix(coef(lmfit)[colind])
    betamat <- as.matrix(coef(lmfit))
  } else {
    #betamat <- cfs[colind,,drop=FALSE]
    #}
    betamat <- cfs
  }
  
  if (inherits(cov.unscaled, "try-error")) {
    stop("fit_contrasts: error computing contrast covariance")
    #browser()
  }
  
  ct <- as.vector(t(cmat) %*% betamat)
  
  rss <- colSums(as.matrix(lmfit$residuals^2))
  rdf <- lmfit$df.residual
  resvar <- rss/rdf
  sigma <- sqrt(resvar)
  

  vc <- map_dbl(1:ncol(betamat), function(i) {
    vcv <- cov.unscaled * sigma[i]^2
    sqrt(diag(t(cmat) %*% vcv %*% cmat))
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