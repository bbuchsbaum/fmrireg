#' @keywords internal
#' @noRd
qr.lm <- getFromNamespace("qr.lm", "stats")

#' @keywords internal
#' @noRd
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



#' Beta Statistics for Linear Model
#'
#' @description
#' This function calculates beta statistics for a linear model fit.
#'
#' @param lmfit Fitted linear model object.
#' @param varnames Vector of variable names.
#' @param se Logical flag indicating whether to calculate standard errors (default: TRUE).
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{\code{estimate}: Estimated beta coefficients.}
#'   \item{\code{stat}: t-statistics of the beta coefficients (if \code{se = TRUE}).}
#'   \item{\code{se}: Standard errors of the beta coefficients (if \code{se = TRUE}).}
#'   \item{\code{prob}: Probabilities associated with the t-statistics (if \code{se = TRUE}).}
#'   \item{\code{stat_type}: Type of the statistics calculated.}
#' }
#' @examples
#' data(mtcars)
#' lm_fit <- lm(mpg ~ wt + qsec + am, data = mtcars)
#' beta_stats(lm_fit, c("Intercept", "wt", "qsec", "am"))
#' @noRd
#' @autoglobal
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
    
    
    ## assumes se computed
    ret <- dplyr::tibble(
      type="beta",
      name="parameter_estimates",
      stat_type="tstat",
      df.residual=rdf,
      conmat=list(NULL),
      colind=list(NULL),
      data=list(dplyr::tibble(
        estimate=list(betamat),
        se=list(vc),
        stat=list(betamat/vc),
        prob=list(2 * (1 - pt(abs(betamat/vc), lmfit$df.residual))),
        sigma=list(sigma)))
    )
  } else {
    stop("unimplemented")
    betamat <- t(betamat)
    colnames(betamat) <- varnames
    ret <- dplyr::tibble(
      type="beta",
      name="parameter_estimates",
      stat_type="effect",
      df.residual=rdf,
      conmat=list(NULL),
      colind=list(NULL),
      data=list(tibble(
        estimate=list(betamat),
        se=list(NULL),
        stat=list(NULL),
        prob=list(NULL),
        sigma=list(NULL)))
    )
  }
}

#' Fit F-contrasts for Linear Model
#'
#' @description
#' This function calculates F-contrasts for a fitted linear model.
#'
#' @param lmfit Fitted linear model object.
#' @param conmat Contrast matrix.
#' @param colind Column indices corresponding to the variables in the contrast matrix.
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{\code{estimate}: Estimated contrasts.}
#'   \item{\code{se}: Residual variance.}
#'   \item{\code{stat}: F-statistics for the contrasts.}
#'   \item{\code{prob}: Probabilities associated with the F-statistics.}
#'   \item{\code{stat_type}: Type of the statistics calculated.}
#' }
#' @noRd
#' @keywords internal
fit_Fcontrasts <- function(lmfit, conmat, colind) {
  Qr <- qr.lm(lmfit)

  cov.unscaled <- try(chol2inv(Qr$qr))
  
  cmat <- matrix(0, ncol(conmat), ncol(cov.unscaled))
  cmat[,colind] <- t(conmat)
  
  cfs <- coef(lmfit)
  
  betamat <- if (is.vector(cfs)) {
    as.matrix(coef(lmfit))
  } else {
    cfs
  }
  
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
  
  estimate <- purrr::map_dbl(1:ncol(betamat), function(i) {
    b <- betamat[,i]
    cb <- cmat %*% b
    t(cb) %*% cm %*% (cb) / r 
  })
  
  se <- sigma2
  Fstat <- estimate/se
  

  return(
    structure(list(
      conmat=conmat,
      estimate=estimate,
      se=sigma2,
      df.residual=rdf,
      stat=Fstat,
      prob=1-pf(Fstat,r,rdf),
      stat_type="Fstat")
      ,class=c("Fstat", "result_stat"))
  )

}


#' @keywords internal
#' @noRd
#' @autoglobal
estimate_contrast.contrast <- function(x, fit, colind, ...) {
  ret <- fit_contrasts(fit, x$weights, colind, se=TRUE)
  
  ## assumes se computed
  tibble(
    type="contrast",
    name=x$name,
    stat_type=ret$stat_type,
    df.residual=ret$df.residual,
    conmat=list(x$weights),
    colind=list(colind),
    data=list(tibble(
      estimate=ret$estimate,
      se=ret$se,
      stat=ret$stat,
      prob=ret$prob,
      sigma=ret$sigma)))
    
}


#' @export
#' @autoglobal
#' @keywords internal
#' @noRd
estimate_contrast.Fcontrast <- function(x, fit, colind, ...) {
  ret <- fit_Fcontrasts(fit, x$weights, colind)
  tibble(
    type="Fcontrast",
    name=x$name,
    stat_type=ret$stat_type,
    df.residual=ret$df.residual,
    conmat=list(x$weights),
    colind=list(colind),
    data=list(tibble(
      estimate=ret$estimate,
      se=ret$se,
      stat=ret$stat,
      prob=ret$prob)))
}

#' Fit Contrasts for Linear Model
#'
#' @description
#' This function calculates contrasts for a fitted linear model.
#'
#' @param lmfit The fitted linear model object.
#' @param conmat The contrast matrix or contrast vector.
#' @param colind The subset column indices in the design associated with the contrast.
#' @param se Whether to compute standard errors, t-statistics, and p-values (default: TRUE).
#'
#' @return A list containing the following elements:
#' \itemize{
#'   \item{\code{conmat}: Contrast matrix.}
#'   \item{\code{sigma}: Residual standard error.}
#'   \item{\code{df.residual}: Degrees of freedom for residuals.}
#'   \item{\code{estimate}: Estimated contrasts.}
#'   \item{\code{se}: Standard errors of the contrasts (if \code{se = TRUE}).}
#'   \item{\code{stat}: t-statistics for the contrasts (if \code{se = TRUE}).}
#'   \item{\code{prob}: Probabilities associated with the t-statistics (if \code{se = TRUE}).}
#'   \item{\code{stat_type}: Type of the statistics calculated.}
#' }
#' @noRd
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
  
    vc <- purrr::map_dbl(1:ncol(betamat), function(i) {
      vcv <- cov.unscaled * sigma[i]^2
      sqrt(diag(t(cmat) %*% vcv %*% cmat))
    })
  
    return(
      structure(list(
        conmat=cmat,
        sigma=sigma,
        df.residual=lmfit$df.residual,
        estimate=ct,
        se=vc,
        stat=ct/vc,
        prob=2 * (1 - pt(abs(ct/vc), lmfit$df.residual)),
        stat_type="tstat"),
        class=c("tstat", "result_stat"))
    )
  } else {
    return(
      structure(list(
        conmat=cmat,
        sigma=sigma,
        df.residual=lmfit$df.residual,
        estimate=ct,
        stat_type="effects"),
        class=c("effect", "result_stat"))
    )
  }
}

