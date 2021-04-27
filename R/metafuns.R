

#' @keywords internal
meta_stouffer <- function(pval, se) {
  if (any(pval == 0)) {
    pval[pval == 0] <- .Machine$double.eps
  }
  inv_var <- 1/(se^2)
  wts <- inv_var/rowSums(inv_var)
  zscore <- qnorm(1-pval)
  wzscore <- zscore * wts
  wzscore <- rowSums(wzscore)
  
  return(
    list(
      estimate=wzscore,
      se=sqrt(rowSums(wts*wts*(se^2))),
      stat=wzscore,
      prob=1-pnorm(wzscore),
      stat_type="zfstat"
    )
  )
}


#' @keywords internal
meta_fixef <- function(beta,se, weighting=c("inv_var", "equal")) {
  weighting <- match.arg(weighting)
  #browser()
  if (weighting == "inv_var") {
    inv_var <- 1/(se^2)
    wts <- inv_var/rowSums(inv_var)
    wbeta <- beta * wts
    wbeta <- rowSums(wbeta)
    pooledse <- sqrt(rowSums(wts*wts*(se^2)))
  } else {
    wbeta <- rowSums(wbeta)
    pooledse <- sqrt(rowSums((se^2)))
  }
  
  return(
    list(
      estimate=wbeta,
      se==pooledse,
      stat=wbeta/pooledse,
      prob=1-pchisq((wbeta/pooledse)^2,1),
      stat_type="zstat")
  )
  
}


#' @keywords internal
meta_Fcontrasts <- function(fres) {
  
  ncon <- length(fres[[1]])
  res <- lapply(seq(1,ncon), function(i) {
    pval <- do.call(cbind, lapply(fres, function(x) as.vector(x[[i]]$prob)))
    se <- do.call(cbind, lapply(fres, function(x) x[[i]]$se))
    
    meta_stouffer(pval,se)
  })
  names(res) <- names(fres)
  res
}


#' @keywords internal
meta_contrasts <- function(cres) {
  ncon <- length(cres[[1]])
  if (ncon > 0) {
    res <- lapply(1:ncon, function(i) {
      beta <- do.call(cbind, lapply(cres, function(x) as.vector(x[[i]]$estimate)))
      se <- do.call(cbind, lapply(cres, function(x) x[[i]]$se))
      meta_fixef(beta,se)
    })
  } else {
    stop("there are no contrasts for this model.")
  }
  
  return(
    list(
      estimate=do.call(cbind, lapply(res, function(x) x$estimate)),
      se=do.call(cbind, lapply(res, function(x) x$se)),
      stat=do.call(cbind, lapply(res, function(x) x$stat)),
      prob=do.call(cbind, lapply(res, function(x) x$prob)),
      stat_type="zstat"
    )
  )
}


#' @keywords internal
meta_betas <- function(bstats, colind) {
  
  len <- length(colind)
  
  res <- lapply(colind, function(i) {
    #print(i)
    beta <- do.call(cbind, lapply(bstats, function(x) x$estimate[,i]))
    se <- do.call(cbind, lapply(bstats, function(x) x$se[,i]))
    meta_fixef(beta,se)
  })
  
  return(
    list(
      estimate=do.call(cbind, lapply(res, function(x) x$estimate)),
      se=do.call(cbind, lapply(res, function(x) x$se)),
      stat=do.call(cbind, lapply(res, function(x) x$stat)),
      prob=do.call(cbind, lapply(res, function(x) x$prob)),
      stat_type="zstat"
    )
  )
}