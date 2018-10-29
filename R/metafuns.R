meta_stouffer <- function(pval, se) {
  inv_var <- 1/(se^2)
  wts <- inv_var/rowSums(inv_var)
  zscore <- qnorm(1-pval)
  wzscore <- zscore * wts
  wzscore <- rowSums(wzscore)
  
  return(
    list(
      estimate=function() wzscore,
      se=function() sqrt(rowSums(wts*wts*(se^2))),
      stat=function() wzscore,
      prob=function() 1-pnorm(wzscore),
      stat_type="zfstat"
    )
  )
}

meta_fixef <- function(beta,se) {
  inv_var <- 1/(se^2)
  
  wts <- inv_var/rowSums(inv_var)
  wbeta <- beta * wts
  wbeta <- rowSums(wbeta)
  pooledse <- sqrt(rowSums(wts*wts*(se^2)))
  
  return(
    list(
      estimate=function() wbeta,
      se=function() pooledse,
      stat=function() wbeta/pooledse,
      prob=function() 1-pchisq((wbeta/pooledse)^2,1),
      stat_type="zstat")
  )
  
}

meta_Fcontrasts <- function(fres) {
  ncon <- length(fres[[1]])
  res <- lapply(seq(1,ncon), function(i) {
    pval <- do.call(cbind, lapply(fres, function(x) as.vector(x[[i]]$prob())))
    se <- do.call(cbind, lapply(fres, function(x) x[[i]]$se()))
    
    meta_stouffer(pval,se)
  })
  names(res) <- names(fres)
  res
}

meta_contrasts <- function(cres) {
  ncon <- length(cres[[1]])
  if (ncon > 0) {
    res <- lapply(1:ncon, function(i) {
      beta <- do.call(cbind, lapply(cres, function(x) as.vector(x[[i]]$estimate())))
      se <- do.call(cbind, lapply(cres, function(x) x[[i]]$se()))
      meta_fixef(beta,se)
    })
  } else {
    stop("there are no contrasts for this model.")
  }
  
  return(
    list(
      estimate=function() do.call(cbind, lapply(res, function(x) x$estimate())),
      se=function() do.call(cbind, lapply(res, function(x) x$se())),
      stat=function() do.call(cbind, lapply(res, function(x) x$stat())),
      prob=function() do.call(cbind, lapply(res, function(x) x$prob())),
      stat_type="zstat"
    )
  )
}

meta_betas <- function(bstats, colind) {
  
  len <- length(colind)
  
  res <- lapply(colind, function(i) {
    print(i)
    beta <- do.call(cbind, lapply(bstats, function(x) x$estimate()[,i]))
    se <- do.call(cbind, lapply(bstats, function(x) x$se()[,i]))
    meta_fixef(beta,se)
  })
  
  return(
    list(
      estimate=function() do.call(cbind, lapply(res, function(x) x$estimate())),
      se=function() do.call(cbind, lapply(res, function(x) x$se())),
      stat=function() do.call(cbind, lapply(res, function(x) x$stat())),
      prob=function() do.call(cbind, lapply(res, function(x) x$prob())),
      stat_type="zstat"
    )
  )
}