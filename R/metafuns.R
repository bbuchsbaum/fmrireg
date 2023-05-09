

#' Meta-analysis using Stouffer's method
#'
#' This function performs a meta-analysis on input p-values and standard errors using Stouffer's method.
#' Stouffer's method combines z-scores by weighting them by their inverse variance.
#' 
#' @param pval A numeric vector of p-values from multiple studies or tests.
#' @param se A numeric vector of standard errors corresponding to the p-values.
#' 
#' @return A list containing the following elements:
#' \itemize{
#'   \item{estimate}{A numeric vector of the weighted z-scores obtained using Stouffer's method.}
#'   \item{se}{A numeric vector of the combined standard errors.}
#'   \item{stat}{Same as the estimate, the weighted z-scores.}
#'   \item{prob}{A numeric vector of the combined p-values.}
#'   \item{stat_type}{A character string indicating the type of statistic used, in this case "zfstat".}
#' }
#' 
#' @keywords internal
#' @examples
#' pval <- c(0.05, 0.01, 0.03)
#' se <- c(0.2, 0.15, 0.25)
#' result <- meta_stouffer(pval, se)
#' print(result)
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


#' Fixed-effects meta-analysis
#'
#' This function performs a fixed-effects meta-analysis on input beta estimates and standard errors.
#' It supports two types of weighting schemes: inverse-variance weighting and equal weighting.
#' 
#' @param beta A numeric matrix of beta estimates from multiple studies or tests.
#' @param se A numeric matrix of standard errors corresponding to the beta estimates.
#' @param weighting A character string specifying the weighting scheme to use. 
#'   Options are "inv_var" (default) for inverse-variance weighting and "equal" for equal weighting.
#' 
#' @return A list containing the following elements:
#' \itemize{
#'   \item{estimate}{A numeric vector of the combined beta estimates.}
#'   \item{se}{A numeric vector of the pooled standard errors.}
#'   \item{stat}{A numeric vector of the combined z-scores.}
#'   \item{prob}{A numeric vector of the combined p-values.}
#'   \item{stat_type}{A character string indicating the type of statistic used, in this case "zstat".}
#' }
#' 
#' @keywords internal
#' @examples
#' beta <- matrix(c(0.5, 1, 0.75), nrow=1)
#' se <- matrix(c(0.2, 0.15, 0.25), nrow=1)
#' result <- meta_fixef(beta, se, weighting = "inv_var")
#' print(result)
meta_fixef <- function(beta,se, weighting=c("inv_var", "equal")) {
  weighting <- match.arg(weighting)
  
  if (weighting == "inv_var") {
    inv_var <- 1/(se^2)
    wts <- inv_var/rowSums(inv_var)
    wbeta <- beta * wts
    wbeta <- rowSums(wbeta)
    pooledse <- sqrt(rowSums(wts*wts*(se^2)))
  } else {
    wbeta <- rowSums(beta)
    pooledse <- sqrt(rowSums((se^2)))
  }
  
  return(
    list(
      estimate=wbeta,
      se=pooledse,
      stat=wbeta/pooledse,
      prob=1-pchisq((wbeta/pooledse)^2,1),
      stat_type="zstat")
  )
  
}


#' Meta-analysis of F-contrasts
#'
#' This function performs a meta-analysis on a list of F-contrasts results obtained from different studies or tests.
#'
#' @param fres A list of F-contrast results, where each element of the list contains the results for a particular study or test.
#'
#' @return A list containing meta-analysis results for each contrast. Each element of the list includes:
#' \itemize{
#'   \item{estimate}{A numeric vector of the combined z-scores.}
#'   \item{se}{A numeric vector of the pooled standard errors.}
#'   \item{stat}{A numeric vector of the combined z-scores.}
#'   \item{prob}{A numeric vector of the combined p-values.}
#'   \item{stat_type}{A character string indicating the type of statistic used, in this case "zfstat".}
#' }
#' 
#' @keywords internal
#' @examples
#' # Assuming `fres` is a list of F-contrast results obtained from different studies or tests
#' # meta_results <- meta_Fcontrasts(fres)
#' # print(meta_results)
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


#' Meta-analysis of contrasts
#'
#' This function performs a meta-analysis on a list of contrast results obtained from different studies or tests.
#'
#' @param cres A list of contrast results, where each element of the list contains the results for a particular study or test.
#' @param weighting A character string specifying the weighting method to use for the meta-analysis. Options are "inv_var" (inverse variance weighting) or "equal" (equal weighting). Default is "inv_var".
#'
#' @return A list containing meta-analysis results for each contrast, including:
#' \itemize{
#'   \item{estimate}{A matrix of the combined contrast estimates.}
#'   \item{se}{A matrix of the pooled standard errors.}
#'   \item{stat}{A matrix of the combined z-scores.}
#'   \item{prob}{A matrix of the combined p-values.}
#'   \item{stat_type}{A character string indicating the type of statistic used, in this case "zstat".}
#' }
#' 
#' @keywords internal
#' @examples
#' # Assuming `cres` is a list of contrast results obtained from different studies or tests
#' # meta_results <- meta_contrasts(cres)
#' # print(meta_results)
meta_contrasts <- function(cres, weighting=c("inv_var", "equal")) {
  weighting <- match.arg(weighting)
  ncon <- length(cres[[1]])
  if (ncon > 0) {
    res <- lapply(1:ncon, function(i) {
      beta <- do.call(cbind, lapply(cres, function(x) as.vector(x[[i]]$estimate)))
      se <- do.call(cbind, lapply(cres, function(x) x[[i]]$se))
      meta_fixef(beta,se, weighting)
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
meta_betas <- function(bstats, colind, weighting=c("inv_var", "equal")) {
  weighting <- match.arg(weighting)
  
  len <- length(colind)
  
  res <- lapply(colind, function(i) {
    #print(i)
    beta <- do.call(cbind, lapply(bstats, function(x) x$estimate[,i]))
    se <- do.call(cbind, lapply(bstats, function(x) x$se[,i]))
    meta_fixef(beta,se,weighting)
  })
  
  #browser()
  
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