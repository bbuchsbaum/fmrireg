
.setup_model <- function(dataset, formula, block_formula, baseline_model=NULL, contrasts=NULL) {
  assert_that(inherits(dataset, "fmri_dataset"))
  
  if (is.null(baseline_model)) {
    baseline_model <- baseline_model(basis="bs", 
                                     degree=ceiling(median(dataset$sampling_frame$blocklens)/100), 
                                     sframe=dataset$sampling_frame)
  }
  
  
  ev_model <- event_model(formula, block_formula, data=dataset$event_table, 
                          sampling_frame=dataset$sampling_frame, contrasts=contrasts)
  model <- fmri_model(ev_model, baseline_model)
  
  conlist <- contrast_weights(ev_model)
  fcon <- Fcontrasts(ev_model)
  
  list(baseline_model=baseline_model, ev_model=ev_model, conlist=conlist,fcon=fcon)
}

#' fmri_lm
#' 
#' @param formula the model furmula for experimental events
#' @param block_formula the model formula for block structure
#' @param baseline_model the \code{baseline_model} object
#' @param dataset an object derived from \code{fmri_dataset} containing the time-series data
#' @param durations a vector of event durations
#' @param drop_empty whether to remove factor levels with size of zero
#' @param contrasts a set of contrasts
#' @param strategy the data splitting strategy
#' 
#' @examples 
#' etab <- data.frame(onset=c(1,30,15,25), fac=factor(c("A", "B", "A", "B")), run=c(1,1,2,2))
#' etab2 <- data.frame(onset=c(1,30,65,75), fac=factor(c("A", "B", "A", "B")), run=c(1,1,1,1))
#' mat <- matrix(rnorm(100*100), 100,100)
#' dset <- matrix_dataset(mat, TR=1, run_length=c(50,50),event_table=etab)
#' dset2 <- matrix_dataset(mat, TR=1, run_length=c(100),event_table=etab2)
#' lm.1 <- fmri_lm(onset ~ hrf(fac), block_formula= ~ run,dataset=dset)
#' lm.2 <- fmri_lm(onset ~ hrf(fac), block_formula= ~ run,dataset=dset2)
#' @export
fmri_lm <- function(formula, block_formula, baseline_model=NULL, dataset, 
                     durations, drop_empty=TRUE, contrasts=NULL, 
                     strategy=c("runwise", "slicewise", "all")) {
  
 
  strategy <- match.arg(strategy)
  
  assert_that(inherits(dataset, "fmri_dataset"))
  
  if (is.null(baseline_model)) {
    baseline_model <- baseline_model(basis="bs", 
                                     degree=ceiling(median(dataset$sampling_frame$blocklens)/100), 
                                     sframe=dataset$sampling_frame)
  }
  

  ev_model <- event_model(formula, block_formula, data=dataset$event_table, sampling_frame=dataset$sampling_frame, contrasts=contrasts)
  model <- fmri_model(ev_model, baseline_model)
  
  conlist <- contrast_weights(ev_model)
  fcon <- Fcontrasts(ev_model)
  
  result <- if (strategy == "runwise") {
    runwise_lm(dataset, model, conlist, fcon)
  } else {
    stop()
  }
  
  ret <- list(
    result=result,
    model=model,
    contrasts=contrasts,
    strategy=strategy)
  
  class(ret) <- "fmri_lm"
  
  ret
}


#' @export
coef.fmri_lm <- function(x) {
 x$result$betas$estimate()
}

#' @export
stats.fmri_lm <- function(x) {
  x$result$betas$stat()
}

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
  res <- lapply(1:ncon, function(i) {
    pval <- do.call(cbind, lapply(fres, function(x) as.vector(x[[i]]$prob())))
    se <- do.call(cbind, lapply(fres, function(x) x[[i]]$se()))
    
    meta_stouffer(pval,se)
  })
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

combine_baseline_betas <- function(bstats, colind) {
  list(
    estimate=function() do.call(rbind, lapply(bstats, function(x) x$estimate()[colind,])),
    se=function() do.call(rbind, lapply(bstats, function(x) x$se()[colind,])),
    stat=function() do.call(rbind, lapply(bstats, function(x) x$stat()[colind,])),
    prob=function() do.call(rbind, lapply(bstats, function(x) x$prob()[colind,])),
    stat_type="tstat"
  )
}

multiresponse_lm <- function(form, data_env, conlist, vnames, fcon) {
  lm.1 <- lm(as.formula(form), data=data_env)
  
  conres <- lapply(conlist, function(con) fit_contrasts(lm.1, con, attr(con, "term_indices")))
  names(conres) <- names(conlist)
  
  Fres <- lapply(fcon, function(con) fit_Fcontrasts(lm.1, t(con), attr(con, "term_indices")))
  
  bstats <- beta_stats(lm.1, vnames)
  #list(conres=conres, Fres=Fres, bstats=bstats, event_indices=eterm_indices, baseline_indices=bterm_indices)
  list(conres=conres, Fres=Fres, bstats=bstats)
}

#' @importFrom foreach foreach %do% %dopar%
runwise_lm <- function(dset, model, conlist, fcon) {
    chunks <- exec_strategy("runwise")(dset)
    
    term_names <- names(terms(model))
    form <- paste(".y ~ ", paste(term_names, collapse = " + "), "-1")

    cres <- foreach( ym = chunks) %do% {
      
      ## get event model for the nth run
      eterm_matrices <- lapply(event_terms(model), function(x) as.matrix(design_matrix(x, ym$chunk_num)))
      
      ## get baseline model for the nth run
      bterm_matrices <- lapply(baseline_terms(model), function(x) as.matrix(design_matrix(x, ym$chunk_num)))
      
      ## column indices of event regressors
      eterm_indices <- 1:sum(sapply(eterm_matrices, ncol))
      start <- length(eterm_indices) +1
      
      ## column indices of baseline regressors
      bterm_indices <- start:(start+sum(sapply(bterm_matrices, ncol)))
     
      term_matrices <- c(eterm_matrices, bterm_matrices)
      names(term_matrices) <- term_names
      
      
      data_env <- list2env(term_matrices)
      data_env[[".y"]] <- ym$data
      
      vnames <- unlist(lapply(term_matrices, colnames))
      
      ret <- multiresponse_lm(form, data_env, conlist, vnames, fcon)
      list(conres=ret$conres, Fres=ret$Fres, bstats=ret$bstats, event_indices=eterm_indices, baseline_indices=bterm_indices)
      
    }
    
    bstats <- lapply(cres, "[[", "bstats")
    conres <- lapply(cres, "[[", "conres")
    Fres <- lapply(cres, "[[", "Fres")

    if (length(cres) > 1) {
      meta_con <- if (length(conres[[1]]) > 0) meta_contrasts(conres) else list()
      meta_beta <- meta_betas(bstats, cres[[1]]$event_indices)
      meta_F <- meta_Fcontrasts(Fres)
      list(contrasts=meta_con, betas=meta_beta, Fcontrasts=meta_F)
    } else {
      list(contrasts=conres[[1]], betas=bstats[[1]], Fcontrasts=Fres[[1]])
    }
    
}
  
    
    



