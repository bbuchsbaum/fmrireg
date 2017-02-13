

#' fmri_glm
#' 
#' @param formula
#' @param block_formula
#' @param baseline_model
#' @param dataset
#' @param durations
#' @param nuisance_matrix
#' @param drop_empty
#' @param contrasts
#' @param strategy
#' @param 
#' @export
fmri_lm <- function(formula, block_formula, baseline_model=NULL, dataset, 
                     durations, drop_empty=TRUE, contrasts=NULL, 
                     strategy=c("runwise", "slicewise", "all")) {
  
 
  strategy <- match.arg(strategy)
  
  assert_that(inherits(dataset, "fmri_dataset"))
  
  if (is.null(baseline_model)) {
    baseline_model <- baseline_model(basis="bs", degree=ceiling(median(dataset$sampling_frame$blocklens)/100), 
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
  
  browser()
  
  
 
  
  model
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
  res <- lapply(1:ncon, function(i) {
    beta <- do.call(cbind, lapply(cres, function(x) as.vector(x[[i]]$estimate())))
    se <- do.call(cbind, lapply(cres, function(x) x[[i]]$se()))
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

meta_betas <- function(bstats, colind) {
  len <- length(colind)

  res <- lapply(colind, function(i) {
    print(i)
    beta <- do.call(cbind, lapply(bstats, function(x) x$estimate()[i,]))
    se <- do.call(cbind, lapply(bstats, function(x) x$se()[i,]))
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



#' @importFrom foreach foreach %do% %dopar%
runwise_lm <- function(dset, model, conlist, fcon) {
    chunks <- exec_strategy("runwise")(dset)
    
    term_names <- names(terms(model))
    form <- as.formula(paste("y ~ ", paste(term_names, collapse = " + "), "-1"))
    
    cres <- foreach( ym = chunks) %do% {
      
      ## replace with function called "baseline_indices"
      eterm_matrices <- lapply(event_terms(model), function(x) as.matrix(design_matrix(x, ym$chunk_num)))
      bterm_matrices <- lapply(baseline_terms(model), function(x) as.matrix(design_matrix(x, ym$chunk_num)))
      
      eterm_indices <- 1:sum(sapply(eterm_matrices, ncol))
      start <- length(eterm_indices) +1
      bterm_indices <- start:(start+sum(sapply(bterm_matrices, ncol)))
     
      term_matrices <- c(eterm_matrices, bterm_matrices)
      names(term_matrices) <- term_names
      
      data_env <- list2env(term_matrices)
      
      y <- ym$data
      lm.1 <- lm(form, data=data_env)
    
      conres <- lapply(conlist, function(con) fit_contrasts(lm.1, con, attr(con, "term_indices")))
      names(conres) <- names(conlist)
      
      Fres <- lapply(fcon, function(con) fit_Fcontrasts(lm.1, t(con), attr(con, "term_indices")))
      
      bstats <- beta_stats(lm.1)
      list(conres=conres, Fres=Fres, bstats=bstats, event_indices=eterm_indices, baseline_indices=bterm_indices)
    }
    
    bstats <- lapply(cres, "[[", "bstats")
    conres <- lapply(cres, "[[", "conres")
    Fres <- lapply(cres, "[[", "Fres")
    
    if (length(cres) > 1) {
      meta_con <- meta_contrasts(conres)
      meta_beta <- meta_betas(bstats, cres[[1]]$event_indices)
      meta_F <- meta_Fcontrasts(Fres)
      list(contrasts=meta_con, betas=meta_beta, Fcontrasts=meta_F)
    } else {
      list(contrasts=conres[[1]], betas=bstats[[1]], Fcontrasts=Fres[[1]])
    }
    
}
  
    
    



