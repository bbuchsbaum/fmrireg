

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
  
  if (strategy == "runwise") {
    runwise_lm(dataset, model, conlist, fcon)
  } else {
    stop()
  }
 
  
  model
}

meta_fixed <- function(beta, vars, vwts) {
  summ<-sum(vwts*beta)/sum(vwts)
  varsum<-sum(vwts*vwts*vars)/(sum(vwts)^2)
  summtest<-summ/sqrt(varsum)
  c(b=summ, se=sqrt(varsum), z=summtest)
}


meta_fixef <- function(beta,se) {
  inv_var <- 1/(se^2)
  
  wts <- inv_var/rowSums(inv_var)
  wbeta <- beta * wts
  wbeta <- rowSums(wbeta)
  pooledse <- sqrt(rowSums(wts*wts*(se^2)))
  tibble::data_frame(b=wbeta, se=pooledse, z=wbeta/pooledse)
}

meta_contrasts <- function(cres) {
  ncon <- length(cres[[1]])
  res <- lapply(1:ncon, function(i) {
    beta <- do.call(cbind, lapply(cres, function(x) x[[i]]$estimate))
    se <- do.call(cbind, lapply(cres, function(x) x[[i]]$se))
    meta_fixef(beta,se)
  })
}

meta_betas <- function(bstats, colind) {
  len <- length(colind)

  res <- lapply(colind, function(i) {
    print(i)
    beta <- do.call(cbind, lapply(bstats, function(x) x$beta[i,]))
    se <- do.call(cbind, lapply(bstats, function(x) x$se[i,]))
    meta_fixef(beta,se)
  })
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
      browser()
      list(contrasts=meta_con, betas=meta_beta)
    } else {
      list(contrasts=conres[[1]], betas=bstats[[1]], Fstats=Fres[[1]])
    }
    
}
  
    
    



