

#' fmri_rlm
#' 
#' robust linear modeling of fmri data
#' 
#' @inheritParams fmri_lm
#' @examples 
#' etab <- data.frame(onset=c(1,30,15,25), fac=factor(c("A", "B", "A", "B")), run=c(1,1,2,2))
#' etab2 <- data.frame(onset=c(1,30,65,75), fac=factor(c("A", "B", "A", "B")), run=c(1,1,1,1))
#' mat <- matrix(rnorm(100*100), 100,100)
#' dset <- matrix_dataset(mat, TR=1, run_length=c(50,50),event_table=etab)
#' dset2 <- matrix_dataset(mat, TR=1, run_length=c(100),event_table=etab2)
#' lm.1 <- fmri_rlm(onset ~ hrf(fac), block_formula= ~ run,dataset=dset)
#' lm.2 <- fmri_rlm(onset ~ hrf(fac), block_formula= ~ run,dataset=dset2)
#' @export
fmri_rlm <- function(formula, block_formula, baseline_model=NULL, dataset, 
                     durations, drop_empty=TRUE, contrasts=NULL, 
                     strategy=c("runwise", "slicewise", "all")) {
  
 
  strategy <- match.arg(strategy)
  
  assert_that(inherits(dataset, "fmri_dataset"))
  fobj <- .setup_model(dataset, formula, block_formula, baseline_model, contrasts)

  result <- if (strategy == "runwise") {
    runwise_rlm(dataset, fobj$model, fobj$conlist, fobj$fcon)
  } else {
    stop()
  }
  
  ret <- list(
    result=result,
    model=fobj$model,
    contrasts=contrasts,
    strategy=strategy)
  
  class(ret) <- c("fmri_rlm", "fmri_lm")
  
  ret

}

#' @importFrom foreach foreach %do% %dopar%
runwise_rlm <- function(dset, model, conlist, fcon) {
  
  ## get an iterator of data chunks
  chunks <- exec_strategy("runwise")(dset)
  
  term_names <- names(terms(model))
  form <- paste(".y ~ ", paste(term_names, collapse = " + "), "-1")
  
  ## iterate over each data chunk
  cres <- foreach( ym = chunks) %do% {
    
    ## get event model for the nth run
    eterm_matrices <- lapply(event_terms(model), 
                             function(x) as.matrix(design_matrix(x, ym$chunk_num)))
    
    ## get baseline model for the nth run
    bterm_matrices <- lapply(baseline_terms(model), 
                             function(x) as.matrix(design_matrix(x, ym$chunk_num)))
    
    ## column indices of event regressors
    eterm_indices <- 1:sum(sapply(eterm_matrices, ncol))
    start <- length(eterm_indices) +1
    
    ## column indices of baseline regressors
    bterm_indices <- start:(start+sum(sapply(bterm_matrices, ncol)))
    
    term_matrices <- c(eterm_matrices, bterm_matrices)
    names(term_matrices) <- term_names
    
    
    data_env <- list2env(term_matrices)
    vnames <- unlist(lapply(term_matrices, colnames))
    
    lapply(1:ncol(ym$data), function(i) {
      data_env[[".y"]] <- ym$data[,i]
      rlm.1 <- lmrob(as.formula(form), data=data_env)
    
      conres <- lapply(conlist, function(con) fit_contrasts(rlm.1, con, attr(con, "term_indices")))
      names(conres) <- names(conlist)
    
      Fres <- lapply(fcon, function(con) fit_Fcontrasts(rlm.1, t(con), attr(con, "term_indices")))
    
      bstats <- beta_stats(rlm.1, vnames)
      list(conres=conres, Fres=Fres, bstats=bstats)
    })
    
  }
  
}