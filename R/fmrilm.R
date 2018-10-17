
get_model_formula.fmri_model <- function(x) {
  term_names <- names(terms(x))
  form <- paste(".y ~ ", paste(term_names, collapse = " + "), "-1")
}

term_matrices.fmri_model <- function(x, blocknum=NULL) {
  eterms <- lapply(event_terms(x), 
                   function(x) as.matrix(design_matrix(x, blocknum)))
  
  bterms <- lapply(baseline_terms(x), 
                   function(x) as.matrix(design_matrix(x, blocknum)))
  
  if (is.null(blocknum)) {
    blocknum <- unique(fmodel$event_model$blockids)
  }
  
  start <- 1
  eterm_indices <- 1:sum(map_int(eterms, ncol))
  start <- length(eterm_indices) +1
  bterm_indices <- start:(start+sum(map_int(bterms, ncol)))
  
  term_matrices <- c(tmats$event_terms, tmats$baseline_terms)
  names(term_matrices) <- names(terms(x))
  
  vnames <- unlist(lapply(term_matrices, colnames))
  
  attr(term_matrices, "event_term_indices") <- eterm_indices
  attr(term_matrices, "baseline_term_indices") <- bterm_indices
  attr(term_matrices, "blocknum") <- blocknum 
  attr(term_matrices, "varnames") <- vnames
  term_matrices
}


# fmri_lm_model <- function(dataset, formula, block_formula, baseline_model=NULL, contrasts=NULL) {
#   assert_that(inherits(dataset, "fmri_dataset"))
#   
#   if (is.null(baseline_model)) {
#     baseline_model <- baseline_model(basis="bs", 
#                                      degree=ceiling(median(dataset$sampling_frame$blocklens)/100), 
#                                      sframe=dataset$sampling_frame)
#   }
#   
#   
#   ev_model <- event_model(formula, block_formula, data=dataset$event_table, 
#                           sampling_frame=dataset$sampling_frame, contrasts=contrasts)
#   
#   model <- fmri_model(ev_model, baseline_model)
#   
#   conlist <- contrast_weights(ev_model)
#   fcon <- Fcontrasts(ev_model)
#   
#   ret <- list(model=model, 
#               conlist=conlist,
#               fcon=fcon)
#   
#   class(ret) <- "fmri_lm_model"
#   ret
# }

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
                     strategy=c("runwise", "vectorwise", "chunkwise", "all")) {
  
 
  strategy <- match.arg(strategy)
  
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
  fcons <- Fcontrasts(ev_model)

  result <- if (strategy == "runwise") {
    runwise_lm(dataset, model, conlist,fcons)
  } else {
    stop("only currently implemented strategy is 'runwise'")
  }
  
  ret <- list(
    result=result,
    model=fobj$model,
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
  
  conres <- lapply(conlist, function(con) {
    if (!is.null(con)) {
      fit_contrasts(lm.1, con, attr(con, "term_indices"))
    } else {
      NULL
    }
  })
  
  names(conres) <- names(conlist)

  Fres <- lapply(fcon, function(con) fit_Fcontrasts(lm.1, t(con), attr(con, "term_indices")))
  
  bstats <- beta_stats(lm.1, vnames)
  #list(conres=conres, Fres=Fres, bstats=bstats, event_indices=eterm_indices, baseline_indices=bterm_indices)
  list(conres=conres, Fres=Fres, bstats=bstats)
}

#' @importFrom foreach foreach %do% %dopar%
runwise_lm <- function(dset, model, conlist, fcon) {
  
    ## get an iterator of data chunks
    chunks <- exec_strategy("runwise")(dset)

    form <- get_formula(model)
   
    ## iterate over each data chunk
    cres <- foreach( ym = chunks, .verbose=TRUE) %do% {

      ## get event model for the nth run
      tmats <- term_matrices(model, ym$chunk_num)
      
      data_env <- list2env(tmats)
      data_env[[".y"]] <- ym$data
      
      
      ret <- multiresponse_lm(form, data_env, conlist, attr(tmats,"varnames"), fcon)
      
      list(conres=ret$conres, Fres=ret$Fres, bstats=ret$bstats, 
           event_indices=attr(tmats, "eterm_indices"), baseline_indices=bterm_indices)
      
    }
    
    bstats <- lapply(cres, "[[", "bstats")
    conres <- lapply(cres, "[[", "conres")
    Fres <- lapply(cres, "[[", "Fres")

    if (length(cres) > 1) {
      hascon <- !sapply(conres[[1]], is.null)
      meta_con <- if (any(hascon)) {
        meta_con <- meta_contrasts(conres[hascon]) 
      } else {
        list()
      }
      meta_beta <- meta_betas(bstats, cres[[1]]$event_indices)
      meta_F <- meta_Fcontrasts(Fres)
      list(contrasts=meta_con, betas=meta_beta, Fcontrasts=meta_F)
    } else {
      list(contrasts=conres[[1]], betas=bstats[[1]], Fcontrasts=Fres[[1]])
    }
    
}
  
    
    



