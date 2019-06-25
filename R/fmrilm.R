

#' @export
get_formula.fmri_model <- function(x) {
  term_names <- names(terms(x))
  form <- paste(".y ~ ", paste(term_names, collapse = " + "), "-1")
}


#' @export
term_matrices.fmri_model <- function(x, blocknum=NULL) {
  eterms <- lapply(event_terms(x), 
                   function(x) as.matrix(design_matrix(x, blocknum)))
  
  bterms <- lapply(baseline_terms(x), 
                   function(x) as.matrix(design_matrix(x, blocknum)))
  
  if (is.null(blocknum)) {
    blocknum <- sort(unique(x$event_model$blockids))
  }
  
  start <- 1
  eterm_indices <- 1:sum(map_int(eterms, ncol))
  start <- length(eterm_indices) +1
  bterm_indices <- start:(start+sum(map_int(bterms, ncol)))
  #browser()
  term_matrices <- c(eterms, bterms)
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
#' @param block the model formula for block structure
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
#' con  <- pair_contrast(~ fac == "A", ~ fac == "B", name="A_min_B")
#' lm.1 <- fmri_lm(onset ~ hrf(fac, contrasts=con), block= ~ run, dataset=dset)
#' lm.2 <- fmri_lm(onset ~ hrf(fac, contrasts=con), block= ~ run, dataset=dset2)
#' lm.3 <- fmri_lm(onset ~ hrf(fac, contrasts=con), block= ~ run, dataset=dset2, strategy="chunkwise")
#' @export
fmri_lm <- function(formula, block, baseline_model=NULL, dataset, 
                     durations, drop_empty=TRUE, contrasts=NULL, 
                     strategy=c("runwise", "chunkwise"), nchunks=10) {
  
 
  strategy <- match.arg(strategy)
  
  assert_that(inherits(dataset, "fmri_dataset"))
  
 
  if (is.null(baseline_model)) {
    baseline_model <- baseline_model(basis="bs", 
                                    degree=ceiling(median(dataset$sampling_frame$blocklens)/100), 
                                    sframe=dataset$sampling_frame)
  }
 
  ev_model <- event_model(formula, block, data=dataset$event_table, 
                         sampling_frame=dataset$sampling_frame, contrasts=contrasts)
     
  model <- fmri_model(ev_model, baseline_model)
  conlist <- unlist(contrast_weights(ev_model), recursive=FALSE)

  fcons <- Fcontrasts(ev_model)

  result <- if (strategy == "runwise") {
    runwise_lm(dataset, model, conlist, fcons)
  } else if (strategy == "chunkwise") {
    chunkwise_lm(dataset, model, conlist,fcons, nchunks)
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


#' @keywords internal
combine_baseline_betas <- function(bstats, colind) {
  list(
    estimate=function() do.call(rbind, lapply(bstats, function(x) x$estimate()[colind,])),
    se=function() do.call(rbind, lapply(bstats, function(x) x$se()[colind,])),
    stat=function() do.call(rbind, lapply(bstats, function(x) x$stat()[colind,])),
    prob=function() do.call(rbind, lapply(bstats, function(x) x$prob()[colind,])),
    stat_type="tstat"
  )
}

#' @keywords internal
multiresponse_lm <- function(form, data_env, conlist, vnames, fcon, modmat=NULL) {
  lm.1 <- if (is.null(modmat)) {
    lm(as.formula(form), data=data_env)
  } else {
    lm.fit(modmat, data_env$.y)
  }
  
  conres <- if (!is.null(conlist)) {
    ret <- lapply(conlist, function(con) {
      fit_contrasts(lm.1, con$weights, attr(con, "term_indices"))
    })
  
    names(ret) <- names(conlist)
    ret
  } 

  Fres <- lapply(fcon, function(con) fit_Fcontrasts(lm.1, t(con), attr(con, "term_indices")))
  names(Fres) <- names(fcon)
  
  bstats <- beta_stats(lm.1, vnames)
  #list(conres=conres, Fres=Fres, bstats=bstats, event_indices=eterm_indices, baseline_indices=bterm_indices)
  list(conres=conres, Fres=Fres, bstats=bstats)
}

#' @keywords internal
wrap_chunked_lm_results <- function(cres) {
  
  extract <- function(l, el, item) {
    do.call(rbind, lapply(l, function(x) x[[el]][[item]]()))
  }
  
  extract2 <- function(l, el, item, i) {
    do.call(rbind, lapply(l, function(x) x[[el]][[i]][[item]]()))
  }
  
  do_extract <- function(l, el, items, efun,...) {
    res <- lapply(items, function(item) {
      function() { efun(l, el, item,...) }
    })
    names(res) <- items
    res
  }
  
  bstats <- c(do_extract(cres, "bstats", c("estimate", "stat", "se", "prob"), extract), 
              list(stat_type=cres[[1]]$bstats$stat_type))
 
  ncon <- length(cres[[1]]$conres)
  
  conres <- if (ncon >= 1) {
    ret <- lapply(1:ncon, function(i)  {
      x <- do_extract(cres, "conres", c("estimate", "stat", "se", "prob"),extract2, i)
      c(x, list(stat_type=cres[[1]]$conres[[i]]$stat_type))
    })
    
    names(ret) <- names(cres[[1]]$conres)
    ret
  }
    

  nf <- length(cres[[1]]$Fres)
  
  Fres <- if (nf >= 1) {
    ret <- lapply(1:nf, function(i)  {
      x <- do_extract(cres, "Fres", c("estimate", "stat", "se", "prob"),extract2, i)
      c(x, list(stat_type=cres[[1]]$Fres[[i]]$stat_type))
    })
    names(ret) <- names(cres[[1]]$Fres)
    ret
  }
  
  list(betas=bstats, contrasts=conres, Fcontrasts=Fres)
    
      
}


#' @keywords internal
chunkwise_lm <- function(dset, model, conlist, fcon, nchunks) {
  chunks <- exec_strategy("chunkwise", nchunks)(dset)
  form <- get_formula(model)
  tmats <- term_matrices(model)
  data_env <- list2env(tmats)
  data_env[[".y"]] <- rep(0, nrow(tmats[[1]]))
  modmat <- model.matrix(as.formula(form), data_env)
  cres <- foreach( ym = chunks, .verbose=TRUE) %dopar% {
    data_env[[".y"]] <- ym$data
    ret <- multiresponse_lm(form, data_env, conlist, attr(tmats,"varnames"), fcon, modmat=modmat)
  }
  
  wrap_chunked_lm_results(cres)

}

#' @importFrom foreach foreach %do% %dopar%
runwise_lm <- function(dset, model, conlist, fcon) {
  
    ## get an iterator of data chunks
    chunks <- exec_strategy("runwise")(dset)

    form <- get_formula(model)
   
    ## iterate over each data chunk
    cres <- foreach( ym = chunks, .verbose=TRUE) %dopar% {
      
      ## get event model for the nth run
      tmats <- term_matrices(model, ym$chunk_num)
      
      data_env <- list2env(tmats)
      data_env[[".y"]] <- ym$data
      
      ret <- multiresponse_lm(form, data_env, conlist, attr(tmats,"varnames"), fcon)
  
      list(conres=ret$conres, Fres=ret$Fres, bstats=ret$bstats, 
           event_indices=attr(tmats, "event_term_indices"), 
           baseline_indices=attr(tmats, "baseline_term_indices") )
      
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
  
    
    



