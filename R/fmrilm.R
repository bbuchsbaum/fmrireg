

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
#' con2  <- unit_contrast(~ fac, name="A_min_baseline")
#' lm.1 <- fmri_lm(onset ~ hrf(fac, contrasts=contrast_set(con,con2)), block= ~ run, dataset=dset)
#' lm.2 <- fmri_lm(onset ~ hrf(fac, contrasts=con), block= ~ run, dataset=dset2)
#' lm.2a <- fmri_lm(onset ~ hrf(fac, contrasts=con), block= ~ run, robust=TRUE, dataset=dset2)
#' lm.3 <- fmri_lm(onset ~ hrf(fac, contrasts=con), block= ~ run, dataset=dset2, strategy="runwise")
#' lm.3a <- fmri_lm(onset ~ hrf(fac, contrasts=con), block= ~ run, robust=TRUE, dataset=dset2, strategy="runwise")
#' @export
fmri_lm <- function(formula, block, baseline_model=NULL, dataset, 
                     durations, drop_empty=TRUE, contrasts=NULL, robust=FALSE,
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
    runwise_lm(dataset, model, conlist, fcons, robust=robust)
  } else if (strategy == "chunkwise") {
    chunkwise_lm(dataset, model, conlist,fcons, nchunks, robust=robust)
  }
  
  ret <- list(
    result=result,
    model=model,
    contrasts=contrasts,
    strategy=strategy,
    fcons=fcons,
    bcons=conlist)
  
  class(ret) <- "fmri_lm"
  
  ret
}


#' @export
coef.fmri_lm <- function(x, type=c("estimates", "contrasts")) {
  type <- match.arg(type)
  if (type == "estimates") {
    ret <- x$result$betas$estimate()
    colnames(ret) <- shortnames(x$model$event_model)#conditions(x$model$event_model)
    as_tibble(ret)
  } else if (type == "contrasts") {
    ret <- x$result$contrasts$estimate()
    colnames(ret) <- names(x$bcons)
    as_tibble(ret)
  } #else if (type == "baseline") {
    #ret <- x$result$contrasts$estimate()
  #}
}

#' @export
stats.fmri_lm <- function(x, type=c("estimates", "contrasts", "F")) {
  type <- match.arg(type)
  if (type == "estimates") {
    ret <- x$result$betas$stat()
    colnames(ret) <- shortnames(x$model$event_model)
    as_tibble(ret)
  } else if (type == "contrasts") {
    if (length(x$result$contrasts) == 0) {
      stop("no computed contrasts for this model.")
    }
    ret <- x$result$contrasts$stat()
    colnames(ret) <- names(x$bcons)
    as_tibble(ret)
  } else if (type == "F") {
    ret <- x$result$Fcontrasts
    ret <- do.call(cbind, lapply(ret, function(f) f$stat()))
    colnames(ret) <- names(x$fcons)
    as_tibble(ret)
  }
}

#' @export
standard_error.fmri_lm <- function(x, type=c("estimates", "contrasts")) {
  type <- match.arg(type)
  if (type == "estimates") {
    ret <- x$result$betas$se()
    colnames(ret) <- shortnames(x$model$event_model)
    as_tibble(ret)
  } else if (type == "contrasts") {
    if (length(x$result$contrasts) == 0) {
      stop("no computed contrasts for this model.")
    }
    ret <- x$result$contrasts$se()
    colnames(ret) <- names(x$bcons)
    as_tibble(ret)
  } else if (type == "F") {
    ret <- x$result$Fcontrasts
    ret <- do.call(cbind, lapply(ret, function(f) f$se()))
    colnames(ret) <- names(x$fcons)
    as_tibble(ret)
  }
}
  


print.fmri_lm <- function(x) {
  cat("fmri_lm model: \n", as.character(x$model$event_model$model_spec$formula), "\n")
  cat("  baseline parameters: ", ncol(design_matrix(x$model$baseline_model)), "\n")
  cat("  design parameters: ", ncol(design_matrix(x$model$event_model)), "\n")
  cat("  contrasts: ", names(x$bcons), "\n")
  cat("  F-contrasts: ", names(x$fcons), "\n")
}

# summary.fmri_lm <- function(x, type=c("coef", "contrasts", "Fcontrasts")) {
#   type <- match.arg(type)
#   if (type == "coef") {
#     betas=x$result$betas
#     list(
#       estimate=betas$estimate(),
#       se=betas$se(),
#       stat=betas$stat(),
#       prob=betas$prob())
#   } else if (type == "contrasts") {
#     x$result$contrasts
#   } else if (type == "Fcontrasts") {
#     x$result$Fcontrasts
#   } else {
#     stop()
#   }
# }


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

fit_lm_contrasts <- function(fit, conlist, fcon, vnames) {
  conres <- if (!is.null(conlist)) {
    ret <- lapply(conlist, function(con) {
      fit_contrasts(fit, con$weights, attr(con, "term_indices"))
    })
    
    names(ret) <- names(conlist)
    ret
  } 
  
  Fres <- lapply(fcon, function(con) fit_Fcontrasts(fit, t(con), attr(con, "term_indices")))
  names(Fres) <- names(fcon)
  
  bstats <- beta_stats(fit, vnames)
  #list(conres=conres, Fres=Fres, bstats=bstats, event_indices=eterm_indices, baseline_indices=bterm_indices)
  list(conres=conres, Fres=Fres, bstats=bstats)
}

#' @keywords internal
multiresponse_lm <- function(form, data_env, conlist, vnames, fcon, modmat=NULL) {
  lm.1 <- if (is.null(modmat)) {
    lm(as.formula(form), data=data_env)
  } else {
    lm.fit(modmat, data_env$.y)
  }
  
  fit_lm_contrasts(lm.1, conlist, fcon, vnames)
}


#' a beastly function that unravels multiple chunk-wise model fits results ...
#' @keywords internal
wrap_chunked_lm_results <- function(cres, event_indices=NULL) {
  if (!is.null(event_indices)) {
    force(event_indices)
  }
  
  standard_cols <- c("estimate", "stat", "se", "prob")
  
  extract <- function(l, el, item) {
    ret <- do.call(rbind, lapply(l, function(x) x[[el]][[item]]()))
    if (is.null(event_indices)) {
      ret
    } else {
      ret[,event_indices]
    }
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
  
  bstats <- c(do_extract(cres, "bstats", standard_cols, extract), 
              list(stat_type=cres[[1]]$bstats$stat_type))
 
  ncon <- length(cres[[1]]$conres)
  
  conres <- if (ncon >= 1) {
    ret <- lapply(1:ncon, function(i)  {
      x <- do_extract(cres, "conres", standard_cols, extract2, i)
      c(x, list(stat_type=cres[[1]]$conres[[i]]$stat_type))
    })
    
    names(ret) <- names(cres[[1]]$conres)
    ret
  }
    

  nf <- length(cres[[1]]$Fres)
  
  Fres <- if (nf >= 1) {
    ret <- lapply(1:nf, function(i)  {
      x <- do_extract(cres, "Fres", standard_cols,extract2, i)
      c(x, list(stat_type=cres[[1]]$Fres[[i]]$stat_type))
    })
    names(ret) <- names(cres[[1]]$Fres)
    ret
  }
  
  list(betas=bstats, contrasts=conres, Fcontrasts=Fres)
    
      
}



#' Run glm for each data chunk (responses split horizontally) and then concatenate chunks
#' @keywords internal
chunkwise_lm <- function(dset, model, conlist, fcon, nchunks, robust=FALSE, verbose=FALSE) {
  chunks <- exec_strategy("chunkwise", nchunks)(dset)
  form <- get_formula(model)
  tmats <- term_matrices(model)
  data_env <- list2env(tmats)
  data_env[[".y"]] <- rep(0, nrow(tmats[[1]]))
  modmat <- model.matrix(as.formula(form), data_env)
  
  lmfun <- if (robust) multiresponse_rlm else multiresponse_lm
  cres <- foreach( ym = chunks, .verbose=verbose) %dopar% {
    data_env[[".y"]] <- ym$data
    ret <- lmfun(form, data_env, conlist, attr(tmats,"varnames"), fcon, modmat=modmat)
  }
  
  event_indices=attr(tmats, "event_term_indices")
  baseline_indices=attr(tmats, "baseline_term_indices")
  
  wrap_chunked_lm_results(cres, event_indices)

}


#' Run glm for each data run (responses split vertically) and then combine over runs via meta-analysis
#' 
#' @importFrom foreach foreach %do% %dopar%
runwise_lm <- function(dset, model, conlist, fcon, robust=FALSE, verbose=FALSE) {
    #method <- match.arg(method)
    
    lmfun <- if (robust) {
      multiresponse_rlm
    } else {
      multiresponse_lm
    }
  
    ## get an iterator of data chunks
    chunks <- exec_strategy("runwise")(dset)

    form <- get_formula(model)
   
    ## iterate over each data chunk
    cres <- foreach( ym = chunks, .verbose=verbose) %dopar% {
      
      ## get event model for the nth run
      tmats <- term_matrices(model, ym$chunk_num)
      
      data_env <- list2env(tmats)
      data_env[[".y"]] <- as.matrix(ym$data)
      ret <- lmfun(form, data_env, conlist, attr(tmats,"varnames"), fcon)

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
        meta_contrasts(conres[hascon]) 
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
  
    
    



