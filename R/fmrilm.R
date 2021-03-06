

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


#' fmri_lm
#' 
#' @param formula the model formula for experimental events
#' @param block the model formula for block structure
#' @param baseline_model the \code{baseline_model} object
#' @param dataset an object derived from \code{fmri_dataset} containing the time-series data
#' @param durations a vector of event durations
#' @param drop_empty whether to remove factor levels with size of zero
#' @param contrasts a set of contrasts
#' @param strategy the data splitting strategy
#' @param robust whether to use robust fitting (TRUE or FALSE)
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
                     durations, drop_empty=TRUE, robust=FALSE,
                     strategy=c("runwise", "chunkwise"), nchunks=10) {
  
 
  strategy <- match.arg(strategy)
  
  assert_that(inherits(dataset, "fmri_dataset"))
  
  if (is.null(baseline_model)) {
    baseline_model <- baseline_model(basis="bs", 
                                    degree=max(ceiling(median(dataset$sampling_frame$blocklens)/100),3), 
                                    sframe=dataset$sampling_frame)
  } else {
    assert_that(inherits(baseline_model, "baseline_model"), msg="'baseline_model' arg must have the class 'baseline_model'")
  }
 
  ev_model <- event_model(formula, block, data=dataset$event_table, 
                         sampling_frame=dataset$sampling_frame)
     
  model <- fmri_model(ev_model, baseline_model)
  
  ret <- fmri_lm_fit(model, dataset, strategy, robust, contrasts, nchunks)
  ret
}


# @export  ## do weed to export?
#' fmri_lm_fit
#' 
#' @inheritParams fmri_lm
#' @param fmrimod an object of type \code{fmri_model}
fmri_lm_fit <- function(fmrimod, dataset, strategy=c("chunkwise", "runwise"), robust=FALSE, contrasts=NULL, nchunks=10) {
  strategy <- match.arg(strategy)
  
  conlist <- unlist(contrast_weights(fmrimod$event_model), recursive=FALSE)
  fcons <- Fcontrasts(fmrimod$event_model)
  
  
  result <- if (strategy == "runwise") {
    runwise_lm(dataset, fmrimod, conlist, fcons, robust=robust)
  } else if (strategy == "chunkwise") {
    #message("chunkwise_lm with ", nchunks)
    chunkwise_lm(dataset, fmrimod, conlist,fcons, nchunks, robust=robust)
  }
  
  ret <- list(
    result=result,
    model=fmrimod,
    contrasts=contrasts,
    strategy=strategy,
    fcons=fcons,
    bcons=conlist,
    dataset=dataset)
  
  class(ret) <- "fmri_lm"
  
  ret
  
  
}


#' @export
coef.fmri_lm <- function(x, type=c("estimates", "contrasts"), recon=FALSE) {
  type <- match.arg(type)
  res <- if (type == "estimates") {
    ret <- x$result$betas$estimate()
    colnames(ret) <- conditions(x$model$event_model)
    #shortnames(x$model$event_model)#conditions(x$model$event_model)
    as_tibble(ret, .name_repair="check_unique")
  } else if (type == "contrasts") {
    ret <- x$result$contrasts$estimate()
    colnames(ret) <- names(x$bcons)
    as_tibble(ret, .name_repair="check_unique")
  } #else if (type == "baseline") {
    #ret <- x$result$contrasts$estimate()
  #}
  
  res
  
  # if (recon && inherits(x$dataset, "fmri_dataset")) {
  #   m <- get_mask(x$dataset)
  #   sp <- space(m)
  #   SparseNeuroVec(as.matrix(res), neuroim2::add_dim(sp, ncol(res)), mask=m)
  # } else {
  #   res
  # }
}

#' @export
stats.fmri_lm <- function(x, type=c("estimates", "contrasts", "F")) {
  type <- match.arg(type)
  if (type == "estimates") {
    ret <- x$result$betas$stat()
    #colnames(ret) <- shortnames(x$model$event_model)
    colnames(ret) <- conditions(x$model$event_model)
    as_tibble(ret,.name_repair="check_unique")
  } else if (type == "contrasts") {
    if (length(x$result$contrasts) == 0) {
      stop("no computed contrasts for this model.")
    }
    ret <- x$result$contrasts$stat()
    colnames(ret) <- names(x$bcons)
    as_tibble(ret, .name_repair="check_unique")
  } else if (type == "F") {
    ret <- x$result$Fcontrasts
    ret <- do.call(cbind, lapply(ret, function(f) f$stat()))
    colnames(ret) <- names(x$fcons)
    as_tibble(ret, .name_repair="check_unique")
  }
}

#' @export
standard_error.fmri_lm <- function(x, type=c("estimates", "contrasts")) {
  type <- match.arg(type)
  if (type == "estimates") {
    ret <- x$result$betas$se()
    #colnames(ret) <- shortnames(x$model$event_model)
    colnames(ret) <- conditions(x$model$event_model)
    as_tibble(ret, .name_repair="check_unique")
  } else if (type == "contrasts") {
    if (length(x$result$contrasts) == 0) {
      stop("no computed contrasts for this model.")
    }
    ret <- x$result$contrasts$se()
    colnames(ret) <- names(x$bcons)
    as_tibble(ret,.name_repair="check_unique")
  } else if (type == "F") {
    ret <- x$result$Fcontrasts
    ret <- do.call(cbind, lapply(ret, function(f) f$se()))
    colnames(ret) <- names(x$fcons)
    as_tibble(ret,.name_repair="check_unique")
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

fit_lm_contrasts <- function(fit, conlist, fcon, vnames,dofpen=0) {
  conres <- if (!is.null(conlist)) {
    ret <- lapply(conlist, function(con) {
      fit_contrasts(fit, con$weights, attr(con, "term_indices"))
    })
    
    names(ret) <- names(conlist)
    ret
  } 
  
  Fres <- lapply(fcon, function(con) fit_Fcontrasts(fit, t(con), attr(con, "term_indices")))
  names(Fres) <- names(fcon)
  
  bstats <- beta_stats(fit, vnames, dofpen)
  #list(conres=conres, Fres=Fres, bstats=bstats, event_indices=eterm_indices, baseline_indices=bterm_indices)
  list(conres=conres, Fres=Fres, bstats=bstats)
}

#' @keywords internal
multiresponse_lm <- function(form, data_env, conlist, vnames, fcon, modmat=NULL, dofpen=0) {
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
      ret[,event_indices,drop=FALSE]
    }
  }
  
  extract2 <- function(l, el, item, i) {
    #do.call(rbind, lapply(l, function(x) x[[el]][[i]][[item]]()))
    lapply(l, function(x) x[[el]][[i]][[item]]())
  }
  
  do_extract <- function(l, el, items, efun,...) {
    res <- lapply(items, function(item) {
      function() { efun(l, el, item,...) }
    })
  
    #if (is.vector(res)) {
    #  res <- matrix(res, ncol=length(res))
    #}
    names(res) <- items
    res
  }
  
  bstats <- c(do_extract(cres, "bstats", standard_cols, extract), 
              list(stat_type=cres[[1]]$bstats$stat_type))
  
  #cstats <- c(do_extract(cres, "conres", standard_cols, extract2), 
  #            list(stat_type=cres[[1]]$bstats$stat_type))
 

  ncon <- length(cres[[1]]$conres)
  
  conres <- if (ncon >= 1) {
    cdat <- lapply(1:ncon, function(i)  {
      x <- do_extract(cres, "conres", standard_cols, extract2, i)
      c(x, list(stat_type=cres[[1]]$conres[[i]]$stat_type))
    })
    
    names(cdat) <- names(cres[[1]]$conres)
    force(cdat)
    
    format_mat <- function(x, nam) {
      ret <- if (is.list(x)) {
        do.call(cbind, x)
      } else if (nrow(x[[1]]) == 1) {
        do.call(cbind, x)
      } else {
        do.call(cbind, x)
      }
      colnames(ret) <- nam
      ret
    }
    
    ## want format with columns as contrasts, rows as voxels
    ## needs to work when nchunks > 1, when there is 1 series per chunk, when there is 1 series and 1 chunk
    list(
      estimate=function() {
        #x <- do.call(rbind, lapply(cdat, function(x) t(x$estimate())))
        format_mat(lapply(cdat, function(x) unlist(x$estimate())), names(cres[[1]]$conres))
      },
      prob=function() {
        #x <- do.call(rbind, lapply(cdat, function(x) t(x$prob())))
        format_mat(lapply(cdat, function(x) unlist(x$prob())), names(cres[[1]]$conres))
      },
      se=function() {
        format_mat(lapply(cdat, function(x) unlist(x$se())), names(cres[[1]]$conres))
      },
      stat=function() {
        format_mat(lapply(cdat, function(x) unlist(x$stat())), names(cres[[1]]$conres))
      },
      stat_type=cdat[[1]]$stat_type
    )
  }
  
  
  #browser()
    

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
#' @importFrom iterators icount
chunkwise_lm <- function(dset, model, conlist, fcon, nchunks, robust=FALSE, verbose=FALSE, dofpen=0) {
  
  
  chunks <- exec_strategy("chunkwise", nchunks=nchunks)(dset)
  form <- get_formula(model)
  tmats <- term_matrices(model)
  data_env <- list2env(tmats)
  data_env[[".y"]] <- rep(0, nrow(tmats[[1]]))
  modmat <- model.matrix(as.formula(form), data_env)
  
  lmfun <- if (robust) multiresponse_rlm else multiresponse_lm
  
  #browser()
  cres <- foreach( ym = chunks, i = icount(), .verbose=verbose) %dopar% {
    message("processing chunk ", i)
    data_env[[".y"]] <- as.matrix(ym$data)
    ret <- lmfun(form, data_env, conlist, attr(tmats,"varnames"), fcon, modmat=modmat,dofpen=dofpen)
  }
  
  event_indices=attr(tmats, "event_term_indices")
  baseline_indices=attr(tmats, "baseline_term_indices")
  
  wrap_chunked_lm_results(cres, event_indices)

}


#' Run glm for each data run (responses split vertically) and then combine over runs via meta-analysis
#' 
#' @importFrom foreach foreach %do% %dopar%
#' @keywords internal
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
        ##meta_contrasts(conres[hascon]) 
        meta_contrasts(conres) 
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
  
    
    



