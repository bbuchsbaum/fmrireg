

#' @export
get_formula.fmri_model <- function(x) {
  term_names <- names(terms(x))
  form <- paste(".y ~ ", paste(term_names, collapse = " + "), "-1")
}


#' Construct term matrices for an fMRI model
#'
#' This function constructs the term matrices for an fMRI model, which consists of event-related terms
#' and baseline-related terms. The term matrices are used for building the design matrix in fMRI data analysis.
#'
#' @param x An object of class "fmri_model" containing the event and baseline models.
#' @param blocknum (Optional) A numeric vector specifying the block numbers to be included in the term matrices.
#'                 By default, all unique block numbers in the event model are included.
#' @return A named list of term matrices, with event terms followed by baseline terms.
#'         Attributes "event_term_indices" and "baseline_term_indices" store the indices of event and baseline terms,
#'         "blocknum" stores the block numbers, and "varnames" stores the variable names.
#' @export
#' @seealso fmri_model
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

#' @keywords internal
create_fmri_model <- function(formula, block, baseline_model=NULL, dataset, 
                  durations, drop_empty=TRUE) {
  if (is.null(baseline_model)) {
    baseline_model <- baseline_model(basis="bs", 
                                     degree=max(ceiling(median(dataset$sampling_frame$blocklens)/100),3), 
                                     sframe=dataset$sampling_frame)
  } else {
    assert_that(inherits(baseline_model, "baseline_model"), msg="'baseline_model' arg must have the class 'baseline_model'")
  }
  
  ev_model <- event_model(formula, block, data=dataset$event_table, 
                          sampling_frame=dataset$sampling_frame)
  
  fmri_model(ev_model, baseline_model)
  
}



#' Fit a linear regression model for fMRI data analysis
#'
#' This function fits a linear regression model for fMRI data analysis using the specified model formula,
#' block structure, and dataset. The model can be fit using either a runwise or chunkwise data splitting strategy,
#' and robust fitting can be enabled if desired.
#'
#' @param formula The model formula for experimental events.
#' @param block The model formula for block structure.
#' @param baseline_model (Optional) The \code{baseline_model} object. Default is NULL.
#' @param dataset An object derived from \code{fmri_dataset} containing the time-series data.
#' @param durations A vector of event durations.
#' @param drop_empty Whether to remove factor levels with a size of zero. Default is TRUE.
#' @param robust Whether to use robust fitting. Default is FALSE.
#' @param strategy The data splitting strategy, either "runwise" or "chunkwise". Default is "runwise".
#' @param nchunks Number of data chunks when strategy is `chunkwise`. Default is 10.
#' @param ... Extra arguments.
#' @return A fitted linear regression model for fMRI data analysis.
#' @examples
#' # Example usage of fmri_lm function
#' # ...
#' @export
#' @seealso fmri_dataset, fmri_lm_fit
fmri_lm <- function(formula, block, baseline_model=NULL, dataset, 
                     durations, drop_empty=TRUE, robust=FALSE,
                     strategy=c("runwise", "chunkwise"), nchunks=10, ...) {
  
 
  strategy <- match.arg(strategy)
  assert_that(inherits(dataset, "fmri_dataset"))

  model <- create_fmri_model(formula, block, baseline_model,dataset, durations, drop_empty)
  ret <- fmri_lm_fit(model, dataset, strategy, robust, nchunks)
  ret
}


#' Fit an fMRI linear regression model with a specified fitting strategy
#'
#' This function fits an fMRI linear regression model using the specified fmri_model object, dataset,
#' and data splitting strategy (either "runwise" or "chunkwise"). It is primarily an internal function
#' used by the fmri_lm function.
#'
#' @param fmrimod An object of type \code{fmri_model}.
#' @param dataset An object derived from \code{fmri_dataset} containing the time-series data.
#' @param strategy The data splitting strategy, either "chunkwise" or "runwise". Default is "chunkwise".
#' @param robust Whether to use robust fitting. Default is FALSE.
#' @param nchunks Number of data chunks when strategy is `chunkwise`. Default is 10.
#' @param ... Extra arguments.
#' @return A fitted fMRI linear regression model with the specified fitting strategy.
#' @keywords internal
#' @seealso fmri_lm, fmri_model, fmri_dataset
fmri_lm_fit <- function(fmrimod, dataset, strategy=c("chunkwise", "runwise"), 
                        robust=FALSE, nchunks=10,...) {
  strategy <- match.arg(strategy)
  
  conlist <- unlist(contrast_weights(fmrimod$event_model), recursive=FALSE)
  fcons <- Fcontrasts(fmrimod$event_model)
  
  result <- if (strategy == "runwise") {
    runwise_lm(dataset, fmrimod, conlist, fcons, robust=robust,...)
  } else if (strategy == "chunkwise") {
    chunkwise_lm(dataset, fmrimod, conlist,fcons, nchunks, robust=robust,...)
  }
  
  #browser()
  
  ret <- list(
    result=result,
    model=fmrimod,
    #contrasts=contrasts,
    strategy=strategy,
    fcons=fcons,
    bcons=conlist,
    dataset=dataset)
  
  class(ret) <- "fmri_lm"
  
  ret
  
  
}


#' @export
coef.fmri_lm <- function(object, type=c("estimates", "contrasts"), recon=FALSE, ...) {
  type <- match.arg(type)
  res <- if (type == "estimates") {
    ret <- object$result$betas$estimate[,object$result$event_indices,drop=FALSE]
    colnames(ret) <- conditions(object$model$event_model)
    #shortnames(x$model$event_model)#conditions(x$model$event_model)
    suppressMessages(as_tibble(ret, .name_repair="check_unique"))
  } else if (type == "contrasts") {
    if (!is.null(object$result$contrasts$estimate)) {
      ret <- object$result$contrasts$estimate
      colnames(ret) <- names(object$bcons)
      suppressMessages(as_tibble(ret, .name_repair="check_unique"))
    } else {
      stop("no contrasts for this model.")
    }
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
stats.fmri_lm <- function(x, type=c("estimates", "contrasts", "F"), ...) {
  type <- match.arg(type)
  if (type == "estimates") {
    ret <- x$result$betas$estimate[,x$result$event_indices,drop=FALSE]/x$result$betas$se[,x$result$event_indices,drop=FALSE]
    #colnames(ret) <- shortnames(x$model$event_model)
    colnames(ret) <- conditions(x$model$event_model)
    suppressMessages(as_tibble(ret,.name_repair="check_unique"))
  } else if (type == "contrasts") {
    if (length(x$result$contrasts) == 0) {
      stop("no computed contrasts for this model.")
    }
    ret <- x$result$contrasts$estimate/x$result$contrasts$se
    colnames(ret) <- names(x$bcons)
    suppressMessages(as_tibble(ret, .name_repair="check_unique"))
  } else if (type == "F") {
    if (is.null(x$result$Fcontrasts)) {
      stop("no computed F contrasts for this model.")
    }
    ret <- x$result$Fcontrasts$stat
    #ret <- do.call(cbind, lapply(ret, function(f) f$stat()))
    colnames(ret) <- names(x$fcons)
    suppressMessages(as_tibble(ret, .name_repair="check_unique"))
  }
}

#' @export
standard_error.fmri_lm <- function(x, type=c("estimates", "contrasts")) {
  type <- match.arg(type)
  if (type == "estimates") {
    ret <- x$result$betas$se
    #colnames(ret) <- shortnames(x$model$event_model)
    colnames(ret) <- conditions(x$model$event_model)
    suppressMessages(as_tibble(ret, .name_repair="check_unique"))
  } else if (type == "contrasts") {
    if (length(x$result$contrasts) == 0) {
      stop("no computed contrasts for this model.")
    }
    ret <- x$result$contrasts$se
    colnames(ret) <- names(x$bcons)
    suppressMessages(as_tibble(ret,.name_repair="check_unique"))
  } else if (type == "F") {
    ret <- x$result$Fcontrasts$se
    #ret <- do.call(cbind, lapply(ret, function(f) f$se()))
    colnames(ret) <- names(x$fcons)
    suppressMessages(as_tibble(ret,.name_repair="check_unique"))
  }
}


  

#' @export
print.fmri_lm <- function(x,...) {
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


#' @keywords internal
fit_lm_contrasts <- function(fit, conlist, fcon, vnames, se=TRUE) {
  conres <- if (!is.null(conlist)) {
    ret <- lapply(conlist, function(con) {
      fit_contrasts(fit, con$weights, attr(con, "term_indices"), se=se)
    })
    
    names(ret) <- names(conlist)
    ret
  } 
  
  Fres <- lapply(fcon, function(con) fit_Fcontrasts(fit, t(con), attr(con, "term_indices")))
  names(Fres) <- names(fcon)
  
  bstats <- beta_stats(fit, vnames,se=se)
  #list(conres=conres, Fres=Fres, bstats=bstats, event_indices=eterm_indices, baseline_indices=bterm_indices)
  list(contrasts=conres, Fres=Fres, bstats=bstats, fit=fit)
  #list(conres=conres, Fres=Fres, bstats=bstats)
}




#' Multiresponse Linear Model
#'
#' This function fits a linear model to multiple responses in an fMRI dataset.
#'
#' @param form The formula used to define the linear model.
#' @param data_env The environment containing the data to be used in the linear model.
#' @param conlist The list of contrasts used in the analysis.
#' @param vnames The names of the variables used in the linear model.
#' @param fcon The F-contrasts used in the analysis.
#' @param modmat The model matrix (default is NULL, which will calculate the model matrix using the formula).
#' @return A list containing the results from the multiresponse linear model analysis.
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
  
  
  nf <- length(cres[[1]]$Fres)
  
  Fres <- if (nf >= 1) {
    ret <- lapply(1:nf, function(i)  {
      x <- do_extract(cres, "Fres", standard_cols,extract2, i)
      c(x, list(stat_type=cres[[1]]$Fres[[i]]$stat_type))
    })
    names(ret) <- names(cres[[1]]$Fres)
    ret
  }
  
  conmats <- lapply(cres[[1]]$conres, function(x) x$conmat)
  list(betas=bstats, contrasts=conres, Fcontrasts=Fres, conmats=conmats)
    
      
}

#' Unpack Chunkwise Results
#'
#' This function processes and unpacks the results of chunkwise analysis in fMRI
#' data. It extracts and reorganizes the effects, standard errors, contrasts, 
#' F-contrasts, and contrast matrices from the chunkwise analysis results.
#'
#' @param cres The results of the chunkwise analysis.
#' @param event_indices The indices of the event-related effects.
#' @param baseline_indices The indices of the baseline-related effects.
#' @return A list containing the betas, contrasts, F-contrasts, and contrast matrices
#'         for the unpacked chunkwise results.
#' @keywords internal
#' @noRd
unpack_chunkwise <- function(cres, event_indices, baseline_indices) {
 
  effects <- do.call(rbind, lapply(cres, function(x) x$bstats$estimate))
  effects_se <- do.call(rbind, lapply(cres, function(x) x$bstats$se))
  
  betas <- list(
    estimate=effects[,event_indices,drop=FALSE],
    se=effects_se[,event_indices,drop=FALSE]
  )
  
  ncon <- length(cres[[1]]$contrasts)
  nf <- length(cres[[1]]$Fres)
  
  ## helper for Fres
  extract_values <- function(Fnames, name) {
    #### used to be 'rbind', why? ret <- do.call(rbind, lapply(Fnames, function(fn) {
    ret <- do.call(cbind, lapply(Fnames, function(fn) {
      do.call(rbind, lapply(cres, function(x) as.matrix(x$Fres[[fn]][[name]])))
    }))
    colnames(ret) <- Fnames
    ret
  }
  
  
  contrasts <- if (ncon > 0) {
    counter= 1
    con_effects <- do.call(rbind, lapply(cres, function(x) {
      #print(x)
      do.call(cbind, lapply(x$contrasts, function(z) {
        z$estimate
      }))
    }))
    
  
    con_se <- do.call(rbind, lapply(cres, function(x) {
      do.call(cbind, lapply(x$contrasts, function(z) z$se))
    }))
    
    
    list(estimate=con_effects,
         se=con_se)
  }
  
  #browser()
  #browser()
  Fres <- if (nf >= 1) {
    Fnames <- names(cres[[1]]$Fres)
    list(estimate=extract_values(Fnames, "estimate"),
         se=extract_values(Fnames, "se"),
         stat=extract_values(Fnames, "stat"),
         prob=extract_values(Fnames, "prob"))
  }
  
  
  conmats <- lapply(cres[[1]]$contrasts, function(x) x$conmat)
  #print("made past conmats")
  list(betas=betas, contrasts=contrasts, Fcontrasts=Fres, conmats=conmats)
}


#' Chunkwise Linear Model for fMRI Dataset
#'
#' This function performs a chunkwise linear model analysis on fMRI dataset, 
#' splitting the dataset into chunks and running the linear model on each chunk.
#'
#' @param dset An object of class \code{fmri_dataset}.
#' @param model The fMRI model used for the analysis.
#' @param conlist The list of contrasts used in the analysis.
#' @param fcon The F-contrasts used in the analysis.
#' @param nchunks The number of chunks to divide the dataset into.
#' @param robust Whether to use robust linear modeling (default is FALSE).
#' @param verbose Whether to display progress messages (default is FALSE).
#' @return A list containing the unpacked chunkwise results.
#' @keywords internal
#' @importFrom iterators icount
#' @autoglobal 
chunkwise_lm.fmri_dataset <- function(dset, model, conlist, fcon, nchunks, robust=FALSE, verbose=FALSE) {
  chunks <- exec_strategy("chunkwise", nchunks=nchunks)(dset)
  form <- get_formula(model)
  tmats <- term_matrices(model)
  data_env <- list2env(tmats)
  data_env[[".y"]] <- rep(0, nrow(tmats[[1]]))
  modmat <- model.matrix(as.formula(form), data_env)
  
  lmfun <- if (robust) multiresponse_rlm else multiresponse_lm

  ym <- NULL
  cres <- foreach( ym = chunks, i = icount(), .verbose=verbose) %dopar% {
    message("processing chunk ", i)
    data_env[[".y"]] <- as.matrix(ym$data)
    ret <- lmfun(form, data_env, conlist, attr(tmats,"varnames"), fcon, modmat=modmat)
  }
  
  event_indices=attr(tmats, "event_term_indices")
  baseline_indices=attr(tmats, "baseline_term_indices")
  
  out <- unpack_chunkwise(cres,event_indices,baseline_indices) %>% purrr::list_modify(event_indices=event_indices,
                                                                                  baseline_indices=baseline_indices)
  out
  
}


#' Runwise Linear Model for fMRI Dataset
#'
#' This function performs a runwise linear model analysis on an fMRI dataset by 
#' running the linear model for each data run (responses split vertically) and 
#' then combines the results over runs via meta-analysis.
#'
#' @param dset An object of class \code{fmri_dataset}.
#' @param model The fMRI model used for the analysis.
#' @param conlist The list of contrasts used in the analysis.
#' @param fcon The F-contrasts used in the analysis.
#' @param robust Whether to use robust linear modeling (default is FALSE).
#' @param verbose Whether to display progress messages (default is FALSE).
#' @return A list containing the combined results from runwise linear model analysis.
#' @importFrom foreach foreach %do% %dopar%
#' @keywords internal
#' @autoglobal
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

      list(conres=ret$contrasts, Fres=ret$Fres, bstats=ret$bstats, 
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
      list(contrasts=meta_con, betas=meta_beta, Fcontrasts=meta_F,
           event_indices=cres[[1]]$event_indices, baseline_indices=cres[[1]]$baseline_indices)
    } else {
      list(contrasts=conres[[1]], betas=bstats[[1]], Fcontrasts=Fres[[1]],
           event_indices=cres[[1]]$event_indices, baseline_indices=cres[[1]]$baseline_indices)
    }
    
}
  
    
    



