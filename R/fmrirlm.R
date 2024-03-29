

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
#' #lm.1 <- fmri_rlm(onset ~ hrf(fac), block= ~ run,dataset=dset)
#' #lm.2 <- fmri_rlm(onset ~ hrf(fac), block= ~ run,dataset=dset2)
#' @keywords experimental
fmri_rlm <- function(formula, block, baseline_model=NULL, dataset, 
                     durations, drop_empty=TRUE, 
                     nchunks=10,
                     strategy=c("runwise", "slicewise", "all")) {
  
 
  strategy <- match.arg(strategy)
  
  
  assert_that(inherits(dataset, "fmri_dataset"))
  ##fobj <- .setup_model(dataset, formula, block, baseline_model, contrasts)
 
  model <- create_fmri_model(formula, block, baseline_model, dataset, durations, drop_empty)
  ret <- fmri_lm_fit(model, dataset, strategy, robust=TRUE, nchunks)
  ret
  
  

}


#' @importFrom foreach foreach %do% %dopar%
#' @autoglobal
runwise_rlm <- function(dset, model, conlist, fcon) {
  
  ## get an iterator of data chunks
  chunks <- exec_strategy("runwise")(dset)
  
  term_names <- names(terms(model))
  form <- paste(".y ~ ", paste(term_names, collapse = " + "), "-1")
  
  ctrl <- lmrob.control(k.max=1000, maxit.scale=1000, max.it=1000)
  ## iterate over each data chunk
  cres <- foreach( ym = chunks) %do% {
    
    ## get event model for the nth run
    eterm_matrices <- lapply(event_terms(model), 
                             function(x) as.matrix(design_matrix(x, ym$chunk_num)))
    
    ## get baseline model for the nth run
    bterm_matrices <- lapply(baseline_terms(model), 
                             function(x) as.matrix(design_matrix(x, ym$chunk_num)))
    
    ## column indices of event regressors
    eterm_indices <- 1:sum(map_int(eterm_matrices, ncol))
    start <- length(eterm_indices) +1
    
    ## column indices of baseline regressors
    bterm_indices <- start:(start+sum(map_int(bterm_matrices, ncol)))
    
    term_matrices <- c(eterm_matrices, bterm_matrices)
    names(term_matrices) <- term_names
    
    data_env <- list2env(term_matrices)
    vnames <- unlist(lapply(term_matrices, colnames))
    
    lapply(1:ncol(ym$data), function(i) {
      data_env[[".y"]] <- ym$data[,i]
      rlm.1 <- lmrob(as.formula(form), data=data_env, control=ctrl)
    
      conres <- lapply(conlist, function(con) fit_contrasts(rlm.1, con, attr(con, "term_indices")))
      names(conres) <- names(conlist)
      
      Fres <- lapply(fcon, function(con) fit_Fcontrasts(rlm.1, t(con), attr(con, "term_indices")))
    
      bstats <- beta_stats(rlm.1, vnames)
      list(contrasts=conres, Fres=Fres, bstats=bstats)
    })
    
  }
  
}


#' Fit a multiresponse ARMA model
#'
#' This function fits a multiresponse ARMA model for fMRI data with specified autocorrelation
#' structure using the provided formula, dataset environment, contrast list, variable names,
#' F-contrast, design matrix, and block ids.
#'
#' @param form The model formula.
#' @param data_env The dataset environment containing the time-series data.
#' @param conlist The contrast list.
#' @param vnames The variable names.
#' @param fcon The F-contrast.
#' @param modmat The design matrix.
#' @param blockids The block ids.
#' @param autocor The autocorrelation structure, either "ar1", "ar2", "arma", or "auto".
#' @return A list containing the contrasts, F-contrasts, and beta statistics for the fitted ARMA model.
#' @importFrom forecast auto.arima
#' @keywords internal
#' @seealso lm.fit, Arima, auto.arima
multiresponse_arma <- function(form, data_env, conlist, vnames, fcon, modmat, 
                              blockids,
                              autocor=c("ar1", "ar2","arma", "auto")) {
  autocor <- match.arg(autocor)
  Y <- data_env$.y
  
  whitened_lm <- function(y, modmat, sq_inv) {
    modmat_wh <- t(sq_inv) %*% modmat
    y_wh <- (y %*% sq_inv)[1,]
    lm.fit(as.matrix(modmat_wh), y_wh)
  }
  

  rlens <- table(blockids)
  i <- NULL
  ret <- foreach(i = 1:ncol(Y)) %dopar% {
    #print(i)
    yi <- Y[,i]
    lm.1 <- lm(yi ~ modmat-1)
    yresid <- resid(lm.1)
    yrsplit <- split(yresid, blockids)
    yresid_padded <- unlist(lapply(yrsplit, function(z) c(z, rep(NA,5))))
    
    fit <- if (autocor == "ar1") {
      afit <- Arima(yresid_padded, order=c(1,0,0), seasonal=FALSE)
      sq_inv <- sq_inv_ar1_by_run(afit$coef[1], rlens)
      whitened_lm(yi, modmat, sq_inv)
    } else if (autocor == "ar2") {
      afit <- Arima(yresid_padded, order=c(2,0,0), seasonal=FALSE)
      #aresid <- resid(afit)
      sq_inv <- sq_inv_arma_by_run(afit$coef[1:2], 0, rlens)
      whitened_lm(yi, modmat, sq_inv)
    } else if (autocor == "arma") {
      afit <- Arima(yresid_padded, order=c(1,0,1), seasonal=FALSE)
      #aresid <- resid(afit)
      sq_inv <- sq_inv_arma_by_run(afit$coef[1], afit$coef[2], rlens)
      whitened_lm(yi, modmat, sq_inv)
    } else if (autocor == "auto") {
      afit <- auto.arima(yresid_padded, max.d=0, max.D=0, seasonal=FALSE)
      nar <- afit$arma[1]
      nma <- afit$arma[2]
      sq_inv <- if (nar == 1 && nma == 0) {
        sq_inv_ar1_by_run(afit$coef[1], rlens)
      } else if (nar == 0 && nma == 1) {
        sq_inv_arma_by_run(0, afit$coef[1], rlens)
      } else if (nar == 0 && nma == 0) {
        Matrix::Diagonal(n=length(yresid))
      } else {
        sq_inv_arma(afit$coef[1:nar], afit$coef[(nar+1):(nar+nma)], rlens)
      }
      
      whitened_lm(yi, modmat, sq_inv)
    }
  
    fit_lm_contrasts(fit, conlist, fcon, vnames)
  }
  

  #ret <- unpack_chunkwise(ret)
  ret
  #conres=ret$conres, Fres=ret$Fres, bstats=ret$bstats
  ## needless name-changes translation...fixme
  #list(contrasts=ret$contrasts, bstats=ret$betas, Fres=ret$Fcontrasts)
}


## could use heavyLm ...
#' @importFrom robustbase lmrob lmrob.control
#' @keywords internal
multiresponse_rlm <- function(form, data_env, conlist, vnames, fcon, modmat=NULL) {
  Y <- data_env$.y
  rcontrol <- lmrob.control(k.max=500, maxit.scale=500)
  ret <- lapply(1:ncol(Y), function(i) {
    data_env[[".y"]] <- Y[,i]
    rlm.1 <- lmrob(as.formula(form), data=data_env,
                   control=rcontrol)
    fit_lm_contrasts(rlm.1, conlist, fcon, vnames)
  })
  ret <- wrap_chunked_lm_results(ret)
  #conres=ret$conres, Fres=ret$Fres, bstats=ret$bstats
  ## needless name-changes translation...fixme
  list(contrasts=ret$contrasts, bstats=ret$betas, Fres=ret$Fcontrasts)
}
