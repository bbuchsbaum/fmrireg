#' Fit a Robust Linear Model for fMRI Data Analysis
#'
#' This function fits a robust linear regression model for fMRI data analysis using the specified model formula,
#' block structure, and dataset. The model can be fit using either a runwise or chunkwise data splitting strategy.
#'
#' @inheritParams fmri_lm
#' @param nchunks Number of data chunks when strategy is "chunkwise". Default is 10.
#' @param meta_weighting Method for combining results across runs/chunks. Either "inv_var" for inverse variance weighting or "equal" for equal weighting. Default is "inv_var".
#' @return A fitted robust linear regression model for fMRI data analysis.
#' @export
#' @examples
#' etab <- data.frame(onset=c(1,30,15,25), fac=factor(c("A", "B", "A", "B")), run=c(1,1,2,2))
#' etab2 <- data.frame(onset=c(1,30,65,75), fac=factor(c("A", "B", "A", "B")), run=c(1,1,1,1))
#' mat <- matrix(rnorm(100*100), 100,100)
#' dset <- matrix_dataset(mat, TR=1, run_length=c(50,50),event_table=etab)
#' dset2 <- matrix_dataset(mat, TR=1, run_length=c(100),event_table=etab2)
#' lm.1 <- fmri_rlm(onset ~ hrf(fac), block= ~ run, dataset=dset)
#' lm.2 <- fmri_rlm(onset ~ hrf(fac), block= ~ run, dataset=dset2)
fmri_rlm <- function(formula, block, baseline_model=NULL, dataset, 
                     durations=0, drop_empty=TRUE,
                     strategy=c("runwise", "chunkwise"), 
                     nchunks=10,
                     meta_weighting=c("inv_var", "equal"),
                     ...) {
  
  strategy <- match.arg(strategy)
  meta_weighting <- match.arg(meta_weighting)
  
  # Error checking
  assert_that(is.formula(formula), msg = "'formula' must be a formula")
  assert_that(is.formula(block), msg = "'block' must be a formula")
  assert_that(inherits(dataset, "fmri_dataset"), msg = "'dataset' must be an 'fmri_dataset'")
  assert_that(is.numeric(durations), msg = "'durations' must be numeric")
  if (strategy == "chunkwise") {
    assert_that(is.numeric(nchunks) && nchunks > 0, msg = "'nchunks' must be a positive number")
  }
  
  model <- create_fmri_model(formula, block, baseline_model, dataset, durations, drop_empty)
  ret <- fmri_lm_fit(model, dataset, strategy, robust=TRUE, nchunks=nchunks, 
                     meta_weighting=meta_weighting, ...)
  
  class(ret) <- c("fmri_rlm", "fmri_lm")
  return(ret)
}

#' Perform Chunkwise Robust Linear Modeling on fMRI Dataset
#'
#' @inheritParams chunkwise_lm.fmri_dataset
#' @keywords internal
chunkwise_rlm <- function(dset, model, conlist, fcon, nchunks, verbose = FALSE) {
  chunks <- exec_strategy("chunkwise", nchunks = nchunks)(dset)
  form <- get_formula(model)
  tmats <- term_matrices(model)
  data_env <- list2env(tmats)
  modmat <- model.matrix(as.formula(form), data_env)
  
  ctrl <- robustbase::lmrob.control(k.max=500, maxit.scale=500)
  
  cres <- foreach(ym = chunks, .verbose = verbose) %dopar% {
    if (verbose) message("Processing chunk ", ym$chunk_num)
    
    results <- lapply(1:ncol(ym$data), function(i) {
      data_env$.y <- ym$data[,i]
      fit <- robustbase::lmrob(as.formula(form), data = data_env, control = ctrl)
      fit_lm_contrasts(fit, conlist, fcon, attr(tmats, "varnames"))
    })
    
    # Combine results for this chunk
    bstats <- lapply(results, function(x) x$bstats) %>% dplyr::bind_rows()
    contrasts <- lapply(results, function(x) x$contrasts) %>% dplyr::bind_rows()
    
    list(
      bstats = bstats,
      contrasts = contrasts,
      event_indices = attr(tmats, "event_term_indices"),
      baseline_indices = attr(tmats, "baseline_term_indices")
    )
  }
  
  unpack_chunkwise(cres, attr(tmats, "event_term_indices"), 
                   attr(tmats, "baseline_term_indices"))
}

#' Perform Runwise Robust Linear Modeling on fMRI Dataset
#'
#' @inheritParams runwise_lm
#' @keywords internal
runwise_rlm <- function(dataset, model, conlist, fcons,
                       meta_weighting=c("inv_var", "equal"), 
                       verbose=FALSE, ...) {
  # Set up robust control parameters
  ctrl <- robustbase::lmrob.control(k.max=500, maxit.scale=500)
  
  runwise_fit(dataset, model, conlist, fcons,
             fit_fun=multiresponse_rlm,
             meta_weighting=meta_weighting,
             verbose=verbose,
             control=ctrl, ...)
}

#' Fit a multiresponse robust linear model
#'
#' This function fits a robust linear model to multiple responses in an fMRI dataset.
#' Unlike the standard multiresponse_lm, this processes one response at a time due
#' to limitations of robust regression implementations.
#'
#' @param form The formula used to define the linear model.
#' @param data_env The environment containing the data to be used in the linear model.
#' @param conlist The list of contrasts used in the analysis.
#' @param vnames The names of the variables used in the linear model.
#' @param fcon The F-contrasts used in the analysis.
#' @param modmat The model matrix (default is NULL, which will calculate the model matrix using the formula).
#' @return A list containing the results from the multiresponse robust linear model analysis.
#' @keywords internal
#' @noRd
multiresponse_rlm <- function(form, data_env, conlist, vnames, fcon, modmat = NULL) {
  with_package("robustbase")
  Y <- data_env$.y
  ctrl <- robustbase::lmrob.control(k.max=500, maxit.scale=500)
  
  # Process each response variable separately
  results <- lapply(1:ncol(Y), function(i) {
    data_env$.y <- Y[,i]
    
    # Fit robust model
    fit <- if (is.null(modmat)) {
      robustbase::lmrob(as.formula(form), data = data_env, control = ctrl)
    } else {
      robustbase::lmrob(modmat, Y[,i], control = ctrl)
    }
    
    # Get contrasts and stats
    fit_lm_contrasts(fit, conlist, fcon, vnames)
  })
  
  # Combine results
  bstats <- lapply(results, function(x) x$bstats) %>% dplyr::bind_rows()
  contrasts <- lapply(results, function(x) x$contrasts) %>% dplyr::bind_rows()
  
  list(
    contrasts = contrasts,
    bstats = bstats
  )
}

#' Print method for fmri_rlm objects
#' 
#' @param x An fmri_rlm object
#' @param ... Additional arguments passed to print.fmri_lm
#' @export
print.fmri_rlm <- function(x, ...) {
  cat("\n═══ Robust fMRI Linear Model ═══\n")
  NextMethod()
}
