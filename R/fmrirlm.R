#' Fit a Robust Linear Model for fMRI Data Analysis
#'
#' This function fits a robust linear regression model for fMRI data analysis using the specified model formula,
#' block structure, and dataset. The model can be fit using either a runwise or chunkwise data splitting strategy.
#'
#' @inheritParams fmri_lm
#' @param nchunks Number of data chunks when strategy is "chunkwise". Default is 10.
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
fmri_rlm <- function(formula, block, baseline_model = NULL, dataset,
                     durations = 0, drop_empty = TRUE,
                     strategy = c("runwise", "chunkwise"),
                     nchunks = 10,
                     ...) {

  res <- fmri_lm(formula, block, baseline_model = baseline_model, dataset = dataset,
                 durations = durations, drop_empty = drop_empty,
                 strategy = strategy, nchunks = nchunks,
                 robust = TRUE, ...)
  class(res) <- c("fmri_rlm", class(res))
  res
}

#' Perform Chunkwise Robust Linear Modeling on fMRI Dataset
#'
#' @inheritParams chunkwise_lm.fmri_dataset
#' @keywords internal
chunkwise_rlm <- function(dset, model, contrast_objects, nchunks,
                          verbose = FALSE, use_fast_path = FALSE,
                          progress = FALSE) {
  chunkwise_lm.fmri_dataset(dset, model, contrast_objects, nchunks,
                            robust = TRUE, verbose = verbose,
                            use_fast_path = use_fast_path, progress = progress)
}

#' Perform Runwise Robust Linear Modeling on fMRI Dataset
#'
#' @inheritParams runwise_lm
#' @keywords internal
runwise_rlm <- function(dset, model, contrast_objects,
                       verbose = FALSE, use_fast_path = FALSE,
                       progress = FALSE) {
  runwise_lm(dset, model, contrast_objects, robust = TRUE,
             verbose = verbose, use_fast_path = use_fast_path,
             progress = progress)
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
