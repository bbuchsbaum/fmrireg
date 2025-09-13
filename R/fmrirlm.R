#' Fit a Robust Linear Model for fMRI Data Analysis
#'
#' This function fits a robust linear regression model for fMRI data analysis using the specified model formula,
#' block structure, and dataset. The model can be fit using either a runwise or chunkwise data splitting strategy.
#'
#' @inheritParams fmri_lm
#' @param nchunks Number of data chunks when strategy is "chunkwise". Default is 10.
#' @param cor_struct Correlation structure: "iid", "ar1", "ar2", or "arp". Default is "iid".
#' @param cor_iter Number of iterations for AR parameter estimation. Default is 1.
#' @param cor_global Whether to use global AR parameters. Default is FALSE.
#' @param ar_p Order of autoregressive model for "arp" structure. Default is NULL.
#' @param ar1_exact_first Whether to use exact AR(1) for first iteration. Default is FALSE.
#' @param robust_psi Robust psi function: "huber" or "bisquare". Default is "huber".
#' @param robust_k_huber Tuning constant for Huber's psi. Default is 1.345.
#' @param robust_c_tukey Tuning constant for Tukey's bisquare. Default is 4.685.
#' @param robust_max_iter Maximum iterations for robust fitting. Default is 2.
#' @param robust_scale_scope Scope for robust scale estimation: "run" or "global". Default is "run".
#' @param ... Additional arguments passed to fmri_lm
#' @return A fitted robust linear regression model for fMRI data analysis.
#' @export
#' @examples
#' etab <- data.frame(onset=c(1,30,15,25), fac=factor(c("A", "B", "A", "B")), run=c(1,1,2,2))
#' etab2 <- data.frame(onset=c(1,30,65,75), fac=factor(c("A", "B", "A", "B")), run=c(1,1,1,1))
#' mat <- matrix(rnorm(100*100), 100,100)
#' dset <- fmridataset::matrix_dataset(mat, TR=1, run_length=c(50,50),event_table=etab)
#' dset2 <- fmridataset::matrix_dataset(mat, TR=1, run_length=c(100),event_table=etab2)
#' lm.1 <- fmri_rlm(onset ~ hrf(fac), block= ~ run, dataset=dset)
#' lm.2 <- fmri_rlm(onset ~ hrf(fac), block= ~ run, dataset=dset2)
fmri_rlm <- function(formula, block, baseline_model = NULL, dataset,
                     durations = 0, drop_empty = TRUE,
                     strategy = c("runwise", "chunkwise"),
                     nchunks = 10,
                     cor_struct = c("iid", "ar1", "ar2", "arp"), cor_iter = 1L,
                     cor_global = FALSE, ar_p = NULL, ar1_exact_first = FALSE,
                     robust_psi = c("huber", "bisquare"), robust_k_huber = 1.345,
                     robust_c_tukey = 4.685, robust_max_iter = 2L,
                     robust_scale_scope = c("run", "global"), ...) {

  # Create robust options
  robust_options <- list(
    type = match.arg(robust_psi),
    k_huber = robust_k_huber,
    c_tukey = robust_c_tukey,
    max_iter = robust_max_iter,
    scale_scope = match.arg(robust_scale_scope)
  )
  
  # Create AR options
  ar_options <- list(
    cor_struct = match.arg(cor_struct),
    p = ar_p,
    iter = cor_iter,
    global = cor_global,
    exact_first = ar1_exact_first
  )

  res <- fmri_lm(formula, block, baseline_model = baseline_model, dataset = dataset,
                 durations = durations, drop_empty = drop_empty,
                 strategy = strategy, nchunks = nchunks,
                 robust = TRUE,
                 robust_options = robust_options,
                 ar_options = ar_options,
                 ...)
  class(res) <- c("fmri_rlm", class(res))
  res
}

# Note: chunkwise_rlm is no longer needed with the new interface
# The fmri_lm function handles robust fitting directly

# Note: runwise_rlm is no longer needed with the new interface
# The fmri_lm function handles robust fitting directly

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
#' @rdname print
print.fmri_rlm <- function(x, ...) {
  cat("\n=== Robust fMRI Linear Model ===\n")
  NextMethod()
}
