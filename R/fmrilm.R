

#' fmri_glm
#' @param formula
#' @param block_formula
#' @param baseline_model
#' @param dataset
#' @param durations
#' @param nuisance_matrix
#' @param drop_empty
#' @param contrasts
#' @param strategy
#' @param 
#' @export
#' @importFrom foreach foreach
fmri_glm <- function(formula, block_formula, baseline_model=NULL, dataset, 
                     durations, drop_empty=TRUE, contrasts=NULL, 
                     strategy=c("runwise", "slicewise", "all")) {
  
 
  strategy <- match.arg(strategy)
  
  assert_that(inherits(dataset, "fmri_dataset"))
  
  if (is.null(baseline_model)) {
    baseline_model <- baseline_model(basis="bs", degree=ceiling(median(dataset$sampling_frame$blocklens)/100), 
                                     sframe=dataset$sampling_frame)
  }
  

  ev_model <- event_model(formula, block_formula, data=dataset$event_table, sampling_frame=dataset$sampling_frame, contrasts=contrasts)
  
  model <- fmri_model(ev_model, baseline_model)
  
  term_names <- names(terms(model))
  term_matrices <- lapply(terms(model), function(x) as.matrix(design_matrix(x)))
  names(term_matrices) <- term_names
  
  form <- as.formula(paste("y ~ ", paste(term_names, collapse = " + "), "-1"))
  
  conlist <- contrast_weights(ev_model)
  
  if (strategy == "runwise") {
    runwise_lm(form, conlist, term_matrices, dset)
  } else {
    stop()
  }
 
  
  model
}
  
runwise_lm <- function(form, conlist, term_matrices, dset) {
    chunks <- exec_strategy("runwise")(dset)
    
    browser()
    cres <- foreach( ym = chunks) %do% {
      data_env <- list2env(lapply(term_matrices, function(x) x[ym$row_ind,]))
      y <- ym$data
      lm.1 <- lm(form, data=data_env)
      
      colind <- attr(conlist, "term_indices")
      conres <- lapply(conlist, function(con) fit_contrasts(lm.1, con, attr(con, "term_indices")))
      names(conres) <- nams(conlist)
      browser()
      fres <- fit_Ftests(lm.1)
      list(fres=fres, conres=conres)
    }
}
  
    
    



