
#' @export
fmri_glm <- function(formula, dataset, basis=HRF_SPMG1, durations, drop_empty=TRUE, nchunks=1) {
  
  
  model <- fmri_model(formula, dataset$event_table, basis=basis, durations=durations, 
                      blockids=dataset$blockids, blocklens=dataset$blocklens, TR=dataset$TR, 
                      aux_data=dataset$aux_data, drop_empty=drop_empty)
  
  
 
 
  
  term_names <- .sanitizeName(names(terms(model)))
  term_matrices <- lapply(terms(model), design_matrix)
  names(term_matrices) <- term_names
  
  form <- as.formula(paste("y ~ ", paste(term_names, collapse = " + ")))
  
  chunks <- data_chunks(dataset, nchunks)
  
  
  
  model
}

