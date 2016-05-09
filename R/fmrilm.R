
#' @export
#' @importFrom foreach foreach
fmri_glm <- function(formula, dataset, basis=HRF_SPMG1, durations, drop_empty=TRUE, nchunks=1) {
  
  
  model <- fmri_model(formula, dataset$event_table, basis=basis, durations=durations, 
                      blockids=dataset$blockids, blocklens=dataset$blocklens, TR=dataset$TR, 
                      aux_data=dataset$aux_data, drop_empty=drop_empty)
  
  
 

  term_names <- names(terms(model))
  term_matrices <- lapply(terms(model), design_matrix)
  term_matrices <- lapply(term_matrices, as.matrix)
  names(term_matrices) <- term_names
  
  form <- as.formula(paste("ym ~ ", paste(term_names, collapse = " + "), "-1"))
  
  conlist <- contrast_weights(model)
  
  chunks <- data_chunks(dataset, nchunks)
  browser()
  cres <- foreach( ym = chunks) %dopar% {
    lm.1 <- lm(form, data=term_matrices)
    conres <- lapply(conlist, function(con) fit_contrasts(lm.1, con))
    fres <- anova(lm.1)
  }
  
  
  model
}

