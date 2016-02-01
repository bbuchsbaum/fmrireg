
#' @export
fmri_glm <- function(formula, dataset, basis=HRF_SPMG1, durations, drop_empty=TRUE) {
  
  
  model <- fmri_model(formula, dataset$event_table, basis=basis, durations=durations, 
                      blockids=dataset$blockids, blocklens=dataset$blocklens, TR=dataset$TR, 
                      aux_data=dataset$aux_data, drop_empty=drop_empty)
  
  ret <- list(
    model_spec=model,
    convolved_model=construct(model),
    dataset=dataset
  )
  
  class(ret) <- c("fmri_glm_spec", "list")
  ret
}

