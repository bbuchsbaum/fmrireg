

#' @param formula
#' @param baseline
#' @param event_table
#' @param aux_table
#' @param basis
#' @param durations
#' @param sampling_frame
#' @export
fmri_model <- function(formula, baseline, sampling_frame, event_table, event_block_ids, 
                       aux_table=data.frame(), basis=HRF_SPMG1,
                       drop_empty=TRUE, durations=0) {
  
  stopifnot(inherits(formula, "formula"))
  
  assert_that(all(sampling_frame$blocklens>0))
  
  assert_that(length(sampling_frame$TR) == 1)
  
  assert_that(sampling_frame$TR > 0)
  
  assert_that(is.data.frame(aux_table))
  
  assert_that(length(event_block_ids) == nrow(event_table))
  
  #if (is.null(resp)) {
  #  stop("need to provide onset vector on left side of formula, e.g. Onsets ~  a + b")
  #}
  
  if (missing(durations)) {
    ## assume zero-duration impulse for all events
    durations <- rep(0, nrow(event_table))
  }
  
  formspec <- function(formula, table) {
    vterms <- extract_terms(formula, table)
    resp <- attr(vterms, "response")
    variables <- extract_variables(formula, table)
    
    if (resp != 0) {
      lhs <- variables[[resp]]
    } else {
      lhs <- NULL
    }
    
    rhs <- variables[(resp+1):length(variables)]
    vclass <- sapply(rhs, class)
    
    ret <- list(vterms=vterms, resp=resp, variables=variables, lhs=lhs, rhs=rhs,vclass=vclass)
    
    class(ret) <- "formula_extraction"
    ret
  }
  
  
  event_spec <- formspec(formula, event_table)
  
  # if (missing(event_block_ids) || is.null(event_block_ids)) {
  #   which_block <- which(sapply(event_spec$rhs, function(x) inherits(x, "blockspec")))
  #   block_spec <- event_spec$rhs[which_block]
  #   event_spec$rhs <- event_spec$rhs[-which_block]
  # } else {
  #   blockspec <- block(event_block_ids)
  # }
  # 
  # block_term <- construct(blockspec)
  
  baseline_spec <- formspec(baseline, aux_table)
  
  model_spec <- list(formula=formula, baseline_spec=baseline_spec, event_table=event_table, aux_table=aux_table,
                     onsets=event_spec$lhs, event_spec=event_spec, event_block_ids=event_block_ids, baseline_spec=baseline_spec, 
                     durations=durations, sampling_frame=sampling_frame, drop_empty=drop_empty)
  
  class(model_spec) <- c("model_spec", "list")
  
  fmodel <- construct_model(model_spec)
  fmodel
}



construct_model <- function(x) {
  
  terms <- lapply(x$event_spec$rhs, function(m) construct(m,x))
  term_names <- sapply(x$varspec, "[[", "label")
  term_names <- .sanitizeName(term_names)
  names(terms) <- term_names
  
  term_lens <- sapply(lapply(terms, conditions), length)
  spans <- c(0, cumsum(term_lens))
  
  term_indices <- lapply(1:(length(spans)-1), function(i) {
    seq(spans[i]+1, spans[i+1])
  })
  
  names(term_indices) <- term_names
  
  ret <- list(
    term_indices=term_indices,
    terms=terms,
    model_spec=x
  )
  
  class(ret) <- c("fmri_model", "list")
  ret
}

