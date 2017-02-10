

#' @importFrom lazyeval f_eval f_rhs f_lhs
#' @param formula
#' @param data
#' @param block
#' @param sampling_frame
#' @param drop_empty
#' @param durations
#' @param contrasts
#' @export
event_model <- function(formula, data, block, sampling_frame, drop_empty=TRUE, durations=0, contrasts=NULL) {
  
  stopifnot(inherits(formula, "formula"))

  if (lazyeval::is_formula(block)) {
    block_rhs <- lazyeval::f_rhs(block)
    blockids <- lazyeval::f_eval_rhs(block, data)
  } else {
    blockids <- block
  }
  
  assert_that(length(blockids) == nrow(data))
  
  if (missing(durations)) {
    ## assume zero-duration impulse for all events
    durations <- rep(0, nrow(data))
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
  
  
  event_spec <- formspec(formula, data)
  
 
  model_spec <- list(formula=formula, event_table=data, onsets=event_spec$lhs, event_spec=event_spec, 
                     blockids=blockids, 
                     durations=durations, 
                     sampling_frame=sampling_frame,
                     drop_empty=drop_empty,
                     contrasts=contrasts)
  
  class(model_spec) <- c("event_model_spec", "list")
  
  fmodel <- construct_model(model_spec)
  fmodel
}


construct_model <- function(x) {

  terms <- lapply(x$event_spec$rhs, function(m) construct(m,x))
  term_names <- sapply(x$event_spec$rhs, "[[", "label")
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
    blockids=x$blockids,
    sampling_frame=x$sampling_frame,
    model_spec=x
  )
  
  class(ret) <- c("event_model", "list")
  ret
}

extract_terms <- function(formula, data) {
  if (!inherits(formula, "terms")) {
    terms(formula, data = data)
  } else {
    formula
  }	
}

extract_covariates <- function(.terms, variables, resp, etab) {
  vars <- attr(.terms, "variables") 
  varnames <- sapply(vars, deparse, width.cutoff = 500)[-1]
  ind.vars <- varnames[-resp] 
  orig.covar.names <- ind.vars[which(sapply(variables[-resp], function(obj) is.numeric(obj) || .is.parametric.basis(obj)))]
  new.covar.names <- names(events(etab))[sapply(events(etab), is_continuous)]
  covar.names <- as.list(orig.covar.names)
  names(covar.names) <- new.covar.names
  covar.names
}

is_parametric_basis <- function(obj) { inherits(obj, "ParametricBasis") }


#' @importFrom formula.tools rhs lhs op op.type
extract_variables <- function(form, data) {
 
  .terms <- extract_terms(form,data)
  env <- environment(.terms)
  #env <- environment(.terms)
  varnames <- sapply(attr(.terms, "variables") , deparse, width.cutoff = 500)[-1]
  #varnames <- c(formula.tools::get.vars(formula.tools::lhs(form)), formula.tools::get.vars(formula.tools::rhs(form)))
  #variables <- as.data.frame(do.call(cbind, lapply(varnames, function(vn) eval(parse(text=vn), data,environment(terms(form))))))
  variables <- eval(attr(.terms, "variables"), data, env) 
  names(variables) <- varnames
  variables
}

# has_block_variable <- function(form) {
#   ret <- op(rhs(form)) == "|"
#   length(ret) > 0 && ret
# }


parse_term <- function(vars, ttype) {
  dim <- length(vars) # number of variables
  term <- deparse(vars[[1]],backtick=TRUE) # first covariate
  if (dim>1) # then deal with further covariates
    for (i in 2:dim) term[i]<-deparse(vars[[i]],backtick=TRUE)
  for (i in 1:dim) term[i] <- attr(terms(reformulate(term[i])),"term.labels")
  
  full.call<-paste(ttype, "(",term[1],sep="")
  if (dim > 1) for (i in 2:dim) full.call<-paste(full.call,",",term[i],sep="")
  label <- paste(full.call,")",sep="")   # label for parameters of this term  
  
  list(term=term, label=label)
}



#' @importFrom tibble as_tibble
design_matrix.event_model_spec <- function(x) {
  termlist <- lapply(x$varspec, function(m) construct(m,x))
  ret <- lapply(termlist, design_matrix)
  vnames <- unlist(lapply(ret, colnames))
  dmat <- tibble::as_tibble(do.call(cbind, ret))
  names(dmat) <- vnames
  dmat
}

#' @importFrom tibble as_tibble
#' @rdname design_matrix
design_matrix.event_model <- function(x, blockid=NULL) {
  ret <- lapply(x$terms, design_matrix, blockid)
  vnames <- unlist(lapply(ret, names))
  dmat <- tibble::as_tibble(do.call(cbind, ret))
  names(dmat) <- vnames
  dmat
}



#' @export
terms.event_model <- function(x) {
  x$terms
}

#' @export
#' @rdname conditions
conditions.event_model <- function(x) {
  unlist(lapply(terms(x), function(t) conditions(t)), use.names=FALSE)
}

#' @export
#' @rdname contrast_weights
contrast_weights.convolved_term <- function(x) {
  lapply(x$contrasts, function(cspec) {
    if (!is.null(cspec))
      contrast_weights(cspec,x)
  })
}

#' @export
#' @rdname contrast_weights
FContrasts.convolved_term <- function(x) {
  Fcontrasts(x$evterm)
}


#' @export
#' @rdname contrast_weights
contrast_weights.fmri_term <- function(x) { stop("unimplemented") }

contrast_weights.fmri_model <- function(x) { contrast_weights(x$event_model) }

#' @export
#' @rdname contrast_weights
contrast_weights.event_model <- function(x) {
  tind <- x$term_indices
  len <- length(conditions(x))
  tnames <- names(terms(x))
  ret <- unlist(lapply(seq_along(terms(x)), function(i) {
    cwlist <- contrast_weights(terms(x)[[i]])
    if (!is.null(cwlist)) {
      ret <- lapply(cwlist, function(cw) {
        out <- numeric(len)
        out[tind[[i]]] <- as.vector(cw$weights)
        attr(out, "term_indices") <- as.vector(tind[[i]])
        out
      })
      
      cnames <- sapply(cwlist, function(cw) cw$name)
    
      prefix <- tnames[i]
      names(ret) <- paste0(prefix, "#", cnames)
      
      ret
    }
  }), recursive=FALSE)

  
  ret
}

Fcontrasts.event_model <- function(x) {
  tind <- x$term_indices
  len <- length(conditions(x))
  tnames <- names(terms(x))
  ret <- unlist(lapply(seq_along(terms(x)), function(i) {
    cwlist <- Fcontrasts(terms(x)[[i]]$evterm)
    if (!is.null(cwlist)) {
      ret <- lapply(cwlist, function(cw) {
        out <- matrix(0, len, ncol(cw))
        out[tind[[i]],] <- cw
        attr(out, "term_indices") <- as.vector(tind[[i]])
        row.names(out) <- row.names(cw)
        out
      })
      
      cnames <- names(cwlist)
      prefix <- tnames[i]
      names(ret) <- paste0(prefix, "#", cnames)
      ret
    }
  }), recursive=FALSE)
}
  
  
#' @export
#' @rdname contrast_weights
design_matrix.convolved_term <- function(x, blockid=NULL) {
  if (is.null(blockid)) {
    x$design_matrix
  } else {
    keep <- blockids(x$sampling_frame) %in% blockid
    x$design_matrix[keep,]
  } 
}

#' @importFrom tibble as_tibble
#' @export
#' @rdname contrast_weights
matrix_term <- function(varname, mat) {
  stopifnot(is.matrix(mat))
  ret <- list(varname=varname, design_matrix=tibble::as_tibble(mat))
  class(ret) <- c("matrix_term", "fmri_term", "list")
  ret
}


#' @export
#' @rdname design_matrix
design_matrix.matrix_term <- function(x) {
  if (is.null(names(x$design_matrix))) {
    cnames <- paste0(x$varname, "_", 1:ncol(x$design_matrix))
    names(x$design_matrix) <- cnames
  }
  
  x$design_matrix

}

#' @export
#' @rdname event_table
event_table.convolved_term <- function(x) event_table(x$evterm)

#' @export
#' @rdname nbasis
nbasis.convolved_term <- function(x) nbasis(x$hrf)

#' @export
#' @rdname longnames
longnames.convolved_term <- function(x) {
  # ignores exclude.basis
  term.cells <- cells(x)
  # ignores exclude.basis
  apply(as.matrix(sapply(1:ncol(term.cells), 
    function(i) {
      paste0(names(term.cells)[i], "#", term.cells[[i]], sep="")
  })), 1, paste, collapse=":")

}

#' @export
#' @rdname longnames
shortnames.convolved_term <- function(x) {
  # ignores exclude.basis
  term.cells <- cells(x)
  # ignores exclude.basis
  apply(as.matrix(sapply(1:ncol(term.cells), 
                         function(i) {
                           term.cells[[i]]
                         })), 1, paste, collapse=":")
  
}


#' @export
#' @rdname longnames
longnames.matrix_term <- function(x) {
  paste0(x$name, "#", colnames(design_matrix(x)))
}


#' @export
print.event_model <- function(object) {
  cat("event_model", "\n")
  cat(" ", Reduce(paste, deparse(object$model_spec$formula)), "\n")
  cat(" ", "Num Terms", length(terms(object)), "\n")
  cat(" ", "Num Events: ", nrow(object$model_spec$event_table), "\n")
  cat(" ", "Num Columns: ", length(conditions(object)), "\n")
  cat(" ", "Num Blocks: ", length(unique(object$blockids)), "\n")
  #cat(" ", "Length of Blocks: ", paste(object$model_spec$blocklens, collapse=", "), "\n")
  for (i in 1:length(terms(object))) {
    cat("\n")
    t <- terms(object)[[i]]
    cat("Term:", i, " ")
    print(t)
    cat("\n")
  }
  
}

#' @importFrom ggplot2 ggplot aes_string geom_line facet_wrap xlab theme_bw
#' @importFrom tidyr gather
#' @export
plot.event_model <- function(x, term_name=NULL, longnames=TRUE) {
  all_terms <- terms(x)
  term_names <- names(all_terms)
  
  sframe <- x$sampling_frame
  
  dflist <- lapply(all_terms, function(term) {
    dm1 <- tibble::as_tibble(design_matrix(term))
    if (!longnames) {
      names(dm1) <- shortnames(term)
    }
    dm1$.block <- sframe$blockids
    dm1$.time <- sframe$time
   
    
    tidyr::gather(dm1, condition, value, -.time, -.block)
  })
  
  names(dflist) <- term_names
  
  if (is.null(term_name)) {
    tn <- term_names[1]
    dfx <-  dflist[[tn]]
  } else {
    dfx <-  dflist[[term_name]]
  }
  
  ggplot2::ggplot(dfx, aes_string(x=".time", y="value", colour="condition")) + geom_line() + facet_wrap(~ .block, ncol=1) +
    xlab("Time") + theme_bw(14)
  
}
  
