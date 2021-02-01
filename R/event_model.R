
#' @export
event_model.list <- function(x, data, block, sampling_frame, drop_empty=TRUE, durations=0, precision=.3) {
  assert_that(all(sapply(x, function(a) inherits(a, "hrfspec"))), msg="event_model.list: all `x` elements must be of type `hrfspec")
  
  if (lazyeval::is_formula(block)) {
    ## TODO check for existence of block in data
    ## TODO warn when onset are way wrong
    block_rhs <- lazyeval::f_rhs(block)
    blockvals <- lazyeval::f_eval_rhs(block, data)
  } else {
    blockvals <- block
  }
  
  assert_that(is.increasing(blockvals), msg="'blockvals' must consist of strictly increasing integers")
  assert_that(length(blockvals) == nrow(data))
  
  blocks <-unique(blockvals)
  ranked_blocks <- rank(blocks)
  blockids <- ranked_blocks[match(blockvals, blocks)]
  
  if (missing(durations)) {
    ## assume zero-duration impulse for all events
    durations <- rep(0, nrow(data))
  }
  
  form <- paste("~", paste(sapply(x, function(z) z$label), collapse="+"))
  
  model_spec <- list(formula=form,
                     event_table=data, 
                     onsets=x[[1]]$onsets, 
                     event_spec=x, 
                     blockvals=blockvals,
                     blockids=blockids, 
                     durations=durations, 
                     sampling_frame=sampling_frame,
                     drop_empty=drop_empty,
                     precision=.3)
  
  class(model_spec) <- c("event_model_spec", "list")
  fmodel <- construct_model(model_spec)
  
}


#' @export
event_model.formula <- function(x, data, block, sampling_frame, drop_empty=TRUE, durations=0, precision=.3) {
  formula <- x
  stopifnot(inherits(formula, "formula"))
  assert_that(inherits(data, "data.frame"), msg="`data` must be a `data.frame`")

  if (lazyeval::is_formula(block)) {
    ## TODO check for existence of block in data
    ## TODO warn when onset are way wrong
    block_rhs <- lazyeval::f_rhs(block)
    blockvals <- lazyeval::f_eval_rhs(block, data)
  } else {
    blockvals <- block
  }
  
  assert_that(is.increasing(blockvals), msg="'blockvals' must consist of strictly increasing integers")
  assert_that(length(blockvals) == nrow(data))
  
  blocks <-unique(blockvals)
  ranked_blocks <- rank(blocks)
  blockids <- ranked_blocks[match(blockvals, blocks)]
  
  if (missing(durations)) {
    ## assume zero-duration impulse for all events
    durations <- rep(0, nrow(data))
  }
  
  formspec <- function(formula, table) {
    vterms <- extract_terms(formula, table)
    resp <- attr(vterms, "response")
    variables <- extract_variables(formula, table, vterms)
    
    if (resp != 0) {
      lhs <- variables[[resp]]
    } else {
      lhs <- NULL
    }
    
    rhs <- variables[(resp+1):length(variables)]
    
    ## vclass <- unlist(lapply(rhs, function(x) class(x)[1]))
    #ret <- list(vterms=vterms, resp=resp, variables=variables, lhs=lhs, rhs=rhs,vclass=vclass)
    ret <- list(lhs=lhs, rhs=rhs)
    class(ret) <- "formula_extraction"
    ret
  }
  
  event_spec <- formspec(formula, data)
  assertthat::assert_that(all(map_lgl(event_spec$rhs, inherits, "hrfspec")),
                          msg="all terms on right hand side must be 'hrf' terms")
  
  model_spec <- list(formula=formula, 
                     event_table=data, 
                     onsets=event_spec$lhs, 
                     event_spec=event_spec$rhs, 
                     blockvals=blockvals,
                     blockids=blockids, 
                     durations=durations, 
                     sampling_frame=sampling_frame,
                     drop_empty=drop_empty,
                     contrasts=contrasts,
                     precision=precision)
  
  class(model_spec) <- c("event_model_spec", "list")
  fmodel <- construct_model(model_spec)
  fmodel
}

#' @export
blocklens.event_model <- function(x) {
  blocklens(x$sampling_frame)
}

#' @keywords internal
construct_model <- function(x) {
  ## what does this function actually need?
  ## it neds a list of `hrfspec` objects as `rhd`
  assert_that(is.numeric(x$onsets))
  #term_names <- sapply(x$event_spec$rhs, "[[", "id")
  #browser()
  
  ## term_names <- sapply(x$event_spec$rhs, "[[", "name")
  term_names <- sapply(x$event_spec, "[[", "name")
  term_names <- .sanitizeName(term_names)
  
  dups <- sum(duplicated(term_names)) > 0
  
  if (dups) {
    dup_ids <- ave(term_names, term_names, FUN=seq_along)
    term_names <- paste0(term_names, "_", dup_ids)
  } 
  

  #terms <- lapply(x$event_spec$rhs, function(m) construct(m,x))
  terms <- lapply(x$event_spec, function(m) construct(m,x))
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
    contrasts=x$contrasts,
    model_spec=x
  )
  
  class(ret) <- c("event_model", "list")
  ret
}

#' @keywords internal
extract_terms <- function(formula, data) {
  if (!inherits(formula, "terms")) {
    terms(formula, data = data)
  } else {
    formula
  }	
}

#' @keywords internal
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

#' @keywords internal
is_parametric_basis <- function(obj) { inherits(obj, "ParametricBasis") }

#' @keywords internal
extract_variables <- function(form, data, .terms=NULL) {
  if (is.null(.terms)) {
    .terms <- extract_terms(form,data)
  }
  
 
  env <- environment(.terms)
  #env <- environment(.terms)
  varnames <- sapply(attr(.terms, "variables") , deparse, width.cutoff = 500)[-1]
  #varnames <- c(formula.tools::get.vars(formula.tools::lhs(form)), formula.tools::get.vars(formula.tools::rhs(form)))
  #variables <- as.data.frame(do.call(cbind, lapply(varnames, function(vn) eval(parse(text=vn), data,environment(terms(form))))))
  #variables <- eval(attr(.terms, "variables"), data, env) 
  variables <- eval(attr(.terms, "variables"), data) 
  names(variables) <- varnames
  variables
}

# has_block_variable <- function(form) {
#   ret <- op(rhs(form)) == "|"
#   length(ret) > 0 && ret
# }

#' @keywords internal
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


#' @export
#' @importFrom tibble as_tibble
design_matrix.event_model_spec <- function(x) {
  termlist <- lapply(x$varspec, function(m) construct(m,x))
  ret <- lapply(termlist, design_matrix)
  vnames <- unlist(lapply(ret, colnames))
  dmat <- tibble::as_tibble(do.call(cbind, ret),.name_repair="check_unique")
  names(dmat) <- vnames
  dmat
}

#' @importFrom tibble as_tibble
#' @export
#' @rdname design_matrix
design_matrix.event_model <- function(x, blockid=NULL) {
  ret <- lapply(x$terms, design_matrix, blockid)
  vnames <- unlist(lapply(ret, names))
  dmat <- tibble::as_tibble(do.call(cbind, ret),.name_repair="check_unique")
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
#' @rdname Fcontrasts
Fcontrasts.convolved_term <- function(x) {
  Fcontrasts(x$evterm)
}


# @export
# @rdname contrast_weights
#contrast_weights.fmri_term <- function(x) { stop("unimplemented") }


#' @export
#' @rdname contrast_weights
contrast_weights.fmri_model <- function(x) { 
  contrast_weights(x$event_model) 
}

#' @export
term_names.event_model <- function(x) {
  xt <- terms(x)
  unlist(lapply(xt, function(term) term$varname))
}


#' @export
#' @rdname contrast_weights
contrast_weights.event_model <- function(x) {
  tnames <- term_names(x)
  tind <- x$term_indices
  ncond <- length(conditions(x))
  ret <- lapply(seq_along(terms(x)), function(i) {
    cwlist <- contrast_weights(terms(x)[[i]])
    len <- length(conditions(terms(x)[[i]]))
    
    if (!is.null(cwlist) && length(cwlist) > 0) {
      ret <- lapply(cwlist, function(cw) {
        out <- matrix(0, ncond, ncol(cw$weights))
        out[tind[[i]],] <- cw$weights
        attr(cw, "term_indices") <- as.vector(tind[[i]])
        attr(cw, "offset_weights") <- out
        cw
      })
      
      cnames <- sapply(cwlist, function(cw) cw$name)

      #prefix <- tnames[i]
      #names(ret) <- paste0(prefix, "#", cnames)
      names(ret) <- cnames
      
      ret
    }
  })

  names(ret) <- tnames
  ret
}


#' @export
Fcontrasts.event_model <- function(x) {
  tind <- x$term_indices
  #len <- length(conditions(x))
  tnames <- names(terms(x))
  ret <- unlist(lapply(seq_along(terms(x)), function(i) {
    eterm <- terms(x)[[i]]$evterm
    len <- length(conditions(eterm))
    cwlist <- Fcontrasts(eterm)
    if (!is.null(cwlist)) {
      
      ret <- lapply(cwlist, function(cw) {
        
        out <- matrix(0, len, ncol(cw))
        #rownames(out) <- rep("C", nrow(out))
        ti <- tind[[i]]
        #out[ti,] <- cw
        out <- as.matrix(cw)
        attr(out, "term_indices") <- as.vector(ti)
        #row.names(out)[ti] <- row.names(cw)
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
#' @rdname design_matrix
design_matrix.convolved_term <- function(x, blockid=NULL) {
  if (is.null(blockid)) {
    x$design_matrix
  } else {
    keep <- blockids(x$sampling_frame) %in% blockid
    x$design_matrix[keep,]
  } 
}

design_matrix.afni_hrf_convolved_term <- function(x, blockid=NULL) {
  stop("afni_hrf_convolved_term delegates design matrix construction to AFNI")
}


#' matrix_term
#' 
#' A set of regression variables stored as a numeric matrix
#' 
#' @importFrom tibble as_tibble
#' @examples 
#' mat <- matrix(rnorm(100*10), 100 ,10)
#' mterm <- matrix_term("mterm", mat)
#' @export
#' @rdname matrix_term
matrix_term <- function(varname, mat) {
  stopifnot(is.matrix(mat))
  ret <- list(varname=varname, design_matrix=tibble::as_tibble(mat,.name_repair="check_unique"))
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
longnames.afni_hrf_convolved_term <- function(x) {
  # do not include basis term
  term.cells <- cells(x, exclude_basis=TRUE)
  # ignores exclude.basis
  apply(as.matrix(sapply(1:ncol(term.cells), 
                         function(i) {
                           paste0(names(term.cells)[i], "#", term.cells[[i]], sep="")
                         })), 1, paste, collapse=":")
}

#' @export
#' @rdname longnames
longnames.event_model <- function(x) {
  unlist(lapply(terms(x), longnames))
 
}

#' @export
#' @rdname longnames
longnames.event_term <- function(x) {
  # ignores exclude.basis
  term.cells <- cells(x)
  # ignores exclude.basis
  apply(as.matrix(sapply(1:ncol(term.cells), 
                         function(i) {
                           paste0(names(term.cells)[i], "#", term.cells[[i]], sep="")
                         })), 1, paste, collapse=":")
}


#' @export
#' @rdname shortnames
shortnames.event_model <- function(x) {
  unlist(lapply(terms(x), shortnames))
}

#' @export
#' @rdname shortnames
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
#' @rdname shortnames
shortnames.matrix_term <- function(x) {
  colnames(x$design_matrix)
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
  term_names <- sapply(all_terms, "[[", "varname")
  
  sframe <- x$sampling_frame
  
  condition = NULL; value = NULL; .time = NULL; .block=NULL
  
  dflist <- lapply(all_terms, function(term) {
    dm1 <- tibble::as_tibble(design_matrix(term),.name_repair="check_unique")
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
  
  p <- ggplot2::ggplot(dfx, ggplot2::aes_string(x=".time", y="value", colour="condition")) + 
    ggplot2::geom_line() + ggplot2::facet_wrap(~ .block, ncol=1) +
    ggplot2::xlab("Time") + ggplot2::theme_bw(14) 
  
  if (length(unique(dfx$condition)) > 25) {
    p <- p + ggplot2::guides(colour=FALSE)
  }
    
  p
}
  
