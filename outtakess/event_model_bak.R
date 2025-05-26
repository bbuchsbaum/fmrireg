
#' event_model.list
#'
#' Construct an event model from a list of HRF specifications.
#'
#' @description
#' This function constructs an event model from a list of hemodynamic response function
#' (HRF) specifications `x`, along with other parameters like data, block, sampling frame,
#' and durations. It is specifically designed to handle lists of HRF specifications.
#'
#' @param x A list of HRF specifications. All elements in `x` must be of type `hrfspec`.
#' @param data The dataset containing the event information.
#' @param block The block variable, either a formula or a vector of block values.
#' @param sampling_frame The time series grid over which to sample the function.
#' @param drop_empty Logical indicating whether to drop empty events. Default is `TRUE`.
#' @param durations A numeric vector of event durations. Default is 0 for all events.
#' @param precision Numeric value indicating the precision hrf sampling. Default is 0.3.
#' @param ... Additional arguments.
#' @return An event model object.
#' @export
event_model.list <- function(x, data, block, sampling_frame, drop_empty=TRUE, durations=0, precision=.3, ...) {
  assert_that(all(sapply(x, function(a) inherits(a, "hrfspec"))), msg="event_model.list: all `x` elements must be of type `hrfspec")
  
  if (lazyeval::is_formula(block)) {
    ## TODO check for existence of block in data
    ## TODO warn when onset are way wrong
    block_rhs <- lazyeval::f_rhs(block)
    blockvals <- lazyeval::f_eval_rhs(block, data)
  } else {
    blockvals <- block
  }
  if (is.factor(blockids)) {
    blockvals <- as.integer(as.character(blockvals))
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


#' event_model.formula
#'
#' Construct an event model from a formula and data.
#'
#' @description
#' This function constructs an event model using a formula `x`, a dataset `data`,
#' a block variable, and other parameters like sampling frame, drop_empty, durations,
#' and precision. It is specifically designed to handle models specified as formulas.
#'
#' @param x A formula specifying the event model.
#' @param data The dataset containing the event information.
#' @param block The block variable, either a formula or a vector of block values.
#' @param sampling_frame The time series grid over which to sample the function.
#' @param drop_empty Logical indicating whether to drop empty events. Default is `TRUE`.
#' @param durations A numeric vector of event durations. Default is 0 for all events.
#' @param precision Numeric value indicating the precision of the model. Default is 0.3.
#' @param ... Additional arguments.
#' @return An event model object.
#' @export
event_model.formula <- function(x, data, block, sampling_frame, drop_empty=TRUE, durations=0, precision=.3, ...) {
  formula <- x
  stopifnot(inherits(formula, "formula"))
  assert_that(inherits(data, "data.frame"), msg="`data` must be a `data.frame`")

  if (lazyeval::is_formula(block)) {
    ## TODO check for existence of block in data
    ## TODO warn when onset are way off
    block_rhs <- lazyeval::f_rhs(block)
    blockvals <- lazyeval::f_eval_rhs(block, data)
  } else {
    blockvals <- block
  }
  
  if (is.factor(blockids)) {
    blockvals <- as.integer(as.character(blockvals))
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
blocklens.event_model <- function(x,...) {
  blocklens(x$sampling_frame)
}

#' @export
blockids.event_model <- function(x) {
  x$blockids
}


#' @keywords internal
#' @noRd
construct_model <- function(x) {
  ## what does this function actually need?
  ## it neds a list of `hrfspec` objects as `rhd`
  assert_that(is.numeric(x$onsets))
  ##term_names <- sapply(x$event_spec$rhs, "[[", "id")
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
#' @noRd
extract_terms <- function(formula, data) {
  if (!inherits(formula, "terms")) {
    terms(formula, data = data)
  } else {
    formula
  }	
}

# @keywords internal
# extract_covariates <- function(.terms, variables, resp, etab) {
#   vars <- attr(.terms, "variables") 
#   varnames <- sapply(vars, deparse, width.cutoff = 500)[-1]
#   ind.vars <- varnames[-resp] 
#   
#   
#   orig.covar.names <- ind.vars[which(sapply(variables[-resp], function(obj) is.numeric(obj) || is_parametric_basis(obj)))]
#   new.covar.names <- names(events(etab))[sapply(events(etab), is_continuous)]
#   covar.names <- as.list(orig.covar.names)
#   names(covar.names) <- new.covar.names
#   covar.names
# }

#' @keywords internal
#' @noRd
is_parametric_basis <- function(obj) { inherits(obj, "ParametricBasis") }

#' @keywords internal
#' @noRd
extract_variables <- function(form, data, .terms=NULL) {
  if (is.null(.terms)) {
    .terms <- extract_terms(form,data)
  }
  
 
  env <- environment(.terms)
  #env <- environment(.terms)
  varnames <- sapply(attr(.terms, "variables") , deparse, width.cutoff = 500)[-1]
  print(varnames)
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
#' @noRd
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
design_matrix.event_model_spec <- function(x, ...) {
  termlist <- lapply(x$varspec, function(m) construct(m,x))
  ret <- lapply(termlist, design_matrix)
  vnames <- unlist(lapply(ret, colnames))
  dmat <- suppressMessages(tibble::as_tibble(do.call(cbind, ret),.name_repair="check_unique"))
  names(dmat) <- vnames
  dmat
}

#' @importFrom tibble as_tibble
#' @export
#' @rdname design_matrix
design_matrix.event_model <- function(x, blockid=NULL, ...) {
  ret <- lapply(x$terms, design_matrix, blockid)
  vnames <- unlist(lapply(ret, names))
  dmat <- suppressMessages(tibble::as_tibble(do.call(cbind, ret),.name_repair="check_unique"))
  names(dmat) <- vnames
  dmat
}



#' @export
terms.event_model <- function(x,...) {
  x$terms
}

#' @export
conditions.event_model <- function(x,...) {
  unlist(lapply(terms(x), function(t) conditions(t)), use.names=FALSE)
}

#' @export
#' @rdname contrast_weights
contrast_weights.convolved_term <- function(x,...) {
  lapply(x$contrasts, function(cspec) {
    if (!is.null(cspec))
      contrast_weights(cspec,x)
  })
}

#' @export
#' @rdname Fcontrasts
Fcontrasts.convolved_term <- function(x,...) {
  #browser()
  #cellcount <- attr(cells(x, drop.empty=FALSE), "count")
  Fcontrasts(x$evterm)
}


# @export
# @rdname contrast_weights
#contrast_weights.fmri_term <- function(x) { stop("unimplemented") }


#' @export
#' @rdname contrast_weights
contrast_weights.fmri_model <- function(x,...) { 
  contrast_weights(x$event_model) 
}

#' @export
term_names.event_model <- function(x) {
  xt <- terms(x)
  unlist(lapply(xt, function(term) term$varname))
}


#' @export
#' @rdname contrast_weights
contrast_weights.event_model <- function(x,...) {
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
Fcontrasts.event_model <- function(x,...) {
  tind <- x$term_indices
  #len <- length(conditions(x))
  tnames <- names(terms(x))
  ret <- unlist(lapply(seq_along(terms(x)), function(i) {
    eterm <- terms(x)[[i]]$evterm
    len <- length(conditions(eterm))
    cwlist <- Fcontrasts(terms(x)[[i]])
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
design_matrix.convolved_term <- function(x, blockid=NULL, ...) {
  if (is.null(blockid)) {
    x$design_matrix
  } else {
    keep <- blockids(x$sampling_frame) %in% blockid
    x$design_matrix[keep,]
  } 
}

#' @export
design_matrix.afni_hrf_convolved_term <- function(x, blockid=NULL, ...) {
  stop("afni_hrf_convolved_term delegates design matrix construction to AFNI")
}


#' matrix_term
#'
#' Create a matrix_term object.
#'
#' @description
#' This function creates a matrix_term object, which is a set of regression
#' variables stored as a numeric matrix.
#'
#' @param varname The name of the variable.
#' @param mat The matrix of values.
#' @return A matrix_term object.
#' @importFrom tibble as_tibble
#' @examples 
#' mat <- matrix(rnorm(100 * 10), 100, 10)
#' mterm <- matrix_term("mterm", mat)
#' @export
matrix_term <- function(varname, mat) {
  stopifnot(is.matrix(mat))
  ret <- list(varname=varname, design_matrix=suppressMessages(tibble::as_tibble(mat,.name_repair="minimal")))
  class(ret) <- c("matrix_term", "fmri_term", "list")
  ret
}


#' @export
#' @rdname design_matrix
design_matrix.matrix_term <- function(x, ...) {
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
nbasis.convolved_term <- function(x,...) nbasis(x$hrf)

#' @export
cells.event_model <- function(x, ...) {
  ret <- lapply(x$terms, function(z) cells(z, ...))
  tnames = names(ret)
  out <- do.call(rbind, lapply(1:length(tnames), function(i) {
    dplyr::tibble(term=tnames[i], level=ret[[i]][,1], basis=ret[[i]][,2])
  })) %>% dplyr::mutate(index=1:dplyr::n())
}

#' @export
#' @rdname longnames
longnames.convolved_term <- function(x, ...) {
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
longnames.afni_hrf_convolved_term <- function(x, ...) {
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
longnames.event_model <- function(x, ...) {
  unlist(lapply(terms(x), longnames))
 
}

#' @export
#' @rdname longnames
longnames.event_term <- function(x, ...) {
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
shortnames.event_model <- function(x, ...) {
  unlist(lapply(terms(x), shortnames))
}

#' @export
#' @rdname shortnames
shortnames.convolved_term <- function(x, ...) {
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
shortnames.matrix_term <- function(x, ...) {
  colnames(x$design_matrix)
}


#' @export
#' @rdname longnames
longnames.matrix_term <- function(x, ...) {
  paste0(x$name, "#", colnames(design_matrix(x)))
}


#' @export
print.event_model <- function(x, ...) {
  # Header with fancy border
  cat("\n╔══════════════════════════════════════════╗")
  cat("\n║           fMRI Event Model               ║")
  cat("\n╚══════════════════════════════════════════╝\n")
  
  # Model Formula Section
  cat("\n Model Formula:\n")
  cat("  ", crayon::cyan(Reduce(paste, deparse(x$model_spec$formula))), "\n")
  
  # Summary Statistics
  cat("\n📈 Model Summary:\n")
  cat("  • Number of Terms:", crayon::yellow(length(terms(x))), "\n")
  cat("  • Total Events:", crayon::yellow(nrow(x$model_spec$event_table)), "\n")
  cat("  • Design Matrix Columns:", crayon::yellow(length(conditions(x))), "\n")
  cat("  • Number of Blocks:", crayon::yellow(length(unique(x$blockids))), "\n")
  
  # Event Table Preview
  if (nrow(x$model_spec$event_table) > 0) {
    cat("\n📋 Event Table Preview:\n")
    cat("  • Variables:", crayon::green(paste(names(x$model_spec$event_table), collapse=", ")), "\n")
    if (nrow(x$model_spec$event_table) > 3) {
      print(head(x$model_spec$event_table, 3))
      cat("  ... (", nrow(x$model_spec$event_table) - 3, " more rows )\n")
    } else {
      print(x$model_spec$event_table)
    }
  }
  
  # Individual Terms Section
  cat("\n🔍 Model Terms:\n")
  for (i in seq_along(terms(x))) {
    term <- terms(x)[[i]]
    cat("\n  Term", crayon::blue(i), ":", crayon::magenta(names(terms(x))[i]), "\n")
    
    # Get conditions for this term
    term_conditions <- conditions(term)
    if (length(term_conditions) > 0) {
      cat("    • Conditions:", crayon::green(paste(term_conditions, collapse=", ")), "\n")
    }
    
    # Get HRF information if available
    if (!is.null(term$hrfspec)) {
      cat("    • HRF Type:", crayon::yellow(attr(term$hrfspec$hrf, "name")), "\n")
      if (!is.null(attr(term$hrfspec$hrf, "nbasis"))) {
        cat("    • Basis Functions:", crayon::yellow(attr(term$hrfspec$hrf, "nbasis")), "\n")
      }
    }
    
    # Get contrast information if available
    if (!is.null(term$contrasts) && length(term$contrasts) > 0) {
      cat("    • Contrasts:", crayon::cyan(paste(names(term$contrasts), collapse=", ")), "\n")
    }
  }
  
  cat("\n") # End with newline for cleaner output
}


#' Plot an event_model object
#'
#' @description
#' This function creates a plot of an event_model object, visualizing the design matrix.
#'
#' @param x An event_model object.
#' @param y Unused.
#' @param term_name Name of the term to plot (optional).
#' @param longnames Logical flag; if TRUE, use long names for the plot (default: TRUE).
#' @param ... Additional arguments to be passed to the ggplot2 function.
#' @importFrom ggplot2 ggplot aes_string geom_line facet_wrap xlab theme_bw
#' @importFrom tidyr gather
#' @return A ggplot2 plot object.
#' @export
plot.event_model <- function(x, y, term_name=NULL, longnames=TRUE, ...) {
  all_terms <- terms(x)
  term_names <- sapply(all_terms, "[[", "varname")
  
  sframe <- x$sampling_frame
  
  condition = NULL; value = NULL; .time = NULL; .block=NULL
  
  dflist <- lapply(all_terms, function(term) {
    dm1 <- suppressMessages(tibble::as_tibble(design_matrix(term),.name_repair="check_unique"))
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
  

#' Create an event model directly from components, accepting variables as strings or expressions
#'
#' @description
#' This function constructs an event model directly from its components, allowing variables
#' to be specified as character strings or expressions, similar to the formula interface.
#'
#' @param event_terms A list of event term specifications. Each term is a list with components:
#'        - `variables`: Character vector of variable names or expressions (captured using `rlang::exprs`)
#'        - `hrf`: HRF specification for the term
#' @param events A data frame of event variables. Variables involved in terms must be present here.
#' @param onsets A numeric vector of event onset times in seconds
#' @param block A vector of block IDs (must be strictly increasing integers)
#' @param sampling_frame The time series grid over which to sample the function
#' @param durations A numeric vector of event durations. Default is 0 for all events.
#' @param drop_empty Logical indicating whether to drop empty events. Default is TRUE.
#' @param precision Numeric value indicating the precision of HRF sampling. Default is 0.3.
#'
#' @return An event_model object
#'
#' @examples
#' library(fmrireg)
#' library(rlang)
#'
#' # Example with variables specified as character strings
#' event_terms <- list(
#'   list(
#'     variables = c("x", "y"),
#'     hrf = "spmg1"
#'   )
#' )
#'
#' # Example with variables specified as expressions
#' event_terms <- list(
#'   list(
#'     variables = exprs(x, Poly(y, 2)),
#'     hrf = "spmg1"
#'   )
#' )
#'
#' @export
create_event_model <- function(event_terms, events, onsets, block, sampling_frame,
                               durations=0, drop_empty=TRUE, precision=0.3) {
  
  # Input validation
  assert_that(is.list(event_terms), msg="event_terms must be a list of term specifications")
  assert_that(length(onsets) == nrow(events), msg="onsets length must match number of events")
  
  # Process block IDs
  if (is.factor(block)) {
    block <- as.integer(as.character(block))
  }
  assert_that(is.increasing(block), msg="'block' must consist of strictly increasing integers")
  assert_that(length(block) == nrow(events))
  
  blocks <- unique(block)
  ranked_blocks <- rank(blocks)
  blockids <- ranked_blocks[match(block, blocks)]
  
  # Process durations
  if (length(durations) == 1) {
    durations <- rep(durations, length(onsets))
  }
  
  # Capture event data in a data environment
  data_env <- list2env(as.list(events), parent = as.environment("package:fmrireg"))
  
  # Create hrfspecs
  hrfspecs <- lapply(event_terms, function(term) {
    variables <- term$variables
    hrf_spec <- term$hrf
    
    # Prepare variables, accepting both strings and expressions
    if (is.character(variables)) {
      # Convert character strings to symbols
      variables <- lapply(variables, rlang::sym)
    } else if(is.expression(variables) || is_quosure(variables)) {
      # Single expression, wrap in list
      variables <- list(variables)
    } else if (!is.list(variables)) {
      stop("variables in event_terms must be character vector or expressions captured by rlang::exprs")
    }
    
    # Prepare HRF basis
    if (is.character(hrf_spec)) {
      # Built-in HRF
      hrf <- getHRF(hrf_spec)
    } else if (inherits(hrf_spec, "HRF")) {
      # Direct HRF object
      hrf <- hrf_spec
    } else if (is.list(hrf_spec)) {
      # HRF with parameters
      hrf <- do.call(getHRF, c(list(hrf_spec$hrf), hrf_spec$parameters))
    } else {
      stop("Invalid HRF specification")
    }
    
    # Create hrfspec object
    hrfspec(
      vars = variables,
      basis = hrf,
      onsets = onsets,
      durations = durations,
      precision = precision,
      data_env = data_env  # Pass the data environment
    )
  })
  
  # Create model specification
  # Construct formula from components
  formula_terms <- lapply(event_terms, function(term) {
    variables <- term$variables
    hrf_spec <- term$hrf
    
    # Prepare variables for formula term string
    var_strings <- sapply(variables, function(var) {
      if (is.symbol(var) || is.language(var)) {
        deparse(var)
      } else if (is.character(var)) {
        var
      } else {
        stop("variables must be symbols, expressions, or character strings")
      }
    })
    
    # Get HRF type for the term
    hrf_type <- if (is.character(hrf_spec)) {
      hrf_spec
    } else if (inherits(hrf_spec, "HRF")) {
      attr(hrf_spec, "name")
    } else if (is.list(hrf_spec)) {
      hrf_spec$hrf
    }
    
    # Construct term string
    var_string <- paste(var_strings, collapse=", ")
    paste0("hrf(", var_string, ", basis='", hrf_type, "')")
  })
  
  # Combine terms with + and create formula
  formula_str <- paste("onsets ~", paste(formula_terms, collapse=" + "))
  model_formula <- as.formula(formula_str)
  
  model_spec <- list(
    formula = model_formula,
    event_table = events,
    onsets = onsets,
    event_spec = hrfspecs,
    blockvals = block,
    blockids = blockids,
    durations = durations,
    sampling_frame = sampling_frame,
    drop_empty = drop_empty,
    precision = precision
  )
  
  class(model_spec) <- c("event_model_spec", "list")
  
  # Construct the model
  construct_model(model_spec)
}
  
  