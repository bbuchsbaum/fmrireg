
#' construct a covariate term that is not convolved with a hemodynamic response function
#' 
#' @param ... a variable argument set of covariate names
#' @param data a data.frame containing the variables
#' @export
#' @inheritParams hrf
#' 
#' @examples
#' df1 <- data.frame(x=rnorm(100), y=rnorm(100))
#' cv <- covariate(x,y,data=df1)
covariate <- function(..., data, id=NULL, prefix=NULL) {
  vars <- as.list(substitute(list(...)))[-1] 
  parsed <- parse_term(vars, "covariate")
  term <- parsed$term
  label <- parsed$label
  
  varnames <- if (!is.null(prefix)) {
    paste0(prefix, "_", term)
  } else {
    term
  }
  
  termname <- paste0(varnames, collapse="::")
  
  if (is.null(id)) {
    id <- termname
  }  

  ret <- list(
    data=data,
    name=termname, ## covariate(x,y), where termname = "x::y"
    id=id, ## covariate(x), id by default is "x::y"
    varnames=varnames, ## list of all variables (e.g. list(x,y))
    vars=term, ## list of unparsed vars
    label=label, ## "covariate(x)" the full expression
    subset=substitute(subset))
  
  class(ret) <- c("covariatespec", "hrfspec", "list")
  ret
}

#' @keywords internal
covariate_term <- function(varname, mat) {
  stopifnot(is.matrix(mat))
  ret <- list(varname=varname, design_matrix=tibble::as_tibble(mat))
  class(ret) <- c("covariate_term", "matrix_term", "fmri_term", "list")
  ret
}

#' @export
construct.covariatespec <- function(x, model_spec, sampling_frame=NULL) {
  mat <- do.call(cbind, lapply(x$vars, function(v) {
    eval(parse(text=v), envir=x$data)
  }))
  
  colnames(mat) <- x$varnames
  
  cterm <- covariate_term(x$name, mat)
  
  sframe <- if (is.null(sampling_frame)) {
    model_spec$sampling_frame
  } else {
    sampling_frame
  }
  
  ret <- list(
    varname=x$name,
    spec=x,
    evterm=cterm,
    design_matrix=cterm$design_matrix,
    sampling_frame=sframe,
    id=if(!is.null(x$id)) x$id else x$varname
  )
  
  class(ret) <- c("covariate_convolved_term", "convolved_term", "fmri_term", "list") 
  ret
}

#' @export
event_table.covariate_convolved_term <- function(x) {
  cnames <- colnames(x$design_matrix)
  ret <- do.call(cbind, lapply(cnames, function(tname) {
    rep(.sanitizeName(tname), nrow(x$design_matrix))
  }))
  
  colnames(ret) <- cnames
  as_tibble(ret,.name_repair="check_unique")
  
}

#' @export
nbasis.covariate_convolved_term <- function(x) {
  ncol(x$design_matrix)
}
