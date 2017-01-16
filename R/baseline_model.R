

#' @param baseline_model 
#' @param basis
#' @param sampling_frame
#' @importFrom lazyeval f_eval f_rhs f_lhs
#' @export
baseline_model <- function(basis="bs", degree=5, sampling_frame, nuisance_matrix=NULL) {
  drift_spec <- baseline(degree=degree, basis=basis, constant=FALSE)
  drift_term <- construct(drift_spec, sampling_frame)
  
  ## automatically add constant term
  block_term <- construct_block_term("constant", sampling_frame)
  
  nuisance_term <- if (!is.null(nuisance_matrix)) {
    assertthat::assert_that(nrow(nuisance_matrix) == length(sampling_frame$blockids))
    if (is.null(colnames(nuisance_matrix))) {
      colnames(nuisance_matrix) <- paste0("nuisance#", 1:ncol(nuisance_matrix))
    }
    matrix_term("nuisance", nuisance_matrix)
  } 
  
  ret <- list(drift_term=drift_term, drift_spec=drift_spec, block_term=block_term, nuisance_term=nuisance_term)
  class(ret) <- c("baseline_model", "list")
  ret
}

#' baseline
#' 
#' A matrix of polynomial regressors for modeling low-frequency drift in fmri time series.
#' 
#' @importFrom splines bs ns
#' @param number of polynomial terms for each image block
#' @param basis the type of polynomial basis.
#' @param name the name of the term
#' @export
baseline <- function(degree=5, basis=c("bs", "poly", "ns"), name=paste0("Baseline_", basis, "_", degree), constant=FALSE) {
  basis <- match.arg(basis)
  bfun <- switch(basis,
                 bs=bs,
                 ns=ns,
                 poly=poly)
  
  ret <- list(
    degree=degree,
    basis=basis,
    fun=bfun,
    constant=constant,
    name=name
  )
  
  class(ret) <- c("baselinespec", "nuisancespec")
  ret
}


#' @export
design_matrix.baseline_model <- function(x) {
  if (is.null(x$nuisance_term)) {
    tibble::as_tibble(cbind(design_matrix(x$block_term), design_matrix(x$drift_term)))
  } else {
    tibble::as_tibble(cbind(design_matrix(x$block_term), design_matrix(x$drift_term), design_matrix(x$nuisance_term)))
  }
}


#' @export
terms.baseline_model <- function(x) {
  if (is.null(x$nuisance_term)) {
    list(x$block_term, x$drift_term)
  } else {
    list(x$block_term, x$drift_term, x$nuisance_term)
  }
}


#' a block variable, which is constant over the span of a scanning run
#' 
#' @value an instance of a class \code{blockspec}
#' @export
block <- function(x) {
  varname <- substitute(x)
  pterm <- parse_term(as.list(substitute(x)), "block")
  
  ret <- list(
    name=varname,
    label=pterm$label
  )
  
  class(ret) <- "blockspec"
  ret
  
}



#' @export  
construct.baselinespec <- function(x, sampling_frame) {
  ret <- if (x$constant) {
    lapply(sampling_frame$blocklens, function(bl) cbind(rep(1, bl), x$fun(seq(1, bl), x$degree)))
  } else {
    lapply(sampling_frame$blocklens, function(bl) x$fun(seq(1, bl), x$degree))
  }
  
  
  nc <- sum(sapply(ret, ncol))
  nc_per_block <- ncol(ret[[1]])
  mat <- matrix(0, sum(sampling_frame$blocklens), nc)
  
  
  for (i in seq_along(ret)) {
    rowstart <- sum(sampling_frame$blocklens[1:i]) - sampling_frame$blocklens[1] + 1
    rowend <- rowstart + sampling_frame$blocklens[i] -1
    colstart <- (nc_per_block)*(i-1) + 1
    colend <- colstart + nc_per_block - 1
    mat[rowstart:rowend, colstart:colend] <- ret[[i]]
  }
  
  
  cnames <- apply(expand.grid(paste0("base_", x$basis, "#", 1:length(sampling_frame$blocklens)), paste0("block#", 1:nc_per_block)), 1, paste, collapse="_")
  colnames(mat) <- cnames
  matrix_term(x$name, mat)	
}

construct_block_term <- function(vname, sampling_frame) {
  blockids <- sampling_frame$blockids
  blockord <- sort(unique(blockids))
  
  expanded_blockids <- ordered(rep(blockord, sampling_frame$blocklens))
  
  mat <- if (length(levels(blockids)) == 1) {
    matrix(rep(1, length(expanded_blockids)), length(expanded_blockids), 1)
  } else {
    mat <- model.matrix(~ expanded_blockids - 1)
  }
  
  colnames(mat) <- paste0(vname, "#", blockord)
  block_term(vname, blockids=blockids, expanded_blockids=expanded_blockids, mat)
  
}

#' nuisance
#' 
#' a 'nuisance' term that consists of an arbitrary numeric matrix with the same number of rows as image time points.
#' 
#' @export
#' @param x a \code{matrix} 
#' @return a class of type \code{nuisancespec}
nuisance <- function(x) {
  varname <- substitute(x)
  
  ret <- list(
    name=varname
  )
  
  class(ret) <- "nuisancespec"
  ret
}

#' @export
construct.nuisancespec <- function(x, model_spec) {
  mat <- base::eval(parse(text=x$varname), envir=model_spec$aux_data, enclos=parent.frame())
  
  matrix_term(x$name, mat)
}  



#' @export
construct.blockspec <- function(x, model_spec) {
  construct_block_term(x$name, model_spec$sampling_frame)
}

print.baseline_model <- function(x) {
  cat("baseline_model", "\n")
  cat("  ", "name: ", x$drift$varname, "\n")
  cat("  ", "basis type: ", x$drift_spec$basis, "\n")
  cat("  ", "degree: ", x$drift_spec$degree, "\n")
  cat("  ", "drift columns: ", ncol(design_matrix(x$drift_term)), "\n")
  cat("  ", "constant columns: ", ncol(design_matrix(x$block_term)), "\n")
  cat("  ", "nuisance columns: ", ncol(design_matrix(x$nuisance_term)), "\n")
  cat("  ", "total columns: ", ncol(design_matrix(x)), "\n")
  cat("  ", "design_matrix: ", "\n")
  print(design_matrix(x))
   
}



