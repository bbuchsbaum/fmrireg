
#' baseline_model
#' 
#' @param baseline_model 
#' @param basis
#' @param degree
#' @param sampling_frame
#' @param nuisance_list a list of nusiance matrices, one per block
#' @importFrom lazyeval f_eval f_rhs f_lhs
#' @export
baseline_model <- function(basis="bs", degree=5, sampling_frame, nuisance_list=NULL) {
  drift_spec <- baseline(degree=degree, basis=basis, constant=FALSE)
  drift_term <- construct(drift_spec, sampling_frame)
  
  ## automatically add constant term
  block_term <- construct_block_term("constant", sampling_frame)
  
  nuisance_term <- if (!is.null(nuisance_list)) {
    total_len <- sum(sapply(nuisance_list, nrow))
    assertthat::assert_that(total_len == length(sampling_frame$blockids))
    assertthat::assert_that(length(nuisance_list) == length(sampling_frame$blocklens))
    
    colind <- vector(mode="list")
    if (is.null(colnames(nuisance_mat))) {
      colnames(nuisance_mat) <- paste0("nuisance#", 1:ncol(nuisance_mat))
    }
    
    matrix_term("nuisance", nuisance_mat)
  } 
  
  ret <- list(drift_term=drift_term, drift_spec=drift_spec, 
              block_term=block_term, nuisance_term=nuisance_term, 
              sampling_frame=sampling_frame)
  
  class(ret) <- c("baseline_model", "list")
  ret
}

#' baseline
#' 
#' A matrix of polynomial regressors for modeling low-frequency drift in fmri time series.
#' 
#' @importFrom splines bs ns
#' @param degree number of polynomial terms for each image block
#' @param basis the type of polynomial basis.
#' @param name the name of the term
#' @export
baseline <- function(degree=5, basis=c("bs", "poly", "ns"), name=paste0("baseline_", basis, "_", degree), constant=FALSE) {
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
  ret <- if (is.null(x$nuisance_term)) {
    list(x$block_term, x$drift_term)
  } else {
    list(x$block_term, x$drift_term, x$nuisance_term)
  }
  
  names(ret) <- sapply(ret, "[[", "varname")
  ret
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
  
  colind <- vector(length(ret), mode="list")
  rowind <- vector(length(ret), mode="list")
  
  for (i in seq_along(ret)) {
    rowstart <- sum(sampling_frame$blocklens[1:i]) - sampling_frame$blocklens[1] + 1
    rowend <- rowstart + sampling_frame$blocklens[i] -1
    colstart <- (nc_per_block)*(i-1) + 1
    colend <- colstart + nc_per_block - 1
    mat[rowstart:rowend, colstart:colend] <- ret[[i]]
    colind[[i]] <- colstart:colend
    rowind[[i]] <- rowstart:rowend
  }
  
  cnames <- apply(expand.grid(paste0("base_", x$basis, "#", 1:length(sampling_frame$blocklens)), paste0("block#", 1:nc_per_block)), 1, paste, collapse="_")
  colnames(mat) <- cnames
  baseline_term(x$name, mat, colind, rowind)	
}

#' baseline_term
#' @importFrom tibble as_tibble
#' @export
baseline_term <- function(varname, mat, colind, rowind) {
  stopifnot(is.matrix(mat))
  ret <- list(varname=varname, design_matrix=tibble::as_tibble(mat), colind=colind, rowind=rowind)
  class(ret) <- c("baseline_term", "matrix_term", "fmri_term", "list")
  ret
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

#' @importFrom tibble as_tibble
#' @export
block_term <- function(varname, blockids, expanded_blockids, mat) {
  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(nrow(mat) == length(blockids))
  assertthat::assert_that(ncol(mat) == length(unique(blockids)))
  
  ret <- list(varname=varname, blockids=blockids, expanded_blockids, design_matrix=tibble::as_tibble(mat), nblocks=ncol(mat))
  class(ret) <- c("block_term", "matrix_term", "fmri_term", "list")
  ret
}


design_matrix.block_term <- function(x, block=NULL) {
  if (is.null(block)) {
    x$design_matrix
  } else {
    x$design_matrix[, block]
  }
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



