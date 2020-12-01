

#' @importFrom purrr map_int
get_col_inds <- function(Xlist) {
  ncols <- purrr::map_int(Xlist, ncol)
  csum <- cumsum(ncols)
  csum1 <- c(0, csum[-length(csum)])
  m <- as.matrix(cbind(csum1+1, csum))
  
  lapply(1:nrow(m), function(i) seq(m[i,1], m[i,2]))
  #apply(m, 1, function(x) seq(x[1], x[2]))
}




#' construct a baseline model to model noise and other non-event-related sources of variance
#' 
#' @param basis the basis function type
#' @param degree the degree of the spline function.
#' @param sframe sframe a \code{sampling_frame} object
#' @param intercept whether to include an intercept for each block. Automatically set to \code{FALSE} when basis == "constant". 
#' @param nuisance_list a list of nuisance matrices, one matrix per fMRI block
#' @importFrom lazyeval f_eval f_rhs f_lhs
#' 
#' @examples 
#' 
#' ## bspline basis with degree = 3. This will produce a design matrix with three splines regressors and a constant intercept.
#' sframe <- sampling_frame(blocklens=c(100,100), TR=2)
#' bmod <- baseline_model(basis="bs", degree=3, sframe=sframe)
#' bmod_global <- baseline_model(basis="bs", degree=3, sframe=sframe, intercept="global")
#' bmod_nointercept <- baseline_model(basis="bs", degree=3, sframe=sframe, intercept="none")
#' stopifnot(ncol(design_matrix(bmod)) == 8)
#' stopifnot(ncol(design_matrix(bmod_global)) == 7)
#' stopifnot(ncol(design_matrix(bmod_nointercept)) == 6)
#' 
#' ## polynomial with no intercept term
#' bmod <- baseline_model(basis="poly", degree=3, sframe=sframe, intercept="none")
#' 
#' ## a baseline model that only has dummy-coded intercept terms, one per block, i.e. to model runwise mean shifts only.
#' bmod <- baseline_model(basis="constant", degree=1, sframe=sframe)
#' 
#' ## global intercept only
#' bmod <- baseline_model(basis="constant", degree=1, sframe=sframe, intercept="global")
#' 
#' 
#' ## add an arbitrary nuisance matrix with two columns, i.e. motion regressors, physiological noise, etc.
#' nuismat <- matrix(rnorm(100*2), 100, 2)
#' bmod <- baseline_model(basis="bs", degree=3, sframe=sframe, nuisance_list=list(nuismat, nuismat))
#' stopifnot(ncol(design_matrix(bmod)) == 12)
#' @export
baseline_model <- function(basis=c("constant", "poly", "bs", "ns"), degree=1, sframe, 
                           intercept=c("runwise", "global", "none"), nuisance_list=NULL) {
  
  basis <- match.arg(basis)
  intercept <- match.arg(intercept)
  
  if (basis %in% c("bs", "ns")) {
    assert_that(degree > 2, msg ="'bs' and 'ns' bases must have degree >= 3")
  }
  
  drift_spec <- baseline(degree=degree, basis=basis, intercept=intercept)
  drift_term <- construct(drift_spec, sframe)
  
  block_term <- if (intercept != "none" && basis != "constant") {
    construct_block_term("constant", sframe, intercept)
  }
  
  nuisance_term <- if (!is.null(nuisance_list)) {
   
    total_len <- sum(purrr::map_int(nuisance_list, nrow))
    assertthat::assert_that(total_len == length(blockids(sframe)), msg="number of rows of `nuisance_list' must match number of scans in sampling_frame")
    assertthat::assert_that(length(nuisance_list) == length(blocklens(sframe)), msg="number of `nuisance_list` elements must match number of blocks in sampling_frame")
    
    colind <- get_col_inds(nuisance_list)
    rowind <- split(1:length(blockids(sframe)), blockids(sframe))
   
    for (i in 1:length(nuisance_list)) {
      nmat <- nuisance_list[[i]]
      colnames(nmat) <- paste0("nuisance#", i, "#", 1:ncol(nmat))
      ## 'unclass' is used below because Matrix::bdiag for some reason doesn't work with class of type c("poly", "matrix")
      nuisance_list[[i]] <- unclass(as.matrix(nmat))
    }
    
    
    #baseline_term("nuisance", Matrix::bdiag(lapply(nuisance_list, unclass)), colind,rowind)
    cnames <- unlist(lapply(nuisance_list, colnames))
    nmat <- as.matrix(Matrix::bdiag(nuisance_list))
    colnames(nmat) <- cnames
    baseline_term("nuisance", nmat, colind,rowind)
  } 
  
  ret <- list(drift_term=drift_term, 
              drift_spec=drift_spec, 
              block_term=block_term, 
              nuisance_term=nuisance_term, 
              sampling_frame=sframe)
  
  class(ret) <- c("baseline_model", "list")
  ret
}


#' Create a model specification for modeling low-frequency drift in fmri time series.
#' 
#' @importFrom splines bs ns
#' 
#' @param degree number of basis terms for each image block (ignored for 'constant' basis)
#' @param basis the type of polynomial basis.
#' @param name the name of the term
#' @param constant whether to include an intercept term
#' @param intercept the type of intercept to include
#' @export
baseline <- function(degree=1, basis=c("constant", "poly", "bs", "ns"), name=NULL, 
                     intercept=c("runwise", "global", "none")) {
  
  basis <- match.arg(basis)
  intercept <- match.arg(intercept)
  
  if (basis == "constant") {
    degree <- 1
  }
  
  bfun <- switch(basis,
                 bs=splines::bs,
                 ns=splines::ns,
                 poly=poly,
                 constant=function(x, degree) { matrix(rep(1, length(x))) })
  
  if (is.null(name)) {
    name <- paste0("baseline_", basis, "_", degree)
  }
  
  ret <- list(
    degree=degree,
    basis=basis,
    fun=bfun,
    intercept=intercept,
    name=name
  )
  
  class(ret) <- c("baselinespec", "nuisancespec")
  ret
}


#' @export
design_matrix.baseline_model <- function(x, blockid=NULL, allrows=FALSE) {
  if (is.null(x$nuisance_term) && is.null(x$block_term)) {
    ## just drift term
    tibble::as_tibble(design_matrix(x$drift_term, blockid,allrows),.name_repair="check_unique")
  } else if (is.null(x$nuisance_term) && !is.null(x$block_term)) {
    ## drift and block term but no nuisance term
    tibble::as_tibble(cbind(design_matrix(x$block_term, blockid, allrows), 
                            design_matrix(x$drift_term, blockid, allrows)), .name_repair="check_unique")
  } else if (!is.null(x$nuisance_term) && is.null(x$block_term)) {
    ## drift and  nuisance term
    tibble::as_tibble(cbind(design_matrix(x$drift_term, blockid, allrows),
                            design_matrix(x$nuisance_term, blockid, allrows)), .name_repair="check_unique") 
  } else {
    ## all three 
    tibble::as_tibble(cbind(design_matrix(x$block_term, blockid, allrows), 
                            design_matrix(x$drift_term, blockid, allrows), 
                            design_matrix(x$nuisance_term, blockid, allrows)),.name_repair="check_unique")
  }
}



#' @export
terms.baseline_model <- function(x) {
  ret <- list(x$block_term, x$drift_term, x$nuisance_term)
  ret <- ret[!map_lgl(ret, is.null)]
  
  names(ret) <- unlist(lapply(ret, "[[", "varname"))
  ret
}


#' a block variable, which is constant over the span of a scanning run
#' 
#' @param x the block variable
#' @return an instance of a class \code{blockspec}
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
#' @importFrom purrr map_chr
construct.baselinespec <- function(x, sampling_frame) {
  ## ret <- if (x$constant && x$basis != "constant") {
  ## add intercept term per run
  ##lapply(sampling_frame$blocklens, function(bl) cbind(rep(1, bl), x$fun(seq(1, bl), degree=x$degree)))
  ## } else {
  ## lapply(sampling_frame$blocklens, function(bl) x$fun(seq(1, bl), degree=x$degree))
  ##}
  
  ret <- lapply(sampling_frame$blocklens, function(bl) x$fun(seq(1, bl), degree=x$degree))
  
  if (x$basis == "constant" && x$intercept == "global") {
    colind <- lapply(ret, function(x) 1:1)
    rowind <- lapply(1:length(ret), function(i) {
      rowstart <- sum(sampling_frame$blocklens[1:i]) - sampling_frame$blocklens[i] + 1
      rowend <- rowstart + sampling_frame$blocklens[i] -1
      rowstart:rowend
    })
    ret <- do.call(rbind, ret)
    cnames <- paste0("base_", x$basis)
    colnames(ret) <- cnames
    baseline_term(x$name, ret,colind, rowind)
  } else {
    nc <- sum(map_int(ret, ncol))
    nc_per_block <- ncol(ret[[1]])
    mat <- matrix(0, sum(sampling_frame$blocklens), nc)
  
    colind <- vector(length(ret), mode="list")
    rowind <- vector(length(ret), mode="list")
  
    for (i in seq_along(ret)) {
      rowstart <- sum(sampling_frame$blocklens[1:i]) - sampling_frame$blocklens[i] + 1
      rowend <- rowstart + sampling_frame$blocklens[i] -1
      colstart <- (nc_per_block)*(i-1) + 1
      colend <- colstart + nc_per_block - 1
      mat[rowstart:rowend, colstart:colend] <- ret[[i]]
      colind[[i]] <- colstart:colend
      rowind[[i]] <- rowstart:rowend
    }

    cnames <- apply(expand.grid(paste0("base_", x$basis, "", 1:nc_per_block), paste0("block_", 1:length(sampling_frame$blocklens))), 1, paste, collapse="_")
    colnames(mat) <- cnames
    baseline_term(x$name, mat, colind, rowind)	
  }
}

#' baseline_term
#' 
#' @param varname the name of the term
#' @param mat the matrix of covariates
#' @param colind the column indices
#' @param rowind the row indices
#' @importFrom tibble as_tibble
#' @import Matrix
#' @export
baseline_term <- function(varname, mat, colind, rowind) {
  #print(paste(varname, ":", class(mat)))
  
  stopifnot(inherits(mat, "matrix") || is.data.frame(mat) || inherits(mat, "Matrix"))

  ret <- list(varname=varname, 
              design_matrix=tibble::as_tibble(as.matrix(mat),.name_repair="minimal"), 
              colind=colind, 
              rowind=rowind)
  class(ret) <- c("baseline_term", "matrix_term", "fmri_term", "list")
  ret
}


#' @export
design_matrix.baseline_term <- function(x, blockid=NULL, allrows=FALSE) {
  if (is.null(blockid)) {
    x$design_matrix
  } else {
    if (!allrows) {
      x$design_matrix[unlist(x$rowind[blockid]), unlist(x$colind[blockid]), drop=FALSE]
    } else {
      x$design_matrix[, unlist(x$colind[blockid]), drop=FALSE]
    }
  }
}


#' @keywords internal
construct_block_term <- function(vname, sframe, intercept=c("global", "runwise")) {
  intercept <- match.arg(intercept)
  blockids <- blockids(sframe)
  blockord <- sort(unique(blockids))
  
  expanded_blockids <- ordered(rep(blockord, blocklens(sframe)))
  
  if (length(blockord) == 1 || intercept == "global") {
    mat <- matrix(rep(1, length(expanded_blockids)), length(expanded_blockids), 1)
    colnames(mat) <- paste0(vname, "_", "global")
  } else {
    mat <- model.matrix(~ expanded_blockids - 1)
    colnames(mat) <- paste0(vname, "_", blockord)
  }
  
  
  block_term(vname, blockids=blockids, expanded_blockids=expanded_blockids, mat, type=intercept)
  
}

#' block_term
#' 
#' construct a constant term
#' 
#' @param varname the name of the term
#' @param blockids the ordered sequence of blockids
#' @param expanded_blockids the vector of blockids expanded by run
#' @param mat the \code{matrix} of covariates
#' @param type the block variable type: 'runwise' or 'global'
#' @importFrom tibble as_tibble
#' @export
block_term <- function(varname, blockids, expanded_blockids, mat, type=c("runwise", "global")) {
  type <- match.arg(type)
  assertthat::assert_that(is.matrix(mat))
  assertthat::assert_that(nrow(mat) == length(blockids))
  
  if (type == "runwise") {
    assertthat::assert_that(ncol(mat) == length(unique(blockids)))
  } else {
    assertthat::assert_that(ncol(mat) == 1)
  }
  
  ret <- list(varname=varname, blockids=blockids, expanded_blockids=expanded_blockids, design_matrix=tibble::as_tibble(mat,.name_repair="check_unique"), nblocks=ncol(mat),
              colind=as.list(1:ncol(mat)), rowind=split(1:length(blockids), blockids))
  class(ret) <- c("block_term", "matrix_term", "fmri_term", "list")
  ret
}


#' @export
design_matrix.block_term <- function(x, blockid=NULL, allrows=FALSE) {
  if (is.null(blockid)) {
    x$design_matrix
  } else {
    if (!allrows) {
      x$design_matrix[unlist(x$rowind[blockid]), unlist(x$colind[blockid]), drop=FALSE]
    } else {
      x$design_matrix[, unlist(x$colind[blockid]), drop=FALSE]
    }
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

#' @export
print.baseline_model <- function(x) {
  cat("baseline_model", "\n")
  cat("  ", "name: ", x$drift_term$varname, "\n")
  cat("  ", "basis type: ", x$drift_spec$basis, "\n")
  cat("  ", "degree: ", x$drift_spec$degree, "\n")
  cat("  ", "drift columns: ", ncol(design_matrix(x$drift_term)), "\n")
  if (!is.null(x$block_term)) {
    cat("  ", "constant columns: ", ncol(design_matrix(x$block_term)), "\n")
  } else {
    cat("  ", "constant columns: 0 ", "\n")
  }
  cat("  ", "nuisance columns: ", ifelse(is.null(x$nuisance_term), 0, ncol(design_matrix(x$nuisance_term))), "\n")
  cat("  ", "total columns: ", ncol(design_matrix(x)), "\n")
  cat("  ", "design_matrix: ", "\n")
  print(design_matrix(x))
   
}

#' @importFrom ggplot2 ggplot aes_string geom_line facet_wrap xlab theme_bw
#' @importFrom tidyr gather
#' @export
plot.baseline_model <- function(x, term_name=NULL) {

  all_terms <- terms(x)
  term_names <- names(all_terms)
  
  cidx <- grep("constant", term_names)
  if (length(cidx) > 0) {
    all_terms <- all_terms[-cidx]
    term_names <- term_names[-cidx]
  }
  
  sframe <- x$sampling_frame
  
  condition = NULL; value = NULL; .time = NULL; .block=NULL
  dflist <- lapply(all_terms, function(term) {
    dm1 <- tibble::as_tibble(design_matrix(term),.name_repair="check_unique")
    dm1$.block <- sframe$blockids
    dm1$.time <- sframe$time
    
    tidyr::gather(dm1, condition, value, -.time, -.block)
  })
  
  names(dflist) <- term_names
  
  
  dfx <- if (is.null(term_name)) {
    tn <- term_names[1]
    dflist[[tn]]
  } else {
    dflist[[term_name]]
  }
  
  
  ggplot2::ggplot(dfx, aes_string(x=".time", y="value", colour="condition")) + geom_line() + facet_wrap(~ .block, ncol=1) +
    xlab("Time") + theme_bw(14)
  
}



