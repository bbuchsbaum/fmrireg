#' fmrireg: regression tools for fMRI data
#'
#' fmrireg provides functions for generating experimental design matrices appropriate for analyzing fMRI data with regression.
#' 
#' 
#' @keywords internal
"_PACKAGE"

#' @useDynLib fmrireg
#' @importFrom Rcpp evalCpp
#' @import assertthat
#' @importFrom stats as.formula
NULL


#if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))