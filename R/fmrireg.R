#' fmrireg: regresssion tools for fMRI data
#'
#' fmrireg provides functions for generating experimental design matrices appropriate for analyzing fMRI data with regression.
#' 
#' 
#' @useDynLib fmrireg
#' @importFrom Rcpp evalCpp
#' @docType package
#' @name fmrireg
#' @import assertthat
#' @import stats
NULL


#if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))