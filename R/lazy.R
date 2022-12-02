
#' #' @keywords internal
#' as_vectors.matrix <- function(x, mask=NULL) {
#'   if (!is.null(mask)) {
#'     assert_that(length(mask) == ncol(x))
#'     x <- x[,as.logical(mask)]
#'   }
#'   
#'   f <- function(i) {
#'     x[,i]
#'   }
#'   
#'   neuroim2::deferred_list(replicate(ncol(x), f))
#' }
#' 
#' #' @keywords internal
#' as_vectors.NeuroVec <- function(x, mask=NULL) {
#'   if (!is.null(mask)) {
#'     assert_that(length(mask) == dim(x)[4])
#'     mask <- which(as.logical(mask))
#'   } else {
#'     mask <- 1:(dim(x)[4])
#'   }
#'   
#'   f <- function(i) {
#'     drop(series(x, mask[i]))
#'   }
#'   
#'   neuroim2::deferred_list(replicate(length(mask), f))
#' }
#' 
#' 
#' #' @keywords internal
#' lazy_series <- function(bvec, i) {
#'   structure(list(
#'     bvec=bvec,
#'     i=i),
#'     class="lazy_series")
#' }
#' 
#' #' @keywords internal
#' as.vector.lazy_series <- function(x) {
#'   series(x$bvec, x$i)
#' }
#' 
#' 
#' #' @keywords internal
#' to_tibble.fmri_mem_dataset <- function(x) {
#'   stop("not implemented")
#' }
