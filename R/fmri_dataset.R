


#' @param scans a vector of file names for the images comprising the dataset
#' @param TR the repetition time in seconds of the scan-to-scan interval.
#' @param blocklens the number of scans in each block
#' @param blockids the id of each block, must be unique for each block and be non-decreasing
#' @param event_table a \code{data.frame} containing the event onsets and experimental variables.
#' @param aux_data a \code{list} of auxilliary data such as nuisance variables that may enter a 
#' regression model unconvolved with hemodynamic response function.
#' @export
fmri_dataset <- function(scans, TR, blocklens, blockids=rep(1:length(blocklens), blocklens), event_table=data.frame(), aux_data=list()) {
  
  ret <- list(
    scans=scans,
    TR=TR,
    blocklens=blocklens,
    blockids=blockids,
    event_table=event_table,
    aux_data=aux_data
  )
  
  class(ret) <- c("fmri_dataset", "list")
  ret
}