#' @param blocklens
#' @param TR
#' @param start_time
#' @export
sampling_frame <- function(blocklens, TR, start_time=TR/2, precision=.1) {
  assert_that(TR > 0)
  assert_that(all(blocklens > 0))
  
  blockids <- rep(1:length(blocklens), blocklens)
  scan_time <- unlist(lapply(1:length(blocklens), function(i) seq(start_time, by=TR, length.out=blocklens[i])))
  ret <- list(blocklens=blocklens,
              TR=TR,
              start_time=start_time,
              blockids=blockids,
              time=scan_time,
              precision=precision)
  
  class(ret) <- c("sampling_frame", "list")
  ret
}

#' @export
samples.sampling_frame <- function(x, blocknum=NULL, global=FALSE) {
  if (is.null(blocknum)) {
    blocknum <- seq(1, length(x$blocklens))
  }
  
  if (!global) {
    x$time
  } else {
    unlist(lapply(blocknum, function(b) {
      start <- if (b > 1) sum(x$blocklens[1:(b-1)])*x$TR + x$start_time else x$start_time
      seq(start, by=x$TR, length.out=x$blocklens[b])
    }))
  }
}


#' @export
global_onsets.sampling_frame <- function(x, onsets, blockids) {
  
  ids <- rep(1:length(unique(blockids)), table(blockids))
  
  if (max(ids) > length(x$blocklens)) {
    stop("there are more block ids than block lengths, cannot compute global onsets")
  }
  
  sapply(1:length(onsets),function(i) {
    blocknum <- ids[i]
    offset <- (sum(x$blocklens[1:blocknum]) - x$blocklens[blocknum])*x$TR
    if (onsets[i] > x$blocklens[blocknum]*x$TR) {
      NA
    } else {
      onsets[i] + offset
    }
    
  })
}  

split_by_block.sampling_frame <- function(x, vals) {
  split(vals, x$blockids)
}
