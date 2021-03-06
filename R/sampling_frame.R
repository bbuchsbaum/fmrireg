


#' construct a \code{sampling_frame}
#' 
#' A \code{sampling_frame} describes the block structure and temporal sampling of an fMRI paradigm. 
#' 
#' 
#' @param blocklens the number of scans in each block, a \code{vector}
#' @param TR the repetition time in seconds; i.e. the spacing between consecutive image acquisitions.
#' @param start_time the offset of first scan of each block (default is \code{TR/2})
#' @param precision the discrete sampling interval used for convolution with hemodynamic response function.
#' @examples 
#' 
#' frame <- sampling_frame(blocklens=c(100,100, 100), TR=2, precision=.5)
#' 
#' ## the relative time (with respect to the last block) in seconds of each sample/acquisition
#' sam <- samples(frame)
#' ## the global time (with respect to first block) of each sample/acquisition
#' gsam <- samples(frame, global=TRUE)
#' @export
sampling_frame <- function(blocklens, TR, start_time=TR/2, precision=.1) {
  assert_that(all(TR > 0))
  assert_that(all(blocklens > 0))
  
  if (length(TR) == 1) {
    TR <- rep(TR, length(blocklens))
  }
  
  if (length(start_time) == 1) {
    start_time <- rep(start_time, length(blocklens))
  }
  
  blockids <- rep(1:length(blocklens), blocklens)
  scan_time <- unlist(lapply(1:length(blocklens), function(i) seq(start_time[i], by=TR[i], length.out=blocklens[i])))
  
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
## TODO screwy things happen when blockids don't start at 1.
samples.sampling_frame <- function(x, blockids=NULL, global=FALSE) {
  if (is.null(blockids)) {
    blockids <- seq(1, length(x$blocklens))
  }
  
  if (!global) {
    unlist(lapply(blockids, function(b) {
      seq(x$start_time[b], by=x$TR[b], length.out=x$blocklens[b])
    }))
  } else {
    unlist(lapply(blockids, function(b) {
      start <- if (b > 1) {
        sum(x$blocklens[1:(b-1)])*x$TR[b] + x$start_time[b] 
      } else {
        x$start_time[b]
      }
      
      seq(start, by=x$TR[b], length.out=x$blocklens[b])
    }))
  }
}


#' @export
global_onsets.sampling_frame <- function(x, onsets, blockids) {
  
  ids <- rep(1:length(unique(blockids)), table(blockids))
  
  if (max(ids) > length(x$blocklens)) {
    stop("there are more block ids than block lengths, cannot compute global onsets")
  }
  
  map_dbl(1:length(onsets),function(i) {
    blocknum <- ids[i]
    offset <- (sum(x$blocklens[1:blocknum]) - x$blocklens[blocknum])*x$TR[blocknum]
    if (onsets[i] > x$blocklens[blocknum]*x$TR[blocknum]) {
      NA
    } else {
      onsets[i] + offset
    }
    
  })
}  

#' @export
split_by_block.sampling_frame <- function(x, vals) {
  split(vals, x$blockids)
}

#' @export
blockids.sampling_frame <- function(x) {
  x$blockids
}

#' @export
blocklens.sampling_frame <- function(x) {
  x$blocklens
}

#' @export
print.sampling_frame <- function(x) {
  cat("sampling_frame: \n")
  cat("  number of blocks:", length(x$blocklens), "\n")
  cat("  blocklens: ", paste(x$blocklens, collapse=", "), "\n")
  cat("  TR: ", paste0(x$TR, "s"), "\n")
  cat("  start_time: ", paste0(x$start_time, "s"), "\n")
}
