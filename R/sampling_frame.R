#' Construct a sampling_frame
#'
#' A \code{sampling_frame} describes the block structure and temporal sampling of an fMRI paradigm.
#'
#' @param blocklens A numeric vector representing the number of scans in each block.
#' @param TR A numeric value or vector representing the repetition time in seconds (i.e., the spacing between consecutive image acquisitions).
#' @param start_time A numeric value or vector representing the offset of the first scan of each block (default is \code{TR/2}).
#' @param precision A numeric value representing the discrete sampling interval used for convolution with the hemodynamic response function (default is 0.1).
#'
#' @examples
#' frame <- sampling_frame(blocklens = c(100, 100, 100), TR = 2, precision = 0.5)
#'
#' # The relative time (with respect to the last block) in seconds of each sample/acquisition
#' sam <- samples(frame)
#' # The global time (with respect to the first block) of each sample/acquisition
#' gsam <- samples(frame, global = TRUE)
#'
#' @return A list with class "sampling_frame" describing the block structure and temporal sampling of an fMRI paradigm.
#' @export
sampling_frame <- function(blocklens, TR, start_time=TR/2, precision=.1) {
  # Add comprehensive input validation
  assert_that(all(TR > 0), msg="TR (repetition time) must be positive")
  assert_that(all(blocklens > 0), msg="Block lengths must be positive integers")
  assert_that(precision > 0 && precision < min(TR), 
              msg="Precision must be positive and less than the minimum TR")
  assert_that(all(start_time >= 0), 
              msg="Start time must be non-negative")
  
  # Add warnings for potential issues
  if (length(blocklens) > 50) {
    warning("There are more than 50 blocks in sampling_frame")
  }
  
  total_samples <- sum(blocklens)
  if (total_samples > 1e6) {
    warning("Large number of samples (>1M) may impact performance")
  }
  
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

#' Extract samples from a sampling_frame
#'
#' This function extracts the relative or global time of each sample/acquisition from a \code{sampling_frame}.
#'
#' @param x A sampling_frame object.
#' @param blockids A numeric vector of block IDs to extract the samples from. If NULL (default), all block IDs are used.
#' @param global A logical value. If TRUE, the global time (with respect to the first block) of each sample/acquisition is returned. If FALSE (default), the relative time (with respect to the last block) of each sample/acquisition is returned.
#' @param ... Additional arguments (currently unused).
#'
#' @examples
#' frame <- sampling_frame(blocklens = c(100, 100, 100), TR = 2, precision = 0.5)
#' # The relative time (with respect to the last block) in seconds of each sample/acquisition
#' sam <- samples(frame)
#' # The global time (with respect to the first block) of each sample/acquisition
#' gsam <- samples(frame, global = TRUE)
#'
#' @return A numeric vector of sample times extracted from the specified \code{sampling_frame}.
#' @export
samples.sampling_frame <- memoise::memoise(function(x, blockids=NULL, global=FALSE, ...) {
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
})

#' Compute global onsets from a sampling_frame
#'
#' This function computes the global onsets (with respect to the first block) for a given \code{sampling_frame}.
#'
#' @param x A sampling_frame object.
#' @param onsets A numeric vector of onsets within the specified blocks.
#' @param blockids A numeric vector of block IDs corresponding to the onsets.
#' @param ... Additional arguments (currently unused).
#'
#' @examples
#' frame <- sampling_frame(blocklens = c(100, 100, 100), TR = 2, precision = 0.5)
#' onsets <- c(10, 20, 30)
#' blockids <- c(1, 2, 3)
#' global_onsets(frame, onsets, blockids)
#'
#' @return A numeric vector of global onsets computed from the specified \code{sampling_frame}.
#' @export
#' @family global_onsets
global_onsets.sampling_frame <- memoise::memoise(function(x, onsets, blockids,...) {
  
  ids <- rep(1:length(unique(blockids)), table(blockids))
  
  if (max(ids) > length(x$blocklens)) {
    stop("there are more block ids than block lengths, cannot compute global onsets")
  }
  
  purrr::map_dbl(1:length(onsets),function(i) {
    blocknum <- ids[i]
    offset <- (sum(x$blocklens[1:blocknum]) - x$blocklens[blocknum])*x$TR[blocknum]
    if (onsets[i] > x$blocklens[blocknum]*x$TR[blocknum]) {
      NA
    } else {
      onsets[i] + offset
    }
    
  })
})  

#' @export
split_by_block.sampling_frame <- function(x, vals, ...) {
  split(vals, x$blockids)
}

#' @export
blockids.sampling_frame <- function(x) {
  x$blockids
}

#' @export
blocklens.sampling_frame <- function(x,...) {
  x$blocklens
}

#' @export
print.sampling_frame <- function(x, ...) {
  # Header
  cat("\n═══ Sampling Frame ═══\n")
  
  # Structure information
  cat("\n📊 Structure:\n")
  cat(crayon::blue("  • Blocks:"), length(x$blocklens), "\n")
  cat(crayon::blue("  • Total timepoints:"), sum(x$blocklens), "\n")
  cat(crayon::blue("  • Block lengths:"), paste(x$blocklens, collapse=", "), "\n")
  
  # Timing information
  cat("\n⏱️  Timing:\n")
  cat(crayon::blue("  • TR:"), paste(unique(x$TR), collapse=", "), "seconds\n")
  cat(crayon::blue("  • Start times:"), paste(unique(x$start_time), collapse=", "), "seconds\n")
  cat(crayon::blue("  • Precision:"), x$precision, "seconds\n")
  
  # Duration information
  cat("\n📈 Duration:\n")
  cat(crayon::blue("  • Total duration:"), sum(x$blocklens * x$TR), "seconds\n")
  
  cat("\n")
}
