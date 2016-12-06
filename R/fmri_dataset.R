

#' fmri_dataset
#' 
#' @param scans a vector of file names of the images comprising the dataset
#' @param mask name of the binary mask file indicating the voxels to include in analysis.
#' @param TR the repetition time in seconds of the scan-to-scan interval.
#' @param blocklens the number of scans in each block.
#' @param blockids the id of each block, must be unique for each block and be non-decreasing.
#' @param event_table a \code{data.frame} containing the event onsets and experimental variables.
#' @param aux_data a \code{list} of auxilliary data such as nuisance variables that may enter a 
#' regression model and are not convolved with hemodynamic response function.
#' @export
fmri_dataset <- function(scans, mask, TR, 
                         blocklens, blockids=rep(1:length(blocklens), blocklens), 
                         event_table=data.frame(), 
                         aux_data=data.frame()) {
  assert_that(length(unique(blockids)) == length(blocklens))
  
  ret <- list(
    scans=scans,
    mask=mask,
    TR=TR,
    blocklens=blocklens,
    blockids=blockids,
    event_table=event_table,
    aux_data=aux_data
  )
  
  class(ret) <- c("fmri_dataset", "list")
  ret
}



gen_chunk_iter <- function(x, maskSeq) {
  nchunks <- length(maskSeq)
  chunkNum <- 1
  
  nextEl <- function() {
    if (chunkNum > nchunks) {
      stop("StopIteration")
    } else {
      bvecs <- lapply(x$scans, function(scan) neuroim::loadVector(scan, mask=maskSeq[[chunkNum]]))
      chunkNum <<- chunkNum + 1
      do.call(rbind, lapply(bvecs, function(bv) bv@data))
    }
  }
  
  iter <- list(nextElem=nextEl)
  class(iter) <- c("chunkiter", "abstractiter", "iter")
  iter
}
  
arbitrary_chunks <- function(x, nchunks) {
  mask <- neuroim::loadVolume(x$mask)
  indices <- which(mask != 0)
  chsize <- as.integer(length(indices)/nchunks)
  chunkids <- rep(1:nchunks, each=chsize, length.out=length(indices))
  
  maskSeq <- lapply(1:nchunks, function(i) {
    m <- mask
    m[] <- 0
    m[indices[chunkids==i]] <- 1
    m
  })
  
  gen_chunk_iter(x, maskSeq)
  
}

slicewise_chunks <- function(x) {
  mask <- neuroim::loadVolume(x$mask)
  template <- BrainVolume(array(0, dim(mask)), space(mask))
  nchunks <- dim(mask)[3]
  
  maskSeq <- lapply(1:nchunks, function(i) {
    m <- template
    m[,,i] <- 1
    m
  })
  
  gen_chunk_iter(x, maskSeq)
  
}

one_chunk <- function(x) {
  mask <- neuroim::loadVolume(x$mask)
  gen_chunk_iter(x, list(mask))

}

#' @import neuroim
data_chunks.fmri_dataset <- function(x, nchunks=1, lazy=TRUE) {
  
  mask <- neuroim::loadVolume(x$mask)
  if (nchunks == 1) {
    one_chunk(x)
  } else if (nchunks == dim(mask)[3]) {
    slicewise_chunks(x)
  } else {
    arbitrary_chunks(x, nchunks)
  }
}


