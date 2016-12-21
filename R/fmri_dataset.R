#' read_fmri_config
#' 
#' @param file_name name of configuration file
#' @importFrom assertthat assert_that
#' @importFrom tibble as_data_frame
#' @export
read_fmri_config <- function(file_name) {
  env <- new.env()
  source(file_name, env)
  
  if (is.null(env$base_path)) {
    env$base_path = ""
  }
  
  if (is.null(env$output_dir)) {
    env$output_dir = "stat_out"
  }
  
  
  
  assert_that(!is.null(env$scans))
  assert_that(!is.null(env$TR))
  assert_that(!is.null(env$mask))
  assert_that(!is.null(env$run_length))
  assert_that(!is.null(env$event_model))
  assert_that(!is.null(env$design))
  assert_that(file.exists(file.path(env$base_path,env$design)))
  
  
  #env$mask <- neuroim::loadVolume(file.path(env$base_path, env$mask))
  env$design <- tibble::as_data_frame(read.table(file.path(env$base_path,env$design), header=TRUE))
  if (is.null(env$aux_data)) {
    env$aux_data=tibble::as_data_frame()
  } else {
    env$aux_data <- tibble::as_data_frame(read.table(file.path(env$base_path,env$aux_data), header=TRUE))
  }
  
  out <- as.list(env)
  class(out) <- c("fmri_config", "list")
  out
}



#' fmri_dataset
#' 
#' @param scans a vector of file names of the images comprising the dataset
#' @param mask name of the binary mask file indicating the voxels to include in analysis.
#' @param TR the repetition time in seconds of the scan-to-scan interval.
#' @param run_length the number of scans in each run.
#' @param runids the id of each block, must be unique for each run and be non-decreasing.
#' @param event_table a \code{data.frame} containing the event onsets and experimental variables.
#' @param aux_data a \code{list} of auxilliary data such as nuisance variables that may enter a 
#' regression model and are not convolved with hemodynamic response function.
#' @export
fmri_dataset <- function(scans, mask, TR, 
                         run_length, 
                         event_table=data.frame(), 
                         aux_data=tibble::as_tibble(list()), base_path="") {
  
  if (length(run_length) == 1) {
    run_length <- rep(run_length, length(scans))
  }
  
  assert_that(length(run_length) == length(scans))
  runids <- rep(1:length(scans), run_length)
  
  ret <- list(
    scans=scans,
    mask_file=mask,
    mask=neuroim::loadVolume(file.path(base_path, mask)),
    TR=TR,
    run_length=run_length,
    runids=runids,
    nruns=length(scans),
    event_table=event_table,
    aux_data=aux_data,
    base_path=base_path
  )
  
  class(ret) <- c("fmri_dataset", "list")
  ret
}


data_chunk <- function(mat, voxel_ind, row_ind) {
  ret <- list(data=mat,
       voxel_ind=voxel_ind,
       row_ind=row_ind)
  class(ret) <- c("data_chunk", "list")
  ret
}

gen_chunk_iter <- function(x, maskSeq) {
  nchunks <- length(maskSeq)
  chunkNum <- 1
  
  nextEl <- function() {
    if (chunkNum > nchunks) {
      stop("StopIteration")
    } else {
      bvecs <- lapply(x$scans, function(scan) neuroim::loadVector(file.path(x$base_path, scan), mask=maskSeq[[chunkNum]]))
      ret <- data_chunk(do.call(rbind, lapply(bvecs, function(bv) bv@data)), voxel_ind=maskSeq[[chunkNum]], row_ind=1:nrow(x$event_table))
      chunkNum <<- chunkNum + 1
      ret
    }
  }
  
  iter <- list(nextElem=nextEl)
  class(iter) <- c("chunkiter", "abstractiter", "iter")
  iter
}
  
arbitrary_chunks <- function(x, nchunks) {
  mask <- x$mask
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

runwise_chunks <- function(x) {
  nchunks <- length(x$scans)
  chunkNum <- 1
  
  nextEl <- function() {
    if (chunkNum > nchunks) {
      stop("StopIteration")
    } else {
      bvec <- neuroim::loadVector(file.path(x$base_path, x$scans[chunkNum]), mask=x$mask)
      ret <- data_chunk(bvec@data, voxel_ind=which(x$mask>0), row_ind=which(x$runids==chunkNum))
      chunkNum <<- chunkNum + 1
      ret
    }
  }
  
  iter <- list(nextElem=nextEl)
  class(iter) <- c("chunkiter", "abstractiter", "iter")
  iter
}

slicewise_chunks <- function(x) {
  mask <- x$mask
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
  mask <- x$mask
  gen_chunk_iter(x, list(mask))

}

#' @import neuroim
data_chunks.fmri_dataset <- function(x, nchunks=1,runwise=FALSE) {
  
  mask <- x$mask
  
  if (runwise) {
    runwise_chunks(x)
  } else if (nchunks == 1) {
    one_chunk(x)
  } else if (nchunks == dim(mask)[3]) {
    slicewise_chunks(x)
  } else {
    arbitrary_chunks(x, nchunks)
  }
}




