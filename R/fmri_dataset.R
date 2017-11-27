


#' read_fmri_config
#' 
#' @param file_name name of configuration file
#' @importFrom assertthat assert_that
#' @importFrom tibble as_data_frame
#' @export
read_fmri_config <- function(file_name, base_path=NULL) {
  print(file_name)
  env <- new.env()
  source(file_name, env)
  
  if (is.null(env$base_path) && is.null(base_path)) {
    env$base_path = "."
  } else {
    env$base_path <- base_path
  }
  
  if (is.null(env$output_dir)) {
    env$output_dir = "stat_out"
  }
  

  assert_that(!is.null(env$scans))
  assert_that(!is.null(env$TR))
  assert_that(!is.null(env$mask))
  assert_that(!is.null(env$run_length))
  assert_that(!is.null(env$event_model))
  assert_that(!is.null(env$event_table))
  assert_that(!is.null(env$block_column))
  #assert_that(file.exists(file.path(env$base_path,env$event_table)))
  
  
  #env$mask <- neuroim::loadVolume(file.path(env$base_path, env$mask))
  print(env$base_path)
 
  dname <- file.path(env$base_path, env$event_table)
  print(dname)
  assert_that(file.exists(dname))
  env$design <- tibble::as_data_frame(read.table(dname, header=TRUE))

  # if (is.null(env$aux_table)) {
  #   env$aux_table=tibble::as_data_frame()
  # } else {
  #   env$aux_table <- tibble::as_data_frame(read.table(file.path(env$base_path,env$aux_table), header=TRUE))
  # }
  
  out <- as.list(env)
  class(out) <- c("fmri_config", "list")
  out
}


#' @export
matrix_dataset <- function(datamat, TR, run_length, event_table=data.frame()) {
  assert_that(sum(run_length) == nrow(datamat))
  
  frame <- sampling_frame(run_length, TR)
  
  ret <- list(
    datamat=datamat,
    nruns=length(run_length),
    event_table=event_table,
    sampling_frame=frame
  )
  
  class(ret) <- c("matrix_dataset", "fmri_dataset", "list")
  ret
  
}

#' fmri_mem_dataset
#' 
#' @param scans a vector of objects of class \code{\linkS4class{BrainVector}}
#' @param mask a object of type \code{\linkS4class{LogicalBrainVolume}} indicating the voxels to include in analyses.
#' @param TR the repetition time in seconds of the scan-to-scan interval.
#' @param event_table a \code{data.frame} containing the event onsets and experimental variables.
#' @export
fmri_mem_dataset <- function(scans, mask, TR, 
                         event_table=data.frame(), 
                         base_path=".") {
  
  assert_that(all(sapply(scans, function(x) inherits(x, "BrainVector"))))
  
  
  run_length <- sapply(scans, function(x) dim(x)[4])

  frame <- sampling_frame(run_length, TR)
  
  assert_that(length(run_length) == length(scans))
  
  ret <- list(
    scans=scans,
    mask=mask,
    nruns=length(scans),
    event_table=event_table,
    base_path=base_path,
    sampling_frame=frame
  )
  
  class(ret) <- c("fmri_mem_dataset", "volumetric_dataset", "fmri_dataset", "list")
  ret
}


#' fmri_dataset
#' 
#' @param scans a vector of file names of the images comprising the dataset
#' @param mask name of the binary mask file indicating the voxels to include in analysis.
#' @param TR the repetition time in seconds of the scan-to-scan interval.
#' @param run_length the number of scans in each run.
#' @param event_table a \code{data.frame} containing the event onsets and experimental variables.
#' @export
fmri_dataset <- function(scans, mask, TR, 
                         run_length, 
                         event_table=data.frame(), 
                         base_path=".") {
  
  if (length(run_length) == 1) {
    run_length <- rep(run_length, length(scans))
  }
  
  frame <- sampling_frame(run_length, TR)
  assert_that(length(run_length) == length(scans))
  
  ret <- list(
    scans=scans,
    mask_file=mask,
    nruns=length(scans),
    event_table=event_table,
    base_path=base_path,
    sampling_frame=frame
  )
  
  class(ret) <- c("fmri_file_dataset", "volumetric_dataset", "fmri_dataset", "list")
  ret
}

#' @keywords internal
data_chunk <- function(mat, voxel_ind, row_ind, chunk_num) {
  ret <- list(data=mat,
       voxel_ind=voxel_ind,
       row_ind=row_ind,
       chunk_num=chunk_num)
  class(ret) <- c("data_chunk", "list")
  ret
}

#' @keywords internal
chunk_iter <- function(x, nchunks, get_chunk) {
  chunk_num <- 1
  
  nextEl <- function() {
    if (chunk_num > nchunks) {
      stop("StopIteration")
    } else {
      ret <- get_chunk(chunk_num)
      chunk_num <<- chunk_num + 1
      ret
    }
  }
  
  iter <- list(nextElem=nextEl)
  class(iter) <- c("chunkiter", "abstractiter", "iter")
  iter
}


#' @importFrom neuroim series
#' @export
data_chunks.fmri_mem_dataset <- function(x, nchunks=1,runwise=FALSE) {
  
  mask <- x$mask
  
  get_run_chunk <- function(chunk_num) {
    bvec <- x$scans[[chunk_num]]
    voxel_ind <- which(x$mask>0)
    row_ind=which(x$sampling_frame$blockids == chunk_num)
    ret <- data_chunk(neuroim::series(bvec,voxel_ind), 
                      voxel_ind=voxel_ind, 
                      row_ind=row_ind, 
                      chunk_num=chunk_num)
  }
  
  get_seq_chunk <- function(chunk_num) {
    bvecs <- x$scans
    voxel_ind <- maskSeq[[chunk_num]]
    ret <- data_chunk(do.call(rbind, lapply(bvecs, function(bv) neuroim::series(bv, voxel_ind))), 
                      voxel_ind=voxel_ind, 
                      row_ind=1:nrow(x$event_table))
    
  }
  
  if (runwise) {
    chunk_iter(x, length(x$scans), get_run_chunk)
  } else if (nchunks == 1) {
    maskSeq <- one_chunk()
    chunk_iter(x, 1, get_seq_chunk)
  } else if (nchunks == dim(mask)[3]) {
    maskSeq <- slicewise_chunks(x)
    chunk_iter(x, length(maskSeq), get_seq_chunk)
  } else {
    maskSeq <- arbitrary_chunks(x, nchunks)
    chunk_iter(x, length(maskSeq), get_seq_chunk)
  }
  
}



#' @import neuroim
#' @export
data_chunks.fmri_dataset <- function(x, nchunks=1,runwise=FALSE) {
  
  mask <- x$mask
  
  get_run_chunk <- function(chunk_num) {
    bvec <- neuroim::loadVector(file.path(x$base_path, x$scans[chunk_num]), mask=x$mask)
    ret <- data_chunk(bvec@data, voxel_ind=which(x$mask>0), row_ind=which(x$sampling_frame$blockids == chunk_num), chunk_num=chunk_num)
  }
  
  get_seq_chunk <- function(chunk_num) {
    bvecs <- lapply(x$scans, function(scan) neuroim::loadVector(file.path(x$base_path, scan), mask=maskSeq[[chunk_num]]))
    ret <- data_chunk(do.call(rbind, lapply(bvecs, function(bv) bv@data)), voxel_ind=maskSeq[[chunk_num]], 
                      row_ind=1:nrow(x$event_table))
    
  }
  
  if (runwise) {
    chunk_iter(x, length(x$scans), get_run_chunk)
  } else if (nchunks == 1) {
    maskSeq <- one_chunk()
    chunk_iter(x, 1, get_seq_chunk)
  } else if (nchunks == dim(mask)[3]) {
    maskSeq <- slicewise_chunks(x)
    chunk_iter(x, length(maskSeq), get_seq_chunk)
  } else {
    maskSeq <- arbitrary_chunks(x, nchunks)
    chunk_iter(x, length(maskSeq), get_seq_chunk)
  }
  
}

#' @import neuroim
#' @export
data_chunks.matrix_dataset <- function(x, runwise=TRUE) {
  get_run_chunk <- function(chunk_num) {
    ind <- which(blockids(x$sampling_frame) == chunk_num)
    mat <- x$datamat[ind,]
    data_chunk(mat, voxel_ind=1:ncol(mat), row_ind=ind, chunk_num=chunk_num)
  }
  
  get_one_chunk <- function(chunk_num) {
    data_chunk(x$datamat, voxel_ind=1:ncol(x$datamat), row_ind=1:nrow(x$datamat), chunk_num=chunk_num)
    
  }
  
  if (runwise) {
    chunk_iter(x, length(x$sampling_frame$blocklens), get_run_chunk)
  } else {
    chunk_iter(x, 1, get_one_chunk)
  } 
  
}

#' @export
exec_strategy <- function(strategy=c("all", "slicewise", "runwise")) {
  strategy <- match.arg(strategy)
  function(dset) {
    if (strategy == "runwise") {
      data_chunks(dset, runwise=TRUE)
    } else if (strategy == "slicewise") {
      data_chunks(dset, nchunks = dim(dset$mask)[3], runwise=FALSE)
    } else if (strategy == "all") {
      data_chunks(dset, nchunks = 1, runwise=FALSE)
    }
  }
  
}


#' @keywords internal
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
  
  maskSeq
  
}

#' @keywords internal
slicewise_chunks <- function(x) {
  mask <- x$mask
  template <- BrainVolume(array(0, dim(mask)), space(mask))
  nchunks <- dim(mask)[3]
  
  maskSeq <- lapply(1:nchunks, function(i) {
    m <- template
    m[,,i] <- 1
    m
  })
  
  maskSeq
  
}

one_chunk <- function(x) {
  mask <- x$mask
  list(mask)
}





