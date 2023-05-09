`%dopar%` <- foreach::`%dopar%`
`%do%` <- foreach::`%do%`

#' @keywords internal
default_config <- function() {
  env <- new.env()
  env$cmd_flags <- ""
  env$jobs <- 1
  env
  
}


#' read a basic fMRI configuration file
#' 
#' @param file_name name of configuration file
#' @param base_path the file path to be prepended to relative file names
#' @importFrom assertthat assert_that
#' @importFrom tibble as_tibble
#' @export
read_fmri_config <- function(file_name, base_path=NULL) {
  #print(file_name)
  env <- default_config()
  
  source(file_name, env)
  
  env$base_path <- if (is.null(env$base_path) && is.null(base_path)) {
   "."
  } else if (!is.null(base_path) && is.null(env$base_path)) {
    base_path
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
  assert_that(!is.null(env$baseline_model))
  
  if (!is.null(env$censor_file)) {
    env$censor_file = NULL
  }
  
  if (!is.null(env$contrasts)) {
    env$contrasts = NULL
  }
  
  if (!is.null(env$nuisance)) {
    env$nuisance = NULL
  }
  
  dname <- file.path(env$base_path, env$event_table)
  
  assert_that(file.exists(dname))
  env$design <- suppressMessages(tibble::as_tibble(read.table(dname, header=TRUE),.name_repair="check_unique"))

  out <- as.list(env)
  class(out) <- c("fmri_config", "list")
  out
}

#' Create a Matrix Dataset Object
#'
#' This function creates a matrix dataset object, which is a list containing information about the data matrix, TR, number of runs, event table, sampling frame, and mask.
#'
#' @param datamat A matrix where each column is a voxel time-series.
#' @param TR Repetition time (TR) of the fMRI acquisition.
#' @param run_length A numeric vector specifying the length of each run in the dataset.
#' @param event_table An optional data frame containing event information. Default is an empty data frame.
#'
#' @return A matrix dataset object of class c("matrix_dataset", "fmri_dataset", "list").
#' @export
#'
#' @examples
#' # A matrix with 100 rows and 100 columns (voxels)
#' X <- matrix(rnorm(100*100), 100, 100)
#' dset <- matrix_dataset(X, TR=2, run_length=100)
#'
#' # An iterator with 5 chunks
#' iter <- data_chunks(dset, nchunks=5)
#' `%do%` <- foreach::`%do%`
#' y <- foreach::foreach(chunk = iter) %do% { colMeans(chunk$data) }
#' length(y) == 5
#'
#' # An iterator with 100 chunks
#' iter <- data_chunks(dset, nchunks=100)
#' y <- foreach::foreach(chunk = iter) %do% { colMeans(chunk$data) }
#' length(y) == 100
#'
#' # A matrix_dataset with 200 rows, 100 columns and 2 runs
#' X <- matrix(rnorm(200*100), 200, 100)
#' dset <- matrix_dataset(X, TR=2, run_length=c(100,100))
#'
#' # Get a "runwise" iterator. For every iteration, an entire run's worth of data is returned.
#' iter <- data_chunks(dset, runwise=TRUE)
#' y <- foreach::foreach(chunk = iter)  %do% { colMeans(chunk$data) }
#' length(y) == 2
matrix_dataset <- function(datamat, TR, run_length, event_table=data.frame()) {
  if (is.vector(datamat)) {
    datamat <- as.matrix(datamat)
  }
  assert_that(sum(run_length) == nrow(datamat))
  
  frame <- sampling_frame(run_length, TR)
  
  ret <- list(
    datamat=datamat,
    TR=TR,
    nruns=length(run_length),
    event_table=event_table,
    sampling_frame=frame,
    mask=rep(1,ncol(datamat))
  )
  
  class(ret) <- c("matrix_dataset", "fmri_dataset", "list")
  ret
  
}

#' Create an fMRI Memory Dataset Object
#'
#' This function creates an fMRI memory dataset object, which is a list containing information about the scans, mask, TR, number of runs, event table, base path, sampling frame, and censor.
#'
#' @param scans A list of objects of class \code{\linkS4class{NeuroVec}}.
#' @param mask A binary mask of class \code{\linkS4class{NeuroVol}} indicating the set of voxels to include in analyses.
#' @param TR Repetition time (TR) of the fMRI acquisition.
#' @param run_length A numeric vector specifying the length of each run in the dataset. Default is the length of the scans.
#' @param event_table An optional data frame containing event information. Default is an empty data frame.
#' @param base_path An optional base path for the dataset. Default is "." (current directory).
#' @param censor An optional numeric vector specifying which time points to censor. Default is NULL.
#'
#' @return An fMRI memory dataset object of class c("fmri_mem_dataset", "volumetric_dataset", "fmri_dataset", "list").
#' @export
#'
#' @examples
#' # Create a NeuroVec object
#' d <- c(10, 10, 10, 10)
#' nvec <- neuroim2::NeuroVec(array(rnorm(prod(d)), d), space=neuroim2::NeuroSpace(d))
#'
#' # Create a NeuroVol mask
#' mask <- neuroim2::NeuroVol(array(rnorm(10*10*10), d[1:3]), space=neuroim2::NeuroSpace(d[1:3]))
#' mask[mask < .5] <- 0
#'
#' # Create an fmri_mem_dataset
#' dset <- fmri_mem_dataset(list(nvec), mask, TR=2)
#'
#' # Create an iterator with 100 chunks
#' iter <- data_chunks(dset, nchunks=100)
fmri_mem_dataset <- function(scans, mask, TR, 
                             run_length=sapply(scans, function(x) dim(x)[4]),
                             event_table=data.frame(), 
                             base_path=".",
                             censor=NULL) {
  
  
  
  assert_that(all(map_lgl(scans, function(x) inherits(x, "NeuroVec"))))
  assert_that(inherits(mask, "NeuroVol"))
  assert_that(all(dim(mask) == dim(scans[[1]][1:3])))
  
  ntotscans <- sum(sapply(scans, function(x) dim(x)[4]))
  #run_length <- map_dbl(scans, ~ dim(.)[4])
  assert_that(sum(run_length) == ntotscans)
  
  if (is.null(censor)) {
    censor <- rep(0, sum(run_length))
  }

  frame <- sampling_frame(run_length, TR)
  
  ret <- list(
    scans=scans,
    mask=mask,
    nruns=length(run_length),
    event_table=event_table,
    base_path=base_path,
    sampling_frame=frame,
    censor=censor
  )
  
  class(ret) <- c("fmri_mem_dataset", "volumetric_dataset", "fmri_dataset", "list")
  ret
}

#' Create a Latent Dataset Object
#'
#' This function creates a latent dataset object, which encapsulates a dimension-reduced subspace of "latent variables".
#' The dataset is a list containing information about the latent neuroimaging vector, TR, number of runs, event table, base path, sampling frame, and censor.
#'
#' @param lvec An instance of class \code{LatentNeuroVec}.
#' @param TR Repetition time (TR) of the fMRI acquisition.
#' @param run_length A numeric vector specifying the length of each run in the dataset.
#' @param event_table An optional data frame containing event information. Default is an empty data frame.
#' @param base_path An optional base path for the dataset. Default is "." (current directory).
#' @param censor An optional numeric vector specifying which time points to censor. Default is NULL.
#'
#' @return A latent dataset object of class c("latent_dataset", "fmri_dataset", "list").
#' @export
#'
#' @examples
#' # Create a matrix with 100 rows and 1000 columns (voxels)
#' X <- matrix(rnorm(100*1000), 100, 1000)
#' pres <- prcomp(X)
#' basis <- pres$x[,1:25]
#' loadings <- pres$rotation[,1:25]
#' offset <- colMeans(X)
#'
#' # Create a LatentNeuroVec object
#' lvec <- neuroim2::LatentNeuroVec(basis, loadings, neuroim2::NeuroSpace(c(10,10,10,100)), 
#' mask=rep(TRUE,1000), offset=offset)
#'
#' # Create a latent_dataset
#' dset <- latent_dataset(lvec, TR=2, run_length=100)
#' @export
latent_dataset <- function(lvec, TR, run_length, event_table=data.frame()) {
  assert_that(sum(run_length) == dim(lvec)[4])
  
  frame <- sampling_frame(run_length, TR)
  
  ret <- list(
    lvec=lvec,
    datamat=lvec@basis,
    TR=TR,
    nruns=length(run_length),
    event_table=event_table,
    sampling_frame=frame,
    mask=rep(1,ncol(lvec@basis))
  )
  
  class(ret) <- c("latent_dataset", "matrix_dataset", "fmri_dataset", "list")
  ret
  
}



#' Create an fMRI Dataset Object from a Set of Scans
#'
#' This function creates an fMRI dataset object from a set of scans, design information, and other data. The dataset is a list containing information about the scans, mask, TR, number of runs, event table, base path, sampling frame, censor, mode, and preload.
#'
#' @param scans A vector of one or more file names of the images comprising the dataset.
#' @param mask Name of the binary mask file indicating the voxels to include in the analysis.
#' @param TR The repetition time in seconds of the scan-to-scan interval.
#' @param run_length A vector of one or more integers indicating the number of scans in each run.
#' @param event_table A data.frame containing the event onsets and experimental variables. Default is an empty data.frame.
#' @param base_path The file path to be prepended to relative file names. Default is "." (current directory).
#' @param censor A binary vector indicating which scans to remove. Default is NULL.
#' @param preload Read image scans eagerly rather than on first access. Default is FALSE.
#' @param mode The type of storage mode ('normal', 'bigvec', 'mmap', filebacked'). Default is 'normal'.
#'
#' @return An fMRI dataset object of class c("fmri_file_dataset", "volumetric_dataset", "fmri_dataset", "list").
#' @export
#'
#' @examples
#' # Create an fMRI dataset with 3 scans and a mask
#' dset <- fmri_dataset(c("scan1.nii", "scan2.nii", "scan3.nii"), 
#'   mask="mask.nii", TR=2, run_length=rep(300, 3), 
#'   event_table=data.frame(onsets=c(3, 20, 99, 3, 20, 99, 3, 20, 99), 
#'   run=c(1, 1, 1, 2, 2, 2, 3, 3, 3))
#' )
#'
#' # Create an fMRI dataset with 1 scan and a mask
#' dset <- fmri_dataset("scan1.nii", mask="mask.nii", TR=2, 
#'   run_length=300, 
#'   event_table=data.frame(onsets=c(3, 20, 99), run=rep(1, 3))
#' ) 
fmri_dataset <- function(scans, mask, TR, 
                         run_length, 
                         event_table=data.frame(), 
                         base_path=".",
                         censor=NULL,
                         preload=FALSE,
                         mode=c("normal", "bigvec", "mmap", "filebacked")) {
  
  assert_that(is.character(mask), msg="'mask' should be the file name of the binary mask file")
  mode <- match.arg(mode)
  
  #if (length(run_length) == 1) {
  #  run_length <- rep(run_length, length(scans))
  #}
  
  ## run_length should equal total length of images in scans -- but we can 't confirm that here.
  
  if (is.null(censor)) {
    censor <- rep(0, sum(run_length))
  }
  
  frame <- sampling_frame(run_length, TR)
  
  #assert_that(length(run_length) == length(scans))
  
  maskfile <- paste0(base_path, "/", mask)
  scans=paste0(base_path, "/", scans)

  maskvol <- if (preload) {
    assert_that(file.exists(maskfile))
    message(paste("preloading masks", maskfile))
    neuroim2::read_vol(maskfile)
  }
  
  vec <- if (preload) {
    message(paste("preloading scans", paste(scans, collapse = " ")))
    neuroim2::read_vec(scans, mode=mode,mask=maskvol)
  }
  
  
  ret <- list(
    scans=scans,
    vec=vec,
    mask_file=maskfile,
    mask=maskvol,
    nruns=length(run_length),
    event_table=suppressMessages(as_tibble(event_table,.name_repair="check_unique")),
    base_path=base_path,
    sampling_frame=frame,
    censor=censor,
    mode=mode,
    preload=preload
  )
  
  class(ret) <- c("fmri_file_dataset", "volumetric_dataset", "fmri_dataset", "list")
  ret
}



#' @export
#' @importFrom neuroim2 NeuroVecSeq 
get_data.latent_dataset <- function(x, ...) {
  x$lvec@basis
}

#' @export
#' @importFrom neuroim2 NeuroVecSeq 
get_data.fmri_mem_dataset <- function(x, ...) {
  if (length(x$scans) > 1) {
    do.call(neuroim2::NeuroVecSeq, x$scans)
  } else {
    x$scans[[1]]
  }
}

#' @export
#' @importFrom neuroim2 NeuroVecSeq 
get_data.matrix_dataset <- function(x, ...) {
  x$datamat
}

#' @export
#' @importFrom neuroim2 NeuroVecSeq FileBackedNeuroVec
get_data.fmri_file_dataset <- function(x, ...) {
  if (is.null(x$vec)) {
    get_data_from_file(x,...)
  } else {
    x$vec
  }
}

#' @import memoise
#' @keywords internal
get_data_from_file <- memoise::memoise(function(x, ...) {
  m <- get_mask(x)
  neuroim2::read_vec(x$scans, mask=m, mode=x$mode, ...)
})



#' @export
get_mask.fmri_file_dataset <- function(x, ...) {
  if (is.null(x$mask)) {
    neuroim2::read_vol(x$mask_file)
  } else {
    x$mask
  }
}


#' @export
get_mask.fmri_mem_dataset <- function(x, ...) {
  x$mask
}

#' @export
get_mask.matrix_dataset <- function(x, ...) {
  x$mask
}

#' @export
get_mask.latent_dataset <- function(x, ...) {
  x$lvec@mask
}


#' @keywords internal
data_chunk <- function(mat, voxel_ind, row_ind, chunk_num) {
  ret <- list(
       data=mat,
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
  
  iter <- list(nchunks=nchunks, nextElem=nextEl)
  class(iter) <- c("chunkiter", "abstractiter", "iter")
  iter
}


#' Create Data Chunks for fmri_mem_dataset Objects
#'
#' This function creates data chunks for fmri_mem_dataset objects. It allows for the retrieval of run-wise or sequence-wise data chunks, as well as arbitrary chunks.
#'
#' @param x An object of class 'fmri_mem_dataset'.
#' @param nchunks The number of data chunks to create. Default is 1.
#' @param runwise If TRUE, the data chunks are created run-wise. Default is FALSE.
#' @param ... Additional arguments.
#'
#' @return A list of data chunks, with each chunk containing the data, voxel indices, row indices, and chunk number.
#' @importFrom neuroim2 series
#' @autoglobal
#' @export
#'
#' @examples
#' # Create an fmri_mem_dataset
#' # ... (see example for fmri_mem_dataset)
#'
#' # Create an iterator with 5 chunks
#' iter <- data_chunks(dset, nchunks=5)
#' `%do%` <- foreach::`%do%`
#' y <- foreach::foreach(chunk = iter) %do% { colMeans(chunk$data) }
#' length(y) == 5
#'
#' # Create an iterator with 100 chunks
#' iter <- data_chunks(dset, nchunks=100)
#' y <- foreach::foreach(chunk = iter) %do% { colMeans(chunk$data) }
#' length(y) == 100
#'
#' # Create a "runwise" iterator
#' iter <- data_chunks(dset, runwise=TRUE)
#' y <- foreach::foreach(chunk = iter) %do% { colMeans(chunk$data) }
#' length(y) == 2
data_chunks.fmri_mem_dataset <- function(x, nchunks=1,runwise=FALSE,...) {
  mask <- get_mask(x)
  #print("data chunks")
  #print(nchunks)
  get_run_chunk <- function(chunk_num) {
    bvec <- x$scans[[chunk_num]]
    voxel_ind <- which(mask>0)
    #print(voxel_ind)
    row_ind <- which(x$sampling_frame$blockids == chunk_num)
    ret <- data_chunk(neuroim2::series(bvec,voxel_ind), 
                      voxel_ind=voxel_ind, 
                      row_ind=row_ind, 
                      chunk_num=chunk_num)
  }
  
  get_seq_chunk <- function(chunk_num) {
    bvecs <- x$scans
    voxel_ind <- maskSeq[[chunk_num]]
    #print(voxel_ind)

    m <- do.call(rbind, lapply(bvecs, function(bv) neuroim2::series(bv, voxel_ind)))
    ret <- data_chunk(do.call(rbind, lapply(bvecs, function(bv) neuroim2::series(bv, voxel_ind))), 
                      voxel_ind=voxel_ind, 
                      row_ind=1:nrow(m),
                      chunk_num=chunk_num)
    
  }
  
  maskSeq <- NULL
  if (runwise) {
    chunk_iter(x, length(x$scans), get_run_chunk)
  } else if (nchunks == 1) {
    maskSeq <- one_chunk()
    chunk_iter(x, 1, get_seq_chunk)
  #} #else if (nchunks == dim(mask)[3]) {
    #maskSeq <<- slicewise_chunks(x)
    #chunk_iter(x, length(maskSeq), get_seq_chunk)
  } else {
    maskSeq <- arbitrary_chunks(x, nchunks)
    chunk_iter(x, length(maskSeq), get_seq_chunk)
  }
  
}


#' Create Data Chunks for fmri_file_dataset Objects
#'
#' This function creates data chunks for fmri_file_dataset objects. It allows for the retrieval of run-wise or sequence-wise data chunks, as well as arbitrary chunks.
#'
#' @param x An object of class 'fmri_file_dataset'.
#' @param nchunks The number of data chunks to create. Default is 1.
#' @param runwise If TRUE, the data chunks are created run-wise. Default is FALSE.
#' @param ... Additional arguments.
#'
#' @return A list of data chunks, with each chunk containing the data, voxel indices, row indices, and chunk number.
#' @export
#'
#' @examples
#' # Create an fmri_file_dataset
#' # ... (see example for fmri_dataset)
#'
#' # Create an iterator with 5 chunks
#' iter <- data_chunks(dset, nchunks=5)
#' `%do%` <- foreach::`%do%`
#' y <- foreach::foreach(chunk = iter) %do% { colMeans(chunk$data) }
#' length(y) == 5
#'
#' # Create an iterator with 100 chunks
#' iter <- data_chunks(dset, nchunks=100)
#' y <- foreach::foreach(chunk = iter) %do% { colMeans(chunk$data) }
#' length(y) == 100
#'
#' # Create a "runwise" iterator
#' iter <- data_chunks(dset, runwise=TRUE)
#' y <- foreach::foreach(chunk = iter) %do% { colMeans(chunk$data) }
#' length(y) == 2
data_chunks.fmri_file_dataset <- function(x, nchunks=1,runwise=FALSE,...) {
  mask <- get_mask(x)
  #maskSeq <- NULL
  iter <- if (runwise) {
    chunk_iter(x, length(x$scans), get_run_chunk)
  } else if (nchunks == 1) {
    maskSeq <- one_chunk(x)
    chunk_iter(x, 1, get_seq_chunk)
  } else {
    maskSeq <- arbitrary_chunks(x, nchunks)
    #print(maskSeq)
    chunk_iter(x, length(maskSeq), get_seq_chunk)
  }
  
  get_run_chunk <- function(chunk_num) {
    bvec <- neuroim2::read_vec(file.path(x$scans[chunk_num]), mask=mask)
    ret <- data_chunk(bvec@data, voxel_ind=which(x$mask>0), 
                      row_ind=which(x$sampling_frame$blockids == chunk_num), 
                      chunk_num=chunk_num)
  }
  
  get_seq_chunk <- function(chunk_num) {
  
    #v <- x$vec
    v <- get_data(x)
    #bvecs <- lapply(x$scans, function(scan) neuroim2::read_vec(scan, mask=maskSeq[[chunk_num]]))
    vind=maskSeq[[chunk_num]]
    m <- series(v, vind)
    ret <- data_chunk(m, voxel_ind=vind, 
                      row_ind=1:nrow(x$event_table), 
                      chunk_num=chunk_num)
    
  }
  
  iter
  
  ##message("nchunks is ", nchunks)
  

  
}


#' Create Data Chunks for matrix_dataset Objects
#'
#' This function creates data chunks for matrix_dataset objects. It allows for the retrieval of run-wise or sequence-wise data chunks, as well as arbitrary chunks.
#'
#' @param x An object of class 'matrix_dataset'.
#' @param nchunks The number of data chunks to create. Default is 1.
#' @param runwise If TRUE, the data chunks are created run-wise. Default is FALSE.
#' @param ... Additional arguments.
#'
#' @return A list of data chunks, with each chunk containing the data, voxel indices, row indices, and chunk number.
#' @export
#'
#' @examples
#' # Create a matrix_dataset
#' # ... (see example for matrix_dataset)
#'
#' # Create an iterator with 5 chunks
#' iter <- data_chunks(dset, nchunks=5)
#' `%do%` <- foreach::`%do%`
#' y <- foreach::foreach(chunk = iter) %do% { colMeans(chunk$data) }
#' length(y) == 5
#'
#' # Create an iterator with 100 chunks
#' iter <- data_chunks(dset, nchunks=100)
#' y <- foreach::foreach(chunk = iter) %do% { colMeans(chunk$data) }
#' length(y) == 100
#'
#' # Create a "runwise" iterator
#' iter <- data_chunks(dset, runwise=TRUE)
#' y <- foreach::foreach(chunk = iter) %do% { colMeans(chunk$data) }
#' length(y) == 2
data_chunks.matrix_dataset <- function(x, nchunks=1, runwise=FALSE,...) {
  get_run_chunk <- function(chunk_num) {
    ind <- which(blockids(x$sampling_frame) == chunk_num)
    mat <- x$datamat[ind,,drop=FALSE]
    #browser()
    data_chunk(mat, voxel_ind=1:ncol(mat), row_ind=ind, chunk_num=chunk_num)
  }
  
  get_one_chunk <- function(chunk_num) {
    data_chunk(x$datamat, voxel_ind=1:ncol(x$datamat), row_ind=1:nrow(x$datamat), chunk_num=chunk_num)
  }
  
  
  if (runwise) {
    chunk_iter(x, length(x$sampling_frame$blocklens), get_run_chunk)
  } else if (nchunks==1) {
    chunk_iter(x, 1, get_one_chunk)
  } else {
    sidx <- split(1:ncol(x$datamat), sort(rep(1:nchunks, length.out=ncol(x$datamat))))
    get_chunk <- function(chunk_num) {
      data_chunk(x$datamat[,sidx[[chunk_num]], drop=FALSE], 
                 voxel_ind=sidx[[chunk_num]], 
                 row_ind=1:nrow(x$datamat), 
                 chunk_num=chunk_num)
    }
    chunk_iter(x, nchunks, get_chunk)
  }
  
}

#' @keywords internal
exec_strategy <- function(strategy=c("voxelwise", "runwise", "chunkwise"), nchunks=NULL) {
  strategy <- match.arg(strategy)
  
  function(dset) {
    if (strategy == "runwise") {
      data_chunks(dset, runwise=TRUE)
    } else if (strategy == "voxelwise") {
      m <- get_mask(dset)
      data_chunks(dset, nchunks = sum(m), runwise=FALSE)
    } else if (strategy == "chunkwise") {
      m <- get_mask(dset)
      ##message("nchunks is", nchunks)
      assert_that(!is.null(nchunks) && is.numeric(nchunks))
      if (nchunks > sum(m)) {
        warning("requested number of chunks is greater than number of voxels in mask")
        nchunks <- sum(m)
      }
      data_chunks(dset, nchunks = nchunks, runwise=FALSE)
    }
  }
  
}


#' @keywords internal
arbitrary_chunks <- function(x, nchunks) {
  #print("arbitrary chunks")
  #browser()
  mask <- get_mask(x)
  #print(mask)
  indices <- as.integer(which(mask != 0))
  chsize <- round(length(indices)/nchunks)
  #print(indices)
  
  assert_that(chsize > 0)
  chunkids <- sort(rep(1:nchunks, each=chsize, length.out=length(indices)))
  #print(chunkids)
  
  mfun <- function(i) indices[chunkids==i]
  #print(mfun)
  ret <- neuroim2::deferred_list2(mfun, nchunks)
  #print(ret[[1]])
  return(ret)
  
}

#' @keywords internal
slicewise_chunks <- function(x) {
  mask <- x$mask
  template <- neuroim2::NeuroVol(array(0, dim(mask)), neuroim2::space(mask))
  nchunks <- dim(mask)[3]
  
  maskSeq <- lapply(1:nchunks, function(i) {
    m <- template
    m[,,i] <- 1
    m
  })
  
  maskSeq
  
}

#' @keywords internal
one_chunk <- function(x) {
  mask <- get_mask(x)
  list(mask)
}

#' @export
print.fmri_dataset <- function(x, ...) {
  cat("fmri_dataset", "\n")
  cat("  number of runs: ", x$nruns, "\n")
  print(x$sampling_frame)
  cat("  event_table: ", "\n")
  print(x$event_table)
}


#' @export
print.matrix_dataset <- function(x,...) {
  cat("matrix_dataset", "\n")
  cat("  number of runs: ", x$nruns, "\n")
  cat("  number of rows: ", nrow(x$datamat), "\n")
  cat("  number of columns: ", ncol(x$datamat), "\n")
  print(x$sampling_frame)
  cat("  event_table: ", "\n")
  print(x$event_table)
}

#' @export
print.latent_dataset <- function(x,...) {
  cat("latent_dataset", "\n")
  cat("  number of runs: ", x$nruns, "\n")
  cat("  number of rows: ", nrow(x$datamat), "\n")
  cat("  number of latent variables: ", ncol(x$datamat), "\n")
  print(x$sampling_frame)
  cat("  event_table: ", "\n")
  print(x$event_table)
}

#' @export
print.chunkiter <- function(x, ...) {
  cat(paste("chunk iterator with", x$nchunks, " chunks"))
}






