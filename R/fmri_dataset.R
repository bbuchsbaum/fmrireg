
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
  print(file_name)
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
  env$design <- tibble::as_tibble(read.table(dname, header=TRUE))

  out <- as.list(env)
  class(out) <- c("fmri_config", "list")
  out
}

#' matrix_dataset
#' 
#' @inheritParams fmri_dataset
#' @param datamat a \code{matrix} each column is a voxel time-series
#' @export
#' @examples 
#' 
#' ## a matrix with 100 rows and 100 columns (voxels)
#' X <- matrix(rnorm(100*100), 100, 100)
#' dset <- matrix_dataset(X, TR=2, run_length=100)
#' 
#' ## an iterator with 5 chunks
#' iter <- data_chunks(dset, nchunks=5)
#' y <- foreach(chunk = iter) %do% { colMeans(chunk$data) }
#' length(y) == 5
#' 
#' ## an iterator with 100 chunks
#' iter <- data_chunks(dset, nchunks=100)
#' y <- foreach(chunk = iter) %do% { colMeans(chunk$data) }
#' length(y) == 100
#' 
#' ## a matrix_dataset with 200 rows, 100 columns and 2 runs
#' X <- matrix(rnorm(200*100), 200, 100)
#' dset <- matrix_dataset(X, TR=2, run_length=c(100,100))
#' 
#' ## get a "runwise" iterator. For every iteration an entire run's worth of data is returned.
#' iter <- data_chunks(dset, runwise=TRUE)
#' y <- foreach(chunk = iter) %do% { colMeans(chunk$data) }
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

#' fmri_mem_dataset
#' 
#' @inheritParams fmri_dataset
#' @param scans a \code{list} of objects of class \code{\linkS4class{NeuroVec}}
#' @param mask a binary mask of class  \code{\linkS4class{NeuroVol}} indicating the set of voxels to include in analyses
#' @export
#' @examples 
#' 
#' d <- c(10,10,10,10)
#' nvec <- neuroim2::NeuroVec(array(rnorm(prod(d)), d), space=neuroim2::NeuroSpace(d))
#' mask <- neuroim2::NeuroVol(array(rnorm(10*10*10), d[1:3]), space=neuroim2::NeuroSpace(d[1:3]))
#' mask[mask < .5] <- 0
#' dset <- fmri_mem_dataset(list(nvec),mask, TR=2)
#' 
#' iter <- data_chunks(dset, nchunks=100)
fmri_mem_dataset <- function(scans, mask, TR, 
                             event_table=data.frame(), 
                             base_path=".",
                             censor=NULL) {
  
  
  
  assert_that(all(map_lgl(scans, function(x) inherits(x, "NeuroVec"))))
  assert_that(inherits(mask, "NeuroVol"))
  assert_that(all(dim(mask) == dim(scans[[1]][1:3])))
  
 
  run_length <- map_dbl(scans, ~ dim(.)[4])
  assert_that(sum(run_length) > length(scans))
  
  if (is.null(censor)) {
    censor <- rep(0, sum(run_length))
  }

  frame <- sampling_frame(run_length, TR)
  
  ret <- list(
    
    scans=scans,
    mask=mask,
    nruns=length(scans),
    event_table=event_table,
    base_path=base_path,
    sampling_frame=frame,
    censor=censor
  )
  
  class(ret) <- c("fmri_mem_dataset", "volumetric_dataset", "fmri_dataset", "list")
  ret
}

#' A dataset that encapsulates a dimension-reduced subspace of "latent variables". 
#' 
#' 
#' @inheritParams fmri_dataset
#' @param lvec an instance of class \code{LatentNeuroVec}
#' @examples 
#' 
#' a matrix with 100 rows and 1000 columns (voxels)
#' X <- matrix(rnorm(100*1000), 100, 1000)
#' pres <- prcomp(X)
#' basis <- pres$x[,1:25]
#' loadings <- pres$rotation[,1:25]
#' offset <- colMeans(X)
#' lvec <- LatentNeuroVec(basis, loadings, NeuroSpace(c(10,10,10,100)), mask=rep(TRUE,1000), offset=offset)
#' dset <- latent_dataset(lvec, TR=2, run_length=100)
#' 
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



#' An fMRI dataset consisting of a set of scans as files, design information, and other data.
#' 
#' @param scans a vector of one or more file names of the images comprising the dataset
#' @param mask name of the binary mask file indicating the voxels to include in analysis.
#' @param TR the repetition time in seconds of the scan-to-scan interval.
#' @param run_length a \code{vector} of one or more integers indicating the number of scans in each run.
#' @param event_table a \code{data.frame} containing the event onsets and experimental variables.
#' @param base_path the file path to be prepended to relative file names.
#' @param censor a binary vector indicating which scans to remove.
#' @export
#' @importFrom tibble as_tibble
#' @examples 
#' 
#' dset <- fmri_dataset(c("scan1.nii", "scan2.nii", "scan3.nii"), mask="mask.nii", TR=2, run_length=rep(300,3), 
#'         event_table=data.frame(onsets=c(3,20,99,3,20,99,3,20,99), run=c(1,1,1,2,2,2,3,3,3)))
#'         
#' dset <- fmri_dataset("scan1.nii", mask="mask.nii", TR=2, run_length=300, 
#'         event_table=data.frame(onsets=c(3,20,99), run=rep(1,3)))
#' 
#' 
fmri_dataset <- function(scans, mask, TR, 
                         run_length, 
                         event_table=data.frame(), 
                         base_path=".",
                         censor=NULL,
                         preload=FALSE,
                         mode=c("bigvec", "mmap", "filebacked")) {
  
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
  assert_that(file.exists(maskfile))
  maskvol <- neuroim2::read_vol(maskfile)
  
  scans=paste0(base_path, "/", scans)

  vec <- if (preload) {
    message(paste("preloading scans", paste(scans, collapse = " ")))
    read_vec(scans, mode=mode,mask=maskvol)
  }
  
  ret <- list(
    scans=paste0(base_path, "/", scans),
    vec=vec,
    mask_file=maskfile,
    mask=maskvol,
    nruns=length(scans),
    event_table=as_tibble(event_table),
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
  do.call(neuroim2::NeuroVecSeq, x$scans)
}

#' @export
#' @importFrom neuroim2 NeuroVecSeq 
get_data.matrix_dataset <- function(x, ...) {
  x$datamat
}

#' @export
#' @importFrom neuroim2 NeuroVecSeq FileBackedNeuroVec
get_data.fmri_file_dataset <- function(x, ...) {
  ## memoise?
  if (x$preload) {
    x$vec
  } else {
    read_vec(scans, mode=mode,mask=x$maskvol)
  }
}

#' @export
get_mask.fmri_file_dataset <- function(x) {
  x$mask
}


#' @export
get_mask.fmri_mem_dataset <- function(x) {
  x$mask
}

#' @export
get_mask.matrix_dataset <- function(x) {
  x$mask
}

#' @export
get_mask.latent_dataset <- function(x) {
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


#' @importFrom neuroim2 series
#' @export
data_chunks.fmri_mem_dataset <- function(x, nchunks=1,runwise=FALSE) {
  
  mask <- get_mask(x)
  
  get_run_chunk <- function(chunk_num) {
    bvec <- x$scans[[chunk_num]]
    voxel_ind <- which(mask>0)
    row_ind <- which(x$sampling_frame$blockids == chunk_num)
    ret <- data_chunk(neuroim2::series(bvec,voxel_ind), 
                      voxel_ind=voxel_ind, 
                      row_ind=row_ind, 
                      chunk_num=chunk_num)
  }
  
  get_seq_chunk <- function(chunk_num) {
    bvecs <- x$scans
    voxel_ind <- maskSeq[[chunk_num]]
    m <- do.call(rbind, lapply(bvecs, function(bv) neuroim2::series(bv, voxel_ind)))
    ret <- data_chunk(do.call(rbind, lapply(bvecs, function(bv) neuroim2::series(bv, voxel_ind))), 
                      voxel_ind=voxel_ind, 
                      row_ind=1:nrow(m),
                      chunk_num=chunk_num)
    
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



#' @import neuroim2
#' @export
data_chunks.fmri_file_dataset <- function(x, nchunks=1,runwise=FALSE) {
  
  mask <- get_mask(x)
  
  iter <- if (runwise) {
    chunk_iter(x, length(x$scans), get_run_chunk)
  } else if (nchunks == 1) {
    maskSeq <<- one_chunk()
    chunk_iter(x, 1, get_seq_chunk)
  } else {
    maskSeq <<- arbitrary_chunks(x, nchunks)
    chunk_iter(x, length(maskSeq), get_seq_chunk)
  }
  
  get_run_chunk <- function(chunk_num) {
    bvec <- neuroim2::read_vec(file.path(x$scans[chunk_num]), mask=mask)
    ret <- data_chunk(bvec@data, voxel_ind=which(x$mask>0), 
                      row_ind=which(x$sampling_frame$blockids == chunk_num), 
                      chunk_num=chunk_num)
  }
  
  get_seq_chunk <- function(chunk_num) {
  
    v <- x$vec
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

#' @import neuroim2
#' @export
data_chunks.matrix_dataset <- function(x, nchunks=1, runwise=FALSE) {
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
      data_chunk(x$datamat[,sidx[[chunk_num]], drop=FALSE], voxel_ind=sidx[[chunk_num]], row_ind=1:nrow(x$datamat), chunk_num=chunk_num)
    }
    chunk_iter(x, nchunks, get_chunk)
  }
  
}

#' @export
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
  mask <- get_mask(x)
  indices <- as.integer(which(mask != 0))
  chsize <- round(length(indices)/nchunks)
  assert_that(chsize > 0)
  chunkids <- sort(rep(1:nchunks, each=chsize, length.out=length(indices)))
  
  mfun <- function(i) indices[chunkids==i]
  neuroim2::deferred_list2(mfun, nchunks)
  
}

#' @keywords internal
slicewise_chunks <- function(x) {
  mask <- x$mask
  template <- NeuroVol(array(0, dim(mask)), space(mask))
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
  mask <- x$mask
  list(mask)
}

#' @export
print.fmri_dataset <- function(object) {
  cat("fmri_dataset", "\n")
  cat("  number of runs: ", object$nruns, "\n")
  print(object$sampling_frame)
  cat("  event_table: ", "\n")
  print(object$event_table)
}


#' @export
print.matrix_dataset <- function(object) {
  cat("matrix_dataset", "\n")
  cat("  number of runs: ", object$nruns, "\n")
  cat("  number of rows: ", nrow(object$datamat), "\n")
  cat("  number of columns: ", ncol(object$datamat), "\n")
  print(object$sampling_frame)
  cat("  event_table: ", "\n")
  print(object$event_table)
}

#' @export
print.latent_dataset <- function(object) {
  cat("latent_dataset", "\n")
  cat("  number of runs: ", object$nruns, "\n")
  cat("  number of rows: ", nrow(object$datamat), "\n")
  cat("  number of latent variables: ", ncol(object$datamat), "\n")
  print(object$sampling_frame)
  cat("  event_table: ", "\n")
  print(object$event_table)
}

print.chunkiter <- function(object) {
  cat(paste("chunk iterator with", object$nchunks, " chunks"))
}






