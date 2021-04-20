

read_confounds <- function(x,...) UseMethod("read_confounds")

bids_source <- function(bids_path, deriv_folder="derivatives/fmriprep", id, 
                        bold_space, task=NULL, session=NULL, confound_vars="FramewiseDisplacement") {

  scan_map <- paste0(bids_path, "/", "sub-", id, "/sub-", id, "_scans.tsv")
  
  
  if (!file.exists(scan_map)) {
    stop("cannot find scan_map named: ", paste0(bids_path, "/", "sub-", id, "/sub-", id, "_scans.tsv"))
  }
  
  scan_map <- read.table(scan_map, header=TRUE, stringsAsFactors=FALSE, colClasses=list(run="character", scan="character", session="character"))
  
  if (!is.null(session) && session != "") {
    qsession <- dplyr::enquo(session)
    scan_map <- dplyr::filter(scan_map, session == !!qsession)
  }
  
  if (!is.null(task)) {
    qtask <- dplyr::enquo(task)
    scan_map <- scan_map %>% dplyr::filter(task == !!qtask)
  } 
  
  zerostr <- function(vals, ndigits=2) {
    #ndigits <- log(max(vals), 10) + 1
    nzeros <- ndigits - as.integer(log(vals,10)) -1
    prefix <- sapply(nzeros, function(nz) paste(rep("0", times=nz), collapse =""))
    paste(prefix, vals, sep="")  
  }
  
  snum <- zerostr(as.numeric(as.character(scan_map$scan)))
  
  sess <- function(isdir=TRUE) {
    if (is.null(session) || session == "") {
      ""
    } else {
      if (isdir) paste0("/ses-", session, "/") else paste0("_ses-", session) 
    }
  }
  
  confounds <- paste0(bids_path, "/", deriv_folder, "/sub-", id, sess(), "/func/", paste0("sub-", id, sess(FALSE), "_task-", task, "_run-", snum, "_bold_confounds.tsv"))

  raw_scans <- scan_map$func
  event_files <- paste0(bids_path, "/sub-", id, sess(), "/func/", gsub("bold.nii.gz", "events.tsv", basename(raw_scans)))
  
  scans <- paste0("sub-", id, sess(FALSE), "_task-", task, "_run-", snum, "_bold_space-", bold_space, ".nii.gz")
  
  structure(
    list(bids_path=bids_path,
         deriv_folder=deriv_folder,
         func_path=paste0(bids_path, "/", deriv_folder, "/sub-", id, sess(), "/func"),
         scan_map=scan_map,
         raw_scans=raw_scans,
         preproc_scans=scans,
         event_files=event_files,
         confound_files=confounds,
         confound_vars=confound_vars,
         task=task),
    class="bids_source")
}

read_confounds.bids_source <- function(x, replace_na=c("median", "mean"), scale=TRUE, center=TRUE, cvars=NULL, npcs=-1, perc_var=-1) {
  replace_na <- match.arg(replace_na)
  
  cfiles <-x$confound_files
  cvars <- if (is.null(cvars)) x$confound_vars else cvars
  
  lapply(cfiles, function(cf) {
    m <- as.matrix(read.table(cf, header=TRUE, na.strings="n/a") %>% dplyr::select(cvars))
    
    
    if (any(is.na(m))) {
      f <- get(replace_na)
      m <- apply(m, 2, function(v) {
        mu <- f(v, na.rm=TRUE)
        v[is.na(v)] <- mu
        v
      })
    }
    
    sm <- scale(m, center=center,scale=scale)
    
    if (npcs > 0 || perc_var > 0 && (ncol(sm) > 1)) {
      pres <- prcomp(sm, scale=TRUE)
      varexp <- cumsum(pres$sdev^2)/sum(pres$sdev^2) * 100
      if (npcs > 0 && perc_var <= 0) {
        sm <- pres$x[,1:npcs,drop=FALSE]
      } else if (npcs <= 0 && perc_var > 0) {
        keep <- which( (varexp - perc_var) > 0)[1]
        sm <- pres$x[, 1:keep,drop=FALSE]
      } else {
        keep <- which( (varexp - perc_var) > 0)[1]
        npcs <- max(c(keep, npcs))
        if (npcs < 1) {
          npcs <- 1
        }
        sm <- pres$x[, 1:npcs, drop=FALSE]
      }
    }
      
    sm
      
  })
}


 