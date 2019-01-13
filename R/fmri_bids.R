
get_prepped_scans <- function(x, id, task, space,abspath) UseMethod("get_prepped_scans")


bids_source <- function(bids_path, deriv_folder="derivatives/fmriprep", id, bold_space, task=NULL, confound_vars=NULL) {
  scan_map <- paste0(bids_path, "/", "sub-", id, "/sub-", id, "_scans.tsv")
  
  
  if (!file.exists(scan_map)) {
    stop("cannot find scan_map named: ", paste0(bids_path, "/", "sub-", id, "/sub-", id, "_scans.tsv"))
  }
  
  scan_map <- read.table(scan_map, header=TRUE, stringsAsFactors=FALSE)
  
  if (!is.null(task)) {
    qtask <- dplyr::quo(task)
    scan_map <- scan_map %>% filter(task == !!qtask)
  } 
  
  snum <- as.character(scan_map$scan)
  confounds <- paste0(bids_path, "/", deriv_folder, "/", paste0("sub-", id, "_task-", task, "_run-", snum, "_bold_confounds.tsv"))

  raw_scans <- scan_map$func
  event_files <- paste0(bids_path, "/sub-", id, "/func/", gsub("bold.nii.gz", "events.tsv", raw_scans))
  
  scans <- paste0("sub-", id, "_task-", task, "_run-", snum, "_bold_space-", bold_space, "_preproc.nii.gz")
  
  structure(
    list(bids_path=bids_path,
         deriv_folder=deriv_folder,
         func_path=paste0(bids_path, "/", deriv_folder, "/func"),
         scan_map=scan_map,
         scans=scans,
         event_files=event_files,
         confounds=confounds),
    class="bids_source")
}

get_scan_map <- function(x, id, task=NULL) {
  #scan_map <- read.table(paste0(x$bids_path, "/", "sub-", id, "/sub-", id, "_scans.tsv"), header=TRUE)
  
  #if (!file.exists(scan_map)) {
  #  stop("cannot find scan_map named: ", paste0(x$bids_path, "/", "sub-", id, "/sub-", id, "_scans.tsv"))
  #}
  
  if (!is.null(task)) {
    qtask <- dplyr::quo(task)
    x$scan_map %>% filter(task == !!qtask)
  } else {
    x$scan_map
  }
}

get_confound_files <- function(x, id, task, abspath=TRUE) {
  scan_map <- get_scan_map(x,id, task)
  snums <- scan_map$scan
  
  paste0(x$bids_path, "/", x$deriv_folder, "/", paste0("sub-", id, "_task-", task, "_run-", snums, "_bold_confounds.tsv"))
  
}

load_confound_files <- function(x, id, task, abspath=TRUE) {
  cfiles <- get_confound_files(x,id,task,abspath)
  lapply(cfiles, read.table, header=TRUE, na.strings="n/a")
}

get_prepped_scans.bids_query <- function(x, id, task, space,abspath=TRUE) {
  scan_map <- get_scan_map(x,id, task)
  
  snums <- scan_map$scan
  scans <- paste0("sub-", id, "_task-", task, "_run-", snums, "_bold_space-", space, "_preproc.nii.gz")
  
  if (abspath) {
    p <- paste0(x$bids_path, "/", x$deriv_folder, "/", "sub-", id, "/func/")
    scans <- paste0(p, scans)
  }
  
  scans
}

get_event_files <- function(x, id, task, space,abspath=TRUE) {
  
}
  
  