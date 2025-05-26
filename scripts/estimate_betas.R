#!/usr/bin/env Rscript

###############################################################################
# estimate_betas_improved.R
#
# A modern R script that demonstrates single-trial beta estimation using
# the 'fmrireg' library. Replaces older AFNI-based scripts.
#
# Usage example:
#   ./estimate_betas_improved.R --bids_path /path/to/BIDS \
#       --subid 01 --task stroop --mask "sub-01_task-stroop_space-MNI152_desc-mask.nii.gz"
#
# See '--help' for details on arguments.
###############################################################################

library(optparse)
suppressPackageStartupMessages({
  library(fmrireg)
  library(bidser)
  library(tidyverse)
})

log_info <- function(msg) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), msg, "\n")
}

################################################################################
## 1) Define Command-Line Options
################################################################################

option_list <- list(
  make_option(c("-b", "--bids_path"), type="character", default=".",
              help="Path to the BIDS directory [default = %default]."),
  make_option(c("--bids_session"), type="character", default="",
              help="BIDS session label (optional)."),
  make_option(c("-d", "--deriv_folder"), type="character", default="derivatives/fmriprep",
              help="Relative path under BIDS containing preprocessed scans [default = %default]."),
  make_option(c("--outdir"), type="character", default="betas",
              help="Output directory name [default = %default]."),
  make_option(c("-t", "--task"), type="character",
              help="Name of the BIDS task (required)."),
  make_option(c("-s", "--subid"), type="character",
              help="Subject ID (required). Example: '01'."),
  make_option(c("-x", "--space"), type="character", default="MNI152NLin2009cAsym_preproc",
              help="BOLD space name in fMRIPrep filenames [default = %default]."),

  make_option(c("-m", "--mask"), type="character", default=NULL,
              help="Name (or relative path) of the binary mask image (NIfTI). Optional if a mask is auto-generated."),

  make_option(c("-l", "--duration"), type="numeric", default=1,
              help="Duration of each stimulus event (constant). If your event file has a 'duration' column, 
                    this script can ignore or override it. [default=%default]"),

  make_option(c("--milliseconds"), action="store_true", default=FALSE,
              help="If TRUE, treat event onsets as milliseconds and divide by 1000 [default=%default]."),

  make_option(c("--onsets"), type="character", default="onset",
              help="Name of the event onsets column in the .tsv events file [default=%default]."),

  make_option(c("-c", "--confounds"), type="character", default=NULL,
              help="Path to a text file listing confound variable names to read from confounds.tsv (one per line)."),

  make_option(c("-v", "--percvar"), type="numeric", default=95,
              help="Percent of confound variance to retain (for confound PCA) [default=%default]."),

  make_option(c("--method"), type="character", default="lss",
              help="Single-trial beta estimation method, e.g., 'lss', 'ols', 'mixed', etc. [default=%default]."),

  make_option(c("-p", "--polort"), type="integer", default=3,
              help="Number of polynomial regressors for baseline modeling (like 'polort' in AFNI). 
                    3 means up to cubic drift. [default=%default]."),

  make_option(c("-o", "--outstem"), type="character", default="betas",
              help="Base name for output files (NIfTI) [default=%default]."),

  make_option(c("--concatenate"), action="store_true", default=FALSE,
              help="If TRUE, concatenate all run-level betas into one file [default=%default].")
)

parser  <- OptionParser(usage="%prog [options]", option_list=option_list)
opt_par <- parse_args(parser)
args    <- opt_par

################################################################################
## 2) Prepare BIDS and fMRIPrep Data
################################################################################

# BIDS project object from 'bidser'. This will let us easily fetch scans, confounds, events, etc.
session_label <- if (nzchar(args$bids_session)) args$bids_session else NULL

bidsproj <- bidser::bids_project(bids_path = args$bids_path,
                                 fmriprep   = TRUE)  # or set if you have multiple derivatives

sub_id   <- paste0("sub-", args$subid)
taskname <- args$task
if (is.null(taskname) || is.null(args$subid)) {
  stop("ERROR: Must specify --subid and --task at minimum.")
}

# 2A) Get confounds if we have a confounds text file:
confound_vars <- if (!is.null(args$confounds)) {
  conf_vars <- scan(args$confounds, what="", quiet=TRUE)
  message("Using confound variables: ", paste(conf_vars, collapse=", "))
  conf_vars
} else {
  NULL
}

# If we have confound_vars, gather them via `bidser:::read_confounds()`.
nuisanceCovars <- NULL
if (!is.null(confound_vars)) {
  cdat <- bidser:::read_confounds(
    bidsproj,
    perc_var = args$percvar,
    cvars    = confound_vars,
    subid    = args$subid,
    task     = args$task,
    nest     = TRUE  # returns a nested list
  )
  # cdat$data is a list of run-wise confound matrices
  # we store them in nuisanceCovars
  nuisanceCovars <- lapply(cdat$data, function(x) {
    # remove columns that are all NA
    keep <- !apply(x, 2, function(z) all(is.na(z)))
    x[, keep, drop=FALSE]
  })
}

# 2B) Identify event files (the .tsv or .csv with columns onsets, durations, etc.)
event_files <- bidser:::search_files(
  bidsproj,
  subid = args$subid,
  task  = taskname,
  kind  = "events",
  full_path = TRUE
)
if (length(event_files) < 1) {
  stop("No event files found for sub-", args$subid, " task-", taskname)
}

if (length(event_files) != length(scanpaths)) {
  warning("Number of event files (", length(event_files), 
          ") doesn't match number of scans (", length(scanpaths), ").")
  # Don't stop, but log the warning
}

# 2C) Identify preprocessed BOLD scans
scans <- bidser:::preproc_scans(
  bidsproj,
  subid  = args$subid,
  task   = taskname,
  session= session_label,
  space  = args$space
)
scanpaths <- file.path(args$bids_path, scans)  # full path
if (!all(file.exists(scanpaths))) {
  stop("Some bold scans don't exist. Missing: ", paste(scanpaths[!file.exists(scanpaths)], collapse=", "))
}

################################################################################
## 3) Build Output Directory and Mask
################################################################################

outdir <- file.path(dirname(dirname(scanpaths[1])), args$outdir)
dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
message("Output directory: ", outdir)

maskpath <- NULL
if (!is.null(args$mask)) {
  # user-provided relative or absolute path
  maskpath <- if (file.exists(args$mask)) {
    normalizePath(args$mask)
  } else {
    # might be under outdir or something else
    file.path(dirname(dirname(scanpaths[1])), args$mask)
  }
  if (!file.exists(maskpath)) {
    stop("Could not locate mask file: ", maskpath)
  }
  message("Using mask: ", maskpath)
}

################################################################################
## 4) Loop over runs: create an fmri_dataset, read event table, build design
################################################################################

all_beta_files   <- c()  # We'll collect run-level output
method_to_use    <- args$method
polort           <- args$polort
onset_col        <- args$onsets
ms_mode          <- args$milliseconds
fixed_duration   <- as.numeric(args$duration)

# Because fmrireg expects "run_length" per run, we need to see how many volumes each scan has:
library(neuroim2)  # for read_nifti_header
run_lengths <- numeric(length(scanpaths))
for (i in seq_along(scanpaths)) {
  hdr <- neuroim2::read_nifti_header(scanpaths[i])
  run_lengths[i] <- hdr@dim_[5]  # 4th dimension = # time points
}
TRvals <- neuroim2::voxdim(hdr)[4]

if (is.na(TRvals) || TRvals <= 0) {
  stop("Invalid TR value: ", TRvals, ". Please check input scan header.")
}

# For each run, read the events, build an fmri_dataset, then estimate betas
for (i in seq_along(scanpaths)) {
  cat("\n========== RUN #", i, "==========\n")
  cat("  BOLD file:", scanpaths[i], "\n")

  # Attempt to read the event file for this run. 
  # If your naming scheme is sub-01_task-xxx_run-1_events, then we might pick the ith file
  # or find matching. For simplicity, assume 1:1 ordering:
  if (length(event_files) < i) {
    stop("Not enough event files to match the # of scans.")
  }
  evfile  <- event_files[i]
  design  <- read.table(evfile, header=TRUE)

  # Convert onsets from ms to s if needed
  if (!onset_col %in% names(design)) {
    stop("Event file must have a column named '", onset_col, "'")
  }
  if (ms_mode) {
    design[[onset_col]] <- design[[onset_col]] / 1000
  }

  # If we do not trust "duration" column in the events or we want to override it:
  if (!("duration" %in% names(design))) {
    # create one
    design$duration <- fixed_duration
  } else {
    # override existing
    design$duration <- fixed_duration
  }

  # Optional: If we have nuisance confounds for this run
  nuisance_mat <- if (!is.null(nuisanceCovars)) {
    nuisanceCovars[[i]]
  } else {
    NULL
  }

  # Build the fmri_dataset
  ds <- fmri_dataset(
    scans = scanpaths[i],
    mask  = if (!is.null(maskpath)) maskpath else NULL,
    TR    = TRvals,
    run_length = run_lengths[i],
    event_table= design
  )

  # Baseline polynomials: we can easily build a baseline model with polynomials
  # fmrireg includes 'baseline_model("constant", sframe=ds$sampling_frame)' for a constant,
  # or 'baseline_model("poly", degree=polort, sframe=ds$sampling_frame)' for polynomials.

  bmod <- baseline_model("poly", degree=polort, sframe=ds$sampling_frame)

  # Next define single-trial formula: onset ~ trialwise() 
  # That gives 1 regressor per event. 
  # Or if your design is more complicated, define a custom formula
  # for the random part:
  st_formula <- onset ~ trialwise()

  # We do not define 'fixed' for typical single-trial. 
  # So we pass fixed=NULL, block= ~1 (if single run) 
  # or block=~ run if multi-run.

  # Estimate betas
  # The 'method' might be "lss", "ols", "mixed", "pls", etc.
  # example: method="lss"
  bet_result <- estimate_betas(
    ds,
    fixed = NULL,
    ran   = st_formula,
    block = ~1, 
    method= method_to_use,
    basemod = bmod
  )

  # bet_result$betas_ran is dimension (#events x #voxels) 
  # if 'matrix_dataset' or a NeuroVec if fmri_file_dataset

  # For a file-based dataset, we can store the betas as a NIfTI 
  #  (or for matrix_dataset, we can create an image from it).
  # Let's see if ds is fmri_file_dataset => bet_result$betas_ran is a NeuroVec
  # or if ds is matrix_dataset => it's a 2D matrix.

  # We'll write out the random betas as NIfTI in the outdir
  outfile_nii <- file.path(
    outdir,
    paste0(args$outstem, "_run-", i, "_", method_to_use, "_betas.nii.gz")
  )

  if (inherits(bet_result$betas_ran, "NeuroVec")) {
    # direct write
    neuroim2::write_nifti(bet_result$betas_ran, outfile_nii)
  } else {
    # It's a matrix => we need to create a "dummy" NeuroVec
    refmask <- get_mask(ds)  # numeric vector or NeuroVol
    if (is.vector(refmask)) {
      log_info("WARNING: Matrix dataset without proper mask information.")
      # Create a simple output file with matrix data
      # Save as a CSV instead
      outfile_csv <- gsub("\\.nii\\.gz$", ".csv", outfile_nii)
      log_info(paste("Saving matrix betas to CSV instead:", outfile_csv))
      write.csv(bet_result$betas_ran, outfile_csv)
    } else {
      # refmask is a NeuroVol => we can properly convert to NIfTI
      sp <- neuroim2::space(refmask)
      outdim <- c(dim(sp)[1:3], nrow(bet_result$betas_ran))
      arr <- array(0, outdim)
      voxel_inds <- which(refmask > 0)
      for (ev in seq_len(nrow(bet_result$betas_ran))) {
        arr_1d <- bet_result$betas_ran[ev, ]
        arr[,, , ev][voxel_inds] <- arr_1d
      }
      outimg <- neuroim2::NeuroVol4D(arr, sp)
      neuroim2::write_nifti(outimg, outfile_nii)
    }
  }

  message("Wrote betas for run ", i, " => ", outfile_nii)
  all_beta_files <- c(all_beta_files, outfile_nii)

  if (!all(c(onset_col) %in% names(design))) {
    stop("Required column '", onset_col, "' not found in event file: ", evfile)
  }
}

################################################################################
## 5) (Optional) Concatenate betas across runs
################################################################################

if (args$concatenate && length(all_beta_files) > 1) {
  cat("\n--- Concatenating run-level betas into single 4D file ---\n")
  final_4d   <- file.path(outdir, paste0(args$outstem, "_", method_to_use, "_allruns.nii.gz"))
  cat("Output => ", final_4d, "\n")
  # We can do this with 3dTcat or with neuroim2 (if same spatial dims)
  # Here we show a quick approach in R using neuroim2:
  vols <- lapply(all_beta_files, function(f) neuroim2::read_vol4D(f))
  # check dims
  # combine
  big <- do.call(neuroim2::concat4D, vols)
  neuroim2::write_nifti(big, final_4d)
}

print_summary <- function(all_beta_files, method, subject, task) {
  cat("\n========== BETA ESTIMATION SUMMARY ==========\n")
  cat("Subject:           ", subject, "\n")
  cat("Task:              ", task, "\n")
  cat("Method:            ", method, "\n")
  cat("Number of runs:    ", length(all_beta_files), "\n")
  cat("Output files:      \n")
  for (i in seq_along(all_beta_files)) {
    cat("  ", i, ": ", basename(all_beta_files[i]), "\n")
  }
  if (length(all_beta_files) > 1 && args$concatenate) {
    cat("Concatenated file: ", basename(final_4d), "\n")
  }
  cat("================================================\n")
}

# Call this function before the script ends
print_summary(all_beta_files, method_to_use, args$subid, args$task)

cat("\nDone.\n")