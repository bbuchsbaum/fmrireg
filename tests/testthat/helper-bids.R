# Build a tiny, self-contained, OFFLINE BIDS dataset (raw + fmriprep
# derivatives) in `root`, for from_bids() tests. Verified to be fully indexed by
# bidser::bids_project(root, fmriprep = TRUE). Requires RNifti + jsonlite.
make_bids_fixture <- function(root,
                              subjects = c("01", "02"),
                              task = "stroop",
                              runs = c(1, 2),
                              nvols = 40L,
                              dim_xyz = c(4L, 4L, 4L),
                              TR = 2.0,
                              space = "MNI152NLin2009cAsym") {
  stopifnot(requireNamespace("RNifti", quietly = TRUE),
            requireNamespace("jsonlite", quietly = TRUE))
  dir.create(root, recursive = TRUE, showWarnings = FALSE)

  write_nii <- function(path, nv = nvols) {
    arr <- array(stats::rnorm(prod(dim_xyz) * nv), dim = c(dim_xyz, nv))
    RNifti::writeNifti(arr, path)
  }
  write_json <- function(path, x)
    writeLines(jsonlite::toJSON(x, auto_unbox = TRUE, pretty = TRUE), path)
  write_tsv <- function(path, df)
    utils::write.table(df, path, sep = "\t", quote = FALSE,
                       row.names = FALSE, col.names = TRUE, na = "n/a")

  # a simple 3D brain mask (all-in) for file-backed datasets
  RNifti::writeNifti(array(1L, dim = dim_xyz), file.path(root, "brain_mask.nii.gz"))

  write_json(file.path(root, "dataset_description.json"),
             list(Name = "Synthetic Stroop", BIDSVersion = "1.8.0",
                  DatasetType = "raw"))
  write_tsv(file.path(root, "participants.tsv"),
            data.frame(participant_id = paste0("sub-", subjects),
                       stringsAsFactors = FALSE))

  prep_dir <- file.path(root, "derivatives", "fmriprep")
  dir.create(prep_dir, recursive = TRUE, showWarnings = FALSE)
  write_json(file.path(prep_dir, "dataset_description.json"),
             list(Name = "fMRIPrep - synthetic", BIDSVersion = "1.8.0",
                  DatasetType = "derivative",
                  GeneratedBy = list(list(Name = "fMRIPrep", Version = "23.0.0"))))

  confound_cols <- c("trans_x", "trans_y", "trans_z", "rot_x", "rot_y", "rot_z",
                     "csf", "white_matter", "global_signal", "framewise_displacement")

  for (s in subjects) {
    rawfunc <- file.path(root, paste0("sub-", s), "func")
    depfunc <- file.path(prep_dir, paste0("sub-", s), "func")
    dir.create(rawfunc, recursive = TRUE, showWarnings = FALSE)
    dir.create(depfunc, recursive = TRUE, showWarnings = FALSE)

    for (r in runs) {
      stem <- sprintf("sub-%s_task-%s_run-%d", s, task, r)
      write_nii(file.path(rawfunc, paste0(stem, "_bold.nii.gz")))
      write_json(file.path(rawfunc, paste0(stem, "_bold.json")),
                 list(RepetitionTime = TR, TaskName = task))
      ev <- data.frame(onset = seq(0, by = 8, length.out = 5),
                       duration = rep(1.5, 5),
                       trial_type = rep(c("congruent", "incongruent"), length.out = 5),
                       stringsAsFactors = FALSE)
      write_tsv(file.path(rawfunc, paste0(stem, "_events.tsv")), ev)

      write_nii(file.path(depfunc,
                          sprintf("%s_space-%s_desc-preproc_bold.nii.gz", stem, space)))
      cf <- as.data.frame(matrix(round(stats::rnorm(nvols * length(confound_cols)), 4),
                                 nrow = nvols, dimnames = list(NULL, confound_cols)))
      cf$framewise_displacement[1] <- NA
      write_tsv(file.path(depfunc, paste0(stem, "_desc-confounds_timeseries.tsv")), cf)
    }
  }
  invisible(root)
}
