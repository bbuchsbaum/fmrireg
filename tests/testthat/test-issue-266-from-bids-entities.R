# from_bids() pairs files by exact BIDS entities (bbuchsbaum/fmrireg#266).
# Before the fix, labels were matched as unanchored regexes and scans, events
# and confounds were paired by run number alone.

skip_if_not_installed("bidser")
skip_if_not_installed("RNifti")
skip_if_not_installed("jsonlite")

fixture_266 <- function(..., edit = NULL) {
  root <- withr::local_tempdir(.local_envir = parent.frame())
  suppressMessages(make_bids_fixture(root, ...))
  if (!is.null(edit)) edit(root)
  root
}

open_266 <- function(root) suppressMessages(bidser::bids_project(root, fmriprep = TRUE))

all_files <- function(root, pattern) {
  list.files(root, pattern, recursive = TRUE, full.names = TRUE)
}

rename_all <- function(root, from, to) {
  for (f in all_files(root, from)) {
    file.rename(f, file.path(dirname(f), gsub(from, to, basename(f))))
  }
}

space <- "MNI152NLin2009cAsym"

test_that("subject labels match exactly: sub-1 is not sub-10", {
  root <- fixture_266(subjects = c("1", "10"))
  mani <- from_bids(open_266(root), task = "stroop", space = space, subjects = "1")
  expect_length(mani, 1L)
  b <- mani[[1]]
  expect_equal(b$id, "sub-1")
  expect_length(b$scans, 2L)                         # was 4: sub-10's runs leaked in
  expect_true(all(grepl("^sub-1_", basename(b$scans))))
  expect_equal(b$run_length, c(40L, 40L))
  expect_equal(sort(unique(b$events$run)), 1:2)
  expect_equal(nrow(b$events), 10L)
})

test_that("task labels match exactly: stroop is not stroopx", {
  root <- fixture_266(edit = function(root) {
    for (f in all_files(root, "task-stroop_")) {
      file.copy(f, file.path(dirname(f), sub("task-stroop_", "task-stroopx_", basename(f))))
    }
  })
  mani <- from_bids(open_266(root), task = "stroop", space = space,
                    confounds = "trans_x")
  for (b in mani) {
    expect_length(b$scans, 2L)
    expect_false(any(grepl("stroopx", b$scans)))
    expect_equal(nrow(b$events), 10L)
    expect_length(b$confounds, 2L)
  }
})

test_that("a second resolution is an error unless `res` picks one", {
  root <- fixture_266(edit = function(root) {
    for (f in all_files(root, "desc-preproc_bold")) {
      file.copy(f, file.path(dirname(f), sub("_desc-", "_res-2_desc-", basename(f))))
    }
  })
  proj <- open_266(root)
  expect_error(from_bids(proj, task = "stroop", space = space),
               "several preprocessed series for one run.*res = ")
  mani <- from_bids(proj, task = "stroop", space = space, res = "2")
  expect_true(all(grepl("_res-2_", unlist(lapply(mani, `[[`, "scans")))))
  expect_length(mani[[1]]$scans, 2L)
})

test_that("sessions are runs of their own, paired with their own events and confounds", {
  root <- fixture_266(subjects = "01", edit = function(root) {
    for (base in c(file.path(root, "sub-01"),
                   file.path(root, "derivatives", "fmriprep", "sub-01"))) {
      for (ses in c("1", "2")) {
        dest <- file.path(base, paste0("ses-", ses))
        dir.create(dest)
        file.copy(file.path(base, "func"), dest, recursive = TRUE)
        rename_all(dest, "^sub-01_", paste0("sub-01_ses-", ses, "_"))
      }
      unlink(file.path(base, "func"), recursive = TRUE)
    }
    # Mark session 2's events and confounds so a mis-pairing shows.
    for (f in all_files(root, "ses-2_.*_events\\.tsv$")) {
      ev <- utils::read.delim(f)
      ev$onset <- ev$onset + 1
      utils::write.table(ev, f, sep = "\t", quote = FALSE, row.names = FALSE)
    }
    for (f in all_files(root, "ses-2_.*_desc-confounds_timeseries\\.tsv$")) {
      cf <- utils::read.delim(f, na.strings = "n/a")
      cf$trans_x <- 100 + seq_len(nrow(cf))
      utils::write.table(cf, f, sep = "\t", quote = FALSE, row.names = FALSE, na = "n/a")
    }
  })
  mani <- from_bids(open_266(root), task = "stroop", space = space,
                    confounds = "trans_x")
  b <- mani[[1]]
  expect_length(b$scans, 4L)
  expect_equal(b$meta$runs$session, c("1", "1", "2", "2"))
  expect_equal(b$meta$runs$run, c(1L, 2L, 1L, 2L))
  expect_true(all(mapply(grepl, c("ses-1_.*run-1", "ses-1_.*run-2",
                                  "ses-2_.*run-1", "ses-2_.*run-2"),
                         basename(b$scans))))
  expect_equal(sort(unique(b$events$run)), 1:4)
  first_onset <- as.vector(tapply(b$events$onset, b$events$run, min))
  expect_equal(first_onset, c(0, 0, 1, 1))
  expect_true(all(b$confounds[[1]] < 100) && all(b$confounds[[2]] < 100))
  expect_true(all(b$confounds[[3]] > 100) && all(b$confounds[[4]] > 100))

  only2 <- from_bids(open_266(root), task = "stroop", space = space, sessions = "2")
  expect_equal(only2[[1]]$meta$runs$session, c("2", "2"))
})

test_that("output sidecars under derivatives/ never supply the TR", {
  root <- fixture_266(edit = function(root) {
    out <- file.path(root, "derivatives", "fmrireg", "sub-01")
    dir.create(out, recursive = TRUE)
    writeLines('{"Name": "fmrireg", "BIDSVersion": "1.8.0", "DatasetType": "derivative", "GeneratedBy": [{"Name": "fmrireg"}]}',
               file.path(dirname(out), "dataset_description.json"))
    stem <- sprintf("sub-01_task-stroop_space-%s", space)
    writeLines('{"RepetitionTime": 99, "ContrastOrder": ["a", "b"]}',
               file.path(out, paste0(stem, "_desc-beta_bold.json")))
    writeLines('{"RepetitionTime": 98, "ContrastOrder": "a"}',
               file.path(out, paste0(stem, "_desc-tstat_bold.json")))
  })
  mani <- from_bids(open_266(root), task = "stroop", space = space)
  expect_equal(vapply(mani, `[[`, numeric(1), "TR"), c(2, 2))
})

test_that("each run's functional mask is paired; anatomical masks stay out", {
  root <- fixture_266(subjects = "01", edit = function(root) {
    prep <- file.path(root, "derivatives", "fmriprep", "sub-01")
    for (f in all_files(prep, "desc-preproc_bold\\.nii\\.gz$")) {
      RNifti::writeNifti(array(1L, c(4, 4, 4)),
                         sub("_desc-preproc_bold", "_desc-brain_mask", f))
    }
    dir.create(file.path(prep, "anat"))
    RNifti::writeNifti(array(1L, c(8, 8, 8)),
                       file.path(prep, "anat", sprintf("sub-01_space-%s_desc-brain_mask.nii.gz", space)))
  })
  b <- from_bids(open_266(root), task = "stroop", space = space)[[1]]
  expect_length(b$meta$runs$mask, 2L)
  expect_true(all(file.exists(b$meta$runs$mask)))
  expect_true(all(mapply(grepl, c("run-1_space-.*_desc-brain_mask",
                                  "run-2_space-.*_desc-brain_mask"),
                         basename(b$meta$runs$mask))))
  expect_false(any(grepl("/anat/", b$meta$runs$mask)))
})

test_that("duplicate events files outside the selection do not block it", {
  root <- fixture_266(edit = function(root) {
    f <- all_files(file.path(root, "sub-02"), "_run-1_events\\.tsv$")
    file.copy(f, sub("_run-1_", "_acq-x_run-1_", f))
  })
  proj <- open_266(root)
  mani <- from_bids(proj, task = "stroop", space = space, subjects = "01")
  expect_length(mani[[1]]$scans, 2L)
  expect_error(from_bids(proj, task = "stroop", space = space, subjects = "02"),
               "several events files for one run")
})
