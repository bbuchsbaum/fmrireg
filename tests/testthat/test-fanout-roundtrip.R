# M8: the fan-out -> fan-in round trip. Per-subject jobs write BIDS-keyed maps to
# disk (the part this feature owns), then collect_results() reads them back into
# a group_data object that fmri_meta() meta-analyses -- the canonical,
# non-deprecated path.
#
# The group read uses the canonical group_data(format="nifti"), which requires
# fmrigds with BIDS-aware subject keying (so BIDS-named beta/se files pair).
# Environments with an older fmrigds skip the group step.

skip_on_cran()
skip_if_not_installed("bidser")
skip_if_not_installed("RNifti")
skip_if_not_installed("jsonlite")
skip_if_not_installed("fmrigds")

# Probe: does the installed fmrigds pair BIDS-named subject files?
fmrigds_pairs_bids <- function() {
  isTRUE(tryCatch(
    identical(fmrigds:::.nifti_subject_key("sub-01_task-x_desc-beta_bold.nii.gz"),
              "sub-01"),
    error = function(e) FALSE
  ))
}

test_that("fan-out writes BIDS-keyed per-subject maps, and the group loop closes", {
  withr::local_options(fmrireg.suppress_deprecation = TRUE)
  root <- withr::local_tempdir()
  suppressMessages(make_bids_fixture(root))
  proj <- suppressMessages(bidser::bids_project(root, fmriprep = TRUE))
  mask <- file.path(proj$path, "brain_mask.nii.gz")

  mani <- suppressWarnings(from_bids(
    proj, task = "stroop", space = "MNI152NLin2009cAsym",
    confounds = c("trans_x", "trans_y", "rot_z"), mask = mask
  ))

  # A formula-embedded contrast makes write_results emit beta AND se maps
  # (required to meta-analyse: the group step needs per-subject variance).
  con <- contrast_set(pair_contrast(~ trial_type == "incongruent",
                                    ~ trial_type == "congruent",
                                    name = "incong_gt_cong"))
  outdir <- file.path(root, "derivatives", "glm")
  tmpl <- fmri_template(
    onset ~ hrf(trial_type, contrasts = con), ~ run,
    baseline = baseline_spec(degree = 3),
    reducer  = reduce_write_results(format = "nifti", stats = c("beta", "se"),
                                    path = outdir)
  )

  res <- suppressWarnings(run_jobs(instantiate(tmpl, mani)))
  expect_true(all(res$ok), info = paste(batch_errors(res), collapse = "; "))

  # --- the fan-out contract this feature owns: correct BIDS-keyed maps on disk ---
  beta_paths <- sort(list.files(outdir, pattern = "desc-beta_bold\\.nii\\.gz$",
                                full.names = TRUE))
  se_paths <- sort(list.files(outdir, pattern = "desc-se_bold\\.nii\\.gz$",
                              full.names = TRUE))
  expect_length(beta_paths, 2)        # one per subject
  expect_length(se_paths, 2)
  expect_match(basename(beta_paths[1]), "^sub-01_task-stroop_space-MNI152")

  # --- the loop closes: collect_results -> canonical group_data -> fmri_meta ---
  skip_if_not(fmrigds_pairs_bids(),
              "installed fmrigds lacks BIDS-aware subject keying")

  gd <- collect_results(outdir, space = "MNI152NLin2009cAsym")
  expect_s3_class(gd, "group_data")
  expect_equal(n_subjects(gd), 2L)
  expect_setequal(get_subjects(gd), c("sub-01", "sub-02"))

  fm <- fmri_meta(gd, formula = ~ 1, method = "fe", verbose = FALSE)
  expect_s3_class(fm, "fmri_meta")
  expect_equal(dim(fm$coefficients), c(64L, 1L))      # 4x4x4 mask, intercept-only
  expect_equal(colnames(fm$coefficients), "(Intercept)")
  expect_true(all(is.finite(fm$coefficients)))
  expect_true(all(is.finite(fm$se)))
})
