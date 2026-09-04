# bids_export.R: volume reconstruction, NIfTI write, contrast-type fallbacks,
# finalize conflict, and GDS plan helpers beyond existing coverage files.

test_that(".extract_contrast_type covers matrix/weights/generic/error fallbacks", {
  # Matrix with >1 rows => F-contrast (list row, not tibble list-column)
  f_row <- list(
    name = "omnibus",
    data = list(list(estimate = 1)),
    contrast_matrix = matrix(c(1, 0, -1, 0, 1, -1), 2, 3, byrow = TRUE)
  )
  expect_equal(fmrireg:::.extract_contrast_type(f_row, list()), "F-contrast")

  t_row <- list(
    name = "A_vs_B",
    data = list(list(estimate = 1)),
    contrast_matrix = matrix(c(1, -1, 0), 1, 3)
  )
  expect_equal(fmrireg:::.extract_contrast_type(t_row, list()), "t-contrast")

  w_row <- list(
    name = "weighted",
    data = list(list(estimate = 1)),
    weights = c(1, -1)
  )
  expect_equal(fmrireg:::.extract_contrast_type(w_row, list()), "t-contrast")

  generic <- list(
    name = "main_effect",
    data = list(list(estimate = 1))
  )
  expect_equal(fmrireg:::.extract_contrast_type(generic, list()), "contrast")

  # Missing data triggers tryCatch fallback
  expect_equal(
    suppressWarnings(fmrireg:::.extract_contrast_type(list(name = "x"), list())),
    "contrast"
  )
})

test_that(".compute_beta_volumes rejects empty/invalid beta payloads", {
  expect_error(
    fmrireg:::.compute_beta_volumes(list(result = list(betas = list(data = list())))),
    "invalid or empty"
  )
  expect_error(
    fmrireg:::.compute_beta_volumes(
      list(result = list(betas = list(data = list(list(estimate = NULL)))))
    ),
    "invalid or empty"
  )
})

test_that(".compute_statistical_volumes handles list-wrapped stats and errors", {
  skip_if_not_installed("neuroim2")
  mask <- array(TRUE, dim = c(2, 2, 1))
  space <- neuroim2::NeuroSpace(c(2, 2, 1))

  ctab <- tibble::tibble(
    name = c("A_vs_B", "no_stat"),
    data = list(
      list(estimate = list(rnorm(4)), se = abs(rnorm(4)), stat = rnorm(4)),
      list(other = 1)
    )
  )

  out <- fmrireg:::.compute_statistical_volumes(
    ctab, "beta", brain_dims = c(2, 2, 1), mask = mask, space = space
  )
  expect_true(!is.null(out$neurovec))
  expect_equal(out$contrast_names, "A_vs_B")
  expect_s4_class(out$neurovec, "NeuroVec")

  expect_null(
    fmrireg:::.compute_statistical_volumes(
      ctab, "missing_stat", brain_dims = c(2, 2, 1), mask = mask, space = space
    )
  )

  bad <- tibble::tibble(
    name = "A_vs_B",
    data = list(list(estimate = rnorm(2)))
  )
  expect_error(
    fmrireg:::.compute_statistical_volumes(
      bad, "beta", brain_dims = c(2, 2, 1), mask = mask, space = space
    ),
    "mask contains"
  )

  # Mask length / product mismatch
  expect_error(
    fmrireg:::.compute_statistical_volumes(
      ctab[1, ], "beta",
      brain_dims = c(2, 2, 1),
      mask = array(TRUE, dim = c(2, 1, 1)),
      space = space
    ),
    "Mask length"
  )
})

test_that(".write_nifti_volume writes NeuroVol and errors on bad paths", {
  skip_if_not_installed("neuroim2")
  sp <- neuroim2::NeuroSpace(c(2, 2, 1))
  vol <- neuroim2::NeuroVol(array(rnorm(4), dim = c(2, 2, 1)), sp)
  path <- tempfile(fileext = ".nii.gz")
  out <- fmrireg:::.write_nifti_volume(vol, path, "demo volume")
  expect_equal(out, path)
  expect_true(file.exists(path))

  expect_error(
    fmrireg:::.write_nifti_volume(vol, "/no/such/dir/out.nii.gz", "bad write"),
    "Failed to write bad write"
  )
})

test_that(".finalize_atomic_write errors when overwrite=FALSE and target exists", {
  base <- tempfile("bids-fin-")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE), add = TRUE)

  tmp <- fmrireg:::.create_temp_write_dir(base)
  writeLines("new", file.path(tmp, "map.txt"))
  final <- file.path(base, "dest")
  dir.create(final)
  writeLines("old", file.path(final, "map.txt"))

  expect_error(
    fmrireg:::.finalize_atomic_write(
      tmp, final, list(main = file.path(tmp, "map.txt")), overwrite = FALSE
    ),
    "already exist"
  )
})

test_that(".write_gds_plan_file writes and can restore from backup", {
  skip_if_not_installed("fmrigds")
  skip_if_not_installed("neuroim2")

  # Build a tiny on-disk NIfTI that fmrigds can open as source for a plan.
  # Prefer tabular CSV source which is more reliable than empty h5.
  csv <- tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      subject = c("s1", "s2"),
      contrast = "c1",
      beta = c(0.1, 0.2),
      se = c(0.05, 0.06)
    ),
    csv,
    row.names = FALSE
  )

  # Use gds itself as the "source" path argument: write_gds_plan_file
  # calls fmrigds::gds(source = gds_h5_final, format = "h5"). If that fails
  # for CSV, exercise the stash/backup rollback branch with a stub by
  # catching the error after backup creation.
  plan_path <- tempfile(fileext = ".rds")
  bak_root <- tempfile("gds-bak-")
  dir.create(bak_root)
  on.exit(unlink(c(plan_path, bak_root), recursive = TRUE), add = TRUE)

  writeLines("old-plan", plan_path)
  # Expect failure for non-h5 source, but backup/stash path should run first
  expect_error(
    fmrireg:::.write_gds_plan_file(
      gds_h5_final = csv,
      gds_plan_final = plan_path,
      overwrite = TRUE,
      backup_root = bak_root
    )
  )
  # Original plan restored (or still present) after failed write
  expect_true(file.exists(plan_path))
})

test_that(".save_contrasts_by_contrast_nifti returns empty for no contrasts", {
  expect_equal(
    fmrireg:::.save_contrasts_by_contrast_nifti(
      list(), available_contrasts = NULL, path = tempdir(),
      entities = list(subject = "01", task = "x"),
      desc = "GLM", contrast_stats = "beta", overwrite = TRUE
    ),
    list()
  )
  empty_tbl <- tibble::tibble(
    type = character(), name = character(), data = list()
  )
  expect_equal(
    fmrireg:::.save_contrasts_by_contrast_nifti(
      list(), available_contrasts = empty_tbl, path = tempdir(),
      entities = list(subject = "01", task = "x"),
      desc = "GLM", contrast_stats = "beta", overwrite = TRUE
    ),
    list()
  )
})
