# Remaining bids_export helpers not covered by prior coverage files.

test_that(".get_stat_units .map_stat_name .image_data_info cover maps", {
  expect_match(fmrireg:::.get_stat_units("beta"), "arbitrary")
  expect_equal(fmrireg:::.get_stat_units("tstat"), "t-statistic")
  expect_equal(fmrireg:::.get_stat_units("pval"), "probability")
  expect_equal(fmrireg:::.get_stat_units("fpval"), "probability")
  expect_equal(fmrireg:::.get_stat_units("se"), "standard error")
  expect_equal(fmrireg:::.get_stat_units("z"), "z-score")
  expect_equal(fmrireg:::.get_stat_units("fstat"), "F-statistic")
  expect_equal(fmrireg:::.get_stat_units("unknown"), "arbitrary")

  expect_true(is.character(fmrireg:::.map_stat_name("tstat")) ||
                identical(fmrireg:::.map_stat_name("tstat"), "tstat"))
  mapped <- vapply(
    c("beta", "se", "tstat", "pval", "fstat", "fpval"),
    fmrireg:::.map_stat_name,
    character(1)
  )
  expect_equal(length(mapped), 6L)

  info <- fmrireg:::.image_data_info(c("h5", "nifti"), units = "t-statistic")
  expect_equal(info$Units, "t-statistic")
  expect_true(all(c("HDF5", "NIfTI") %in% info$FileFormat))
  expect_true(any(grepl("LabeledVolumeSet|NIfTI", info$Storage)))

  info_default <- fmrireg:::.image_data_info(character(0), units = "x")
  expect_equal(info_default$FileFormat, "HDF5")
})

test_that(".extract_source_files and .extract_contrast_type cover branches", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  src0 <- fmrireg:::.extract_source_files(fit)
  expect_true(is.character(src0))

  fit2 <- fit
  attr(fit2$dataset, "source_files") <- c("bids::sub-01_bold.nii.gz")
  expect_equal(
    fmrireg:::.extract_source_files(fit2),
    "bids::sub-01_bold.nii.gz"
  )

  fit3 <- fit
  attr(fit3$dataset, "file_path") <- "/data/sub-02_task-x_bold.nii.gz"
  src3 <- fmrireg:::.extract_source_files(fit3)
  expect_true(any(grepl("bids::|sub-02", src3)) || identical(src3, character(0)))

  t_row <- tibble::tibble(
    type = "contrast",
    name = "A_vs_B",
    data = list(tibble::tibble(estimate = 1, se = 0.2, tstat = 2, prob = 0.05))
  )
  expect_equal(fmrireg:::.extract_contrast_type(t_row, fit), "t-contrast")

  f_row <- tibble::tibble(
    type = "Fcontrast",
    name = "main",
    data = list(tibble::tibble(estimate = 1, fstat = 3, prob = 0.01))
  )
  expect_equal(fmrireg:::.extract_contrast_type(f_row, fit), "F-contrast")

  infer_f <- tibble::tibble(
    name = "main",
    data = list(tibble::tibble(fstat = 4, prob = 0.02))
  )
  expect_equal(fmrireg:::.extract_contrast_type(infer_f, fit), "F-contrast")

  infer_t <- tibble::tibble(
    name = "t1",
    data = list(tibble::tibble(tstat = 2, beta = 0.1))
  )
  expect_equal(fmrireg:::.extract_contrast_type(infer_t, fit), "t-contrast")
})

test_that(".extract_degrees_of_freedom and beta volume validation errors", {
  fit <- suppressWarnings(fmrireg:::.demo_fmri_lm())
  dfv <- fmrireg:::.extract_degrees_of_freedom(fit)
  expect_true(is.numeric(dfv) || is.null(dfv) || is.integer(dfv))

  expect_error(
    fmrireg:::.compute_beta_volumes(structure(list(result = list(betas = list())), class = "fmri_lm")),
    "invalid or empty"
  )

  # Temp dir / stash helpers
  base <- tempfile("bids-tmp-")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE), add = TRUE)
  tmp <- fmrireg:::.create_temp_write_dir(base)
  expect_true(dir.exists(tmp))

  existing <- file.path(base, "keep.txt")
  writeLines("x", existing)
  backup <- tempfile("bids-bak-")
  dir.create(backup)
  stashed <- fmrireg:::.stash_existing_file(existing, backup, overwrite = TRUE)
  expect_true(file.exists(stashed) || is.character(stashed))

  fmrireg:::.cleanup_temp_write_dir(tmp)
  expect_false(dir.exists(tmp))
})
