# More bids_export.R helpers: atomic write, metadata, stat mapping.

test_that(".map_stat_name / .get_stat_units / .image_data_info cover mappings", {
  expect_equal(fmrireg:::.map_stat_name("beta"), "estimate")
  expect_equal(fmrireg:::.map_stat_name("tstat"), "stat")
  expect_equal(fmrireg:::.map_stat_name("fstat"), "stat")
  expect_equal(fmrireg:::.map_stat_name("pval"), "prob")
  expect_equal(fmrireg:::.map_stat_name("fpval"), "prob")
  expect_equal(fmrireg:::.map_stat_name("se"), "se")
  expect_equal(fmrireg:::.map_stat_name("sigma"), "sigma")
  expect_equal(fmrireg:::.map_stat_name("custom_stat"), "custom_stat")

  expect_match(fmrireg:::.get_stat_units("beta"), "arbitrary")
  expect_equal(fmrireg:::.get_stat_units("tstat"), "t-statistic")
  expect_equal(fmrireg:::.get_stat_units("pval"), "probability")
  expect_equal(fmrireg:::.get_stat_units("fpval"), "probability")
  expect_equal(fmrireg:::.get_stat_units("se"), "standard error")
  expect_equal(fmrireg:::.get_stat_units("z"), "z-score")
  expect_equal(fmrireg:::.get_stat_units("f"), "F-statistic")
  expect_equal(fmrireg:::.get_stat_units("fstat"), "F-statistic")
  expect_equal(fmrireg:::.get_stat_units("unknown"), "arbitrary")

  info <- fmrireg:::.image_data_info(c("h5", "nifti"), "t-statistic")
  expect_equal(info$Units, "t-statistic")
  expect_true(all(c("HDF5", "NIfTI") %in% info$FileFormat))

  info_empty <- fmrireg:::.image_data_info(character(0), "x")
  expect_equal(info_empty$FileFormat, "HDF5")
})

test_that(".extract_contrast_type and .extract_degrees_of_freedom branches", {
  row_t <- tibble::tibble(
    type = "contrast",
    name = "A_vs_B",
    data = list(tibble::tibble(estimate = 1, se = 0.1, stat = 2, prob = 0.05))
  )
  expect_equal(fmrireg:::.extract_contrast_type(row_t, list()), "t-contrast")

  row_f <- tibble::tibble(
    type = "Fcontrast",
    name = "main",
    data = list(tibble::tibble(estimate = 1, se = 1, stat = 3, prob = 0.01))
  )
  expect_equal(fmrireg:::.extract_contrast_type(row_f, list()), "F-contrast")

  # Infer from data fields when type missing
  row_infer_f <- tibble::tibble(
    name = "omnibus",
    data = list(list(fstat = 2.5, estimate = 1))
  )
  expect_equal(fmrireg:::.extract_contrast_type(row_infer_f, list()), "F-contrast")

  row_infer_t <- tibble::tibble(
    name = "A_vs_B",
    data = list(list(tstat = 2, estimate = 1))
  )
  expect_equal(fmrireg:::.extract_contrast_type(row_infer_t, list()), "t-contrast")

  row_name <- tibble::tibble(
    name = "Faces_vs_Places",
    data = list(list(estimate = 1))
  )
  expect_equal(fmrireg:::.extract_contrast_type(row_name, list()), "t-contrast")

  expect_equal(
    fmrireg:::.extract_degrees_of_freedom(list(result = list(df = list(inference = 40)))),
    40
  )
  # Prefer betas$df.residual when result$df is absent. Avoid a bare
  # result$df.residual alone: result$df partial-matches it and errors into NULL.
  expect_equal(
    fmrireg:::.extract_degrees_of_freedom(
      list(result = list(betas = data.frame(df.residual = 22)))
    ),
    22
  )
  expect_equal(
    fmrireg:::.extract_degrees_of_freedom(
      list(result = list(contrasts = data.frame(df.residual = 18)))
    ),
    18
  )
  expect_null(fmrireg:::.extract_degrees_of_freedom(list(result = list())))
})

test_that("atomic write helpers create, stash, finalize, cleanup, rollback", {
  base <- tempfile("bids-atomic-")
  dir.create(base)
  on.exit(unlink(base, recursive = TRUE), add = TRUE)

  tmp <- fmrireg:::.create_temp_write_dir(base)
  expect_true(dir.exists(tmp))
  expect_match(basename(tmp), "^\\.tmp_write_")

  staged <- file.path(tmp, "out.txt")
  writeLines("v1", staged)

  # Relative paths must stay under staging root
  rel <- fmrireg:::.relative_staged_paths(staged, tmp)
  expect_equal(rel, "out.txt")
  expect_error(
    fmrireg:::.relative_staged_paths("/tmp/outside.txt", tmp),
    "outside the staged"
  )

  final_path <- file.path(base, "final")
  dir.create(final_path)
  created <- list(main = staged)
  finalized <- fmrireg:::.finalize_atomic_write(
    tmp, final_path, created, overwrite = FALSE
  )
  expect_true(file.exists(file.path(final_path, "out.txt")))
  expect_equal(finalized$files$main, file.path(final_path, "out.txt"))

  # Stash existing + overwrite path
  tmp2 <- fmrireg:::.create_temp_write_dir(base)
  writeLines("v2", file.path(tmp2, "out.txt"))
  # Also create nested
  dir.create(file.path(tmp2, "sub"), recursive = TRUE)
  writeLines("nested", file.path(tmp2, "sub", "x.txt"))
  created2 <- list(
    main = file.path(tmp2, "out.txt"),
    nested = list(a = file.path(tmp2, "sub", "x.txt"))
  )
  finalized2 <- fmrireg:::.finalize_atomic_write(
    tmp2, final_path, created2, overwrite = TRUE
  )
  expect_equal(readLines(file.path(final_path, "out.txt")), "v2")
  expect_true(file.exists(file.path(final_path, "sub", "x.txt")))
  expect_equal(finalized2$files$nested$a, file.path(final_path, "sub", "x.txt"))

  # stash without overwrite errors
  existing <- file.path(final_path, "out.txt")
  expect_error(
    fmrireg:::.stash_existing_file(existing, file.path(base, "bak"), overwrite = FALSE),
    "already exists"
  )
  expect_null(fmrireg:::.stash_existing_file(
    file.path(base, "missing.txt"), file.path(base, "bak"), overwrite = TRUE
  ))

  # Rollback transaction: remove moved file, restore backup
  bak_root <- file.path(base, "rollback_bak")
  dir.create(bak_root)
  target <- file.path(final_path, "rollback_me.txt")
  writeLines("current", target)
  backup <- file.path(bak_root, "rollback_me.txt")
  writeLines("old", backup)
  fmrireg:::.rollback_finalized_write(list(
    moved_files = target,
    backup_map = stats::setNames(backup, target)
  ))
  expect_equal(readLines(target), "old")

  # Cleanup
  leftover <- fmrireg:::.create_temp_write_dir(base)
  writeLines("x", file.path(leftover, "z.txt"))
  fmrireg:::.cleanup_temp_write_dir(leftover)
  expect_false(dir.exists(leftover))

  # Path rewriter
  updated <- fmrireg:::.update_file_paths_in_results(
    list(a = file.path(tmp2, "ghost.txt"), b = list(c = "/abs/keep.txt")),
    temp_dir = tmp2,
    final_path = final_path
  )
  expect_equal(updated$b$c, "/abs/keep.txt")
})

test_that(".fmri_lm_bids_inference_metadata and .release_unused_hdf5_handles", {
  expect_null(fmrireg:::.fmri_lm_bids_inference_metadata(list(result = list())))

  vm <- structure(
    list(
      method = "model",
      covariance_scope = "run",
      df_inference = c(40, 41, 39),
      df_nominal = 40,
      metadata = list(
        df_method = "residual",
        estimation_scope = "joint",
        noise = list(struct = "ar1"),
        robust = list(type = "none"),
        selected_max_lag = 2L,
        taper = "none"
      )
    ),
    class = "fmri_lm_variance_model"
  )
  meta <- fmrireg:::.fmri_lm_bids_inference_metadata(list(result = list(variance_model = vm)))
  expect_equal(meta$VarianceMethod, "model")
  expect_equal(meta$InferenceDegreesOfFreedom$Min, 39)
  expect_equal(meta$InferenceDegreesOfFreedom$Max, 41)
  expect_equal(meta$NoiseStructure, "ar1")

  expect_true(fmrireg:::.release_unused_hdf5_handles(collector = function() NULL))
})

test_that(".extract_source_files reads dataset attributes", {
  expect_equal(fmrireg:::.extract_source_files(list()), character(0))

  obj <- list(dataset = structure(list(), source_files = c("a.nii", "b.nii")))
  attr(obj$dataset, "source_files") <- c("a.nii", "b.nii")
  expect_equal(fmrireg:::.extract_source_files(obj), c("a.nii", "b.nii"))

  obj2 <- list(dataset = structure(list(), class = "fmri_dataset"))
  attr(obj2$dataset, "file_path") <- "/data/sub-01_bold.nii.gz"
  src <- fmrireg:::.extract_source_files(obj2)
  expect_match(src, "bids::")
  expect_match(src, "sub-01_bold")
})
