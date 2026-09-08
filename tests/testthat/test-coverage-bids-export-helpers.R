# Coverage for bids_export.R helpers / validation branches beyond test_write_results.R.

make_contrast_tbl <- function() {
  tibble::tibble(
    type = c("contrast", "contrast", "Fcontrast"),
    name = c("A_vs_B", "A_vs_C", "main_effect"),
    stat_type = c("tstat", "tstat", "Fstat")
  )
}

test_that(".validate_fmrilm_object and entity helpers reject incomplete inputs", {
  expect_error(fmrireg:::.validate_fmrilm_object(list()), "fmri_lm")
  expect_error(
    fmrireg:::.validate_fmrilm_object(structure(list(), class = "fmri_lm")),
    "not been fitted"
  )
  expect_error(
    fmrireg:::.validate_fmrilm_object(
      structure(list(result = list()), class = "fmri_lm")
    ),
    "beta estimates"
  )
  expect_error(
    fmrireg:::.validate_fmrilm_object(
      structure(list(result = list(betas = 1)), class = "fmri_lm")
    ),
    "model specification"
  )
  expect_error(
    fmrireg:::.validate_fmrilm_object(
      structure(
        list(result = list(betas = 1), model = list()),
        class = "fmri_lm"
      )
    ),
    "dataset information"
  )

  expect_error(
    fmrireg:::.validate_required_entities(list(task = "nback")),
    "Subject identifier"
  )
  expect_error(
    fmrireg:::.validate_required_entities(list(subject = "01")),
    "Task identifier"
  )
  expect_warning(
    fmrireg:::.validate_required_entities(list(subject = "01", task = "nback")),
    "Spatial reference"
  )

  entities <- fmrireg:::.create_bids_entities("sub-01!", "my-task", "MNI152")
  expect_equal(entities$subject, "sub01")
  expect_equal(entities$task, "mytask")
  expect_equal(entities$space, "MNI152")

  expect_null(fmrireg:::.sanitize_label(NULL))
  expect_warning(fmrireg:::.sanitize_label("!!!"), "empty string")
  expect_null(suppressWarnings(fmrireg:::.sanitize_label("!!!")))

  fname <- fmrireg:::.generate_bids_filename(
    entities,
    desc = "GLM",
    contrast = "A_vs_B",
    stat = "tstat",
    suffix = "statmap",
    extension = "h5"
  )
  expect_match(fname, "^sub-sub01_task-mytask_space-MNI152")
  expect_match(fname, "contrast-AvsB")
  expect_match(fname, "\\.h5$")

  expect_equal(fmrireg:::.stat_desc("GLM", "tstat"), "tstat")
  expect_equal(fmrireg:::.image_output_formats("h5"), "h5")
  expect_equal(
    fmrireg:::.image_output_formats(c("nifti", "h5")),
    c("h5", "nifti")
  )
  expect_equal(fmrireg:::.image_output_formats(c("gds")), character(0))
})

test_that(".select_contrasts covers exact/auto/regex and invalid regex", {
  ctab <- make_contrast_tbl()

  expect_equal(nrow(fmrireg:::.select_contrasts(ctab, NULL)), 3L)

  exact <- fmrireg:::.select_contrasts(ctab, "A_vs_B", contrast_match = "exact")
  expect_equal(exact$name, "A_vs_B")

  # exact match mode ignores partial regex-like selectors
  none <- fmrireg:::.select_contrasts(ctab, "A_vs", contrast_match = "exact")
  expect_equal(nrow(none), 0L)

  auto <- fmrireg:::.select_contrasts(ctab, "A_vs_B", contrast_match = "auto")
  expect_equal(auto$name, "A_vs_B")

  regex <- fmrireg:::.select_contrasts(ctab, "A_vs_", contrast_match = "regex")
  expect_equal(sort(regex$name), c("A_vs_B", "A_vs_C"))

  expect_error(
    fmrireg:::.select_contrasts(ctab, "(unclosed", contrast_match = "regex"),
    "Invalid contrast regex"
  )
})

test_that(".resolve_contrast_families and .merge_created_files behave", {
  ctab <- make_contrast_tbl()
  families <- fmrireg:::.resolve_contrast_families(
    ctab,
    contrast_stats = c("beta", "tstat", "pval", "se", "fstat", "fpval")
  )
  expect_true(length(families) >= 1L)
  types <- vapply(families, function(f) unique(f$contrasts$type), character(1))
  expect_true("contrast" %in% types || "Fcontrast" %in% types)

  # F-only requested stats drop t-family when only F stats remain applicable
  f_only <- fmrireg:::.resolve_contrast_families(ctab, contrast_stats = c("fstat", "fpval"))
  f_types <- unlist(lapply(f_only, function(f) unique(as.character(f$contrasts$type))))
  expect_true(all(f_types == "Fcontrast") || length(f_only) >= 1L)

  # Empty / inapplicable stats yield empty family list
  empty <- fmrireg:::.resolve_contrast_families(ctab[0, ], c("beta", "tstat"))
  expect_equal(length(empty), 0L)

  merged <- fmrireg:::.merge_created_files(
    list(betas = "a.h5", contrasts = list(c0 = "old.h5")),
    list(betas = "b.h5", contrasts = list(c1 = "c.h5"))
  )
  expect_equal(merged$betas, "b.h5")
  expect_equal(merged$contrasts$c0, "old.h5")
  expect_equal(merged$contrasts$c1, "c.h5")

  expect_equal(
    fmrireg:::.merge_created_files(list(a = 1), list()),
    list(a = 1)
  )
  merged2 <- fmrireg:::.merge_created_files(list(), list(x = 1:2))
  expect_equal(merged2$x, 1:2)
})

test_that(".predict_output_files builds beta and contrast predictions", {
  ctab <- tibble::tibble(
    type = c("contrast", "Fcontrast"),
    name = c("A_vs_B", "main"),
    data = list(
      tibble::tibble(estimate = 1, se = 0.2, stat = 2, prob = 0.05),
      tibble::tibble(estimate = 1, se = 0.2, stat = 3, prob = 0.01)
    )
  )
  fake_lm <- structure(
    list(
      result = list(
        betas = TRUE,
        contrasts = ctab
      ),
      model = list(),
      dataset = list()
    ),
    class = "fmri_lm"
  )

  entities <- list(subject = "01", task = "test", space = "MNI")
  out_dir <- tempfile("bids-pred-")
  predicted <- fmrireg:::.predict_output_files(
    path = out_dir,
    entities = entities,
    desc = "GLM",
    format = "h5",
    strategy = "by_stat",
    save_betas = TRUE,
    contrasts = NULL,
    contrast_match = "auto",
    contrast_stats = c("beta", "tstat", "fstat"),
    fmrilm_obj = fake_lm
  )
  expect_true(is.character(predicted))
  expect_true(any(grepl("betas\\.h5$", predicted)))
  expect_true(any(grepl("betas\\.json$", predicted)))
  expect_true(any(grepl("bold\\.(h5|json)$", predicted)))

  predicted2 <- fmrireg:::.predict_output_files(
    path = out_dir,
    entities = entities,
    desc = "GLM",
    format = c("h5", "nifti"),
    strategy = "by_contrast",
    save_betas = FALSE,
    contrasts = "A_vs_B",
    contrast_match = "exact",
    contrast_stats = c("beta", "tstat"),
    fmrilm_obj = fake_lm
  )
  expect_true(any(grepl("contrast-AvsB", predicted2)))
  expect_true(any(grepl("\\.nii\\.gz$", predicted2)))

  predicted_gds <- fmrireg:::.predict_output_files(
    path = out_dir,
    entities = entities,
    desc = "GLM",
    format = "gds",
    strategy = "by_stat",
    save_betas = FALSE,
    contrasts = NULL,
    contrast_match = "auto",
    contrast_stats = "beta",
    fmrilm_obj = fake_lm
  )
  expect_true(any(grepl("gds\\.h5$", predicted_gds)))
  expect_true(any(grepl("gds\\.rds$", predicted_gds)))
})
