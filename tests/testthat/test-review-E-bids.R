# Review 2026-10, stream E: write_results() hygiene (E3, E4, E6).

.review_e_bids_fit <- function(con = NULL) {
  skip_if_not_installed("neuroim2")
  set.seed(123)
  dims <- c(3, 3, 2); nt <- 50
  scans <- lapply(1:2, function(run) {
    arr <- array(rnorm(prod(dims) * nt), c(dims, nt))
    neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(dim = c(dims, nt)))
  })
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, dims), neuroim2::NeuroSpace(dim = dims))
  ev <- data.frame(
    onset = rep(c(5, 15, 25, 35, 45), 2),
    condition = factor(rep(c("A", "B", "A", "B", "A"), 2)),
    run = rep(1:2, each = 5)
  )
  dset <- neurovec_frame(scans = scans, mask = mask, TR = 1.5, event_table = ev)
  con <- con %||% contrast_set(
    pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B")
  )
  fmri_lm(onset ~ hrf(condition, contrasts = con), block = ~run, dataset = dset,
          durations = 0, strategy = "chunkwise", nchunks = 2, progress = FALSE)
}

.review_e_tmpdir <- function(env = parent.frame()) {
  d <- tempfile("review_e_bids_")
  dir.create(d)
  withr::defer(unlink(d, recursive = TRUE), envir = env)
  d
}

test_that("E3: write_results does not consume the global RNG stream", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("withr")
  fit <- .review_e_bids_fit()
  out <- .review_e_tmpdir()
  set.seed(77)
  before <- get(".Random.seed", envir = globalenv())
  write_results(fit, path = out, subject = "01", task = "t", space = "MNI",
                format = "nifti", contrast_stats = "beta")
  expect_true(identical(get(".Random.seed", envir = globalenv()), before))

  # Staging directories stay unique and are removed afterwards.
  d1 <- fmrireg:::.create_temp_write_dir(out)
  d2 <- fmrireg:::.create_temp_write_dir(out)
  expect_false(identical(d1, d2))
  expect_true(dir.exists(d1) && dir.exists(d2))
  expect_identical(normalizePath(dirname(d1), winslash = "/"),
                   normalizePath(out, winslash = "/"))
  unlink(c(d1, d2), recursive = TRUE)
  expect_false(any(grepl("^\\.tmp_write_", list.files(out, all.files = TRUE))))
})

test_that("E4: NIfTI export works without fmristore; h5 errors clearly", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("withr")
  fit <- .review_e_bids_fit()
  local_mocked_bindings(.has_fmristore = function() FALSE)

  out <- .review_e_tmpdir()
  res <- write_results(fit, path = out, subject = "01", task = "t", space = "MNI",
                       format = "nifti", contrast_stats = c("beta", "tstat"))
  expect_true(file.exists(res$betas$nifti))
  js <- jsonlite::read_json(res$betas$json)
  expect_null(js$SoftwareVersions$fmristore)
  expect_true(file.exists(res$tstat$nifti))

  out2 <- .review_e_tmpdir()
  expect_error(
    write_results(fit, path = out2, subject = "01", task = "t", space = "MNI",
                  format = "h5"),
    "fmristore"
  )
  # Nothing staged or written for the refused request.
  expect_length(list.files(out2, all.files = TRUE, no.. = TRUE), 0L)
})

test_that("E6: DegreesOfFreedom is a scalar when constant, a range otherwise", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("withr")
  fit <- .review_e_bids_fit()
  V <- 18L

  fit$result$df$inference <- rep(97, V)
  out <- .review_e_tmpdir()
  res <- write_results(fit, path = out, subject = "01", task = "t", space = "MNI",
                       format = "nifti", contrast_stats = "tstat")
  js <- jsonlite::read_json(res$tstat$json)
  expect_identical(js$DegreesOfFreedom, 97L)

  fit$result$df$inference <- seq(80, 97, length.out = V)
  out2 <- .review_e_tmpdir()
  res2 <- write_results(fit, path = out2, subject = "01", task = "t", space = "MNI",
                        format = "nifti", strategy = "by_contrast",
                        contrast_stats = "tstat")
  js2 <- jsonlite::read_json(res2$A_vs_B$json)
  expect_equal(js2$DegreesOfFreedom, list(Min = 80, Max = 97))

  expect_identical(fmrireg:::.summarize_df(c(5, 5, 5)), 5)
  expect_equal(fmrireg:::.summarize_df(c(5, NA, 7)), list(Min = 5, Max = 7))
  expect_null(fmrireg:::.summarize_df(NULL))
})

test_that("E6: contrast labels that sanitise to the same string are an error", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("withr")
  con <- contrast_set(
    pair_contrast(~ condition == "A", ~ condition == "B", name = "A_vs_B"),
    pair_contrast(~ condition == "B", ~ condition == "A", name = "A-vs-B")
  )
  fit <- .review_e_bids_fit(con)
  expect_setequal(fit$result$contrasts$name, c("A_vs_B", "A-vs-B"))

  out <- .review_e_tmpdir()
  expect_error(
    write_results(fit, path = out, subject = "01", task = "t", space = "MNI",
                  format = "nifti", strategy = "by_contrast", contrast_stats = "beta"),
    "A_vs_B.*A-vs-B|A-vs-B.*A_vs_B"
  )
  expect_length(list.files(out, all.files = TRUE, no.. = TRUE, recursive = TRUE), 0L)

  # by_stat stores the original names as volume labels: no file-name clash.
  res <- write_results(fit, path = out, subject = "01", task = "t", space = "MNI",
                       format = "nifti", strategy = "by_stat", contrast_stats = "beta")
  expect_true(file.exists(res$beta$nifti))
})

test_that("E6: desc is kept (sanitised) in statistic and contrast file names", {
  skip_if_not_installed("RNifti")
  skip_if_not_installed("withr")
  fit <- .review_e_bids_fit()
  out <- .review_e_tmpdir()

  a <- write_results(fit, path = out, subject = "01", task = "t", space = "MNI",
                     desc = "run-A", format = "nifti", strategy = "by_stat",
                     save_betas = FALSE, contrast_stats = "tstat")
  b <- write_results(fit, path = out, subject = "01", task = "t", space = "MNI",
                     desc = "runB", format = "nifti", strategy = "by_stat",
                     save_betas = FALSE, contrast_stats = "tstat")
  expect_identical(basename(a$tstat$nifti), "sub-01_task-t_space-MNI_desc-runAtstat_bold.nii.gz")
  expect_identical(basename(b$tstat$nifti), "sub-01_task-t_space-MNI_desc-runBtstat_bold.nii.gz")
  expect_true(file.exists(a$tstat$nifti) && file.exists(b$tstat$nifti))

  c1 <- write_results(fit, path = out, subject = "01", task = "t", space = "MNI",
                      desc = "my_model", format = "nifti", strategy = "by_contrast",
                      save_betas = FALSE, contrast_stats = "beta")
  expect_identical(basename(c1$A_vs_B$nifti),
                   "sub-01_task-t_space-MNI_contrast-AvsB_desc-mymodelbeta_bold.nii.gz")

  # betas: desc sanitised too
  d1 <- write_results(fit, path = out, subject = "01", task = "t", space = "MNI",
                      desc = "my_model", format = "nifti", contrast_stats = character(0))
  expect_identical(basename(d1$betas$nifti), "sub-01_task-t_space-MNI_desc-mymodel_betas.nii.gz")
})
