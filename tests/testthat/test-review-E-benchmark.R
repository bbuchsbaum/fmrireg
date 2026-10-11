# Review 2026-10, stream E5: benchmark data must not hard-fail when the
# installed fmrihrf differs from the version recorded in the HRF recipe.

test_that("E5: fmrihrf version mismatch warns and still reconstructs the HRF", {
  recipe <- load_benchmark_dataset("metadata")$hrf_recipes$variant1
  ref <- fmrireg:::.reconstruct_hrf_object(recipe)
  grid <- seq(0, 30, by = 0.1)

  local_mocked_bindings(.installed_fmrihrf_version = function() "99.0.0")
  expect_warning(
    hrf <- fmrireg:::.reconstruct_hrf_object(recipe),
    "recorded with fmrihrf .*installed version is 99.0.0"
  )
  expect_equal(fmrihrf::evaluate(hrf, grid), fmrihrf::evaluate(ref, grid))

  # Semantic recipe checks still fail closed.
  bad <- recipe
  bad$lag <- recipe$lag + 1
  expect_error(suppressWarnings(fmrireg:::.reconstruct_hrf_object(bad)),
               "inconsistent fields")
})

test_that("E5: loading a dataset under a mismatched fmrihrf version works", {
  local_mocked_bindings(.installed_fmrihrf_version = function() "99.0.0")
  env <- fmrireg:::.benchmark_data_env
  had <- exists("fmri_benchmark_datasets_processed", envir = env)
  saved <- if (had) get("fmri_benchmark_datasets_processed", envir = env) else NULL
  if (had) rm("fmri_benchmark_datasets_processed", envir = env)
  withr::defer({
    if (exists("fmri_benchmark_datasets_processed", envir = env)) {
      rm("fmri_benchmark_datasets_processed", envir = env)
    }
    if (had) assign("fmri_benchmark_datasets_processed", saved, envir = env)
  })

  ds <- suppressWarnings(load_benchmark_dataset("BM_Canonical_HighSNR"))
  expect_true(is.matrix(ds$Y_noisy))
  expect_warning(
    {
      rm("fmri_benchmark_datasets_processed", envir = env)
      load_benchmark_dataset("BM_Canonical_HighSNR")
    },
    "fmrihrf"
  )
})
