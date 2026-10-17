# Simplification 2026-10: checks added while merging the parallel streams.

test_that("a failed AR estimate warns instead of silently dropping whitening", {
  local_mocked_bindings(.estimate_ar_via_fmriAR = function(...) stop("boom"))
  set.seed(1)
  r <- matrix(rnorm(200), 50, 4)
  expect_warning(
    phi <- .estimate_ar_parameters_routed(r, 1L),
    class = "fmrireg_ar_estimation_failed"
  )
  expect_identical(phi, 0)
})

test_that(".fmrireg_lapply gives the same results sequentially and through a future plan", {
  skip_if_not_installed("future")
  f <- function(i) sum(seq_len(i)) + stats::runif(1) * 0
  seq_out <- fmrireg:::.fmrireg_lapply(1:6, f)
  old_plan <- future::plan(future::multisession, workers = 2)
  on.exit(future::plan(old_plan), add = TRUE)
  par_out <- fmrireg:::.fmrireg_lapply(1:6, f, parallel = TRUE)
  expect_identical(par_out, seq_out)
  expect_identical(unlist(seq_out), cumsum(1:6) * 1.0)
})

test_that("plotly, bslib and furrr are not hard dependencies", {
  imports <- packageDescription("fmrireg", fields = "Imports")
  if (is.na(imports)) imports <- read.dcf(test_path("../../DESCRIPTION"), fields = "Imports")[1, 1]
  for (pkg in c("plotly", "bslib", "furrr", "RcppParallel")) {
    expect_false(grepl(paste0("\\b", pkg, "\\b"), imports), info = pkg)
  }
})
