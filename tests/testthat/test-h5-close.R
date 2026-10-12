# hdf5r's H5File$close_all() starts with a full gc() (to finalise forgotten
# object handles). When only the file itself is open -- writers that close
# what they create -- the file can be closed directly. .h5_close() does that,
# and falls back to close_all() whenever anything else is still open.

.h5_tmp <- function(leave_dataset_open) {
  skip_if_not_installed("hdf5r")
  f <- tempfile(fileext = ".h5")
  h5 <- hdf5r::H5File$new(f, mode = "w")
  d <- h5$create_dataset("x", robj = as.numeric(1:5))
  if (!leave_dataset_open) d$close()
  list(file = f, h5 = h5, dataset = d)
}

.h5_reopen_ok <- function(f) {
  h <- hdf5r::H5File$new(f, mode = "r")
  on.exit(h$close_all(), add = TRUE)
  as.numeric(h[["x"]]$read())
}

test_that(".h5_close closes the file directly when nothing else is open (no gc)", {
  tmp <- .h5_tmp(leave_dataset_open = FALSE)
  on.exit(unlink(tmp$file), add = TRUE)
  closed_all <- 0L
  trace("gc", quote(closed_all <<- closed_all + 1L), print = FALSE, where = baseenv())
  on.exit(untrace("gc", where = baseenv()), add = TRUE)
  fmrireg:::.h5_close(tmp$h5)
  untrace("gc", where = baseenv())
  expect_identical(closed_all, 0L)
  expect_false(tmp$h5$is_valid)
  expect_equal(.h5_reopen_ok(tmp$file), 1:5)
})

test_that(".h5_close falls back to close_all() when objects are still open", {
  tmp <- .h5_tmp(leave_dataset_open = TRUE)
  on.exit(unlink(tmp$file), add = TRUE)
  fmrireg:::.h5_close(tmp$h5)
  expect_false(tmp$h5$is_valid)
  expect_false(tmp$dataset$is_valid)
  expect_equal(.h5_reopen_ok(tmp$file), 1:5)
})

test_that("write_results closes every H5 file it writes through .h5_close", {
  skip_if_not_installed("fmristore")
  calls <- 0L
  real <- fmrireg:::.h5_close
  testthat::local_mocked_bindings(
    .h5_close = function(h5) { calls <<- calls + 1L; real(h5) },
    .package = "fmrireg")
  skip_if_not_installed("neuroim2")
  set.seed(123)
  dims <- c(3, 3, 2); nt <- 50
  scans <- lapply(1:2, function(run) {
    neuroim2::NeuroVec(array(rnorm(prod(dims) * nt), c(dims, nt)),
                       neuroim2::NeuroSpace(dim = c(dims, nt)))
  })
  mask <- neuroim2::LogicalNeuroVol(array(TRUE, dims), neuroim2::NeuroSpace(dim = dims))
  ev <- data.frame(onset = rep(c(5, 15, 25, 35, 45), 2),
                   condition = factor(rep(c("A", "B", "A", "B", "A"), 2)),
                   run = rep(1:2, each = 5))
  dset <- neurovec_frame(scans = scans, mask = mask, TR = 1.5, event_table = ev)
  con <- contrast_set(pair_contrast(~ condition == "A", ~ condition == "B",
                                    name = "A_vs_B"))
  fit <- fmri_lm(onset ~ hrf(condition, contrasts = con), block = ~ run,
                 dataset = dset, durations = 0, progress = FALSE)
  out <- withr::local_tempdir()
  files <- suppressMessages(write_results(fit, path = out, subject = "01",
                                          task = "test", space = "MNI152NLin2009cAsym",
                                          format = "h5",
                                          strategy = "by_contrast"))
  h5_files <- list.files(out, "\\.h5$", recursive = TRUE)
  expect_gt(length(h5_files), 0L)
  expect_identical(calls, length(h5_files))
})
