# Performance pass 2, P2-4: the timecourse / residual-ACF plots read only the
# voxels they draw. When that column read failed, the helpers silently fell
# back to loading the whole dataset; now the error surfaces instead.

.p2_plot_frames <- function() {
  set.seed(8)
  n <- 40
  etab <- data.frame(onset = c(3, 13, 23, 33), cond = factor(c("a", "b", "a", "b")),
                     run = 1L)
  Y <- matrix(rnorm(n * 7), n, 7)
  mf <- matrix_frame(Y, TR = 2, run_length = n, event_table = etab)
  dims <- c(3, 2, 2)
  arr <- array(rnorm(prod(dims) * n), c(dims, n))
  vec <- neuroim2::NeuroVec(arr, neuroim2::NeuroSpace(dim = c(dims, n)))
  m <- array(TRUE, dims)
  m[c(2, 7)] <- FALSE
  mask <- neuroim2::LogicalNeuroVol(m, neuroim2::NeuroSpace(dim = dims))
  nf <- neurovec_frame(scans = list(vec), mask = mask, TR = 2, event_table = etab)
  list(matrix = mf, volume = nf)
}

test_that("P2-4: voxel readers return the requested columns of the data", {
  frames <- .p2_plot_frames()
  for (ds in frames) {
    full <- as.matrix(fmrireg:::.dset_data_matrix(ds))
    for (v in c(1L, 4L, ncol(full))) {
      expect_identical(fmrireg:::.voxel_series(ds, v), as.numeric(full[, v]))
    }
    vv <- c(5L, 2L, 3L)
    expect_equal(unname(fmrireg:::.voxel_matrix(ds, vv)), unname(full[, vv, drop = FALSE]),
                 tolerance = 0)
  }
})

test_that("P2-4 guard: a failed voxel read is an error, not a full-data load", {
  frames <- .p2_plot_frames()
  ds <- frames$volume
  full_reads <- 0L
  real_dm <- get(".dset_data_matrix", envir = asNamespace("fmrireg"))
  real_read <- get(".dset_read", envir = asNamespace("fmrireg"))
  # Simulate a backend whose column-subset read fails: every read of fewer
  # than all features errors; reads of the whole dataset are counted.
  testthat::local_mocked_bindings(
    .dset_data_matrix = function(x, ...) {
      if (identical(x, ds)) {
        full_reads <<- full_reads + 1L
        return(real_dm(x, ...))
      }
      stop("simulated subset read failure")
    },
    .dset_read = function(x, rows = NULL, cols = NULL) {
      if (!is.null(cols)) stop("simulated subset read failure")
      full_reads <<- full_reads + 1L
      real_read(x, rows = rows, cols = cols)
    }
  )
  expect_error(fmrireg:::.voxel_series(ds, 2L), "simulated subset read failure")
  expect_error(fmrireg:::.voxel_matrix(ds, c(1L, 3L)), "simulated subset read failure")
  expect_identical(full_reads, 0L)
})
