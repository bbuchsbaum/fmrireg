# Simplification 2026-10, stream A: regression tests for the fitting pipeline.

# B11. The reference backend's joint (chunkwise) AR fit estimated phi inside
# every voxel chunk from at most 64 of its voxels, while the matrix backend
# estimates phi once per run from all voxels. Reference results therefore
# changed with `voxel_chunks`. Both backends now estimate per run on the full
# data and hold phi fixed across chunks.
.b11_fixture <- function(n = 120L, V = 100L, seed = 1101L) {
  set.seed(seed)
  ons <- seq(6, 2 * n - 20, by = 16)
  ev <- do.call(rbind, lapply(1:2, function(r) {
    data.frame(onset = ons,
               condition = factor(rep(c("A", "B"), length.out = length(ons))),
               run = r)
  }))
  Y <- matrix(0, 2 * n, V)
  for (v in seq_len(V)) {
    # Voxel-specific AR(1) strength, so a voxel subsample gives a different
    # shared phi than all voxels.
    phi_v <- 0.1 + 0.6 * (v - 1) / (V - 1)
    for (r in 1:2) {
      Y[(r - 1L) * n + seq_len(n), v] <-
        as.numeric(arima.sim(list(ar = phi_v), n)) + 20 + r
    }
  }
  matrix_frame(Y, TR = 2, run_length = rep(n, 2), event_table = ev)
}

.b11_fit <- function(dset, backend, chunks) {
  fmri_lm(onset ~ hrf(condition, contrasts = pair_contrast(
            ~ condition == "A", ~ condition == "B", name = "AvB")),
          block = ~run, dataset = dset, durations = 0,
          control = fmri_lm_control(estimation = estimation_spec("joint"),
                                    noise = noise_spec("ar1")),
          compute = compute_spec(voxel_chunks = chunks, backend = backend))
}

.b11_get <- function(fit, what) {
  as.matrix(fit$result$betas$data[[1]][[what]][[1]])
}

test_that("B11: reference joint AR fit does not depend on voxel_chunks", {
  dset <- .b11_fixture()
  one <- .b11_fit(dset, "reference", 1L)
  four <- .b11_fit(dset, "reference", 4L)

  expect_equal(.b11_get(four, "estimate"), .b11_get(one, "estimate"), tolerance = 1e-12)
  expect_equal(.b11_get(four, "se"), .b11_get(one, "se"), tolerance = 1e-12)
  expect_equal(four$result$contrasts$data[[1]]$se,
               one$result$contrasts$data[[1]]$se, tolerance = 1e-12)
})

test_that("B11: reference joint AR fit matches the matrix backend", {
  dset <- .b11_fixture()
  ref <- .b11_fit(dset, "reference", 4L)
  mat <- .b11_fit(dset, "matrix", 4L)

  # Same per-run phi ...
  expect_equal(unname(unlist(ref$result$ar_coef[[1]])), unname(unlist(mat$result$ar_coef)),
               tolerance = 1e-10)
  # ... hence the same GLS estimates and standard errors.
  expect_equal(.b11_get(ref, "estimate"), .b11_get(mat, "estimate"), tolerance = 1e-8)
  expect_equal(.b11_get(ref, "se"), .b11_get(mat, "se"), tolerance = 1e-8)
})
