# Review 2026-10, stream D, item P-D2: lazy chunking.
# Results must not depend on the number of voxel chunks, and the engines must
# read each voxel's data a bounded number of times (not the whole dataset
# once per helper).

.ck_fixture <- function(nrun = 2L, n = 80L, V = 21L, seed = 5L) {
  set.seed(seed)
  ons <- seq(8, 2 * n - 24, by = 14)
  ev <- do.call(rbind, lapply(seq_len(nrun), function(r) {
    data.frame(onset = ons,
               condition = factor(rep(c("A", "B"), length.out = length(ons))),
               run = r)
  }))
  Y <- matrix(0, n * nrun, V)
  for (r in seq_len(nrun)) {
    for (v in seq_len(V)) {
      Y[(r - 1L) * n + seq_len(n), v] <- as.numeric(arima.sim(list(ar = 0.4), n)) + 30
    }
  }
  Y[c(5, 90), ] <- Y[c(5, 90), ] + 12
  matrix_frame(Y, TR = 2, run_length = rep(n, nrun), event_table = ev)
}

.ck_fit <- function(dset, noise = "iid", robust = "none", chunks = 1L,
                    variance = NULL, scope = "joint") {
  fmri_lm(onset ~ hrf(condition), block = ~run, dataset = dset, durations = 0,
          control = fmri_lm_control(estimation = estimation_spec(scope),
                                    noise = noise_spec(noise),
                                    robust = robust_spec(robust),
                                    variance = variance),
          compute = compute_spec(voxel_chunks = chunks))
}

.ck_stats <- function(fit) {
  d <- fit$result$betas$data[[1]]
  list(est = d$estimate[[1]], se = d$se[[1]], stat = d$stat[[1]],
       sigma = fit$result$sigma, df = fit$result$df$inference)
}

test_that("joint fits are identical across 1, 3 and 7 voxel chunks", {
  dset <- .ck_fixture()
  configs <- list(c("iid", "none"), c("ar1", "none"), c("iid", "bisquare"),
                  c("ar1", "huber"))
  for (cfg in configs) {
    ref <- .ck_stats(.ck_fit(dset, cfg[1], cfg[2], chunks = 1L))
    for (k in c(3L, 7L)) {
      got <- .ck_stats(.ck_fit(dset, cfg[1], cfg[2], chunks = k))
      expect_equal(got, ref, tolerance = 1e-10,
                   info = paste(cfg[1], cfg[2], k))
    }
  }
})

test_that("each voxel is read a bounded number of times", {
  dset <- .ck_fixture()
  total <- prod(dim(as.matrix(fmridataset::collect_assay(dset))))
  cells <- 0
  orig_read <- fmrireg:::.dset_read
  testthat::local_mocked_bindings(
    .dset_read = function(x, rows = NULL, cols = NULL) {
      out <- orig_read(x, rows, cols)
      cells <<- cells + length(out)
      out
    }
  )
  passes <- function(...) {
    cells <<- 0
    fit <- .ck_fit(dset, ...)
    cells / total
  }
  # OLS: one pass over the data, however it is chunked.
  expect_equal(passes("iid", "none", chunks = 1L), 1)
  expect_equal(passes("iid", "none", chunks = 7L), 1)
  # AR / robust: one pass to estimate the noise model, one to fit the chunks.
  expect_equal(passes("ar1", "none", chunks = 7L), 2)
  expect_equal(passes("iid", "bisquare", chunks = 7L), 2)
  expect_equal(passes("ar1", "huber", chunks = 7L), 2)
  # Runwise: each run is read once.
  expect_equal(passes("ar1", "none", scope = "runwise_meta"), 1)
  # HAC inference reuses the chunk residuals; no extra pass.
  expect_equal(passes("ar1", "none", chunks = 7L,
                      variance = variance_spec("hac", max_lag = 2)), 2)
})

test_that("chunk descriptors are lazy and carry no response matrix", {
  dset <- .ck_fixture()
  chunks <- fmrireg:::.dset_chunks(dset, nchunks = 3)
  expect_true(all(vapply(chunks, function(ch) is.null(ch[["data"]]), logical(1))))
  # The data are still available on request, and correct.
  Y <- as.matrix(fmridataset::collect_assay(dset))
  expect_equal(unname(chunks[[2]]$data), unname(Y[, chunks[[2]]$voxel_ind]))
  run_chunks <- fmrireg:::.dset_run_chunks(dset)
  expect_true(all(vapply(run_chunks, function(ch) is.null(ch[["data"]]), logical(1))))

  model <- create_fmri_model(onset ~ hrf(condition), block = ~run,
                             dataset = dset, durations = 0)
  cfg <- fmri_lm_control(estimation = estimation_spec("joint"),
                         noise = noise_spec("ar1"))
  precomp <- fmrireg:::prepare_chunkwise_matrices(model, dset, cfg)
  expect_false(any(vapply(precomp$run_info, function(ri) "Y_orig" %in% names(ri),
                          logical(1))))
})
