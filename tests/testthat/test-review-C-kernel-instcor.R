# Review 2026-10, stream C.
# C6: hrf_smoothing_kernel() normalised by K / diag(K) (row scaling: not
#     symmetric), its "cosine" variant normalised the trial columns of X rather
#     than the scans, and buffer_scans = 0 gave NaN (a scan with no signal).
# C9: instantaneous_correlation_rcpp(fill = "na") let one missing sample
#     poison every later value, and N == 0 read y[-1].

test_that("C6: normalised kernels are symmetric correlation matrices", {
  form <- onset ~ trialwise(basis = "gaussian")
  for (method in c("gram", "cosine")) {
    for (buf in c(0L, 3L)) {
      K <- hrf_smoothing_kernel(30, TR = 1.5, form = form, method = method,
                                buffer_scans = buf, normalise = TRUE)
      info <- paste(method, "buffer", buf)
      expect_true(all(is.finite(K)), info = info)
      expect_equal(K, t(K), tolerance = 1e-12, info = info)
      expect_equal(diag(K), rep(1, 30), tolerance = 1e-12, info = info)
      expect_true(all(abs(K) <= 1 + 1e-12), info = info)
    }
  }
})

test_that("C6: gram kernel equals X X' with symmetric scaling; cosine is the row cosine", {
  form <- onset ~ trialwise(basis = "gaussian")
  len <- 25L; buf <- 2L; TR <- 2
  sf <- fmrihrf::sampling_frame(len + 2 * buf, TR)
  dfx <- data.frame(onset = fmrihrf::samples(sf), block = 1L)
  X <- as.matrix(design_matrix(event_model(form, data = dfx, block = ~ block,
                                           sampling_frame = sf)))
  keep <- (buf + 1):(buf + len)

  G <- X %*% t(X)
  Kg <- hrf_smoothing_kernel(len, TR = TR, form = form, buffer_scans = buf,
                             normalise = FALSE, method = "gram")
  expect_equal(Kg, G[keep, keep], tolerance = 1e-12)

  d <- sqrt(diag(G))
  ref <- G / outer(d, d)
  Kn <- hrf_smoothing_kernel(len, TR = TR, form = form, buffer_scans = buf,
                             normalise = TRUE, method = "gram")
  expect_equal(Kn, ref[keep, keep], tolerance = 1e-12)

  # Cosine similarity between scans (rows of X).
  Kc <- hrf_smoothing_kernel(len, TR = TR, form = form, buffer_scans = buf,
                             normalise = FALSE, method = "cosine")
  expect_equal(Kc, ref[keep, keep], tolerance = 1e-12)
})

test_that("C9: fill = 'na' leaves later correlations finite and matches the 'zero' path inside the data", {
  set.seed(4)
  x <- rnorm(80)
  y <- 0.6 * x + rnorm(80, sd = 0.5)
  # offset = 2: y[k - 2] is outside the series for k = 0, 1
  r_na <- instantaneous_correlation_rcpp(x, y, eta = 0.2, offset = 2L,
                                         warmup = 0L, fill = "na")
  expect_length(r_na, 80L)
  expect_true(all(is.na(r_na[1:2])))
  expect_true(all(is.finite(r_na[4:80])))

  # The 'na' path skips the missing samples: it equals running the recursion
  # on the aligned pairs (x[k], y[k - 2]) for k >= 2.
  r_ref <- instantaneous_correlation_rcpp(x[3:80], y[1:78], eta = 0.2,
                                          offset = 0L, warmup = 0L, fill = "zero")
  expect_equal(r_na[3:80], r_ref, tolerance = 1e-12)

  # Negative offset: missing values at the end only.
  r_neg <- instantaneous_correlation_rcpp(x, y, eta = 0.2, offset = -3L,
                                          warmup = 0L, fill = "na")
  expect_true(all(is.na(r_neg[78:80])))
  expect_true(all(is.finite(r_neg[2:77])))
})

test_that("C9: empty input returns an empty vector", {
  expect_identical(instantaneous_correlation_rcpp(numeric(0), numeric(0), eta = 0.1),
                   numeric(0))
  expect_identical(instantaneous_correlation_rcpp(numeric(0), numeric(0), eta = 0.1,
                                                  fill = "last"),
                   numeric(0))
})
