# Review 2026-10, stream C.
# C6: hrf_smoothing_kernel() normalised by K / diag(K) (row scaling: not
#     symmetric), its "cosine" variant normalised the trial columns of X rather
#     than the scans, and buffer_scans = 0 gave NaN (a scan with no signal).

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

