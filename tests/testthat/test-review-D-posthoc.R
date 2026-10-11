# Review 2026-10, stream D, item C5: the stored covariance must match the
# stored scale, so post-hoc contrasts (fit_contrasts.fmri_lm) reproduce the
# contrasts computed during the fit for AR and robust fits.

.ph_fixture <- function(nrun = 2L, n = 100L, V = 6L, seed = 2L) {
  set.seed(seed)
  ons <- seq(10, 2 * n - 30, by = 16)
  ev <- do.call(rbind, lapply(seq_len(nrun), function(r) {
    data.frame(onset = ons,
               condition = factor(rep(c("A", "B"), length.out = length(ons))),
               run = r)
  }))
  Y <- matrix(0, n * nrun, V)
  for (r in seq_len(nrun)) {
    for (v in seq_len(V)) {
      Y[(r - 1L) * n + seq_len(n), v] <- as.numeric(arima.sim(list(ar = 0.5), n)) + 100
    }
  }
  Y[c(7, 60), ] <- Y[c(7, 60), ] + 15
  matrix_frame(Y, TR = 2, run_length = rep(n, nrun), event_table = ev)
}

.ph_check <- function(dset, scope, noise, robust, backend, chunks) {
  fit <- fmri_lm(
    onset ~ hrf(condition, contrasts = pair_contrast(~condition == "A",
                                                     ~condition == "B",
                                                     name = "AvB")),
    block = ~run, dataset = dset, durations = 0,
    control = fmri_lm_control(estimation = estimation_spec(scope),
                              noise = noise_spec(noise),
                              robust = robust_spec(robust)),
    compute = compute_spec(voxel_chunks = chunks, backend = backend)
  )
  infit <- fit$result$contrasts
  infit <- infit$data[[which(infit$name == "AvB")]]
  ph <- fit_contrasts(fit, list(AvB = pair_contrast(~condition == "A",
                                                    ~condition == "B",
                                                    name = "AvB")))[[1]]
  label <- paste(scope, noise, robust, backend, chunks)
  expect_equal(ph$estimate, infit$estimate, tolerance = 1e-8, info = label)
  expect_equal(ph$se, infit$se, tolerance = 1e-8, info = label)
  expect_equal(ph$stat, infit$stat, tolerance = 1e-8, info = label)
}

test_that("post-hoc contrasts equal in-fit contrasts for joint AR fits", {
  dset <- .ph_fixture()
  for (backend in c("matrix", "reference")) {
    for (chunks in c(1L, 3L)) {
      .ph_check(dset, "joint", "ar1", "none", backend, chunks)
    }
  }
})

test_that("post-hoc contrasts equal in-fit contrasts for joint robust fits", {
  dset <- .ph_fixture()
  for (backend in c("matrix", "reference")) {
    for (chunks in c(1L, 3L)) {
      .ph_check(dset, "joint", "iid", "huber", backend, chunks)
      .ph_check(dset, "joint", "ar1", "bisquare", backend, chunks)
    }
  }
})

test_that("post-hoc contrasts equal in-fit contrasts for one-run runwise AR fits", {
  dset <- .ph_fixture(nrun = 1L)
  .ph_check(dset, "runwise_meta", "ar1", "none", "matrix", 1L)
  .ph_check(dset, "runwise_meta", "ar1", "huber", "matrix", 1L)
})
