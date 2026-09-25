# term_matrices() used to request every run by id, and fmridesign returns a
# global intercept once per requested run. The fit then carried duplicated
# intercept columns and baseline_indices that ran past the end of the design.

global_intercept_fixture <- function(intercept = "global") {
  withr::local_seed(1)
  nscan <- 60L
  sf <- fmrihrf::sampling_frame(blocklens = rep(nscan, 3), TR = 2)
  ev <- do.call(rbind, lapply(1:3, function(r) {
    data.frame(run = r, onset = seq(6, 100, by = 12),
               condition = factor(rep(c("a", "b"), length.out = 8)))
  }))
  Y <- matrix(stats::rnorm(3 * nscan * 40), 3 * nscan) + 50
  ds <- matrix_frame(Y, TR = 2, run_length = rep(nscan, 3), event_table = ev)
  bm <- fmridesign::baseline_model(basis = "poly", degree = 2, sframe = sf,
                                   intercept = intercept)
  em <- fmridesign::event_model(onset ~ hrf(condition), data = ev, block = ~run,
                                sampling_frame = sf)
  list(Y = Y, ds = ds, bm = bm, em = em, sf = sf, ev = ev)
}

test_that("term_matrices() keeps a global intercept as one column", {
  fx <- global_intercept_fixture()
  fm <- fmri_model(fx$em, fx$bm, fx$ds)
  n_design <- ncol(as.matrix(design_matrix(fm)))
  n_base <- ncol(as.matrix(design_matrix(fx$bm)))

  expect_no_message(tm <- term_matrices(fm))
  expect_equal(sum(vapply(tm, ncol, integer(1))), n_design)
  bi <- attr(tm, "baseline_term_indices")
  expect_length(bi, n_base)
  expect_equal(max(bi), n_design)
  expect_false(anyDuplicated(attr(tm, "varnames")) > 0)
})

test_that("a single run still gets that run's columns", {
  fx <- global_intercept_fixture()
  fm <- fmri_model(fx$em, fx$bm, fx$ds)
  tm2 <- term_matrices(fm, blocknum = 2)
  expect_equal(nrow(tm2[[1]]), 60L)
  expect_equal(length(attr(tm2, "baseline_term_indices")),
               ncol(as.matrix(design_matrix(fx$bm, blockid = 2))))
})

test_that("fits with a global intercept index the design correctly", {
  fx <- global_intercept_fixture()
  for (scope in c("joint", "runwise_meta")) {
    fit <- expect_no_warning(fmri_lm(
      onset ~ hrf(condition), block = ~run, dataset = fx$ds, baseline_model = fx$bm,
      control = fmri_lm_control(estimation = estimation_spec(scope))
    ))
    n_design <- ncol(as.matrix(design_matrix(fit$model)))
    expect_equal(max(fit$result$baseline_indices), n_design)
    expect_length(fit$result$baseline_indices, ncol(as.matrix(design_matrix(fx$bm))))
  }

  # joint estimates agree with ordinary least squares on the same design
  fit <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = fx$ds,
                 baseline_model = fx$bm,
                 control = fmri_lm_control(estimation = estimation_spec("joint")))
  X <- as.matrix(design_matrix(fit$model))
  ref <- stats::lm.fit(X, fx$Y[, 1])$coefficients[1:2]
  est <- as.matrix(stats(fit))[1, ] * as.matrix(standard_error(fit))[1, ]
  expect_equal(unname(est), unname(ref), tolerance = 1e-6)
})

test_that("runwise intercepts are unchanged", {
  fx <- global_intercept_fixture(intercept = "runwise")
  fm <- fmri_model(fx$em, fx$bm, fx$ds)
  tm <- term_matrices(fm)
  expect_length(attr(tm, "baseline_term_indices"), ncol(as.matrix(design_matrix(fx$bm))))
  expect_equal(sum(vapply(tm, ncol, integer(1))), ncol(as.matrix(design_matrix(fm))))
})
