# Regression tests: term_matrices() must not duplicate a global intercept
# when the whole design (all runs) is requested.

.tm_fixture <- function(runs = c(60L, 60L, 60L), intercept = "global",
                        nvox = 7L, TR = 2, seed = 11) {
  set.seed(seed)
  nr <- length(runs)
  ev <- do.call(rbind, lapply(seq_len(nr), function(r) {
    on <- seq(6, runs[r] * TR - 20, by = 12)
    data.frame(
      onset = on,
      condition = factor(rep(c("A", "B"), length.out = length(on))),
      run = r
    )
  }))
  sf <- fmrihrf::sampling_frame(blocklens = runs, TR = TR)
  emod <- fmridesign::event_model(onset ~ hrf(condition), data = ev,
                                  block = ~run, sampling_frame = sf)
  bmod <- fmridesign::baseline_model(basis = "poly", degree = 1, sframe = sf,
                                     intercept = intercept)
  # True design: event columns + baseline columns, each built without blockid
  X <- cbind(as.matrix(fmridesign::design_matrix(emod)),
             as.matrix(fmridesign::design_matrix(bmod)))
  beta <- matrix(rnorm(ncol(X) * nvox, sd = 2), ncol(X), nvox)
  beta[seq_len(2), ] <- matrix(c(3, -1.5), 2, nvox)
  Y <- X %*% beta + matrix(rnorm(nrow(X) * nvox, sd = 1e-3), nrow(X), nvox)
  ds <- matrix_frame(Y, TR = TR, run_length = runs, event_table = ev)
  list(ds = ds, bmod = bmod, emod = emod, X = X, beta = beta)
}

.joint_ctrl <- function() {
  fmri_lm_control(estimation = estimation_spec("joint"))
}

# Orient a beta matrix as predictors x voxels regardless of coef() convention.
.as_pred_by_vox <- function(b, npred, nvox) {
  b <- as.matrix(b)
  if (nrow(b) == npred && ncol(b) == nvox) b else t(b)
}

test_that("term_matrices() keeps a single global intercept across runs", {
  fx <- .tm_fixture()
  fm <- fmri_model(fx$emod, fx$bmod, dataset = fx$ds)
  dm <- design_matrix(fm)

  expect_no_message(tm <- term_matrices(fm))
  bidx <- attr(tm, "baseline_term_indices")
  eidx <- attr(tm, "event_term_indices")

  expect_equal(ncol(dm), ncol(fx$X))
  expect_true(all(bidx <= ncol(dm)))
  expect_equal(c(eidx, bidx), seq_len(ncol(dm)))
  expect_equal(sum(vapply(tm, ncol, integer(1))), ncol(dm))
  expect_equal(unname(do.call(cbind, tm)), unname(fx$X))

  vn <- attr(tm, "varnames")
  expect_false(anyDuplicated(vn) > 0)
  expect_equal(sum(vn == "constant_global"), 1L)
  expect_equal(vn[bidx], colnames(fmridesign::design_matrix(fx$bmod)))
  expect_equal(qr(do.call(cbind, tm))$rank, ncol(dm))
})

test_that("design_matrix.fmri_model() with all runs matches the full design", {
  fx <- .tm_fixture()
  fm <- fmri_model(fx$emod, fx$bmod, dataset = fx$ds)
  expect_no_message(dm_all <- design_matrix(fm, blockid = 1:3))
  expect_equal(as.matrix(dm_all), as.matrix(design_matrix(fm)))
})

test_that("term_matrices() on a run subset is unchanged", {
  fx <- .tm_fixture()
  fm <- fmri_model(fx$emod, fx$bmod, dataset = fx$ds)
  tm2 <- term_matrices(fm, 2)
  expect_equal(attr(tm2, "blocknum"), 2)
  expect_equal(nrow(tm2[[1]]), 60L)
  rows <- 61:120
  expect_equal(unname(do.call(cbind, tm2)),
               unname(fx$X[rows, c(1:2, 4, 6)]))
  tm23 <- term_matrices(fm, 2:3)
  expect_equal(nrow(tm23[[1]]), 120L)
})

test_that("joint fit with a global intercept recovers known betas", {
  fx <- .tm_fixture()
  expect_no_warning(expect_no_message(
    fit <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = fx$ds,
                   baseline_model = fx$bmod, control = .joint_ctrl(),
                   compute = compute_spec(voxel_chunks = 1L))
  ))
  nb <- ncol(fx$X)
  expect_true(all(fit$result$baseline_indices <= nb))
  expect_equal(c(fit$result$event_indices, fit$result$baseline_indices),
               seq_len(nb))
  b <- .as_pred_by_vox(coef(fit, include_baseline = TRUE), nb, 7L)
  expect_equal(unname(b), unname(fx$beta), tolerance = 1e-3)
})

test_that("runwise-scope fit with a global intercept recovers event betas", {
  fx <- .tm_fixture()
  expect_no_warning(expect_no_message(
    fit <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = fx$ds,
                   baseline_model = fx$bmod)
  ))
  b <- .as_pred_by_vox(coef(fit), 2L, 7L)
  expect_equal(unname(b), unname(fx$beta[1:2, ]), tolerance = 1e-3)
})

test_that("runwise-intercept multi-run and single-run fits are unchanged", {
  fx <- .tm_fixture(intercept = "runwise")
  fm <- fmri_model(fx$emod, fx$bmod, dataset = fx$ds)
  tm <- term_matrices(fm)
  expect_equal(unname(do.call(cbind, tm)), unname(fx$X))
  expect_equal(attr(tm, "baseline_term_indices"), 3:ncol(fx$X))
  expect_no_warning(
    fit <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = fx$ds,
                   baseline_model = fx$bmod, control = .joint_ctrl(),
                   compute = compute_spec(voxel_chunks = 1L))
  )
  b <- .as_pred_by_vox(coef(fit, include_baseline = TRUE), ncol(fx$X), 7L)
  expect_equal(unname(b), unname(fx$beta), tolerance = 1e-3)

  fx1 <- .tm_fixture(runs = 90L, intercept = "global")
  fm1 <- fmri_model(fx1$emod, fx1$bmod, dataset = fx1$ds)
  expect_no_message(tm1 <- term_matrices(fm1))
  expect_equal(unname(do.call(cbind, tm1)), unname(fx1$X))
  expect_equal(attr(tm1, "baseline_term_indices"), 3:ncol(fx1$X))
  fit1 <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = fx1$ds,
                  baseline_model = fx1$bmod, control = .joint_ctrl(),
                  compute = compute_spec(voxel_chunks = 1L))
  b1 <- .as_pred_by_vox(coef(fit1, include_baseline = TRUE), ncol(fx1$X), 7L)
  expect_equal(unname(b1), unname(fx1$beta), tolerance = 1e-3)
})

test_that("a run without events still counts as a run when all are requested", {
  # Run ids come from the sampling frame. The event model's per-event ids
  # omit a run with no events (run 3 here), so requesting every run
  # explicitly would not be recognised as the whole design and, with an
  # fmridesign that replicates a global intercept per requested run, the
  # intercept would be duplicated again. Current fmridesign no longer
  # replicates it (fmridesign#39), so this guards the invariant rather than
  # reproducing the old failure.
  set.seed(5)
  runs <- c(60L, 60L, 60L)
  ev <- do.call(rbind, lapply(c(1L, 2L), function(r) {
    data.frame(onset = seq(6, 100, by = 12),
               condition = factor(rep(c("A", "B"), length.out = 8)), run = r)
  }))
  sf <- fmrihrf::sampling_frame(blocklens = runs, TR = 2)
  emod <- fmridesign::event_model(onset ~ hrf(condition), data = ev,
                                  block = ~run, sampling_frame = sf)
  bmod <- fmridesign::baseline_model(basis = "poly", degree = 1, sframe = sf,
                                     intercept = "global")
  Y <- matrix(rnorm(sum(runs) * 3), sum(runs), 3)
  ds <- matrix_frame(Y, TR = 2, run_length = runs, event_table = ev)
  fm <- fmri_model(emod, bmod, dataset = ds)
  n_design <- ncol(as.matrix(design_matrix(fm)))

  for (bn in list(NULL, 1:3)) {
    tm <- term_matrices(fm, blocknum = bn)
    expect_equal(sum(vapply(tm, ncol, integer(1))), n_design)
    expect_equal(max(attr(tm, "baseline_term_indices")), n_design)
    expect_false(anyDuplicated(attr(tm, "varnames")) > 0)
  }
  expect_equal(ncol(as.matrix(design_matrix(fm, blockid = 1:3))), n_design)
})
