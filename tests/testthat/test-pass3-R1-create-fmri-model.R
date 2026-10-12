# Pass 3, R1-3: create_fmri_model() had three definitions (fmri_model.R,
# fmri_model_utils.R, fmrilm.R). Collation made fmrilm.R's run while the
# exported help described the fmri_model_utils.R one ("poly", degree 3,
# global intercept, block given as a vector). One definition remains; these
# tests pin the defaults it actually applies.

.p3_cfm_frame <- function(run_length = c(450L, 450L)) {
  set.seed(1)
  nrun <- length(run_length)
  ev <- data.frame(onset = rep(c(10, 40, 70), nrun),
                   condition = factor(rep(c("a", "b", "a"), nrun)),
                   run = rep(seq_len(nrun), each = 3))
  Y <- matrix(rnorm(sum(run_length) * 2), sum(run_length), 2)
  matrix_frame(Y, TR = 1, run_length = run_length, event_table = ev)
}

test_that("create_fmri_model has exactly one definition", {
  src <- list.files(test_path("..", "..", "R"), pattern = "[.]R$", full.names = TRUE)
  skip_if(length(src) == 0L, "package sources not available")
  hits <- unlist(lapply(src, function(f) {
    grep("^create_fmri_model\\s*<-", readLines(f, warn = FALSE), value = TRUE)
  }))
  expect_length(hits, 1L)
})

test_that("default baseline: B-spline drift, degree max(ceiling(median run length / 100), 3), runwise intercept", {
  dset <- .p3_cfm_frame(c(450L, 450L))
  m <- create_fmri_model(onset ~ hrf(condition), block = ~run, dataset = dset)
  spec <- m$baseline_model$drift_spec
  expect_identical(spec$basis, "bs")
  expect_equal(spec$degree, 5)                 # ceiling(450 / 100) = 5
  expect_identical(spec$intercept, "runwise")
  nm <- colnames(design_matrix(m$baseline_model))
  expect_identical(sum(grepl("^constant_", nm)), 2L)   # one intercept per run

  # Short runs: degree floors at 3.
  m_short <- create_fmri_model(onset ~ hrf(condition), block = ~run,
                               dataset = .p3_cfm_frame(c(120L, 120L)))
  expect_equal(m_short$baseline_model$drift_spec$degree, 3)
  expect_identical(m_short$baseline_model$drift_spec$basis, "bs")
})

test_that("a supplied baseline model is used as is; block must be a formula", {
  dset <- .p3_cfm_frame(c(120L, 120L))
  bm <- baseline_model(basis = "poly", degree = 2, sframe = frame_sframe(dset))
  m <- create_fmri_model(onset ~ hrf(condition), block = ~run,
                         baseline_model = bm, dataset = dset)
  expect_identical(m$baseline_model, bm)
  expect_error(create_fmri_model(onset ~ hrf(condition), block = "run",
                                 dataset = dset), "formula")
  expect_error(create_fmri_model(onset ~ hrf(condition), block = ~run,
                                 baseline_model = list(), dataset = dset),
               "baseline_model")
})

test_that("durations: a scalar applies to every event", {
  dset <- .p3_cfm_frame(c(120L, 120L))
  m0 <- create_fmri_model(onset ~ hrf(condition), block = ~run, dataset = dset)
  m2 <- create_fmri_model(onset ~ hrf(condition), block = ~run, dataset = dset,
                          durations = 2)
  m2v <- create_fmri_model(onset ~ hrf(condition), block = ~run, dataset = dset,
                           durations = rep(2, 6))
  X0 <- as.matrix(design_matrix(m0$event_model))
  X2 <- as.matrix(design_matrix(m2$event_model))
  expect_false(isTRUE(all.equal(X0, X2)))
  expect_equal(X2, as.matrix(design_matrix(m2v$event_model)))
})
