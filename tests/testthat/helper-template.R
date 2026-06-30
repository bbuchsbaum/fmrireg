# Shared fixtures for the model-template / fan-out tests (M2, M3, M5, ...).

# A small in-memory matrix_dataset with a 2-run, 2-condition design.
make_test_matrix_dataset <- function(nvox = 4L, runs = c(40L, 40L), TR = 2) {
  set.seed(42)
  total <- sum(runs)
  datamat <- matrix(rnorm(total * nvox), total, nvox)
  ev <- do.call(rbind, lapply(seq_along(runs), function(r) {
    on <- seq(4, runs[r] * TR - 12, by = 8)
    data.frame(
      onset = on,
      condition = factor(rep(c("A", "B"), length.out = length(on))),
      run = r
    )
  }))
  matrix_dataset(datamat, TR = TR, run_length = runs, event_table = ev)
}

# A minimal fitted fmri_lm on the matrix_dataset fixture.
make_test_fit <- function(ds = make_test_matrix_dataset()) {
  fmri_lm(onset ~ hrf(condition), block = ~ run, dataset = ds, strategy = "runwise")
}
