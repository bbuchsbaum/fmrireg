# Regression tests for issue #183: a rank-deficient design (e.g. duplicate /
# collinear nuisance regressors, as produced by fMRIPrep confound tables) must
# NOT silently NA-poison contrasts that do not load on the aliased columns.
# The aliased coefficients come back NA in betamat; because 0 * NA = NA in R,
# the slow-path contrast estimate t(cmat) %*% betamat used to propagate NA into
# every contrast, even estimable ones. The fast and slow paths must agree.

make_rank_deficient_dataset <- function() {
  set.seed(7)
  ev <- data.frame(
    onset = sort(runif(24, 0, 230)), run = 1,
    cond = factor(rep(c("A", "B"), 12)), m = rnorm(24)
  )
  ev$m2 <- ev$m # exact duplicate modulator -> aliased / rank-deficient design
  fmridataset::matrix_dataset(
    matrix(rnorm(120 * 5), 120, 5),
    TR = 2, run_length = 120, event_table = ev
  )
}

test_that("estimable contrast is not NA-poisoned by an aliased nuisance column (#183)", {
  dset <- make_rank_deficient_dataset()
  con <- contrast_set(pair_contrast(~ cond == "A", ~ cond == "B", name = "A_vs_B"))
  form <- onset ~ hrf(cond, contrasts = con) + hrf(m) + hrf(m2)

  fit_slow <- suppressWarnings(fmri_lm(form,
    block = ~run, dataset = dset,
    strategy = "runwise", use_fast_path = FALSE, progress = FALSE
  ))
  fit_fast <- suppressWarnings(fmri_lm(form,
    block = ~run, dataset = dset,
    strategy = "runwise", use_fast_path = TRUE, progress = FALSE
  ))

  v_slow <- as.numeric(unlist(stats(fit_slow, type = "contrasts")[["A_vs_B"]]))
  v_fast <- as.numeric(unlist(stats(fit_fast, type = "contrasts")[["A_vs_B"]]))

  # A_vs_B loads only on `cond`, never on the aliased modulator -> must be finite.
  expect_true(all(is.finite(v_slow)),
    info = "estimable contrast NA-poisoned by aliased nuisance column (slow path)"
  )
  expect_true(all(is.finite(v_fast)))
  # The slow path must agree with the (already-correct) fast path.
  expect_equal(v_slow, v_fast, tolerance = 1e-6)
})
