context("Robust argument passing")

# simple synthetic dataset
etab <- data.frame(onset=c(1,10), repnum=factor(c("A","B")), run=c(1,1))
Y <- matrix(rnorm(20*5), 20, 5)
dset <- matrix_dataset(Y, TR=1, run_length=20, event_table=etab)

# accept valid robust arguments
test_that("fmri_lm accepts robust arguments", {
  # Robust fitting requires the fast engine (the runwise slow path rejects it).
  fit <- fmri_lm(onset ~ hrf(repnum), block = ~ run, dataset = dset,
                 robust = TRUE, robust_psi = "huber", robust_max_iter = 2,
                 robust_scale_scope = "run", use_fast_path = TRUE)
  cfg <- attr(fit, "executed_config")
  expect_equal(cfg$robust$type, "huber")
  expect_equal(cfg$robust$max_iter, 2)
  expect_equal(cfg$robust$scale_scope, "run")

  # robust_psi selects the psi function on its own. This previously had no
  # effect at all -- `robust` defaulted to FALSE and always claimed `type`
  # first -- so the argument was silently dead.
  fit2 <- fmri_lm(onset ~ hrf(repnum), block = ~ run, dataset = dset,
                  robust_psi = "bisquare", robust_max_iter = 1,
                  robust_scale_scope = "global", use_fast_path = TRUE)
  cfg2 <- attr(fit2, "executed_config")
  expect_equal(cfg2$robust$type, "bisquare")
  expect_equal(cfg2$robust$scale_scope, "global")
})

test_that("robust = TRUE and robust_psi = 'bisquare' is a conflict, not a silent win", {
  # This combination used to be accepted and quietly fit Huber. It is a genuine
  # contradiction: TRUE canonicalises to huber.
  expect_error(
    fmri_lm(onset ~ hrf(repnum), block = ~ run, dataset = dset,
            robust = TRUE, robust_psi = "bisquare", use_fast_path = TRUE),
    "Conflicting robust settings"
  )
})

# validation of arguments
test_that("robust argument validation works", {
  # Test invalid robust_psi through robust_options
  expect_error(
    fmri_lm(onset ~ hrf(repnum), block = ~ run, dataset = dset,
            robust_options = list(type = "bogus")),
    "robust_options\\$type.*TRUE.*FALSE.*huber.*bisquare"
  )

  expect_error(
    fmri_lm(onset ~ hrf(repnum), block = ~ run, dataset = dset,
            robust = TRUE, robust_max_iter = 0),
    "max_iter.*>= 1"
  )
})
