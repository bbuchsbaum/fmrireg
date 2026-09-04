# Coverage for fmrilm AR routing helpers and get_formula.

test_that("get_formula.fmri_model builds response formula from terms", {
  fx <- {
    etab <- data.frame(
      onset = c(5, 15, 25, 35),
      condition = factor(c("A", "B", "A", "B")),
      run = 1L
    )
    Y <- matrix(rnorm(50 * 3), 50, 3)
    dset <- matrix_dataset(Y, TR = 1, run_length = 50, event_table = etab)
    emod <- event_model(
      onset ~ hrf(condition), data = etab, block = ~ run,
      sampling_frame = dset$sampling_frame
    )
    bmod <- baseline_model(basis = "poly", degree = 1, sframe = dset$sampling_frame)
    fmri_model(emod, bmod, dset)
  }
  form <- fmrireg:::get_formula.fmri_model(fx)
  expect_true(inherits(form, "formula"))
  expect_match(deparse(form), "\\.y ~")
  expect_error(fmrireg:::get_formula.fmri_model(list()), "fmri_model")
})

test_that("estimate_ar_parameters_routed covers edge orders and fallbacks", {
  set.seed(3)
  n <- 80
  x <- rnorm(n)
  # order <= 0
  expect_equal(fmrireg:::.estimate_ar_parameters_routed(x, 0L), numeric(0))
  expect_equal(fmrireg:::.estimate_ar_parameters_routed(x, NULL), numeric(0))

  # short series
  short <- fmrireg:::.estimate_ar_parameters_routed(rnorm(2), 2L)
  expect_equal(length(short), 2)

  # AR1 / AR2 / arp paths with OLS residuals vs design
  X <- cbind(1, rnorm(n))
  y <- as.numeric(X %*% c(0.5, 1) + arima.sim(list(ar = 0.4), n = n))
  resid <- y - X %*% solve(crossprod(X), crossprod(X, y))

  phi1 <- fmrireg:::.estimate_ar_parameters_routed(
    resid, 1L, design = X, ar_options = list(cor_struct = "ar1")
  )
  expect_equal(length(phi1), 1)

  phi2 <- fmrireg:::.estimate_ar_parameters_routed(
    resid, 2L, design = X, ar_options = list(exact_first = TRUE)
  )
  expect_equal(length(phi2), 2)

  phi3 <- fmrireg:::.estimate_ar_parameters_routed(
    resid, 3L, design = X
  )
  expect_equal(length(phi3), 3)

  # Multi-run
  run_idx <- list(1:40, 41:80)
  phi_runs <- fmrireg:::.estimate_ar_parameters_routed(
    resid, 1L, run_indices = run_idx, design = X,
    ar_options = list(global = FALSE)
  )
  expect_true(is.numeric(phi_runs) || is.list(phi_runs))
})
