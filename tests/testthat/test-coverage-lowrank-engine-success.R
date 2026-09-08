# Successful .run_lowrank_engine paths across sketch methods + plugin dispatch.

make_lowrank_fixture <- function(n = 80L, V = 10L, seed = 19) {
  set.seed(seed)
  etab <- data.frame(
    onset = c(5, 15, 25, 35, 45, 55),
    condition = factor(c("A", "B", "A", "B", "A", "B")),
    run = 1L
  )
  Y <- matrix(rnorm(n * V), n, V)
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = etab)
  emod <- event_model(
    onset ~ hrf(condition), data = etab, block = ~ run,
    sampling_frame = dset$sampling_frame
  )
  bmod <- baseline_model(basis = "poly", degree = 1, sframe = dset$sampling_frame)
  list(
    model = fmri_model(emod, bmod, dset),
    dataset = dset,
    cfg = fmri_lm_control()
  )
}

expect_lowrank_fit <- function(fit) {
  expect_s3_class(fit, "fmri_lm")
  expect_true(!is.null(fit$result$betas))
  expect_true(!is.null(fit$result$sigma) || !is.null(fit$sigma2))
  expect_equal(attr(fit, "strategy"), "sketch")
}

test_that(".run_lowrank_engine succeeds for gaussian/srht/countsketch (IID)", {
  fx <- make_lowrank_fixture()

  fit_g <- fmrireg:::.run_lowrank_engine(
    fx$model, fx$dataset,
    lowrank = list(time_sketch = list(method = "gaussian", m = 24L)),
    cfg = fx$cfg
  )
  expect_lowrank_fit(fit_g)

  fit_s <- fmrireg:::.run_lowrank_engine(
    fx$model, fx$dataset,
    lowrank = list(time_sketch = list(method = "srht", m = 24L)),
    cfg = fx$cfg
  )
  expect_lowrank_fit(fit_s)

  fit_c <- fmrireg:::.run_lowrank_engine(
    fx$model, fx$dataset,
    lowrank = list(time_sketch = list(method = "countsketch", m = 24L)),
    cfg = fx$cfg
  )
  expect_lowrank_fit(fit_c)
})

test_that(".run_lowrank_engine ihs sketch path succeeds without landmarks", {
  fx <- make_lowrank_fixture(n = 70L, V = 8L, seed = 21)
  fit <- fmrireg:::.run_lowrank_engine(
    fx$model, fx$dataset,
    lowrank = list(time_sketch = list(method = "ihs", m = 20L, iters = 2L)),
    cfg = fx$cfg
  )
  expect_lowrank_fit(fit)
})

test_that(".run_lowrank_engine global AR whitening path (ar1) succeeds", {
  skip_if_not_installed("fmriAR")
  fx <- make_lowrank_fixture(n = 90L, V = 8L, seed = 23)
  cfg <- fmri_lm_control(ar_options = list(struct = "ar1"))
  fit <- fmrireg:::.run_lowrank_engine(
    fx$model, fx$dataset,
    lowrank = list(time_sketch = list(method = "gaussian", m = 24L)),
    cfg = cfg
  )
  expect_lowrank_fit(fit)
  expect_true(!is.null(fit$ar_coef) || !is.null(fit$result$ar_coef))
})

test_that(".fit_lowrank_engine_plugin and fmri_lm_lowrank_dispatch work", {
  fx <- make_lowrank_fixture(n = 60L, V = 6L, seed = 25)
  low <- list(time_sketch = list(method = "gaussian", m = 16L))

  fit_plugin <- fmrireg:::.fit_lowrank_engine_plugin(
    fx$model, fx$dataset,
    args = list(lowrank = low),
    cfg = fx$cfg
  )
  expect_lowrank_fit(fit_plugin)

  expect_null(fmrireg:::fmri_lm_lowrank_dispatch(fx$model, fx$dataset, engine = NULL))

  fit_disp <- fmrireg:::fmri_lm_lowrank_dispatch(
    fx$model, fx$dataset,
    engine = "sketch",
    lowrank = low,
    cfg = fx$cfg
  )
  expect_lowrank_fit(fit_disp)

  # Formula path through dispatch
  etab <- fx$dataset$event_table
  fit_form <- fmrireg:::fmri_lm_lowrank_dispatch(
    onset ~ hrf(condition),
    dataset = fx$dataset,
    engine = "latent_sketch",
    lowrank = low,
    block = ~ run,
    cfg = fx$cfg
  )
  expect_lowrank_fit(fit_form)
})
