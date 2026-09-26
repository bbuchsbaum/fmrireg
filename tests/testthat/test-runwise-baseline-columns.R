# Runwise fits pool betas by *global* design column. Each run's design holds
# the event columns plus only that run's baseline columns, so run-local column
# j is not global column j. Before this was fixed, a two-run fit stored only
# the run-local width (6 of 10 columns) and coef(include_baseline = TRUE)
# labelled the pooled local drift/intercept columns with the first global
# names (e.g. "base_bs1_block_2" on what was run 1's intercept).

.runwise_fixture <- function(baseline = NULL, noise_sd = 1e-6, scale = 100, seed = 11) {
  set.seed(seed)
  ev <- data.frame(
    onset = rep(seq(10, 90, by = 20), 2),
    cond = factor(rep(c("A", "B"), 5)),
    run = rep(1:2, each = 5)
  )
  mk <- function(Y) matrix_frame(Y, TR = 2, run_length = c(60, 60), event_table = ev)
  args <- list(onset ~ hrf(cond), block = ~run,
               control = fmri_lm_control(estimation = estimation_spec("runwise")))
  if (!is.null(baseline)) args$baseline_model <- baseline
  probe <- do.call(fmri_lm, c(args, list(dataset = mk(matrix(rnorm(360), 120, 3)))))
  X <- as.matrix(design_matrix(probe$model))
  B <- matrix(rnorm(3 * ncol(X), sd = scale), 3, ncol(X), dimnames = list(NULL, colnames(X)))
  Y <- X %*% t(B) + matrix(rnorm(120 * 3, sd = noise_sd), 120, 3)
  fit <- do.call(fmri_lm, c(args, list(dataset = mk(Y))))
  list(fit = fit, X = X, B = B)
}

test_that("runwise include_baseline labels and recovers per-run baseline columns", {
  s <- .runwise_fixture()
  expect_identical(attr(s$fit, "strategy"), "runwise")
  bb <- coef(s$fit, include_baseline = TRUE)
  expect_identical(dim(bb), c(3L, ncol(s$X)))
  expect_identical(colnames(bb), colnames(s$X))
  expect_true(all(c("base_bs1_block_2", "constant_1", "constant_2") %in% colnames(bb)))
  # Every column, including each run's own drift and intercept, is the known beta.
  expect_equal(unname(bb), unname(s$B), tolerance = 1e-6)
  expect_equal(bb[, "constant_2"], unname(s$B[, "constant_2"]), tolerance = 1e-6)
})

test_that("runwise pooled columns shared by all runs recover the shared beta", {
  sf <- fmrihrf::sampling_frame(blocklens = c(60, 60), TR = 2)
  bm <- baseline_model(basis = "poly", degree = 1, sframe = sf, intercept = "global")
  s <- .runwise_fixture(baseline = bm, seed = 12)
  bb <- coef(s$fit, include_baseline = TRUE)
  expect_identical(colnames(bb), colnames(s$X))
  # The global intercept is estimated in both runs and pooled across them.
  shared <- grep("global", colnames(s$X), value = TRUE)
  expect_length(shared, 1L)
  expect_equal(unname(bb), unname(s$B), tolerance = 1e-6)
})

test_that("runwise estimates and standard errors stay finite at tiny noise", {
  for (sd in c(1e-6, 1e-10, 0)) {
    s <- .runwise_fixture(noise_sd = sd, seed = 13)
    b <- coef(s$fit)
    expect_true(all(is.finite(b)), info = paste("noise sd", sd))
    expect_equal(unname(b), unname(s$B[, 1:2]), tolerance = 1e-6)
    se <- as.matrix(standard_error(s$fit))
    expect_true(all(is.finite(se)), info = paste("noise sd", sd))
    if (sd > 0) expect_true(all(se > 0), info = paste("noise sd", sd))
  }
})

test_that("do_fixef takes the zero-variance limit instead of returning NaN", {
  beta <- cbind(c(1, 2), c(3, 4))
  se <- cbind(c(0, 1), c(1, 1))
  out <- fmrireg:::do_fixef(se, beta, "inv_var")
  # Row 1: run 1 has infinite precision and gets all the weight.
  expect_equal(out$estimate[1], 1)
  expect_equal(out$se[1], 0)
  # Row 2: ordinary inverse-variance pooling.
  expect_equal(out$estimate[2], 3)
  expect_equal(out$se[2], sqrt(0.5))
  both_zero <- fmrireg:::do_fixef(cbind(0, 0), cbind(2, 4), "inv_var")
  expect_equal(both_zero$estimate, 3)
})

test_that("the lean RSS path is accurate for near-exact fits", {
  set.seed(14)
  X <- cbind(1, rnorm(50))
  e <- rnorm(50, sd = 1e-7)
  Y <- X %*% c(1e3, -2e3) + e
  proj <- fmrireg:::.fast_preproject(X)
  lean <- fmrireg:::solve_glm_core(fmrireg:::glm_context(X = X, Y = Y, proj = proj))
  exact <- sum(stats::lm.fit(X, Y)$residuals^2)
  expect_gt(lean$rss, 0)
  expect_equal(lean$rss, exact, tolerance = 1e-6)
})
