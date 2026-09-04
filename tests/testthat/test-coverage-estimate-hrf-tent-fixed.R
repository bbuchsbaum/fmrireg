# estimate_hrf success with tent/bspline + fixed formula and validation gates.

make_hrf_cov_fixture <- function(n = 140L, n_voxels = 3L, seed = 181L) {
  set.seed(seed)
  onsets <- seq(8, n - 30, by = 12)
  events <- data.frame(
    onset = onsets,
    condition = factor(rep(c("A", "B"), length.out = length(onsets))),
    nuisance = sin(seq_along(onsets) * 0.7),
    run = 1L
  )
  Y <- matrix(rnorm(n * n_voxels), n, n_voxels)
  colnames(Y) <- paste0("vox", seq_len(n_voxels))
  list(
    dataset = matrix_dataset(Y, TR = 1, run_length = n, event_table = events),
    events = events
  )
}

test_that("estimate_hrf tent and bspline succeed with fixed nuisance formula", {
  fx <- make_hrf_cov_fixture()

  tent <- estimate_hrf(
    onset ~ hrf(condition),
    fixed = onset ~ hrf(nuisance),
    block = ~ run,
    dataset = fx$dataset,
    basis = "tent",
    k = 6L,
    lambda = 0,
    rsam = seq(0, 18, by = 1),
    progress = FALSE
  )
  expect_s3_class(tent, "fmri_hrf_estimate")
  expect_identical(tent$basis, "tent")
  expect_equal(dim(tent$estimate)[2:3], c(2L, 3L))
  expect_true(all(is.finite(tent$estimate)))
  expect_true(all(is.finite(tent$std.error)))

  bsp <- estimate_hrf(
    onset ~ hrf(condition),
    fixed = onset ~ hrf(nuisance),
    block = ~ run,
    dataset = fx$dataset,
    basis = "bspline",
    k = 6L,
    lambda = 0.5,
    progress = FALSE
  )
  expect_identical(bsp$basis, "bspline")
  expect_equal(length(bsp$curves), 2L)
  expect_gt(bsp$df.residual, 0)

  # Deprecated bs= maps tent correctly
  expect_warning(
    dep <- estimate_hrf(
      onset ~ hrf(condition),
      block = ~ run,
      dataset = fx$dataset,
      bs = "tent",
      k = 5L,
      lambda = 0.1,
      progress = FALSE
    ),
    "bs is deprecated"
  )
  expect_identical(dep$basis, "tent")

  # fx deprecated forces unpenalized when lambda missing
  expect_warning(
    unpen <- estimate_hrf(
      onset ~ hrf(condition),
      block = ~ run,
      dataset = fx$dataset,
      fx = TRUE,
      k = 5L,
      progress = FALSE
    ),
    "fx is deprecated"
  )
  expect_equal(unpen$lambda, 0)
})

test_that("estimate_hrf validation and helper boundary errors", {
  fx <- make_hrf_cov_fixture(n = 100L)

  expect_error(
    estimate_hrf(~ hrf(condition), block = ~ run, dataset = fx$dataset, progress = FALSE),
    "two-sided"
  )
  expect_error(
    estimate_hrf(onset ~ factor(condition), block = ~ run, dataset = fx$dataset, progress = FALSE),
    "hrf\\(\\)"
  )
  expect_error(
    estimate_hrf(onset ~ hrf(condition), block = ~ run, dataset = list(), progress = FALSE),
    "fmri_dataset"
  )
  expect_error(
    estimate_hrf(
      onset ~ hrf(condition), block = ~ run, dataset = fx$dataset,
      basis = "bspline", k = 3L, progress = FALSE
    ),
    "at least 4"
  )
  expect_error(
    estimate_hrf(
      onset ~ hrf(condition), block = ~ run, dataset = fx$dataset,
      basis = "tent", k = 1L, progress = FALSE
    ),
    "at least 2"
  )
  expect_error(
    estimate_hrf(
      onset ~ hrf(condition), block = ~ run, dataset = fx$dataset,
      rsam = c(1, 2, 3), progress = FALSE
    ),
    "start at zero"
  )
  expect_error(
    estimate_hrf(
      onset ~ hrf(condition), block = ~ run, dataset = fx$dataset,
      rsam = c(0, 2, 1), progress = FALSE
    ),
    "strictly increasing"
  )
  expect_error(
    estimate_hrf(
      onset ~ hrf(condition), block = ~ run, dataset = fx$dataset,
      fx = "nope", progress = FALSE
    ),
    "TRUE, or FALSE"
  )
  expect_error(
    estimate_hrf(
      onset ~ hrf(condition), block = ~ run, dataset = fx$dataset,
      ci_level = 1.5, progress = FALSE
    ),
    "ci_level"
  )
  expect_error(
    estimate_hrf(
      onset ~ hrf(condition), block = ~ run, dataset = fx$dataset,
      lambda = "aic", progress = FALSE
    ),
    "gcv"
  )

  # Helper: order/block validation
  expect_error(
    fmrireg:::.order_hrf_event_data(fx$events, bad ~ hrf(condition), ~ run),
    "onset column"
  )
  expect_error(
    fmrireg:::.order_hrf_event_data(fx$events, onset ~ hrf(condition), run ~ 1),
    "one-sided"
  )

  # Lambda selector numeric path
  XtX <- diag(4)
  XtY <- matrix(rnorm(4 * 2), 4, 2)
  pen <- diag(4)
  sel <- fmrireg:::.select_hrf_lambda(
    XtX, XtY, response_ss = c(10, 12), penalty = pen,
    n_effective = 80, lambda = 0.25, lambda_grid = 1, progress = FALSE
  )
  expect_equal(sel$lambda, 0.25)

  expect_null(fmrireg:::.solve_hrf_system(-diag(3), matrix(1, 3, 1), diag(3), 0))
})
