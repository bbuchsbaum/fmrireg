# Twelfth wave: estimate_hrf prepare/partial-out and tidy curve/voxel filters.

test_that(".prepare_hrf_estimation builds designs with and without fixed", {
  set.seed(261)
  n <- 100L
  events <- data.frame(
    onset = seq(8, n - 20, by = 12),
    condition = factor(rep(c("A", "B"), length.out = 7)),
    nuisance = rnorm(7),
    run = 1L
  )
  Y <- matrix(rnorm(n * 2), n, 2)
  colnames(Y) <- c("v1", "v2")
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = events)
  spec <- fmrireg:::.new_hrf_basis_spec("tent", k = 5L, span = 24)

  prep <- fmrireg:::.prepare_hrf_estimation(
    form = onset ~ hrf(condition),
    fixed = NULL,
    block = ~ run,
    dataset = dset,
    basemod = NULL,
    basis_spec = spec
  )
  expect_true(is.matrix(prep$event_design))
  expect_true(is.matrix(prep$nuisance_design))
  expect_true(length(prep$curve_indices) >= 1L)
  expect_s3_class(prep$event_model, "event_model")
  expect_null(prep$fixed_model)

  prep_fix <- fmrireg:::.prepare_hrf_estimation(
    form = onset ~ hrf(condition),
    fixed = onset ~ hrf(nuisance),
    block = ~ run,
    dataset = dset,
    basemod = baseline_model("constant", sframe = dset$sampling_frame),
    basis_spec = spec
  )
  expect_true(!is.null(prep_fix$fixed_model))
  expect_true(ncol(prep_fix$nuisance_design) > ncol(prep$nuisance_design))

  expect_error(
    fmrireg:::.prepare_hrf_estimation(
      form = onset ~ hrf(condition),
      fixed = "not-a-formula",
      block = ~ run,
      dataset = dset,
      basemod = NULL,
      basis_spec = spec
    ),
    "fixed must be NULL or an event-model formula"
  )
})

test_that("tidy.fmri_hrf_estimate filters curves/voxels and validates bounds", {
  set.seed(262)
  n <- 110L
  events <- data.frame(
    onset = seq(8, n - 24, by = 14),
    condition = factor(rep(c("A", "B"), length.out = 6)),
    run = 1L
  )
  Y <- matrix(rnorm(n * 3), n, 3)
  colnames(Y) <- paste0("vox", 1:3)
  dset <- matrix_dataset(Y, TR = 1, run_length = n, event_table = events)
  fit <- estimate_hrf(
    onset ~ hrf(condition),
    block = ~ run,
    dataset = dset,
    basis = "tent",
    k = 5L,
    lambda = 0.1,
    progress = FALSE
  )

  tid_all <- tidy(fit)
  expect_true(nrow(tid_all) > 0)
  expect_true(all(c("time", "curve", "voxel", "estimate", "std.error") %in% names(tid_all)))

  tid_curve <- tidy(fit, curve = fit$curves[[1]])
  expect_equal(unique(as.character(tid_curve$curve)), fit$curves[[1]])

  tid_vox <- tidy(fit, voxel = 2L)
  expect_equal(unique(as.character(tid_vox$voxel)), fit$voxels[[2]])

  expect_error(tidy(fit, curve = 99L), "out-of-range")
  expect_error(tidy(fit, voxel = "missing_vox"), "out-of-range")
})

test_that("validate_estimate_hrf_formula and order helpers reject bad input", {
  expect_error(fmrireg:::.validate_estimate_hrf_formula(~ hrf(x)), "two-sided")
  expect_error(
    fmrireg:::.validate_estimate_hrf_formula(onset ~ trialwise()),
    "trialwise"
  )
  expect_error(
    fmrireg:::.validate_estimate_hrf_formula(onset ~ poly(x)),
    "hrf\\(\\)"
  )
  expect_true(fmrireg:::.validate_estimate_hrf_formula(onset ~ hrf(condition)))

  events <- data.frame(onset = 1:4, condition = factor(c("A", "B", "A", "B")), run = 1L)
  expect_error(
    fmrireg:::.order_hrf_event_data(events, wrong ~ hrf(condition), ~ run),
    "onset column"
  )
  expect_error(
    fmrireg:::.order_hrf_event_data(events, onset ~ hrf(condition), run ~ x),
    "one-sided"
  )
  ordered <- fmrireg:::.order_hrf_event_data(
    events[c(4, 2, 3, 1), ], onset ~ hrf(condition), ~ run
  )
  expect_equal(ordered$onset, sort(events$onset))
})
