# Regression tests for #247: supplying fmri_lm_control() for an unrelated
# section (robust, noise, ...) must not silently switch fmri_lm() from its
# runwise meta-estimation default to a single joint fit.

.est_default_dset <- function() {
  set.seed(1)
  events <- data.frame(
    onset = rep(seq(5, 89, by = 12), 2),
    cond = factor(rep(c("a", "b"), 8)),
    run = rep(1:2, each = 8)
  )
  matrix_frame(
    matrix(rnorm(120 * 4), 120, 4), TR = 2,
    run_length = c(60, 60), event_table = events
  )
}

test_that("fmri_lm_control() without `estimation` defaults to runwise_meta", {
  expect_identical(fmri_lm_control()$estimation$scope, "runwise_meta")
  expect_identical(
    fmri_lm_control(robust = robust_spec())$estimation$scope,
    "runwise_meta"
  )
  expect_identical(
    fmri_lm_control(noise = noise_spec("ar1"))$estimation$scope,
    "runwise_meta"
  )
  # An explicit choice is always honoured.
  expect_identical(
    fmri_lm_control(estimation = estimation_spec("joint"))$estimation$scope,
    "joint"
  )
  expect_false(fmrireg:::.fmri_lm_estimation_supplied(fmri_lm_control()))
  expect_true(fmrireg:::.fmri_lm_estimation_supplied(
    fmri_lm_control(estimation = estimation_spec("runwise_meta"))
  ))
})

test_that("a control object for robust fitting keeps fmri_lm()'s runwise strategy", {
  dset <- .est_default_dset()

  fit_default <- fmri_lm(onset ~ hrf(cond), block = ~run, dataset = dset)
  fit_robust <- fmri_lm(
    onset ~ hrf(cond), block = ~run, dataset = dset,
    control = fmri_lm_control(robust = robust_spec()),
    compute = compute_spec(voxel_chunks = 1)
  )
  fit_ar <- fmri_lm(
    onset ~ hrf(cond), block = ~run, dataset = dset,
    control = fmri_lm_control(noise = noise_spec("ar1"))
  )

  expect_identical(fit_default$strategy, "runwise")
  expect_identical(fit_robust$strategy, "runwise")
  expect_identical(fit_ar$strategy, "runwise")
})

test_that("an explicit iid control reproduces the default fit exactly", {
  dset <- .est_default_dset()
  fit_default <- fmri_lm(onset ~ hrf(cond), block = ~run, dataset = dset)
  fit_ctrl <- fmri_lm(
    onset ~ hrf(cond), block = ~run, dataset = dset,
    control = fmri_lm_control(noise = noise_spec("iid"))
  )
  expect_identical(fit_ctrl$strategy, fit_default$strategy)
  expect_equal(coef(fit_ctrl), coef(fit_default))
  expect_equal(standard_error(fit_ctrl), standard_error(fit_default))
})

test_that("joint estimation still runs when explicitly requested", {
  dset <- .est_default_dset()
  fit_joint <- fmri_lm(
    onset ~ hrf(cond), block = ~run, dataset = dset,
    control = fmri_lm_control(
      robust = robust_spec(),
      estimation = estimation_spec("joint")
    ),
    compute = compute_spec(voxel_chunks = 1)
  )
  expect_identical(fit_joint$strategy, "chunkwise")
})

test_that("MA noise no longer needs an explicit runwise estimation spec", {
  ctrl <- fmri_lm_control(noise = noise_spec("ar1", q = 1L))
  expect_identical(ctrl$estimation$scope, "runwise_meta")
  expect_error(
    fmri_lm_control(
      noise = noise_spec("ar1", q = 1L),
      estimation = estimation_spec("joint")
    ),
    "MA terms currently require"
  )
})

test_that("engines keep their joint default unless estimation is explicit", {
  engine_name <- "estimation_default_probe_engine"
  captured <- NULL
  register_engine(
    engine_name,
    fit = function(model, dataset, args, cfg) {
      captured <<- cfg$estimation$scope
      stop("probe-engine-stop", call. = FALSE)
    }
  )
  on.exit(rm(list = engine_name, envir = fmrireg:::.fmrireg_engine_registry),
          add = TRUE)
  dset <- .est_default_dset()
  run_probe <- function(...) {
    expect_error(
      fmri_lm(onset ~ hrf(cond), block = ~run, dataset = dset,
              engine = engine_name, ...),
      "probe-engine-stop"
    )
    captured
  }

  expect_identical(run_probe(), "joint")
  expect_identical(
    run_probe(control = fmri_lm_control(noise = noise_spec("iid"))),
    "joint"
  )
  expect_identical(
    run_probe(control = fmri_lm_control(
      estimation = estimation_spec("runwise_meta")
    )),
    "runwise_meta"
  )
})

test_that("joint-only inference resolves an omitted scope to joint, with a message", {
  for (v in list(variance_spec("sandwich"), variance_spec("hac"),
                 variance_spec(df = "satterthwaite"))) {
    expect_message(
      ctrl <- fmri_lm_control(variance = v),
      "require joint estimation"
    )
    expect_identical(ctrl$estimation$scope, "joint")
  }
  # An explicit scope is never overridden and produces no message.
  expect_no_message(
    ctrl <- fmri_lm_control(
      variance = variance_spec("sandwich"),
      estimation = estimation_spec("joint")
    )
  )
  expect_identical(ctrl$estimation$scope, "joint")
  expect_no_message(fmri_lm_control(variance = variance_spec("model")))
})

test_that("engine dispatch resolves an omitted scope to joint", {
  dset <- .est_default_dset()
  model <- create_fmri_model(onset ~ hrf(cond), block = ~run, dataset = dset)
  # rrr_gls supports only the joint scope; an unspecified control must reach it
  # as joint rather than as the built-in runwise default.
  fit <- fmrireg:::.fmri_lm_dispatch_engine(
    model = model, dataset = dset, engine = "rrr_gls",
    cfg = fmri_lm_control(),
    engine_args = list(rank = 1L, bootstrap_n = 5L)
  )
  expect_s3_class(fit, "fmri_lm")
  expect_error(
    fmrireg:::.fmri_lm_dispatch_engine(
      model = model, dataset = dset, engine = "rrr_gls",
      cfg = fmri_lm_control(estimation = estimation_spec("runwise_meta")),
      engine_args = list(rank = 1L, bootstrap_n = 5L)
    ),
    "does not support estimation scope 'runwise_meta'"
  )
})
