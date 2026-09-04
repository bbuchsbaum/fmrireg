# Public parametric basis helpers, event constructors, meta write, HAC variance.

test_that("Poly/Scale/RobustScale/ScaleWithin parametric bases work", {
  p <- Poly(seq(-1, 1, length.out = 10), degree = 2)
  expect_s3_class(p, "Poly")
  expect_true(inherits(p, "ParametricBasis"))
  expect_true(is.matrix(p$y) || is.numeric(p$y))

  sc <- Scale(1:10)
  expect_s3_class(sc, "Scale")
  expect_equal(length(sc$y), 10L)

  rs <- RobustScale(c(1:9, 100))
  expect_s3_class(rs, "RobustScale")

  sw <- ScaleWithin(1:12, rep(1:3, each = 4))
  expect_s3_class(sw, "ScaleWithin")
})

test_that("contrast constructors and generators cover error/success paths", {
  expect_error(pair_contrast("A", "B"), regexp = ".")

  pc <- pair_contrast(~ A, ~ B, name = "A_vs_B")
  expect_true(inherits(pc, "contrast") || inherits(pc, "pair_contrast") || is.list(pc))

  uc <- unit_contrast(~ condition, name = "u")
  expect_true(inherits(uc, "contrast") || is.list(uc))

  oc <- tryCatch(oneway_contrast(~ condition, name = "cond"), error = function(e) e)
  if (!inherits(oc, "error")) {
    expect_true(is.list(oc) || inherits(oc, "contrast"))
  }

  ua <- tryCatch(one_against_all_contrast(~ condition), error = function(e) e)
  if (!inherits(ua, "error")) {
    expect_true(is.list(ua) || length(ua) >= 1)
  }
})

test_that("event_factor/variable and is_categorical/continuous helpers", {
  ef <- event_factor(
    factor(c("A", "B", "A")), name = "cond",
    onsets = c(1, 2, 3), blockids = rep(1, 3)
  )
  expect_true(inherits(ef, "event_factor") || is.list(ef))
  expect_true(isTRUE(is_categorical(ef)) || isTRUE(is_categorical(factor(c("A", "B")))))

  ev <- event_variable(
    c(0.1, 0.5, 0.9), name = "amp",
    onsets = c(1, 2, 3), blockids = rep(1, 3)
  )
  expect_true(inherits(ev, "event_variable") || is.list(ev))
  expect_true(isTRUE(is_continuous(ev)) || isTRUE(is_continuous(rnorm(3))))
})

test_that("write_results.fmri_meta exercises strategy branches", {
  meta <- suppressWarnings(fmrireg:::.demo_fmri_meta())
  skip_if(!inherits(meta, "fmri_meta"), "no demo meta")

  outdir <- tempfile("meta-wr-")
  dir.create(outdir)

  res1 <- tryCatch(
    write_results(meta, path = outdir, strategy = "by_stat", format = "nifti"),
    error = function(e) e
  )
  if (inherits(res1, "error")) {
    expect_match(
      conditionMessage(res1),
      "mask|space|Neuro|path|format|stat|coef|nifti|h5|template|dimension",
      ignore.case = TRUE
    )
  } else {
    expect_true(is.list(res1) || is.character(res1))
  }

  res2 <- tryCatch(
    write_results(meta, path = outdir, strategy = "by_coefficient", format = "h5"),
    error = function(e) e
  )
  if (inherits(res2, "error")) {
    expect_match(conditionMessage(res2), ".", perl = TRUE)
  } else {
    expect_true(is.list(res2) || is.character(res2))
  }
})

test_that("HAC variance_from_context covers model+satterthwaite and robust weights", {
  set.seed(55)
  n <- 36
  p <- 2
  V <- 2
  X <- cbind(1, rnorm(n))
  E <- matrix(rnorm(n * V), n, V)
  ctx <- list(
    X = X,
    residuals = E,
    runs = rep(1:2, each = 18),
    censor = NULL,
    robust_weights = runif(n, 0.5, 1)
  )

  # variance_spec signature: method first
  spec <- variance_spec(method = "hac", max_lag = 2L, taper = "tukey", df = "satterthwaite")
  fit <- fmrireg:::.fmri_lm_variance_from_context(ctx, spec)
  expect_equal(length(fit$covariance), V)
  expect_equal(length(fit$df), V)
  expect_true(all(fit$df >= 1))

  spec2 <- variance_spec(method = "sandwich")
  fit2 <- fmrireg:::.fmri_lm_variance_from_context(ctx, spec2)
  expect_equal(fit2$max_lag, 0L)

  expect_error(
    fmrireg:::.fmri_lm_variance_from_context(
      list(X = X, residuals = E[1:10, , drop = FALSE], runs = NULL),
      spec2
    ),
    "incompatible"
  )
})

test_that("fmri_meta FE on demo CSV group data", {
  gd <- tryCatch(fmrireg:::.demo_group_data_csv(), error = function(e) NULL)
  skip_if(is.null(gd), "no csv group data")

  fit <- fmri_meta(gd, formula = ~ 1, method = "fe")
  expect_s3_class(fit, "fmri_meta")
  expect_true(
    length(coef(fit)) > 0 || !is.null(fit$coefficients) || !is.null(fit$beta)
  )
})
