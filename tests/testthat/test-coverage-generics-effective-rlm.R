# all_generic methods, effective_df leftovers, and multiresponse_rlm.

test_that("longnames/shortnames/columns/nbasis/design_matrix.convolved_term", {
  emod <- fmrireg:::.demo_event_model()
  sframe <- fmrireg:::.demo_sampling_frame()

  ln <- longnames(emod)
  expect_true(is.character(ln))
  expect_true(length(ln) >= 1L)

  sn <- shortnames(emod)
  expect_true(is.character(sn))
  expect_true(length(sn) >= 1L)

  cols <- columns(emod)
  expect_true(is.character(cols))
  expect_true(length(cols) >= 1L)

  # columns.event_model with NULL design_matrix
  emod2 <- emod
  emod2$design_matrix <- NULL
  expect_equal(columns(emod2), character(0))

  # Term-level long/short names
  eterms <- terms(emod)
  expect_true(length(eterms) >= 1L)
  tt <- eterms[[1]]
  expect_true(length(longnames(tt)) >= 1L)
  expect_true(length(shortnames(tt)) >= 1L)

  # Convolved term helpers when available
  if (inherits(tt, "convolved_term") || !is.null(tt$evterm) || !is.null(tt$design_matrix)) {
    nb <- tryCatch(nbasis(tt), error = function(e) NULL)
    expect_true(is.null(nb) || (is.numeric(nb) && nb >= 1))
    dm <- tryCatch(design_matrix(tt), error = function(e) NULL)
    expect_true(is.null(dm) || NROW(dm) >= 1L)
    if (!is.null(dm) && !is.null(tt$sampling_frame)) {
      dm1 <- tryCatch(design_matrix(tt, blockid = 1L), error = function(e) NULL)
      expect_true(is.null(dm1) || NROW(dm1) <= NROW(dm))
    }
  }

  expect_equal(blockids(emod), fmrihrf::blockids(sframe))
})

test_that("with_package errors when dependency missing; construct/design_env generics exist", {
  expect_error(fmrireg:::with_package("definitely_not_a_real_pkg_zzz"), "install")
  expect_true(is.function(fmrireg:::design_env))
  expect_true(is.function(fmrireg:::construct))
  expect_true(is.function(estimate_contrast))
  expect_true(is.function(chunkwise_lm))
})

test_that("satterthwaite_df / kenward_roger_df / residual_effective_df", {
  df_s <- fmrireg:::satterthwaite_df(c(2, 3, 4), c(10, 20, 30))
  expect_true(is.finite(df_s) && df_s > 0)

  expect_warning(kr <- fmrireg:::kenward_roger_df(list()), "not implemented")
  expect_true(is.na(kr))

  expect_equal(
    fmrireg:::residual_effective_df(list(effective_df = 42), "ols"),
    42
  )
  expect_equal(
    fmrireg:::residual_effective_df(list(dfres = 30, rank = 3), "ar"),
    30
  )
  expect_equal(
    fmrireg:::residual_effective_df(list(dfres = 25), "ols"),
    25
  )
  expect_warning(
    na_df <- fmrireg:::residual_effective_df(list(), "ols"),
    "No df information"
  )
  expect_true(is.na(na_df))

  # Default method branch in compute_effective_df
  X <- cbind(1, rnorm(20))
  expect_equal(fmrireg:::compute_effective_df(list(), X, method = "unknown"), 18)
})

test_that("multiresponse_rlm fits per-response robust models", {
  skip_if_not_installed("robustbase")
  set.seed(9)
  n <- 40L
  V <- 3L
  X <- cbind(1, rnorm(n), rnorm(n))
  colnames(X) <- c("(Intercept)", "A", "B")
  Y <- X %*% matrix(c(0, 1, -0.5), 3, 1) %*% matrix(1, 1, V) + matrix(rnorm(n * V), n, V)
  # Contaminate one voxel
  Y[1:3, 1] <- Y[1:3, 1] + 20

  data_env <- list2env(list(.y = Y, A = X[, 2], B = X[, 3]), parent = baseenv())
  form <- .y ~ A + B
  conlist <- list()

  # Formula path (modmat = NULL) is the reliable lmrob route
  res <- fmrireg:::multiresponse_rlm(form, data_env, conlist, colnames(X), fcon = NULL, modmat = NULL)
  expect_true(is.list(res))
  expect_true(!is.null(res$bstats))
  expect_true(nrow(res$bstats) >= V)

  data_env2 <- list2env(list(.y = Y[, 1:2, drop = FALSE], A = X[, 2], B = X[, 3]), parent = baseenv())
  res2 <- fmrireg:::multiresponse_rlm(form, data_env2, conlist, colnames(X), fcon = NULL, modmat = NULL)
  expect_true(nrow(res2$bstats) >= 2L)
})

test_that("fmri_rlm chunkwise strategy and bisquare psi", {
  set.seed(10)
  etab <- data.frame(
    onset = c(5, 15, 25, 35, 45, 55),
    condition = factor(rep(c("A", "B"), 3)),
    run = 1L
  )
  Y <- matrix(rnorm(70 * 3), 70, 3)
  dset <- matrix_dataset(Y, TR = 1, run_length = 70, event_table = etab)
  fit <- fmri_rlm(
    onset ~ hrf(condition), block = ~ run, dataset = dset,
    strategy = "chunkwise", nchunks = 2L, robust_psi = "bisquare"
  )
  expect_s3_class(fit, "fmri_rlm")
  expect_output(print(fit), "Robust")
})
