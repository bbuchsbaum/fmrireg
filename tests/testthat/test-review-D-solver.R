# Review 2026-10, stream D: core solver (P-D4), lazy pseudo-inverse (P-D1)
# and rank-deficient residual df (D5).

# ---- P-D4: solve_glm_core uses one QR pass for betas and exact RSS ---------

.rd_near_collinear <- function(n = 200, seed = 11) {
  set.seed(seed)
  x1 <- rnorm(n)
  # Columns 2 and 3 are nearly collinear (condition number ~ 2e7) and the
  # response sits on a large offset, so yTy - b'X'y cancels badly (relative
  # RSS error ~1e-7 with the normal-equation identity).
  X <- cbind(1, x1, x1 + 1e-7 * rnorm(n), rnorm(n))
  colnames(X) <- c("int", "a", "a_near", "b")
  Y <- 1000 + X %*% cbind(c(0, 1, 2, -1), c(0, -1, 0.5, 3)) +
    matrix(rnorm(2 * n, sd = 1), n, 2)
  list(X = X, Y = Y)
}

test_that("solve_glm_core betas and RSS match lm on a near-collinear design", {
  d <- .rd_near_collinear()
  proj <- fmrireg:::.fast_preproject(d$X)
  for (ret_fitted in c(FALSE, TRUE)) {
    res <- fmrireg:::solve_glm_core(
      fmrireg:::glm_context(X = d$X, Y = d$Y, proj = proj),
      return_fitted = ret_fitted
    )
    for (v in 1:2) {
      ref <- lm(d$Y[, v] ~ d$X - 1)
      expect_equal(unname(res$betas[, v]), unname(coef(ref)), tolerance = 1e-6)
      rss_ref <- sum(residuals(ref)^2)
      expect_equal(res$rss[v], rss_ref, tolerance = 1e-10)
      expect_equal(res$sigma2[v], rss_ref / ref$df.residual, tolerance = 1e-10)
    }
  }
})

test_that("solve_glm_core matches lm on an exactly aliased design", {
  set.seed(3)
  n <- 60
  X <- cbind(int = 1, a = rnorm(n), dup = 0, b = rnorm(n))
  X[, "dup"] <- 2 * X[, "a"]
  Y <- matrix(rnorm(n * 3, mean = 50), n, 3)
  proj <- suppressWarnings(fmrireg:::.fast_preproject(X))
  res <- fmrireg:::solve_glm_core(fmrireg:::glm_context(X = X, Y = Y, proj = proj))
  for (v in 1:3) {
    ref <- lm(Y[, v] ~ X - 1)
    cf <- coef(ref)
    expect_equal(unname(res$betas[!is.na(cf), v]), unname(cf[!is.na(cf)]),
                 tolerance = 1e-10)
    expect_equal(res$rss[v], sum(residuals(ref)^2), tolerance = 1e-10)
    expect_equal(res$dfres, ref$df.residual)
  }
})

# ---- P-D1: the n x n pseudo-inverse is formed only on request -------------

test_that(".fast_preproject does not form an n x n identity solve", {
  set.seed(1)
  n <- 150
  X <- cbind(1, matrix(rnorm(n * 4), n, 4))
  big_solves <- 0L
  trace("qr.coef", where = baseenv(), print = FALSE,
        tracer = bquote({
          yy <- as.matrix(y)
          if (nrow(yy) == .(n) && ncol(yy) == .(n)) {
            assign("big_solves", get("big_solves", envir = .(environment())) + 1L,
                   envir = .(environment()))
          }
        }))
  on.exit(untrace("qr.coef", where = baseenv()), add = TRUE)

  proj <- fmrireg:::.fast_preproject(X)
  Y <- matrix(rnorm(n * 3), n, 3)
  res <- fmrireg:::solve_glm_core(fmrireg:::glm_context(X = X, Y = Y, proj = proj))
  expect_equal(big_solves, 0L)

  # The pseudo-inverse is still available on request and is exact.
  expect_equal(c(proj$Pinv), c(MASS::ginv(X)), tolerance = 1e-10)
  expect_equal(dim(proj$Pinv), c(ncol(X), n))
  expect_equal(c(proj$Pinv %*% Y), c(res$betas), tolerance = 1e-10)
})

test_that("lazy Pinv equals the eager one for rank-deficient designs", {
  set.seed(2)
  n <- 40
  X <- cbind(int = 1, a = rnorm(n), dup = 0)
  X[, "dup"] <- X[, "a"]
  proj <- suppressWarnings(fmrireg:::.fast_preproject(X))
  eager <- qr.coef(proj$qr, diag(n))
  eager[!is.finite(eager)] <- 0
  eager[proj$aliased, ] <- 0
  expect_equal(c(proj$Pinv), c(eager), tolerance = 1e-10)
  expect_equal(dim(proj$Pinv), dim(eager))
  expect_identical(rownames(proj$Pinv), colnames(X))
  expect_identical(attr(proj$Pinv, "aliased"), proj$aliased)
})

test_that("voxelwise AR fitting never builds a pseudo-inverse", {
  set.seed(5)
  n <- 60
  ev <- data.frame(onset = seq(4, 100, by = 12),
                   condition = factor(rep(c("A", "B"), length.out = 9)),
                   run = 1)
  Y <- matrix(rnorm(n * 6), n, 6)
  dset <- matrix_frame(Y, TR = 2, run_length = n, event_table = ev)
  calls <- 0L
  orig <- fmrireg:::.preproject_pinv
  testthat::local_mocked_bindings(
    .preproject_pinv = function(...) {
      calls <<- calls + 1L
      orig(...)
    }
  )
  fit <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = dset,
                 durations = 0,
                 control = fmri_lm_control(
                   estimation = estimation_spec("runwise_meta"),
                   noise = noise_spec("ar1", voxelwise = TRUE)
                 ))
  ols <- fmri_lm(onset ~ hrf(condition), block = ~run, dataset = dset,
                 durations = 0,
                 control = fmri_lm_control(estimation = estimation_spec("joint")),
                 compute = compute_spec(voxel_chunks = 2))
  expect_equal(calls, 0L)
})

# ---- D5: residual df on rank-deficient designs -----------------------------

.rd_rankdef_fixture <- function(seed = 9) {
  set.seed(seed)
  sframe <- fmrihrf::sampling_frame(blocklens = c(60, 60), TR = 2)
  etab <- data.frame(onset = rep(seq(4, 100, by = 16), 2),
                     condition = factor(rep(c("A", "B"), length.out = 14)),
                     run = rep(1:2, each = 7))
  Y <- matrix(rnorm(120 * 3), 120, 3)
  for (r in 1:2) {
    for (v in 1:3) {
      Y[(r - 1) * 60 + 1:60, v] <- as.numeric(arima.sim(list(ar = 0.4), 60))
    }
  }
  dset <- matrix_frame(Y, TR = 2, run_length = c(60, 60), event_table = etab)
  # The nuisance column duplicates each run's intercept: one aliased column
  # per run.
  bspec <- suppressWarnings(baseline_model(basis = "poly", degree = 1, sframe = sframe,
                          intercept = "runwise",
                          nuisance_list = list(matrix(1, 60, 1), matrix(1, 60, 1))))
  list(dset = dset, bspec = bspec, Y = Y)
}

test_that("voxelwise AR residual df uses the design rank, not its width", {
  fx <- .rd_rankdef_fixture()
  fit <- suppressWarnings(fmri_lm(
    onset ~ hrf(condition), block = ~run, dataset = fx$dset,
    baseline_model = fx$bspec, durations = 0,
    control = fmri_lm_control(
      estimation = estimation_spec("runwise_meta"),
      noise = noise_spec("ar1", voxelwise = TRUE)
    )
  ))
  model <- fit$model
  X1 <- as.matrix(design_matrix(model))[1:60, ]
  X1 <- X1[, colSums(abs(X1)) > 0, drop = FALSE]
  rank1 <- qr(X1)$rank
  expect_lt(rank1, ncol(X1))
  # Each run contributes n - rank; the pooled df is their sum.
  expect_equal(fit$result$rdf, 2 * (60 - rank1))
})

test_that("rank-deficient AR fits report lm's residual df", {
  fx <- .rd_rankdef_fixture()
  fit <- suppressWarnings(fmri_lm(
    onset ~ hrf(condition), block = ~run, dataset = fx$dset,
    baseline_model = fx$bspec, durations = 0,
    control = fmri_lm_control(estimation = estimation_spec("joint"),
                              noise = noise_spec("ar1")),
    compute = compute_spec(voxel_chunks = 1)
  ))
  X <- as.matrix(design_matrix(fit$model))
  ref <- lm(fx$Y[, 1] ~ X - 1)
  expect_equal(fit$result$rdf, ref$df.residual)
  expect_equal(unique(fit$result$df$inference), ref$df.residual)
})
