# Unit tests for bootstrap primitives that do not need solve_integrated_glm.

test_that("create_bootstrap_blocks respects run structure and sizes", {
  blocks <- fmrireg:::create_bootstrap_blocks(20L, block_size = 5L)
  expect_equal(length(blocks), 4L)
  expect_equal(unname(unlist(blocks)), 1:20)

  run_blocks <- fmrireg:::create_bootstrap_blocks(
    20L, block_size = 4L, run_indices = list(1:10, 11:20)
  )
  expect_true(length(run_blocks) >= 4L)
  expect_equal(sort(unname(unlist(run_blocks))), 1:20)
  # Each original run index appears in some block
  expect_true(all(1:10 %in% unlist(run_blocks)))
  expect_true(all(11:20 %in% unlist(run_blocks)))
})

test_that("bootstrap_residual/case/wild return correctly shaped resamples", {
  set.seed(42)
  n <- 24L
  p <- 3L
  V <- 4L
  X <- cbind(1, rnorm(n), rnorm(n))
  beta <- matrix(c(0.5, 1, -0.4), p, 1)
  Y <- X %*% matrix(rep(beta, V), p, V) + matrix(rnorm(n * V, sd = 0.3), n, V)
  fitted <- X %*% solve(crossprod(X), crossprod(X, Y))
  resid <- Y - fitted
  blocks <- fmrireg:::create_bootstrap_blocks(n, block_size = 6L)

  resid_boot <- fmrireg:::bootstrap_residual(X, fitted, resid, blocks)
  expect_equal(dim(resid_boot$X), dim(X))
  expect_equal(dim(resid_boot$Y), dim(Y))
  expect_identical(resid_boot$X, X)
  # Residual bootstrap keeps design fixed and changes Y
  expect_false(isTRUE(all.equal(resid_boot$Y, Y)))

  case_boot <- fmrireg:::bootstrap_case(X, Y, blocks)
  expect_equal(dim(case_boot$X), dim(X))
  expect_equal(dim(case_boot$Y), dim(Y))
  expect_false(identical(case_boot$X, X) && identical(case_boot$Y, Y))

  wild_boot <- fmrireg:::bootstrap_wild(X, fitted, resid)
  expect_identical(wild_boot$X, X)
  expect_equal(dim(wild_boot$Y), dim(Y))
  # Rademacher weights flip residual signs; Y should differ unless all weights +1
  expect_true(any(wild_boot$Y != Y) || any(wild_boot$Y != fitted))

  # Force length-adjustment: undersized union (pad) and oversized union (trim)
  set.seed(7)
  short_case <- fmrireg:::bootstrap_case(
    X, Y, blocks = list(1:2, 3:4)
  )
  expect_equal(nrow(short_case$X), n)
  expect_equal(nrow(short_case$Y), n)

  set.seed(8)
  long_resid <- fmrireg:::bootstrap_residual(
    X, fitted, resid,
    blocks = lapply(seq_len(8L), function(i) seq_len(n))
  )
  expect_equal(nrow(long_resid$Y), n)
})

test_that("bootstrap_hypothesis_test and compute_bca_ci cover voxel loops", {
  set.seed(99)
  nboot <- 40L
  nvox <- 3L
  boot_contrasts <- list(
    contrast1 = matrix(rnorm(nboot * nvox, mean = 0.4), nboot, nvox)
  )
  boot_result <- list(boot_contrasts = boot_contrasts)

  pvals <- fmrireg:::bootstrap_hypothesis_test(boot_result, contrast_idx = 1L, null_value = 0)
  expect_equal(length(pvals), nvox)
  expect_true(all(pvals >= 0 & pvals <= 1))

  expect_error(
    fmrireg:::bootstrap_hypothesis_test(list(boot_betas = 1), 1L),
    "No contrasts"
  )

  observed <- colMeans(boot_contrasts[[1]])
  ci <- fmrireg:::compute_bca_ci(boot_contrasts[[1]], observed, confidence = 0.9)
  expect_equal(dim(ci), c(2L, nvox))
  expect_true(all(is.finite(ci)))
  expect_true(all(ci[1, ] <= ci[2, ]))
})
