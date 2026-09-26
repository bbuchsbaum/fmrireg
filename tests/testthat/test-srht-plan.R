test_that("SRHT plan yields correct dims and deterministic sketches", {
  skip_on_cran()

  set.seed(11)
  Tlen <- 150L  # not a power of two
  p <- 5L; k <- 3L
  X <- matrix(rnorm(Tlen * p), Tlen, p)
  Z <- matrix(rnorm(Tlen * k), Tlen, k)
  m <- 4L * p

  plan <- fmrireg:::make_srht_plan(Tlen, m)
  Xs1 <- fmrireg:::srht_apply(X, plan)
  Zs1 <- fmrireg:::srht_apply(Z, plan)

  # same plan -> identical output
  Xs2 <- fmrireg:::srht_apply(X, plan)
  Zs2 <- fmrireg:::srht_apply(Z, plan)

  expect_equal(dim(Xs1), c(m, p))
  expect_equal(dim(Zs1), c(m, k))
  expect_equal(Xs1, Xs2)
  expect_equal(Zs1, Zs2)

  # G = Xs'Xs PSD and sized
  G <- crossprod(Xs1)
  ev <- eigen(G, symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(ev > -1e-8))
})


test_that("SRHT kernels are thread-count invariant and consistent with dense S", {
  # Guards the column-parallel FWHT: results must not depend on the thread
  # count, apply must be linear column by column (S X = S I X), and the
  # adjoint must be the exact transpose, including non-power-of-two T.
  for (Tlen in c(7L, 150L, 256L)) {
    set.seed(Tlen)
    m <- max(1L, Tlen %/% 3L)
    plan <- fmrireg:::make_srht_plan(Tlen, m)
    X <- matrix(rnorm(Tlen * 17), Tlen, 17)
    B <- matrix(rnorm(m * 5), m, 5)
    run <- function(threads, f) withr::with_options(
      list(fmrireg.num_threads = threads), f())
    Xs1 <- run(1L, function() fmrireg:::srht_apply(X, plan))
    Xs4 <- run(4L, function() fmrireg:::srht_apply(X, plan))
    expect_identical(Xs1, Xs4)
    S <- run(1L, function() fmrireg:::srht_apply(diag(Tlen), plan))
    expect_equal(Xs4, S %*% X, tolerance = 1e-12)
    At1 <- run(1L, function() fmrireg:::srht_adjoint(B, plan))
    At4 <- run(4L, function() fmrireg:::srht_adjoint(B, plan))
    expect_identical(At1, At4)
    expect_equal(At4, t(S) %*% B, tolerance = 1e-12)
  }
})
