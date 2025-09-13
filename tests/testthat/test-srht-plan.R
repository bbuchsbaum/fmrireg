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

