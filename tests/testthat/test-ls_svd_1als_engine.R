library(testthat)

context("ls_svd_1als_engine")

simple_ls_svd_data <- function() {
  set.seed(123)
  n <- 50
  d <- 3
  k <- 2
  v <- 4
  h_true <- matrix(rnorm(d * v), d, v)
  beta_true <- matrix(rnorm(k * v), k, v)
  X_list <- lapply(seq_len(k), function(i) matrix(rnorm(n * d), n, d))
  Xbig <- do.call(cbind, X_list)
  Y <- Xbig %*% as.vector(matrix(h_true, d, v) %*% t(beta_true))
  Y <- matrix(Y, n, v)
  list(X_list = X_list, Y = Y, d = d, k = k)
}

test_that("ls_svd_1als_engine returns matrices with correct dimensions", {
  dat <- simple_ls_svd_data()
  res <- ls_svd_1als_engine(dat$X_list, dat$Y,
                            lambda_init = 0,
                            lambda_b = 0.1,
                            lambda_h = 0.1)
  expect_equal(dim(res$h), c(dat$d, ncol(dat$Y)))
  expect_equal(dim(res$beta), c(dat$k, ncol(dat$Y)))
  expect_equal(dim(res$h_ls_svd), c(dat$d, ncol(dat$Y)))
  expect_equal(dim(res$beta_ls_svd), c(dat$k, ncol(dat$Y)))
})

single_condition_data <- function() {
  set.seed(42)
  n <- 40
  d <- 2
  v <- 3
  X <- matrix(rnorm(n * d), n, d)
  h_true <- matrix(rnorm(d * v), d, v)
  b_true <- rnorm(v)
  Y <- (X %*% h_true) * matrix(rep(b_true, each = n), n, v)
  list(X_list = list(X), Y = Y)
}

test_that("fullXtX_flag has no effect for single condition", {
  dat <- single_condition_data()
  res_diag <- ls_svd_1als_engine(dat$X_list, dat$Y,
                                 lambda_init = 0,
                                 lambda_b = 0.1,
                                 lambda_h = 0.1,
                                 fullXtX_flag = FALSE)
  res_full <- ls_svd_1als_engine(dat$X_list, dat$Y,
                                 lambda_init = 0,
                                 lambda_b = 0.1,
                                 lambda_h = 0.1,
                                 fullXtX_flag = TRUE)
  expect_equal(res_diag$h, res_full$h)
  expect_equal(res_diag$beta, res_full$beta)
})

correlated_data <- function() {
  set.seed(99)
  n <- 40
  d <- 2
  k <- 2
  v <- 2
  baseX <- matrix(rnorm(n * d), n, d)
  X_list <- list(baseX, baseX + matrix(rnorm(n * d, sd = 0.2), n, d))
  h_true <- matrix(rnorm(d * v), d, v)
  beta_true <- matrix(rnorm(k * v), k, v)
  Xbig <- do.call(cbind, X_list)
  Y <- Xbig %*% as.vector(matrix(h_true, d, v) %*% t(beta_true))
  Y <- matrix(Y, n, v)
  list(X_list = X_list, Y = Y)
}

test_that("fullXtX_flag influences estimates when conditions correlate", {
  dat <- correlated_data()
  res_diag <- ls_svd_1als_engine(dat$X_list, dat$Y,
                                 lambda_init = 0,
                                 lambda_b = 0.1,
                                 lambda_h = 0.1,
                                 fullXtX_flag = FALSE)
  res_full <- ls_svd_1als_engine(dat$X_list, dat$Y,
                                 lambda_init = 0,
                                 lambda_b = 0.1,
                                 lambda_h = 0.1,
                                 fullXtX_flag = TRUE)
  expect_false(isTRUE(all.equal(res_diag$h, res_full$h)))
})
