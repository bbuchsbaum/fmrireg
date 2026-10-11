# Performance pass 2, stream P3-5: GCV lambda selection for the soft
# projection uses the closed form
#   RSS(lambda) = ||Y||^2 - sum_k (2 s_k - s_k^2) ||u_k' Y||^2
# from precomputed summaries (O(k) per lambda) instead of forming two n x V
# matrices per evaluation; the SVD of N skips the unused right vectors.

# Dense reference: the GCV criterion evaluated on the full residual matrix.
p3_gcv_dense <- function(N, Y) {
  s <- svd(N)
  d <- s$d; U <- s$u; n <- nrow(N)
  UY <- crossprod(U, Y)
  score <- function(log_lambda) {
    sh <- d^2 / (d^2 + exp(log_lambda))
    denom <- (1 - sum(sh) / n)^2
    if (denom < 1e-10) return(Inf)
    sum((Y - U %*% (sh * UY))^2) / denom
  }
  list(score = score,
       lambda = exp(stats::optimize(score, log(c(min(d^2) / 100, max(d^2) * 100)))$minimum))
}

p3_soft_data <- function(n = 80L, k = 12L, V = 60L, seed = 4L) {
  set.seed(seed)
  N <- matrix(stats::rnorm(n * k), n, k) %*% diag(seq(3, 0.3, length.out = k))
  Y <- N %*% matrix(stats::rnorm(k * V, sd = 0.3), k, V) + matrix(stats::rnorm(n * V), n, V)
  list(N = N, Y = Y)
}

test_that("P3-5: closed-form GCV selects the same lambda and projection as the dense criterion", {
  for (seed in 1:4) {
    dat <- p3_soft_data(seed = seed)
    ref <- p3_gcv_dense(dat$N, dat$Y)
    proj <- soft_projection(dat$N, lambda = "gcv", Y = dat$Y)
    expect_equal(proj$lambda, ref$lambda, tolerance = 1e-8)
    expect_equal(proj$method, "gcv")
    # projected data equal the dense projection at the reference lambda
    s <- svd(dat$N)
    sh <- s$d^2 / (s$d^2 + ref$lambda)
    dense_proj <- dat$Y - s$u %*% (sh * crossprod(s$u, dat$Y))
    expect_equal(proj$P_lambda(dat$Y), dense_proj, tolerance = 1e-8)
    expect_equal(proj$effective_df, sum(sh), tolerance = 1e-8)
  }
})

test_that("P3-5: the closed-form criterion equals the dense GCV score", {
  dat <- p3_soft_data(n = 50L, k = 8L, V = 30L, seed = 9L)
  ref <- p3_gcv_dense(dat$N, dat$Y)
  s <- svd(dat$N, nv = 0)
  summ <- fmrireg:::.gcv_summaries(s$u, dat$Y)
  for (ll in seq(-4, 4, by = 0.5)) {
    expect_equal(fmrireg:::.gcv_score_closed_form(ll, s$d, summ, nrow(dat$N)),
                 ref$score(ll), tolerance = 1e-10)
  }
})

test_that("P3-5: wide nuisance matrices (more columns than timepoints) still work", {
  set.seed(3)
  N <- matrix(stats::rnorm(30 * 50), 30, 50)
  Y <- matrix(stats::rnorm(30 * 10), 30, 10)
  proj <- soft_projection(N, lambda = "gcv", Y = Y)
  expect_equal(proj$lambda, p3_gcv_dense(N, Y)$lambda, tolerance = 1e-8)
  proj_auto <- soft_projection(N)
  expect_equal(proj_auto$lambda, stats::median(svd(N)$d^2))
})

test_that("P3-5 guard: GCV touches Y once (summaries), not per lambda evaluation", {
  dat <- p3_soft_data()
  n_ops <- 0L
  registerS3method("Ops", "p3_counted", function(e1, e2) {
    n_ops <<- n_ops + 1L
    e1 <- unclass(e1)
    if (!missing(e2)) e2 <- unclass(e2)
    get(.Generic)(e1, e2)
  }, envir = environment())
  Yc <- structure(dat$Y, class = "p3_counted")
  lam <- fmrireg:::.select_lambda_gcv(dat$N, Yc, svd(dat$N, nu = ncol(dat$N), nv = 0))
  expect_equal(lam, p3_gcv_dense(dat$N, dat$Y)$lambda, tolerance = 1e-8)
  # one elementwise pass over Y (its squared norm); the dense criterion did
  # one n x V subtraction per optimize() evaluation
  expect_lte(n_ops, 1L)
})

test_that("P3-5 guard: soft_projection does not keep right singular vectors", {
  dat <- p3_soft_data()
  ref <- svd(dat$N)
  rec <- new.env()
  suppressMessages(trace("svd", bquote(assign("nv", c(get0("nv", envir = .(rec)), nv),
                                              envir = .(rec))),
                         print = FALSE, where = baseenv()))
  withr::defer(suppressMessages(untrace("svd", where = baseenv())))
  proj <- soft_projection(dat$N, lambda = "auto")
  expect_equal(rec$nv, 0)
  expect_null(proj$svd$v)
  expect_equal(proj$svd$d, ref$d)
  expect_equal(proj$svd$u, ref$u)
})
