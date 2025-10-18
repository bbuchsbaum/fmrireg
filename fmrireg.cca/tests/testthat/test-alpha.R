test_that("alpha-mix enforces plausible HRF envelope (cca2)", {
  skip_on_cran()
  skip_if_not_installed("fmrireg.cca")
  library(fmrireg.cca)

  tr <- 2; span <- 30
  basis <- fmrireg.cca:::cca_basis(alpha = 0.3, TR = tr, span = span)
  tt <- seq(0, span, by = tr)
  vals <- fmrihrf::evaluate(basis, tt)
  expect_true(is.matrix(vals) && ncol(vals) >= 2)
  y1t <- vals[,1]; y2t <- vals[,2]

  alpha <- 0.3
  y1 <- 0.5 * (y1t + y2t)
  y2 <- (y1t - y2t) / (2*alpha)
  y1 <- y1 / sqrt(sum(y1^2)); y2 <- y2 / sqrt(sum(y2^2))
  for (i in 1:50) {
    w <- runif(2)
    y <- w[1]*y1t + w[2]*y2t
    coef_y1 <- sum(y * y1); coef_y2 <- sum(y * y2)
    expect_lte(abs(coef_y2), alpha * coef_y1 + 1e-6)
  }
})

