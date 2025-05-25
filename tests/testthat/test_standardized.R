test_that("Standardized assigns column name", {
  x <- c(1, 2, 3)
  b <- Standardized(x)
  expect_equal(colnames(b$y), b$name)
  p <- predict(b, c(4, 5))
  expect_equal(colnames(p), b$name)
})

test_that("Standardized handles zero variance", {
  x <- rep(1, 5)
  b <- Standardized(x)
  expect_equal(b$sd, 1e-6)
  expect_equal(as.numeric(b$y[,1]), rep(0, 5))
})
