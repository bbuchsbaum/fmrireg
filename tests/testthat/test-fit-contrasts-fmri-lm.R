test_that("fit_contrasts.fmri_lm computes non-placeholder t and F stats", {
  # p x V coefficients
  B <- matrix(c(1, 2,
                3, 4), nrow = 2, byrow = FALSE)
  # Stored as V x p in beta stats table
  betas_vxp <- t(B)

  fake <- list(
    result = list(
      betas = tibble::tibble(
        data = list(tibble::tibble(estimate = list(betas_vxp)))
      ),
      sigma = c(2, 3),
      rdf = 10,
      cov.unscaled = diag(c(4, 9))
    )
  )
  class(fake) <- "fmri_lm"

  cres <- fit_contrasts(
    fake,
    list(
      c1 = c(1, -1),
      f1 = rbind(c(1, 0), c(0, 1))
    )
  )

  expect_true(all(c("c1", "f1") %in% names(cres)))

  c1 <- cres$c1
  expect_equal(c1$type, "contrast")
  expect_equal(c1$estimate, c(-1, -1), tolerance = 1e-12)
  expect_true(all(c1$se > 0))
  expect_true(all(is.finite(c1$stat)))
  expect_true(all(c1$prob > 0 & c1$prob < 1))
  expect_false(all(abs(c1$se - 1) < 1e-12))
  expect_false(all(abs(c1$prob - 0.05) < 1e-12))

  f1 <- cres$f1
  expect_equal(f1$type, "Fcontrast")
  expect_true(all(f1$se > 0))
  expect_true(all(is.finite(f1$stat)))
  expect_true(all(f1$prob > 0 & f1$prob < 1))
})
