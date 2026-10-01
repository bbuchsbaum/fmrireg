test_that("contrast accessors return named values without name-repair noise", {
  block <- tibble::tibble(
    estimate = c(2, 4), se = c(0.5, 1),
    stat = c(4, 4), prob = c(0.01, 0.02)
  )
  fit <- structure(list(result = list(contrasts = tibble::tibble(
    type = c("contrast", "contrast", "Fcontrast"),
    name = c("A_vs_B", "B_vs_C", "overall"),
    data = list(block, block * 2, block * 3)
  ))), class = "fmri_lm")

  accessors <- list(
    estimate = function(x) coef(x, type = "contrasts"),
    se = function(x) standard_error(x, type = "contrasts"),
    stat = function(x) stats(x, type = "contrasts"),
    prob = function(x) p_values(x, type = "contrasts")
  )
  for (element in names(accessors)) {
    expect_silent(out <- accessors[[element]](fit))
    expect_equal(out, tibble::tibble(
      A_vs_B = block[[element]], B_vs_C = block[[element]] * 2
    ))
    expect_silent(revised <- fmrireg:::pull_stat_revised(fit, "contrasts", element))
    expect_equal(revised, out)
  }

  expect_silent(out <- stats(fit, type = "F"))
  expect_equal(out, tibble::tibble(overall = block$stat * 3))
  expect_silent(revised <- fmrireg:::pull_stat_revised(fit, "F", "stat"))
  expect_equal(revised, out)

  fit$result$contrasts <- fit$result$contrasts[1, ]
  expect_silent(out <- stats(fit, type = "contrasts"))
  expect_equal(out, tibble::tibble(A_vs_B = block$stat))
})
