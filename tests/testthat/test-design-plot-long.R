# design_plot() reshapes design matrices with an internal helper instead of
# tidyr::pivot_longer(). These tests pin that the long data are unchanged.

tidyr_design_long <- function(dm) {
  tidyr::pivot_longer(
    dm,
    cols = -c(".time", ".block"),
    names_to = "condition",
    values_to = "value"
  )
}

test_that(".pivot_design_long matches tidyr::pivot_longer", {
  skip_if_not_installed("tidyr")

  dm <- tibble::tibble(
    `cond_A` = c(0.1, 0.2, 0.3, 0.4),
    `cond B` = c(1, 2, 3, 4),
    poly = c(-1, 0, NA, 1),
    .block = c(1L, 1L, 2L, 2L),
    .time = c(0, 2, 4, 6)
  )
  expect_identical(fmrireg:::.pivot_design_long(dm), tidyr_design_long(dm))

  # Id columns keep their position relative to each other.
  dm2 <- dm[, c(".time", "poly", ".block", "cond_A")]
  expect_identical(fmrireg:::.pivot_design_long(dm2), tidyr_design_long(dm2))

  # A single regressor and a zero-row design.
  dm3 <- dm[, c("poly", ".block", ".time")]
  expect_identical(fmrireg:::.pivot_design_long(dm3), tidyr_design_long(dm3))
  expect_identical(
    fmrireg:::.pivot_design_long(dm[0, ]),
    tidyr_design_long(dm[0, ])
  )

  # Mixed integer/double regressors take the common double type.
  dm4 <- dm
  dm4$cond_A <- 1:4
  expect_identical(fmrireg:::.pivot_design_long(dm4), tidyr_design_long(dm4))
})

test_that("design_plot long data are identical to the tidyr reshape", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("plotly")
  skip_if_not_installed("bslib")
  skip_if_not_installed("tidyr")

  etab <- data.frame(
    onset = c(5, 15, 25, 5, 15, 25),
    condition = factor(c("A", "B", "A", "B", "A", "B")),
    modulator = c(0.5, -1, 1.5, 0.2, -0.3, 1),
    run = rep(1:2, each = 3)
  )
  Y <- matrix(rnorm(80 * 2), 80, 2)
  dset <- matrix_frame(Y, TR = 1, run_length = c(40, 40), event_table = etab)
  emod <- event_model(
    onset ~ hrf(condition) + hrf(modulator), data = etab, block = ~run,
    sampling_frame = frame_sframe(dset)
  )
  bmod <- baseline_model(basis = "poly", degree = 2, sframe = frame_sframe(dset))
  fmod <- fmri_model(emod, bmod, dset)

  design_long <- function(longnames) {
    app <- design_plot(fmod, longnames = longnames)
    get("df_long", envir = environment(app$serverFuncSource()))
  }

  for (longnames in c(TRUE, FALSE)) {
    new <- design_long(longnames)
    old <- with_mocked_bindings(
      design_long(longnames),
      .pivot_design_long = tidyr_design_long
    )

    expect_identical(names(new), names(old))
    expect_gte(length(new), 3L)
    for (term in names(old)) {
      expect_identical(names(new[[term]]), c(".block", ".time", "condition", "value"))
      expect_identical(new[[term]], old[[term]])
    }
  }
})
