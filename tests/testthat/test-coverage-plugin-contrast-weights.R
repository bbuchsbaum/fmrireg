# aaa_plugin_api.R: contrast offset weight resolution branches.

test_that(".resolve_contrast_fit_weights covers vector/matrix and errors", {
  event_design <- cbind(A = c(1, 0, 1, 0), B = c(0, 1, 0, 1))
  full_design <- cbind(event_design, poly1 = rnorm(4))

  expect_null(
    fmrireg:::.resolve_contrast_fit_weights(
      list(name = "c1", offset_weights = NULL),
      event_design, full_design
    )
  )

  # Unnamed vector matching event columns
  w <- fmrireg:::.resolve_contrast_fit_weights(
    list(name = "c1", offset_weights = c(1, -1)),
    event_design, full_design
  )
  expect_equal(as.numeric(w$weights), c(1, -1))
  expect_equal(as.integer(w$colind), c(1L, 2L))

  # Named vector reordered to event column order
  w2 <- fmrireg:::.resolve_contrast_fit_weights(
    list(name = "c1", offset_weights = c(B = -1, A = 1)),
    event_design, full_design
  )
  expect_equal(as.numeric(w2$weights), c(1, -1))

  expect_error(
    fmrireg:::.resolve_contrast_fit_weights(
      list(name = "c1", offset_weights = c(1, -1, 0)),
      event_design, full_design
    ),
    "offset_weights length"
  )
  expect_error(
    fmrireg:::.resolve_contrast_fit_weights(
      list(name = "c1", offset_weights = c(A = 1, Z = -1)),
      event_design, full_design
    ),
    "do not cover every event-design column"
  )

  # Matrix path with rownames
  mat <- rbind(B = c(1, 0), A = c(0, 1))
  wmat <- fmrireg:::.resolve_contrast_fit_weights(
    list(name = "F1", offset_weights = mat),
    event_design, full_design
  )
  expect_true(is.matrix(wmat$weights) || is.numeric(wmat$weights))
  expect_true(length(wmat$colind) >= 1L)

  expect_error(
    fmrireg:::.resolve_contrast_fit_weights(
      list(name = "F1", offset_weights = matrix(1, 3, 2)),
      event_design, full_design
    ),
    "offset_weights rows"
  )
  expect_error(
    fmrireg:::.resolve_contrast_fit_weights(
      list(name = "F1", offset_weights = rbind(A = c(1, 0), Z = c(0, 1))),
      event_design, full_design
    ),
    "do not cover every event-design column"
  )
})
