test_that("covariate data length must match sampling frame", {
  sframe <- sampling_frame(blocklens = c(50, 50), TR = 1)

  # create covariate data that is too short
  bad_dat <- data.frame(x = rnorm(75), y = rnorm(75))

  expect_error(
    event_model(onset ~ covariate(x, y, data = bad_dat),
                data = data.frame(onset = seq_len(100), run = rep(1:2, each = 50)),
                block = ~ run, sampling_frame = sframe),
    "sampling_frame expects"
  )
})
