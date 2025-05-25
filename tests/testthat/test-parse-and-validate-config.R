library(testthat)

AD <- fmrireg:::parse_and_validate_config

write_yaml <- function(lst, path) {
  yaml::write_yaml(lst, path)
}

context("parse_and_validate_config dataset validation")

test_that("valid minimal config passes", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(),
    terms = list(),
    models = list()
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  res <- AD(tf)
  expect_equal(res$dataset$relpath, "func")
  expect_equal(res$events$onset_column, "onset")
})

test_that("missing top-level section fails", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(events = list(onset_column = "onset", duration_column = "duration", block_column = "run"))
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  expect_error(AD(tf), "Configuration validation failed")
})

test_that("invalid subject id pattern is caught", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data", subjects = list(include = c("bad1", "sub-01"))),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(),
    terms = list(),
    models = list()
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  expect_error(AD(tf), "Configuration validation failed")
})

test_that("missing events field fails", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(duration_column = "duration", block_column = "run"),
    variables = list(),
    terms = list(),
    models = list()
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  expect_error(AD(tf), "Configuration validation failed")
})

test_that("hrfs block validates entries", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    hrfs = list(custom = list(type = "GammaFunction")),
    variables = list(),
    terms = list(),
    models = list()
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  res <- AD(tf)
  expect_equal(res$hrfs$custom$type, "GammaFunction")
})

test_that("invalid hrfs entry fails", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    hrfs = list(bad = list(type = "BadType")),
    variables = list(),
    terms = list(),
    models = list()
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  expect_error(AD(tf), "Configuration validation failed")
})

test_that("variables block validates entries", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(
      rt = list(bids_column = "reaction_time", role = "Numeric")
    ),
    terms = list(),
    models = list()
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  res <- AD(tf)
  expect_equal(res$variables$rt$role, "Numeric")
})

test_that("invalid variables entry fails", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(bad = list(bids_column = "col")),
    terms = list(),
    models = list()
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  expect_error(AD(tf), "Configuration validation failed")
})

test_that("transformations block validates ops", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(run = list(bids_column = "run", role = "Factor"), rt = list(bids_column = "rt", role = "Numeric")),
    transformations = list(
      rt_scaled = list(
        source_variable = "rt",
        ops = list("center", list(type = "scale-within-group", group_by_variable = "run"))
      )
    ),
    terms = list(),
    models = list()
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  res <- AD(tf)
  expect_equal(length(res$transformations$rt_scaled$ops), 2)
})

test_that("invalid transformations entry fails", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(x = list(bids_column = "x", role = "Numeric")),
    transformations = list(bad = list(source_variable = "x", ops = list(list(group_by_variable = "run")))),
    terms = list(),
    models = list()
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  expect_error(AD(tf), "Configuration validation failed")
})
