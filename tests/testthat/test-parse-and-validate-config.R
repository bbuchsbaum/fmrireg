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
