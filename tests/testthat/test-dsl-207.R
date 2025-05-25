library(testthat)

AD <- fmrireg:::parse_and_validate_config

write_yaml <- function(lst, path) yaml::write_yaml(lst, path)

context("DSL-207 role/type compatibility")

test_that("parametric modulation requires numeric mod_var", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(
      cond = list(bids_column = "cond", role = "Factor"),
      mod  = list(bids_column = "rt", role = "Factor")
    ),
    terms = list(pm = list(type = "ParametricModulation", selector_vars = list("cond"), mod_var = "mod")),
    models = list(list(name = "m", terms = list("pm")))
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  expect_error(AD(tf), "mod_var")
})

test_that("modulator_basis with non-numeric mod_var errors", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(
      cond = list(bids_column = "cond", role = "Factor"),
      mod  = list(bids_column = "rt", role = "Factor")
    ),
    terms = list(pm = list(type = "ParametricModulation", selector_vars = list("cond"), mod_var = "mod",
                           modulator_basis = list(type = "Polynomial", parameters = list(degree = 2)))),
    models = list(list(name = "m", terms = list("pm")))
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  expect_error(AD(tf), "modulator_basis")
})

test_that("event variable with wrong role errors", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(cond = list(bids_column = "cond", role = "NuisanceSource")),
    terms = list(t1 = list(type = "EventRelated", event_variables = list("cond"))),
    models = list(list(name = "m", terms = list("t1")))
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  expect_error(AD(tf), "event_variables")
})

test_that("nuisance regressors require nuisance role", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(cond = list(bids_column = "cond", role = "Factor")),
    terms = list(nu = list(type = "NuisanceRegressors", nuisance_source_variables = list("cond"))),
    models = list(list(name = "m", terms = list("nu")))
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  expect_error(AD(tf), "nuisance_source_variables")
})
