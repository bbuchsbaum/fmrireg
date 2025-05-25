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

test_that("missing hrfs block uses canonical default", {
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
  expect_equal(res$hrfs$canonical$type, "SPMCanonical")
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

test_that("confound_groups block validates entries", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(),
    confound_groups = list(motion = list(select_by_pattern = c("trans.*", "rot.*"))),
    terms = list(),
    models = list()
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  res <- AD(tf)
  expect_true("motion" %in% names(res$confound_groups))
})

test_that("confound_groups entry missing selectors fails", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(),
    confound_groups = list(bad = list()),
    terms = list(),
    models = list()
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  expect_error(AD(tf), "Configuration validation failed")
})

test_that("terms block validates conditional fields", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(cond = list(bids_column = "condition", role = "Factor")),
    terms = list(myterm = list(type = "EventRelated", event_variables = list("cond"))),
    models = list()
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  res <- AD(tf)
  expect_equal(res$terms$myterm$type, "EventRelated")
})

test_that("parametric modulation term missing mod_var fails", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(cond = list(bids_column = "condition", role = "Factor"),
                     rt = list(bids_column = "rt", role = "Numeric")),
    terms = list(pm = list(type = "ParametricModulation", selector_vars = list("cond"))),
    models = list()
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  expect_error(AD(tf), "Configuration validation failed")
})

test_that("models block validates entries", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(cond = list(bids_column = "condition", role = "Factor")),
    terms = list(t1 = list(type = "EventRelated", event_variables = list("cond"))),
    models = list(list(name = "m1", terms = list("t1")))
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  res <- AD(tf)
  expect_equal(res$models[[1]]$name, "m1")
  expect_equal(res$models[[1]]$baseline$intercept, "PerRun")
})

test_that("invalid model entry fails", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(cond = list(bids_column = "condition", role = "Factor")),
    terms = list(t1 = list(type = "EventRelated", event_variables = list("cond"))),
    models = list(list(terms = list("t1")))
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  expect_error(AD(tf), "Configuration validation failed")
})

test_that("contrasts and settings validate", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(cond = list(bids_column = "condition", role = "Factor")),
    terms = list(t1 = list(type = "EventRelated", event_variables = list("cond"))),
    contrasts = list(c1 = list(type = "Formula", expression = "t1")),
    models = list(list(name = "m1", terms = list("t1"))),
    default_model = "m1",
    validation_settings = list(cross_references = "Warn", bids_content_checks = "Off")
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  res <- AD(tf)
  expect_equal(res$contrasts$c1$type, "Formula")
  expect_equal(res$default_model, "m1")
  expect_equal(res$validation_settings$cross_references, "Warn")
})

test_that("undefined HRF reference triggers error", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(cond = list(bids_column = "condition", role = "Factor")),
    terms = list(t1 = list(type = "EventRelated", event_variables = list("cond"), hrf = "missing")),
    models = list(list(name = "m1", terms = list("t1")))
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  expect_error(AD(tf), "Configuration validation failed")
})

test_that("missing term reference warns when cross_references Warn", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./data"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(cond = list(bids_column = "condition", role = "Factor")),
    terms = list(t1 = list(type = "EventRelated", event_variables = list("cond"))),
    models = list(list(name = "m1", terms = list("t1", "t2"))),
    validation_settings = list(cross_references = "Warn")
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  expect_warning(res <- AD(tf))
  expect_equal(res$models[[1]]$name, "m1")
})
