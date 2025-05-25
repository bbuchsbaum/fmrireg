library(testthat)
AD <- fmrireg:::parse_and_validate_config
BC <- fmrireg:::build_config_from_ior
LF <- fmrireg::load_fmri_config

write_yaml <- function(lst, path) yaml::write_yaml(lst, path)

test_that("nonexistent dataset path fails", {
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = "./nonexistent"),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(),
    terms = list(),
    models = list()
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf))

  expect_error(LF(tf), "BIDS directory does not exist")
})
