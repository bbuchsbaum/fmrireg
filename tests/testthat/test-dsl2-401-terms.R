skip_if_not_installed("bidser")

LF <- fmrireg::load_fmri_config
LD <- fmrireg:::load_and_prepare_subject_data
PV <- fmrireg:::process_variables_and_transformations
CT <- fmrireg:::convert_terms_to_specs

write_yaml <- function(lst, path) yaml::write_yaml(lst, path)

# Helper dataset from DSL2-301 test
create_temp_bids401 <- function() {
  participants_df <- tibble::tibble(participant_id = "01", age = 30)
  file_structure_df <- tibble::tribble(
    ~subid, ~session, ~datatype, ~task, ~run, ~suffix,                    ~fmriprep,
    "01",   NA,       "func",    "taskA", "01", "bold",                   FALSE,
    "01",   NA,       "func",    "taskA", "01", "events",                 FALSE,
    "01",   NA,       "func",    "taskA", "01", "desc-confounds_timeseries", TRUE
  )

  ev_file <- bidser:::generate_bids_filename(subid="01", task="taskA", run="01", suffix="events.tsv")
  ev_path <- file.path("sub-01", "func", ev_file)
  event_data <- list()
  event_data[[ev_path]] <- tibble::tibble(
    onset = c(1, 3),
    duration = c(1, 1),
    cond = c("a", "b"),
    rt = c(1.2, 2.2),
    run = c(1, 1)
  )

  conf_file <- bidser:::generate_bids_filename(subid="01", task="taskA", run="01", suffix="desc-confounds_timeseries.tsv")
  conf_path <- file.path("sub-01", "func", conf_file)
  conf_data <- list()
  conf_data[[conf_path]] <- tibble::tibble(motion_x = c(0.1,0.2), motion_y = c(0.2,0.3))

  tmp_dir <- tempfile("bids_")
  dir.create(tmp_dir, recursive = TRUE)
  bidser::create_mock_bids(
    project_name = "Mock",
    participants = participants_df,
    file_structure = file_structure_df,
    event_data = event_data,
    confound_data = conf_data,
    create_stub = TRUE,
    stub_path = tmp_dir
  )
  tmp_dir
}

# Parametric modulation term

test_that("PM term references basis-expanded variable", {
  bids_dir <- create_temp_bids401()
  on.exit(unlink(bids_dir, recursive = TRUE), add = TRUE)
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = bids_dir, subjects = list(include = "sub-01"), tasks = c("taskA")),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(
      cond = list(bids_column = "cond", role = "Factor"),
      rt   = list(bids_column = "rt", role = "Numeric")
    ),
    transformations = list(rt_c = list(source_variable="rt", ops=list("center"))),
    terms = list(pm = list(type="ParametricModulation", selector_vars=list("cond"), mod_var="rt_c",
                           modulator_basis=list(type="Polynomial", parameters=list(degree=2)))),
    models = list(list(name="m1", terms=list("pm")))
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf), add = TRUE)

  config <- LF(tf)
  sdat <- LD(config, "01")
  env <- PV(config, sdat)
  specs <- CT(config, config$spec$models[[1]], env)
  expect_named(specs, "pm")
  expect_s3_class(specs$pm, "hrfspec")
  vars <- vapply(specs$pm$vars, rlang::as_label, character(1))
  expect_true(any(grepl("rt_c_polynomial_deg2", vars)))
})

# NuisanceRegressor term

test_that("nuisance term creates covariatespec", {
  bids_dir <- create_temp_bids401()
  on.exit(unlink(bids_dir, recursive = TRUE), add = TRUE)
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = bids_dir, subjects = list(include = "sub-01"), tasks = c("taskA")),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(
      motion_x = list(bids_column = "motion_x", role = "NuisanceSource"),
      motion_y = list(bids_column = "motion_y", role = "NuisanceSource")
    ),
    terms = list(nu = list(type="NuisanceRegressors", nuisance_source_variables=list("motion_x","motion_y"))),
    models = list(list(name="m1", terms=list("nu")))
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf), add = TRUE)

  config <- LF(tf)
  sdat <- LD(config, "01")
  env <- PV(config, sdat)
  specs <- CT(config, config$spec$models[[1]], env)
  expect_named(specs, "nu")
  expect_s3_class(specs$nu, "covariatespec")
  expect_equal(specs$nu$varnames, c("motion_x","motion_y"))
})
