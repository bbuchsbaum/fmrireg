library(testthat)

skip_if_not_installed("bidser")

AD <- fmrireg:::parse_and_validate_config
BC <- fmrireg:::build_config_from_ior
LF <- fmrireg::load_fmri_config

create_temp_bids <- function() {
  participants_df <- tibble::tibble(participant_id = "01", age = 30)
  file_structure_df <- tibble::tribble(
    ~subid, ~session, ~datatype, ~task,   ~run, ~suffix,                    ~fmriprep,
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
    run = c(1, 1)
  )

  conf_file <- bidser:::generate_bids_filename(subid="01", task="taskA", run="01", suffix="desc-confounds_timeseries.tsv")
  conf_path <- file.path("sub-01", "func", conf_file)
  conf_data <- list()
  conf_data[[conf_path]] <- tibble::tibble(motion_x = c(0.1, 0.2))

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

write_yaml <- function(lst, path) yaml::write_yaml(lst, path)


test_that("variable mapping passes when columns exist", {
  bids_dir <- create_temp_bids()
  on.exit(unlink(bids_dir, recursive = TRUE), add = TRUE)
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = bids_dir, subjects = list(include = "sub-01"), tasks = c("taskA")),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(cond = list(bids_column = "cond", role = "Factor"), motion = list(bids_column = "motion_x", role = "NuisanceSource")),
    terms = list(t1 = list(type = "EventRelated", event_variables = list("cond"))),
    models = list(list(name = "m1", terms = list("t1")))
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf), add = TRUE)

  expect_no_error(LF(tf))
})


test_that("variable mapping missing columns triggers error", {
  bids_dir <- create_temp_bids()
  on.exit(unlink(bids_dir, recursive = TRUE), add = TRUE)
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = bids_dir, subjects = list(include = "sub-01"), tasks = c("taskA")),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(cond = list(bids_column = "cond_missing", role = "Factor")),
    terms = list(t1 = list(type = "EventRelated", event_variables = list("cond"))),
    models = list(list(name = "m1", terms = list("t1")))
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf), add = TRUE)

  expect_error(LF(tf), "not found")
})

