skip_if_not_installed("bidser")

LF <- fmrireg::load_fmri_config
LD <- fmrireg:::load_and_prepare_subject_data
BM <- fmrireg:::build_baseline_model_from_dsl

write_yaml <- function(lst, path) yaml::write_yaml(lst, path)

create_temp_bids501 <- function() {
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
    run = c(1, 1)
  )

  conf_file <- bidser:::generate_bids_filename(subid="01", task="taskA", run="01", suffix="desc-confounds_timeseries.tsv")
  conf_path <- file.path("sub-01", "func", conf_file)
  conf_data <- list()
  conf_data[[conf_path]] <- tibble::tibble(motion_x = c(0.1,0.2), motion_y = c(0.2,0.3), run = c(1,1))

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


test_that("baseline_model built with confound groups and shorthand basis", {
  bids_dir <- create_temp_bids501()
  on.exit(unlink(bids_dir, recursive = TRUE), add = TRUE)
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = bids_dir,
                   subjects = list(include = "sub-01"),
                   tasks = c("taskA"),
                   scan_params = list(TR = 2, run_lengths = list(taskA = 2))),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(cond = list(bids_column = "cond", role = "Factor")),
    confound_groups = list(motion = list(select_by_pattern = c("motion_.*"))),
    terms = list(t1 = list(type = "EventRelated", event_variables = list("cond"))),
    models = list(list(name="m1", terms=list("t1"), baseline=list(basis="BSpline(3)", include_confound_groups=list("motion"))))
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf), add = TRUE)

  config <- LF(tf)
  sdat <- LD(config, "01")
  blmod <- BM(config, config$spec$models[[1]], sdat)
  expect_s3_class(blmod, "baseline_model")
  expect_true(!is.null(blmod$terms$nuisance))
})
