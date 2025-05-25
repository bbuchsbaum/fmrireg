skip_if_not_installed("bidser")

LF <- fmrireg::load_fmri_config
BF <- fmrireg::build_fmri_model_from_config

write_yaml <- function(lst, path) yaml::write_yaml(lst, path)

create_temp_bids504 <- function() {
  participants_df <- tibble::tibble(participant_id = "01", age = 30)
  file_structure_df <- tibble::tribble(
    ~subid, ~session, ~datatype, ~task, ~run, ~suffix, ~fmriprep,
    "01",   NA,       "func",    "taskA", "01", "bold",  FALSE,
    "01",   NA,       "func",    "taskA", "01", "events", FALSE
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

  tmp_dir <- tempfile("bids_")
  dir.create(tmp_dir, recursive = TRUE)
  bidser::create_mock_bids(
    project_name = "Mock",
    participants = participants_df,
    file_structure = file_structure_df,
    event_data = event_data,
    create_stub = TRUE,
    stub_path = tmp_dir
  )
  tmp_dir
}


test_that("fmri_model built from config", {
  bids_dir <- create_temp_bids504()
  on.exit(unlink(bids_dir, recursive = TRUE), add = TRUE)
  tf <- tempfile(fileext = ".yml")
  cfg <- list(
    dataset = list(path = bids_dir,
                   subjects = list(include = "sub-01"),
                   tasks = c("taskA"),
                   scan_params = list(TR = 2, run_lengths = list(taskA = 2))),
    events = list(onset_column = "onset", duration_column = "duration", block_column = "run"),
    variables = list(cond = list(bids_column = "cond", role = "Factor")),
    terms = list(t1 = list(type = "EventRelated", event_variables = list("cond"))),
    models = list(list(name="m1", terms=list("t1"), baseline=list(basis="Constant")))
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf), add = TRUE)

  config <- LF(tf)
  fm <- BF(config, "01")
  expect_s3_class(fm, "fmri_model")
})
