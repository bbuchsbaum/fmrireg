library(testthat)

skip_if_not_installed("bidser")

LF <- fmrireg::load_fmri_config
LD <- fmrireg:::load_and_prepare_subject_data
PV <- fmrireg:::process_variables_and_transformations

create_temp_bids301b <- function() {
  participants_df <- tibble::tibble(participant_id = "01", age = 30)
  file_structure_df <- tibble::tribble(
    ~subid, ~session, ~datatype, ~task,   ~run, ~suffix,                    ~fmriprep,
    "01",   NA,       "func",    "taskA", "01", "bold",                   FALSE,
    "01",   NA,       "func",    "taskA", "01", "events",                 FALSE
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

write_yaml <- function(lst, path) yaml::write_yaml(lst, path)


test_that("variable processing generates expected env", {
  bids_dir <- create_temp_bids301b()
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
    terms = list(t1 = list(type="ParametricModulation", selector_vars=list("cond"), mod_var="rt_c",
                           modulator_basis=list(type="Polynomial", parameters=list(degree=2)))),
    models = list(list(name="m1", terms=list("t1")))
  )
  write_yaml(cfg, tf)
  on.exit(unlink(tf), add = TRUE)

  config <- LF(tf)
  sdat <- LD(config, "01")
  env <- PV(config, sdat)
  expect_true(all(c("cond","rt","rt_c") %in% ls(env)))
  expect_equal(round(env$rt_c,1), c(-0.5,0.5))
  poly_var <- grep("rt_c_", ls(env), value=TRUE)
  expect_true(length(poly_var)==1)
  expect_equal(ncol(env[[poly_var]]), 2)
})
