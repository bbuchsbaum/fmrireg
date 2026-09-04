# group_data_h5 helpers: mocked read_h5_metadata, detect type, subject ids.

test_that("read_h5_metadata covers list and S4 mocked handles", {
  skip_if_not_installed("fmristore")
  skip_if_not_installed("methods")

  fake_list <- list(
    dim = c(3L, 3L, 1L, 3L),
    labels = c("beta", "se", "tstat"),
    mask_dim = c(3L, 3L, 1L),
    mask = TRUE,
    close_all = function() invisible(NULL),
    close = function() invisible(NULL)
  )

  meta_list <- with_mocked_bindings(
    {
      fmrireg:::read_h5_metadata("mock-list.h5")
    },
    read_labeled_vec = function(path) fake_list,
    .package = "fmristore"
  )
  expect_equal(meta_list$dim, c(3L, 3L, 1L, 3L))
  expect_equal(meta_list$labels, c("beta", "se", "tstat"))
  expect_true(isTRUE(meta_list$has_mask))

  methods::setClass(
    "CovFakeLabeledVol",
    slots = c(labels = "character", mask = "ANY", dim = "integer"),
    where = topenv()
  )
  on.exit(methods::removeClass("CovFakeLabeledVol", where = topenv()), add = TRUE)

  fake_s4 <- methods::new(
    "CovFakeLabeledVol",
    labels = c("FacesVsPlaces", "GoVsNoGo"),
    mask = array(TRUE, dim = c(2, 2, 1)),
    dim = c(2L, 2L, 1L, 2L)
  )

  meta_s4 <- with_mocked_bindings(
    {
      fmrireg:::read_h5_metadata("mock-s4.h5")
    },
    read_labeled_vec = function(path) fake_s4,
    .package = "fmristore"
  )
  expect_equal(as.integer(meta_s4$dim[1:3]), c(2L, 2L, 1L))
  expect_equal(meta_s4$labels, c("FacesVsPlaces", "GoVsNoGo"))
  expect_true(isTRUE(meta_s4$has_mask))

  expect_error(
    with_mocked_bindings(
      fmrireg:::read_h5_metadata("bad.h5"),
      read_labeled_vec = function(path) structure(1:3, class = "weird"),
      .package = "fmristore"
    ),
    "Failed to read HDF5 metadata|Unsupported"
  )
})

test_that("detect_h5_file_type and extract_subject_ids cover branches", {
  expect_identical(
    fmrireg:::detect_h5_file_type(list(labels = c("beta", "se", "pval"))),
    "by_stat"
  )
  expect_identical(
    fmrireg:::detect_h5_file_type(list(labels = c("beta_01", "reg_face", "Intercept"))),
    "betas"
  )
  expect_identical(
    fmrireg:::detect_h5_file_type(list(labels = c("FaceVsScene", "LeftVsRight"))),
    "by_contrast"
  )

  bids_ids <- fmrireg:::extract_subject_ids_from_paths(c(
    "/proj/sub-01/func/sub-01_task-nback_desc-GLMstatmap_bold.h5",
    "/proj/sub-99/func/sub-99_task-nback_desc-GLMstatmap_bold.h5"
  ))
  expect_equal(bids_ids, c("sub-01", "sub-99"))

  plain <- fmrireg:::extract_subject_ids_from_paths(c("alpha.nii.gz", "beta.nii.gz"))
  expect_equal(plain, c("alpha", "beta"))

  mixed <- fmrireg:::extract_subject_ids_from_paths(c("no_prefix_here.h5", "sub-07_x.h5"))
  expect_equal(mixed[[2]], "sub-07")
})

test_that("group_data_from_h5 reaches metadata/contrast gates with mocked reader", {
  skip_if_not_installed("fmristore")

  fake <- list(
    dim = c(2L, 2L, 1L, 2L),
    labels = c("beta", "se"),
    mask_dim = c(2L, 2L, 1L),
    mask = TRUE,
    close_all = function() invisible(NULL),
    close = function() invisible(NULL)
  )

  paths <- c(tempfile(fileext = ".h5"), tempfile(fileext = ".h5"))
  file.create(paths)
  on.exit(unlink(paths), add = TRUE)

  gd <- with_mocked_bindings(
    {
      suppressWarnings(group_data_from_h5(
        paths,
        subjects = c("sub-01", "sub-02"),
        covariates = data.frame(age = c(20, 30)),
        contrast = NULL,
        stat = c("beta", "se"),
        validate = TRUE
      ))
    },
    read_labeled_vec = function(path) fake,
    .package = "fmristore"
  )
  expect_true(inherits(gd, "group_data_h5") || inherits(gd, "group_data"))
  expect_equal(gd$subjects, c("sub-01", "sub-02"))
  expect_identical(gd$file_type, "by_stat")

  expect_error(
    with_mocked_bindings(
      {
        suppressWarnings(group_data_from_h5(
          paths,
          subjects = c("sub-01", "sub-02"),
          contrast = "MissingContrast",
          validate = TRUE
        ))
      },
      read_labeled_vec = function(path) fake,
      .package = "fmristore"
    ),
    "not found|Available"
  )
})
