# Review 2026-10, pass 3, stream R3-1: optional dependencies.
#
# plotly, bslib and fmrigds moved from Imports to Suggests. design_plot()
# needs shiny + plotly + bslib and group_data() needs fmrigds; both must fail
# with an error that names the missing package. Before this pass,
# with_package() handed a character vector to requireNamespace(), which only
# checked the first element, so a missing plotly/bslib went unnoticed until the
# app was built.

.r3_fmri_model <- function() {
  etab <- data.frame(
    onset = c(5, 25, 45, 65),
    condition = factor(c("A", "B", "A", "B")),
    run = 1L
  )
  sframe <- fmrihrf::sampling_frame(blocklens = 80, TR = 1)
  Y <- matrix(rnorm(80 * 2), 80, 2)
  dset <- matrix_frame(Y, TR = 1, run_length = 80, event_table = etab)
  emod <- event_model(onset ~ hrf(condition), data = etab, block = ~ run,
                      sampling_frame = sframe)
  bmod <- baseline_model(basis = "constant", sframe = sframe)
  fmri_model(emod, bmod, dset)
}

test_that("R3-1: with_package checks every package it is given", {
  expect_silent(fmrireg:::with_package(c("stats", "utils")))
  expect_error(
    fmrireg:::with_package(c("stats", "fmrireg.no.such.pkg.r3")),
    "fmrireg.no.such.pkg.r3"
  )
  err <- tryCatch(
    fmrireg:::with_package(c("fmrireg.no.such.a", "stats", "fmrireg.no.such.b")),
    error = conditionMessage
  )
  expect_match(err, "fmrireg.no.such.a", fixed = TRUE)
  expect_match(err, "fmrireg.no.such.b", fixed = TRUE)
  expect_no_match(err, "stats", fixed = TRUE)
})

test_that("R3-1: design_plot names a missing Suggests package", {
  fmod <- .r3_fmri_model()
  for (missing_pkg in c("shiny", "plotly", "bslib")) {
    local({
      local_mocked_bindings(
        .pkg_available = function(pkg) pkg != missing_pkg,
        .package = "fmrireg"
      )
      expect_error(design_plot(fmod), paste0("'", missing_pkg, "'"),
                   fixed = TRUE)
      expect_error(design_plot(fmod), "design_plot", fixed = TRUE)
    })
  }
})

test_that("R3-1: design_plot still builds the app when the packages exist", {
  skip_if_not_installed("shiny")
  skip_if_not_installed("plotly")
  skip_if_not_installed("bslib")
  app <- design_plot(.r3_fmri_model())
  expect_s3_class(app, "shiny.appobj")
})

test_that("R3-1: design_plot draws lines with linewidth, not the deprecated size", {
  body_txt <- paste(deparse(body(design_plot)), collapse = "\n")
  expect_match(body_txt, "geom_line(linewidth = line_size", fixed = TRUE)
  expect_no_match(body_txt, "geom_line(size", fixed = TRUE)
})

test_that("R3-1: group_data names fmrigds when it is unavailable", {
  local_mocked_bindings(
    .pkg_available = function(pkg) pkg != "fmrigds",
    .package = "fmrireg"
  )
  expect_error(
    group_data(data.frame(subject = "s1", beta = 1, se = 1), format = "csv"),
    "'fmrigds'", fixed = TRUE
  )
})

test_that("R3-1: optional packages are Suggests, not Imports", {
  desc <- read.dcf(system.file("DESCRIPTION", package = "fmrireg"))
  field_pkgs <- function(field) {
    if (!field %in% colnames(desc)) return(character(0))
    x <- strsplit(desc[1, field], ",")[[1]]
    trimws(sub("\\(.*$", "", trimws(x)))
  }
  imports <- field_pkgs("Imports")
  suggests <- field_pkgs("Suggests")
  for (p in c("fmrigds", "plotly", "bslib")) {
    expect_false(p %in% imports, info = p)
    expect_true(p %in% suggests, info = p)
  }
  expect_true("generics" %in% imports)
  # nothing imported from the moved packages, nor the unused aes_string
  ns_imports <- getNamespaceImports("fmrireg")
  expect_false(any(c("plotly", "bslib", "fmrigds") %in% names(ns_imports)))
  expect_false("aes_string" %in% unlist(ns_imports[names(ns_imports) == "ggplot2"]))
})
