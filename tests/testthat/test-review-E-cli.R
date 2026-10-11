# Review 2026-10, stream E8: install_cli() must not default to writing into
# the user's home directory (CRAN policy).

test_that("E8: install_cli requires an explicit dest_dir", {
  no_default <- identical(formals(install_cli)$dest_dir, quote(expr = ))
  expect_true(no_default)
  # Only exercise the no-argument call once it can no longer write to ~.
  if (no_default) {
    expect_error(install_cli(), "dest_dir")
  }
  expect_error(install_cli(dest_dir = NULL), "dest_dir")
  expect_error(install_cli(dest_dir = c("a", "b")), "dest_dir")
  expect_error(install_cli(dest_dir = ""), "dest_dir")
})

test_that("E8: install_cli still installs into an explicit directory", {
  bin_dir <- tempfile("bin_")
  on.exit(unlink(bin_dir, recursive = TRUE), add = TRUE)
  paths <- suppressMessages(install_cli(bin_dir))
  expect_length(paths, 1L)
  expect_true(file.exists(file.path(bin_dir, "fmrireg")))
})
