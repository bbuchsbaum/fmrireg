test_that("staged paths are made relative on POSIX and Windows", {
  expect_identical(
    fmrireg:::.relative_staged_paths(
      c("/tmp/output/.tmp_write_1/a.h5", "/tmp/output/.tmp_write_1/sub/b.json"),
      "/tmp/output/.tmp_write_1"
    ),
    c("a.h5", "sub/b.json")
  )

  expect_identical(
    fmrireg:::.relative_staged_paths(
      c(
        "C:\\Users\\runner\\out\\.tmp_write_1\\a.h5",
        "C:\\Users\\runner\\out\\.tmp_write_1\\sub\\b.json"
      ),
      "C:\\Users\\runner\\out\\.tmp_write_1"
    ),
    c("a.h5", "sub/b.json")
  )
})

test_that("staged paths outside the write directory are rejected", {
  expect_error(
    fmrireg:::.relative_staged_paths(
      "C:\\Users\\runner\\out\\.tmp_write_10\\a.h5",
      "C:\\Users\\runner\\out\\.tmp_write_1"
    ),
    "outside the staged write directory"
  )
})

test_that("result paths are updated without regular-expression matching", {
  files <- list(
    h5 = "C:\\Users\\runner\\out\\.tmp_write_1\\a.h5",
    nested = list(json = "C:\\Users\\runner\\out\\.tmp_write_1\\sub\\b.json"),
    external = "C:\\Users\\runner\\other\\keep.txt"
  )

  updated <- fmrireg:::.update_file_paths_in_results(
    files,
    "C:\\Users\\runner\\out\\.tmp_write_1",
    "C:\\Users\\runner\\out"
  )

  expect_identical(updated$h5, file.path("C:\\Users\\runner\\out", "a.h5"))
  expect_identical(
    updated$nested$json,
    file.path("C:\\Users\\runner\\out", "sub/b.json")
  )
  expect_identical(updated$external, files$external)
})
