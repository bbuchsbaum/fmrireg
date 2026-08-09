validation_script <- function(file) {
  path <- test_path("..", "..", "inst", "validation", file)
  if (!file.exists(path)) {
    path <- system.file("validation", file, package = "fmrireg")
  }
  if (!nzchar(path) || !file.exists(path)) {
    stop(sprintf("Validation script not found: %s", file), call. = FALSE)
  }
  path
}
