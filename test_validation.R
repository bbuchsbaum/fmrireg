# Test validation fix
devtools::load_all('.')

# Test the sanitization fix
cat("Testing sanitization:\n")
test1 <- fmrireg:::.sanitize_label("")
cat("Empty string sanitizes to:", test1, "is.null:", is.null(test1), "\n")

test2 <- fmrireg:::.sanitize_label("@#$")
cat("Invalid chars sanitize to:", test2, "is.null:", is.null(test2), "\n")

# Test validation order
cat("Testing validation order:\n")

entities1 <- fmrireg:::.create_bids_entities(subject = NULL, task = "test", space = NULL)
cat("Subject=NULL, Task=test:\n")
cat("  subject:", entities1$subject, "is.null:", is.null(entities1$subject), "\n")
cat("  task:", entities1$task, "is.null:", is.null(entities1$task), "\n")

entities2 <- fmrireg:::.create_bids_entities(subject = "01", task = NULL, space = NULL)
cat("Subject=01, Task=NULL:\n")
cat("  subject:", entities2$subject, "is.null:", is.null(entities2$subject), "\n")
cat("  task:", entities2$task, "is.null:", is.null(entities2$task), "\n")

# Test actual validation
tryCatch({
  fmrireg:::.validate_required_entities(entities1)
}, error = function(e) {
  cat("Validation error for entities1:", e$message, "\n")
})

tryCatch({
  fmrireg:::.validate_required_entities(entities2)
}, error = function(e) {
  cat("Validation error for entities2:", e$message, "\n")
})