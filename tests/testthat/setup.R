# Test setup file for fmrireg
# This file is automatically sourced by testthat before running tests

# Suppress deprecation warnings during tests
# Tests intentionally use deprecated functions to verify parity with new implementations
options(fmrireg.suppress_deprecation = TRUE)

# Restore option on exit
withr::defer(options(fmrireg.suppress_deprecation = FALSE), teardown_env())
