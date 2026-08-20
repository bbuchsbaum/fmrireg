test_that("performance guardrail checks partition and retained-state contracts", {
  source(validation_script("fmri_lm_performance_guardrail.R"), local = TRUE)
  out <- run_fmri_lm_performance_guardrail(
    n = 70L, nvox = 12L, chunks = c(1L, 2L), seed = 8111L,
    max_fit_to_data_ratio = 100,
    max_chunk_slowdown = 100
  )

  expect_s3_class(out, "fmri_lm_performance_guardrail")
  expect_named(out$checks, c("partition_invariant", "transient_state_absent",
                             "retained_size_bounded", "chunk_slowdown_bounded"))
  expect_true(out$checks[["partition_invariant"]])
  expect_true(out$checks[["transient_state_absent"]])
})
