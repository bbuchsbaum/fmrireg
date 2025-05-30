# Benchmark script for evaluate.regressor C++ methods (fft vs conv)

# --- Libraries ---
library(fmrireg)
library(microbenchmark)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)

# Increase print output width for tibbles
options(width = 120)

# --- Parameters ---
param_grid <- expand.grid(
  len = c(500, 1000, 2000, 4000), # Increased max length
  num_events = c(50, 200, 500, 1000), # Increased max events
  duration_type = c("constant_0", "constant_1", "variable"), # Event duration types
  nbasis = c(1, 2, 3),           # Number of HRF basis functions
  stringsAsFactors = FALSE
)

# Filter out impossible combinations (more likely now with increased num_events)
param_grid <- param_grid |>
  filter(num_events <= len)

# Add HRF objects based on nbasis
hrf_list <- list(HRF_SPMG1, HRF_SPMG2, HRF_SPMG3)
param_grid <- param_grid |>
  mutate(hrf = map(nbasis, ~hrf_list[[.x]]))

# TR for sampling frame
TR <- 2

# Number of benchmark repetitions
BENCH_TIMES <- 10

# --- Benchmark Function ---
run_benchmark <- function(len, num_events, duration_type, nbasis, hrf, ...) {
  cat(sprintf("Running: len=%d, events=%d, duration=%s, nbasis=%d\n",
              len, num_events, duration_type, nbasis))

  # Create sampling frame
  sf <- sampling_frame(len, TR = TR)
  grid <- samples(sf)

  # Generate onsets (ensure they are within reasonable bounds)
  max_time <- len * TR
  # Reduce max onset time further to avoid issues with span near the end
  # Use attr() to access S3 attribute 'span'
  hrf_span_val <- attr(hrf, "span") 
  if (is.null(hrf_span_val)) { 
      # Fallback or default if span attribute isn't found (should not happen with HRF_SPMG*)
      hrf_span_val <- 24 
      warning("Could not find 'span' attribute for HRF, using default 24s.")
  }
  max_onset_time <- max(0, max_time - (hrf_span_val * 2)) 
  if (max_onset_time == 0) {
      warning(sprintf("Series length too short for span (%d * TR = %d vs span = %f). Skipping.", len, max_time, hrf_span_val))
      return(NULL)
  }
  onsets <- sort(runif(num_events, 0, max_onset_time))

  # Generate durations
  durations <- switch(duration_type,
    "constant_0" = rep(0, num_events),
    "constant_1" = rep(1, num_events),
    "variable" = runif(num_events, 0, 5),
    stop("Invalid duration_type")
  )

  # Generate amplitudes
  amplitudes <- rep(1, num_events)

  
  # Create regressor
  reg <- tryCatch({
      regressor(onsets = onsets, hrf = hrf, duration = durations, amplitude = amplitudes)
  }, error = function(e) {
      warning(sprintf("Failed to create regressor for params: %s. Skipping. Error: %s",
                      paste(list(len=len, events=num_events, duration=duration_type, nbasis=nbasis), collapse=", "),
                      e$message))
      return(NULL)
  })

  if (is.null(reg)) {
      return(NULL) # Skip if regressor creation failed
  }

  # Clear the C++ FFT cache before this benchmark run
  # clear_fft_cache() # Removed as cache was removed from C++

  # Run microbenchmark - Only fft and conv
  mb_results <- tryCatch({
      microbenchmark(
          fft = evaluate(reg, grid, method = "fft"),
          conv = evaluate(reg, grid, method = "conv"),
          # Rconv = evaluate(reg, grid, method = "Rconv"), # Removed
          # loop = evaluate(reg, grid, method = "loop"),   # Removed
          times = BENCH_TIMES,
          unit = "ms"
      )
  }, error = function(e) {
       warning(sprintf("Microbenchmark failed for params: %s. Skipping. Error: %s",
                      paste(list(len=len, events=num_events, duration=duration_type, nbasis=nbasis), collapse=", "),
                      e$message))
       return(NULL) # Skip if benchmark itself fails
  })

  if (is.null(mb_results)) {
      return(NULL)
  }

  # Augment results with parameters
  mb_summary <- summary(mb_results) |>
      mutate(
          len = len,
          num_events = num_events,
          duration_type = duration_type,
          nbasis = nbasis
          # rconv_eligible removed
      ) |>
      rename(method = expr) |> # Rename expr column to method
      droplevels() # Precaution: drop unused factor levels if any exist

  return(mb_summary)
}

# --- Run Benchmarks ---
cat("Starting benchmarks...\n")
all_results <- pmap_dfr(param_grid, run_benchmark)
cat("Benchmarks finished.\n")

# --- Process Results ---

# --- Debugging Step: Inspect unique methods found ---
cat("\nUnique methods found in raw results:\n")
print(unique(all_results$method))
cat("---------------------------------------------\n\n")

# Filter out any potential NA methods before final processing
all_results_filtered <- all_results |> filter(!is.na(method))

if (nrow(all_results_filtered) == 0) {
    stop("No valid benchmark results were generated after filtering NAs. Check warnings.")
}

# Calculate median time in ms
results_processed <- all_results_filtered |>
  select(len, num_events, duration_type, nbasis, method, median) |>
  mutate(method = factor(method, levels = c("fft", "conv"))) # Order methods

if (any(is.na(results_processed$method))) {
    warning("NA values were introduced in the 'method' column during factor conversion. This indicates unexpected method names were present in the results.")
    print("Problematic rows (before final factor conversion):")
    print(all_results_filtered |> filter(!(method %in% c("fft", "conv"))))
}

if (nrow(results_processed) == 0) {
    stop("No benchmark results remaining after processing. Check warnings.")
}

# --- Plotting ---
cat("Generating plot...\n")

# Plot: Time vs Number of Events, faceted by Length and Duration Type, colored by Method
p_events <- ggplot(results_processed, aes(x = num_events, y = median, color = method, group = method)) +
  geom_line() +
  geom_point() + # Removed shape aesthetic
  # scale_shape_manual removed
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000), # Extended breaks
                labels = scales::comma) +
  scale_x_continuous(breaks = unique(results_processed$num_events)) +
  facet_grid(nbasis ~ len + duration_type, labeller = label_both, scales = "free_y") +
  labs(
    title = "C++ evaluate.regressor Benchmark (fft vs conv): Time vs Number of Events",
    subtitle = "Faceted by HRF Basis Functions, Series Length, and Duration Type",
    x = "Number of Events",
    y = "Median Execution Time (ms, log scale)",
    color = "Method"
  ) +
  theme_bw(base_size = 11) + # Slightly larger base size
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.background = element_rect(fill = "grey90"),
    panel.spacing = unit(1, "lines")
  )

# Plot: Time vs Length, faceted by Num Events and Duration Type, colored by Method
p_length <- ggplot(results_processed, aes(x = len, y = median, color = method, group = method)) +
  geom_line() +
  geom_point() + # Removed shape aesthetic
  # scale_shape_manual removed
  scale_y_log10(breaks = c(1, 10, 100, 1000, 10000, 100000, 1000000), # Extended breaks
                labels = scales::comma) +
  scale_x_continuous(breaks = unique(results_processed$len)) +
  facet_grid(nbasis ~ num_events + duration_type, labeller = label_both, scales = "free_y") +
  labs(
    title = "C++ evaluate.regressor Benchmark (fft vs conv): Time vs Series Length",
    subtitle = "Faceted by HRF Basis Functions, Number of Events, and Duration Type",
    x = "Series Length (Number of Scans)",
    y = "Median Execution Time (ms, log scale)",
    color = "Method"
  ) +
  theme_bw(base_size = 11) + # Slightly larger base size
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    strip.background = element_rect(fill = "grey90"),
    panel.spacing = unit(1, "lines")
  )

# --- Save Plot ---
plot_file <- "bench/evaluate_regressor_cpp_benchmark.pdf" # Changed filename slightly
cat(sprintf("Saving plot to %s...\n", plot_file))
pdf(plot_file, width = 18, height = 10) # Increased width slightly for more facets
print(p_events)
print(p_length)
dev.off()

cat("Benchmark script complete.\n") 