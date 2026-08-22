# Benchmark Datasets

## Introduction

How do you know whether your HRF fitting or beta estimation method
actually works? You need data where the ground truth is known — specific
effect sizes, noise levels, and HRF shapes — so you can measure how well
your method recovers them.

The `fmrireg` package ships a set of benchmark datasets designed for
exactly this purpose. The canonical high- and low-SNR datasets include
complete noiseless condition-beta oracles. The remaining datasets expose
narrower truths—such as HRF-group membership or trial amplitudes—and
record their limits in `oracle_contract` rather than pretending every
scenario supports the same exact check.

## Available Benchmark Datasets

Let’s start by exploring what benchmark datasets are available:

``` r

# List all available datasets
datasets_info <- list_benchmark_datasets()
print(datasets_info)
#>                                                         Dataset
#> BM_Canonical_HighSNR                       BM_Canonical_HighSNR
#> BM_Canonical_LowSNR                         BM_Canonical_LowSNR
#> BM_HRF_Variability_AcrossVoxels BM_HRF_Variability_AcrossVoxels
#> BM_Trial_Amplitude_Variability   BM_Trial_Amplitude_Variability
#> BM_Complex_Realistic                       BM_Complex_Realistic
#>                                                                                                                        Description
#> BM_Canonical_HighSNR                                 Canonical HRF (SPMG1), high SNR, 3 conditions, fixed amplitudes per condition
#> BM_Canonical_LowSNR                                   Canonical HRF (SPMG1), low SNR, 3 conditions, fixed amplitudes per condition
#> BM_HRF_Variability_AcrossVoxels                                         HRF varies across voxel groups, 2 conditions, moderate SNR
#> BM_Trial_Amplitude_Variability                              Single condition with significant trial-to-trial amplitude variability
#> BM_Complex_Realistic            Complex realistic scenario: 3 HRF groups, 3 conditions, variable durations/amplitudes, AR(2) noise
```

## Loading and Exploring a Dataset

Let’s load the high SNR canonical dataset and explore its structure:

``` r

# Load the high SNR dataset
data <- load_benchmark_dataset("BM_Canonical_HighSNR")

# Get a summary of the dataset
summary_info <- get_benchmark_summary("BM_Canonical_HighSNR")
print(summary_info$dimensions)
#> $n_timepoints
#> [1] 150
#> 
#> $n_voxels
#> [1] 100
#> 
#> $n_events
#> [1] 45
#> 
#> $n_conditions
#> [1] 3
print(summary_info$experimental_design)
#> $conditions
#> [1] "Cond1" "Cond2" "Cond3"
#> 
#> $events_per_condition
#> $events_per_condition$Cond1
#> [1] 15
#> 
#> $events_per_condition$Cond2
#> [1] 15
#> 
#> $events_per_condition$Cond3
#> [1] 15
#> 
#> 
#> $TR
#> [1] 2
#> 
#> $total_time
#> [1] 300
#> 
#> $run_length
#> [1] 150
#> 
#> $target_snr
#> [1] 4
```

## Examining the Data Structure

Each benchmark dataset is a list. Key components include:

- `description`: A text summary.
- `Y_noisy`: The matrix of noisy BOLD time series (time points x
  voxels).
- `Y_clean`: (When available) The BOLD signal without noise.
- `event_onsets`: Vector of event start times.
- `condition_labels`: Vector of condition names for each event.
- `event_durations`: Vector of event durations.
- `true_betas_condition`: Ground truth beta values for each condition.
- `true_hrf_parameters`: The reconstructed HRF and its complete
  generation recipe (base HRF, lag, width, normalization, and `fmrihrf`
  version).
- `oracle_contract`: Which stored objects form a valid numerical oracle,
  or why only partial ground truth is available.
- `TR`, `total_time`, `run_length`: Scan parameters.

``` r

# Look at the BOLD time series dimensions and event structure
cat("Y_noisy BOLD data dimensions:", dim(data$Y_noisy), "\n")
#> Y_noisy BOLD data dimensions: 150 100
cat("Number of events:", length(data$event_onsets), "\n")
#> Number of events: 45
cat("Conditions:", unique(data$condition_labels), "\n")
#> Conditions: Cond1 Cond2 Cond3
cat("Events per condition:", table(data$condition_labels), "\n")
#> Events per condition: 15 15 15
cat("TR:", data$TR, "\n")
#> TR: 2
cat("Run length:", data$run_length, "\n")
#> Run length: 150
```

## Visualizing the Data

Let’s visualize some aspects of the benchmark dataset:

``` r

# Plot the first few voxels' time series
n_timepoints <- nrow(data$Y_noisy)
time_points <- seq(0, by = data$TR, length.out = n_timepoints)

# Create a data frame for plotting
plot_data <- data.frame(
  Time = rep(time_points, 3),
  BOLD = c(data$Y_noisy[, 1], data$Y_noisy[, 2], data$Y_noisy[, 3]),
  Voxel = rep(paste("Voxel", 1:3), each = n_timepoints)
)

# Build a separate event carpet so dense event markers do not obscure the
# response traces.
event_data <- data.frame(
  Time = data$event_onsets,
  Condition = data$condition_labels
)

p_signal <- ggplot(plot_data, aes(x = Time, y = BOLD)) +
  geom_line(linewidth = 0.35) +
  facet_wrap(~Voxel, scales = "free_y") +
  labs(title = "BOLD Time Series",
       x = "Time (seconds)", y = "BOLD Signal") +
  theme_minimal()

p_events <- ggplot(event_data, aes(x = Time, y = Condition, color = Condition)) +
  geom_point(shape = "|", size = 7, stroke = 1.1) +
  scale_x_continuous(limits = range(time_points), expand = c(0, 0)) +
  labs(title = "Event carpet", x = "Time (seconds)", y = NULL) +
  theme_minimal() +
  theme(legend.position = "none")

gridExtra::grid.arrange(p_signal, p_events, ncol = 1, heights = c(4, 1.3))
```

![Three representative BOLD traces aligned above an event carpet whose
colored rows identify Cond1, Cond2, and Cond3 onsets across the
run.](benchmark_datasets_files/figure-html/visualize_data-1.png)

## Creating Design Matrices

One of the key features is the ability to create design matrices with
different HRF assumptions:

``` r

# Reconstruct the exact generating HRF from the stored recipe.
true_hrf <- data$true_hrf_parameters$hrf_object
X_true <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", true_hrf)

# Create design matrix with a different HRF (e.g., a Gaussian HRF instead of canonical)
X_wrong <- create_design_matrix_from_benchmark("BM_Canonical_HighSNR", fmrihrf::HRF_GAUSSIAN)

cat("True HRF design matrix dimensions:", dim(X_true), "\n")
#> True HRF design matrix dimensions: 150 4
cat("Alternative HRF design matrix dimensions:", dim(X_wrong), "\n")
#> Alternative HRF design matrix dimensions: 150 4
```

``` r

recipe <- data$true_hrf_parameters$hrf_recipe
knitr::kable(data.frame(
  Base_HRF = recipe$base_hrf_name,
  Lag = recipe$lag,
  Width = recipe$width,
  Normalize = recipe$normalize,
  fmrihrf_version = recipe$fmrihrf_version
), caption = "Reconstructable HRF provenance for this benchmark")
```

| Base_HRF  | Lag | Width | Normalize | fmrihrf_version |
|:----------|----:|------:|:----------|:----------------|
| HRF_SPMG1 |   0 |     0 | FALSE     | 0.4.0           |

Reconstructable HRF provenance for this benchmark {.table}

## Does the Benchmark Recover Its Own Ground Truth?

A benchmark should first pass an exact oracle: with no noise and the HRF
used to generate the signal, ordinary least squares should recover the
declared betas. This check catches design scaling, time-grid, and
ground-truth errors before the dataset is used to compare methods.

``` r

# Exact noiseless oracle
betas_oracle <- qr.solve(X_true, data$Y_clean)
oracle_max_error <- max(abs(betas_oracle[-1, ] - data$true_betas_condition))

# Fit ordinary least squares with the correct HRF
betas_correct <- qr.solve(X_true, data$Y_noisy)

# Fit OLS with the wrong HRF assumption
betas_wrong <- qr.solve(X_wrong, data$Y_noisy)

# Evaluate performance (remove intercept for comparison)
performance_correct <- evaluate_method_performance("BM_Canonical_HighSNR", 
                                                   betas_correct[-1, ], 
                                                   "OLS_Correct_HRF")

performance_wrong <- evaluate_method_performance("BM_Canonical_HighSNR",
                                                 betas_wrong[-1, ], 
                                                 "OLS_Wrong_HRF")

performance_table <- data.frame(
  HRF = c("Correct (SPMG1)", "Wrong (Gaussian)"),
  Correlation = c(
    performance_correct$overall_metrics$correlation,
    performance_wrong$overall_metrics$correlation
  ),
  RMSE = c(
    performance_correct$overall_metrics$rmse,
    performance_wrong$overall_metrics$rmse
  ),
  MAE = c(
    performance_correct$overall_metrics$mae,
    performance_wrong$overall_metrics$mae
  )
)
```

``` r

cat("Noiseless oracle maximum absolute beta error:",
    format(oracle_max_error, scientific = TRUE, digits = 3), "\n")
#> Noiseless oracle maximum absolute beta error: 3.55e-15

knitr::kable(
  performance_table,
  digits = 3,
  caption = "Noisy-data beta recovery under correct and misspecified HRFs"
)
```

| HRF              | Correlation |  RMSE |   MAE |
|:-----------------|------------:|------:|------:|
| Correct (SPMG1)  |        0.99 | 0.042 | 0.034 |
| Wrong (Gaussian) |        0.99 | 8.397 | 7.999 |

Noisy-data beta recovery under correct and misspecified HRFs {.table}

The noiseless error is at floating-point precision, establishing that
the stored design, response, and beta truth use the same scale. On noisy
data, the Gaussian HRF can retain a high correlation while producing
much larger RMSE and MAE. Correlation preserves ordering but is
insensitive to absolute scale; the scale-sensitive errors reveal the
misspecification.

## Comparing True vs Estimated Betas

``` r

# Get true betas
true_betas <- data$true_betas_condition

# Create comparison plots
comparison_data <- dplyr::bind_rows(lapply(seq_len(nrow(true_betas)), function(i) {
  data.frame(
    True = true_betas[i, ],
    Estimated_Correct = betas_correct[-1, ][i, ],
    Estimated_Wrong = betas_wrong[-1, ][i, ],
    Condition = paste("Condition", i)
  )
}))

# Plot true vs estimated (correct HRF)
p1 <- ggplot(comparison_data, aes(x = True, y = Estimated_Correct, color = Condition)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Correct HRF", x = "True Beta", y = "Estimated Beta") +
  theme_minimal()

# Plot true vs estimated (wrong HRF)
p2 <- ggplot(comparison_data, aes(x = True, y = Estimated_Wrong, color = Condition)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(title = "Wrong HRF", x = "True Beta", y = "Estimated Beta") +
  theme_minimal()

# Display plots side by side
gridExtra::grid.arrange(p1, p2, ncol = 2)
```

![Side-by-side true-versus-estimated beta plots: correct-HRF estimates
lie near the identity line, whereas Gaussian-HRF estimates are severely
inflated in
scale.](benchmark_datasets_files/figure-html/compare_betas-1.png)

## Testing Different Datasets

Let’s compare performance across different benchmark scenarios:

``` r

# Test on different datasets
datasets_to_test <- c("BM_Canonical_HighSNR", "BM_Canonical_LowSNR")
results <- list()

for (dataset_name in datasets_to_test) {
  # Load dataset and create design matrix
  data_test <- load_benchmark_dataset(dataset_name)
  X <- create_design_matrix_from_benchmark(
    dataset_name,
    data_test$true_hrf_parameters$hrf_object
  )
  
  # Fit model
  betas <- qr.solve(X, data_test$Y_noisy)
  oracle_betas <- qr.solve(X, data_test$Y_clean)[-1, , drop = FALSE]
  oracle_error <- max(abs(oracle_betas - data_test$true_betas_condition))
  
  # Evaluate performance
  perf <- evaluate_method_performance(dataset_name, betas[-1, ], "OLS")
  
  results[[dataset_name]] <- list(
    correlation = perf$overall_metrics$correlation,
    rmse = perf$overall_metrics$rmse,
    oracle_error = oracle_error,
    target_snr = data_test$target_snr
  )
}

# Display results
results_df <- data.frame(
  Dataset = names(results),
  Correlation = sapply(results, function(x) round(x$correlation, 3)),
  RMSE = sapply(results, function(x) round(x$rmse, 3)),
  Noiseless_Max_Error = sapply(results, function(x) x$oracle_error),
  Target_SNR = sapply(results, function(x) x$target_snr)
)

print(results_df)
#>                                   Dataset Correlation  RMSE Noiseless_Max_Error
#> BM_Canonical_HighSNR BM_Canonical_HighSNR       0.990 0.042        3.552714e-15
#> BM_Canonical_LowSNR   BM_Canonical_LowSNR       0.617 0.383        7.549517e-15
#>                      Target_SNR
#> BM_Canonical_HighSNR        4.0
#> BM_Canonical_LowSNR         0.5
```

## HRF Variability Dataset

Let’s explore the dataset with HRF variability across voxels:

``` r

# Load the HRF variability dataset
hrf_data <- load_benchmark_dataset("BM_HRF_Variability_AcrossVoxels")

# Examine the HRF group assignments
cat("HRF group assignments:", table(hrf_data$true_hrf_group_assignment), "\n")
#> HRF group assignments: 50 50

hrf_group_recipes <- do.call(rbind, lapply(
  names(hrf_data$true_hrf_parameters),
  function(group_name) {
    group_recipe <- hrf_data$true_hrf_parameters[[group_name]]$hrf_recipe
    data.frame(
      Group = group_name,
      Recipe = group_recipe$recipe_name,
      Base_HRF = group_recipe$base_hrf_name,
      Lag = group_recipe$lag,
      Width = group_recipe$width,
      Normalize = group_recipe$normalize
    )
  }
))
knitr::kable(hrf_group_recipes, caption = "Exact HRF recipe by voxel group")
```

| Group  | Recipe    | Base_HRF  | Lag | Width | Normalize |
|:-------|:----------|:----------|----:|------:|:----------|
| group1 | HRF_SPMG1 | HRF_SPMG1 |   0 |   0.0 | FALSE     |
| group2 | variant1  | HRF_SPMG1 |   1 |   1.2 | TRUE      |

Exact HRF recipe by voxel group {.table}

This dataset records exact HRF recipes, group assignments, and nominal
condition effects, but it does not ship a clean response/design pair.
Its `oracle_contract` therefore does not claim exact noiseless OLS
recovery.

## Trial Amplitude Variability

Let’s examine the trial-to-trial variability dataset:

``` r

# Load the trial variability dataset
trial_data <- load_benchmark_dataset("BM_Trial_Amplitude_Variability")

# Look at the trial-wise amplitudes
true_trial_amps <- trial_data$true_amplitudes_trial

# Plot amplitude variability across trials for first few voxels
amp_plot_data <- data.frame(
  Trial = rep(1:nrow(true_trial_amps), 3),
  Amplitude = c(true_trial_amps[, 1], true_trial_amps[, 2], true_trial_amps[, 3]),
  Voxel = rep(paste("Voxel", 1:3), each = nrow(true_trial_amps))
)

ggplot(amp_plot_data, aes(x = Trial, y = Amplitude)) +
  geom_line() +
  geom_point() +
  facet_wrap(~Voxel) +
  labs(title = "Trial-to-Trial Amplitude Variability",
       x = "Trial Number", y = "True Amplitude") +
  theme_minimal()
```

![Three voxel facets showing that the stored true amplitude fluctuates
from trial to trial around
one.](benchmark_datasets_files/figure-html/trial_variability-1.png)

## Summary

The package currently ships five synthetic cases with different truth
levels:

1.  **Canonical high and low SNR**: complete clean-response, design, and
    condition-beta oracles
2.  **HRF variability**: HRF recipes and voxel-group assignments,
    without a complete clean-response/design oracle
3.  **Trial-amplitude variability**: per-trial amplitude truth for LSS
    checks
4.  **Complex scenario**: HRF-group, duration, and amplitude metadata
    under AR(2) noise, with the limitations stated in its contract

The exercised contracts are concrete:

- **Explicit oracle scope**: Complete condition-beta recovery is
  available for the two canonical datasets; other datasets state which
  partial truths are valid
- **Exact HRF provenance**: Every declared HRF group has a fail-closed,
  reconstructable recipe
- **Declared noise models**: AR(1) and AR(2) parameters are stored with
  the data
- **Checked dimensions**: responses, stored designs, and declared truth
  matrices agree on time-point and voxel dimensions
- **Deterministic generation**: fixed seeds reproduce the packaged
  artifact

Use only the truth level declared by each dataset’s `oracle_contract`.
That distinction is what makes comparisons reproducible without
overstating what a particular synthetic scenario can prove.

## Next

- [`vignette("group_analysis", package = "fmrireg")`](https://bbuchsbaum.github.io/fmrireg/articles/group_analysis.md)
  — Group analysis
