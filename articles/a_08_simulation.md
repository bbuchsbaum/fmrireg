# Simulating fMRI Data

## Introduction to fMRI Data Simulation

When checking an fMRI analysis pipeline, it is useful to have data where
the ground truth is known – real fMRI data do not come with known effect
sizes or noise parameters. Simulation supplies controlled examples. The
`fmrireg` package offers several functions with varying levels of
complexity:

1.  **`simulate_bold_signal`**: Simulates clean BOLD responses for
    multiple experimental conditions
2.  **`simulate_noise_vector`**: Generates structured, fMRI-like noise
    with temporal autocorrelation, drift, and simplified physiological
    components
3.  **`simulate_simple_dataset`**: Combines signal and noise for a
    complete dataset based on SNR
4.  **`simulate_fmri_matrix`**: Creates multiple time series with shared
    event timing but column-specific variation in parameters

This vignette demonstrates their contracts and assumptions. The outputs
are teaching simulations, not empirically validated generative models of
scanner data.

## Simulating Clean BOLD Signals

Let’s start with `simulate_bold_signal`, which generates a clean
hemodynamic response signal for multiple experimental conditions.

``` r

# Simulate 3 conditions with different amplitudes
sim <- simulate_bold_signal(ncond = 3, amps = c(1.0, 1.5, 2.0), TR = 2)

# Extract the data
time <- sim$mat[,1]  # First column contains time
signals <- sim$mat[,-1]  # Other columns contain condition signals
```

We reshape the data into long format so ggplot2 can map each condition
to a separate color.

``` r

df <- data.frame(
  Time = time,
  Cond1 = signals[,1],
  Cond2 = signals[,2],
  Cond3 = signals[,3]
)

df_long <- tidyr::pivot_longer(df, cols = c(Cond1, Cond2, Cond3),
                        names_to = "Condition",
                        values_to = "Response")
```

``` r

ggplot(df_long, aes(x = Time, y = Response, color = Condition)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Simulated BOLD Signals for Three Conditions",
       x = "Time (seconds)",
       y = "BOLD Response",
       color = "Condition") +
  scale_color_brewer(palette = "Set1")
```

![Three colored clean BOLD-response curves with distinct condition
amplitudes and event timings over a 180-second
run.](a_08_simulation_files/figure-html/plot-bold-signals-1.png)

Adding dashed vertical lines at each event onset makes the relationship
between stimulus timing and the hemodynamic response easier to see.

``` r

ggplot(df_long, aes(x = Time, y = Response, color = Condition)) +
  geom_line() +
  geom_vline(xintercept = sim$onset, linetype = "dashed", alpha = 0.3) +
  theme_minimal() +
  labs(title = "Simulated BOLD Signals with Event Onsets",
       x = "Time (seconds)",
       y = "BOLD Response",
       color = "Condition") +
  scale_color_brewer(palette = "Set1")
```

![The three clean BOLD-response curves with dashed vertical markers
showing every event
onset.](a_08_simulation_files/figure-html/plot-bold-with-onsets-1.png)

The function returns a list containing:

- **`onset`**: Event onset times
- **`condition`**: Condition labels for each event
- **`mat`**: Matrix with time points and BOLD responses for each
  condition

You can control:

- Number of conditions (`ncond`)
- Number of repetitions per condition (`nreps`)
- HRF shape (`hrf`)
- Amplitudes for each condition (`amps`)
- Inter-stimulus interval range (`isi`)
- Amplitude variability (`ampsd`)

## Simulating Structured fMRI-like Noise

The `simulate_noise_vector` function combines stylized components
motivated by common fMRI noise sources. To understand exactly what it
generates, we simulate and visualize each component separately and then
combine them.

``` r

n_timepoints <- 200
TR <- 2
time <- seq(0, (n_timepoints - 1) * TR, by = TR)
set.seed(123)
```

First, pure white noise – random fluctuations with no temporal
structure.

``` r

noise_white <- simulate_noise_vector(n_timepoints, TR = TR,
                                  ar = numeric(0), ma = numeric(0),
                                  drift_amplitude = 0, physio = FALSE, sd = 1)
```

Next, ARMA noise introduces temporal autocorrelation, which makes the
signal “smoother” than white noise.

``` r

noise_arma <- simulate_noise_vector(n_timepoints, TR = TR,
                                 ar = c(0.6), ma = c(0.3),
                                 drift_amplitude = 0, physio = FALSE, sd = 1)
```

Scanner drift is a slow oscillation that shifts the baseline over the
course of a run.

``` r

drift_freq <- 1/128
drift_amplitude <- 2
noise_drift <- drift_amplitude * sin(2 * pi * drift_freq * time)
```

Physiological noise originates in continuous time, while an fMRI scanner
records only one sample every TR. Here respiration is 0.3 Hz and cardiac
activity is 1.1 Hz. Both exceed the 0.25 Hz Nyquist frequency at TR = 2
seconds, so they appear at lower frequencies after sampling: 0.2 Hz and
0.1 Hz, respectively. Generating the continuous signals first makes that
aliasing explicit.

``` r

physio_time <- seq(0, max(time), by = 0.05)
cardiac_continuous <- 0.5 * sin(2 * pi * 1.1 * physio_time)
respiratory_continuous <- 0.8 * sin(2 * pi * 0.3 * physio_time)

noise_cardiac <- approx(physio_time, cardiac_continuous, xout = time)$y
noise_respiratory <- approx(physio_time, respiratory_continuous, xout = time)$y
noise_physio <- noise_cardiac + noise_respiratory
```

Finally, we combine the ARMA, drift, and sinusoidal physiological
components into the declared structured teaching model. This sum
illustrates separable noise sources; it is not an empirically validated
scanner-noise distribution.

``` r

noise_combined <- noise_arma + noise_drift + noise_physio
```

``` r

noise_df <- data.frame(
  Time = time,
  White_Noise = noise_white,
  ARMA_Noise = noise_arma,
  Drift_Component = noise_drift,
  Physiological_Noise = noise_physio,
  Combined_Noise = noise_combined
)

noise_long <- tidyr::pivot_longer(noise_df,
                                  cols = -Time,
                                  names_to = "NoiseType",
                                  values_to = "Signal")

noise_long$NoiseType <- factor(noise_long$NoiseType,
                              levels = c("White_Noise", "ARMA_Noise",
                                        "Drift_Component", "Physiological_Noise",
                                        "Combined_Noise"),
                              labels = c("White Noise", "ARMA Noise",
                                        "Drift Component", "Physiological Noise",
                                        "Combined Noise"))
```

``` r

ggplot(noise_long, aes(x = Time, y = Signal, color = NoiseType)) +
  geom_line() +
  facet_wrap(~NoiseType, ncol = 1, scales = "free_y") +
  theme_minimal() +
  theme(legend.position = "none") +
  labs(title = "Simulated fMRI Noise Components: Time Courses",
       x = "Time (seconds)",
       y = "Signal Amplitude")
```

![Five stacked time-course panels showing white, ARMA, sinusoidal drift,
sampled physiological, and combined structured-noise
components.](a_08_simulation_files/figure-html/plot-noise-timecourses-1.png)

### Power Spectrum Analysis

Power spectra reveal the frequency content of each noise type. White
noise has a flat spectrum, while ARMA noise concentrates power at lower
frequencies.

``` r

spec_white <- calculate_spectrum(noise_white, TR)
spec_arma <- calculate_spectrum(noise_arma, TR)
spec_drift <- calculate_spectrum(noise_drift, TR)
spec_physio <- calculate_spectrum(noise_physio, TR)
spec_combined <- calculate_spectrum(noise_combined, TR)
```

``` r

spec_white$NoiseType <- "White Noise"
spec_arma$NoiseType <- "ARMA Noise"
spec_drift$NoiseType <- "Drift Component"
spec_physio$NoiseType <- "Physiological Noise"
spec_combined$NoiseType <- "Combined Noise"

spec_df <- rbind(spec_white, spec_arma, spec_drift, spec_physio, spec_combined)

spec_df$NoiseType <- factor(spec_df$NoiseType,
                           levels = c("White Noise", "ARMA Noise",
                                     "Drift Component", "Physiological Noise",
                                     "Combined Noise"))
```

We also compute a high-resolution drift spectrum from a longer signal to
better resolve the low-frequency peak.

``` r

n_long <- 1024
time_long <- seq(0, (n_long - 1) * TR, by = TR)
drift_long <- drift_amplitude * sin(2 * pi * drift_freq * time_long)
spec_drift_long <- calculate_spectrum(drift_long, TR)
spec_drift_long$NoiseType <- factor(
  "Drift Component",
  levels = levels(spec_df$NoiseType)
)
```

Each panel below shows the power spectrum of a different noise
component. The dashed green line in the drift panel shows the
high-resolution version.

``` r

ggplot() +
  geom_line(data = spec_df, aes(x = Frequency, y = Power, color = NoiseType)) +
  geom_line(data = subset(spec_drift_long, Frequency <= 0.05),
           aes(x = Frequency, y = Power), color = "darkgreen", linetype = "dashed") +
  theme_minimal() +
  facet_wrap(~NoiseType, ncol = 1, scales = "free_y") +
  theme(legend.position = "none") +
  labs(title = "Power Spectra of Different Noise Components",
       x = "Frequency (Hz)",
       y = "Power") +
  scale_color_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(0, 0.25))
```

![Five stacked power spectra: white noise is broad, ARMA power favors
low frequencies, drift peaks near 0.008 hertz, and physiological aliases
peak near 0.1 and 0.2
hertz.](a_08_simulation_files/figure-html/plot-spectra-faceted-1.png)

Plotting all spectra on a single log-scaled axis makes it easy to
compare their relative magnitudes.

``` r

ggplot(spec_df, aes(x = Frequency, y = Power, color = NoiseType)) +
  geom_line() +
  theme_minimal() +
  labs(title = "Power Spectra of Different Noise Components (Log Scale)",
       x = "Frequency (Hz)",
       y = "Power (log scale)") +
  scale_y_log10() +
  scale_color_brewer(palette = "Set1") +
  guides(color = guide_legend(nrow = 2, byrow = TRUE)) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 8)
  ) +
  coord_cartesian(xlim = c(0, 0.25))
```

![Overlaid log-scale power spectra for all five noise components, with
the physiological alias peaks visible near 0.1 and 0.2
hertz.](a_08_simulation_files/figure-html/plot-spectra-log-1.png)

Zooming in on the very low-frequency region confirms the drift component
peaks at the expected frequency (marked by the red dotted line).

``` r

drift_freq_idx <- which.min(abs(spec_df$Frequency - drift_freq))

ggplot() +
  geom_line(data = subset(spec_df, NoiseType == "Drift Component"),
           aes(x = Frequency, y = Power), color = "darkgreen") +
  geom_line(data = subset(spec_drift_long, Frequency <= 0.05),
           aes(x = Frequency, y = Power), color = "darkgreen", linetype = "dashed") +
  geom_vline(xintercept = drift_freq, linetype = "dotted", color = "red") +
  theme_minimal() +
  labs(title = "Zoomed Power Spectrum of Drift Component",
       subtitle = paste("Expected peak at", drift_freq, "Hz"),
       x = "Frequency (Hz)",
       y = "Power") +
  coord_cartesian(xlim = c(0, 0.05))
```

![A low-frequency zoom of the drift spectrum, with the resolved dashed
curve peaking beside the expected 0.0078125-hertz
marker.](a_08_simulation_files/figure-html/plot-drift-zoom-1.png)

The simulation shows five distinct types of noise components:

1.  **White Noise**: Random fluctuations with equal power across all
    frequencies (flat power spectrum).

2.  **ARMA Noise**: Temporally autocorrelated noise with a
    characteristic “smoothed” appearance. The power spectrum shows more
    power at lower frequencies.

3.  **Drift Component**: A very slow oscillation typical of scanner
    drift or physiological trends. The power spectrum shows a dominant
    peak at a very low frequency.

4.  **Physiological Noise**: Simplified respiratory and cardiac
    sinusoids at 0.3 Hz and 1.1 Hz are above the 0.25 Hz Nyquist limit.
    At TR = 2 seconds they appear at distinct aliases of 0.2 Hz and 0.1
    Hz. These components demonstrate sampling and aliasing; they are not
    a validated physiological-noise model, nor evidence that the
    underlying physiology occurs at the aliased frequencies.

5.  **Combined Noise**: The sum of the stylized components. Its power
    spectrum shows features from each contributor.

Adding these components to task-related signals creates controlled,
structured time series for demonstrations and numerical checks.

## Creating a Complete Dataset with Signal and Noise

The `simulate_simple_dataset` function combines clean signals and noise
to create a complete fMRI dataset with a specified signal-to-noise ratio
(SNR).

``` r

snr_values <- c(1.0, 0.5, 0.2)
snr_datasets <- lapply(snr_values, function(snr) {
  simulate_simple_dataset(ncond = 3, TR = 2, snr = snr, seed = 42)
})
names(snr_datasets) <- paste0("SNR = ", snr_values)

data_snr_1.0 <- snr_datasets[["SNR = 1"]]
data_snr_0.5 <- snr_datasets[["SNR = 0.5"]]
data_snr_0.2 <- snr_datasets[["SNR = 0.2"]]
```

The same seed deliberately preserves the event schedule, clean signal,
and standardized noise realization. Only the scale of that noise
changes, so the panels isolate SNR rather than comparing three unrelated
simulations.

``` r

plot_df <- rbind(
  create_plot_df(data_snr_1.0, "SNR = 1.0"),
  create_plot_df(data_snr_0.5, "SNR = 0.5"),
  create_plot_df(data_snr_0.2, "SNR = 0.2")
)

plot_df_long <- tidyr::pivot_longer(plot_df,
                                    cols = -c(Time, SNR),
                                    names_to = "Type",
                                    values_to = "Signal")
```

The overlay plot shows the clean signal (red) against the noisy
measurement (blue) at each SNR level. As SNR decreases, the noise
increasingly obscures the true signal.

``` r

ggplot(plot_df_long, aes(x = Time, y = Signal, color = Type)) +
  geom_line() +
  facet_wrap(~SNR, ncol = 1) +
  theme_minimal() +
  labs(title = "Clean vs. Noisy Signals at Different SNR Levels",
       subtitle = "Single condition shown for clarity",
       x = "Time (seconds)",
       y = "Signal",
       color = "Signal Type") +
  scale_color_manual(values = c("Clean" = "darkred", "Noisy" = "steelblue")) +
  theme(legend.position = "top")
```

![Three stacked overlays compare one shared clean signal with noisy
observations at SNR 0.2, 0.5, and 1.0; masking decreases as SNR
rises.](a_08_simulation_files/figure-html/plot-snr-overlay-1.png)

The decomposition plots below separate each SNR level into its signal,
noise, and combined components. At SNR = 1.0 the signal is clearly
visible; by SNR = 0.2 the noise dominates.

``` r

plot_faceted(plot_df, "SNR = 1.0")
```

![Three stacked panels decompose the SNR 1.0 simulation into signal plus
noise, clean signal, and
noise.](a_08_simulation_files/figure-html/plot-decomposition-snr1-1.png)

``` r

plot_faceted(plot_df, "SNR = 0.5")
```

![Three stacked panels decompose the SNR 0.5 simulation into signal plus
noise, clean signal, and
noise.](a_08_simulation_files/figure-html/plot-decomposition-snr05-1.png)

``` r

plot_faceted(plot_df, "SNR = 0.2")
```

![Three stacked panels decompose the SNR 0.2 simulation into signal plus
noise, clean signal, and dominant
noise.](a_08_simulation_files/figure-html/plot-decomposition-snr02-1.png)

``` r

snr_stats <- data.frame(
  Requested_SNR = snr_values,
  Empirical_SNR = unname(empirical_snr),
  Signal_SD = vapply(snr_datasets, function(dataset) {
    sd(as.vector(dataset$clean$mat[, -1]))
  }, numeric(1)),
  Noise_SD = vapply(snr_datasets, function(dataset) {
    sd(as.vector(dataset$noise))
  }, numeric(1))
)
```

``` r

knitr::kable(
  snr_stats,
  digits = 3,
  caption = "Requested and empirically realized SNR from a controlled simulation"
)
```

|           | Requested_SNR | Empirical_SNR | Signal_SD | Noise_SD |
|:----------|--------------:|--------------:|----------:|---------:|
| SNR = 1   |           1.0 |           1.0 |     0.861 |    0.861 |
| SNR = 0.5 |           0.5 |           0.5 |     0.861 |    1.722 |
| SNR = 0.2 |           0.2 |           0.2 |     0.861 |    4.306 |

Requested and empirically realized SNR from a controlled simulation
{.table}

This visualization shows how different SNR levels affect the fMRI time
series. The lower the SNR, the more the noise dominates the signal. For
each SNR level, we show:

1.  **The original clean signal** (red line): The true underlying BOLD
    response
2.  **The noisy signal** (blue line): The synthetic observation produced
    by this model
3.  **Signal decomposition**: Visualization of how the signal, noise,
    and combined signal relate at each SNR level

With SNR = 1.0, the signal pattern remains visible despite the noise. At
SNR = 0.5, some features are obscured; at SNR = 0.2, the noise
substantially masks the signal in this controlled example.

The function returns:

- **`clean`**: The simulated signals without noise
- **`noisy`**: The signals with added noise
- **`noise`**: The simulated noise component
- **`onsets`**: Trial onset times
- **`conditions`**: Condition labels for each trial

## Simulating Matrix Time Series with Column-Specific Variation

The most advanced simulation function, `simulate_fmri_matrix`, generates
multiple time series (columns) that share the same event timing but have
column-specific variations in amplitudes and durations. This is
particularly useful for simulating multiple voxels or regions with
related but slightly different response profiles.

``` r

requested_onsets <- seq(12, 156, length.out = 10)
sim_matrix <- simulate_fmri_matrix(
  n = 5,                  # 5 voxels/regions
  total_time = 200,       # 200 seconds of scan time
  TR = 2,                 # TR = 2 seconds
  onsets = requested_onsets,
  amplitudes = 1,         # Base amplitude = 1
  amplitude_sd = 0.3,     # Amplitude variability
  durations = 2,          # Base duration = 2 seconds
  duration_sd = 0.5,      # Duration variability
  noise_type = "ar1",     # AR(1) noise
  noise_sd = 0.5          # Noise standard deviation
)
```

When you supply `onsets`, those values are authoritative: the function
skips ISI generation and derives the event count from their length.

``` r

ts_data <- sim_matrix$time_series
matrix_data <- ts_data$datamat

time_points <- seq(0, by = 2, length.out = nrow(matrix_data))
plot_data <- data.frame(Time = time_points)

for(i in 1:ncol(matrix_data)) {
  plot_data[[paste0("Voxel", i)]] <- matrix_data[, i]
}

plot_data_long <- tidyr::pivot_longer(
  plot_data,
  cols = starts_with("Voxel"),
  names_to = "Voxel",
  values_to = "Signal"
)
```

All five simulated voxels share the same event timing, but their
response amplitudes and durations vary independently.

``` r

ggplot(plot_data_long, aes(x = Time, y = Signal, color = Voxel)) +
  geom_line(alpha = 0.8) +
  theme_minimal() +
  labs(title = "Multiple Simulated Time Series with Column-Specific Variation",
       x = "Time (seconds)",
       y = "Signal",
       color = "Voxel") +
  scale_color_brewer(palette = "Set2")
```

![Five colored simulated voxel time series share event timing but differ
in response amplitude, duration, and independent
noise.](a_08_simulation_files/figure-html/plot-matrix-timeseries-1.png)

``` r

amp_df <- as.data.frame(sim_matrix$ampmat)
colnames(amp_df) <- paste0("Voxel", 1:ncol(amp_df))
amp_df$Event <- 1:nrow(amp_df)

dur_df <- as.data.frame(sim_matrix$durmat)
colnames(dur_df) <- paste0("Voxel", 1:ncol(dur_df))
dur_df$Event <- 1:nrow(dur_df)

amp_long <- tidyr::pivot_longer(
  amp_df,
  cols = starts_with("Voxel"),
  names_to = "Voxel",
  values_to = "Amplitude"
)

dur_long <- tidyr::pivot_longer(
  dur_df,
  cols = starts_with("Voxel"),
  names_to = "Voxel",
  values_to = "Duration"
)
```

The amplitude plot shows how each voxel’s response magnitude varies from
event to event around the base amplitude of 1.

``` r

ggplot(amp_long, aes(x = Event, y = Amplitude, color = Voxel, group = Voxel)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Amplitude Variation Across Events and Voxels",
       x = "Event Number",
       y = "Amplitude",
       color = "Voxel") +
  scale_color_brewer(palette = "Set2")
```

![Five colored lines show independently sampled response amplitudes
across ten events for five
voxels.](a_08_simulation_files/figure-html/plot-amplitude-variation-1.png)

Duration variation follows the same pattern – each voxel draws its event
durations independently around the base of 2 seconds.

``` r

ggplot(dur_long, aes(x = Event, y = Duration, color = Voxel, group = Voxel)) +
  geom_line() +
  geom_point() +
  theme_minimal() +
  labs(title = "Duration Variation Across Events and Voxels",
       x = "Event Number",
       y = "Duration (seconds)",
       color = "Voxel") +
  scale_color_brewer(palette = "Set2")
```

![Five colored lines show independently sampled response durations
across ten events for five
voxels.](a_08_simulation_files/figure-html/plot-duration-variation-1.png)

This function exposes the following controls for multiple related time
series:

- **Shared event timing** but individual variation in:
  - **Amplitude** (per event, per column)
  - **Duration** (per event, per column)
- **Independent noise** generation for each column
- **Complex output** including:
  - Time series matrix
  - Amplitude and duration matrices
  - HRF and noise parameter information

## Summary and Comparison

The four simulation functions in `fmrireg` serve different purposes and
offer increasing levels of complexity:

1.  **`simulate_bold_signal`**: Generate clean BOLD signals for multiple
    conditions
2.  **`simulate_noise_vector`**: Create structured, fMRI-like noise with
    temporal structure
3.  **`simulate_simple_dataset`**: Combine signal and noise with a
    specific SNR
4.  **`simulate_fmri_matrix`**: Create multiple time series with
    trial-by-trial, column-specific parameter variation

Choose the appropriate function based on your simulation needs:

- For basic signal generation: use `simulate_bold_signal`
- For structured teaching noise: use `simulate_noise_vector`
- For a complete dataset with controlled SNR: use
  `simulate_simple_dataset`
- For simulating multiple voxels/regions with shared timing but response
  variation: use `simulate_fmri_matrix`

These functions provide controlled inputs for implementation checks and
teaching. Their explicit stochastic and signal assumptions should be
inspected before using them as a benchmark for a scientific method.

## Next

- [`vignette("a_09_linear_model", package = "fmrireg")`](https://bbuchsbaum.github.io/fmrireg/articles/a_09_linear_model.md)
  — Fit and interpret an fMRI linear model
