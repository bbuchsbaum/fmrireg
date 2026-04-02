# Design Plot for fMRI Model

Generates an interactive Shiny app that plots the design matrix for a
given fMRI model. The design matrix is first converted into a
long-format tibble and then plotted over time, faceted by block. Several
customization options allow you to adjust the title, axis labels, line
thickness, color palette, and more.

## Usage

``` r
design_plot(
  fmrimod,
  term_name = NULL,
  longnames = FALSE,
  plot_title = NULL,
  x_label = "Time (s)",
  y_label = "Amplitude",
  line_size = 1,
  color_palette = "viridis",
  facet_ncol = 2,
  theme_custom = ggplot2::theme_minimal(base_size = 15) + ggplot2::theme(panel.spacing =
    ggplot2::unit(1, "lines")),
  legend_threshold = 30,
  ...
)
```

## Arguments

- fmrimod:

  An `fmri_model` object.

- term_name:

  Optional: Name of the term to plot. If `NULL` (the default), the first
  term is used.

- longnames:

  Logical; if TRUE, use long condition names in the legend. Default is
  FALSE.

- plot_title:

  Optional plot title. If `NULL`, a default title is generated.

- x_label:

  Label for the x-axis. Default is "Time".

- y_label:

  Label for the y-axis. Default is "Value".

- line_size:

  Numeric; line thickness for the plot. Default is 1.

- color_palette:

  Character; name of a ColorBrewer palette to use (e.g., "Set1").
  Default is "Set1".

- facet_ncol:

  Number of columns for facet_wrap. Default is 1.

- theme_custom:

  A ggplot2 theme to apply. Default is `theme_bw(base_size = 14)`.

- legend_threshold:

  Numeric; if the number of unique conditions exceeds this value, the
  legend is hidden. Default is 25.

- ...:

  Additional arguments passed to ggplot2::geom_line().

## Value

A Shiny app that displays the design plot.

## Examples

``` r
if (interactive()) {
  ## --- Construct a sampling frame ---
  sframe <- fmrihrf::sampling_frame(blocklens = c(100, 100), TR = 2, precision = 0.5)

  ## --- Create a dummy event table ---
  set.seed(123)
  event_table <- data.frame(
    onset = seq(10, 190, length.out = 20),
    x = rnorm(20),
    y = rnorm(20),
    run = rep(1:2, each = 10)
  )

  ## --- Construct a baseline model ---
  base_mod <- baseline_model(basis = "bs", degree = 3, sframe = sframe, intercept = "runwise")

  ## --- Construct an event model using a formula ---
  ev_mod <- event_model(x = onset ~ hrf(x) + hrf(y), data = event_table,
                        block = ~ run, sampling_frame = sframe,
                        drop_empty = TRUE, durations = rep(0, nrow(event_table)))

  ## --- Combine into an fMRI model ---
  fmri_mod <- fmri_model(ev_mod, base_mod)

  ## --- Launch the interactive design plot ---
  design_plot(fmrimod = fmri_mod,
              term_name = NULL,
              longnames = TRUE,
              plot_title = "fMRI Design Matrix",
              x_label = "Time (s)",
              y_label = "Signal",
              line_size = 1.5,
              color_palette = "Set2",
              facet_ncol = 1)
}
```
