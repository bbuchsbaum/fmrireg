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
