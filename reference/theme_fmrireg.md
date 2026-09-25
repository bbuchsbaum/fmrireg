# The fmrireg ggplot2 theme

A restrained theme used by all fmrireg plots: hairline major grid only,
left-aligned titles, and muted axis ink so the data carry the contrast.
Add it to any ggplot, or override it by adding your own theme
afterwards.

## Usage

``` r
theme_fmrireg(
  base_size = 11,
  base_family = "",
  grid = c("y", "x", "xy", "none")
)
```

## Arguments

- base_size:

  Base font size in points.

- base_family:

  Base font family.

- grid:

  Which major gridlines to keep: `"y"`, `"x"`, `"xy"`, or `"none"`.

## Value

A
[`ggplot2::theme()`](https://ggplot2.tidyverse.org/reference/theme.html)
object.

## Details

fmrireg plots use this theme by default. To use a different base theme
for every fmrireg plot, set `options(fmrireg.plot_theme = <theme>)`, for
example `options(fmrireg.plot_theme = ggplot2::theme_bw())`.

## Examples

``` r
library(ggplot2)
ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_fmrireg()
```
