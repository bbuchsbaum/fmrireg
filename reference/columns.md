# Extract Column Names or Identifiers

Extract column names or identifiers from an object. For parametric basis
objects, this returns tokens representing the type of variables
(categorical or continuous).

## Usage

``` r
columns(x, ...)

# S3 method for class 'event_model'
columns(x, ...)
```

## Arguments

- x:

  An object (typically a ParametricBasis)

- ...:

  Additional arguments passed to method-specific implementations.

## Value

A character vector of column identifiers

## Examples

``` r
dat <- data.frame(
  onsets = c(0, 4, 8),
  condition = factor(c("A", "B", "A")),
  run = 1
)
ev <- event_model(
  onsets ~ hrf(condition),
  data = dat,
  block = ~ run,
  sampling_frame = fmrihrf::sampling_frame(blocklens = 12, TR = 2)
)
columns(ev)
#> [1] "condition_condition.A" "condition_condition.B"
```
