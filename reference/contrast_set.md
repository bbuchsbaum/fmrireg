# Create a contrast set

Create a contrast set

## Usage

``` r
contrast_set(...)
```

## Arguments

- ...:

  contrast specifications

## Value

A list of class "contrast_set" containing the specified contrasts

## Examples

``` r
cs <- contrast_set(
  fmridesign::pair_contrast(~condition == "A", ~condition == "B", name = "A_vs_B")
)
```
