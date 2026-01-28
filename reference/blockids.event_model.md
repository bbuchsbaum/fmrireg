# Block IDs for event_model

Return the run/block IDs associated with an event_model's sampling
frame.

## Usage

``` r
# S3 method for class 'event_model'
blockids(x, ...)
```

## Arguments

- x:

  An event_model object

- ...:

  Additional arguments passed through

## Value

Integer vector of block IDs

## Examples

``` r
ev <- fmrireg:::.demo_event_model()
blockids(ev)
#> [1] 1 1 1 1 2 2 2 2
```
