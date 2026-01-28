# Register an HRF/basis constructor

Register an HRF/basis constructor

## Usage

``` r
register_basis(name, constructor)
```

## Arguments

- name:

  Character scalar used in formulas (e.g. `basis = "friman2"`).

- constructor:

  Function returning an object understood by `fmrihrf` (typically an
  `HRF`). Additional arguments from the original
  [`hrf()`](https://bbuchsbaum.github.io/fmridesign/reference/hrf.html)
  call are forwarded when the basis is constructed.

## Value

Invisibly, `TRUE`.
