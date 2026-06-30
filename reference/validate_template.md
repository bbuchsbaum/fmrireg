# Validate an fmri_template

Re-checks the structural invariants of a template. Useful in preflight
before fanning a model out over many subjects.

## Usage

``` r
validate_template(x)
```

## Arguments

- x:

  An
  [fmri_template](https://bbuchsbaum.github.io/fmrireg/reference/fmri_template.md).

## Value

`TRUE` invisibly if valid; otherwise an error is raised.
