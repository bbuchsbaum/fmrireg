# Resolve a registered basis function by name

This function is used internally by the formula processing system to
lazily resolve basis names (e.g., "cca2", "cca3", "spmg1") into actual
basis objects. It's part of the plugin API that allows extension
packages to register custom basis functions.

## Usage

``` r
resolve_basis(name, ...)
```

## Arguments

- name:

  Character string naming a registered basis function

- ...:

  Additional arguments passed to the basis constructor

## Value

An HRF basis object

## Examples

``` r
if (FALSE) { # \dontrun{
# Resolve the cca3 basis with specific parameters
basis <- resolve_basis("cca3", span = 30, TR = 2)
} # }
```
