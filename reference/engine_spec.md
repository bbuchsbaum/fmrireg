# Inspect a Registered Engine Specification

Returns a read-only description of a registered engine, including its
normalized capabilities, source (`"builtin"` vs `"plugin"`), aliases,
and dispatch strategy. This is intended for extension authors and
diagnostic tooling; it does not expose the underlying fit/preflight
functions.

## Usage

``` r
engine_spec(name)
```

## Arguments

- name:

  Engine name, such as `"rrr_gls"` or `"latent_sketch"`.

## Value

A list of class `fmrireg_engine_spec`, or `NULL` if the engine is not
registered.
