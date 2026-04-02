# Advanced Plugin Development

## Why use the plugin API?

The plugin API is for advanced users who want to extend `fmrireg`
without forking the package. Two extension points are available:

- [`register_basis()`](https://bbuchsbaum.github.io/fmrireg/reference/register_basis.md)
  adds new basis names that can be used inside `hrf(...)`.
- [`register_engine()`](https://bbuchsbaum.github.io/fmrireg/reference/register_engine.md)
  adds new fitting engines that still return standard `fmri_lm` objects.

This vignette shows how to register both kinds of extensions, how engine
capabilities interact with global
[`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md)
options, and how to inspect the engine-spec metadata that drives
dispatch.

## How do you build a minimal dataset?

We will work with a small matrix-backed dataset so every example is fast
and fully reproducible.

``` r
n_time <- 80L
n_voxels <- 24L

events <- data.frame(
  onsets = c(8, 20, 32, 44, 56),
  condition = factor(c("A", "B", "A", "B", "A")),
  run = 1L
)

Y <- matrix(rnorm(n_time * n_voxels), nrow = n_time, ncol = n_voxels)
dset <- fmridataset::matrix_dataset(
  Y,
  TR = 2,
  run_length = n_time,
  event_table = events
)

stopifnot(inherits(dset, "matrix_dataset"))
dim(fmridataset::get_data_matrix(dset))
#> [1] 80 24
```

## How do you register a custom basis?

[`register_basis()`](https://bbuchsbaum.github.io/fmrireg/reference/register_basis.md)
maps a string name to a constructor. That name can then be used inside
`hrf(...)` just like built-in bases.

``` r
basis_name <- "vignette_bspline_basis"

register_basis(
  basis_name,
  function(span = 18, ...) {
    fmrihrf::gen_hrf(fmrihrf::hrf_bspline, N = 5, span = span)
  }
)

model_with_basis <- create_fmri_model(
  formula = onsets ~ hrf(condition, basis = basis_name, span = 16),
  block = ~run,
  dataset = dset
)

stopifnot(inherits(model_with_basis, "fmri_model"))
head(colnames(design_matrix(model_with_basis)))
#> [1] "condition_condition.A_b01" "condition_condition.B_b01"
#> [3] "condition_condition.A_b02" "condition_condition.B_b02"
#> [5] "condition_condition.A_b03" "condition_condition.B_b03"
```

The important rule is that the constructor must return an object that
`fmrihrf` understands as a valid basis.

## How do you register a custom engine?

An engine receives `(model, dataset, args, cfg)` and should return an
`fmri_lm` object. In many cases the easiest path is:

1.  extract or transform the response matrix,
2.  call one of the helper constructors,
3.  return the resulting `fmri_lm` object.

The example below centers each voxel time series before fitting and then
delegates to
[`fit_glm_on_transformed_series()`](https://bbuchsbaum.github.io/fmrireg/reference/fit_glm_on_transformed_series.md).

``` r
engine_name <- "vignette_centered_engine"

register_engine(
  name = engine_name,
  capabilities = list(
    robust = FALSE,
    preprocessing = FALSE
  ),
  fit = function(model, dataset, args, cfg) {
    Y_raw <- as.matrix(fmridataset::get_data_matrix(dataset))
    Y_centered <- if (isTRUE(args$center)) {
      sweep(Y_raw, 2, colMeans(Y_raw), FUN = "-")
    } else {
      Y_raw
    }

    fit_glm_on_transformed_series(
      model,
      Y_centered,
      cfg = cfg,
      dataset = dataset,
      engine = engine_name,
      strategy = "engine"
    )
  }
)

fit_plugin <- fmri_lm(
  onsets ~ hrf(condition, basis = basis_name, span = 16),
  block = ~run,
  dataset = dset,
  engine = engine_name,
  engine_args = list(center = TRUE),
  robust = FALSE,
  robust_options = list(max_iter = 10L)
)

stopifnot(inherits(fit_plugin, "fmri_lm"))
stopifnot(identical(attr(fit_plugin, "engine"), engine_name))
stopifnot(identical(attr(fit_plugin, "strategy"), "engine"))

list(
  engine = attr(fit_plugin, "engine"),
  strategy = attr(fit_plugin, "strategy"),
  n_beta_rows = nrow(fit_plugin$result$betas)
)
#> $engine
#> [1] "vignette_centered_engine"
#> 
#> $strategy
#> [1] "engine"
#> 
#> $n_beta_rows
#> [1] 1
```

## Which helper should you use?

Use
[`fit_glm_on_transformed_series()`](https://bbuchsbaum.github.io/fmrireg/reference/fit_glm_on_transformed_series.md)
when your engine has already produced a time-by-voxel matrix and you
want standard OLS-style inference.

Use
[`fit_glm_with_config()`](https://bbuchsbaum.github.io/fmrireg/reference/fit_glm_with_config.md)
when your engine still wants the integrated AR/robust solver applied to
a transformed response matrix.

Use
[`fit_glm_from_suffstats()`](https://bbuchsbaum.github.io/fmrireg/reference/fit_glm_from_suffstats.md)
when your engine works from streamed cross-products and never
materializes the transformed series directly.

## How do capabilities affect the executed configuration?

The central dispatcher validates unsupported options before engine code
runs. It also keeps two configurations on the fitted object:

- `requested_config`: what the caller supplied to
  [`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md)
- `executed_config`: the engine-scoped configuration actually passed
  into the engine after unsupported sections were normalized away. In
  the current implementation, subordinate robust-tuning values such as
  `max_iter` may be reset even when the caller already requested
  `robust = FALSE`.

``` r
requested_cfg <- attr(fit_plugin, "requested_config")
executed_cfg <- attr(fit_plugin, "executed_config")

cfg_value <- function(cfg, section, field) {
  if (is.null(cfg) || is.null(cfg[[section]]) || is.null(cfg[[section]][[field]])) {
    return(NA)
  }
  cfg[[section]][[field]]
}

requested_max_iter <- cfg_value(requested_cfg, "robust", "max_iter")
executed_max_iter <- cfg_value(executed_cfg, "robust", "max_iter")

data.frame(
  config = c("requested", "executed"),
  robust_type = c(
    as.character(cfg_value(requested_cfg, "robust", "type")),
    as.character(cfg_value(executed_cfg, "robust", "type"))
  ),
  robust_max_iter = c(
    requested_max_iter,
    executed_max_iter
  ),
  normalized = c(
    FALSE,
    !isTRUE(all.equal(requested_max_iter, executed_max_iter))
  )
)
#>      config robust_type robust_max_iter normalized
#> 1 requested       FALSE              10      FALSE
#> 2  executed       FALSE               2       TRUE
```

If a caller enables an unsupported feature such as `robust = TRUE`, the
engine is rejected before `fit()` is called.

``` r
unsupported_message
#> [1] "vignette_centered_engine does not support robust fitting; set robust = FALSE"
```

## How can you inspect engine specs?

Built-in and plugin engines now share the same spec object. The
read-only accessors
[`engine_spec()`](https://bbuchsbaum.github.io/fmrireg/reference/engine_spec.md)
and
[`engine_specs()`](https://bbuchsbaum.github.io/fmrireg/reference/engine_specs.md)
expose that metadata without requiring access to internal registries. If
you are working against an older installed build where those helpers are
not yet available, you can still fall back to the registry names for a
lightweight diagnostic.

``` r
has_public_specs <- exists("engine_spec", envir = asNamespace("fmrireg"), inherits = FALSE)

if (has_public_specs) {
  plugin_spec <- get("engine_spec", envir = asNamespace("fmrireg"))(engine_name)
  builtin_spec <- get("engine_spec", envir = asNamespace("fmrireg"))("rrr_gls")

  print(plugin_spec)
  print(builtin_spec)
} else {
  registry <- get(".fmrireg_engine_registry", envir = asNamespace("fmrireg"))
  data.frame(
    name = c(engine_name, "rrr_gls"),
    registered = c(
      exists(engine_name, envir = registry, inherits = FALSE),
      exists("rrr_gls", envir = registry, inherits = FALSE)
    )
  )
}
#> <fmrireg_engine_spec>
#> name: vignette_centered_engine
#> source: plugin | strategy: engine
#> aliases: <none>
#> capabilities: robust=FALSE, preprocessing=FALSE, ar_voxelwise=TRUE, ar_by_cluster=TRUE
#> <fmrireg_engine_spec>
#> name: rrr_gls
#> source: builtin | strategy: engine
#> aliases: <none>
#> capabilities: robust=FALSE, preprocessing=FALSE, ar_voxelwise=FALSE, ar_by_cluster=FALSE
#> requires: event regressors
```

The full registry is also available as a list of specs.

``` r
if (exists("engine_specs", envir = asNamespace("fmrireg"), inherits = FALSE)) {
  spec_names <- names(get("engine_specs", envir = asNamespace("fmrireg"))())
} else {
  registry <- get(".fmrireg_engine_registry", envir = asNamespace("fmrireg"))
  spec_names <- sort(ls(registry, all.names = TRUE))
}

spec_names
#> [1] "latent_sketch"            "rrr_gls"                 
#> [3] "vignette_centered_engine"
```

## What should extension authors remember?

- Keep engine inputs small and explicit: `model`, `dataset`, `args`,
  `cfg`.
- Declare unsupported global features in `capabilities` so the
  dispatcher can reject them before your engine runs.
- Use `requested_config` versus `executed_config` when you need to
  explain what was requested versus what your engine actually used.
- Prefer the helper constructors over hand-building `fmri_lm` objects
  unless your engine truly needs custom result assembly.
