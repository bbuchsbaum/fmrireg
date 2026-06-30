# Collect per-subject result maps into a group_data object

Globs the per-subject statistical maps written by
[`reduce_write_results()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_write_results.md)
under `dir`, pairs them by subject, and builds a
[`group_data()`](https://bbuchsbaum.github.io/fmrireg/reference/group_data.md)
object for
[`fmri_meta()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_meta.md).
When both beta and se maps are present the result carries variance
(suitable for inverse-variance meta-analysis); a beta-only result is
returned otherwise.

## Usage

``` r
collect_results(
  dir,
  space = NULL,
  beta_desc = "beta",
  se_desc = "se",
  format = c("nifti")
)
```

## Arguments

- dir:

  Directory (searched recursively) containing the written maps.

- space:

  Optional BIDS space label to filter on (e.g. `"MNI152NLin2009cAsym"`);
  `NULL` takes all.

- beta_desc, se_desc:

  The BIDS `desc-` labels of the effect and standard-error maps
  (defaults `"beta"` / `"se"`).

- format:

  Group-data backend format (currently `"nifti"`).

## Value

A `group_data` object.

## Note

Pairing by BIDS subject requires fmrigds with BIDS-aware subject keying.
Multi-file `"h5"` aggregation is not yet supported upstream.

## See also

[`reduce_write_results()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_write_results.md),
[`run_jobs()`](https://bbuchsbaum.github.io/fmrireg/reference/run_jobs.md),
[`fmri_meta()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_meta.md)

## Examples

``` r
if (FALSE) { # \dontrun{
tmpl <- fmri_template(onset ~ hrf(trial_type, contrasts = con), ~ run,
                      reducer = reduce_write_results(format = "nifti",
                                                     stats = c("beta", "se"),
                                                     path = "study/glm"))
run_jobs(instantiate(tmpl, manifest))
gd <- collect_results("study/glm", space = "MNI152NLin2009cAsym")
fm <- fmri_meta(gd, ~ 1, method = "fe")
} # }
```
