# Build a manifest from a BIDS project

Uses bidser to discover preprocessed BOLD scans, events, confounds, TR,
and run lengths for each subject of a BIDS dataset, returning an
`fmri_manifest` (a list of per-subject bindings) suitable for
[`instantiate()`](https://bbuchsbaum.github.io/fmrireg/reference/instantiate.md).

## Usage

``` r
from_bids(
  proj,
  task,
  space = NULL,
  confounds = NULL,
  mask = NULL,
  desc = "preproc",
  subjects = NULL,
  ...
)
```

## Arguments

- proj:

  A
  [`bidser::bids_project`](https://bbuchsbaum.github.io/bidser/reference/bids_project.html)
  (opened with `fmriprep = TRUE` for derivative discovery).

- task:

  BIDS task label.

- space:

  BIDS space for preprocessed scans (e.g. `"MNI152NLin2009cAsym"`);
  `NULL` matches any.

- confounds:

  Confound selection passed to
  [`bidser::read_confounds()`](https://bbuchsbaum.github.io/bidser/reference/read_confounds.html)
  as `cvars` (character vector or a `confound_set`); `NULL` reads no
  confounds.

- mask:

  A brain-mask path (length-1 character) applied to every subject's
  file-backed dataset. Required for file-backed BIDS data, which needs a
  mask (auto-resolution of the fmriprep brain mask is a planned
  enhancement).

- desc:

  BIDS `desc-` label of the preprocessed scans (default `"preproc"`).

- subjects:

  Optional character vector of subject labels (bare ids, e.g. `"01"`);
  default is all participants.

- ...:

  Forwarded to
  [`bidser::read_confounds()`](https://bbuchsbaum.github.io/bidser/reference/read_confounds.html)
  (e.g. `npcs`, `clean`, `na_action`).

## Value

An `fmri_manifest`: a list of bindings, one per subject.

## Details

The design formula's variables must match the columns of the BIDS
`events.tsv` (e.g. `trial_type`); a `run` column is added for the block
structure. Confound columns are selected by `confounds` (a character
vector or
[`bidser::confound_set()`](https://bbuchsbaum.github.io/bidser/reference/confound_set.html)
result) passed straight to
[`bidser::read_confounds()`](https://bbuchsbaum.github.io/bidser/reference/read_confounds.html).

## See also

[`as_manifest()`](https://bbuchsbaum.github.io/fmrireg/reference/as_manifest.md),
[`instantiate()`](https://bbuchsbaum.github.io/fmrireg/reference/instantiate.md),
[`preflight()`](https://bbuchsbaum.github.io/fmrireg/reference/preflight.md)

## Examples

``` r
if (FALSE) { # \dontrun{
proj <- bidser::bids_project("study/", fmriprep = TRUE)
mani <- from_bids(proj, task = "stroop", space = "MNI152NLin2009cAsym",
                  confounds = bidser::confound_set("motion6"),
                  mask = "study/derivatives/.../space-MNI..._desc-brain_mask.nii.gz")
jobs <- instantiate(template, mani)
} # }
```
