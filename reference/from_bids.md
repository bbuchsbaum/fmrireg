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
  sessions = NULL,
  res = NULL,
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
  enhancement). The functional brain masks paired with each run are
  listed in each binding's `meta$runs$mask`.

- desc:

  BIDS `desc-` label of the preprocessed scans (default `"preproc"`).

- subjects:

  Optional character vector of subject labels (bare ids, e.g. `"01"`);
  default is all participants.

- sessions:

  Optional character vector of session labels (bare, e.g. `"1"`);
  default is every session.

- res:

  Optional BIDS `res-` label, for derivatives written at several
  resolutions; `NULL` accepts a single resolution per run.

- ...:

  Forwarded to
  [`bidser::read_confounds()`](https://bbuchsbaum.github.io/bidser/reference/read_confounds.html)
  (e.g. `npcs`, `clean`, `na_action`).

## Value

An `fmri_manifest`: a list of bindings, one per subject. Each binding's
`meta$runs` is a data frame with one row per run (`session`, `run`,
`scan`, `events`, `mask`).

## Details

Files are paired by their BIDS entities, parsed from the file names: a
preprocessed series, its events file and its confounds file belong to
the same run only when `sub`, `ses`, `task` and `run` are equal, and
labels are compared exactly (`subjects = "1"` never selects `sub-10`;
`task = "faces"` never selects `task-facesloc`). Runs are ordered by
session, then run number, and numbered `1, 2, ...` in that order; each
session's runs are blocks of their own. Two series for one run (another
resolution, acquisition or echo) are an error unless `res` picks one of
them.

The design formula's variables must match the columns of the BIDS
`events.tsv` (e.g. `trial_type`); a `run` column holding the block
number is added for the block structure. Confound columns are selected
by `confounds` (a character vector or
[`bidser::confound_set()`](https://bbuchsbaum.github.io/bidser/reference/confound_set.html)
result) passed straight to
[`bidser::read_confounds()`](https://bbuchsbaum.github.io/bidser/reference/read_confounds.html).

The TR is the `RepetitionTime` of the series' own JSON sidecar, else of
the raw series' sidecar (or an inherited `task-<label>_bold.json`), else
the NIfTI header's. Sidecars of other files under `derivatives/` (for
example maps written by
[`write_results()`](https://bbuchsbaum.github.io/fmrireg/reference/write_results.md))
are never read.

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
