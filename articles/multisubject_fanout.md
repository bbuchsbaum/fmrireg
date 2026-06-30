# Running one model across many subjects

``` r

library(fmrireg)
```

## The problem

A common study runs the **same model** across many subjects. The design
is constant – the same conditions, HRF, baseline, contrasts, and fitting
options – but the data changes per subject: different files, different
event order, different run lengths. Looping
[`fmri_lm()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_lm.md)
by hand is tedious and does not lend itself to parallel or cluster
execution.

`fmrireg` separates the part that is **invariant** (an `fmri_template`)
from the part that **varies per subject** (a *binding*: that subject’s
data). Combining the two produces a small, serializable `fmri_job` – a
recipe that holds no voxel data and can be shipped to a worker (a
`future`, or a node in an array job) that reconstructs the dataset,
fits, and reduces.

The four pieces:

| Piece | What it is |
|----|----|
| [`fmri_template()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_template.md) | the invariant model: formula, baseline spec, contrasts, fit control, reducer |
| a *binding* | one subject’s data: scans, run lengths, TR, events, confounds |
| `fmri_job` | template + one binding – a serializable recipe (built by [`instantiate()`](https://bbuchsbaum.github.io/fmrireg/reference/instantiate.md)) |
| a *reducer* | what each fitted subject becomes (a tidy table, or files on disk) |

This closes a loop: BIDS in (via `bidser`) → per-subject fit →
BIDS-keyed maps out (`write_results`) → group level (`collect_results` →
`fmri_meta`).

## A minimal, runnable example

Define the template once, describe two subjects as in-memory bindings,
then fan the model out. (Real studies use file paths; in-memory
`matrix_dataset` bindings keep this example fast.)

``` r

# one subject's data as a binding
make_subject <- function(id, seed) {
  set.seed(seed)
  runs <- c(40L, 40L); TR <- 2
  Y <- matrix(rnorm(sum(runs) * 5), sum(runs), 5)        # timepoints x voxels
  ev <- do.call(rbind, lapply(1:2, function(r)
    data.frame(onset = seq(4, 60, by = 8),
               condition = factor(rep(c("A", "B"), length.out = 8)),
               run = r)))
  list(id = id, scans = Y, TR = TR, run_length = runs, events = ev)
}
subjects <- list(make_subject("sub-01", 1), make_subject("sub-02", 2))

# the invariant model, defined once
tmpl <- fmri_template(
  onset ~ hrf(condition), ~ run,
  baseline = baseline_spec(degree = 3),
  reducer  = reduce_betas()        # each subject -> a tidy beta table
)

# bind data -> jobs -> run
jobs <- instantiate(tmpl, subjects)
res  <- run_jobs(jobs)
res
#> <fmri_batch_result> 2 job(s): 2 ok, 0 failed

values <- batch_values(res)        # named by job id
head(values[["sub-01"]])
#>   job_id                  term voxel    estimate        se        stat
#> 1 sub-01 condition_condition.A     1 -0.01645972 0.2494796 -0.06597621
#> 2 sub-01 condition_condition.B     1 -0.16257429 0.2485769 -0.65402017
#> 3 sub-01 condition_condition.A     2  0.02754297 0.2478427  0.11113087
#> 4 sub-01 condition_condition.B     2 -0.15066662 0.2469459 -0.61011997
#> 5 sub-01 condition_condition.A     3 -0.72182835 0.2720789 -2.65301068
#> 6 sub-01 condition_condition.B     3 -0.47638960 0.2710944 -1.75728281
```

[`run_jobs()`](https://bbuchsbaum.github.io/fmrireg/reference/run_jobs.md)
isolates per-job failures: a subject that errors is recorded in the
result rather than aborting the batch (see
[`batch_errors()`](https://bbuchsbaum.github.io/fmrireg/reference/batch_errors.md)).

## Reducers: what each subject becomes

A reducer runs on the worker, right after the fit, so only the reduced
output crosses the worker→driver (or worker→disk) boundary – not a whole
fitted model.

| Reducer | Output |
|----|----|
| [`reduce_identity()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_identity.md) | the entire `fmri_lm` object (largest) |
| [`reduce_betas()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_betas.md) | tidy data frame: `job_id, term, voxel, estimate, se, stat` |
| [`reduce_contrasts()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_contrasts.md) | tidy data frame of fitted contrasts |
| [`reduce_write_results()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_write_results.md) | writes BIDS-keyed maps to disk, returns the paths |

You can also pass your own `function(fit, job)`; keep it a
top-level/package function so it serializes to a worker.

## Discovering subjects from BIDS

For BIDS-formatted data,
[`from_bids()`](https://bbuchsbaum.github.io/fmrireg/reference/from_bids.md)
(which uses the `bidser` package) builds the per-subject bindings –
scans, events, confounds, TR, and run lengths – so you do not assemble
them by hand.

``` r

proj <- bidser::bids_project("study/", fmriprep = TRUE)
mani <- from_bids(
  proj, task = "stroop", space = "MNI152NLin2009cAsym",
  confounds = bidser::confound_set("motion6"),
  mask = "study/derivatives/.../space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz"
)
jobs <- instantiate(tmpl, mani)
```

The design formula’s variables must match the columns of the BIDS
`events.tsv` (e.g. `trial_type`); a `run` column is added for the block
structure.

## Preflight before you fan out

[`preflight()`](https://bbuchsbaum.github.io/fmrireg/reference/preflight.md)
validates jobs on the driver – design columns present, TR and run
lengths consistent, confound dimensions correct – so problems surface in
seconds rather than on a compute node after the queue drains.

``` r

preflight(jobs)
```

## Running: locally, in parallel, or on a cluster

Sequential is the default. For parallelism, set a `future` plan;
[`run_jobs()`](https://bbuchsbaum.github.io/fmrireg/reference/run_jobs.md)
dispatches through it – including `future.batchtools` cluster plans –
with **no scheduler-specific code** in `fmrireg`.

``` r

library(future)
plan(multisession, workers = 4)        # or future.batchtools::batchtools_slurm
res <- run_jobs(jobs, parallel = TRUE)
```

For an array scheduler, export the jobs to disk and drive them with
whatever system you have:

``` r

export_jobs(jobs, "study/jobs")
# study/jobs/manifest.rds + a backend-agnostic run_one.R
# then, per array task:
#   SLURM: Rscript study/jobs/run_one.R $SLURM_ARRAY_TASK_ID
#   local: for i in $(seq 1 N); do Rscript study/jobs/run_one.R $i; done
```

`run_one.R` reads the manifest, reconstructs job `i`, fits, and writes
its reduced output. The job is a portable recipe, so the worker need
only have `fmrireg` installed.

## Closing the loop: the group level

With
[`reduce_write_results()`](https://bbuchsbaum.github.io/fmrireg/reference/reduce_write_results.md),
each worker writes BIDS-keyed statistical maps.
[`collect_results()`](https://bbuchsbaum.github.io/fmrireg/reference/collect_results.md)
gathers them back into a `group_data` object that
[`fmri_meta()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_meta.md)
meta-analyses. (Writing both `beta` and `se` requires a contrast in the
model, so the group step has per-subject variance.)

``` r

con <- contrast_set(pair_contrast(~ trial_type == "incongruent",
                                  ~ trial_type == "congruent",
                                  name = "incong_gt_cong"))
tmpl <- fmri_template(
  onset ~ hrf(trial_type, contrasts = con), ~ run,
  baseline = baseline_spec(degree = 3, confounds = bidser::confound_set("motion6")),
  reducer  = reduce_write_results(format = "nifti", stats = c("beta", "se"),
                                  path = "study/glm")
)
run_jobs(instantiate(tmpl, mani))

gd <- collect_results("study/glm", space = "MNI152NLin2009cAsym")
fm <- fmri_meta(gd, ~ 1, method = "fe")
```

## Summary

- Define the model once with
  [`fmri_template()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_template.md);
  vary only the data.
- [`instantiate()`](https://bbuchsbaum.github.io/fmrireg/reference/instantiate.md)
  produces serializable `fmri_job` recipes;
  [`from_bids()`](https://bbuchsbaum.github.io/fmrireg/reference/from_bids.md)
  populates them from a BIDS dataset.
- [`preflight()`](https://bbuchsbaum.github.io/fmrireg/reference/preflight.md)
  catches problems before fan-out.
- [`run_jobs()`](https://bbuchsbaum.github.io/fmrireg/reference/run_jobs.md)
  runs sequentially, in parallel via `future`, or through a custom
  backend;
  [`export_jobs()`](https://bbuchsbaum.github.io/fmrireg/reference/export_jobs.md)
  emits a scheduler-agnostic runner for array jobs.
- Reducers keep per-subject output compact;
  [`collect_results()`](https://bbuchsbaum.github.io/fmrireg/reference/collect_results.md)
  →
  [`fmri_meta()`](https://bbuchsbaum.github.io/fmrireg/reference/fmri_meta.md)
  close the loop at the group level.
