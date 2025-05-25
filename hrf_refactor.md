# HRF Module – Surgical Refactor Plan

> **Overall mission** — collapse the duplicate HRF code‑paths into a single, composable pipeline built from tiny orthogonal verbs.  Every downstream consumer (evaluate(), Toeplitz, regressors, contrasts, …) continues to work unchanged.

---

## 0  Guiding idea

All historic helpers (`gen_hrf*()`, `HRF_*` constructors, _etc._) repeat the same three stages:

| Stage | Today | Refactor principle |
|-------|-------|--------------------|
| **Construction** | wrap a raw function and attach attributes | **one** tiny constructor: `as_hrf()` |
| **Transform** | shift, block‑convolve, normalise, bind basis sets via ad‑hoc helpers | orthogonal decorators (`lag_hrf()`, `block_hrf()`, `normalise_hrf()`, `bind_basis()`) |
| **Execution** | several bespoke evaluators; `convolve_block()` duplicates logic | **one** generic: `evaluate.HRF()` |

With exactly **one** generic for each stage, the forest of helpers collapses into a 12‑line pipeline.

---

## 1  One tiny constructor

```r
#' Turn any f(t) into an HRF object
as_hrf <- function(f, name = deparse(substitute(f)), nbasis = 1L, span = 24,
                   params = list()) {
  structure(
    f,
    class        = c("HRF", "function"),
    name         = name,
    nbasis       = as.integer(nbasis),
    span         = span,
    param_names  = names(params),
    params       = params
  )
}
```

All canonical HRFs become one‑liners:

```r
HRF_GAMMA    <- as_hrf(hrf_gamma,    "gamma",    params = c("shape", "rate"))
HRF_SPMG1    <- as_hrf(hrf_spmg1,    "spmg1")
HRF_GAUSSIAN <- as_hrf(hrf_gaussian, "gaussian", params = c("mean",  "sd"))
```

Need the SPMG tuple with derivatives?  Build it with `bind_basis()` (see § 2).

---

## 2  Four orthogonal decorators

```r
# shift in time ---------------------------------------------------------------
lag_hrf <- function(hrf, lag) {
  force(hrf)
  as_hrf(function(t) hrf(t - lag),
         name  = paste0(attr(hrf,"name"), "_lag", lag),
         nbasis = attr(hrf,"nbasis"),
         span   = attr(hrf,"span") + lag)
}

# block‑convolution (partial boxcar) -----------------------------------------
block_hrf <- function(hrf, width, precision = .1, half_life = Inf,
                      summate   = TRUE, normalise = FALSE) {
  force(hrf)
  fun <- function(t) {
    offs <- seq(0, width, by = precision)
    hmat <- vapply(offs, \(o) hrf(t - o) * exp(-o / half_life), numeric(length(t)))
    out  <- if (summate) rowSums(hmat)
            else         hmat[cbind(seq_along(t), apply(hmat,1,which.max))]
    if (normalise) out / max(abs(out)) else out
  }
  as_hrf(fun,
         name  = paste0(attr(hrf,"name"), "_block", width),
         nbasis = attr(hrf,"nbasis"),
         span   = attr(hrf,"span") + width)
}

# bind several HRFs into a basis set -----------------------------------------
bind_basis <- function(...) {
  xs <- list(...)
  as_hrf(
    function(t) do.call(cbind, lapply(xs, function(f) f(t))),
    name   = paste(sapply(xs, attr, "name"), collapse = "+"),
    nbasis = sum(vapply(xs, attr, 0L, "nbasis")),
    span   = max(vapply(xs, attr, 0,  "span"))
  )
}

# rescale (peak = 1) ----------------------------------------------------------
normalise_hrf <- function(hrf) {
  as_hrf(function(t) { y <- hrf(t); y / max(abs(y)) },
         name  = paste0(attr(hrf,"name"), "_norm"),
         nbasis = attr(hrf,"nbasis"),
         span   = attr(hrf,"span"))
}
```

**Example usage** – replaces five old helpers in one pipeline:

```r
HRF_SPMG1 |>
  lag_hrf(3)        |>  # temporal shift
  block_hrf(2)      |>  # boxcar of width 2 s
  normalise_hrf()      # peak = 1
```

---

## 3  Single evaluate generic

```r
evaluate.HRF <- function(x, grid, amplitude = 1, duration = 0,
                         precision = .2, summate = TRUE, normalise = FALSE, ...) {
  base <- function(g) amplitude * x(g)

  out <- if (duration < precision) {
    base(grid)
  } else {
    offs <- seq(0, duration, by = precision)
    hmat <- vapply(offs, \(o) base(grid - o), numeric(length(grid)))
    if (summate) rowSums(hmat) else hmat[cbind(seq_along(grid), apply(hmat,1,which.max))]
  }
  if (normalise) out / max(abs(out)) else out
}
```

`convolve_block()` / `gen_hrf_blocked()` disappear—every caller now hits this method.

---

## 4  Optional registry helper

```r
.HRF_LIBRARY <- list(
  spmg1    = HRF_SPMG1,
  gamma    = HRF_GAMMA,
  gaussian = HRF_GAUSSIAN
)

get_hrf <- function(name) .HRF_LIBRARY[[match.arg(name, names(.HRF_LIBRARY))]]
```

Old `getHRF()` shrinks to four lines.

---

## 5  Why `convolve_block()` can be deleted

| Old helper                                | New pipeline                                  |
|-------------------------------------------|-----------------------------------------------|
| `convolve_block(t, hrf_spmg1, width = 6)` | `block_hrf(HRF_SPMG1, width = 6)(t)` |

**Key facts**
* _Exact_ maths reproduced inside the decorator.
* Returned object still **inherits from "HRF"**, so any downstream call keeps working.
* Duplicated evaluation logic removed → single source of truth.

Feature parity table

| Capability                         | Old | New |
|------------------------------------|-----|-----|
| variable *width*                   | ✅ | ✅ |
| custom *precision*                 | ✅ | ✅ |
| `half_life` attenuation            | ✅ | ✅ |
| `summate = FALSE` (max‑pick)       | ✅ | ✅ |
| `normalise = TRUE`                 | ✅ | ✅ |
| multi‑basis compatible             | ✅ | ✅ |

---

## 6  LOC impact

| Area                                  | Before | After |
|---------------------------------------|--------|-------|
| `gen_hrf*`, `hrf_*` wrappers          | ~310   | 0 |
| individual `HRF_*` constructors       | ~120   | 15 |
| block helpers & duplicate evaluators  | ~140   | 35 |
| **Total saved**                       | **> 500** |

---

## 7  Migration & deprecation plan

1. **Introduce decorators** – add code + tests.
2. **Gen HRF** – replace `width` branch with `hrf <- block_hrf(hrf,width, …)`.
3. **Soft‑deprecate** `convolve_block()` / `gen_hrf_blocked()` using `.Deprecated("block_hrf")`.
4. **Update vignettes / docs** to show verb chaining.
5. **Remove legacy helpers** in next major version after grace period.

_No breaking change for user code that never called `convolve_block()` directly._

---

### TL;DR

> *We swap 500+ lines for ~120 by introducing **one** constructor, **four** decorators and a **single** `evaluate.HRF()`.  Users compose HRFs with pipe‑friendly verbs instead of memorising helper names; all S3 APIs and downstream analyses keep working unchanged.*

Operation Maximal  Elegance & Efficiency — HRF subsystem TODO

Goal: collapse the current 700  +  LOC forest of ad‑hoc helpers into a single mini‑DSL
as_hrf() |> lag_hrf() |> block_hrf() |> normalise_hrf() |> bind_basis()
while preserving every public behaviour, S3 signature, vignette example and unit test.

⸻

0 ▪ Ground rules

Principle Concretely means
No hard breaks All exported objects keep the same names & classes until the next major version.
One new thing at a time Constructor → decorators → evaluate → removal stubs, each merged only after tests pass.
Deprecate, don't delete Obsolete helpers stay as shims that .Deprecated()‑ping and call the new pipeline.
100 % test parity first Only after the old and new code paths produce byte‑identical results do we rip out legacy.



⸻

1 ▪ Scaffolding phase (➜ no user‑facing change)

# Task Files / functions Notes
1.1 Add the canonical constructor as_hrf() new‑file hrf-core.R Include attributes name, nbasis, span, param_names, params.
1.2 Port evaluate logic → re‑implemented evaluate.HRF() replace body in hrf-evaluate.R Must accept all old args (duration, summate, …) and dispatch to new pipeline.
1.3 Implement four decoratorslag_hrf(), block_hrf(), normalise_hrf(), bind_basis() same file Keep argument names identical to old helpers (width, precision, …).
1.4 Internal helpers• .assert_hrf()• .inherit_attrs() same file Utility for decorators to preserve metadata.
1.5 Add hidden registry .HRF_LIBRARY, get_hrf() replace old getHRF() guts later.



⸻

2 ▪ Parallellising phase (add new pipe, keep old)

# Task Files / functions Notes
2.1 Rewrite canonical objects (HRF_SPMG1, HRF_GAMMA, …) as one‑liners that call as_hrf() hrf-canon.R Unit tests must confirm they are still "HRF" "function" and numerical output identical.
2.2 Re‑implement gen_hrf_*, hrf_* wrappers as calls to the decorators hrf-legacy.R Do not remove exports yet; each wrapper becomes:function(...) { .Deprecated("lag_hrf/block_hrf"); lag_hrf(...); … }.
2.3 Re‑route getHRF() to the registry & wrappers ditto Old argument list untouched.
2.4 Update internal callers (regressor, hrfspec, Toeplitz, etc.) to use the new pipeline via old names → zero downstream diff. whole code‑base search Keep both versions side‑by‑side until Step 4.



⸻

3 ▪ Parity tests & benchmarks

# Task Coverage
3.1 Golden‑file tests — for each exported HRF object and each old helper call pattern, store a hash of evaluate(*, grid = seq(0,32,.1)).
3.2 Decorator combinatorics — random combos of lag, block, normalise, basis binding vs. old gen_hrf() permutations.
3.3 Performance bench — micro‑benchmark new vs. old on 10 k grids; assert ≤ 5 % slower or flag for optimisation.



⸻

4 ▪ Deprecation sweep (≋ minor version bump)

# Action
4.1 Mark gen_hrf*, hrf_blocked, convolve_block, gen_hrf_lagged, gen_hrf_set, HRF_* constructors (duplicated ones) with .Deprecated() notes pointing to the pipeline verbs.
4.2 Export a NEWS.md entry + vignette "Composing HRFs the new way".
4.3 Increment package minor version + set Deprecation: next‑x.y header for each function.



⸻

5 ▪ Removal sweep (≋ next major version)

# Remove Keep / move
5.1 All shim functions after grace period —
5.2 Redundant evaluate.hrfspec, convolve_block, duplicated make_hrf logic replaced by the core evaluate & decorators
5.3 LOC audit — expected ≥ 500 lines gone.



⸻

6 ▪ Documentation & examples
 • Update inline examples to use pipe‑style composition.
 • Add cheatsheet table: Old helper → New pipeline expression.
 • Enhance print.HRF to mention chain history (e.g. "spmg1 ▹ lag 2 ▹ block 6 ▹ norm").

⸻

7 ▪ Risk & compatibility checklist
 • S3 classes remain ("HRF","function") ⇒ all downstream generics untouched.
 • nbasis(), span, scale_factor attributes preserved.
 • Formula interface hrf() still works because it calls getHRF() which now returns new objects.
 • CRAN reverse‑dependency check on packages that import these helpers.

⸻

8 ▪ Timeline suggestion

Week Milestone
1 Steps 1.1 – 1.4 merged, unit tests green.
2 Steps 2.x live behind feature flag (options(hrf.refactor=TRUE)).
3 Parity tests pass → enable by default.
4 Deprecation release to CRAN (v Z.Y+1).
+6 mo Shim removal, major version bump.



⸻

Deliverables snapshot
 • R/hrf-core.R (≤ 250 LOC) — constructor, decorators, evaluate.
 • tests/testthat/test-hrf-refactor.R — gold‑file & speed tests.
 • vignettes/hrf-pipeline.Rmd — new user docs.
 • NEWS.md / lifecyle badges — deprecation guidance.

⸻

Next action:
Spin up branch refactor/hrf-pipeline, add as_hrf() & the four decorators, push for CI.


---
## Folding the "odds‑and‑ends" into the one‑page HRF core

Below I show how each of the three helpers melts away once we have the canonical
constructor as_hrf() and the composable decorator bind_basis().

⸻

1 Empirical HRF  →  3‑line wrapper

# old ------------------------------------------------------
gen_empirical_hrf <- function(t, y, name = "empirical_hrf") {
  f <- approxfun(t, y, yright = 0, yleft = 0)
  HRF(f, name = name, nbasis = 1)
}

# new ------------------------------------------------------
empirical_hrf <- function(t, y, name = "empirical_hrf") {
  as_hrf(approxfun(t, y, yright = 0, yleft = 0),
         name = name, nbasis = 1, span = max(t))
}

Exactly the same public behaviour, just written with the single constructor.
If we want to keep the old name for backward compatibility:

gen_empirical_hrf <- function(...) {
  .Deprecated("empirical_hrf()")
  empirical_hrf(...)
}



⸻

2 Basis‑set combiner  →  bind_basis()

# old ------------------------------------------------------
gen_hrf_set <- function(..., span = 32, name = "hrf_set") { … }

# new ------------------------------------------------------
hrf_set <- function(..., name = "hrf_set") {
  bind_basis(...) |>
    (\(h) { attr(h, "name") <- name; h })()
}

All the busywork (checking each element is a function, cbind‑eval, computing
nbasis and span) already lives inside bind_basis().
Again, deprecate the alias if you wish:

gen_hrf_set <- function(...) {
  .Deprecated("hrf_set()")
  hrf_set(...)
}



⸻

3 Parameter‑grid library  →  one tidy‑map call

# old ------------------------------------------------------
gen_hrf_library <- function(fun, pgrid, ...) { … }

# new ------------------------------------------------------
hrf_library <- function(fun, pgrid, ...) {
  purrr::pmap(pgrid, \(...) fun |> purrr::partial(...) |> as_hrf(...)) |>
    (\(fs) do.call(bind_basis, fs))()
}

Differences / improvements

Aspect	Old code	New one‑liner
Handles arbitrary parameter names?	yes	yes (pmap)
Carries each HRF's own attributes?	partly	yes (as_hrf() inside the loop)
Combines into basis set?	ad‑hoc gen_hrf_set()	bind_basis()
Lines of code	~25	5

Deprecation wrapper:

gen_hrf_library <- function(...) {
  .Deprecated("hrf_library()")
  hrf_library(...)
}



⸻

Minimal impact checklist ✔︎

Down‑stream function	Still works?	Why
evaluate()	✔︎	All objects inherit function + "HRF", attributes unchanged.
nbasis(), print.*	✔︎	Attributes set by as_hrf() / bind_basis().
Formula DSL hrf(cond, basis = …)	✔︎	make_hrf() calls as_hrf()‑based constructors.
Legacy code	✔︎ w/ message	Compatibility shims using .Deprecated().



⸻

Take‑away

Those "odds and ends" collapse to ≤ 20 tidy lines once the single‑source
primitives (as_hrf(), bind_basis()) exist.
This keeps the public API intact and enforces the one constructor / one
decorator rule we set for maximal elegance and maintainability.