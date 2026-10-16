# Low-rank / sketch controls for fast GLM

Control object to enable the optional sketched GLM engine.

## Usage

``` r
lowrank_control(
  parcels = NULL,
  landmarks = NULL,
  k_neighbors = 16L,
  time_sketch = list(method = "gaussian", m = NULL),
  ncomp = NULL,
  noise_pcs = 0L
)
```

## Arguments

- parcels:

  Optional parceling, e.g., a neuroim2::ClusteredNeuroVol or integer
  vector (length = number of voxels in mask).

- landmarks:

  Optional integer; number of landmark voxels for optional Nyström
  extension (NULL = off).

- k_neighbors:

  Integer; k for k-NN in Nyström extension.

- time_sketch:

  List controlling the temporal sketch, with elements

  `method`

  :   `"gaussian"` (default), `"countsketch"` or `"srht"` (subsampled
      randomized Hadamard transform). `"ihs"` (iterative Hessian sketch)
      is deprecated: it now computes the exact OLS fit and emits a
      once-per-session message.

  `m`

  :   Number of sketch rows; `NULL` uses `min(8 * p, T)` for a design
      with `p` columns and `T` scans. The sketch-and-solve methods
      require `m > p`; `"ihs"` ignores `m`.

  The former `"ihs"` controls `iters` and `tol` are accepted and
  ignored.

  All sketches are normalised (\\E\\Sr\\^2 = \\r\\^2\\), so residual
  variances, `rss` and covariances are reported on the data scale.
  `"gaussian"`, `"countsketch"` and `"srht"` are sketch-and-solve
  estimators: they fit the model to `m` sketched rows, and their
  standard errors are the conditional-on-sketch standard errors of that
  estimator (roughly `sqrt(T / m)` times the full-data OLS ones), with
  Satterthwaite residual degrees of freedom (about `m - p`). `"ihs"`
  returns exact OLS coefficients, covariance, residual variance and
  `T - rank` degrees of freedom. With many response columns and few
  design columns, one exact least-squares pass costs a single \\X'Z\\
  product, which no iterative sketch can undercut; `"ihs"` was never
  faster than exact OLS here.

  For speed, prefer `"countsketch"` or `"gaussian"`. With `p = 20`,
  `m = 8p` and `r = 5e4` response columns, a sketch-and-solve fit took
  0.04-0.06 s for `T` from 400 to 1000 scans, against 0.05-0.14 s for
  exact OLS (Apple M3 Max, Accelerate BLAS). `"srht"` took 0.07-0.09 s:
  its fast Walsh-Hadamard transform costs \\O(T \log T)\\ per response
  column, CountSketch touches each scan once (\\O(T)\\), and the
  Gaussian sketch's \\O(mT)\\ dense product runs at BLAS speed, so
  `"srht"` loses to both at these sizes and overtakes exact OLS only for
  longer series. Its advantage is an orthonormal-row sketch (\\SS' =
  (T/m)I\\ when `T` is a power of two). The SRHT kernel runs in parallel
  over response columns with `getOption("fmrireg.num_threads")` threads
  when the package is built with OpenMP.

  For sketch-and-solve fits, `result$rss` is the residual sum of squares
  of the `m` sketched rows. Its expectation is \\\sigma^2 \kappa\\, with
  \\\kappa = tr(PSS')\\ stored in `result$sketch$kappa`, so the residual
  variance is `rss / kappa`, not `rss / rdf`. For `"ihs"`, `kappa`
  equals the residual degrees of freedom. With `landmarks`, `rss` is
  interpolated from the landmark voxels, like the residual variance.

  Estimability is judged on the full (whitened) design with the same
  pivoted-QR rule that exact fits use. Aliased coefficients are reported
  as `NA`, and contrasts that load on them are `NA` with a warning. If
  the sketched design restricted to the estimable columns is numerically
  rank deficient, `m` is too small for the design and the fit stops with
  an error rather than regularising.

- ncomp:

  Optional integer; number of latent components within parcels (PCA).

- noise_pcs:

  Integer; optional GLMdenoise-style PCs from low-R2 parcels.

## Value

A list with class "lowrank_control".
