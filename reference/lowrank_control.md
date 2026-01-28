# Low-rank / sketch controls for fast GLM

Control object to enable the optional sketched GLM engine.

## Usage

``` r
lowrank_control(
  parcels = NULL,
  landmarks = NULL,
  k_neighbors = 16L,
  time_sketch = list(method = "gaussian", m = NULL, iters = 0L),
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

  List(method = "gaussian" \| "countsketch", m = NULL, iters = 0L).

- ncomp:

  Optional integer; number of latent components within parcels (PCA).

- noise_pcs:

  Integer; optional GLMdenoise-style PCs from low-R2 parcels.

## Value

A list with class "lowrank_control".
