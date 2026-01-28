# Apply Nyström / barycentric extension: B (p x V) from B_L (p x L)

Apply Nyström / barycentric extension: B (p x V) from B_L (p x L)

## Usage

``` r
extend_betas_landmarks(BL, W)
```

## Arguments

- BL:

  p x L dense matrix of betas on landmarks

- W:

  V x L sparse weight matrix (rows sum to 1)

## Value

p x V dense matrix of extended betas
