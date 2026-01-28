# Create 3D Blocks for Voxelwise Analysis

Creates spatial blocks from a 3D mask for use with spatial_fdr. This is
useful for voxelwise analyses where you want to group nearby voxels
together for more powerful multiple comparisons correction.

## Usage

``` r
create_3d_blocks(mask, block_size = c(10, 10, 10), connectivity = 26)
```

## Arguments

- mask:

  Numeric 3D array or NeuroVol object defining the brain mask. Non-zero
  values indicate voxels to include.

- block_size:

  Integer vector of length 3 specifying block dimensions in voxels
  (default: c(10, 10, 10))

- connectivity:

  Integer scalar; type of connectivity for neighbors: 6 (face
  connectivity) or 26 (face, edge, and corner connectivity). Default: 26

## Value

List with components:

- group_id:

  Integer vector of group IDs for each voxel in the mask

- neighbors:

  List of length n_groups where element i contains integer vector of
  1-based neighbor IDs for group i

- n_groups:

  Integer scalar; total number of groups created

- block_size:

  Block size used

## Examples

``` r
mask <- array(c(1L, 1L, 0L,
                1L, 1L, 0L,
                0L, 0L, 0L), dim = c(3, 3, 1))
blocks <- create_3d_blocks(mask, block_size = c(2, 2, 1))
blocks$n_groups
#> [1] 1
```
