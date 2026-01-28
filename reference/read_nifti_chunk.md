# Read Data Chunk from NIfTI Files

Reads a subset of voxels from all subjects' NIfTI files

## Usage

``` r
read_nifti_chunk(gd, voxel_indices)
```

## Arguments

- gd:

  A group_data_nifti object

- voxel_indices:

  Integer vector of voxel indices to read

## Value

Matrix with dimensions (subjects, voxels)
