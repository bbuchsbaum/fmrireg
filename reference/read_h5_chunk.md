# Read Data Chunk from HDF5 Files

Reads a subset of voxels from all subjects' HDF5 files

## Usage

``` r
read_h5_chunk(gd, voxel_indices, stat = NULL)
```

## Arguments

- gd:

  A group_data_h5 object

- voxel_indices:

  Integer vector of voxel indices to read

- stat:

  Character vector of statistics to extract

## Value

List with one element per subject, each containing extracted data
