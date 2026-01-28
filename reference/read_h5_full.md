# Read All Data from HDF5 Files

Reads complete data from all subjects' HDF5 files. Warning: This can use
a lot of memory for whole-brain data.

## Usage

``` r
read_h5_full(gd, stat = NULL)
```

## Arguments

- gd:

  A group_data_h5 object

- stat:

  Character vector of statistics to extract

## Value

Array with dimensions (voxels, subjects, statistics)
