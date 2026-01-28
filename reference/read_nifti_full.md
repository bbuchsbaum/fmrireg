# Read All Data from NIfTI Files

Reads complete data from all subjects' NIfTI files using memory mapping
when possible.

## Usage

``` r
read_nifti_full(gd, use_mask = NULL)
```

## Arguments

- gd:

  A group_data_nifti object

- use_mask:

  Logical. Apply mask to data (default: TRUE if mask exists)

## Value

List with data matrices
