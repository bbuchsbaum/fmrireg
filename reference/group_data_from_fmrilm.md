# Create Group Dataset from fmri_lm Objects

Creates a group dataset directly from a list of fitted fmri_lm objects

## Usage

``` r
group_data_from_fmrilm(
  lm_list,
  contrast = NULL,
  stat = c("beta", "se"),
  subjects = NULL,
  covariates = NULL
)
```

## Arguments

- lm_list:

  List of fmri_lm objects

- contrast:

  Character string specifying which contrast to extract

- stat:

  Character vector of statistics to extract

- subjects:

  Character vector of subject identifiers

- covariates:

  Data frame of subject-level covariates

## Value

A group_data_h5 object (in-memory variant)
