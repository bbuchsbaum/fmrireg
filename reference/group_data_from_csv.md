# Create Group Dataset from CSV File or Data Frame

Creates a group dataset from tabular data containing pre-extracted
statistics such as ROI means, effect sizes, and standard errors. This
format is useful for ROI-based analyses or when working with summary
statistics.

## Usage

``` r
group_data_from_csv(
  data,
  effect_cols,
  subject_col = "subject",
  roi_col = NULL,
  contrast_col = NULL,
  covariate_cols = NULL,
  wide_format = FALSE
)
```

## Arguments

- data:

  Either a path to a CSV file or a data frame containing the data

- effect_cols:

  Named vector or list specifying column names for effect statistics.
  E.g., c(beta = "mean_activation", se = "std_error") or c(t = "t_stat",
  df = "df")

- subject_col:

  Character string specifying the column containing subject IDs

- roi_col:

  Character string specifying the column containing ROI names (optional)

- contrast_col:

  Character string specifying the column containing contrast names
  (optional)

- covariate_cols:

  Character vector of column names to use as covariates (optional)

- wide_format:

  Logical. If TRUE, expects wide format with ROIs as columns (default:
  FALSE)

## Value

A group_data_csv object

## Examples

``` r
if (FALSE) { # \dontrun{
# Long format: one row per subject-ROI combination
gd <- group_data_from_csv(
  "roi_statistics.csv",
  effect_cols = c(beta = "mean_beta", se = "se_beta"),
  subject_col = "participant_id",
  roi_col = "roi_name",
  covariate_cols = c("age", "sex", "group")
)

# Wide format: one row per subject, ROIs as columns
gd <- group_data_from_csv(
  "subject_summary.csv",
  effect_cols = c(beta = "roi_"),  # Prefix for ROI columns
  subject_col = "subject",
  wide_format = TRUE
)

# From data frame with multiple contrasts
df <- read.csv("contrast_results.csv")
gd <- group_data_from_csv(
  df,
  effect_cols = c(beta = "estimate", se = "std_error", t = "t_value"),
  subject_col = "subject_id",
  contrast_col = "contrast_name",
  roi_col = "region"
)
} # }
```
