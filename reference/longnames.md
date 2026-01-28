# Extract Long Names of Variable Levels

Get the extended names of variable levels, which include the term prefix
and any basis function information. Long names provide the complete
specification of each condition in the model. For example, if a term has
conditions "level1" and "level2" with basis functions "basis1" and
"basis2", the long names would be "term#level1:basis1",
"term#level1:basis2", "term#level2:basis1", "term#level2:basis2".

## Usage

``` r
longnames(x, ...)

# S3 method for class 'event_term'
longnames(x, ...)

# S3 method for class 'event_seq'
longnames(x, ...)

# S3 method for class 'convolved_term'
longnames(x, ...)

# S3 method for class 'event_model'
longnames(x, ...)
```

## Arguments

- x:

  The object to extract names from (typically an event_term,
  event_model, or convolved_term)

- ...:

  Additional arguments passed to methods. Common arguments include:

  exclude_basis

  :   Logical; if TRUE, exclude basis function labels from names

  drop_empty

  :   Logical; if TRUE, drop empty condition levels

## Value

A character vector containing the full condition names with term
prefixes and basis functions

## See also

[`shortnames()`](https://bbuchsbaum.github.io/fmrireg/reference/shortnames.md),
[`event_model()`](https://bbuchsbaum.github.io/fmridesign/reference/event_model.html),
[`event_term()`](https://bbuchsbaum.github.io/fmridesign/reference/event_term.html)

## Examples

``` r
# Create example data with multiple conditions
event_data <- data.frame(
  condition = factor(c("A", "B", "C", "A", "B", "C")),
  rt = c(0.8, 1.2, 0.9, 1.1, 0.7, 1.3),
  onsets = c(1, 10, 20, 30, 40, 50),
  run = c(1, 1, 1, 1, 1, 1)
)

# Create sampling frame
sframe <- sampling_frame(blocklens = 60, TR = 2)

# Create event model with multiple basis functions
evmodel <- event_model(
  onsets ~ hrf(condition, basis = "fourier", nbasis = 2),
  data = event_data,
  block = ~run,
  sampling_frame = sframe
)

# Get long names including basis functions
lnames <- longnames(evmodel)
# Returns: c("condition#A:basis1", "condition#A:basis2",
#           "condition#B:basis1", "condition#B:basis2",
#           "condition#C:basis1", "condition#C:basis2")

# Create simple event term
eterm <- event_term(
  list(condition = event_data$condition),
  onsets = event_data$onsets,
  blockids = event_data$run
)

# Get long names for term
term_names <- longnames(eterm)
# Returns: c("condition#A", "condition#B", "condition#C")
```
