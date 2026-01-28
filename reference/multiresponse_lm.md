# Fit Multiresponse Linear Model

This function fits a linear model to multiple responses in an fMRI
dataset.

## Usage

``` r
multiresponse_lm(form, data_env, conlist, vnames, fcon, modmat = NULL)
```

## Arguments

- form:

  The formula used to define the linear model.

- data_env:

  The environment containing the data to be used in the linear model.

- conlist:

  The list of contrasts used in the analysis.

- vnames:

  The names of the variables used in the linear model.

- fcon:

  The F-contrasts used in the analysis.

- modmat:

  The model matrix (default is `NULL`, which will calculate the model
  matrix using the formula).

## Value

A list containing the results from the multiresponse linear model
analysis.
