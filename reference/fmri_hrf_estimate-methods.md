# Methods for smooth HRF estimates

Methods for smooth HRF estimates

## Usage

``` r
# S3 method for class 'fmri_hrf_estimate'
print(x, ...)

# S3 method for class 'fmri_hrf_estimate'
coef(object, ...)

# S3 method for class 'fmri_hrf_estimate'
as.matrix(x, curve = NULL, what = c("estimate", "std.error"), ...)

# S3 method for class 'fmri_hrf_estimate'
predict(object, newdata = object$time, se.fit = FALSE, ...)
```

## Arguments

- x, object:

  An `fmri_hrf_estimate` object.

- ...:

  Unused.

- curve:

  A curve name or one-based curve index. Required by
  [`as.matrix()`](https://rdrr.io/r/base/matrix.html) when the object
  contains multiple curves.

- what:

  Return estimated curves or their standard errors.

- newdata:

  Numeric post-stimulus times within the fitted span.

- se.fit:

  Return prediction standard errors with fitted curves.
