# Coerce a Gaussian assimilation path to a data frame

Coerce a Gaussian assimilation path to a data frame

## Usage

``` r
# S3 method for class 'da_path_gaussian'
as.data.frame(x, ...)
```

## Arguments

- x:

  A `da_path_gaussian` object.

- ...:

  Ignored, for consistency with the generic.

## Value

A data frame in long form with one row per time and hidden component,
carrying the mean and its marginal variance.
