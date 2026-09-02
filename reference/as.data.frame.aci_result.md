# aci() result to data.frame

aci() result to data.frame

## Usage

``` r
# S3 method for class 'aci_result'
as.data.frame(x, ...)
```

## Arguments

- x:

  An `aci_result` object.

- ...:

  Ignored, for consistency with the generic.

## Value

A data frame with one row per time, carrying the metric and, when
retained, its signal and dispersion parts.
