# Coerce a lag table to a data frame

Coerce a lag table to a data frame

## Usage

``` r
# S3 method for class 'lag_table'
as.data.frame(x, ...)
```

## Arguments

- x:

  A `lag_table` object.

- ...:

  Ignored, for consistency with the generic.

## Value

A data frame in long form with one row per retained cell, carrying the
anchor time, the lag and the divergence.
