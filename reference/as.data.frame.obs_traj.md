# Coerce an observed trajectory to a data frame

Coerce an observed trajectory to a data frame

## Usage

``` r
# S3 method for class 'obs_traj'
as.data.frame(x, ...)
```

## Arguments

- x:

  An `obs_traj` object.

- ...:

  Ignored, for consistency with the generic.

## Value

A data frame in long form with columns `t`, `var` and `value`.
