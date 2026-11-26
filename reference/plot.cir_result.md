# Plot a causal influence range result

Plot a causal influence range result

## Usage

``` r
# S3 method for class 'cir_result'
plot(x, ...)
```

## Arguments

- x:

  A `cir_result` object.

- ...:

  Further arguments to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html). A `type`,
  `xlab`, `ylab` or `main` given here replaces the default. A range with
  a single anchor is drawn as a title only, so only its `main`, `sub`,
  `xlab` and `ylab` are used, by
  [`title()`](https://rdrr.io/r/graphics/title.html). A range with more
  than one anchor whose values are all `NA`, as when every anchor is
  masked, needs a `ylim` that is not `NULL`.

## Value

`x`, invisibly; called for the plot it draws.
