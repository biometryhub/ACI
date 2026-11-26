# Plot an ACI result

Plot an ACI result

## Usage

``` r
# S3 method for class 'aci_result'
plot(x, decompose = TRUE, ...)
```

## Arguments

- x:

  An `aci_result` object.

- decompose:

  `TRUE` to draw the signal and dispersion parts alongside the total.

- ...:

  Further arguments to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), which draws
  the total. A `type`, `xlab`, `ylab` or `main` given here replaces the
  default.

## Value

`x`, invisibly; called for the plot it draws.
