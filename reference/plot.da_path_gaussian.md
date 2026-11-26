# Plot a Gaussian assimilation path

Plot a Gaussian assimilation path

## Usage

``` r
# S3 method for class 'da_path_gaussian'
plot(x, component = 1, truth = NULL, ...)
```

## Arguments

- x:

  A `da_path_gaussian` object.

- component:

  Integer index of the hidden component to draw.

- truth:

  Optional numeric vector of true hidden values to overlay.

- ...:

  Further arguments to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html), which draws
  the mean. A `type`, `xlab`, `ylab`, `ylim` or `main` given here
  replaces the default.

## Value

`x`, invisibly; called for the plot it draws.
