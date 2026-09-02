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

  Passed to the underlying base-graphics calls.

## Value

`x`, invisibly; called for the plot it draws.
