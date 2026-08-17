# Plot an assimilative causal inference result

Draws the observed signal and the causal-information metric on a shared
time axis. The upper panel shows bursts in the observed signal. The
lower panel shows the steps at which the future of that signal sharpens
the estimate of the unobserved component.

## Usage

``` r
# S3 method for class 'aci'
plot(x, y, ...)
```

## Arguments

- x:

  An `aci` object, as returned by
  [`aci()`](https://biometryhub.github.io/ACI/reference/aci.md).

- y:

  Ignored, for compatibility with
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

- ...:

  Further graphical parameters passed to the underlying
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) calls.

## Value

The result `x`, invisibly. Called for the plot it draws.

## Details

The method uses base graphics so that plotting costs the package no
dependency. The colours are chosen to remain distinguishable in the
common forms of colour vision deficiency. For a publication-grade
figure, take
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) of the
result and draw it with the grammar of graphics; the *Assimilative
causal inference on the nonlinear dyad model* vignette shows the recipe.

## See also

[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md),
[`as.data.frame.aci()`](https://biometryhub.github.io/ACI/reference/as.data.frame.aci.md)

## Examples

``` r
model <- aci_dyad_model()
fit <- aci(aci_simulate(model, n = 5000, seed = 333)$x, model)
plot(fit)

```
