# Plot a causal influence range

Draws two panels. The upper panel shows the objective range against
time, with censored times drawn hollow so that a lower bound is not
mistaken for a measurement. The lower panel shows the peak divergence.

## Usage

``` r
# S3 method for class 'aci_cir'
plot(x, y, ...)
```

## Arguments

- x:

  An object of class `aci_cir`.

- y:

  Ignored.

- ...:

  Passed to
  [`graphics::plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

`x`, invisibly.

## Examples

``` r
model <- aci_dyad_model()
sim <- aci_simulate(model, n = 900, seed = 1)
comp <- aci_dyad_components(sim$x, model$parameters)
rng <- aci_cir(sim$x, comp, dt = 0.001, window = 50:150,
               mu0 = model$y0, R0 = 0.1)
plot(rng)

```
