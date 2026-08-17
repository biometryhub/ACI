# Coerce a causal influence range to a data frame

The per-time quantities only. The subjective matrix is left behind
because it is threshold-by-time and does not belong in the same
rectangle.

## Usage

``` r
# S3 method for class 'aci_cir'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An object of class `aci_cir`.

- row.names, optional:

  Passed through to the data frame.

- ...:

  Ignored.

## Value

A data frame with one row per reported time.

## Details

`monotone` is carried because it is the condition under which
`objective` and `objective_exact` are the same functional, and it is
what
[`summary.aci_cir()`](https://biometryhub.github.io/ACI/reference/summary.aci_cir.md)
directs the reader to. A diagnostic reachable only through a print
method cannot be used programmatically.

## Examples

``` r
model <- aci_dyad_model()
sim <- aci_simulate(model, n = 900, seed = 1)
comp <- aci_dyad_components(sim$x, model$parameters)
head(as.data.frame(aci_cir(sim$x, comp, dt = 0.001, window = 50:150,
                           mu0 = model$y0, R0 = 0.1)))
#>    time index  objective objective_exact      peak   status monotone saturated
#> 1 0.049    50 0.05294826      0.09658520 0.2906256 resolved    FALSE     FALSE
#> 2 0.050    51 0.05266212      0.09470754 0.2937984 resolved    FALSE     FALSE
#> 3 0.051    52 0.05263072      0.09522146 0.2969950 resolved    FALSE     FALSE
#> 4 0.052    53 0.05232731      0.09339877 0.3002157 resolved    FALSE     FALSE
#> 5 0.053    54 0.05230943      0.09274205 0.3034609 resolved    FALSE     FALSE
#> 6 0.054    55 0.05200586      0.09496806 0.3067306 resolved    FALSE     FALSE
```
