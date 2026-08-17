# Summarise a causal influence range

Reports the range alongside the two diagnostics a reader needs in order
to interpret it. The first is how many reported times are censored
rather than resolved. The second is how often the divergence sequence is
monotone, which matters because `objective` and `objective_exact` are
the same functional only when it is.

## Usage

``` r
# S3 method for class 'aci_cir'
summary(object, ...)

# S3 method for class 'summary.aci_cir'
print(x, ...)
```

## Arguments

- object:

  An object of class `aci_cir`.

- ...:

  Ignored.

- x:

  An object of class `summary.aci_cir`.

## Value

An object of class `summary.aci_cir`.

## Examples

``` r
model <- aci_dyad_model()
sim <- aci_simulate(model, n = 900, seed = 1)
comp <- aci_dyad_components(sim$x, model$parameters)
summary(aci_cir(sim$x, comp, dt = 0.001, window = 50:150,
                mu0 = model$y0, R0 = 0.1))
#> Causal influence range, summary
#> 
#>   reported times : 101 over [0.049, 0.149]
#>   thresholds     : 129
#> 
#>   resolution
#>     resolved         101
#>     censored         0
#>     below_threshold  0
#>     insufficient     0
#>     censored cells   5.8% of subjective cells
#> 
#>   ranges
#>     objective        0.02497  (median 0.04139, max 0.05295)
#>     objective_exact  median 0.06507
#>     peak divergence  median 0.394
#> 
#>   divergence monotone at 0% of times
#>   Where it is not, the two objective ranges are different
#>   functionals rather than two quadratures of one: the range is
#>   measured as a last exit, which exceeds the superlevel measure.
```
