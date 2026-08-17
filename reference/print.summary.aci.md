# Print a summary of an assimilative causal inference result

Print a summary of an assimilative causal inference result

## Usage

``` r
# S3 method for class 'summary.aci'
print(x, ...)
```

## Arguments

- x:

  A `summary.aci` object, as returned by
  [`summary.aci()`](https://biometryhub.github.io/ACI/reference/summary.aci.md).

- ...:

  Ignored, for compatibility with
  [`print()`](https://rdrr.io/r/base/print.html).

## Value

The summary `x`, invisibly.

## See also

[`summary.aci()`](https://biometryhub.github.io/ACI/reference/summary.aci.md)

## Examples

``` r
model <- aci_dyad_model()
print(summary(aci(aci_simulate(model, n = 2000, seed = 1)$x, model)))
#> <summary.aci> assimilative causal inference
#>   model: nonlinear dyad model with intermittent extreme events
#>   steps: 2000, dt: 0.001, time span: [0, 1.999]
#> 
#>   causal-information metric:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.00000 0.09512 0.18530 0.32462 0.41822 1.60235 
#>   peak 1.60235 at time 1.204
#> 
#>   diagnostics:
#>     smallest covariance: filter 0.1, smoother 0.0641812
#>     terminal identity residual: 0 (zero analytically)
#>     round-off clamps to zero: 0 of 2000 steps
```
