# Print an assimilative causal inference result

Prints a compact summary of an `aci` object, covering the model it was
run for, the number of steps, the integration step and time span, and a
five-number summary of the causal-information metric.

## Usage

``` r
# S3 method for class 'aci'
print(x, ...)
```

## Arguments

- x:

  An `aci` object, as returned by
  [`aci()`](https://biometryhub.github.io/ACI/reference/aci.md).

- ...:

  Ignored, for compatibility with
  [`print()`](https://rdrr.io/r/base/print.html).

## Value

The result `x`, invisibly.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications*, 17, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

## See also

[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md),
[`summary.aci()`](https://biometryhub.github.io/ACI/reference/summary.aci.md)

## Examples

``` r
model <- aci_dyad_model()
fit <- aci(aci_simulate(model, n = 2000, seed = 1)$x, model)
print(fit)
#> <aci> assimilative causal inference
#>   model: nonlinear dyad model with intermittent extreme events
#>   steps: 2000, dt: 0.001, time span: [0, 1.999]
#>   causal-information metric:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.00000 0.09512 0.18530 0.32462 0.41822 1.60235 
```
