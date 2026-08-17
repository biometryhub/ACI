# Assimilative causal-information metric

Computes the assimilative causal-information metric at each time step.
The metric is the relative entropy (Kullback-Leibler divergence) of the
smoother posterior of the unobserved component from its filter
posterior, with the smoother posterior as the integrating density. A
larger value marks a step at which conditioning on the future
observations sharpens the estimate of the unobserved component more
strongly, and the metric is non-negative.

## Usage

``` r
aci_metric(filt, smooth)
```

## Arguments

- filt:

  A list with numeric vectors `mean` and `cov`, the filtered mean and
  covariance, as returned by
  [`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md).

- smooth:

  A list with numeric vectors `mean` and `cov`, the smoothed mean and
  covariance, as returned by
  [`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md).

## Value

A numeric vector of the causal-information metric at each time step,
non-negative throughout.

## Details

For scalar conditional Gaussian posteriors the relative entropy
separates into a signal part, driven by the shift between the smoothed
and filtered means relative to the filtered covariance, and a dispersion
part, driven by the ratio of the smoothed to the filtered covariance.

The dispersion part is evaluated in a form that stays accurate when the
two covariances nearly agree. Writing the covariance ratio as `1 + d`,
the direct expression subtracts two nearly equal quantities and loses
precision exactly where the metric is smallest; the form used here does
not. Values that round to a small negative number are clamped to zero,
and anything more negative than round-off is an error rather than a
result. The number of values clamped this way is recorded by
[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md) and
reported by
[`summary.aci()`](https://biometryhub.github.io/ACI/reference/summary.aci.md);
a value that is exactly zero because the posteriors agree, as at the
final step, needs no clamp and is not counted.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications*, 17, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

## See also

[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md),
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md),
[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md)

## Examples

``` r
model <- aci_dyad_model()
sim <- aci_simulate(model, n = 2000, seed = 1)
comp <- aci_dyad_components(sim$x, model$parameters)
filt <- aci_filter(sim$x, comp, dt = 0.001, mu0 = model$y0, R0 = 0.1)
smooth <- aci_smoother(sim$x, comp, dt = 0.001, filt)
summary(aci_metric(filt, smooth))
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.00000 0.09512 0.18530 0.32462 0.41822 1.60235 
```
