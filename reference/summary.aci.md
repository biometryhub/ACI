# Summarise an assimilative causal inference result

Summarises an `aci` object into the quantities a reader needs to judge
the result rather than only to read it. These are the distribution of
the causal-information metric, where and when it peaks, and the
diagnostics that say whether the numerical assumptions held.

## Usage

``` r
# S3 method for class 'aci'
summary(object, ...)
```

## Arguments

- object:

  An `aci` object, as returned by
  [`aci()`](https://biometryhub.github.io/ACI/reference/aci.md).

- ...:

  Ignored, for compatibility with
  [`summary()`](https://rdrr.io/r/base/summary.html).

## Value

An object of class `summary.aci`, a list with the number of steps `n`,
the step `dt`, the time `span`, the model `label`, the metric
five-number summary `metric`, the peak metric value `peak` and its time
`peak_time`, the smallest covariances `min_filter_cov` and
`min_smoother_cov`, the `terminal_residual`, and the number of
round-off-clamped metric values `n_clamped`.

## Details

Three diagnostics are reported alongside the metric. The smallest
filtered and smoothed covariances indicate how close the explicit Euler
discretisation came to losing positivity; a value near zero means the
result is sensitive to `dt`. The terminal residual is the gap between
the smoother and the filter at the final step, which is zero
analytically and so measures accumulated numerical error. The clamp
count reports how many metric values were negative by no more than
round-off and were clamped to zero, which can happen where the two
posteriors agree to near machine precision. The exact zero at the final
step is the terminal identity at work, not a clamp, and is not counted.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications*, 17, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

## See also

[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md),
[`plot.aci()`](https://biometryhub.github.io/ACI/reference/plot.aci.md),
[`as.data.frame.aci()`](https://biometryhub.github.io/ACI/reference/as.data.frame.aci.md)

## Examples

``` r
model <- aci_dyad_model()
fit <- aci(aci_simulate(model, n = 5000, seed = 333)$x, model)
summary(fit)
#> <summary.aci> assimilative causal inference
#>   model: nonlinear dyad model with intermittent extreme events
#>   steps: 5000, dt: 0.001, time span: [0, 4.999]
#> 
#>   causal-information metric:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.00000 0.06962 0.11010 0.22070 0.25589 2.45911 
#>   peak 2.45911 at time 0.382
#> 
#>   diagnostics:
#>     smallest covariance: filter 0.1, smoother 0.0682916
#>     terminal identity residual: 0 (zero analytically)
#>     round-off clamps to zero: 0 of 5000 steps
```
