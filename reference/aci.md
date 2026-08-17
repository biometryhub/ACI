# Run assimilative causal inference

Runs the whole assimilative causal inference workflow on an observed
signal for a given `aci_model`. It assembles the conditional Gaussian
components from the model and the signal, runs the forward filter and
backward smoother, and scores each step with the causal-information
metric. This is the single user-facing entry point that ties the
numerical core together.

## Usage

``` r
aci(x, model, dt = 0.001, mu0 = NULL, R0 = 0.1, time = NULL)
```

## Arguments

- x:

  Numeric vector. The observed signal, one complete finite value per
  time step; at least two.

- model:

  An `aci_model` object; see
  [`aci_cgns_model()`](https://biometryhub.github.io/ACI/reference/aci_cgns_model.md).

- dt:

  Numeric scalar. The integration time step; must be positive. Ignored
  with a check when `time` is supplied.

- mu0:

  Numeric scalar. The initial filtered mean of the unobserved component.
  Defaults to the model's initial unobserved state `y0`.

- R0:

  Numeric scalar. The initial filtered covariance of the unobserved
  component; must be positive.

- time:

  Optional numeric vector. The observed time grid, one value per
  observation, strictly increasing and equally spaced. When supplied,
  `dt` is derived from it.

## Value

An `aci` object, a list with the `model`, the time vector `t`, the
observed signal `x`, the `filter` and `smoother` statistics (each a list
of `mean` and `cov`), the causal-information metric `aci`, the count
`n_clamped` of metric values clamped from a round-off negative to zero,
and the step `dt`. See
[`summary.aci()`](https://biometryhub.github.io/ACI/reference/summary.aci.md),
[`plot.aci()`](https://biometryhub.github.io/ACI/reference/plot.aci.md)
and
[`as.data.frame.aci()`](https://biometryhub.github.io/ACI/reference/as.data.frame.aci.md).

## Details

The workflow assumes the observed signal is complete and sampled on a
regular grid. The closed-form recursions integrate a fixed step, and
there is no contract for missing observations or for an irregular grid.
By default the time vector is constructed from `dt`; supply `time`
instead to have the step derived from an observed grid and the
regularity checked.

The metric this returns is a statement about the supplied model, not
about the world. See the *Assumptions and interpretation* vignette for
what an ACI peak does and does not support.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications*, 17, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

## See also

[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md),
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md),
[`aci_metric()`](https://biometryhub.github.io/ACI/reference/aci_metric.md),
[`summary.aci()`](https://biometryhub.github.io/ACI/reference/summary.aci.md)

## Examples

``` r
model <- aci_dyad_model()
sim <- aci_simulate(model, n = 5000, seed = 333)
fit <- aci(sim$x, model)
fit
#> <aci> assimilative causal inference
#>   model: nonlinear dyad model with intermittent extreme events
#>   steps: 5000, dt: 0.001, time span: [0, 4.999]
#>   causal-information metric:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.00000 0.06962 0.11010 0.22070 0.25589 2.45911 

# An observed time grid may be supplied instead of a step.
fit_time <- aci(sim$x, model, time = sim$t)
identical(fit$aci, fit_time$aci)
#> [1] TRUE
```
