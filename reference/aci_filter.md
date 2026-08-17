# Forward conditional Gaussian filter

Runs the forward filter of a conditional Gaussian nonlinear system.
Given an observed signal and the system components, it returns the
filtered mean and covariance of the unobserved component at each time
step, that is the mean and covariance of the unobserved state
conditional on the observed path up to and including the current step.

## Usage

``` r
aci_filter(x, comp, dt, mu0, R0)
```

## Arguments

- x:

  Numeric vector. The observed signal, one value per time step; at least
  two complete, finite observations.

- comp:

  A conditional Gaussian components list; see
  [aci_components](https://biometryhub.github.io/ACI/reference/aci_components.md).

- dt:

  Numeric scalar. The integration time step; must be positive.

- mu0:

  Numeric scalar. The initial filtered mean of the unobserved component.

- R0:

  Numeric scalar. The initial filtered covariance of the unobserved
  component; must be positive.

## Value

A list with two numeric vectors, `mean` and `cov`, the filtered mean and
covariance of the unobserved component at each time step.

## Details

The recursion is the closed-form conditional Gaussian filter and is
integrated with a first-order (Euler) step of width `dt`. The observed
signal is assumed to be complete and sampled on a regular grid of
spacing `dt`.

The filtered covariance is checked at every step. An explicit Euler
scheme can drive the covariance non-positive when `dt` is too large for
the system, even for a perfectly admissible model, and the relative
entropy that scores the result has no meaning in that state. The
recursion therefore stops at the first offending step and names it,
rather than returning a trajectory that looks like a result.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications*, 17, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

## See also

[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md),
[`aci_metric()`](https://biometryhub.github.io/ACI/reference/aci_metric.md),
[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md)

## Examples

``` r
model <- aci_dyad_model()
sim <- aci_simulate(model, n = 2000, seed = 1)
comp <- aci_dyad_components(sim$x, model$parameters)
filt <- aci_filter(sim$x, comp, dt = 0.001, mu0 = model$y0, R0 = 0.1)
str(filt)
#> List of 2
#>  $ mean: num [1:2000] 2 1.99 1.99 1.98 2 ...
#>  $ cov : num [1:2000] 0.1 0.101 0.101 0.102 0.103 ...
```
