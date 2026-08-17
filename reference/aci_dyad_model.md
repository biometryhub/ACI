# Nonlinear dyad model

Builds the `aci_model` object for the nonlinear dyad model with
intermittent extreme events, the flagship conditional Gaussian nonlinear
system of the package. The observed process `x` couples to the
unobserved process `y` through the state-dependent term `gamma x`, which
makes `x` intermittently extreme, and the unobserved process has
constant self-drift `-d_y`.

## Usage

``` r
aci_dyad_model(
  d_x = 0.5,
  d_y = 0.5,
  gamma = 2,
  F_x = 0.5,
  F_y = 1,
  sigma_x = 0.5,
  sigma_y = 1
)
```

## Arguments

- d_x:

  Numeric scalar. Damping of the observed process.

- d_y:

  Numeric scalar. Damping of the unobserved process.

- gamma:

  Numeric scalar. Strength of the quadratic coupling.

- F_x:

  Numeric scalar. Constant forcing of the observed process.

- F_y:

  Numeric scalar. Constant forcing of the unobserved process.

- sigma_x:

  Numeric scalar. Noise coefficient of the observed process; must be
  non-zero, since its square is the observation-noise covariance.

- sigma_y:

  Numeric scalar. Noise coefficient of the unobserved process.

## Value

An `aci_model` object; see
[`aci_cgns_model()`](https://biometryhub.github.io/ACI/reference/aci_cgns_model.md).

## Details

The two governing equations are
`dx = (- d_x x + gamma x y + F_x) dt + sigma_x dW_x` and
`dy = (- d_y y - gamma x^2 + F_y) dt + sigma_y dW_y`, with independent
noise. The default parameters are those of the reference implementation.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications*, 17, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

## See also

[`aci_cgns_model()`](https://biometryhub.github.io/ACI/reference/aci_cgns_model.md),
[`aci_simulate()`](https://biometryhub.github.io/ACI/reference/aci_simulate.md),
[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md)

## Examples

``` r
model <- aci_dyad_model()
model
#> <aci_model> nonlinear dyad model with intermittent extreme events
#>   unobserved self-drift L_y: -0.5
#>   noise Grammians: S_xoS_x = 0.25, S_yoS_y = 1, S_yoS_x = 0
#>   initial state: x0 = 1, y0 = 2
#>   parameters: d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1, sigma_x = 0.5, sigma_y = 1
sim <- aci_simulate(model, n = 5000, seed = 333)
fit <- aci(sim$x, model)
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
