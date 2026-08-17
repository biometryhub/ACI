# Conditional Gaussian nonlinear system model

Builds an `aci_model` object for a conditional Gaussian nonlinear system
(CGNS). A CGNS pairs an observed process `x` with an unobserved process
`y` whose statistics, conditional on the observed path, are Gaussian in
closed form. The model is described by the coefficients of the pair of
governing stochastic differential equations, and this constructor is the
general entry point for supplying them.

## Usage

``` r
aci_cgns_model(
  L_x,
  f_x,
  L_y,
  f_y,
  S_xoS_x,
  S_yoS_y,
  S_yoS_x = 0,
  x0 = 0,
  y0 = 0,
  label = "conditional Gaussian nonlinear system",
  parameters = NULL
)
```

## Arguments

- L_x:

  Coupling of the unobserved component into the drift of the observed
  process. A vectorised function of the observed signal returning one
  value per observation, or a numeric scalar for a constant coupling.

- f_x:

  Remaining drift of the observed process. A vectorised function of the
  observed signal, or a numeric scalar.

- L_y:

  Linear self-drift of the unobserved component. A vectorised function
  of the observed signal, or a numeric scalar for a self-drift constant
  in time.

- f_y:

  Remaining drift of the unobserved process. A vectorised function of
  the observed signal, or a numeric scalar.

- S_xoS_x:

  Numeric scalar. The observation-noise covariance; must be positive.

- S_yoS_y:

  Numeric scalar. The latent-noise covariance of the unobserved process;
  must be non-negative.

- S_yoS_x:

  Numeric scalar. The latent-to-observation noise cross-covariance. The
  default `0` is independent noise.

- x0:

  Numeric scalar. The initial value of the observed process, used when
  simulating.

- y0:

  Numeric scalar. The initial value of the unobserved process, used when
  simulating and as the default initial filtered mean.

- label:

  Character string. A human-readable name for the model.

- parameters:

  Optional named list of the finite numeric scalars that define the
  model, retained for printing and reproducibility.

## Value

An `aci_model` object, a list of the model coefficients, noise Grammians
and initial state, with class `"aci_model"`.

## Details

The observed process obeys `dx = (L_x(x) y + f_x(x)) dt + dW_x` and the
unobserved process obeys `dy = (L_y(x) y + f_y(x)) dt + dW_y`, where the
drift of the unobserved component is linear in `y`. All four
coefficients may vary with the observed signal and are supplied as
vectorised functions of `x`; a numeric constant is accepted and treated
as constant in time.

A self-drift `L_y` that varies with the observed signal keeps the system
conditionally Gaussian, because the coefficient is measurable with
respect to the observed path. This is what the noisy predator-prey model
needs, where the latent population's growth rate is set by the
population being watched. What would leave the class is a coefficient
depending on the *unobserved* component, and no such system can be built
here.

The noise is described by its Grammians, the entries of the noise
covariance of the pair of processes. `S_xoS_x` and `S_yoS_y` are the
observation-noise and latent-noise covariances, and `S_yoS_x` is the
noise cross-covariance, zero for independent noise. For the scalar
systems this package integrates the cross-covariance is symmetric, so
the transpose `S_xoS_y` of the components schema is derived rather than
supplied.

The noise covariance must be mathematically admissible, with `S_xoS_x`
strictly positive (the filter inverts it), `S_yoS_y` non-negative, and
the joint covariance positive semidefinite, that is
`S_xoS_x * S_yoS_y - S_yoS_x^2 >= 0`. A model that violates any of these
cannot be constructed. A singular system, whose determinant is zero,
describes perfectly correlated noise and is admissible, but it can drive
the filtered covariance toward zero;
[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md)
reports the step at which that happens rather than returning an
uninterpretable result.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications*, 17, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

## See also

[`aci_dyad_model()`](https://biometryhub.github.io/ACI/reference/aci_dyad_model.md),
[`aci_simulate()`](https://biometryhub.github.io/ACI/reference/aci_simulate.md),
[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md),
[aci_components](https://biometryhub.github.io/ACI/reference/aci_components.md)

## Examples

``` r
# A linear Ornstein-Uhlenbeck pair with constant coupling.
model <- aci_cgns_model(
  L_x = 1, f_x = function(x) -0.5 * x, L_y = -0.5, f_y = 0,
  S_xoS_x = 0.25, S_yoS_y = 1, y0 = 0
)
sim <- aci_simulate(model, n = 2000, dt = 0.01, seed = 1)
fit <- aci(sim$x, model, dt = 0.01)
summary(fit)
#> <summary.aci> assimilative causal inference
#>   model: conditional Gaussian nonlinear system
#>   steps: 2000, dt: 0.01, time span: [0, 19.99]
#> 
#>   causal-information metric:
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#> 0.00000 0.07284 0.17307 0.24688 0.33611 1.95375 
#>   peak 1.95375 at time 12.04
#> 
#>   diagnostics:
#>     smallest covariance: filter 0.1, smoother 0.0810466
#>     terminal identity residual: 0 (zero analytically)
#>     round-off clamps to zero: 0 of 2000 steps

# An inadmissible noise covariance is rejected at construction.
try(aci_cgns_model(
  L_x = 1, f_x = 0, L_y = -0.5, f_y = 0,
  S_xoS_x = 1, S_yoS_y = -1
))
#> Error : `S_yoS_y` must be non-negative; it is the latent-noise covariance.
```
