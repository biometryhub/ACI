# Conditional Gaussian components of the nonlinear dyad model

Builds the conditional Gaussian components list for the nonlinear dyad
model with intermittent extreme events, given an observed signal and the
model parameters. The returned list is the `comp` argument consumed by
[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md)
and
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md).

## Usage

``` r
aci_dyad_components(x, p)
```

## Arguments

- x:

  Numeric vector. The observed signal, one value per time step.

- p:

  A named list of dyad parameters with elements `d_x`, `d_y`, `gamma`,
  `F_x`, `F_y`, `sigma_x` and `sigma_y`, each a finite numeric scalar.
  The `parameters` entry of an
  [`aci_dyad_model()`](https://biometryhub.github.io/ACI/reference/aci_dyad_model.md)
  is exactly this list.

## Value

A conditional Gaussian components list; see
[aci_components](https://biometryhub.github.io/ACI/reference/aci_components.md).

## Details

In the dyad model the observed process couples to the unobserved process
through the state-dependent term `gamma * x`, and the unobserved process
has constant linear self-drift `-d_y`. The two noise sources are
independent, so the noise cross-covariances are zero.

Most users should reach for
[`aci_dyad_model()`](https://biometryhub.github.io/ACI/reference/aci_dyad_model.md)
and [`aci()`](https://biometryhub.github.io/ACI/reference/aci.md)
instead. This function is the worked example of the components schema,
and it serves as the template for a conditional Gaussian system the
package supplies no constructor for. See
[aci_components](https://biometryhub.github.io/ACI/reference/aci_components.md)
for the schema itself.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications*, 17, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

## See also

[aci_components](https://biometryhub.github.io/ACI/reference/aci_components.md),
[`aci_dyad_model()`](https://biometryhub.github.io/ACI/reference/aci_dyad_model.md),
[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md)

## Examples

``` r
model <- aci_dyad_model()
sim <- aci_simulate(model, n = 2000, seed = 1)
comp <- aci_dyad_components(sim$x, model$parameters)
str(comp)
#> List of 8
#>  $ L_x    : num [1:2000] 2 1.99 2 1.98 2.04 ...
#>  $ f_x    : num [1:2000] 0 0.002953 -0.000477 0.004169 -0.010324 ...
#>  $ L_y    : num -0.5
#>  $ f_y    : num [1:2000] -1 -0.976 -1.004 -0.967 -1.083 ...
#>  $ S_xoS_x: num 0.25
#>  $ S_yoS_y: num 1
#>  $ S_yoS_x: num 0
#>  $ S_xoS_y: num 0
```
