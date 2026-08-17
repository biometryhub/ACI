# Noisy predator-prey model

Builds the stochastic Lotka-Volterra predator-prey model as an
`aci_model` object, ready for
[`aci_simulate()`](https://biometryhub.github.io/ACI/reference/aci_simulate.md)
and [`aci()`](https://biometryhub.github.io/ACI/reference/aci.md). The
parameter defaults are those of the reference implementation.

## Usage

``` r
aci_predprey_model(
  direction = "prey_to_predator",
  alpha = 0.4,
  beta = 0.1,
  gamma = 1.1,
  delta = 0.4,
  sigma_x = 0.3,
  sigma_y = 0.3,
  x0 = 4,
  y0 = 4
)
```

## Arguments

- direction:

  Character scalar, `"prey_to_predator"` (the default) or
  `"predator_to_prey"`. See
  [`aci_predprey_components()`](https://biometryhub.github.io/ACI/reference/aci_predprey_components.md).

- alpha:

  Numeric scalar. The predator's natural death rate.

- beta:

  Numeric scalar. The effect of prey availability on predator growth.

- gamma:

  Numeric scalar. The prey's natural growth rate.

- delta:

  Numeric scalar. The effect of predator presence on prey decline.

- sigma_x, sigma_y:

  Numeric scalars. The noise amplitudes of the predator and prey
  processes.

- x0, y0:

  Numeric scalars. Initial predator and prey populations.

## Value

An object of class `aci_model`.

## Details

The model is studied in two directions, and `direction` selects which.
The two are genuinely different questions rather than one question and
its mirror. The observed process differs, the latent process differs,
and the causal-information metric is not symmetric between them.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications*, 17, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

## See also

[`aci_dyad_model()`](https://biometryhub.github.io/ACI/reference/aci_dyad_model.md),
[`aci_predprey_components()`](https://biometryhub.github.io/ACI/reference/aci_predprey_components.md),
[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md)

## Examples

``` r
model <- aci_predprey_model()
model
#> <aci_model> noisy predator-prey model (prey to predator)
#>   unobserved self-drift L_y: state-dependent
#>   noise Grammians: S_xoS_x = 0.09, S_yoS_y = 0.09, S_yoS_x = 0
#>   initial state: x0 = 4, y0 = 4
#>   parameters: alpha = 0.4, beta = 0.1, gamma = 1.1, delta = 0.4, sigma_x = 0.3, sigma_y = 0.3
```
