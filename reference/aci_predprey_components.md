# Conditional Gaussian components of the noisy predator-prey model

Builds the conditional Gaussian components list for a stochastic
Lotka-Volterra predator-prey pair, given an observed signal, the model
parameters, and which of the two causal directions is being studied.

## Usage

``` r
aci_predprey_components(x, p, direction)
```

## Arguments

- x:

  Numeric vector. The observed signal, one value per time step. Which
  population this is depends on `direction`.

- p:

  A named list with elements `alpha`, `beta`, `gamma`, `delta`,
  `sigma_x` and `sigma_y`, each a finite numeric scalar. The
  `parameters` entry of an
  [`aci_predprey_model()`](https://biometryhub.github.io/ACI/reference/aci_predprey_model.md)
  is exactly this list.

- direction:

  Character scalar. `"predator_to_prey"` observes the prey and treats
  the predator as latent, asking whether the predator influences the
  future of the prey. `"prey_to_predator"` is the converse.

## Value

A conditional Gaussian components list; see
[aci_components](https://biometryhub.github.io/ACI/reference/aci_components.md).

## Details

The system is \$\$\mathrm{d}x = (\beta x y - \alpha x)\\\mathrm{d}t +
\sigma_x\\\mathrm{d}W_x\$\$ \$\$\mathrm{d}y = (\gamma y - \delta x
y)\\\mathrm{d}t + \sigma_y\\\mathrm{d}W_y\$\$ with `x` the predator and
`y` the prey. Either variable may be taken as the observed one, and the
choice determines which causal question is asked.

In both directions the latent variable's self-drift depends on the
observed state, so `L_y` is a vector rather than a scalar. This is what
distinguishes the model from the dyad, whose latent self-drift is the
constant `-d_y`.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications*, 17, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

## See also

[aci_components](https://biometryhub.github.io/ACI/reference/aci_components.md),
[`aci_predprey_model()`](https://biometryhub.github.io/ACI/reference/aci_predprey_model.md),
[`aci_dyad_components()`](https://biometryhub.github.io/ACI/reference/aci_dyad_components.md)

## Examples

``` r
model <- aci_predprey_model("prey_to_predator")
sim <- aci_simulate(model, n = 2000, dt = 0.005, seed = 1)
comp <- aci_predprey_components(sim$x, model$parameters, model$direction)
str(comp)
#> List of 8
#>  $ L_x    : num [1:2000] 0.4 0.399 0.399 0.397 0.401 ...
#>  $ f_x    : num [1:2000] -1.6 -1.59 -1.6 -1.59 -1.6 ...
#>  $ L_y    : num [1:2000] -0.5 -0.495 -0.496 -0.489 -0.503 ...
#>  $ f_y    : num [1:2000] 0 0 0 0 0 0 0 0 0 0 ...
#>  $ S_xoS_x: num 0.09
#>  $ S_yoS_y: num 0.09
#>  $ S_yoS_x: num 0
#>  $ S_xoS_y: num 0
```
