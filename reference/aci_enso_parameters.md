# Parameters of the stochastic ENSO model

Builds the parameter list of the ENSO model, with the defaults of the
reference implementation. Most parameters are derived from a smaller set
of physical ones, so they are computed here rather than left for a
caller to keep consistent.

## Usage

``` r
aci_enso_parameters(
  factor = 0.65,
  b_0 = 2.5,
  mu = 0.5,
  d_tau = 2,
  lambda = 2/60,
  m = 2
)
```

## Arguments

- factor:

  Numeric scalar. The diversity-modulating parameter that governs how
  many extreme eastern-Pacific events the model produces.

- b_0:

  Numeric scalar. High-end estimate of the thermocline tilt.

- mu:

  Numeric scalar. Relative coupling coefficient.

- d_tau:

  Numeric scalar. Damping of the wind burst.

- lambda:

  Numeric scalar. Damping of the Walker circulation.

- m:

  Numeric scalar. Target equilibrium mean of the Walker circulation.

## Value

A named list of the model's parameters, derived and primitive.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications*, 17, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

## See also

[`aci_enso_components()`](https://biometryhub.github.io/ACI/reference/aci_enso_components.md),
[`aci_enso_model()`](https://biometryhub.github.io/ACI/reference/aci_enso_model.md)

## Examples

``` r
str(aci_enso_parameters())
#> List of 22
#>  $ factor : num 0.65
#>  $ b_0    : num 2.5
#>  $ mu     : num 0.5
#>  $ d_tau  : num 2
#>  $ lambda : num 0.0333
#>  $ m      : num 2
#>  $ alpha_1: num 0.0264
#>  $ alpha_2: num 0.0813
#>  $ delta_u: num 0.033
#>  $ delta_h: num 0.102
#>  $ r      : num 0.163
#>  $ gamma_C: num 0.488
#>  $ gamma_E: num 0.488
#>  $ r_C    : num 0.305
#>  $ r_E    : num 0.914
#>  $ zeta_C : num 0.305
#>  $ zeta_E : num 0.305
#>  $ C_u    : num 0.0195
#>  $ sigma_u: num 0.0322
#>  $ sigma_h: num 0.0161
#>  $ sigma_C: num 0.0322
#>  $ sigma_E: num 0.018
```
