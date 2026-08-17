# Conditional Gaussian components of the stochastic ENSO model

Builds the vector-valued components list of the ENSO model from the
observed paths, for the configuration in which the central- and
eastern-Pacific temperature anomalies and the Walker circulation are
observed and the zonal current, thermocline depth and wind burst are
not.

## Usage

``` r
aci_enso_components(
  T_C,
  T_E,
  I,
  p = aci_enso_parameters(),
  time = NULL,
  dt = 0.005
)
```

## Arguments

- T_C, T_E, I:

  Numeric vectors of equal length. The observed central-Pacific
  temperature anomaly, eastern-Pacific temperature anomaly, and Walker
  circulation strength.

- p:

  A parameter list, as returned by
  [`aci_enso_parameters()`](https://biometryhub.github.io/ACI/reference/aci_enso_parameters.md).

- time:

  Numeric vector or `NULL`. The observation times, needed because the
  coefficients are seasonally modulated. When `NULL`, a regular grid of
  spacing `dt` starting at zero is used.

- dt:

  Numeric scalar. The step used to build the default `time` grid.

## Value

A vector-valued conditional Gaussian components list; see
[aci_components](https://biometryhub.github.io/ACI/reference/aci_components.md).

## Details

The observed state is `x = (T_C, T_E, I)` and the unobserved state is
`y = (u, h_W, tau)`, both three-dimensional. Every coefficient varies
with the observed state, the season, or both.

Two noise covariances vary in time and neither is optional. The Walker
circulation's observation noise is multiplicative in its own state, with
variance `lambda * (4 - I) * I`, which vanishes at the ends of the
interval `[0, 4]` and so keeps the process inside it. The latent noise
of the wind burst depends on the observed central-Pacific temperature
and the season. A filter that read either at its first step alone would
be integrating a different system.

Because the Walker-circulation noise vanishes at `I = 0` and `I = 4`, an
observed path that reaches either endpoint gives a singular
observation-noise covariance. Such a path is rejected rather than
regularised. The filter inverts that covariance, and perturbing it to
restore invertibility would introduce an observation the data does not
support.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications*, 17, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

## See also

[`aci_enso_parameters()`](https://biometryhub.github.io/ACI/reference/aci_enso_parameters.md),
[`aci_enso_model()`](https://biometryhub.github.io/ACI/reference/aci_enso_model.md),
[`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md)

## Examples

``` r
p <- aci_enso_parameters()
n <- 200
tt <- seq(0, by = 0.005, length.out = n)
comp <- aci_enso_components(
  T_C = 0.05 * sin(tt), T_E = 0.04 * cos(tt),
  I = 1.6 + 0.2 * sin(0.5 * tt), p = p, time = tt
)
dim(comp$L_x)
#> [1]   3   3 200
```
