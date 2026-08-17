# Stochastic ENSO model

Returns the ENSO model as a labelled description of its structure and
parameters. Unlike the scalar models, this one carries no simulator. The
six-dimensional system is integrated with a mixed scheme in the
reference implementation, and its wind-burst update contains a
correction this package does not reproduce. See the design note on that
anomaly before simulating the system by any route.

## Usage

``` r
aci_enso_model(p = aci_enso_parameters())
```

## Arguments

- p:

  A parameter list, as returned by
  [`aci_enso_parameters()`](https://biometryhub.github.io/ACI/reference/aci_enso_parameters.md).

## Value

A list describing the model, with its `label`, `parameters`, the names
of its `observed` and `unobserved` components, and the `components`
function that builds a components list from observed paths.

## Details

Use
[`aci_enso_components()`](https://biometryhub.github.io/ACI/reference/aci_enso_components.md)
with observed paths to obtain the components consumed by
[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md),
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md)
and
[`aci_metric()`](https://biometryhub.github.io/ACI/reference/aci_metric.md),
and
[`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md)
to ask the causal question of one observable given the others.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications*, 17, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

## See also

[`aci_enso_components()`](https://biometryhub.github.io/ACI/reference/aci_enso_components.md),
[`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md)

## Examples

``` r
model <- aci_enso_model()
model$observed
#> [1] "T_C" "T_E" "I"  
model$unobserved
#> [1] "u"   "h_W" "tau"
```
