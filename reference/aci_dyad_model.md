# Nonlinear dyad (andreou2026aci eq. 1-2)

Nonlinear dyad (andreou2026aci eq. 1-2)

## Usage

``` r
aci_dyad_model(variant = "p1", observe = "x", params = NULL)
```

## Arguments

- variant:

  Paper-specific parameter preset.

- observe:

  Which dyad component is treated as observed.

- params:

  Optional complete parameter list overriding the preset.

## Value

A `cgns_model` for the dyad, ready for
[`simulate()`](https://rdrr.io/r/stats/simulate.html) and the
assimilation verbs; the observed-x variant carries a sealed realiser
descriptor so that compilation takes the directed route.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications* **17**, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

## Examples

``` r
aci_dyad_model()
#> <cgns_model> 'dyad[p1] y->x': k = 1 observed, l = 1 hidden
```
