# Noisy predator-prey benchmark model

andreou2026aci supplementary model (SI.4.2; ACI_code
noisy_predator_prey) Construct either causal partition of the stochastic
Lotka-Volterra example used in andreou2026aci. The two partitions should
be compared separately because the supplied MATLAB file contains
sequential direction-specific blocks.

## Usage

``` r
aci_predprey_model(hidden = c("prey", "predator"), params = list())
```

## Arguments

- hidden:

  Either `"prey"` or `"predator"`, naming the hidden component.

- params:

  Optional named list overriding `alpha`, `beta`, `gamma`, `delta`,
  `s_x` and `s_y`.

## Value

An object of class `cgns_model`.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications* **17**, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

## See also

[`aci_dyad_model()`](https://biometryhub.github.io/ACI/reference/aci_dyad_model.md),
[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md)

## Examples

``` r
aci_predprey_model(hidden = "prey")
#> <cgns_model> 'predator_prey[prey hidden]': k = 1 observed, l = 1 hidden
```
