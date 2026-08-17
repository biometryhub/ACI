# Print a conditional Gaussian nonlinear system model

Prints a compact, one-field-per-line summary of an `aci_model` object,
covering its label, the self-drift of the unobserved component, the
noise Grammians, the initial state and, when present, the named scalar
parameters that define the model.

## Usage

``` r
# S3 method for class 'aci_model'
print(x, ...)
```

## Arguments

- x:

  An `aci_model` object; see
  [`aci_cgns_model()`](https://biometryhub.github.io/ACI/reference/aci_cgns_model.md).

- ...:

  Ignored, for compatibility with
  [`print()`](https://rdrr.io/r/base/print.html).

## Value

The model `x`, invisibly.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications*, 17, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

## See also

[`aci_cgns_model()`](https://biometryhub.github.io/ACI/reference/aci_cgns_model.md),
[`aci_dyad_model()`](https://biometryhub.github.io/ACI/reference/aci_dyad_model.md)

## Examples

``` r
print(aci_dyad_model())
#> <aci_model> nonlinear dyad model with intermittent extreme events
#>   unobserved self-drift L_y: -0.5
#>   noise Grammians: S_xoS_x = 0.25, S_yoS_y = 1, S_yoS_x = 0
#>   initial state: x0 = 1, y0 = 2
#>   parameters: d_x = 0.5, d_y = 0.5, gamma = 2, F_x = 0.5, F_y = 1, sigma_x = 0.5, sigma_y = 1
```
