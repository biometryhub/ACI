# Coerce an assimilative causal inference result to a data frame

Flattens the posterior trajectories and the causal-information metric
into one tidy row per time step, for plotting, export or joining to
other series.

## Usage

``` r
# S3 method for class 'aci'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)
```

## Arguments

- x:

  An `aci` object, as returned by
  [`aci()`](https://biometryhub.github.io/ACI/reference/aci.md).

- row.names:

  Passed to [`data.frame()`](https://rdrr.io/r/base/data.frame.html).

- optional:

  Passed to [`data.frame()`](https://rdrr.io/r/base/data.frame.html).

- ...:

  Ignored, for compatibility with
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).

## Value

A data frame with one row per time step and the columns `t`, `x`,
`filter_mean`, `filter_cov`, `smoother_mean`, `smoother_cov` and `aci`.

## See also

[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md),
[`summary.aci()`](https://biometryhub.github.io/ACI/reference/summary.aci.md)

## Examples

``` r
model <- aci_dyad_model()
fit <- aci(aci_simulate(model, n = 2000, seed = 1)$x, model)
head(as.data.frame(fit))
#>       t         x filter_mean filter_cov smoother_mean smoother_cov        aci
#> 1 0.000 1.0000000    2.000000  0.1000000      1.972334   0.06583734 0.04200534
#> 2 0.001 0.9940949    1.990076  0.1007400      1.970049   0.06609311 0.04076791
#> 3 0.002 1.0009543    1.990428  0.1014788      1.967868   0.06634390 0.04189209
#> 4 0.003 0.9916614    1.977640  0.1022122      1.965637   0.06658970 0.04069774
#> 5 0.004 1.0206483    1.996005  0.1029456      1.963551   0.06683055 0.04572747
#> 6 0.005 1.0298176    1.998215  0.1036661      1.961148   0.06706643 0.04784566
```
