# Gaussian relative entropy along a pair of paths

Evaluates the relative entropy at every time of two assimilation paths
on a common grid. Gaussian relative entropy is oriented as smoother
relative to filter.

## Usage

``` r
aci_metric(p, q, decompose = TRUE)
```

## Arguments

- p:

  A `da_path_gaussian` object, the integrating distribution.

- q:

  A `da_path_gaussian` object on the same time grid.

- decompose:

  `TRUE` to return the signal and dispersion parts alongside the total.

## Value

A data frame with the time column `t` and either `total` alone or
`total`, `signal` and `dispersion`.

## See also

[`aci_metric_pair()`](https://biometryhub.github.io/ACI/reference/aci_metric_pair.md),
[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md)

## Examples

``` r
m <- aci_dyad_model()
sim <- simulate(m, seed = 1, t_end = 2, dt = 0.01)
ob <- as_obs(sim)
f <- aci_filter(m, ob)
#> Warning: No init$cov supplied; using a diffuse prior. Discard an initial burn-in window when interpreting results.
s <- aci_smoother(m, ob, filter = f)
head(aci_metric(s, f))
#>      t     total    signal dispersion
#> 1 0.00 1.7548143 0.5349939  1.2198205
#> 2 0.01 2.4916834 1.9988182  0.4928653
#> 3 0.02 2.0958385 1.6500602  0.4457784
#> 4 0.03 2.5225164 2.1165660  0.4059504
#> 5 0.04 1.2681616 0.8911746  0.3769870
#> 6 0.05 0.9794228 0.6355366  0.3438862
```
