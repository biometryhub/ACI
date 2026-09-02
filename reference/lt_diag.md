# Diagonal of a lag table

Accesses the zero-lag entries of a table without depending on its
storage representation.

## Usage

``` r
lt_diag(x)
```

## Arguments

- x:

  A `lag_table` object.

## Value

Numeric vector of the zero-lag divergences, one per time.

## See also

[`lt_row()`](https://biometryhub.github.io/ACI/reference/lt_row.md)

## Examples

``` r
m <- aci_dyad_model()
sim <- simulate(m, seed = 1, t_end = 2, dt = 0.01)
ob <- as_obs(sim)
tb <- lag_table(m, ob, mode = "forward")
#> Warning: No init$cov supplied; using a diffuse prior. Discard an initial burn-in window when interpreting results.
head(lt_diag(tb))
#> [1] 1.7308809 2.3989221 2.0003218 2.3959394 1.1828418 0.9013055
```
