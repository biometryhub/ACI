# One row of a lag table

Accesses the divergences at one anchor time across increasing positive
lag, padding the cells that adaptive storage did not retain.

## Usage

``` r
lt_row(x, j, pad = c("zero", "na"))
```

## Arguments

- x:

  A `lag_table` object.

- j:

  Integer index of the anchor time.

- pad:

  Either `"zero"` or `"na"`, the value used for cells the table did not
  retain.

## Value

Numeric vector of divergences at increasing positive lag.

## See also

[`lt_diag()`](https://biometryhub.github.io/ACI/reference/lt_diag.md),
[`aci_range()`](https://biometryhub.github.io/ACI/reference/aci_range.md)

## Examples

``` r
m <- aci_dyad_model()
sim <- simulate(m, seed = 1, t_end = 2, dt = 0.01)
ob <- as_obs(sim)
tb <- lag_table(m, ob, mode = "forward")
#> Warning: No init$cov supplied; using a diffuse prior. Discard an initial burn-in window when interpreting results.
head(lt_row(tb, 1))
#> [1] 1.7308809 2.3055094 1.8721223 2.1967183 1.0153139 0.7241675
```
