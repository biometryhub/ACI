# Heuristic tail estimate of a lag table

The historical `lt_tail_bound()` name is retained for compatibility, but
its value is a heuristic tail estimate, not a certified mathematical
error bound.

## Usage

``` r
lt_tail_bound(x, j = NULL)
```

## Arguments

- x:

  A `lag_table` object.

- j:

  Optional integer index of a single anchor time; `NULL` returns the
  estimate at every time.

## Value

Numeric vector of heuristic tail estimates.

## Examples

``` r
m <- aci_dyad_model()
sim <- simulate(m, seed = 1, t_end = 2, dt = 0.01)
ob <- as_obs(sim)
tb <- lag_table(m, ob, mode = "forward")
#> Warning: No init$cov supplied; using a diffuse prior. Discard an initial burn-in window when interpreting results.
head(lt_tail_bound(tb))
#> [1] 0 0 0 0 0 0
```
