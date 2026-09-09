# Heuristic tail estimate of a lag table

The historical `lt_tail_bound()` name is retained for compatibility, but
its value is a heuristic tail estimate, not a certified mathematical
error bound. What the package computes is a norm-based suffix
accumulation with a multiplier of 1.5; andreou2026smoother eq. 3.19 is a
spectral-radius condition, and individual spectral radii do not
establish contraction of an arbitrary ordered product. The value is
therefore a diagnostic on the retained record, not a guarantee about the
cells the truncation dropped.

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
#> Warning: No init$cov supplied; using a diffuse prior. Its opening steps are prior-dominated; a prior far wider than the hidden state's own scale can also destabilise the explicit step, which is a refusal rather than a window to discard.
head(lt_tail_bound(tb))
#> [1] 0 0 0 0 0 0
```
