# Reduce a model by prescribing the conditioning channels

Rewrites a model so that the conditioning observed channels enter as
known forcing rather than as observations to be assimilated, leaving the
target channels as the only observed process. x_B becomes prescribed
forcing in the (x_A, y) system.

## Usage

``` r
aci_conditional_reduce(model, obs, spec)
```

## Arguments

- model:

  A `cgns_model` object.

- obs:

  An observed trajectory, or anything
  [`as_obs()`](https://biometryhub.github.io/ACI/reference/as_obs.md)
  accepts.

- spec:

  A `aci_conditional_spec` object; see
  [`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md).

## Value

A list with the reduced `model` and the reduced `obs`.

## See also

[`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md)

## Examples

``` r
m2 <- aci_model(
  Lx = function(t, x) matrix(c(1, 0), 2, 1),
  fx = function(t, x) -0.5 * x,
  Ly = function(t, x) matrix(-0.5, 1, 1),
  fy = function(t, x) 0,
  Sx1 = function(t, x) diag(0.5, 2),
  Sy2 = function(t, x) matrix(1, 1, 1),
  k = 2, l = 1)
sim2 <- simulate(m2, seed = 1, t_end = 1, dt = 0.01)
aci_conditional_reduce(m2, as_obs(sim2), aci_conditional(2, "reduce"))
#> $model
#> <cgns_model> 'cgns_model|reduced': k = 1 observed, l = 1 hidden
#> 
#> $obs
#> <obs_traj> k = 1, N+1 = 101, dt = 0.01, span [0, 1]
#> 
#> $map
#> $map$A
#> [1] 1
#> 
#> $map$B
#> [1] 2
#> 
#> 
```
