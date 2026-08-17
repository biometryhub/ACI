# Simulate a conditional Gaussian nonlinear system

Simulates a realisation of the observed and unobserved processes of an
`aci_model` by an Euler-Maruyama integration of its pair of stochastic
differential equations. The simulation starts from the model's initial
state and draws independent standard normal increments at each step.

## Usage

``` r
aci_simulate(
  model,
  n,
  dt = 0.001,
  seed = NULL,
  scheme = c("euler_maruyama", "milstein"),
  sigma_x = NULL,
  d_sigma_x = NULL,
  increments = NULL
)
```

## Arguments

- model:

  An `aci_model` object; see
  [`aci_cgns_model()`](https://biometryhub.github.io/ACI/reference/aci_cgns_model.md).

- n:

  Integer scalar. The number of time steps to simulate, including the
  initial step; at least two.

- dt:

  Numeric scalar. The integration time step; must be positive.

- seed:

  Optional integer. If supplied, seeds a reproducible path and restores
  the caller's generator state on exit.

- scheme:

  Character scalar. The integration scheme, `"euler_maruyama"` (the
  default) or `"milstein"`. The two coincide unless the diffusion varies
  with the state it integrates, so `"milstein"` requires `sigma_x` and
  `d_sigma_x`.

- sigma_x:

  Optional function of the observed state. Supplying it makes the
  observation noise multiplicative, overriding the model's constant
  `S_xoS_x` **for the simulation only**. This is a simulation
  capability. The filter still requires a constant observation-noise
  covariance, so a path generated this way cannot yet be assimilated by
  this package.

- d_sigma_x:

  Optional function of the observed state, the derivative of `sigma_x`.
  Required by the Milstein scheme, and never estimated numerically on
  the caller's behalf.

- increments:

  Optional named list of the standard normal variates to integrate, with
  elements `dW_x` and `dW_y`, each a complete finite numeric vector of
  length `n - 1`. Supplying it drives the integrator with those variates
  instead of drawing any, which is how a path captured from another
  implementation is reproduced here. The Wiener increment applied at a
  step is `sqrt(dt)` times the variate, following the reference
  implementation's own naming. Cannot be combined with `seed`.

## Value

A data frame with `n` rows and columns `t` (time), `x` (the observed
process) and `y` (the unobserved process).

## Details

The random increments come from R's generator, which differs from that
of the reference MATLAB implementation, so a simulated path here does
not reproduce the reference path even at a matching seed. The
independent-oracle grade of the numerical core is therefore run on the
MATLAB signal itself, not on a fresh simulation.

Supplying `seed` makes the path reproducible without disturbing the
calling session. The generator state is saved before the draw and
restored when the function exits, so a seeded call has no effect on the
sequence a caller would otherwise have seen. Leaving `seed` as `NULL`
consumes the global stream in the ordinary way.

Supplied increments provide what a matching seed cannot. `increments`
takes the standard normal variates the integrator would otherwise draw,
one per transition in each of `dW_x` and `dW_y`, and no random numbers
are generated at all. This is what makes the integration scheme
gradeable against another implementation. Euler-Maruyama is invertible,
so the variates a reference run used are recoverable from its own
captured path by subtracting the drift and dividing by the noise
coefficient; feeding them back here reproduces that path itself, rather
than merely a path with the same statistics. Supplying `increments`
together with `seed` is an error, because the seed would govern nothing.

Simulation supports independent noise, that is a model with zero noise
cross-covariance.

## References

Andreou, M., Chen, N. and Bollt, E. (2026). Assimilative causal
inference. *Nature Communications*, 17, 1854.
[doi:10.1038/s41467-026-68568-0](https://doi.org/10.1038/s41467-026-68568-0)

## See also

[`aci_dyad_model()`](https://biometryhub.github.io/ACI/reference/aci_dyad_model.md),
[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md)

## Examples

``` r
model <- aci_dyad_model()
sim <- aci_simulate(model, n = 5000, seed = 333)
head(sim)
#>       t         x        y
#> 1 0.000 1.0000000 2.000000
#> 2 0.001 1.0026906 2.017046
#> 3 0.002 1.0373242 2.002320
#> 4 0.003 1.0090259 2.033431
#> 5 0.004 1.0175164 2.056711
#> 6 0.005 0.9975656 2.059818

# A seeded call leaves the caller's random stream untouched.
set.seed(1)
before <- runif(1)
set.seed(1)
invisible(aci_simulate(model, n = 10, seed = 99))
identical(before, runif(1))
#> [1] TRUE

# Supplied increments make the path a pure function of its input, which is
# how a run captured elsewhere is reproduced here rather than only matched
# in distribution.
set.seed(7)
z <- list(dW_x = rnorm(999), dW_y = rnorm(999))
driven <- aci_simulate(model, n = 1000, increments = z)
identical(driven$x, aci_simulate(model, n = 1000, increments = z)$x)
#> [1] TRUE
```
