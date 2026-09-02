# Coerce to an observed trajectory

Generic converting the supported observation representations to an
`obs_traj` object.

## Usage

``` r
as_obs(x, ...)

# S3 method for class 'obs_traj'
as_obs(x, ...)

# S3 method for class 'matrix'
as_obs(x, dt = 1, t0 = 0, ...)

# S3 method for class 'aci_sim'
as_obs(x, ...)
```

## Arguments

- x:

  Object to coerce: an `obs_traj`, a numeric matrix, or a simulation of
  class `aci_sim`.

- ...:

  Arguments passed to methods.

- dt:

  Positive 1-length numeric step used to build the time grid.

- t0:

  1-length numeric time of the first observation.

## Value

An object of class `obs_traj`; see
[`observed_trajectory()`](https://biometryhub.github.io/ACI/reference/observed_trajectory.md).

## Methods (by class)

- `as_obs(obs_traj)`: Returns the trajectory unchanged.

- `as_obs(matrix)`: Builds a uniform grid from `dt` and `t0` for a
  matrix of observations.

- `as_obs(aci_sim)`: Extracts the observation component of a simulation.

## See also

[`observed_trajectory()`](https://biometryhub.github.io/ACI/reference/observed_trajectory.md)

## Examples

``` r
m <- aci_dyad_model()
sim <- simulate(m, seed = 1, t_end = 2, dt = 0.01)
ob <- as_obs(sim)
as_obs(sim)
#> <obs_traj> k = 1, N+1 = 201, dt = 0.01, span [0, 2]
```
