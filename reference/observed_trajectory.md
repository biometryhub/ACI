# Observed trajectory on a uniform time grid

This release assimilates a complete record of the observed state on a
uniform grid, taken as noise-free. Every observed channel must be
present at every time, the times must be strictly increasing and
uniformly spaced, and the values must be finite; independent sensor
error is refused rather than absorbed. `observed_trajectory()`
constructs the observation object consumed throughout the package and
enforces those requirements.

## Usage

``` r
observed_trajectory(t, x, noise_free = TRUE, names = NULL)
```

## Arguments

- t:

  Numeric vector of observation times, strictly increasing and uniformly
  spaced.

- x:

  Numeric matrix of observations with one row per time and at least one
  column, or a vector coerced to a single column.

- noise_free:

  Logical; must be `TRUE`, since noisy observations are not supported in
  this version.

- names:

  Optional character vector of unique, non-empty column names, one per
  observed channel.

## Value

An object of class `obs_traj`: a list with the time vector `t`, the
observation matrix `x`, the step `dt`, the observed dimension `k` and
the flag `noise_free`.

## Details

The noise-free restriction follows the method's current published scope:
andreou2026cir Section 2.1 leaves noise-contaminated observations to
future work, and its Discussion lists them as an open direction, so this
is a limitation of the framework as published, not a modelling
assumption added by the package.

## See also

[`as_obs()`](https://biometryhub.github.io/ACI/reference/as_obs.md),
[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md)

## Examples

``` r
observed_trajectory(t = seq(0, 1, by = 0.1),
                    x = matrix(rnorm(11), ncol = 1))
#> <obs_traj> k = 1, N+1 = 11, dt = 0.1, span [0, 1]
```
