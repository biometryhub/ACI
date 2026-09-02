# Fixed-lag online data assimilation

Generic reconstructing the hidden state at each time from the observed
record up to a fixed number of steps ahead of it. The closed-form method
is used for a `cgns_model`; a general `stochastic_model` is out of scope
in this release.

## Usage

``` r
aci_online(model, obs, lag, ...)

# S3 method for class 'cgns_model'
aci_online(
  model,
  obs,
  lag,
  filter = NULL,
  conditional = NULL,
  init = NULL,
  regularize = NULL,
  force_validate = FALSE,
  ...
)

# S3 method for class 'stochastic_model'
aci_online(model, obs, lag, ...)
```

## Arguments

- model:

  A `cgns_model` or `stochastic_model` object.

- obs:

  An observed trajectory, or anything
  [`as_obs()`](https://biometryhub.github.io/ACI/reference/as_obs.md)
  accepts.

- lag:

  Number of future steps each estimate may condition on. A non-negative
  whole number, or `Inf` for the whole record. No default: the lag is
  the argument the function exists for, and defaulting it invites the
  full-lag result to be mistaken for
  [`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md).

- ...:

  Arguments passed to methods.

- filter:

  Optional precomputed filter path; recomputed when `NULL`. It must be
  the explicit single-step filter, which is the discretization the
  Theorem 3 recursions are exact for.

- conditional:

  Optional `aci_conditional_spec` selecting a conditional ACI reduction;
  see
  [`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md).

- init:

  Optional list with the initial hidden `mean` and `cov`.

- regularize:

  Covariance policy for this call; see
  [`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md).
  One record covers the whole call, and is returned in
  `meta$regularization`.

- force_validate:

  `FALSE` (the default) lets a `filter` that
  [`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md)
  produced for this same run, unaltered since, skip per-step
  re-validation, as in
  [`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md).

## Value

An assimilation path of kind `"online"`, carrying `meta$lag`, the
per-anchor `meta$lag_effective`, `meta$saturated` and `meta$scheme`. Its
kind is what keeps it out of the places a complete smoother is required:
`lag_table(smoother = )` rejects it with `aci_error_dims`.

## Details

`lag` is the number of future observations each estimate may condition
on: the estimate at index `j` uses the observed record through index
`j + lag`, and saturates at the end of the record. `lag = 0` returns the
filter moments unchanged. `lag = Inf` returns the complete Theorem 3
posterior given the whole record.

## Methods (by class)

- `aci_online(cgns_model)`: Closed-form fixed-lag online Theorem 3
  smoother for a conditional-Gaussian model.

- `aci_online(stochastic_model)`: Classed not-implemented condition for
  a general (non-CGNS) stochastic model.

## Scheme

`aci_online()` computes the **discrete** Theorem 3 posterior: the exact
conditional law of the hidden state given the observed increments on the
sampling grid under the explicit single-step discretization.
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md)
integrates the **continuous** backward smoothing equations with an Euler
step of the same size. These are two discretizations of the same
continuous-time object and they agree only to first order in the step,
so at full lag `aci_online()` does **not** reproduce
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md).
The gap grows with the length of the record; it is not a constant
offset. On the packaged ENSO partition (`l = 3`, `dt = 0.005`) the
smoothed means differ by up to 1.89e-02 against a mean scale of 0.388
over 401 steps, and by up to 9.58e-02 against a scale of 2.22 over 4001
steps. The resulting ACI values differ by up to 0.104 against a scale of
1.093 at 401 steps and 0.482 against 2.347 at 4001.
[`lag_table()`](https://biometryhub.github.io/ACI/reference/lag_table.md)
and the causal influence range estimators use the discrete scheme
throughout, which is why
[`lt_diag()`](https://biometryhub.github.io/ACI/reference/lt_diag.md)
and `aci()$aci` can differ by that same amount. `meta$scheme` records
which scheme produced a path: `"theorem3_discrete"` here and on the lag
table's reference smoother, `"backward_ode_euler"` on
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md).

The one exact boundary is the other end: at `lag = 0` the returned
moments are the filter moments, unchanged value for value.

## References

Andreou, M., Chen, N. and Li, Y. (2026). An adaptive online smoother
with closed-form solutions and information-theoretic lag selection for
conditional Gaussian nonlinear systems. *Journal of Nonlinear Science*
**36**(4), 71.
[doi:10.1007/s00332-026-10271-x](https://doi.org/10.1007/s00332-026-10271-x)

## See also

[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md),
[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md),
[`lag_table()`](https://biometryhub.github.io/ACI/reference/lag_table.md)

## Examples

``` r
m <- aci_dyad_model()
sim <- simulate(m, seed = 1, t_end = 2, dt = 0.01)
ob <- as_obs(sim)
aci_online(m, ob, lag = 5)
#> Warning: No init$cov supplied; using a diffuse prior. Discard an initial burn-in window when interpreting results.
#> <da_path_gaussian> kind = online, l = 1, N+1 = 201
```
