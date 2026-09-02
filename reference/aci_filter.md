# Data assimilation filter

Generic reconstructing the hidden state from the observed record up to
each time. The closed-form method is used for a `cgns_model`; a general
`stochastic_model` is out of scope in this release.

## Usage

``` r
aci_filter(model, obs, ...)

# S3 method for class 'cgns_model'
aci_filter(
  model,
  obs,
  init = NULL,
  conditional = NULL,
  stepper = c("explicit", "implicit"),
  nsub = 1L,
  regularize = NULL,
  loglik = TRUE,
  ...
)

# S3 method for class 'stochastic_model'
aci_filter(model, obs, ...)
```

## Arguments

- model:

  A `cgns_model` or `stochastic_model` object.

- obs:

  An observed trajectory, or anything
  [`as_obs()`](https://biometryhub.github.io/ACI/reference/as_obs.md)
  accepts.

- ...:

  Arguments passed to methods.

- init:

  Optional list with the initial hidden `mean` and `cov`; `NULL` uses a
  diffuse prior and warns.

- conditional:

  Optional `aci_conditional_spec` selecting a conditional ACI reduction;
  see
  [`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md).

- stepper:

  Either `"explicit"` or `"implicit"`; the implicit Riccati step
  preserves positivity.

- nsub:

  Positive whole number of sub-steps taken per observation.

- regularize:

  Covariance policy for this call. `"none"` (the default, and the value
  of `getOption("aci.regularize")` when it is unset) stops with a
  classed `aci_error_covariance_not_spd` naming the site, grid index and
  time as soon as a covariance leaves the positive-definite cone.
  `"floor"` is the previous behaviour: the covariance is projected back
  by
  [`spd_floor()`](https://biometryhub.github.io/ACI/reference/spd_floor.md)
  and every such event is recorded in the result's
  `meta$regularization`.

- loglik:

  `TRUE` (the default) accumulates the predictive log-likelihood into
  `meta$loglik`. `FALSE` skips that work; the filter moments are
  unchanged and `meta$loglik` is `NULL`. The likelihood is not used by
  ACI, so `FALSE` is the cheaper choice when only the state estimate is
  wanted.

## Value

An assimilation path: `da_path_gaussian` for the closed-form engine. Its
`meta$loglik` holds the predictive log-likelihood of the observed
record, or `NULL` when the method was called with `loglik = FALSE`.

## Methods (by class)

- `aci_filter(cgns_model)`: Closed-form filter for a
  conditional-Gaussian model.

- `aci_filter(stochastic_model)`: Classed not-implemented condition for
  a general (non-CGNS) stochastic model.

## See also

[`aci_smoother()`](https://biometryhub.github.io/ACI/reference/aci_smoother.md),
[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md),
[`lag_table()`](https://biometryhub.github.io/ACI/reference/lag_table.md)

## Examples

``` r
m <- aci_dyad_model()
sim <- simulate(m, seed = 1, t_end = 2, dt = 0.01)
ob <- as_obs(sim)
f <- aci_filter(m, ob)
#> Warning: No init$cov supplied; using a diffuse prior. Discard an initial burn-in window when interpreting results.
f
#> <da_path_gaussian> kind = filter, l = 1, N+1 = 201
```
