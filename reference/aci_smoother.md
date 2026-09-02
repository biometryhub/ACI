# Data assimilation smoother

Generic reconstructing the hidden state from the whole observed record.
The closed-form method is used for a `cgns_model`; a general
`stochastic_model` is out of scope in this release.

## Usage

``` r
aci_smoother(model, obs, ...)

# S3 method for class 'cgns_model'
aci_smoother(
  model,
  obs,
  filter = NULL,
  conditional = NULL,
  init = NULL,
  stepper = c("explicit", "implicit"),
  nsub = 1L,
  regularize = NULL,
  force_validate = FALSE,
  ...
)

# S3 method for class 'stochastic_model'
aci_smoother(model, obs, ...)
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

- filter:

  Optional precomputed filter path; recomputed when `NULL`.

- conditional:

  Optional `aci_conditional_spec` selecting a conditional ACI reduction;
  see
  [`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md).

- init:

  Optional list with the initial hidden `mean` and `cov`.

- stepper:

  Either `"explicit"` or `"implicit"`.

- nsub:

  Positive whole number of sub-steps taken per observation.

- regularize:

  Covariance policy for this call; see
  [`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md).
  One record covers the whole call, so a filter recomputed here and the
  backward recursion that consumes it share the `meta$regularization` on
  the returned smoother.

- force_validate:

  `FALSE` (the default) lets a `filter` that
  [`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md)
  produced for this same model, observations and conditional
  specification, and that has not been altered since, skip the per-step
  re-validation of its covariances. Any other supplied path, including
  one that has been through
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html), is validated in
  full as before. `TRUE` validates unconditionally. The smoother result
  is the same either way.

## Value

An assimilation path: `da_path_gaussian` for the closed-form engine.

## Methods (by class)

- `aci_smoother(cgns_model)`: Closed-form backward-ODE smoother for a
  conditional-Gaussian model.

- `aci_smoother(stochastic_model)`: Classed not-implemented condition
  for a general (non-CGNS) stochastic model.

## See also

[`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md),
[`aci()`](https://biometryhub.github.io/ACI/reference/aci.md)

## Examples

``` r
m <- aci_dyad_model()
sim <- simulate(m, seed = 1, t_end = 2, dt = 0.01)
ob <- as_obs(sim)
f <- aci_filter(m, ob)
#> Warning: No init$cov supplied; using a diffuse prior. Discard an initial burn-in window when interpreting results.
aci_smoother(m, ob, filter = f)
#> <da_path_gaussian> kind = smoother, l = 1, N+1 = 201
```
