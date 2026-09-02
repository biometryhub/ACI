# Assimilative causal inference

Runs the filter and smoother for a model and observed record and scores
each time by the relative entropy of the smoother from the filter. A
normal `aci()` call uses the supplied-code backward-ODE headline
smoother, including its correlated-noise correction, independently of
`keep`.
[`lag_table()`](https://biometryhub.github.io/ACI/reference/lag_table.md)
and `aci(table = ...)` instead use the complete online Theorem 3
smoother; their finite-grid diagonal can therefore differ from headline
ACI.

## Usage

``` r
aci(
  model,
  obs,
  engine = c("auto", "cgns"),
  conditional = NULL,
  table = NULL,
  keep = c("paths", "table", "none"),
  decompose = TRUE,
  init = NULL,
  stepper = c("explicit", "implicit"),
  nsub = 1L,
  regularize = NULL,
  loglik = TRUE,
  ...
)
```

## Arguments

- model:

  A `cgns_model` object.

- obs:

  An observed trajectory, or anything
  [`as_obs()`](https://biometryhub.github.io/ACI/reference/as_obs.md)
  accepts.

- engine:

  One of `"auto"` or `"cgns"`; `"auto"` selects the closed-form engine
  for a conditional-Gaussian model.

- conditional:

  Optional `aci_conditional_spec`; see
  [`aci_conditional()`](https://biometryhub.github.io/ACI/reference/aci_conditional.md).

- table:

  Optional precomputed `lag_table`, whose online diagonal is used in
  place of the headline smoother.

- keep:

  One of `"paths"`, `"table"` or `"none"`, selecting which objects are
  retained on the result. It does not select a smoother.

- decompose:

  `TRUE` to retain the signal and dispersion parts.

- init:

  Optional list with the initial hidden `mean` and `cov`.

- stepper:

  Either `"explicit"` or `"implicit"`.

- nsub:

  Positive whole number of sub-steps taken per observation.

- regularize:

  Covariance policy for this call; see
  [`aci_filter()`](https://biometryhub.github.io/ACI/reference/aci_filter.md).
  One record covers the filter, the smoother and any table this call
  builds, and is returned in `meta$regularization`.

- loglik:

  `TRUE` (the default) accumulates the predictive log-likelihood on the
  internal filter, where `keep = "paths"` exposes it as
  `paths$filter$meta$loglik`. ACI itself never uses it, so `FALSE` skips
  that work and leaves `paths$filter$meta$loglik` `NULL`; every ACI
  quantity is unchanged.

- ...:

  Must be empty; unused arguments are an error.

## Value

An object of class `aci_result`.

## See also

[`aci_range()`](https://biometryhub.github.io/ACI/reference/aci_range.md),
[`lag_table()`](https://biometryhub.github.io/ACI/reference/lag_table.md),
[`aci_metric()`](https://biometryhub.github.io/ACI/reference/aci_metric.md)

## Examples

``` r
m <- aci_dyad_model()
sim <- simulate(m, seed = 1, t_end = 2, dt = 0.01)
ob <- as_obs(sim)
a <- aci(m, ob)
#> Warning: No init$cov supplied; using a diffuse prior. Discard an initial burn-in window when interpreting results.
a
#> <aci_result> engine = cgns | peak ACI = 2.523 at t = 0.03
```
