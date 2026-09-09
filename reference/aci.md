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
  builds, and is returned in `meta$regularization`. Flooring changes the
  numerical covariance so that the recursion can continue; it
  establishes nothing about the accuracy of the reconstruction or of the
  resulting information score, and a large finite ACI obtained after a
  floor is a diagnostic, not a result.

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

## Details

ACI is measured under the model, prior and observation record supplied
to it: it scores how much the later observed record sharpens the
hidden-state reconstruction those inputs imply, and does not test
whether they are correct. The per-time value is in nats; a record-length
total is in nats times model time. A positive value alone is not an
empirically identified causal effect, an intervention effect or a
significance statement, and it carries the discretisation. On a
structurally independent Brownian null (`dx = dW1`, `dy = dW2`,
uncoupled in both directions) with prior variance 1, one observed
interval of length `dt`, the observed path `x(t) = 2t`, the explicit
stepper at `nsub = 1` and `regularize = "none"`, the headline backward
smoother returns 1.0135e-4 nats at `dt = 0.1` and 2.4419e-8 at
`dt = 0.0125`. The same null returns much larger values at a step that
is coarse against the prior and process variance scales: at prior
variance 0.1 and `dt = 0.07` it returns 1.4660115 nats, which
`nsub = 20` reduces to 0.0140779. Every one of those runs is positive
definite and finite and records zero floor events, so positive-definite,
finite output with no floor event does not certify that the step
resolved the problem.

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
#> Warning: No init$cov supplied; using a diffuse prior. Its opening steps are prior-dominated; a prior far wider than the hidden state's own scale can also destabilise the explicit step, which is a refusal rather than a window to discard.
a
#> <aci_result> engine = cgns | peak ACI = 2.523 at t = 0.03
```
