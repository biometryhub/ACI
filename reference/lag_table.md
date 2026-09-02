# Finite-lag divergence table

Builds the table of finite-lag divergences consumed by the causal
influence range estimators. A lag table uses the complete online Theorem
3 smoother as its reference. That reference costs O(N) time-point work;
table construction then costs work proportional to the retained lag
cells, with O(N^2) cells for a full table in the worst case.

## Usage

``` r
lag_table(
  model,
  obs,
  mode = c("forward", "full"),
  tol = getOption("aci.default_tol", 1e-08),
  window = 3L,
  max_lag = Inf,
  filter = NULL,
  smoother = NULL,
  conditional = NULL,
  init = NULL,
  stepper = "explicit",
  nsub = 1L,
  regularize = NULL,
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

- mode:

  Either `"forward"` or `"full"`, selecting which cells are retained.

- tol:

  Positive tolerance below which a row is frozen by the adaptive storage
  rule.

- window:

  Number of consecutive steps a row must stay below `tol` before it is
  frozen.

- max_lag:

  Maximum positive lag retained, or `Inf` for no cap.

- filter:

  Optional precomputed filter path.

- smoother:

  Optional precomputed smoother path.

- conditional:

  Optional `aci_conditional_spec`; see
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
  One record covers the filter, the Theorem 3 reference smoother and
  every relative-entropy denominator the table forms, and is returned in
  `meta$regularization`.

- ...:

  Must be empty; unused arguments are an error.

## Value

An object of class `lag_table`.

## References

Andreou, M. and Chen, N. (2026). Bridging prediction and attribution:
identifying forward and backward causal influence ranges using
assimilative causal inference. arXiv:2510.21889v2, 4 August 2026.
[doi:10.48550/arXiv.2510.21889](https://doi.org/10.48550/arXiv.2510.21889)

Andreou, M., Chen, N. and Li, Y. (2026). An adaptive online smoother
with closed-form solutions and information-theoretic lag selection for
conditional Gaussian nonlinear systems. *Journal of Nonlinear Science*
**36**(4), 71.
[doi:10.1007/s00332-026-10271-x](https://doi.org/10.1007/s00332-026-10271-x)

## See also

[`aci_range()`](https://biometryhub.github.io/ACI/reference/aci_range.md),
[`lt_row()`](https://biometryhub.github.io/ACI/reference/lt_row.md)

## Examples

``` r
m <- aci_dyad_model()
sim <- simulate(m, seed = 1, t_end = 2, dt = 0.01)
ob <- as_obs(sim)
lag_table(m, ob, mode = "forward")
#> Warning: No init$cov supplied; using a diffuse prior. Discard an initial burn-in window when interpreting results.
#> <lag_table> mode = forward, N+1 = 201, tol = 1e-08
#>   mean retained lag: 100.0 steps; max heuristic tail estimate: 0.00e+00
```
