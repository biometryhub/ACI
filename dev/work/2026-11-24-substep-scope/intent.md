# Intent: state what `nsub` and `stepper` change, and test it

Author: the `acir` authors. Status: draft 2026-11-24.

## Problem

The help describes `nsub` only as "sub-steps taken per observation" and says
the implicit Riccati step "preserves positivity". Both are incomplete.

- `nsub` sub-steps the forward filter, and the filter then approaches the
  exact solution of its equations over each interval. The smoother's backward
  pass takes the same sub-steps but holds the filter mean and covariance at
  their values at the start of each interval, so it does not approach the
  exact smoother. On a scalar linear model at `dt = 0.04` the largest error
  in the smoothed variance is 0.0044 at `nsub = 1` and 0.0122 at
  `nsub = 512`. On the uncoupled Brownian null in `?aci`, whose exact ACI is
  0, the peak ACI is 1.466 at `nsub = 1`, 0.0141 at `nsub = 20` and 0.0121 at
  `nsub = 10000`.
- `stepper` applies to the forward filter only. The backward pass is
  explicit under either stepper. `aci(stepper = "implicit")` can stop in the
  backward pass with a message that advises `stepper = "implicit"`.
- No row of `inst/evidence/register.csv` covers `nsub > 1` or the implicit
  stepper, and no test compares either with an exact reference. Reversing the
  sign of the innovation in the sub-stepped scalar smoother leaves the whole
  suite passing.

## Proposed outcome

- One help section, "Sub-steps and steppers", on `?aci_filter`, says what is
  sub-stepped, what is held fixed, what a larger `nsub` does and does not
  change, and that `stepper` applies to the forward filter only. `?aci`,
  `?aci_smoother` and `?aci_online` point to it.
- The covariance error message is true whichever stepper the caller used.
- One new test file compares the sub-stepped filter (both steppers) and the
  sub-stepped smoother with closed-form references written in the test, on
  the scalar and the matrix kernels, and fails when the innovation sign of
  either sub-stepped smoother is reversed.
- The register gains rows for these checks, and the README states its true
  size.

## Affected code and users

`acir/R/aci-assimilation.R`, `acir/R/aci-core.R`, `acir/R/aci-utils.R` (one
message), `acir/R/aci-model-library.R` (example comments), their help pages,
`acir/inst/evidence/register.csv`, `acir/README.md`, `acir/NEWS.md`. Users of
`nsub > 1` or `stepper = "implicit"`.

## Constraints

No numerical code changes; no stored output, fixture or tolerance moves. The
condition class of the covariance error stays. No existing `@param` line
changes.

## Open questions

Whether the smoother should sub-step against interpolated filter moments is a
change of method and is not part of this work.
