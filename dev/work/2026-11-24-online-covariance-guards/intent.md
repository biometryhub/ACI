# Intent: apply the covariance policy to the fixed-lag online smoother

Author: the `acir` authors. Status: draft 2026-11-24.

## Problem

`aci_online()` at a finite lag can return a covariance that is not positive
definite, with no condition. On a 13-point record of the linear model
`aci_linear_model(lambda_x = 2, lambda_y = 1, sigma_x = 0.5,
sigma_y = 0.05)` at `dt = 0.09`, with prior variance 0.5,
`aci_online(lag = 8)` returns a variance of -0.00787 at the first anchor
under the default policy: no error, no warning and no event in
`meta$regularization`. The same happens with two hidden components. At
`lag = Inf`, and in `lag_table()`, the same inputs stop with
`aci_error_covariance_not_spd`. The help of `aci_filter()` says the policy
stops a call as soon as a covariance leaves the positive-definite cone; the
finite-lag route only symmetrises the covariance it accumulates.

A second fault was reported on the same routes: with a prior variance below
about one step of hidden noise, the complete Theorem 3 smoother returns a
smoothed variance above the filter variance at the first anchor (0.826
against 0.001 on the dyad at `dt = 0.01`), and `aci(table = )` reports a peak
ACI of 409.8 with no condition.

## Proposed outcome

At a finite lag, each anchor's covariance passes through the same policy as
every other route: under `regularize = "none"` the call stops with
`aci_error_covariance_not_spd` naming the anchor's index and time; under
`regularize = "floor"` the covariance is floored by `spd_floor()`, the event
is recorded in `meta$regularization` and one `aci_warn_regularized` is
raised.

For the tight-prior fault, a guard is added only if a statistic that the
method guarantees (a smoothed covariance cannot exceed the filter covariance
at the same time) separates the accepted inputs from the faulty ones by at
least two orders of magnitude on each side. The margin is measured before
any threshold is chosen.

## Affected code and users

`acir/R/aci-online-smoother.R` (`.online_window_compiled()`,
`.online_at_compiled()`, `.da_online_compiled()`), `acir/R/aci-utils.R` (the
site register), the help of `aci_online()`. Users of `aci_online()` at a
finite lag.

## Constraints

No value changes on any input the package accepts. No change to any fixture,
tolerance or existing test. The Theorem 3 scheme itself is not changed.

## Open questions

None.
