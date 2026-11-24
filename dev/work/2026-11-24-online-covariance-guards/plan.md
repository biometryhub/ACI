# Plan: covariance policy on the fixed-lag online smoother (from intent.md 2026-11-24)

Status: draft 2026-11-24.

## Files that change

- `acir/R/aci-utils.R`: one new entry in `.ACI_COV_SITES`,
  `smoother_window`, for the fixed-lag smoothed covariance.
- `acir/R/aci-online-smoother.R`: `.online_window_compiled()` and
  `.online_at_compiled()` take the recorder as a last argument with the same
  default as `.online_aux_compiled()`, and pass each anchor's covariance
  through `.cov_guard()` in place of `sym()`; `.da_online_compiled()` passes
  its recorder on. The internal help of `.online_window_compiled()` says so.
- `acir/R/aci-assimilation.R`: the `regularize` help of `aci_online()` says
  that each anchor's fixed-lag covariance is checked.
- `acir/man/aci_online.Rd`: regenerated.
- `acir/tests/testthat/test-48-online-covariance-guards.R`: new.
- `acir/NEWS.md`: one bullet under 0.2.2, Bug fixes.

## Order of work

1. Reproduce both faults on the base code.
2. Measure, for the tight-prior fault, the largest generalised eigenvalue of
   the smoothed covariance against the filter covariance, minus one, over
   every anchor, on the records and priors the test suite uses, the library
   models with one, two and three hidden components, the generic closure
   model, long records, strong coupling and both `dt` scales, and on a sweep
   of priors through one step of hidden noise. Add a guard only if a
   threshold sits two orders of magnitude from both sides.
3. Write the tests and confirm they fail on `ba3559b`.
4. Change the window and pair routes, the site register and the help.
5. Run the full suite and confirm that no existing test, fixture or
   tolerance moves.

## Risks

- A guard could change an accepted value. `.cov_guard()` returns
  `(R + t(R)) / 2` when the Cholesky factor exists, the same arithmetic as
  `sym()`, so accepted values are identical to the bit; the tests compare
  them with `identical()`.
- The guard adds one Cholesky factorisation per anchor to the finite-lag
  route. The matrices are `l` by `l`.
- Under `regularize = "floor"` the pair route floors each returned snapshot,
  not the running sum, so that it agrees with the window route.

## Proof

`test-48-online-covariance-guards.R`: the refusal under the default policy
with one and two hidden components, through the window route and the pair
route; the flooring and its record under `regularize = "floor"`; with and
without a supplied `filter =`; a library model and a closure model; and
bit-identical results on accepted inputs at lags next to the failing one.
The full suite at 0 failures.

## Departures

Empty at approval.

### 2026-11-24, tight prior

- **No guard for the tight-prior fault.** Step 2 found no threshold with two
  orders of magnitude on each side. On accepted inputs the excess of the
  smoothed covariance over the filter covariance reaches 2.2e-16 at
  `lag = Inf` and 4.7e-05 at finite lags (ENSO, three hidden components,
  prior 0.01). On the faulty side the excess rises continuously from zero as
  the prior variance falls below one step of hidden noise: 1.4e-08 on the
  linear model at 1.0051 steps, 3.4e-03 on the dyad at 1.01 steps, 826 on
  the reported dyad record. The fault is left unchanged.

### 2026-11-24, one hidden component

- **An inline test comes before the check.** With one hidden component the
  window route tests the variance for finite and positive, as the scalar
  kernels do, and calls the check only when that test fails. Calling the
  check at every anchor added about 30% to a fixed-lag run on a 4,001-point
  scalar record (0.050 s to 0.066 s). With the inline test the time is
  0.050 s again and every value is unchanged.

### 2026-11-25, tight prior: a warning keyed to the step ratio

This replaces "no guard" above. The fault has no usable boundary in the
output and a sharp one in the input.

- **Cause.** To leading order the Theorem 3 update factor is the identity
  minus `Q Rf^-1 dt`, with `Rf` the filter covariance at the start of the
  interval and `Q = gyy - gyx gxx^-1 gxy` the hidden noise covariance given
  the observed increment. The exact factor stays positive. The leading-order
  one changes sign in a direction once the largest eigenvalue of
  `Q Rf^-1 dt`, the step ratio, passes 1.
- **Measured, 168 cases whose filter runs.** The ratio is above 1 in 79, at
  the first interval only in every one. Of the 73 cases with a smoothed variance above the filter
  variance, 72 have a ratio above 1. The other is the three-component
  finite-lag case at 4.7e-05, with a ratio of 0.57, which this change does
  not cover. Seven cases above 1 show no excess; all have a ratio of 1.1 or
  less, and their first-anchor variance is 1.3% to 23% above the exact
  discrete value.
- **Correlated noise.** The boundary follows `Q`, not `gyy`. On the
  two-component closure model `gyy Rf^-1 dt` passes 1 at a prior of 0.0068
  and nothing is wrong up to 1.59; the excess starts between step ratios of
  0.997 and 1.031.
- **Change.** `.thmD1_aux_compiled()` and `.online_aux_scalar()`, the two
  places every Theorem 3 update is formed, compute the step ratio and call a
  new `.aci_step_ratio_warn()`, which raises one `aci_warn_leading_order`
  per public call through a latch on the recorder. With several hidden
  components the trace of `Q Rf^-1 dt` bounds the ratio from above, and the
  eigenvalues are taken only on an interval the trace does not clear. The
  scalar helpers take the recorder as a last argument. The threshold is 1,
  the point where the factor changes sign, and is not tuned.
- **A floored covariance is left out.** Under `regularize = "floor"` an
  interval whose filter covariance the policy has just floored does not
  enter the ratio: the ratio against a floor value says nothing about the
  model, and the flooring is reported by `aci_warn_regularized`.
- **Not changed.** Every returned value. The call does not stop: seven of
  the 79 cases return a variance the method allows. A lag table reused
  through `aci(table = )` forms no update and raises nothing.
- **Files added to the list above.** `acir/R/aci-online-scalar.R`,
  `acir/R/aci-cir.R` (one call site), `acir/man/lag_table.Rd`; a second
  NEWS bullet.
- **Proof added.** In `test-48-online-covariance-guards.R`: one warning per
  call on `aci_online()` at `lag = Inf` and at a finite lag, `lag_table()`
  and `aci_range()`, with one and several hidden components, a library
  model and a closure model, and a filter read back from a saved copy; the
  ratio and index the condition carries; no warning just below 1; and the
  correlated-noise pair on each side of the boundary, which a ratio built
  from `gyy` alone gets wrong.
