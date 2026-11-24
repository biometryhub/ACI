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
