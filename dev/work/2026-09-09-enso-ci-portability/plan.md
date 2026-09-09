# ENSO regularization test portability

## Problem and baseline

The five supplied Windows/Linux check logs fail the same four assertions in
`test-38-regularization-policy.R`. They compare downstream values from the
deliberately unstable, floored ENSO example with macOS snapshots at tolerance
`1e-9`. The warning, event-count and valid-example checks pass. The baseline
is clean commit `790fa9210ec49e7f0eba4f065b742d194f668447` on
`provision-observation-scale-guard`; its focused tests pass locally.

The first floor produces a covariance with condition number about `1.69e11`.
Perturbing one entry by `2.22e-16` changes the downstream smoother eigenvalues
and ACI beyond the old tolerance. Changing R's matrix-product option alone
does not reproduce the runner differences. This establishes sensitivity to
rounding, without identifying the first differing operation on each runner.

## Decision and scope

Correct the reporting test's use of these unstable outputs as portable
numerical references. Do not change package calculations, priors, covariance
floors, defaults, dependencies, workflows, fixtures or registered tolerances.
Changing those calculations to recover a platform-specific snapshot would
change the method for a case whose numerical accuracy is not established.
Increasing the tolerance or adding platform-specific expected values would
retain an unjustified reference.

Only `acir/tests/testthat/test-38-regularization-policy.R` and this plan change.
There is no user-visible package behavior change and no release-note change.
The test change is an evidence-review finding under `REVIEW.md`; it is not a
new independent numerical grade for the unstable ENSO result.

## Implementation

1. Retain every warning, print, site, count, location and stable first-step
   numerical assertion. Add online site counts and locations.
2. Compare each downstream recorded eigenvalue with the strict-policy error
   obtained by replaying that route from the same already-floored filter.
   There is one event at each site, so the first refusal is also the recorded
   worst value. Require finite negative values and the correct error class,
   site, index and time. The strict path stops before the recorder update.
3. Check that the stress result remains finite and grossly inflated (peak
   above `1e20`, an order-of-magnitude diagnostic bound, not an accuracy
   tolerance). Compare all times and metric values returned by `aci_metric()`
   with the saved ACI result, replacing its duplicate peak snapshot.
4. Add a deterministic matrix-guard check with known eigenvalues and repeated
   events. It independently checks recording of minima and first locations;
   replay consistency alone cannot establish the correctness of shared code.

## Verification and limits

Record the untouched source hashes, installation, environment, five supplied
logs and numerical probes under
`outputs/implementation/20260909-enso-ci-portability/` in the review workspace.
Verify the replacement checks under small controlled perturbations of the
first floored covariance, applied only in disposable diagnostic processes.
Use deliberate recorder/metric defects to check that the new assertions fail.

Run the focused tests, the full suite with the retained parity fixtures,
fixture provenance, a fresh build and installation, and `R CMD check --as-cran`.
Compare source hashes and installed baseline/candidate numerical outputs.
Report unavailable checks and any existing advisory benchmark result. Actual
Windows/Linux candidate verification requires a new CI run; local sensitivity
experiments do not replace it. No staging, committing or pushing is authorized
for this work session, so the plan and implementation remain an uncommitted
diff for review.
