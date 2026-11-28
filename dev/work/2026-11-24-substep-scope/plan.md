# Plan: sub-step and stepper scope (from intent.md 2026-11-24)

Status: draft 2026-11-24.

## Files that change

- `acir/R/aci-assimilation.R`: new section "Sub-steps and steppers" in the
  `aci_filter()` block before `@returns`; one pointer sentence each in the
  `aci_smoother()` and `aci_online()` blocks.
- `acir/R/aci-core.R`: one pointer sentence in the `aci()` description.
- `acir/R/aci-utils.R`: the remedy sentence of `.aci_stop_cov()` says the
  implicit stepper acts on the forward filter only. Class unchanged.
- `acir/R/aci-model-library.R`: two example comments in `aci_enso_model()`
  that call a sub-stepped run "refining" are reworded.
- `acir/man/`: `aci_filter.Rd`, `aci_smoother.Rd`, `aci_online.Rd`,
  `aci.Rd`, `aci_enso_model.Rd`, regenerated.
- `acir/tests/testthat/test-50-substep-scope.R`: new.
- `acir/inst/evidence/register.csv`: three rows appended.
- `acir/README.md`: the register size sentence.
- `acir/NEWS.md`: one Documentation entry and one Bug fixes entry under
  0.2.2.

## Order of work

1. Write the tests against closed-form references:
   a. the explicit sub-stepped filter against the exact solution of the
      scalar linear filter equations over each interval, error halving as
      `nsub` doubles, on a library model, a closure-built `aci_model()` and a
      two-component diagonal model (matrix kernel);
   b. the implicit filter against the same reference, same order, scalar and
      matrix;
   c. the sub-stepped smoother against the closed-form solution of its
      backward equations with filter moments held per interval, error halving
      as `nsub` doubles, scalar and matrix kernels, explicit and implicit
      filter input;
   d. on an uncoupled model, where the exact smoother is the exact filter and
      the exact ACI is 0, the smoother's gap and the peak ACI do not shrink
      as `nsub` grows and halve when `dt` halves;
   e. a backward-pass covariance failure under `stepper = "implicit"` keeps
      its class and does not advise `stepper = "implicit"`, scalar and matrix
      smoother.
2. Confirm (e) fails on `ba3559b` and (a) to (d) pass there.
3. Prove (c) fails with the innovation sign reversed in each sub-stepped
   smoother kernel, on a copy of the package.
4. Change the message, write the help, regenerate, append the register rows,
   set the README counts from the file, add NEWS.
5. Run the full suite, the base check and the prose guard.

## Risks

- An order test can pass by accident if the error is dominated by round-off;
  the `nsub` range is kept where the error is far above it.
- The register rows could claim more than the tests check; each row names
  the test file and states the order bound, not a tolerance.

## Proof

`tests/testthat/test-50-substep-scope.R`; the two mutant runs; the full
suite.

## Departures

Empty at approval.

### 2026-11-24, build

- The innovation sign reversed in the matrix smoother kernel is already
  caught by the existing suite (59 failures in eight files), because that
  kernel has no separate single-step branch. The new smoother test catches it
  too.
- The implicit filter is graded only by its order against the exact filter.
  A separate check against the stationary variance would add nothing the
  order test does not already show.
- The error message keeps "raise nsub": on the failing case in the new test,
  `nsub = 16` runs to completion.

### 2026-11-24, vignette

- `acir/vignettes/vignette-2-advanced.Rmd` also changes, in three passages
  that described sub-stepped runs: the uncoupled null, the ENSO sequence
  that settles near 1.813, and the list of remedies after a refusal. Each
  now says that more sub-steps do not bring the smoothed result closer to
  the exact posterior. No code chunk changes.

## Review amendment, 27 November 2026

The numerical examples and tests remain intact. They establish residual
frozen-filter error, not that increasing nsub can never improve accuracy.
The null example itself improves substantially before reaching a nonzero error.

Correct the Sub-steps and steppers section and smoother pointer in
aci-assimilation.R, the ENSO example comment in aci-model-library.R and the
matching aci_filter.Rd, aci_smoother.Rd and aci_enso_model.Rd. Clarify that the
filter refinement holds the observed rate and coefficient values fixed.
Remove general error-rate claims that the selected exact-reference cases do
not establish. Update NEWS consistently.

Apply only vignette review items V01-V04: retain the recorded null values,
explain the backward pass under Steppers, and remove the unsupported ENSO-limit
and categorical no-improvement insertions. No other vignette proposal belongs
in this amendment.

Validate exact prose/help correspondence, unchanged executable R/example/chunk
content and patch application at the PR head. Preserve all existing tests,
fixtures and expectations. Their recorded execution above is not repeated.
No package tests/checks, benchmarks or vignette execution are run at the author's
direction. CI on the amended head and the other author's review remain necessary
before merge.
