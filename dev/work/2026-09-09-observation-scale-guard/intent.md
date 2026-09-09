# Intent: assess a diagnostic for observation-noise scale and time resolution

Author: Max Moldovan. Status: draft.

Follows `2026-09-08-initial-release-corrections`, which added the realised-Gram
check. That check closed the exactly-zero case. A further diagnostic is deferred
past 0.1.0. The release documents and tests the remaining limitation without
changing calculations or deciding that small positive noise is invalid.

## Problem

`.check_gram_path()` tests `rcond(gxx)` against `.ACI_GRAM_RCOND_MIN`. `rcond()`
is invariant under a uniform rescaling of the Gram, which is what makes it the
right test for conditioning and the wrong one for scale. For a one-dimensional
observed state the mathematical reciprocal condition number of every positive
Gram is 1. Small positive Grams can pass, although condition estimation can
also fail at extreme floating-point scales.

On an unresolved time step, small observation noise can cause the explicit
Riccati step to become non-positive and fail. The implicit stepper can instead
return a finite positive covariance with no warning or regularization event.
Neither successful execution nor a large variance drop establishes accuracy.
Weak coupling can also let the explicit route complete without a warning.

Reproducible cases are in `acir/tests/testthat/test-40-observation-scale.R`.
The near-zero case has `Lx = 1`, `Ly = -0.5`, `Sy2 = 0.5`, other drift/noise
terms zero, prior mean 0 and variance 1, and times and observed values
`(0, 0.01, 0.02, 0.03)`. `Sx1` is `1e-8` at time 0.01 and 0.5 elsewhere.
The explicit route refuses; the implicit route returns variance near `1e-14`
at time 0.02 without a warning/event. The corresponding matrix case is tested.

The same file supplies controls that a future diagnostic must distinguish:

- A static hidden state with prior variance 1 and constant observation noise
  `1e-8` has exact covariance `1/(1+t/1e-16)`. Its accurate implicit result
  drops by `1e14` in the first interval. A large drop alone cannot justify refusal.
- With `Sx1 = 1e-3`, `Ly = -0.5`, `Sy2 = 0.5`, and the prior at the exact
  stationary covariance, the implicit variance is about 81% too small at
  `dt = 0.01`, despite a largest one-step drop of only 5.16. Refinement reduces
  that error. A large-drop threshold is not an accuracy certificate.
- Consistent changes of observation, hidden-state and time units preserve the
  hidden posterior after conversion back. A raw absolute Gram cutoff would
  depend on the units.

The original `1e-160` explanation was incorrect: its squared Gram is still
representable as a positive subnormal value. In the local R/LAPACK check,
`rcond` returned zero at that scale; the Gram itself became zero at `1e-162`.
Floating-point cancellation that leaves a small nonzero coefficient is a
separate issue and cannot be diagnosed from its returned magnitude alone.

## Proposed outcome

A diagnostic that records a clearly defined numerical concern. No threshold,
warning class or new refusal policy is selected for this release.

**Outcome side**: measure covariance contraction, with a matrix statistic that
does not rely only on marginal variances. Distinguish a large update from an
inaccurate one. The explicit route does not provide such a diagnostic already.

**Input side**: assess a dimensionless comparison involving the observation
Gram, coupling, state covariance and step size. A highly informative observation
is not by itself a model-contract violation. Correlated noise and conditional
precision need to be covered by the chosen definition.

Historical margin measurements reported with the original proposal:

| Record | Worst one-step variance drop |
|---|---|
| dyad, the authors' graded record | 1.01 |
| ENSO documented configuration, `y1` / `y2` / `y3` | 1.001 / 1.004 / 1.277 |
| the original `1e-8` counterexample | about 1e14 |
| the original `1e-30` counterexample | about 1e58 |

The complete configurations and raw outputs for that sweep were not recorded
here. These values do not establish a safe threshold. The valid static-state
case also drops by `1e14`, and the unresolved stationary case drops by only
5.16. Sweep all retained graded records before proposing a threshold, while
keeping the independent valid and inaccurate controls above.

## Affected code and users

`R/aci-kernels-scalar.R` and `R/aci-kernels-matrix.R` at the covariance guards;
`R/aci-conditional.R` if the test goes on the input side beside
`.check_gram_path()`; the `regularize` documentation in `R/aci-assimilation.R`,
whose current text records the limitation; `NEWS.md`; and the regression files
`tests/testthat/test-35-observation-gram.R` (exactly-zero refusal) and
`tests/testthat/test-40-observation-scale.R` (small positive noise and controls).

Users affected by a future diagnostic depend on its definition and threshold.
A model accepted today and refused afterwards is a breaking change and needs
its NEWS entry.

## Constraints

**The guard may refuse or record. It may never alter a value on an input it
accepts.** Retaining accepted-path arithmetic preserves numerical comparisons
on those inputs; refusing a previously accepted model still changes the supported
input domain. Flooring the Gram, adding jitter, or changing the solver or gain
can change numerical outputs and requires separate review and comparisons.
Nothing here touches `safe_chol()`, whose jitter is a separate published
contract.

Keep every existing fixture and tolerance gate. Confirm that each graded
comparison ran rather than skipped, and report any new warnings/refusals.

The threshold must be invariant under a rescaling of the state and under a
change of time unit. Uniform Gram rescaling leaves the mathematical `rcond`
unchanged; arbitrary multivariate coordinate changes need separate treatment.

## Open questions

1. What does the diagnostic claim to detect: a large update, an unresolved step,
   or a model coefficient inconsistent with its intended scale?
2. What statistic, reference scale and threshold support that claim? A recorded
   advisory warning is preferable to refusal unless inadmissibility is established.
3. Does "step" mean an observation interval or an internal substep, and how does
   changing `nsub` affect the diagnostic?
4. How are scalar/matrix, generic/specialised, cached and conditional routes
   covered, including correlated noise and both covariance policies? How are
   events distinguished from flooring and propagated to downstream results?

## What this work package does not do

It does not repair a degenerate observation model, introduce a pseudoinverse or
a new degenerate-observation interpretation, or change `safe_chol()`. It does
not revisit the exactly-zero Gram check, which stands.
