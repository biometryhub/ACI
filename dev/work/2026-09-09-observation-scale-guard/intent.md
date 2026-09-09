# Intent: refuse or record an observation model that is degenerate in scale, not only in conditioning

Author: Max Moldovan. Status: draft.

Follows `2026-09-08-initial-release-corrections`, which added the realised-Gram
check. That check closed the exactly-zero case. This work package closes the
class the check cannot see, and is deliberately held out of 0.1.0: for the
release the limitation is documented rather than removed.

## Problem

`.check_gram_path()` tests `rcond(gxx)` against `.ACI_GRAM_RCOND_MIN`. `rcond()`
is invariant under a uniform rescaling of the Gram, which is what makes it the
right test for conditioning and the wrong one for scale. For a one-dimensional
observed state the reciprocal condition number of every non-zero Gram is exactly
1, so only a Gram that underflows to zero is refused.

The two steppers then diverge, because the explicit kernel's only covariance
tests are `!is.finite(R)` and `R <= 0` (`R/aci-kernels-scalar.R:365-368`). A
near-zero Gram makes the explicit Riccati step overshoot to negative, so it
refuses. The implicit stepper preserves positivity by construction, so nothing
fires and the filter returns a collapsed posterior.

Measured on a scalar model whose `Sx1` falls at one interval start of a
`dt = 0.01` grid, prior variance 1, all other intervals at `Sx1 = 0.5`:

| `Sx1` at the interval | explicit | implicit |
|---|---|---|
| `1e-4` | refuses, `aci_error_covariance_not_spd` | completes, mean 0.999999, sd 1e-3 |
| `1e-8` | refuses | completes, mean 1, sd 1e-7, finite loglik |
| `1e-30` | refuses | completes, mean 1, sd 1e-29 |
| `1e-160` | refuses, `aci_error_gram_path` | refuses (the Gram underflows to exactly 0) |

The implicit rows are the defect: a posterior pinned to the noiseless limit with
a variance collapsed by fourteen orders of magnitude or more in a single step,
no condition raised, no event recorded in `meta`.

## Proposed outcome

A test that measures scale rather than conditioning, at one of two sites.

**Outcome side** (preferred): give the implicit route the collapse test the
explicit route gets for free, by refusing or recording when the posterior
variance falls by more than a stated ratio in one step with no regularisation
policy in force.

**Input side**: compare the realised `gxx` against the state uncertainty
propagated through the coupling, roughly `Lx R Lx' h`, and refuse where the gain
has saturated at its noiseless limit.

The separation is wide enough that the threshold is not a close call. Worst
single-step variance drop on the graded records:

| Record | Worst one-step variance drop |
|---|---|
| dyad, the authors' graded record | 1.01 |
| ENSO documented configuration, `y1` / `y2` / `y3` | 1.001 / 1.004 / 1.277 |
| the `1e-8` counterexample above | about 1e14 |
| the `1e-30` counterexample above | about 1e58 |

Any threshold between 1e2 and 1e12 separates the two populations, with nothing
of ours near the boundary.

## Affected code and users

`R/aci-kernels-scalar.R` and `R/aci-kernels-matrix.R` at the covariance guards;
`R/aci-conditional.R` if the test goes on the input side beside
`.check_gram_path()`; the `regularize` documentation in `R/aci-assimilation.R`,
whose current text records this limitation and would be rewritten; `NEWS.md`;
`tests/testthat/test-35-observation-gram.R`, which covers both steppers for the
exactly-zero case and would gain the near-zero rows.

Users who notice: anyone running the implicit stepper on a model whose
observation noise approaches zero somewhere in the record. A model that is
accepted today and refused afterwards is a breaking change and needs its NEWS
entry.

## Constraints

**The guard may refuse or record. It may never alter a value on an input it
accepts.** The parity claim is a claim about numbers on the graded surface, so a
refusal costs nothing and a repair -- flooring the Gram, adding jitter, changing
the solver or the gain -- would move numbers on admissible inputs and end it.
Nothing here touches `safe_chol()`, whose jitter is a separate published
contract.

The acceptance evidence is the existing suite: 7,135 assertions including the 28
fixture-backed grades. Green after the change means every graded comparison
still ran and still met its recorded tolerance.

The threshold must be invariant under a rescaling of the state and under a
change of time unit, which `rcond()` supplied for nothing and a scale test does
not.

## Open questions

1. Outcome side or input side, and whether the matrix kernel needs the same
   test: `rcond` sees conditioning there but still not scale.
2. Error or warning. A very precise sensor is a legitimate model, not a contract
   violation, which argues for a recorded warning rather than a refusal.
3. Whether the margin sweep should cover every graded record before a threshold
   is fixed. The dyad and ENSO records are measured above; the predator-prey and
   partition records are not.

## What this work package does not do

It does not repair a degenerate observation model, introduce a pseudoinverse or
a new degenerate-observation interpretation, or change `safe_chol()`. It does
not revisit the exactly-zero Gram check, which stands.
