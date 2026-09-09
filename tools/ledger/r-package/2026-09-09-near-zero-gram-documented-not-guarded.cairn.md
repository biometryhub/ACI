---
id: 2026-09-09-near-zero-gram-documented-not-guarded
schema_version: 1.4
date: 2026-09-09
tier: T3
classification: open
export_status: local_only
consensus_mode: none
domain: r-package
project: aciR
status: accepted
title: "The near-zero observation Gram ships documented for 0.1.0, and the scale guard is queued as its own work package"
tags: [release, observation_model, guards, documentation, parity]
triggers:
  - adr_class_commitment
  - audit_finding_downstream
  - interface_contract
reversal_cost: low
decision_pressure: publication
review_due: null
review_trigger: "when the observation-scale guard work package is implemented, or if a user reports a collapsed posterior on the implicit stepper"
supersedes: []
superseded_by: null
related: [2026-09-02-release-0.1.0-parity-milestone]
---

## Status

Decided 2026-09-09, in the session that verified the 8 September initial-release
action report against the package and reviewed the corrections branch
(`c6481ce`). The report's ADV-01 finding was confirmed in substance and found to
be mischaracterised in route: the realised-Gram check on that branch closes the
exactly-zero case, and a verification sweep found the near-zero case alive on
the implicit stepper. The decision covers the 0.1.0 release only. It is recorded
beside `dev/work/2026-09-09-observation-scale-guard/intent.md` (commit
`bcccba3`), which carries the measurements and the design constraint.

## Alternatives considered

| Option | One-line description |
|---|---|
| A | Extend the guard now with a scale-aware test, on the realised Gram or on the per-step collapse of the posterior variance |
| B (CHOSEN) | Correct the documentation to state the limitation truthfully, pin the behaviour with a test, and queue the guard as a post-release work package |
| C | Record rather than refuse: emit a classed warning and a `meta` event when the posterior variance collapses on the implicit route |
| D | Ship as it stands, leaving the help's claim that such a case fails later as integration instability |

## Rationale for rejection

### Option A

Not wrong, and not now. `rcond()` supplies invariance under a rescaling of the
Gram for nothing, and a scale test does not: the threshold acquires units, and
has to stay invariant under a rescaling of the state and under a change of time
unit. A legitimately precise sensor is a valid model rather than a contract
violation, so whether the test refuses or warns is unresolved. The measured
margin shows the change would be *safe* -- worst single-step variance drop 1.01
on the dyad graded record and 1.277 on the ENSO configuration, against about
1e14 for the counterexample -- but safe is not the same as designed, and a
threshold chosen under release pressure is the kind of constant nobody can
defend a year later. Deferred with its measurements preserved, not dropped.

### Option C

The cheaper half of A, and the better value of the two if the release had room
for one behavioural change: it converts silence into signal without deciding
whether the model is admissible. It was rejected for this release because it
still requires the same collapse threshold, which is the part that needs the
design pass. It remains the likely shape of the queued work package.

### Option D

Rejected outright. The help would have continued to promise that a near-zero
observation noise fails later as integration instability, which is true on the
explicit stepper and false on the implicit one. A false assurance is worse than
the silence it replaced, and a claim that outruns the implementation is the
defect class the audit exists to find.

## Implemented option (B)

The two lines of the `regularize` section that promise a later failure are
replaced with what the routes actually do, including the measured
counterexample, and a test pins the implicit route completing on `Sx1 = 1e-8` so
the sentence and the behaviour stay coupled. The claim drifted in the first
place because nothing exercised it.

Parity is untouched by construction. B changes no code, so no number can move;
the graded surface, its 28 fixture-backed rows and their tolerances are
unchanged, and the branch's suite of 7,135 assertions was verified green before
the decision. The parity claim is a claim about numbers on the graded surface,
which is what makes the documentation-only route sufficient here.

## Forward cost

- 0.1.0 ships a documented route -- the implicit stepper on a near-degenerate
  scalar observation model -- that returns a confident posterior with no
  condition raised and no event recorded in `meta`.
- The truth of the help text now rests on a test that pins a limitation. When
  the guard lands, that test and that sentence move together, or the
  documentation drifts a second time.
- The queued work package inherits a constraint that cannot be traded away: the
  guard may refuse or record, and may never alter a value on an input it
  accepts. A repair-style implementation -- flooring the Gram, adding jitter,
  changing the solver or the gain -- would move numbers on admissible inputs and
  end the parity claim.
- A model accepted today and refused once the guard lands is a breaking change,
  and carries its own NEWS entry then.
- The margin measurement covers the dyad graded record and the ENSO
  configuration. The predator-prey and partition records are unmeasured, and are
  swept before any threshold is fixed.

## References

### Methodological

- Golub, G. H. and Van Loan, C. F. (2013). *Matrix Computations*, 4th edition. Johns Hopkins University Press. The reciprocal condition number is a relative measure, invariant under a uniform rescaling of the matrix, which is why it detects a rank-deficient Gram and not a vanishing one.

### Empirical

- `acir/R/aci-conditional.R` -- `.check_gram_path()`, the realised-Gram test added on the corrections branch, and the threshold it shares with the constructor.
- `acir/R/aci-kernels-scalar.R:365-368` -- the explicit route's only covariance tests, `!is.finite(R)` and `R <= 0`; the implicit route preserves positivity, so neither fires on a collapse.
- Verification runs, 2026-09-09: with `Sx1` falling at one interval start of a `dt = 0.01` grid, the explicit route refuses at 1e-4, 1e-8 and 1e-30 with `aci_error_covariance_not_spd`, while the implicit route completes with the mean at the noiseless limit and sd 1e-7, 1e-29 respectively; both refuse at 1e-160, where the Gram underflows to exactly zero.
- Margin, same runs: worst single-step variance drop 1.01 on the dyad graded record and 1.001 / 1.004 / 1.277 on the ENSO documented configuration, against about 1e14 and 1e58 for the two counterexamples.
- Branch suite `c6481ce`: 7,135 assertions passing, 0 failures, 1 environment-gated skip.

### Operational

- Work package: `dev/work/2026-09-09-observation-scale-guard/intent.md`, commit `bcccba3`.
- Cairn: `2026-09-02-release-0.1.0-parity-milestone.cairn.md` (this project) -- the parity claim this decision protects.
- Branch: `acir-0.1.0-corrections` at `c6481ce`, which added the exactly-zero guard.

<!-- generated_by: ledger add 1.0.0; schema: 1.4; ts: 2026-09-09T02:34:23Z -->
