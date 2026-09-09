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

Decided 2026-09-09 following the initial-release corrections at `c6481ce`.
The realised-Gram check closes the exactly-zero case. Small positive noise can
still produce finite positive implicit results without a warning or event,
including results on an unresolved time step. The decision covers 0.1.0 only.

The documentation and tests implementing option B are included with this update;
they were absent from `f7d2d26`, although this record described them as complete.
`dev/work/2026-09-09-observation-scale-guard/intent.md` records the deferred
diagnostic and the controls its design must address.

## Alternatives considered

| Option | One-line description |
|---|---|
| A | Extend the guard now with a scale-aware test, on the realised Gram or on the per-step collapse of the posterior variance |
| B (CHOSEN) | Correct the documentation to state the limitation truthfully, pin the behaviour with a test, and queue the guard as a post-release work package |
| C | Record rather than refuse: emit a classed warning and a `meta` event when the posterior variance collapses on the implicit route |
| D | Ship as it stands, leaving the help's claim that such a case fails later as integration instability |

## Rationale for rejection

### Option A

Deferred because the proposed variance-drop threshold cannot distinguish valid
information from numerical error. The static-state control in
`acir/tests/testthat/test-40-observation-scale.R` has an accurate `1e14` drop,
while the unresolved stationary-covariance control has about 81% error with a
largest one-step drop of only 5.16. The earlier dyad/ENSO margins do not establish
that a threshold is safe for other valid inputs. A future diagnostic needs a
defined purpose, justified threshold, matrix coverage and invariance under the
relevant changes of units. No new refusal rule is selected.

### Option C

A recorded warning could draw attention to a numerical concern without
declaring the model invalid. It still needs a defined statistic and threshold,
and must not imply that an unflagged result is accurate. That design remains
deferred; it is a more suitable starting point than a mandatory refusal.

### Option D

Rejected because the help promised a later failure that need not occur on
either stepper. The explicit route's behaviour depends on coupling and step
size; implicit positivity also does not establish time resolution or accuracy.

## Implemented option (B)

The `regularize` help in `acir/R/aci-assimilation.R` and generated
`acir/man/aci_filter.Rd` no longer promises a later failure. It explains that
small positive noise can pass the conditioning check, that finite positive
results may have no warning/event, and that time-resolution sensitivity needs
assessment. `acir/NEWS.md` records the correction.

`acir/tests/testthat/test-40-observation-scale.R` covers the implicit near-zero
return for scalar and matrix states, a weak-coupling explicit return, a valid
large variance drop, consistent changes of units, and an analytic covariance
case whose error decreases with refinement. These tests describe a limitation
and its controls; they do not certify the inaccurate output as a reference
result. Existing exactly-zero refusal tests remain in `test-35`.

No executable package expression, default, fixture or existing tolerance changes.
Numerical comparisons retain their existing scope; this does not establish
universal MATLAB equivalence or validate a future guard.

## Forward cost

- Small positive observation noise can still produce inaccurate positive
  covariances without a warning/event on an unresolved implicit step. Valid
  models can also produce small covariances accurately.
- The truth of the help text now rests on a test that pins a limitation. When
  the guard lands, that test and that sentence move together, or the
  documentation drifts a second time.
- A future diagnostic must preserve accepted-path values. Refusal changes the
  supported input domain; flooring the Gram or changing the solver/gain can
  change numbers and requires separate review and renewed comparisons.
- A model accepted today and refused once the guard lands is a breaking change,
  and carries its own NEWS entry then.
- The original margin summary covers dyad and ENSO only and lacks full run
  provenance. All retained graded cases need a reproducible sweep before any
  threshold is selected, alongside independent valid and inaccurate controls.

## References

### Methodological

- Golub, G. H. and Van Loan, C. F. (2013). *Matrix Computations*, 4th edition. Johns Hopkins University Press. The reciprocal condition number is a relative measure, invariant under a uniform rescaling of the matrix, which is why it detects a rank-deficient Gram and not a vanishing one.

### Empirical

- `acir/R/aci-conditional.R` -- `.check_gram_path()`, the realised-Gram test added on the corrections branch, and the threshold it shares with the constructor.
- `acir/tests/testthat/test-40-observation-scale.R` -- complete model inputs and analytic controls. From the repository root: `Rscript -e 'devtools::test("acir", filter = "40-observation-scale")'`.
- Extreme-scale correction: locally, `1e-160` squared is positive (about `1e-320`), but R/LAPACK returns `rcond = 0`; at `1e-162` the squared Gram actually becomes zero. The earlier underflow explanation conflated these cases.
- Historical margin summary: dyad 1.01 and ENSO 1.001 / 1.004 / 1.277. Complete run inputs and raw outputs were not recorded here; these values alone do not establish a safe threshold.
- Historical branch suite `c6481ce`: reported 7,135 passes and one environment-gated skip. The untouched `f7d2d26` review run had 7,155 passes and no skips with `ACI_ORACLE_PARITY_ROOT` pointing at the local evidence. New validation is separate from those baseline runs.

### Operational

- Work package: `dev/work/2026-09-09-observation-scale-guard/intent.md`, commit `bcccba3`.
- Cairn: `2026-09-02-release-0.1.0-parity-milestone.cairn.md` (this project) -- the parity claim this decision protects.
- Branch: `acir-0.1.0-corrections` at `c6481ce`, which added the exactly-zero guard.
