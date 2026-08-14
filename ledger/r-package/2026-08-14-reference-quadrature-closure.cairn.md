---
id: 2026-08-14-reference-quadrature-closure
schema_version: 1.4
date: 2026-08-14
tier: T3
classification: open
export_status: local_only
consensus_mode: none
domain: r-package
project: aciR
status: accepted
title: "The odd-interval Simpson closure follows the reference, trading accuracy for transcription"
tags: [numerics, quadrature, transcription_fidelity, external_review]
triggers:
  - public_api_change
  - methodology_change_downstream
  - audit_finding_downstream
reversal_cost: low
decision_pressure: null
review_due: null
review_trigger: "if the package's charter ever puts accuracy of the integral above reproduction of the authors' number, this reverses"
supersedes: []
superseded_by: null
related:
  - 2026-08-14-mixed-tolerance-comparison-criterion
  - 2026-08-14-censored-range-not-missing
---

## Status

Decided 2026-08-14 on the maintainer's ruling that `objective` may **not**
differ from the reference at 1e-7, following an external review that identified
a residual I had wrongly closed.

`.aci_simpson()` closes an odd interval count with the quadratic through the
last three samples integrated over the last interval, which is the reference's
rule. It previously used a Simpson 3/8 panel over the last three intervals.
`objective` moved from **4.58e-09 to 1.37e-14** against the authors'
numbers -- the maximum over the full reported region of both graded datasets
(751 and 1001 times). An earlier draft quoted 3.57e-15, which was a maximum
over a SUBSET: it was measured at `margin = 0.001`, which censors more times
and so compares fewer of them. Quote the maximum, and say what it is over.

The 3/8 panel is the more accurate rule on an equally spaced grid, measurably:
against the exact integral its error is about a third of the reference's at
n >= 6. That accuracy was the original reason for choosing it. It is being
given up deliberately.

The change also removed a latent defect. The 3/8 panel assumed equal spacing,
and `objective_exact` integrates a **logarithmic** threshold grid. Neither the
129-point default nor the reference's 513 triggered it -- both have an even
interval count -- but any even-length `epsilon` would have. The replacement is
exact for unequal spacing, so the two problems closed as one change.

The new closure is **derived, not transcribed**: integrating the Lagrange basis
through the last three abscissae over the final interval, with a test asserting
it reduces to `h/12 * (-y0 + 8*y1 + 5*y2)` at equal spacing.

## Alternatives considered

**1. Keep the 3/8 panel** and document the 1e-7 as a designed difference.

**2. Offer both** behind `quadrature = c("reference", "simpson38")`.

**3. Keep 3/8 for `objective_exact`** and use the reference's closure only for
`objective`.

**4. Raise the oracle tolerance** to 1e-6 and call the residual closed.

## Rationale for rejection

**1** is defensible and was the reviewer's stated alternative. Rejected by the
maintainer's ruling: this package's claim is fidelity to the authors on the
quantity the method leads with, and an unexplained-looking 1e-7 on the headline
number costs more than the accuracy gain buys. A package that is *better* in a
way nobody asked for and *different* in a way everybody checks has chosen
badly.

**2** was the reviewer's first preference. Rejected as surface for a need
nobody has: two rules mean two sets of numbers to grade, two paths through the
oracle, and a user decision that has one right answer given (1)'s rejection.
Reconsider if a caller ever asks for the more accurate rule.

**3** mixes the two, which is the one thing the reviewer explicitly warned
against. Both quantities are compared against the same reference; splitting the
rule between them makes the residual harder to attribute, not easier.

**4** is forbidden by the mixed-tolerance cairn and would have hidden the
finding rather than closed it.

## Forward cost

`.aci_simpson()` is now exact for quadratics rather than cubics on an odd
interval count, and the quadrature tests say so explicitly instead of asserting
a cubic exactness the rule no longer has. Any future caller passing an odd
interval count inherits third-order behaviour on the interior and second-order
on the closing interval, which is the reference's behaviour.

The lesson that generalises is not about quadrature. I had eliminated the
Simpson hypothesis by argument -- the CIR integrand decays to zero at its
endpoint, so the closure term is multiplied by nothing -- and the argument was
true of a full-record row and false of a truncated one, which is the only case
that mattered. Measured at `j = 408`: the full-record row ends at exactly 0,
the truncated row at 9.25e-06, or 2.7e-04 of the peak. **Substituting an
argument for a measurement at the point where the argument is load-bearing** is
the same failure as a fixture that annihilates the term it grades, and it is
the third instance in this project.

## References

- `aciR/R/aci-quadrature.R` -- `.aci_simpson_closure()` and its derivation
- `aciR/tests/testthat/test-quadrature.R` -- reduction and exactness assertions
- `design/2026-08-14_reviewer_OUTPUT.md` section 3(a)
- `design/artefacts/2026-08-14_cir_residual_rows.csv` -- row `index == 408`
