---
id: 2026-08-14-mixed-tolerance-comparison-criterion
schema_version: 1.4
date: 2026-08-14
tier: T3
classification: open
export_status: local_only
consensus_mode: none
domain: r-package
project: aciR
status: accepted
title: "Parity verdicts use a mixed absolute/relative criterion, with the relative tolerance left untouched"
tags: [numerics, tolerance, parity_harness, honest_reporting]
triggers:
  - adr_class_commitment
  - methodology_change_downstream
reversal_cost: low
decision_pressure: null
review_due: null
review_trigger: "if any quantity ever passes on the absolute term alone at a magnitude where the relative term should have carried it"
supersedes: []
superseded_by: null
related:
  - 2026-08-14-extraction-by-byte-verified-hoist
---

## Status

Decided and implemented 2026-08-14, in response to the parity harness's first
full run: 25 of 26 quantities within tolerance, one exceeding.

The exceedance was `online_mean(lag=Inf)` on the dyad record, at a declared
relative tolerance of 1e-11 and a measured maximum relative error of 4.16e-11.
It was investigated before anything was changed. The maximum occurs at index
2000, where the online smoothed mean is -1.678e-4 -- the only point in the
3001-step series with magnitude below 1e-3, against a series median of 0.864.
The absolute difference there is 6.99e-15, on a quantity formed by cancellation
between operands of order one, so it is roughly thirty units in the last place
of its inputs. The maximum absolute difference over the whole series is
1.53e-14.

The defect was therefore in the criterion, not in either implementation. A
purely relative test divides by the value, so where a series crosses zero it
divides by nothing in particular and reports a large number for a difference
that is small on every scale that matters.

The criterion is now `|a - b| <= atol + rtol * |b|`, elementwise, with
`atol = 1e-12` and every relative tolerance left exactly where it was. Both
`max_abs` and `max_rel` remain in the report alongside a `headroom` column --
the reciprocal of the fraction of the permitted budget the worst element used
-- so a reader can see which term carried each verdict and whether a pass was
comfortable. All 26 quantities pass; the tightest headroom is 143x.

A `SUSPECT-EXACT` verdict was added at the same time, for any quantity whose
maximum absolute difference is exactly zero. Two implementations in two
languages agreeing bit-for-bit is not a triumph, it is usually evidence that
one side did not run. No row currently triggers it, which is itself the
evidence that both sides genuinely execute.

## Alternatives considered

**1. Relax the relative tolerance to 1e-10** so the row passes.

**2. Judge on absolute error alone.**

**3. Exclude near-zero elements from the comparison** by masking `|value|`
below a floor.

**4. Report the exceedance and change nothing**, treating 4e-11 as the honest
measured agreement of that quantity.

## Rationale for rejection

**1** is the move this project must not make and the reason this cairn exists.
It would have produced a green table, and it would have raised the ceiling for
every element of that series -- including the 3000 well-conditioned ones where
1e-11 is the right bar -- to hide the behaviour of one. A tolerance chosen
after seeing the difference describes the difference; it does not test it.

**2** is scale-blind in the other direction. These quantities run over orders
of magnitude between models, and a fixed absolute bar that suits the dyad's
order-one means would be meaningless on the ENSO components.

**3** was rejected because masking removes the element from the grading
entirely, which is the failure this package is built around: a term that is
present, reachable and excluded from the check. The zero crossing is real data
and should be judged, on a criterion that can judge it.

**4** was considered seriously and is the most conservative option. Rejected
because the reported number would then be an artefact of the criterion rather
than a property of the code, and a reader comparing 4e-11 here against 5e-15
elsewhere would draw a conclusion about the online smoother that is not true.

## Forward cost

`atol` becomes a per-quantity declaration as the harness reaches models whose
natural scale is not order one; the default of 1e-12 is stated per comparison
rather than assumed globally. The `headroom` column must be watched rather than
merely printed -- a quantity whose headroom drifts toward 1 over time is
degrading even while it passes, and a quantity whose headroom is enormous
should be checked for reaching the code it claims to.

## References

- `tools/oracle/parity/tools/compare.R` -- the criterion and its justification
- `tools/oracle/parity/reports/parity_scalar.csv` -- 26 of 26, headroom per row
