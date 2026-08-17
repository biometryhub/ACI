---
id: 2026-08-14-censored-range-not-missing
schema_version: 1.4
date: 2026-08-14
tier: T3
classification: open
export_status: local_only
consensus_mode: none
domain: r-package
project: aciR
status: accepted
title: "An unresolved causal influence range is right-censored, and is reported as a bound"
tags: [reporting, censoring, honest_reporting, external_review]
triggers:
  - public_api_change
  - methodology_change_downstream
reversal_cost: low
decision_pressure: null
review_due: null
review_trigger: "if a calibration study across the dyad, predator-prey and ENSO shows margin = 0.1 is the wrong default"
supersedes: []
superseded_by: null
related:
  - 2026-08-13-cir-saturation-margin
  - 2026-08-14-reference-quadrature-closure
---

## Status

Decided 2026-08-14 on the maintainer's ruling, following an external review.

A reported time whose range consumes more than `1 - margin` of the sequence it
was measured against now returns the computed value as a **right-censored lower
bound**, with the caveat in a new four-valued `status` (`"resolved"`,
`"censored"`, `"below_threshold"`, `"insufficient"`). It previously returned
`NA`. `saturated` is retained and equals `status == "censored"`, so callers
testing it keep working; what changed is that there is now a number beside it.

Only `"insufficient"` -- fewer than three later observations, where no
quadrature exists -- still returns `NA`.

Measured on the arbitrary cross-noise dataset at the reference horizon: **1049
of 1051 times now carry a finite bound where all 1049 were holes.**

The reasoning is the reviewer's and is correct. A time whose last-exit consumes
91 to 100 per cent of the available record is not an absent measurement. It is
a censored one, and the record supports the statement that the range is at
least this long. Returning `NA` discarded the only true statement available,
and did so at the end of the record, which is where a user studying a recent
event looks.

`subjective_censored`, a logical matrix, marks the individual thresholds that
ran long, because resolution is judged per quantity and a demanding threshold
outruns the objective range without condemning it.

## Alternatives considered

**1. Keep returning `NA`.** The status quo: an unresolved range is no range.

**2. Return the value with no flag**, as the reference does.

**3. Return the bound and drop `saturated`**, leaving only `status`.

**4. Report a Kaplan-Meier-style curve** of subjective range against threshold
with censoring times marked.

## Rationale for rejection

**1** conflates "we did not measure this" with "we measured a lower bound",
which are different statements, and picks the less informative one. It was
chosen originally as the conservative option; on reflection, discarding
information is not conservatism.

**2** is the reference's behaviour and is the reason a reader of its output
cannot tell a resolved range from a truncated one. The flag is the whole
contribution.

**3** would have been cleaner but breaks every caller testing `saturated`, for
a rename. Retaining it as a derived field costs one line.

**4** is the reviewer's own "that is a paper, not a patch". Agreed.

## Forward cost

`aci_cir()` returns two more fields, so anything consuming the object
positionally breaks -- nothing in the package does. The default `margin` of 0.1
is unchanged and remains a judgement rather than a calibrated constant; the
reviewer's recommendation of a calibration figure across the dyad, the
predator-prey pair and ENSO stands as open work, and is the review trigger on
this cairn.

`summary.aci_cir()` reports the censored fraction, so the caveat is visible at
the prompt rather than only in the object.

## References

- `aciR/R/aci-cir.R` -- the status taxonomy
- `aciR/R/aci-cir-methods.R` -- `summary.aci_cir()` reporting
- `aciR/tests/testthat/test-cir.R`, `test-cir-horizon.R`
- `tools/design/2026-08-14_reviewer_OUTPUT.md` section 3(b), S2
