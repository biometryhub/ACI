---
id: 2026-08-13-cir-saturation-margin
schema_version: 1.4
date: 2026-08-13
tier: T3
classification: open
export_status: local_only
consensus_mode: none
domain: r-package
project: aciR
status: accepted
title: "Unresolvable causal-influence-range times return NA rather than the record-limited value"
tags: [causal_influence_range, public_api, divergence_from_reference, interpretation]
triggers:
  - public_api_change
  - methodology_change_downstream
reversal_cost: medium
decision_pressure: null
review_due: null
review_trigger: "author review of aci_cir(); this is a deliberate divergence from the reference and has not yet been agreed with anyone but its implementer"
supersedes: []
superseded_by: null
related:
  - 2026-08-13-online-smoother-log-reconstruction
---

## Status

Decided and implemented 2026-08-13 while building `aci_cir()`; shipped in
commit `08bedb6` on `main` at version 0.1.0.9000. **This decision has been
made by the implementer alone and is flagged for author review**; it is the
one place in the causal-influence-range work where aciR deliberately behaves
differently from the reference implementation, and it is recorded here
partly so that the divergence is impossible to inherit unknowingly.

The causal influence range at time *t* is measured by comparing the
posterior informed by the whole record against the posterior informed only
up to each later observation. A time near the end of the record is therefore
compared against almost nothing, and returns a small number -- which reads
as a short, confident range rather than as an absence of evidence.

## Alternatives considered

| Option | One-line description |
|---|---|
| A | Report every computed value, as the reference does, and leave the reader to notice |
| B (CHOSEN) | Mark a time whose range consumes more than a retained fraction of the record as saturated and return `NA` |
| C | Trim the reported window by a fixed lookahead, as the reference's plotting code does |
| D | Refuse the computation outright for times too close to the end |

## Rationale for rejection

### Option A

Faithful to the reference and the cheapest to defend on those grounds, but
it exports a known trap. The reference is a set of analysis scripts whose
author knows where the artefact lives; a package is used by people who do
not, and the failure is silent and directional -- it always biases toward
shorter ranges, at exactly the end of the record where a user studying a
recent event is most likely to look. A number that is wrong in a known
direction, with nothing marking it, is worse than a missing number.

Would have been chosen if the goal were bit-faithful reproduction rather
than a reusable instrument. It is still available: `margin` may be set
arbitrarily small, and the oracle test does exactly that to compare like
with like.

### Option C

The reference's own workaround: compute over a window, then plot only the
part with a fixed lookahead of clear record beyond it. Rejected because a
fixed lookahead in time units is not transferable -- the right amount
depends on how fast the divergence decays, which depends on the system and
on the threshold. Encoding one number that happened to suit the dyad would
be a magic constant dressed as a method.

### Option D

Over-strict, and it would have made the boundary a hard error rather than
an observation about the data. The range near the end of a record is not
ill-posed; it is under-determined, which is what `NA` means.

## Implemented option (B)

`margin` (default 0.1) states the fraction of each comparison sequence that
must remain unused for that time's range to count as resolved. Resolution is
judged **per quantity, not per time**: the subjective range at a demanding
threshold outruns the objective range by a wide margin, and condemning the
objective range because the strictest threshold was unresolved would discard
the quantity the method leads with. So the objective range is judged against
the decay of the divergence to the reporting threshold, and each subjective
threshold against its own count.

One implementation note that is easy to get wrong and was got wrong first:
the divergence is exactly zero at the final observation by construction, so
a range can never formally reach the end of its row, and a saturation test
written against the row's end never fires. The test is against the retained
margin.

## Forward cost

- **A defaulted number now sits in the public interface with no external
  warrant.** 0.1 is a judgement, not a derivation. It must stay documented
  as this package's own device, and `API_STABILITY.md` marks `aci_cir()`
  Experimental principally because of it.
- **The oracle grade must stand the margin down.** The reference reports
  everything, so comparing like with like means disabling the margin for the
  grade (`margin = 1e-9`) and grading the margin separately. Two test paths
  where a faithful transcription would have had one.
- **Any downstream comparison with the reference's published figures must
  account for this.** aciR will show gaps where the reference shows values.
  That is intended, and anyone reproducing a figure needs to be told.
- **If the author rejects this, the reversal is cheap but not free**: the
  returned object's shape, and every test asserting `NA` at late times,
  change together.

## References

### Methodological

- Andreou, M., Chen, N. and Bollt, E. (2026). *Assimilative causal
  inference*. Nature Communications, 17, 1854.
  https://doi.org/10.1038/s41467-026-68568-0 -- the causal influence range,
  eqs. (8)-(9).

### Empirical

- `aciR/R/aci-cir.R:181-206` -- the per-quantity resolution test.
- `aciR/tests/testthat/test-cir.R` -- "times the record cannot resolve are
  reported as unresolved"; asserts saturation affects late times and not
  early ones.
- `aciR/inst/extdata/oracle-manifest.yml` -- fixture `cir`,
  `does_not_grade`, which records that the margin has no upstream
  counterpart and is graded only by aciR's own tests.
- `matlab_reference/dyad_interaction_model.m:460-470` -- the reference's
  `lookahead_tolerance = 0.6`, the workaround rejected as option C.

### Operational

- Cairn: `2026-08-13-online-smoother-log-reconstruction.cairn.md` (this
  project)
- Commit: `08bedb6` (main, 0.1.0.9000)
- Review guide, station 7 and station 10 item 2 -- flagged to the author as
  the decision most worth disagreeing with.
