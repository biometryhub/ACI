---
id: 2026-08-14-extraction-by-byte-verified-hoist
schema_version: 1.4
date: 2026-08-14
tier: T3
classification: open
export_status: local_only
consensus_mode: none
domain: infrastructure
project: aciR
status: accepted
title: "Reference kernels become callable by byte-verified hoist, and the generated functions are not committed"
tags: [independent_oracle, parity_harness, provenance, licensing]
triggers:
  - adr_class_commitment
  - test_infrastructure
reversal_cost: low
decision_pressure: null
review_due: null
review_trigger: "if the parity harness moves to a repository that vendors the reference, the non-commit half of this decision should be revisited"
supersedes: []
superseded_by: null
related:
  - 2026-08-14-minimal-substitution-capture-gate
  - 2026-08-13-grading-register-as-build-gate
---

## Status

Decided and implemented 2026-08-14. The reference publishes seven top-to-bottom
scripts and no function library -- three files are callable, and none of them
carries any of the method. `simps.m` and `legendUnq.m` are third-party
FileExchange code; `progress_bar.m` is the authors' own console progress
indicator.

> **Amended 2026-08-16.** The paragraph above previously read "only `simps.m`
> and `legendUnq.m` are callable, and both are third-party". A coverage audit
> against the reference source found `progress_bar.m`, a third callable file,
> written by the authors and called from all six model scripts. The decision
> this cairn records is unaffected: a progress indicator is not a kernel, so
> the reference still offered nothing to grade `aci_filter()` against.
To run the authors' algebra on a dataset other than the one their script
generates, the computational passages have to become functions.

They are hoisted, never retyped. `tools/extract.R` emits each function as a
signature, a byte-exact slice of the reference between `BEGIN VERBATIM` /
`END VERBATIM` markers, and `end`; `check_extraction()` re-reads the slice from
the source and fails on a single byte of drift, and separately refuses any
executable line outside the markers. Seven functions are declared for the dyad,
totalling 343 verbatim lines.

Gate G1 settles empirically what the mechanism only claims: each extracted
function is called with the inputs the script itself used, and every output must
EQUAL -- not approximate -- the value the script itself produced. **58 of 58
outputs across two profiles, maximum absolute difference 0.**

The generated functions are NOT committed. `matlab_reference/` is already
gitignored in this repository as re-clonable external code, and the extracts
reproduce it; `manifest/extracts.dcf` and the generator are committed instead,
so the extraction is reproducible from an upstream clone in one command.

A `ConsumedOnly` assertion is part of the manifest: it declares that a slice
touches an input at one index and nowhere else, and the extractor verifies it
against the slice's code. This is what makes it legitimate to pass a scalar
initial mean where the reference's filter expects the true latent path -- the
filter reads it only as `y(1)`, and that is checked rather than believed.

## Alternatives considered

**1. Keep the existing hand-transcribed harnesses** in `tools/oracle/aci_oracle_*.m`
and call them the MATLAB side.

**2. Re-implement the kernels in MATLAB from the paper**, as clean functions
with designed signatures.

**3. Hoist, but commit the generated `.m` files** with the upstream MIT notice.

**4. Refactor the reference scripts in place** into functions, and run those.

## Rationale for rejection

**1** is the status quo and the reason this work exists. Those harnesses are a
transcription, and the same person transcribed them and wrote the R port; an
error common to both is invisible to every fixture the package owns. Notably
the transcription turned out to be faultless -- see Forward cost -- but that is
a result, not something that could be assumed.

**2** is the worst option available and would have looked the most
professional. A third implementation by the same author grades nothing; it
merely produces a third opportunity for the same misreading.

**3** is defensible -- MIT permits redistribution with the notice, and the
existing harnesses already carry it. Rejected for consistency: this repository
already treats the reference as re-clonable rather than vendored, and a
generated artefact reproducing a third of a file is better regenerated than
stored. The rendered parity document does quote the excerpts side by side with
the notice, which is where a reviewer actually needs to see them.

**4** was rejected because it destroys the thing being verified. Once the
reference file is edited, "byte-identical to the reference" has no referent.

## Forward cost

Each further model costs a manifest record per kernel plus a G1 run; the
generator and checker are written. Upstream re-clones must be followed by
`check_extraction()`, which will fail loudly if line numbering moved -- that is
the intended behaviour, and it is why the reduction preserves line counts.

The first payoff was immediate and is worth recording, because it retires a
standing doubt rather than confirming one. The committed dyad fixture
(`aciR/inst/extdata/dyad_reference.csv`) was re-graded against the authors' own
code at the PUBLISHED N = 30000 -- possible because the setup, signal, filter
and smoother are all O(N) even though the online smoother and CIR are not --
and all eight quantities over 301 rows agreed with maximum absolute difference
**exactly 0**. The transcription introduced no error. The fixture the package
grades against is the authors' output, and that is now a measurement rather
than a hope.

## References

- `tools/oracle/parity/manifest/extracts.dcf` -- declared ranges, inputs, outputs
- `tools/oracle/parity/tools/extract.R` -- generator and byte checker
- `tools/oracle/parity/tools/gate.R`, `tools/oracle/parity/matlab/bitcompare.m` -- gate G1
- `tools/oracle/parity/reports/g1_small.csv`, `g1_tiny.csv` -- 58/58 at difference 0
- `matlab_reference/LICENSE` -- MIT, Copyright (c) 2025 Marios Andreou
