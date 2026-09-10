---
id: 2026-09-10-release-020-guards-and-scoped-claims
schema_version: 1.4
date: 2026-09-10
tier: T3
classification: open
export_status: local_only
consensus_mode: none
domain: r-package
project: aciR
status: accepted
title: "acir 0.2.0 releases the guards and the scoped evidence claims, with the near-zero limitation documented rather than closed"
tags: [release, versioning, guards, evidence, governance]
triggers:
  - adr_class_commitment
  - release_reconciliation
  - breaking_change_with_deprecation
reversal_cost: medium
decision_pressure: publication
review_due: null
review_trigger: "at the next release, or when the observation-scale guard work package lands and changes what is refused"
supersedes: []
superseded_by: null
related: [2026-09-02-release-0.1.0-parity-milestone, 2026-09-09-near-zero-gram-documented-not-guarded]
---

## Status

Released 2026-09-10. Tag `acir-v0.2.0` is annotated, sits on `6461892` (the
merge of the release pull request), and its tree is byte-identical to the
candidate that returned `Status: 1 NOTE` from `R CMD check --as-cran` locally,
where the single note is `New submission`. Ten checks green on that commit, read
back from the host API rather than inferred. `main` then opened at `0.2.0.9000`
and branch protection was enabled and read back the same way.

## Alternatives considered

| Option | One-line description |
|---|---|
| A | Patch release 0.1.1, treating the work as documentation corrections |
| B (CHOSEN) | Minor release 0.2.0, on the ground that three changes stop working calls from working |
| C | Hold the corrections unreleased until the observation-scale guard closes the near-zero case |
| D | Remove the retired `T` argument at 0.2.0, as the expiring deprecation notice promised |

## Rationale for rejection

### Option A

`API_STABILITY.md` states that at 0.x the interface may change with the change
announced under a "Breaking changes" heading. Three of this release's items make
previously working calls fail: the declared R floor moved to 4.1.0 so
installation on 4.0.x is refused; a degenerate observation-noise Gram is refused
wherever it occurs, where before it returned numbers; and overflow in a
recursion raises rather than returning a silent non-finite value. A patch
release asserts that nothing a caller wrote stops working, which would be false.

### Option C

The near-zero case cannot be closed by a threshold chosen under release
pressure -- the reasoning is in the companion cairn. Holding an otherwise
finished release for it would have kept the corrected claims off the public
record for the sake of a limitation that is now stated in the help, in NEWS, in
five tests and in a work package. Correct claims shipped beat correct claims
withheld.

### Option D

The notice said `T` would be accepted "until 0.2.0", so arriving here forced a
choice. Removing it would break callers at the same release that already carries
three breaks, for an alias that costs nothing to keep. The promise was retired
instead of the argument, and NEWS records that no removal version is scheduled.
The cost is a deprecation with no end date, which the 1.0 planning inherits.

## Implemented option (B)

Version, citation metadata and the software-metadata record all read 0.2.0 and
R 4.1.0. `NAMESPACE` is byte-identical to the `acir-v0.1.0` tag, so no exported
surface moved. The fixture and evidence subtrees hash identically to the
pre-release `main`, and the only R file the release touched carries a
deprecation message and its roxygen. No graded quantity changed.

## Forward cost

- The release ships a documented route on which the implicit stepper returns a
  collapsed posterior from a near-degenerate observation model, and a flooring
  policy that can return a value of order 1e21 with a warning. Both are stated;
  neither is closed.
- The `T` alias now has no scheduled removal. Semantic versioning is promised
  from 1.0, so the removal needs a home in the 1.0 plan or the promise erodes.
- Branch protection binds the ordinary workflow from this release onward: each
  author now needs the other's approval, and administrators may bypass. The
  process record was corrected to describe what the host actually enforces.
- The evidence register carries 68 rows of which 40 are behavioural-only, and
  the fixture producer route remains static and unexercised. Both are recorded
  gaps that a reader of the register will meet, and the next release is when
  they are most likely to be read closely.

## References

### Methodological

- Preston-Werner, T. *Semantic Versioning 2.0.0*. https://semver.org/ -- the 0.x clause under which a minor bump carries breaking change before 1.0.

### Empirical

- `acir/NEWS.md`, `# acir 0.2.0` -- the Breaking changes heading and its three items.
- `acir/API_STABILITY.md` -- the 0.x interface promise this release is measured against.
- Local release-candidate check, 2026-09-10: `R CMD build` then `R CMD check --as-cran` returned `Status: 1 NOTE` (`New submission`); tests, examples with `--run-donttest`, and all three vignette rebuilds passed; the installed candidate reported version 0.2.0 and `citation("acir")` followed it.
- Subtree hashes against the pre-release `main`: fixtures and evidence register identical; all of `R/` except `aci-model.R` identical.

### Operational

- Cairn: `2026-09-09-near-zero-gram-documented-not-guarded.cairn.md` -- why the near-zero limitation ships documented.
- Cairn: `2026-09-02-release-0.1.0-parity-milestone.cairn.md` -- the release this one succeeds.
- Work package: `dev/work/2026-09-10-release-0.2.0/plan.md`, and `dev/work/2026-09-09-observation-scale-guard/intent.md` for the deferred guard.
