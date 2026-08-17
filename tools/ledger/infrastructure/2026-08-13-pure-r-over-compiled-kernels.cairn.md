---
id: 2026-08-13-pure-r-over-compiled-kernels
schema_version: 1.4
date: 2026-08-13
tier: T3
classification: open
export_status: local_only
consensus_mode: none
domain: infrastructure
project: aciR
status: accepted
title: "Package stays pure R; compiled kernels deferred to measurement and to the companion repository"
tags: [architecture, cran, dependencies, performance]
triggers:
  - adr_class_commitment
  - new_external_dependency
reversal_cost: medium
decision_pressure: null
review_due: null
review_trigger: "if a package-advertised use case is measured to be unusable in pure R -- particularly the vector CGNS core at ENSO resolution"
supersedes: []
superseded_by: null
related:
  - 2026-08-13-online-smoother-log-reconstruction
---

## Status

Decided 2026-08-13, reversing a decision taken earlier the same day. The
earlier decision adopted `cpp11` kernels for the quadratic online-smoother
recursion and the three-million-step ENSO filter, and was approved on the
scoping estimate. It was reversed when two later constraints landed
together: that the package targets CRAN first, and that reproduction at the
published scale moves to a companion repository rather than living in the
package.

Both halves of the original justification moved out of the package with
that second constraint. What remains inside the package runs at a scale
pure R handles.

## Alternatives considered

| Option | One-line description |
|---|---|
| A | Adopt `cpp11` now, as originally approved, and carry compiled code from the outset |
| B (CHOSEN) | Stay pure R; measure; add compiled kernels only if a package-advertised use case proves unusable |
| C | Ship compiled kernels as an optional accelerator with the R path as the graded reference |
| D | Put the compiled kernels only in the companion repository, never in the package |

## Rationale for rejection

### Option A

The justification did not survive the change of scope. The quadratic
recursion at package scale is on the order of `1e5` inner operations, which
R does comfortably; the three-million-step filter that motivated compiled
code belongs to the reproduction, which is no longer in the package.

The cost, by contrast, is permanent and recurring rather than one-off. A
compiled package on CRAN takes on portability obligations across compiler
versions and sanitiser configurations, and it would end the package's
current property of two base `Imports` and no compiled code -- which is a
real asset for a small methods package whose value proposition is
auditability.

Would be chosen the moment a measurement shows a package-advertised use
case cannot run. That is this cairn's review trigger, and the ordering is
the point: measure, then compile.

### Option C

Attractive on paper and rejected on maintenance grounds. Two
implementations of one recursion are two places for it to drift, and the
cross-check that would prevent the drift is an obligation, not a free
benefit -- it has to be written, kept fast enough to run, and never deleted
for being slow. Reasonable if the speedup were load-bearing; it is not, at
package scale.

### Option D

Nearly the chosen option and differs only in finality: D forecloses
compiled code in the package permanently, B defers it pending evidence.
Deferral is preferred because the vector CGNS core has not been built yet
and its cost is the widest remaining unknown in the roadmap. Choosing D
would be deciding that question before measuring it.

## Forward cost

- **The pure-R ceiling is now an accepted constraint.** The published
  causal-influence-range figure, quadratic in a window of about 20,600
  steps, is a batch computation and not an interactive one. Documented in
  `aci_cir()`; not guarded in code.
- **The companion repository inherits the performance problem.** Whatever
  reproduces the ENSO case study at full resolution must solve it there,
  and cannot assume the package will.
- **The reversal must stay visible.** An approved decision was reversed on
  changed constraints rather than on new information about the original
  question; this record is what prevents that reading being lost.
- **The measurement that would trigger reversal has not been taken.** The
  vector core does not exist yet, so the claim that pure R suffices is
  established for what is built and asserted for what is not.

## References

### Methodological

- Andreou, M., Chen, N. and Li, Y. (2026). *An adaptive online smoother
  with closed-form solutions and information-theoretic lag selection for
  conditional Gaussian nonlinear systems*. Journal of Nonlinear Science.
  https://doi.org/10.48550/arXiv.2411.05870 -- eq. 3.19, the spectral bound
  that makes the pure-R path affordable at all.

### Empirical

- `aciR/DESCRIPTION` -- `Imports: graphics, stats`; no `LinkingTo`, no
  `src/`.
- Measured 2026-08-13: full-lag online smoother over the 30,001-step dyad
  signal, 11.3 s in pure R with an effective lag of 15,711.
- Measured 2026-08-13: `R CMD check --as-cran` 0 errors, 0 warnings; whole
  suite of 399 assertions inside the check budget.

### Operational

- Cairn: `2026-08-13-online-smoother-log-reconstruction.cairn.md` (this
  project) -- the algorithmic decision that removed most of the need for
  compiled code.
- Commit: `08bedb6` (main, 0.1.0.9000)
