---
id: 2026-08-13-online-smoother-log-reconstruction
schema_version: 1.4
date: 2026-08-13
tier: T3
classification: open
export_status: local_only
consensus_mode: none
domain: r-package
project: aciR
status: accepted
title: "Online smoother reconstructs update products in logarithms rather than storing the reference triangle"
tags: [online_smoother, numerical_methods, memory, causal_influence_range]
triggers:
  - public_api_change
  - adr_class_commitment
reversal_cost: medium
decision_pressure: null
review_due: null
review_trigger: "when the vector or conditional CGNS core lands -- ordered products of matrices do not commute, and the logarithmic reconstruction does not carry over unchanged"
supersedes: []
superseded_by: null
related:
  - 2026-08-13-cir-saturation-margin
  - 2026-08-13-pure-r-over-compiled-kernels
---

## Status

Decided and implemented 2026-08-13 while building `aci_online_smoother()`,
the estimator the causal influence range is constructed from. Shipped in
commit `08bedb6` on `main` at version 0.1.0.9000. The decision was forced
by a resource wall rather than chosen freely: a faithful transcription of
the reference algorithm cannot run at the scale of the published figure on
ordinary hardware.

The reference implementation (`dyad_interaction_model.m`, upstream commit
`733c49fc`) stores the online smoother as nested staggered arrays -- one row
per observation, each row holding an estimate at every earlier time instant.
At the published figure's scale that triangle, plus the divergence matrix
built from it, needs roughly fourteen gigabytes.

## Alternatives considered

| Option | One-line description |
|---|---|
| A | Transcribe the reference faithfully; accept the quadratic memory |
| B (CHOSEN) | Accumulate the ordered update products as cumulative logarithms and truncate them once they fall below tolerance |
| C | Keep the quadratic algorithm but move it to compiled code so the constant factor is affordable |
| D | Expose only a short fixed lag, so the triangle never grows |

## Rationale for rejection

### Option A

Correct, and the only option needing no argument for its equivalence, but it
does not run. Fourteen gigabytes of intermediates to produce quantities that
are immediately reduced to a handful of scalars per time step is not a
tuning problem; it is the wrong shape. It would also have made the package
unusable inside `R CMD check`, where memory is constrained, and so
unshippable to CRAN with any example exercising it.

Would have been chosen if the quantities the triangle holds were themselves
wanted as output. They are not: every row is consumed by a reduction.

### Option C

Moves the wall rather than removing it. Compiled code buys a constant
factor, perhaps an order of magnitude, against a quadratic; the published
figure would still need gigabytes. It also carries a permanent cost --
a compiled package on CRAN takes on portability obligations across
compilers and sanitiser configurations -- for a benefit the logarithmic
reconstruction obtains for free. Rejected alongside the broader
pure-R posture recorded separately.

### Option D

Would have been honest but crippling: the fixed-lag family is the whole
point of the estimator, and the range is defined over the full sweep of
lags. Capping the lag short enough to bound memory would have capped it
below the range being measured, which is measuring the instrument rather
than the system.

## Implemented option (B)

The update applied to an earlier step when a new observation arrives is
damped by an ordered product of the per-step auxiliary matrices over a
contiguous range. The source paper (Andreou, Chen and Li, Theorem 3 and
eq. 3.19) bounds the spectral radius of every factor below one, so the
products decay geometrically and the influence of an observation on distant
earlier steps falls away exponentially.

Two consequences are exploited. First, a contiguous ordered product is a
difference of cumulative logarithms, so any product is recoverable in
constant time from two `O(n)` arrays -- a cumulative log-magnitude and a
cumulative sign. Second, once a product falls below tolerance the remaining
terms cannot change the result at double precision, so the accumulation
stops there. Work becomes `O(n * lag_effective)` and memory `O(n)`.

The difference is taken in logarithms and only then exponentiated. Forming
the two products separately would underflow; their ratio does not.

## Forward cost

- **The redesign must be graded, not assumed.** A literal transcription of
  the reference triangle is carried permanently in
  `tests/testthat/test-online-smoother.R` and the fast path is checked
  against it (measured 2.442e-15 at n = 400). That test is slow by design
  and may not be deleted for being slow: it is the only thing standing
  between the reconstruction and a silent divergence from the algorithm it
  replaces.
- **The tolerance is now a documented part of the interface.** `tol`
  defaults to 1e-18 and is a numerical tolerance on a converged geometric
  series, not a modelling choice. It must stay documented as such, or users
  will tune it as though it changed the estimand.
- **The argument does not survive the move to vector states.** Ordered
  products of matrices do not commute and cannot be reduced to differences
  of cumulative scalars. When the vector core lands, this file's fast path
  either gains a matrix-specific derivation or reverts to accumulation over
  a bounded lag. Recorded as this cairn's review trigger.
- **A factor of exactly zero must stay handled.** The logarithm carries it
  as negative infinity, and a range beginning past it would otherwise
  evaluate as an indeterminate difference. Such ranges are genuinely zero.

## References

### Methodological

- Andreou, M., Chen, N. and Li, Y. (2026). *An adaptive online smoother with
  closed-form solutions and information-theoretic lag selection for
  conditional Gaussian nonlinear systems*. Journal of Nonlinear Science.
  https://doi.org/10.48550/arXiv.2411.05870 -- Algorithm 1 (the recursion),
  eq. 3.17 (the ordered product), eq. 3.19 (the spectral bound this
  decision rests on), eqs. 3.22-3.23 (the fixed-lag form).
- Andreou, M., Chen, N. and Bollt, E. (2026). *Assimilative causal
  inference*. Nature Communications, 17, 1854.
  https://doi.org/10.1038/s41467-026-68568-0

### Empirical

- `aciR/R/aci-online-smoother.R:200-239` -- `.aci_online_aux()` and
  `.aci_online_product()`, the cumulative-logarithm reconstruction.
- `aciR/tests/testthat/test-online-smoother.R:46-87` -- the literal
  triangle transcription and the agreement test that grades the redesign.
- `aciR/inst/extdata/oracle-manifest.yml` -- `analytic_oracles` entries
  `online_smoother_zero_lag` and `online_smoother_triangle`, with measured
  values.
- Measured 2026-08-13: reconstruction against literal triangle 2.442e-15
  (n = 400) and 5.773e-15 (n = 800); against the reference MATLAB, online
  smoother mean 5.107026e-15 and covariance 7.216450e-16.

### Operational

- Cairn: `2026-08-13-pure-r-over-compiled-kernels.cairn.md` (this project)
- Cairn: `2026-08-13-cir-saturation-margin.cairn.md` (this project)
- Commit: `08bedb6` (main, 0.1.0.9000)
- Harness: `tools/oracle/aci_oracle_cir.m`
