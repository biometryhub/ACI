---
id: 2026-08-13-vector-core-dispatch-and-cross-noise-gap
schema_version: 1.4
date: 2026-08-13
tier: T3
classification: open
export_status: local_only
consensus_mode: none
domain: r-package
project: aciR
status: accepted
title: "Vector core added beside the scalar one, and its cross-noise path grounded analytically"
tags: [vector_states, conditional_aci, oracle_scope, rb22, public_api]
triggers:
  - public_api_change
  - adr_class_commitment
reversal_cost: high
decision_pressure: null
review_due: null
review_trigger: "if an authors'-reference implementation of the matrix cross-noise terms ever appears upstream, the grade should be redone against it and this cairn's central caveat retired"
supersedes: []
superseded_by: null
related:
  - 2026-08-13-l-y-admits-observed-state-dependence
  - 2026-08-13-online-smoother-log-reconstruction
---

## Status

Decided and implemented 2026-08-13 as stage R3, adding vector-valued states and
the conditional construction. Shipped at version 0.1.0.9000.

Two decisions are recorded together because they are one design: where the
vector code lives relative to the scalar code, and what its most exposed path is
graded against. The second is the one that matters, and it is the reason this
cairn exists rather than a commit message.

A spike preceded the work: the matrix core, instantiated at one dimension,
reproduced the dyad oracle at 4.574119e-14 -- bit-identical to the scalar core,
component for component -- which established that the generalisation is exact
before any of it was shipped.

## Alternatives considered

| Option | One-line description |
|---|---|
| A | Replace the scalar core with the general one, since at one dimension they agree |
| B (CHOSEN) | Keep both; the public functions dispatch on the shape of the components |
| C | A separate family of vector-only exports |
| D | Ship the vector core but decline to expose the noise cross-covariance until an upstream grading exists |

## Rationale for rejection

### Option A

The tempting one, and wrong. The two do agree at one dimension -- that is
graded, not assumed -- but the scalar recursion is about thirty times faster
and is the code the package's oldest and strongest oracle grades. Replacing a
validated, fast implementation with a general, slow one that computes the same
numbers spends an asset and buys nothing. The user-visible consequence would
have been a thirty-fold slowdown on the flagship path in exchange for less
evidence.

### Option C

Would have doubled the vocabulary -- `aci_filter()` beside
`aci_filter_matrix()` -- for a distinction the user does not care about. The
question a user asks is the same in both cases; only the shape of their system
differs, and the components list already carries that shape. Dispatching on it
is the honest expression of that.

### Option D

Seriously considered, and rejected on the grounds that concealment is worse
than disclosure. The matrix cross-noise terms have no upstream grading and will
not acquire one, since the reference contains no model that exercises them.
Withholding the capability until an impossible condition is met would have left
users to write the terms themselves, ungraded and undisclosed, which is
strictly worse than shipping them with the gap stated. The gap is recorded in
the manifest's `does_not_grade`, in the validation vignette, and here.

## Implemented option (B)

The discriminator is the shape of the latent-noise covariance: a matrix selects
the vector path, a bare number the scalar one. A one-by-one matrix routes to
the vector path deliberately -- someone who wrote their system in matrices gets
the matrix recursions, and the two agree.

Three things needed more care than they did in one dimension. The covariance is
re-symmetrised at each step, because an explicit Euler step breaks symmetry at
round-off and the asymmetry compounds. Positive-definiteness is tested by
attempting a Cholesky factorisation, which is the definition rather than a
proxy. And the relative entropy is evaluated through Cholesky factors
throughout.

**That last point is not stylistic, and it is the finding worth carrying
forward.** The trace term of the multivariate relative entropy is
`tr(R_f^-1 R_s)`. Written as `tr(c_f^-1 R_s c_f^-1)` it is wrong -- but it
agrees for DIAGONAL covariances, so no scalar test, and no block-diagonal test,
can distinguish the two. The spike's block-diagonal grading passed with the
wrong expression in place. The matrix cross-noise fixture, which carries genuine
off-diagonal structure, failed on its first run and produced a negative relative
entropy, which the existing domain guard refused to return.

The conditional construction is the reference implementation's: inflating a
component's observational uncertainty without bound sends its filter weight to
zero, implemented by supplying an inverse noise Grammian supported only on the
target block. This is why the schema admits `S_xoS_x_inv` directly -- the
object handed in is deliberately not the inverse of any covariance.

## Forward cost

- **The package now exposes a path with no authors'-reference grounding, and
  must keep saying so.** The matrix cross-noise terms are graded by a second
  independent implementation, which refutes a transcription error, and by
  block-diagonal algebraic Riccati fixed points, which depend on neither
  implementation and are the primary grounding. Neither is an authors'
  reference, and no future summary of the package's evidence may imply
  otherwise.
- **Two cores to keep in step.** They agree at one dimension today, and a test
  holds them to it against the dyad fixture. Any change to the scalar
  recursion must be mirrored, and the failure mode is a silent divergence
  between two paths a user cannot tell apart.
- **A diagonal test can no longer be trusted to grade a matrix expression.**
  The trace defect is the standing example, and the regression test that pins
  it asserts explicitly that the wrong expression would have produced a
  different number -- so the test cannot itself decay into grading nothing.
- **The online smoother and the causal influence range remain scalar-only,**
  and the gap is now more visible because everything around them generalised.
  The ordered products the online smoother reduces to cumulative logarithms do
  not commute once they are matrices, so that path needs a fresh derivation
  rather than a widening.
- **Conditioning changes the estimand.** Every downstream summary must carry
  that, or the conditional metric will be read as a sharper answer to the
  unconditional question.

## References

### Methodological

- Andreou, M., Chen, N. and Bollt, E. (2026). *Assimilative causal inference*.
  Nature Communications, 17, 1854.
  https://doi.org/10.1038/s41467-026-68568-0 -- conditional ACI, Supplementary
  Information sections 1.2.1 and 1.4; the matrix CGNS state estimation
  equations, section 2.1.2.

### Empirical

- `aciR/R/aci-core-mv.R` -- the vector filter, smoother and multivariate
  relative entropy.
- `aciR/R/aci-conditional.R` -- the masked inverse Grammian.
- `aciR/tests/testthat/test-identities-mv.R` -- the analytic grounding: the
  block-diagonal Riccati fixed points with non-zero cross-noise, and the
  regression test pinning the trace term.
- `matlab_reference/ENSO_model_cond_ACI_*.m` -- thirty-five occurrences of
  "NOISE CROSS-INTERACTION TERMS ARE ABSENT FROM THIS MODEL", the measured
  basis for the claim that no upstream grading exists.
- Measured 2026-08-13 against `oracle/aci_oracle_mv.m`: maximum absolute error
  5.773160e-15 across all ten compared series, with a cross-covariance
  Frobenius norm of 3.990144e-01.
- Measured 2026-08-13 (spike): the one-dimensional collapse reproduced the dyad
  oracle at 4.574119e-14, and a block-diagonal two-dimensional system
  reproduced two independent scalar runs at 0.00e+00.

### Operational

- Cairn: `2026-08-13-l-y-admits-observed-state-dependence.cairn.md` (this
  project) -- the widening this one builds on.
- Cairn: `2026-08-13-online-smoother-log-reconstruction.cairn.md` (this
  project) -- the path that did NOT generalise here, and why.
- `aciR/inst/extdata/oracle-manifest.yml` -- fixture `mv`, whose
  `does_not_grade` field carries the central caveat.
