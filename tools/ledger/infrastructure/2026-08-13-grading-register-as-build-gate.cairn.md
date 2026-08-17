---
id: 2026-08-13-grading-register-as-build-gate
schema_version: 1.4
date: 2026-08-13
tier: T3
classification: open
export_status: local_only
consensus_mode: none
domain: infrastructure
project: aciR
status: accepted
title: "Oracle scope tracked by a measured register that fails the build on an undeclared gap"
tags: [oracle_scope, rb22, build_gate, mechanism_over_rule]
triggers:
  - adr_class_commitment
reversal_cost: medium
decision_pressure: null
review_due: null
review_trigger: "if a gap is ever found that the register did not catch -- that would falsify the claim this cairn rests on"
supersedes: []
superseded_by: null
related:
  - 2026-08-13-vector-core-dispatch-and-cross-noise-gap
  - 2026-08-13-online-smoother-log-reconstruction
---

## Status

Decided and implemented 2026-08-13, after stages R1 to R4 had roughly tripled
the package's surface. Shipped at version 0.1.0.9000 as
`tests/testthat/test-grading-matrix.R`.

The problem it addresses is this package's founding one, returning at scale. A
term can be publicly reachable, execute on every run, be reported as covered by
line coverage, AND be graded by an oracle -- and still be checked against
nothing, because the fixture that grades it multiplies it by zero. That is what
happened to the scalar noise cross-covariance in 0.1.0, and it is the failure
mode of any oracle whose scope is narrower than the surface it appears to
grade.

Through R1 to R4 the answer to "what is actually graded?" was maintained by
hand, in the manifest's `does_not_grade` fields and in my own head. Two things
made that untenable. The surface grew from three consumers to five and from one
schema to two. And the R3 trace defect proved the hand-tracking was already
wrong once: a block-diagonal test was believed to grade the multivariate metric
and did not, because block-diagonal covariances have diagonal Cholesky factors.

## Alternatives considered

| Option | One-line description |
|---|---|
| A | Continue recording scope by hand in the manifest's `does_not_grade` fields |
| B (CHOSEN) | A register that MEASURES which terms each grading scenario exercises non-degenerately, and fails the build on an ungraded, undeclared cell |
| C | Require a per-term non-degeneracy assertion in each oracle test, without a register |
| D | Raise the coverage target and rely on it |

## Rationale for rejection

### Option A

What was already being done, and it had already failed once. Prose in a
manifest is a record of an intention to check; it cannot detect that a new
fixture annihilates a term, and it goes stale in the direction that flatters
the package -- a scope note written when a gap existed keeps reading as current
after the gap closes, and vice versa.

The deeper objection is the project's own recorded wisdom: rank impossible
above caught above remembered. A hand-maintained scope note is the
*remembered* tier, and it fails under exactly the load that stages R1 to R4
applied.

### Option C

Partly adopted -- every oracle test does now assert its fixture is
non-degenerate before trusting it, and those assertions stay. Rejected as
*sufficient* because they are local: each says "this fixture exercises what it
claims", and none says "every reachable cell is exercised by some fixture". The
union is the property that matters, and only a register can see it.

### Option D

Coverage is the metric that lied here in the first place. The cross-covariance
terms were at 100% line coverage while graded by nothing, because coverage
records that a line ran, not that its result was checked against anything. A
higher target would have made the lie more emphatic.

## Implemented option (B)

The register enumerates the schema's terms against the consumers of each schema
variant -- fifty-six cells at the time of writing. For each grading scenario it
builds the components that scenario's oracle actually uses and records, per
term, whether the value is present AND not identically zero. "Present but zero"
is precisely the state that reads as covered and grades nothing, so it counts
as ungraded.

A cell that is ungraded must appear in a declared list with a written reason,
or the build fails. Three further assertions keep the register from decaying:
a declared gap that has since been graded fails (a stale exemption is worse
than none), the term list is held against the documented schema so a new
component cannot silently escape tracking, and the register must be
non-trivially satisfied so it cannot pass vacuously.

**The register earned itself immediately.** On first run it measured 54 of 56
cells graded, with the two gaps being the online smoother and the causal
influence range against the noise cross-covariance -- both graded only by a
fixture reusing the dyad signal, whose cross-covariance is zero. That was
suspected and disclosed but never measured. A new harness closed both, and the
register now reads 56 of 56 with an empty exemption list.

## Forward cost

- **Every future capability must be graded or declared before it can ship.**
  Adding a consumer, a schema term, or a fixture that annihilates something
  breaks the build. That is the intended cost and the whole value.
- **The register is only as honest as its scenario list.** It reads the
  components each oracle actually grades, so it cannot be satisfied by
  intention -- but a scenario whose `consumers` field overstates what its
  oracle compares would overstate coverage. That field is the one place a
  false claim could enter, and it is small enough to review.
- **Non-degeneracy is defined as "not identically zero", which is a floor
  rather than a standard.** A term that is non-zero but negligible would count
  as graded. Tightening that would need a per-term magnitude threshold, and
  the current fixtures are far from the boundary.
- **An empty exemption list is a state to be maintained, not an achievement to
  be assumed.** It is empty today because a fixture was built to make it so.

## References

### Methodological

- The Independent Oracle Principle, and the RB22 finding it generated: a test
  whose oracle shares the claim's source proves self-consistency, not
  correspondence. Estate sweep, 2026-07: 11 of 21 oracle-bearing packages
  exposed to the same class.

### Empirical

- `aciR/tests/testthat/test-grading-matrix.R` -- the register.
- `tools/oracle/aci_oracle_cir_cross.m` -- the harness the register caused to exist.
- Measured 2026-08-13: 54 of 56 cells graded on first build; 56 of 56 after
  the cross-noise range fixture, with the exemption list empty.
- Measured 2026-08-13: the closed gap grades at 8.881784e-15 (online smoother
  mean), 1.352390e-14 (peak), 6.317863e-14 (objective range).

### Operational

- Cairn: `2026-08-13-vector-core-dispatch-and-cross-noise-gap.cairn.md` (this
  project) -- the trace defect that showed hand-tracking had already failed.
- `aciR/inst/extdata/oracle-manifest.yml` -- fixture `cir_cross`, whose
  `exists_because` field records that the register found the gap it closes.
